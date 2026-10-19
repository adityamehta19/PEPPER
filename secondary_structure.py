"""
secondary_structure.py
======================
Assigns secondary structure to every residue in a motif window.

Three methods are supported in priority order:

1. DSSP (primary)
-----------------
Runs the external ``mkdssp`` (DSSP 2.x, 3.x or 4.x) and reads its classic
output with ``Bio.PDB.DSSP.make_dssp_dict``.  Gold-standard hydrogen-bond-based
algorithm (Kabsch & Sander, 1983).
Install with: conda install -c conda-forge dssp  OR  sudo apt install dssp

mkdssp is looked for on PATH and in conda/mamba installations (e.g.
``<env>/Library/bin/mkdssp.exe`` on Windows).  mkdssp 4.x needs the libcifpp
data files (components.cif, mmcif_pdbx.dic); conda only points to them in an
activated environment, so the tool finds that folder next to the executable
and sets LIBCIFPP_DATA_DIR itself.  DSSP runs on the original file; if mkdssp
rejects it, or the chain IDs do not fit DSSP's one-letter column (e.g. "A-2"
in assembly files), it runs on a clean PDB copy of the model's protein.

2. VMD/STRIDE (secondary - highly recommended fallback)
--------------------------------------------------------
Uses VMD in headless batch mode to run the STRIDE algorithm
(Frishman & Argos, 1995), which is bundled inside every VMD installation.
STRIDE is hydrogen-bond-based like DSSP and gives very similar assignments.

VMD is invoked as:
    vmd -dispdev text -eofexit -e stride_script.tcl     (stdin closed)

The Tcl script loads a clean PDB copy of the model (VMD cannot be relied on to
read every mmCIF or multi-letter chain ID), calls ``mol ssrecalc`` (which runs
STRIDE), and writes the ``structure`` code of each residue to a file.  VMD is
run with a temporary folder whose path has no spaces: VMD passes that path to
STRIDE unquoted, and with a space STRIDE fails and VMD silently reports every
residue as coil.  That failure is detected and treated as an error.

Both programs are accepted only after they assign a built-in test helix
correctly, and the reason for rejecting one is logged (``--check_ss``).

STRIDE eight-state alphabet (same as DSSP):
    H — Alpha helix        G — 3-10 helix       I — Pi helix
    E — Extended strand    B — Isolated bridge   T — Turn
    C — Coil               b — Isolated bridge (lower-case in STRIDE output)

3. Bio.PDB dihedral fallback (last resort)
------------------------------------------
Ramachandran-region classification from backbone phi/psi angles.
Less accurate — cannot detect turns or bends. Used only when neither
DSSP nor VMD is available.

Reference
---------
Frishman, D. & Argos, P. (1995). Proteins, 23, 566-579.
Kabsch, W. & Sander, C. (1983). Biopolymers, 22, 2577-2637.
"""

from __future__ import annotations

import functools
import glob
import logging
import os
import re
import shutil
import string
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from Bio.PDB import is_aa
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

log = logging.getLogger(__name__)

# Eight-state to three-state mapping (shared by DSSP and STRIDE)
_SS_3STATE: dict[str, str] = {
    "H": "Helix",  "G": "Helix",  "I": "Helix",
    "E": "Strand", "B": "Strand", "b": "Strand",
    "T": "Coil",   "S": "Coil",   "C": "Coil",   "-": "Coil", "?": "Coil",
    "P": "Coil",
}

_SS_FULLNAME: dict[str, str] = {
    "H": "Alpha-helix",
    "G": "3-10-helix",
    "I": "Pi-helix",
    "E": "Beta-strand",
    "B": "Beta-bridge",
    "b": "Beta-bridge",
    "T": "Turn",
    "S": "Bend",
    "P": "PPII-helix",      # polyproline II, assigned by DSSP 4.x
    "C": "Coil",
    "-": "Coil",
    "?": "Coil",
}

_AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "U",
}


# ===========================================================================
# ResidueSSRecord dataclass
# ===========================================================================

@dataclass
class ResidueSSRecord:
    """Secondary structure record for one residue in a motif window."""
    position: int = 0
    residue_name: str = ""
    residue_number: str = ""
    one_letter: str = "X"
    ss_eight: str = "?"
    ss_three: str = "Coil"
    ss_fullname: str = "Coil"
    phi: float | None = None
    psi: float | None = None
    method: str = "unknown"


def _aa3to1(resname: str) -> str:
    return _AA3TO1.get(resname.strip().upper(), "X")


# ===========================================================================
# Shared plumbing for the external programs (mkdssp, VMD)
# ===========================================================================

# Explicit executable paths (--dssp_exe / --vmd_exe).  When set, ONLY that path is
# tried, so a wrong path is reported instead of silently using something else.
_EXE_OVERRIDE: dict[str, str | None] = {"dssp": None, "vmd": None}

# Why each candidate executable was rejected by the last detection (shown by
# --check_ss and logged), so "not found" and "found but broken" can be told apart.
_DETECT_ERRORS: dict[str, list[str]] = {"dssp": [], "vmd": []}

# Marker for "this program could not be run on that file" in the per-file caches.
_FAILED = object()

# Results per structure file: {(chain id, resseq, icode): (ss8, phi, psi)} or _FAILED.
# DSSP / VMD run once per file, not once per motif window.
_DSSP_CACHE: dict[tuple[str, object], object] = {}
_STRIDE_CACHE: dict[tuple[str, object], object] = {}

_DSSP_TIMEOUT = 600      # seconds per structure
_VMD_TIMEOUT = 600

# Windows "process could not start" exit codes, with a readable explanation.
_WINDOWS_EXIT_CODES = {
    0xC0000135: "Windows could not find a DLL that the program needs (exit code 0xC0000135)",
    0xC0000139: "a DLL the program needs is the wrong version (exit code 0xC0000139)",
    0xC000007B: "the program or one of its DLLs is for the wrong architecture (exit code 0xC000007B)",
}


def configure_executables(dssp_exe: str | None = None, vmd_exe: str | None = None) -> None:
    """Set explicit paths of the DSSP (mkdssp) and VMD executables, and reset detection."""
    _EXE_OVERRIDE["dssp"] = dssp_exe or None
    _EXE_OVERRIDE["vmd"] = vmd_exe or None
    _dssp_available.cache_clear()
    _vmd_available.cache_clear()


def clear_caches() -> None:
    """Forget cached DSSP / STRIDE results and tool detection (for tests)."""
    _DSSP_CACHE.clear()
    _STRIDE_CACHE.clear()
    _DSSP_TOOLS.clear()
    _dssp_available.cache_clear()
    _vmd_available.cache_clear()


def _resolve_exe(name: str) -> str | None:
    """Full path of *name* (a path, or a program name looked up on PATH), or None."""
    if os.path.isfile(name):
        return os.path.abspath(name)
    return shutil.which(name)


def _short(text: str, limit: int = 400) -> str:
    """The informative part of a program's output, on one line."""
    lines = [l.strip() for l in (text or "").splitlines() if l.strip()]
    joined = " | ".join(lines)
    return joined if len(joined) <= limit else "..." + joined[-limit:]


def _exit_reason(returncode: int) -> str:
    code = returncode & 0xFFFFFFFF
    if code in _WINDOWS_EXIT_CODES:
        return _WINDOWS_EXIT_CODES[code]
    return f"exit code {returncode}" if returncode else ""


def _helix_ok(codes: str) -> bool:
    """True when most residues of the test helix were assigned as helix."""
    return bool(codes) and sum(c in "HGI" for c in codes) >= 0.8 * len(codes)


def _run_on_test_helix(run) -> tuple[bool, str]:
    """Run *run(pdb_path, model)* on the built-in ideal helix and judge the result.

    *run* returns ``(residue map or None, error text)``.  Used during detection, so
    a program is only reported as available when it actually assigns a helix.
    """
    from Bio.PDB import PDBParser
    with tempfile.TemporaryDirectory() as tmp:
        pdb = os.path.join(tmp, "selftest_helix.pdb")
        with open(pdb, "w") as fh:
            fh.write(_SELFTEST_PDB)
        model = PDBParser(QUIET=True).get_structure("T", pdb)[0]
        result, err = run(pdb, model)
    if result is None:
        return False, err or "no result"
    codes = "".join(result.get(("A", n, " "), ("?",))[0] for n in range(5, 15))
    if not _helix_ok(codes):
        return False, f"gave {codes!r} for an ideal alpha-helix (expected mostly H)"
    return True, ""


# ---------------------------------------------------------------------------
# Clean single-model PDB copy, readable by every DSSP / STRIDE version
# ---------------------------------------------------------------------------

# Single-character chain IDs that DSSP and STRIDE both accept.
_CHAIN_LETTERS = string.ascii_uppercase + string.ascii_lowercase + string.digits


def _pdb_atom_name(name: str, element: str) -> str:
    """Atom name in PDB columns 13-16 (one-letter elements start in column 14)."""
    if len(name) >= 4 or len(element) == 2:
        return f"{name:<4}"[:4]
    return f" {name:<3}"


def _write_clean_pdbs(model, directory: str) -> list[tuple[str, dict[str, str]]] | None:
    """Write the protein residues of *model* as plain PDB files.

    The original file can be unreadable for the external programs (mmCIF for
    VMD and old DSSP, chain IDs longer than one character, PDB records that
    mkdssp 4 rejects, alternate locations).  This copy has only ATOM records of
    amino-acid residues (no hydrogens, first alternate location) with
    single-character chain IDs.

    Returns ``[(path, {chain letter in file: original chain id})]`` - one file per
    62 chains, because a PDB file cannot hold more chain IDs - or None when the
    model has no protein or a residue number does not fit the PDB format.
    """
    chains = []
    for chain in model:
        residues = [r for r in chain if is_aa(r) and "CA" in r]
        if residues:
            chains.append((chain.id, residues))
    if not chains:
        return None

    files: list[tuple[str, dict[str, str]]] = []
    batch_size = len(_CHAIN_LETTERS)
    for b in range(0, len(chains), batch_size):
        group = chains[b:b + batch_size]
        keep = {cid for cid, _ in group if len(cid) == 1 and cid in _CHAIN_LETTERS}
        spare = iter(c for c in _CHAIN_LETTERS if c not in keep)
        mapping: dict[str, str] = {}
        lines: list[str] = []
        serial = 0
        for cid, residues in group:
            letter = cid if cid in keep else next(spare)
            mapping[letter] = cid
            for res in residues:
                _, resseq, icode = res.id
                if not -999 <= resseq <= 9999:
                    log.debug("residue number %s does not fit a PDB file", resseq)
                    return None
                resname = res.get_resname().strip()[:3]
                for atom in res:            # a disordered atom yields its selected altloc
                    element = (atom.element or atom.get_id()[:1]).strip().upper()
                    if element in ("H", "D"):
                        continue
                    serial += 1
                    x, y, z = atom.coord
                    lines.append(
                        "ATOM  %5d %s %3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n"
                        % (serial % 100000, _pdb_atom_name(atom.get_id(), element), resname,
                           letter, resseq, icode.strip()[:1] or " ", x, y, z, 1.0, 0.0,
                           element[:2])
                    )
            # A full TER record: mkdssp 4 rejects a bare "TER" line and then loses
            # the first residue of the next chain.
            last = residues[-1]
            serial += 1
            lines.append("TER   %5d      %3s %1s%4d%1s\n" % (
                serial % 100000, last.get_resname().strip()[:3], letter, last.id[1],
                last.id[2].strip()[:1] or " "))
        lines.append("END\n")
        path = os.path.join(directory, f"clean_{len(files) + 1}.pdb")
        with open(path, "w") as fh:
            fh.writelines(lines)
        files.append((path, mapping))
    return files


def _protein_coverage(result: dict, model) -> float:
    """Fraction of the model's protein residues that have an entry in *result*."""
    total = found = 0
    for chain in model:
        for res in chain:
            if is_aa(res) and "CA" in res:
                total += 1
                found += (chain.id, res.id[1], res.id[2]) in result
    return found / total if total else 0.0


def _records_from_map(
    ss_map: dict, window_residues: list[Residue], method: str,
) -> list[ResidueSSRecord]:
    records: list[ResidueSSRecord] = []
    for pos, res in enumerate(window_residues, 1):
        chain_id = res.get_parent().id
        key = (chain_id, res.id[1], res.id[2])
        if key in ss_map:
            ss8, phi, psi = ss_map[key]
        else:
            ss8, phi, psi = "?", None, None
        records.append(ResidueSSRecord(
            position=pos, residue_name=res.get_resname().strip(),
            residue_number=f"{res.id[1]}{res.id[2].strip()}",
            one_letter=_aa3to1(res.get_resname()),
            ss_eight=ss8, ss_three=_SS_3STATE.get(ss8, "Coil"),
            ss_fullname=_SS_FULLNAME.get(ss8, "Coil"),
            phi=phi, psi=psi, method=method,
        ))
    return records


# ===========================================================================
# Method 1: DSSP
# ===========================================================================

@dataclass
class _DsspTool:
    path: str                    # full path of the executable
    version: tuple[int, ...]     # e.g. (4, 6, 1); decides the command line
    env: dict                    # environment it is run with


# Executables that passed detection, by the name _dssp_available() returned.
_DSSP_TOOLS: dict[str, _DsspTool] = {}


def _conda_prefixes() -> list[Path]:
    """Conda / mamba installations and their environments (DSSP is usually installed with conda)."""
    roots: list[Path] = []
    for var in ("CONDA_PREFIX", "MAMBA_ROOT_PREFIX", "CONDA_ROOT"):
        if os.environ.get(var):
            roots.append(Path(os.environ[var]))
    home = Path.home()
    for name in ("anaconda3", "miniconda3", "miniforge3", "mambaforge", "micromamba", ".conda"):
        roots.append(home / name)
    if os.name == "nt":
        for base in (os.environ.get("ProgramData", r"C:\ProgramData"), r"C:\\"):
            for name in ("anaconda3", "miniconda3", "miniforge3", "mambaforge"):
                roots.append(Path(base) / name)
    prefixes: list[Path] = []
    for root in roots:
        if root.is_dir():
            prefixes.append(root)
            envs = root / "envs"
            if envs.is_dir():
                prefixes.extend(sorted(p for p in envs.iterdir() if p.is_dir()))
    return prefixes


def _dssp_candidates() -> list[str]:
    explicit = _EXE_OVERRIDE["dssp"]
    if explicit:
        return [explicit]
    candidates = ["mkdssp", "dssp", "mkdssp.exe", "dssp.exe"]
    for prefix in _conda_prefixes():
        for rel in ("Library/bin/mkdssp.exe", "bin/mkdssp", "Scripts/mkdssp.exe", "bin/mkdssp.exe"):
            p = prefix / rel
            if p.is_file():
                candidates.append(str(p))
    return candidates


def _libcifpp_data_dir(exe_path: str) -> str | None:
    """Folder with the data files mkdssp 4 needs (components.cif, mmcif_pdbx.dic).

    mkdssp 4.x (libcifpp) refuses to read any structure without them. conda sets
    LIBCIFPP_DATA_DIR only inside an *activated* environment, so a mkdssp run
    from outside it fails on every file.  The folder is therefore looked for
    next to the executable as well.
    """
    def has_data(d: Path) -> int:
        if (d / "components.cif").exists() or (d / "components.cif.gz").exists():
            return 2
        return 1 if (d / "mmcif_pdbx.dic").exists() else 0

    current = os.environ.get("LIBCIFPP_DATA_DIR")
    if current and Path(current).is_dir() and has_data(Path(current)):
        return current
    exe_dir = Path(exe_path).resolve().parent
    dirs = [
        exe_dir.parent / "share" / "libcifpp",                 # <prefix>/bin/mkdssp
        exe_dir.parent.parent / "share" / "libcifpp",          # <prefix>\Library\bin\mkdssp.exe (conda, Windows)
        exe_dir.parent / "Library" / "share" / "libcifpp",
        exe_dir.parent.parent / "Library" / "share" / "libcifpp",
    ]
    if os.environ.get("CONDA_PREFIX"):
        dirs.append(Path(os.environ["CONDA_PREFIX"]) / "share" / "libcifpp")
    dirs += [Path(p) for p in ("/usr/share/libcifpp", "/usr/local/share/libcifpp",
                               "/opt/homebrew/share/libcifpp", "/var/cache/libcifpp")]
    scored = [(has_data(d), i, d) for i, d in enumerate(dirs) if d.is_dir()]
    scored = [s for s in scored if s[0]]
    if not scored:
        return None
    best = max(scored, key=lambda s: (s[0], -s[1]))
    return str(best[2])


def _dssp_env(exe_path: str) -> dict:
    """Environment for running mkdssp: its data folder, and its DLL folder on PATH."""
    env = os.environ.copy()
    data_dir = _libcifpp_data_dir(exe_path)
    if data_dir:
        env["LIBCIFPP_DATA_DIR"] = data_dir
    exe_dir = Path(exe_path).resolve().parent
    extra = [str(exe_dir)]
    if exe_dir.name.lower() == "bin" and exe_dir.parent.name.lower() == "library":
        extra.append(str(exe_dir.parent.parent))              # conda prefix (Windows DLLs)
    env["PATH"] = os.pathsep.join(extra + [env.get("PATH", "")])
    return env


def _run_dssp(tool: _DsspTool, in_file: str | Path) -> tuple[dict | None, str]:
    """Run DSSP on *in_file*.

    Returns ``({(chain, resseq, icode): (ss8, phi, psi)}, "")``, or
    ``(None, reason)`` if DSSP failed.
    """
    from Bio.PDB.DSSP import make_dssp_dict
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "out.dssp")
        if tool.version >= (4,):
            cmd = [tool.path, "--output-format", "dssp", str(in_file), out]
        else:                                    # DSSP 2.x / 3.x
            cmd = [tool.path, "-i", str(in_file), "-o", out]
        try:
            proc = subprocess.run(
                cmd, capture_output=True, text=True, encoding="utf-8", errors="replace",
                timeout=_DSSP_TIMEOUT, env=tool.env, stdin=subprocess.DEVNULL,
            )
        except subprocess.TimeoutExpired:
            return None, f"timed out after {_DSSP_TIMEOUT} s"
        except OSError as exc:
            return None, str(exc)
        if not os.path.exists(out) or os.path.getsize(out) == 0:
            reason = " - ".join(x for x in (_exit_reason(proc.returncode),
                                            _short(proc.stderr or proc.stdout)) if x)
            return None, reason or "DSSP wrote no output"
        try:
            dssp_dict, _ = make_dssp_dict(out)
        except Exception as exc:
            return None, f"could not read the DSSP output: {exc}"
    if not dssp_dict:
        return None, "the DSSP output has no residues"
    result = {}
    for (chain, res_id), entry in dssp_dict.items():
        ss8 = entry[1] if entry[1].strip() else "-"
        result[(chain, res_id[1], res_id[2])] = (ss8, _angle_or_none(entry[3]), _angle_or_none(entry[4]))
    return result, ""


def _probe_dssp(candidate: str) -> _DsspTool | None:
    """Return the tool if *candidate* runs and assigns a test helix correctly."""
    errors = _DETECT_ERRORS["dssp"]
    path = _resolve_exe(candidate)
    if path is None:
        errors.append(f"{candidate}: not found")
        return None
    env = _dssp_env(path)
    try:
        proc = subprocess.run(
            [path, "--version"], capture_output=True, text=True, encoding="utf-8",
            errors="replace", timeout=60, env=env, stdin=subprocess.DEVNULL,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        errors.append(f"{path}: could not be started ({exc})")
        return None
    text = (proc.stdout or "") + (proc.stderr or "")
    match = re.search(r"(\d+)\.(\d+)(?:\.(\d+))?", text)
    if proc.returncode != 0 and not match:
        errors.append(f"{path}: " + (" - ".join(
            x for x in (_exit_reason(proc.returncode), _short(text)) if x) or "failed to start"))
        return None
    version = tuple(int(g) for g in match.groups() if g) if match else (4,)
    tool = _DsspTool(path, version, env)
    ok, err = _run_on_test_helix(lambda pdb, model: _run_dssp(tool, pdb))
    if not ok:
        hint = ""
        if version >= (4,) and "LIBCIFPP_DATA_DIR" not in env:
            hint = (" (mkdssp 4 needs its data files components.cif / mmcif_pdbx.dic; set "
                    "LIBCIFPP_DATA_DIR to the folder that holds them, e.g. <conda env>/share/libcifpp)")
        errors.append(f"{path} (version {'.'.join(map(str, version))}): {err}{hint}")
        return None
    return tool


@functools.lru_cache(maxsize=1)
def _dssp_available() -> str | None:
    """Return the name or path of a working DSSP executable, or None.

    Candidates: the path given with ``--dssp_exe`` (alone, if given), otherwise
    ``mkdssp`` / ``dssp`` on PATH and mkdssp in conda installations.  A
    candidate is accepted only if it starts AND assigns a built-in test helix,
    so a broken install (missing DLLs, missing libcifpp data files) is reported
    with the reason instead of failing silently on every structure.  The result
    is cached, so the probe runs once per session.
    """
    _DETECT_ERRORS["dssp"] = []
    seen: set[str] = set()
    for candidate in _dssp_candidates():
        path = _resolve_exe(candidate)
        if path is not None:
            if path.lower() in seen:
                continue
            seen.add(path.lower())
        tool = _probe_dssp(candidate)
        if tool is not None:
            _DSSP_TOOLS[candidate] = tool
            log.debug("DSSP: using %s (version %s, LIBCIFPP_DATA_DIR=%s)", tool.path,
                      ".".join(map(str, tool.version)), tool.env.get("LIBCIFPP_DATA_DIR"))
            return candidate
    broken = [e for e in _DETECT_ERRORS["dssp"] if not e.endswith(": not found")]
    if broken:
        for e in broken:
            log.warning("DSSP found but not usable: %s", e)
    elif _EXE_OVERRIDE["dssp"]:
        log.warning("--dssp_exe %r could not be run; DSSP will not be used", _EXE_OVERRIDE["dssp"])
    return None


def _dssp_tool(dssp_exe: str) -> _DsspTool:
    tool = _DSSP_TOOLS.get(dssp_exe)
    if tool is None:                       # called directly with a name (tests, self-test)
        path = _resolve_exe(dssp_exe) or dssp_exe
        tool = _DsspTool(path, (4,), _dssp_env(path) if os.path.exists(path) else os.environ.copy())
        _DSSP_TOOLS[dssp_exe] = tool
    return tool


def _dssp_map_for(model, pdb_file: str | Path, dssp_exe: str):
    """DSSP residue map for the first model of *pdb_file* (cached per file)."""
    key = (str(Path(pdb_file).resolve()), model.id)
    if key in _DSSP_CACHE:
        return _DSSP_CACHE[key]
    tool = _dssp_tool(dssp_exe)
    result: object = _FAILED
    reason = ""

    # 1) the original file, as long as its chain IDs fit DSSP's one-letter column
    single_letter = all(len(c.id) == 1 and c.id.strip() for c in model)
    old_and_cif = tool.version < (4,) and Path(pdb_file).suffix.lower() in (".cif", ".mmcif")
    if single_letter and not old_and_cif:
        raw, reason = _run_dssp(tool, pdb_file)
        if raw is not None and _protein_coverage(raw, model) >= 0.9:
            result = raw
        elif raw is not None:
            reason = "DSSP residues did not match the parsed structure"

    # 2) a clean PDB copy of the protein (multi-letter chains, records DSSP rejects, ...)
    if result is _FAILED:
        with tempfile.TemporaryDirectory() as tmp:
            files = _write_clean_pdbs(model, tmp)
            merged: dict | None = {} if files else None
            for path, mapping in files or []:
                raw, err = _run_dssp(tool, path)
                if raw is None:
                    merged, reason = None, err
                    break
                for (letter, resseq, icode), value in raw.items():
                    if letter in mapping:
                        merged[(mapping[letter], resseq, icode)] = value
            if merged:
                result = merged
                if reason:
                    log.debug("DSSP could not use %s directly (%s); used a clean copy",
                              pdb_file, reason)
            elif files is None:
                reason = reason or "no protein residues that can be written to a PDB file"

    if result is _FAILED:
        log.warning("DSSP failed on %s: %s", Path(pdb_file).name, reason or "unknown error")
    _DSSP_CACHE[key] = result
    return result


def _angle_or_none(value) -> float | None:
    """DSSP reports 360.0 for an undefined phi/psi; convert to None."""
    if value is None:
        return None
    v = float(value)
    if abs(v - 360.0) < 1e-6:
        return None
    return round(v, 2)


def _assign_dssp(
    structure: Structure,
    pdb_file: str | Path,
    window_residues: list[Residue],
    dssp_exe: str = "mkdssp",
) -> list[ResidueSSRecord] | None:
    """Assign secondary structure with DSSP using the executable *dssp_exe*."""
    model = window_residues[0].get_parent().get_parent()
    dssp_map = _dssp_map_for(model, pdb_file, dssp_exe)
    if dssp_map is _FAILED:
        return None
    return _records_from_map(dssp_map, window_residues, "DSSP")


# ===========================================================================
# Method 2: VMD/STRIDE
# ===========================================================================

# Usual VMD install locations on Windows (tried only if the file exists).  The
# VMD 1.9.4 installer defaults to "C:\Program Files\VMD"; older ones used the
# "University of Illinois" folder.
_VMD_WINDOWS_PATHS = (
    r"C:\Program Files\VMD\vmd.exe",
    r"C:\Program Files (x86)\VMD\vmd.exe",
    r"C:\Program Files\University of Illinois\VMD\vmd.exe",
    r"C:\Program Files (x86)\University of Illinois\VMD\vmd.exe",
)


def _vmd_candidates() -> list[str]:
    explicit = _EXE_OVERRIDE["vmd"]
    if explicit:
        return [explicit]
    candidates = ["vmd", "vmd.exe", "vmd_LINUXAMD64", "vmd_WIN64.exe", "vmd_MACOSXX86_64",
                  "vmd_MACOSXARM64"]
    candidates += [p for p in _VMD_WINDOWS_PATHS if os.path.exists(p)]
    if os.name == "nt":
        for pattern in (r"C:\Program Files*\VMD*\vmd.exe",
                        r"C:\Program Files*\University of Illinois\VMD*\vmd.exe"):
            candidates += sorted(glob.glob(pattern))
    if os.environ.get("VMDDIR"):
        candidates += [os.path.join(os.environ["VMDDIR"], n) for n in ("vmd.exe", "vmd")]
    candidates += [p for p in ("/usr/local/bin/vmd", "/opt/vmd/bin/vmd") if os.path.exists(p)]
    return candidates


def _space_free_tmpdir() -> str | None:
    """A writable temporary folder whose path has no spaces, or None.

    VMD hands STRIDE a file in the temporary folder on its command line without
    quoting it. With a space in that path (e.g. C:\\Users\\First Last\\AppData\\
    Local\\Temp) STRIDE cannot read its input and VMD silently reports every
    residue as coil.
    """
    base = tempfile.gettempdir()
    if " " not in base:
        return base
    options = []
    if os.name == "nt":
        try:
            import ctypes
            buf = ctypes.create_unicode_buffer(1024)
            if ctypes.windll.kernel32.GetShortPathNameW(base, buf, 1024):
                options.append(buf.value)            # 8.3 form, e.g. C:\Users\FIRSTL~1\...
        except Exception:
            pass
        options += [os.path.join(os.environ.get("SystemRoot", r"C:\Windows"), "Temp"),
                    os.path.join(os.environ.get("PUBLIC", r"C:\Users\Public"), "Documents"),
                    r"C:\Temp"]
    else:
        options += ["/tmp", "/var/tmp"]
    for d in options:
        if " " in d or not os.path.isdir(d):
            continue
        try:
            probe = tempfile.mkdtemp(prefix="mmx_", dir=d)
            os.rmdir(probe)
            return d
        except OSError:
            continue
    return None


# VMD Tcl script: load the clean PDB, run STRIDE (mol ssrecalc), write one line
# per residue.  Both paths are in a space-free temporary folder.
_VMD_TCL_TEMPLATE = """\
mol new {{{pdb_file}}} type pdb waitfor all
set molid [molinfo top]
if {{[catch {{mol ssrecalc $molid}} err]}} {{
    puts "MMX_STRIDE_ERROR $err"
}}
set sel [atomselect $molid "name CA"]
set fh [open {{{out_file}}} w]
foreach a [$sel get {{chain resid insertion resname structure phi psi}}] {{
    puts $fh [join $a "\\t"]
}}
close $fh
quit
"""

# VMD's STRIDE codes (H G I E B T C) are already in the DSSP alphabet.
_VMD_TO_SS8: dict[str, str] = {
    "H": "H", "G": "G", "I": "I",
    "E": "E", "B": "B", "b": "B",
    "T": "T", "S": "S", "C": "C", "": "C",
}

# Text VMD prints when its STRIDE call failed (it then reports all residues as coil).
_VMD_STRIDE_FAILURES = ("unable to read output file from stride", "mmx_stride_error",
                        "unable to call stride", "stride failed")


def _run_vmd_stride(model, vmd_exe: str) -> tuple[dict | None, str]:
    """Run STRIDE through VMD on a clean PDB copy of *model*.

    Returns ``({(chain, resseq, icode): (ss8, phi, psi)}, "")`` or ``(None, reason)``.
    """
    tmp_root = _space_free_tmpdir()
    if tmp_root is None:
        return None, ("no temporary folder without spaces in its path was found; VMD's "
                      "STRIDE cannot work there (set TMP to a folder such as C:\\Temp)")
    result: dict = {}
    with tempfile.TemporaryDirectory(prefix="mmx_vmd_", dir=tmp_root) as tmp:
        files = _write_clean_pdbs(model, tmp)
        if not files:
            return None, "no protein residues that can be written to a PDB file"
        env = os.environ.copy()
        for var in ("TMP", "TEMP", "TMPDIR"):
            env[var] = tmp
        for i, (pdb_path, mapping) in enumerate(files, 1):
            out_path = os.path.join(tmp, f"stride_{i}.txt")
            tcl_path = os.path.join(tmp, f"stride_{i}.tcl")
            with open(tcl_path, "w") as fh:
                # Tcl accepts forward slashes on every platform; backslashes of
                # Windows paths would be read as escapes.
                fh.write(_VMD_TCL_TEMPLATE.format(
                    pdb_file=Path(pdb_path).as_posix(), out_file=Path(out_path).as_posix()))
            # The script is passed with -e and stdin is closed: the Windows build
            # of VMD does not read commands piped through stdin.
            try:
                proc = subprocess.run(
                    [vmd_exe, "-dispdev", "text", "-eofexit", "-e", Path(tcl_path).as_posix()],
                    stdin=subprocess.DEVNULL, capture_output=True, text=True,
                    encoding="utf-8", errors="replace", timeout=_VMD_TIMEOUT, env=env, cwd=tmp,
                )
            except subprocess.TimeoutExpired:
                return None, f"timed out after {_VMD_TIMEOUT} s"
            except OSError as exc:
                return None, str(exc)
            output = (proc.stdout or "") + (proc.stderr or "")
            failure = next((l for l in output.splitlines()
                            if any(f in l.lower() for f in _VMD_STRIDE_FAILURES)), None)
            if failure:
                return None, f"VMD's STRIDE run failed: {failure.strip()}"
            if not os.path.exists(out_path):
                return None, "VMD wrote no output - " + (_short(output) or _exit_reason(proc.returncode))
            with open(out_path, encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    parts = line.rstrip("\r\n").split("\t")
                    if len(parts) < 7:
                        continue
                    letter, resid_str, icode, _resname, ss_raw, phi_str, psi_str = parts[:7]
                    letter = letter.strip()
                    try:
                        resseq = int(resid_str)
                    except ValueError:
                        continue
                    if letter not in mapping:
                        continue
                    result[(mapping[letter], resseq, icode.strip() or " ")] = (
                        _VMD_TO_SS8.get(ss_raw.strip(), "C"),
                        _vmd_angle(phi_str), _vmd_angle(psi_str),
                    )
    if not result:
        return None, "VMD reported no residues"
    return result, ""


def _vmd_angle(text: str) -> float | None:
    """VMD reports 0.0 for an undefined phi/psi (chain ends); convert to None."""
    try:
        v = float(text)
    except ValueError:
        return None
    return None if v == 0.0 else round(v, 2)


@functools.lru_cache(maxsize=1)
def _vmd_available() -> str | None:
    """Return the VMD executable (name or path) if it runs STRIDE correctly, else None.

    Tried: the path given with ``--vmd_exe`` (alone, if given), otherwise
    ``vmd`` on PATH, the platform-specific binary names, and the usual install
    folders.  Each candidate is checked by running STRIDE on a built-in helix
    (VMD is started with a script file, never with commands on stdin, which the
    Windows build ignores).
    """
    _DETECT_ERRORS["vmd"] = []
    seen: set[str] = set()
    for candidate in _vmd_candidates():
        path = _resolve_exe(candidate)
        if path is None:
            _DETECT_ERRORS["vmd"].append(f"{candidate}: not found")
            continue
        if path.lower() in seen:
            continue
        seen.add(path.lower())
        ok, err = _run_on_test_helix(lambda pdb, model: _run_vmd_stride(model, path))
        if ok:
            log.debug("VMD/STRIDE: using %s", path)
            return candidate if candidate == _EXE_OVERRIDE["vmd"] else path
        _DETECT_ERRORS["vmd"].append(f"{path}: {err}")
    broken = [e for e in _DETECT_ERRORS["vmd"] if not e.endswith(": not found")]
    if broken:
        for e in broken:
            log.warning("VMD found but STRIDE is not usable: %s", e)
    elif _EXE_OVERRIDE["vmd"]:
        log.warning("--vmd_exe %r could not be run; VMD/STRIDE will not be used", _EXE_OVERRIDE["vmd"])
    return None


def _vmd_stride_map(model, pdb_file: str | Path, vmd_exe: str):
    """Run VMD/STRIDE once for the model of a file and return the residue map (cached)."""
    key = (str(Path(pdb_file).resolve()), model.id)
    if key in _STRIDE_CACHE:
        return _STRIDE_CACHE[key]
    result, reason = _run_vmd_stride(model, vmd_exe)
    if result is None:
        log.warning("VMD/STRIDE failed on %s: %s", Path(pdb_file).name, reason)
        result = _FAILED
    _STRIDE_CACHE[key] = result
    return result


def _assign_vmd_stride(
    pdb_file: str | Path,
    window_residues: list[Residue],
    vmd_exe: str,
) -> list[ResidueSSRecord] | None:
    """Assign secondary structure using VMD's built-in STRIDE algorithm.

    VMD is run in headless batch mode (-dispdev text) with a Tcl script that
    loads a clean PDB copy of the model, triggers STRIDE via ``mol ssrecalc``,
    and writes per-residue structure codes to a temp file.  The result is
    cached per file, so VMD runs once per structure.

    Returns
    -------
    list[ResidueSSRecord] or None
        Per-residue records, or None if VMD fails.
    """
    model = window_residues[0].get_parent().get_parent()
    stride_map = _vmd_stride_map(model, pdb_file, vmd_exe)
    if stride_map is _FAILED:
        return None
    records = _records_from_map(stride_map, window_residues, "VMD_STRIDE")
    log.debug("VMD/STRIDE assigned %d residues in motif window", len(records))
    return records


# ===========================================================================
# Method 3: Bio.PDB dihedral fallback
# ===========================================================================

def _dihedral_angle(p1, p2, p3, p4):
    b1 = p2 - p1; b2 = p3 - p2; b3 = p4 - p3
    n1 = np.cross(b1, b2); n2 = np.cross(b2, b3)
    n1n = np.linalg.norm(n1); n2n = np.linalg.norm(n2)
    if n1n < 1e-6 or n2n < 1e-6:
        return float("nan")
    n1 /= n1n; n2 /= n2n
    m1 = np.cross(n1, b2 / (np.linalg.norm(b2) + 1e-12))
    return float(np.degrees(np.arctan2(np.dot(m1, n2), np.dot(n1, n2))))


def _get_phi_psi(chain_residues: list[Residue]) -> dict[int, tuple]:
    result: dict[int, tuple] = {}
    n = len(chain_residues)

    def _atom(res, name):
        return res[name].get_vector().get_array() if name in res else None

    for i, res in enumerate(chain_residues):
        phi = psi = None
        if i > 0:
            c_prev = _atom(chain_residues[i-1], "C")
            n_cur  = _atom(res, "N")
            ca_cur = _atom(res, "CA")
            c_cur  = _atom(res, "C")
            if all(x is not None for x in [c_prev, n_cur, ca_cur, c_cur]):
                v = _dihedral_angle(c_prev, n_cur, ca_cur, c_cur)
                phi = round(v, 2) if not np.isnan(v) else None
                phi = -phi

        if i < n - 1:
            n_cur  = _atom(res, "N")
            ca_cur = _atom(res, "CA")
            c_cur  = _atom(res, "C")
            n_next = _atom(chain_residues[i+1], "N")
            if all(x is not None for x in [n_cur, ca_cur, c_cur, n_next]):
                v = _dihedral_angle(n_cur, ca_cur, c_cur, n_next)
                psi = round(v, 2) if not np.isnan(v) else None
                psi = -psi

        result[res.id[1]] = (phi, psi)
    return result


def _ss_from_dihedrals(phi, psi):
    if phi is None or psi is None:
        return "?", "Coil", "Coil"
    if -130 <= phi <= -40 and -70 <= psi <= 20:
        return "H", "Helix", "Alpha-helix"
    if -170 <= phi <= -50 and (70 <= psi <= 180 or -180 <= psi <= -130):
        return "E", "Strand", "Beta-strand"
    if -100 <= phi <= -30 and -30 <= psi <= 10:
        return "G", "Helix", "3-10-helix"
    return "C", "Coil", "Coil"


def _assign_biopdb(
    structure: Structure,
    window_residues: list[Residue],
) -> list[ResidueSSRecord]:
    chain_id = window_residues[0].get_parent().id if window_residues else "A"
    chain_residues: list[Residue] = []
    # Use the model the window residues belong to (not always model 0).
    window_model = window_residues[0].get_parent().get_parent() if window_residues else None
    if window_model is not None:
        for chain in window_model:
            if chain.id == chain_id:
                chain_residues = sorted(
                    [r for r in chain if is_aa(r)],
                    key=lambda r: (r.id[1], r.id[2]),
                )
    # phi/psi of a whole chain is computed once and kept on the chain object;
    # without this, every motif window of a long chain recomputed it.
    window_chain = window_residues[0].get_parent()
    phi_psi = getattr(window_chain, "_mmx_phi_psi", None)
    if phi_psi is None:
        phi_psi = _get_phi_psi(chain_residues)
        try:
            window_chain._mmx_phi_psi = phi_psi
        except AttributeError:
            pass
    records: list[ResidueSSRecord] = []
    for pos, res in enumerate(window_residues, 1):
        one_letter = _aa3to1(res.get_resname())
        resnum = f"{res.id[1]}{res.id[2].strip()}"
        phi, psi = phi_psi.get(res.id[1], (None, None))
        ss8, ss3, ssname = _ss_from_dihedrals(phi, psi)
        records.append(ResidueSSRecord(
            position=pos, residue_name=res.get_resname().strip(),
            residue_number=resnum, one_letter=one_letter,
            ss_eight=ss8, ss_three=ss3, ss_fullname=ssname,
            phi=phi, psi=psi, method="BioPDB_dihedral",
        ))
    log.debug("BioPDB dihedral fallback assigned %d residues", len(records))
    return records


# ===========================================================================
# Public entry point
# ===========================================================================

def assign_secondary_structure(
    structure: Structure,
    pdb_file: str | Path,
    window_residues: list[Residue],
    motif_id: str,
) -> list[ResidueSSRecord]:
    """Assign secondary structure to motif window residues.

    Priority order:
    1. DSSP  (mkdssp binary)           — gold standard, hydrogen-bond-based
    2. VMD/STRIDE (vmd binary)         — equally accurate, ships with VMD
    3. Bio.PDB dihedral                — last resort, less accurate

    The method used is recorded in each ResidueSSRecord.method field
    and reported in both output CSV files.

    Parameters
    ----------
    structure : Bio.PDB.Structure
        Full parsed parent structure (needed for correct H-bond context).
    pdb_file : str or Path
        Path to the full PDB file.
    window_residues : list[Residue]
        Motif window residues.
    motif_id : str
        Motif identifier for logging.

    Returns
    -------
    list[ResidueSSRecord]
        One record per residue in the window.
    """
    if not window_residues:
        return []

    # DSSP and STRIDE (run on the file) analyse the first model of a file.
    # Windows from any other model therefore use the dihedral method, so
    # that a window is never labelled with another model's secondary
    # structure.
    window_model = window_residues[0].get_parent().get_parent()
    first_model = structure.get_list()[0]
    if window_model is not first_model:
        log.debug(
            "%s: window is in model %s (not the first model) - using the "
            "dihedral method", motif_id, window_model.id,
        )
        return _assign_biopdb(structure, window_residues)

    # --- Try DSSP first ---
    dssp_exe = _dssp_available()
    if dssp_exe:
        records = _assign_dssp(structure, pdb_file, window_residues, dssp_exe)
        if records is not None:
            log.debug("%s: secondary structure assigned via DSSP (%s)", motif_id, dssp_exe)
            return records
        log.debug("%s: DSSP failed for this structure - trying VMD/STRIDE", motif_id)

    # --- Try VMD/STRIDE ---
    vmd_exe = _vmd_available()
    if vmd_exe:
        records = _assign_vmd_stride(pdb_file, window_residues, vmd_exe)
        if records is not None:
            log.debug("%s: secondary structure assigned via VMD/STRIDE", motif_id)
            return records
        log.debug(
            "%s: VMD/STRIDE failed for this structure - falling back to dihedral method",
            motif_id,
        )
    else:
        log.debug(
            "%s: VMD not found - using Bio.PDB dihedral fallback. "
            "Install VMD (https://www.ks.uiuc.edu/Research/vmd/) for "
            "accurate STRIDE-based assignment without needing mkdssp.",
            motif_id,
        )

    # --- Last resort: backbone dihedrals ---
    return _assign_biopdb(structure, window_residues)


# ===========================================================================
# Self-test: which secondary-structure methods work on this machine?
# ===========================================================================

_SELFTEST_PDB = """HEADER    IDEAL ALPHA HELIX (SELF-TEST)           01-JAN-00   TEST
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1           1
ATOM      1  N   ALA A   1      -0.525   1.362   0.000  1.00  0.00           N  
ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  
ATOM      3  C   ALA A   1       1.520   0.000   0.000  1.00  0.00           C  
ATOM      4  O   ALA A   1       2.144  -0.723   0.775  1.00  0.00           O  
ATOM      5  CB  ALA A   1      -0.507  -0.774  -1.206  1.00  0.00           C  
ATOM      6  N   ALA A   2       2.116   0.811  -0.869  1.00  0.00           N  
ATOM      7  CA  ALA A   2       3.571   0.893  -0.958  1.00  0.00           C  
ATOM      8  C   ALA A   2       4.180   1.269   0.384  1.00  0.00           C  
ATOM      9  O   ALA A   2       5.144   0.647   0.828  1.00  0.00           O  
ATOM     10  CB  ALA A   2       3.992   1.908  -2.008  1.00  0.00           C  
ATOM     11  N   ALA A   3       3.617   2.286   1.030  1.00  0.00           N  
ATOM     12  CA  ALA A   3       4.118   2.734   2.326  1.00  0.00           C  
ATOM     13  C   ALA A   3       4.112   1.599   3.336  1.00  0.00           C  
ATOM     14  O   ALA A   3       5.096   1.385   4.042  1.00  0.00           O  
ATOM     15  CB  ALA A   3       3.286   3.890   2.856  1.00  0.00           C  
ATOM     16  N   ALA A   4       3.001   0.870   3.407  1.00  0.00           N  
ATOM     17  CA  ALA A   4       2.882  -0.246   4.341  1.00  0.00           C  
ATOM     18  C   ALA A   4       3.987  -1.266   4.118  1.00  0.00           C  
ATOM     19  O   ALA A   4       4.625  -1.713   5.069  1.00  0.00           O  
ATOM     20  CB  ALA A   4       1.530  -0.924   4.202  1.00  0.00           C  
ATOM     21  N   ALA A   5       4.214  -1.634   2.860  1.00  0.00           N  
ATOM     22  CA  ALA A   5       5.250  -2.607   2.528  1.00  0.00           C  
ATOM     23  C   ALA A   5       6.609  -2.153   3.036  1.00  0.00           C  
ATOM     24  O   ALA A   5       7.339  -2.932   3.646  1.00  0.00           O  
ATOM     25  CB  ALA A   5       5.319  -2.830   1.026  1.00  0.00           C  
ATOM     26  N   ALA A   6       6.949  -0.892   2.786  1.00  0.00           N  
ATOM     27  CA  ALA A   6       8.230  -0.347   3.227  1.00  0.00           C  
ATOM     28  C   ALA A   6       8.394  -0.488   4.731  1.00  0.00           C  
ATOM     29  O   ALA A   6       9.441  -0.926   5.207  1.00  0.00           O  
ATOM     30  CB  ALA A   6       8.357   1.116   2.837  1.00  0.00           C  
ATOM     31  N   ALA A   7       7.360  -0.118   5.482  1.00  0.00           N  
ATOM     32  CA  ALA A   7       7.406  -0.209   6.938  1.00  0.00           C  
ATOM     33  C   ALA A   7       7.705  -1.631   7.387  1.00  0.00           C  
ATOM     34  O   ALA A   7       8.557  -1.847   8.247  1.00  0.00           O  
ATOM     35  CB  ALA A   7       6.093   0.250   7.549  1.00  0.00           C  
ATOM     36  N   ALA A   8       7.005  -2.600   6.804  1.00  0.00           N  
ATOM     37  CA  ALA A   8       7.208  -4.002   7.157  1.00  0.00           C  
ATOM     38  C   ALA A   8       8.657  -4.412   6.953  1.00  0.00           C  
ATOM     39  O   ALA A   8       9.256  -5.043   7.824  1.00  0.00           O  
ATOM     40  CB  ALA A   8       6.304  -4.902   6.331  1.00  0.00           C  
ATOM     41  N   ALA A   9       9.222  -4.055   5.804  1.00  0.00           N  
ATOM     42  CA  ALA A   9      10.609  -4.395   5.500  1.00  0.00           C  
ATOM     43  C   ALA A   9      11.550  -3.854   6.564  1.00  0.00           C  
ATOM     44  O   ALA A   9      12.426  -4.571   7.046  1.00  0.00           O  
ATOM     45  CB  ALA A   9      11.013  -3.850   4.140  1.00  0.00           C  
ATOM     46  N   ALA A  10      11.370  -2.588   6.931  1.00  0.00           N  
ATOM     47  CA  ALA A  10      12.215  -1.963   7.944  1.00  0.00           C  
ATOM     48  C   ALA A  10      12.171  -2.742   9.249  1.00  0.00           C  
ATOM     49  O   ALA A  10      13.210  -3.017   9.847  1.00  0.00           O  
ATOM     50  CB  ALA A  10      11.784  -0.528   8.196  1.00  0.00           C  
ATOM     51  N   ALA A  11      10.968  -3.098   9.691  1.00  0.00           N  
ATOM     52  CA  ALA A  11      10.805  -3.849  10.932  1.00  0.00           C  
ATOM     53  C   ALA A  11      11.589  -5.150  10.891  1.00  0.00           C  
ATOM     54  O   ALA A  11      12.298  -5.483  11.840  1.00  0.00           O  
ATOM     55  CB  ALA A  11       9.338  -4.149  11.193  1.00  0.00           C  
ATOM     56  N   ALA A  12      11.464  -5.888   9.791  1.00  0.00           N  
ATOM     57  CA  ALA A  12      12.171  -7.156   9.641  1.00  0.00           C  
ATOM     58  C   ALA A  12      13.671  -6.968   9.797  1.00  0.00           C  
ATOM     59  O   ALA A  12      14.324  -7.723  10.515  1.00  0.00           O  
ATOM     60  CB  ALA A  12      11.876  -7.782   8.288  1.00  0.00           C  
ATOM     61  N   ALA A  13      14.218  -5.959   9.124  1.00  0.00           N  
ATOM     62  CA  ALA A  13      15.650  -5.684   9.198  1.00  0.00           C  
ATOM     63  C   ALA A  13      16.089  -5.457  10.636  1.00  0.00           C  
ATOM     64  O   ALA A  13      17.090  -6.018  11.079  1.00  0.00           O  
ATOM     65  CB  ALA A  13      16.010  -4.468   8.360  1.00  0.00           C  
ATOM     66  N   ALA A  14      15.340  -4.635  11.365  1.00  0.00           N  
ATOM     67  CA  ALA A  14      15.666  -4.344  12.758  1.00  0.00           C  
ATOM     68  C   ALA A  14      15.734  -5.619  13.582  1.00  0.00           C  
ATOM     69  O   ALA A  14      16.678  -5.820  14.345  1.00  0.00           O  
ATOM     70  CB  ALA A  14      14.641  -3.401  13.367  1.00  0.00           C  
ATOM     71  N   ALA A  15      14.734  -6.482  13.429  1.00  0.00           N  
ATOM     72  CA  ALA A  15      14.695  -7.740  14.169  1.00  0.00           C  
ATOM     73  C   ALA A  15      15.945  -8.565  13.910  1.00  0.00           C  
ATOM     74  O   ALA A  15      16.559  -9.078  14.845  1.00  0.00           O  
ATOM     75  CB  ALA A  15      13.465  -8.550  13.792  1.00  0.00           C  
ATOM     76  N   ALA A  16      16.322  -8.694  12.642  1.00  0.00           N  
ATOM     77  CA  ALA A  16      17.507  -9.465  12.275  1.00  0.00           C  
ATOM     78  C   ALA A  16      18.744  -8.934  12.982  1.00  0.00           C  
ATOM     79  O   ALA A  16      19.523  -9.705  13.541  1.00  0.00           O  
ATOM     80  CB  ALA A  16      17.729  -9.432  10.772  1.00  0.00           C  
ATOM     81  N   ALA A  17      18.925  -7.616  12.959  1.00  0.00           N  
ATOM     82  CA  ALA A  17      20.078  -6.996  13.605  1.00  0.00           C  
ATOM     83  C   ALA A  17      20.137  -7.356  15.081  1.00  0.00           C  
ATOM     84  O   ALA A  17      21.192  -7.732  15.590  1.00  0.00           O  
ATOM     85  CB  ALA A  17      20.032  -5.484  13.454  1.00  0.00           C  
ATOM     86  N   ALA A  18      19.004  -7.243  15.768  1.00  0.00           N  
ATOM     87  CA  ALA A  18      18.943  -7.562  17.192  1.00  0.00           C  
ATOM     88  C   ALA A  18      19.398  -8.989  17.452  1.00  0.00           C  
ATOM     89  O   ALA A  18      19.415  -9.440  18.596  1.00  0.00           O  
ATOM     90  CB  ALA A  18      17.533  -7.374  17.727  1.00  0.00           C  
END
"""

_INSTALL_HINT = (
    "How to get it:\n"
    "  DSSP (mkdssp)  Windows/macOS/Linux with conda:  conda install -c conda-forge dssp\n"
    "                 Windows without conda: use WSL (Ubuntu) and  sudo apt install dssp\n"
    "                 Ubuntu/Debian:  sudo apt install dssp\n"
    "                 then, if it is not on PATH:  --dssp_exe <full path to mkdssp>\n"
    "                 (conda on Windows: <env>\\Library\\bin\\mkdssp.exe; its data folder\n"
    "                 <env>\\share\\libcifpp is found automatically, or set LIBCIFPP_DATA_DIR)\n"
    "  VMD/STRIDE     install VMD (https://www.ks.uiuc.edu/Research/vmd/), then, if it is\n"
    "                 not on PATH:  --vmd_exe <full path to vmd>  (Windows: C:\\Program Files\\VMD\\vmd.exe)"
)


def available_methods() -> str:
    """One line naming the secondary-structure method the run will use (runs detection)."""
    dssp_exe = _dssp_available()
    if dssp_exe:
        tool = _DSSP_TOOLS.get(dssp_exe)
        version = ".".join(map(str, tool.version)) if tool else "?"
        return f"DSSP ({tool.path if tool else dssp_exe}, version {version})"
    vmd_exe = _vmd_available()
    if vmd_exe:
        return f"VMD/STRIDE ({vmd_exe}); DSSP not usable"
    return "dihedral fallback only - DSSP and VMD/STRIDE not usable (run --check_ss)"


def _problems(tool: str) -> list[str]:
    """Detection errors of executables that exist but do not work."""
    return [e for e in _DETECT_ERRORS[tool] if not e.endswith(": not found")]


def self_test(verbose: bool = True) -> dict:
    """Run every secondary-structure method on a small ideal alpha-helix.

    Returns a dict ``{method: {"available": bool, "ok": bool, "codes": str}}``
    for ``DSSP``, ``VMD_STRIDE`` and ``BioPDB_dihedral``; a method is ``ok`` when
    most interior residues of the helix are assigned as helix.
    """
    import tempfile
    from Bio.PDB import PDBParser

    results: dict[str, dict] = {}
    with tempfile.TemporaryDirectory() as tmp:
        pdb = os.path.join(tmp, "selftest_helix.pdb")
        with open(pdb, "w") as fh:
            fh.write(_SELFTEST_PDB)
        structure = PDBParser(QUIET=True).get_structure("T", pdb)
        residues = [r for r in structure[0]["A"] if is_aa(r)]
        window = residues[4:14]                      # interior residues

        def judge(recs):
            codes = "".join(r.ss_eight for r in recs) if recs else ""
            ok = bool(recs) and sum(c in "HGI" for c in codes) >= 0.8 * len(codes)
            return codes, ok

        clear_caches()
        dssp_exe = _dssp_available()
        recs = _assign_dssp(structure, pdb, window, dssp_exe) if dssp_exe else None
        codes, ok = judge(recs)
        results["DSSP"] = {"available": bool(dssp_exe), "executable": dssp_exe,
                           "ok": ok, "codes": codes, "problems": _problems("dssp")}

        vmd_exe = _vmd_available()
        recs = _assign_vmd_stride(pdb, window, vmd_exe) if vmd_exe else None
        codes, ok = judge(recs)
        results["VMD_STRIDE"] = {"available": bool(vmd_exe), "executable": vmd_exe,
                                 "ok": ok, "codes": codes, "problems": _problems("vmd")}

        recs = _assign_biopdb(structure, window)
        codes, ok = judge(recs)
        results["BioPDB_dihedral"] = {"available": True, "executable": "built in",
                                      "ok": ok, "codes": codes, "problems": []}
        clear_caches()

    if verbose:
        print("Secondary-structure self-test (ideal alpha-helix; expected: H = helix)")
        for name, r in results.items():
            if not r["available"] and r["problems"]:
                status = "FOUND BUT NOT WORKING"
            elif not r["available"]:
                status = "NOT FOUND"
            elif r["ok"]:
                status = "works"
            else:
                status = f"found but gave an unexpected result ({r['codes'] or 'no output'})"
            exe = f" [{r['executable']}]" if r["available"] and r["executable"] != "built in" else ""
            print(f"  {name:<16}: {status}{exe}" + (f"  codes: {r['codes']}" if r["ok"] else ""))
            if not r["available"]:
                for problem in r["problems"]:
                    print(f"      reason: {problem}")
        best = next((n for n in ("DSSP", "VMD_STRIDE") if results[n]["ok"]), None)
        if best:
            print(f"\nThe tool will use {best} (priority: DSSP, then VMD/STRIDE, then dihedral).")
        else:
            print("\nNeither DSSP nor VMD/STRIDE works here, so the tool falls back to the "
                  "less accurate dihedral method.\n" + _INSTALL_HINT)
    return results
