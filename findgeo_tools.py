"""
findgeo_tools.py
================
Tools for comparing MetalMotifExtractor's geometry assignment with FindGeo
(Andreini, Cavallaro & Lorenzini, Bioinformatics 28, 1658-1660, 2012).

Two comparisons are supported, as suggested by the reviewer:

A. **Same donor atoms** (tests the geometry classifier alone).
   MetalMotifExtractor writes, for every motif, a small PDB file with the
   metal and exactly the donor atoms it counted (``write_donor_pdb``).
   FindGeo is run on these files with a threshold large enough to keep every
   atom, so both programs see identical atoms.

B. **Each program defines its own coordination sphere** (tests the whole
   pipeline).  FindGeo is run on the full structure with its own sphere
   definition, and sites are matched by ion identity.  Differences in
   coordination number, donor atoms and geometry are tabulated, so the effect
   of water handling and distance cutoffs is visible.

FindGeo facts this module relies on (from the FindGeo help page, checked
when this module was written; verify them against the version you run)
------------------------------------------------------------------------
* A donor is any atom within a threshold distance (default 2.8 A) of the
  metal, except excluded elements (default C and H).
* For superposition FindGeo sets all metal-donor distances to 3.0 A, so its
  RMSD is on a 3.0 A scale.  MetalMotifExtractor's RMSD uses the mean
  metal-donor distance; use ``rmsd_at_3A`` here to compare numbers.
* Output per metal, in a directory named ``<METAL>_<resnum>_<serial>_<chain>``:
  ``findgeo.input`` (PDB lines; the metal first, then the donor atoms) and
  ``findgeo.out`` (RMSD for every geometry tested, and the best geometry).
* The library has 36 geometries, including "with a vacancy" shapes.

Geometry names
--------------
FindGeo uses three-letter codes.  ``FINDGEO_TO_OURS`` maps the codes whose
shape is also in MetalMotifExtractor's library.  Five of these (trv, tev,
spv, pva, bvp) are vacancy shapes whose equivalence was verified numerically
(vertex removal gives an exact fit, RMSD 0); the others are mapped by
description and assume the two programs' ideal coordinates are the same
shape.  Geometries present in only one library cannot agree and are reported
as "not comparable".
"""

from __future__ import annotations

import logging
import math
import os
import re
import shlex
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# FindGeo library (36 geometries) and the mapping to this tool's library
# ---------------------------------------------------------------------------
FINDGEO_GEOMETRIES: dict[str, tuple[int, str]] = {
    "lin": (2, "Linear"), "trv": (2, "Trigonal plane with a vacancy"),
    "tri": (3, "Trigonal plane"), "tev": (3, "Tetrahedron with a vacancy"),
    "spv": (3, "Square plane with a vacancy"),
    "tet": (4, "Tetrahedron"), "spl": (4, "Square plane"),
    "bva": (4, "Trigonal bipyramid with a vacancy (axial)"),
    "bvp": (4, "Trigonal bipyramid with a vacancy (equatorial)"),
    "pyv": (4, "Square pyramid with a vacancy (equatorial)"),
    "spy": (5, "Square pyramid"), "tbp": (5, "Trigonal bipyramid"),
    "tpv": (5, "Trigonal prism with a vacancy"),
    "oct": (6, "Octahedron"), "tpr": (6, "Trigonal prism"),
    "pva": (6, "Pentagonal bipyramid with a vacancy (axial)"),
    "pvp": (6, "Pentagonal bipyramid with a vacancy (equatorial)"),
    "cof": (6, "Octahedron, face monocapped with a vacancy (capped face)"),
    "con": (6, "Octahedron, face monocapped with a vacancy (non-capped face)"),
    "ctf": (6, "Trigonal prism, square-face monocapped with a vacancy (capped face)"),
    "ctn": (6, "Trigonal prism, square-face monocapped with a vacancy (non-capped edge)"),
    "pbp": (7, "Pentagonal bipyramid"), "coc": (7, "Octahedron, face monocapped"),
    "ctp": (7, "Trigonal prism, square-face monocapped"),
    "hva": (7, "Hexagonal bipyramid with a vacancy (axial)"),
    "hvp": (7, "Hexagonal bipyramid with a vacancy (equatorial)"),
    "cuv": (7, "Cube with a vacancy"), "sav": (7, "Square antiprism with a vacancy"),
    "hbp": (8, "Hexagonal bipyramid"), "cub": (8, "Cube"),
    "sqa": (8, "Square antiprism"), "boc": (8, "Octahedron, trans-bicapped"),
    "bts": (8, "Trigonal prism, square-face bicapped"),
    "btt": (8, "Trigonal prism, triangular-face bicapped"),
    "ttp": (9, "Trigonal prism, square-face tricapped"),
    "csa": (9, "Square antiprism, square-face monocapped"),
}

# FindGeo code -> name used in ideal_polyhedra.py
FINDGEO_TO_OURS: dict[str, str] = {
    "lin": "Linear", "trv": "Angular (120 deg)",
    "tri": "Trigonal planar", "tev": "Trigonal pyramidal", "spv": "T-shaped",
    "tet": "Tetrahedral", "spl": "Square planar", "bvp": "Seesaw",
    "spy": "Square pyramidal", "tbp": "Trigonal bipyramidal",
    "oct": "Octahedral", "tpr": "Trigonal prismatic", "pva": "Pentagonal pyramidal",
    "pbp": "Pentagonal bipyramidal", "coc": "Capped octahedral",
    "ctp": "Capped trigonal prismatic (square-face cap)",
    "hbp": "Hexagonal bipyramidal", "cub": "Cubic", "sqa": "Square antiprismatic",
    "bts": "Bicapped trigonal prismatic (square-face caps)",
    "btt": "Bicapped trigonal prismatic (triangular-face caps)",
    "ttp": "Tricapped trigonal prismatic", "csa": "Capped square antiprismatic",
}
OURS_TO_FINDGEO = {v: k for k, v in FINDGEO_TO_OURS.items()}

_METAL_ELEMENTS = {
    "LI", "NA", "K", "MG", "CA", "SR", "BA", "MN", "FE", "CO", "NI", "CU", "ZN",
    "CD", "HG", "LA", "CE", "ND", "SM", "EU", "GD", "TB", "YB", "LU",
}


# ===========================================================================
# 1. Writing the donor-atom PDB for FindGeo (comparison A)
# ===========================================================================

def _split_resnum(resnum: str) -> tuple[int, str]:
    m = re.match(r"^(-?\d+)([A-Za-z]?)$", str(resnum).strip())
    if not m:
        return 0, " "
    num = int(m.group(1))
    if abs(num) > 9999:           # PDB residue number field is 4 characters
        num = num % 10000
    return num, (m.group(2) or " ")


def _pdb_line(record, serial, name, resname, chain, resnum, xyz, element) -> str:
    num, icode = _split_resnum(resnum)
    element = (element or "").strip().upper()
    name = name.strip()
    atom_name = name.ljust(4) if (len(element) >= 2 or len(name) >= 4) else (" " + name).ljust(4)
    return (
        f"{record:<6s}{serial:5d} {atom_name} {resname:>3s} {(chain or ' ')[:1]}"
        f"{num:4d}{icode}   {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
        f"{1.0:6.2f}{0.0:6.2f}          {element:>2s}"
    )


def write_donor_pdb(path: str | Path, ion, donors) -> Path:
    """Write the metal and exactly the donor atoms counted for a site.

    The metal is the first line, then the donor atoms in the order of the
    donor list.  Running FindGeo on this file with a large enough threshold
    (e.g. 4.5 A) makes FindGeo see precisely the atoms MetalMotifExtractor
    used.  Chain IDs longer than one character are truncated to one.

    Parameters
    ----------
    ion : Bio.PDB Residue
        The metal ion residue.
    donors : list[coordination.DonorAtom]
        The donor atoms (with identity fields filled in).
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    atom = next(iter(ion.get_atoms()))
    element = (atom.element or ion.get_resname()).strip()
    lines = [_pdb_line(
        "HETATM", 1, atom.get_name(), ion.get_resname().strip(), ion.get_parent().id,
        f"{ion.id[1]}{ion.id[2].strip()}", atom.get_coord(), element,
    )]
    for i, d in enumerate(donors, start=2):
        lines.append(_pdb_line(
            "HETATM" if d.hetero else "ATOM", i, d.atom_name, d.resname, d.chain,
            d.resnum, d.coord, d.element,
        ))
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")
    return path


# ===========================================================================
# 2. Parsing FindGeo output
# ===========================================================================

@dataclass
class FindGeoAtom:
    name: str
    resname: str
    chain: str
    resnum: str
    element: str
    xyz: tuple[float, float, float]

    @property
    def key(self) -> tuple[str, str, str, str]:
        return (self.resname, self.chain, self.resnum, self.name)


@dataclass
class FindGeoSite:
    group: str                      # first directory below the run root (motif ID or PDB ID)
    directory: str                  # e.g. ZN_262_2041_A
    metal: str = ""
    resnum: str = ""
    serial: str = ""
    chain: str = ""
    atoms: list[FindGeoAtom] = field(default_factory=list)   # donors only (metal excluded)
    best: str = ""                  # three-letter code
    best_tag: str = ""              # "regular" / "irregular"
    rmsd: dict[str, float] = field(default_factory=dict)

    @property
    def cn(self) -> int:
        return len(self.atoms)

    @property
    def best_rmsd(self) -> float | None:
        return self.rmsd.get(self.best)


_OUT_ROW = re.compile(r"^\s*(\w{3})\s+-\s+(.*?)\s*\|\s*(\w+)\s*\|\s*([-\d.]+)\s*$")
_OUT_BEST = re.compile(r"Best geometry of the site:\s*(\w{3})\s*\(\s*(\w+)\s*\)", re.I)
_DIR = re.compile(r"^([A-Za-z]{1,2})_(-?\d+[A-Za-z]?)_(\d+)_(.+)$")


def parse_findgeo_out(text: str) -> dict:
    """Parse the text of a ``findgeo.out`` file.

    Returns ``{"rmsd": {code: float}, "tags": {code: str}, "best": code,
    "best_tag": str}``; ``best`` is "" if FindGeo reported none.
    """
    rmsd: dict[str, float] = {}
    tags: dict[str, str] = {}
    for line in text.splitlines():
        m = _OUT_ROW.match(line)
        if m:
            rmsd[m.group(1).lower()] = float(m.group(4))
            tags[m.group(1).lower()] = m.group(3).lower()
    best, best_tag = "", ""
    m = _OUT_BEST.search(text)
    if m:
        best, best_tag = m.group(1).lower(), m.group(2).lower()
    return {"rmsd": rmsd, "tags": tags, "best": best, "best_tag": best_tag}


def parse_findgeo_input(text: str) -> tuple[FindGeoAtom | None, list[FindGeoAtom]]:
    """Parse a ``findgeo.input`` file: (metal, donor atoms)."""
    atoms: list[FindGeoAtom] = []
    for line in text.splitlines():
        if not line.startswith(("ATOM", "HETATM")) or len(line) < 54:
            continue
        name = line[12:16].strip()
        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        if not element:
            element = re.sub(r"[^A-Za-z]", "", name)[:1].upper()
        atoms.append(FindGeoAtom(
            name=name, resname=line[17:20].strip(), chain=line[21].strip(),
            resnum=f"{line[22:26].strip()}{line[26].strip()}", element=element,
            xyz=(float(line[30:38]), float(line[38:46]), float(line[46:54])),
        ))
    if not atoms:
        return None, []
    return atoms[0], atoms[1:]


def collect_runs(run_root: str | Path) -> list[FindGeoSite]:
    """Read every FindGeo result below *run_root*.

    Expected layout: ``run_root/<group>/<METAL_resnum_serial_chain>/`` with
    ``findgeo.out`` and ``findgeo.input``.  ``<group>`` is the motif ID
    (comparison A) or the PDB ID (comparison B).
    """
    run_root = Path(run_root)
    sites: list[FindGeoSite] = []
    for out_file in sorted(run_root.rglob("findgeo.out")):
        rel = out_file.parent.relative_to(run_root).parts
        if len(rel) < 2:
            continue
        group, directory = rel[0], rel[1]
        site = FindGeoSite(group=group, directory=directory)
        m = _DIR.match(directory)
        if m:
            site.metal, site.resnum, site.serial, site.chain = m.groups()
        parsed = parse_findgeo_out(out_file.read_text(errors="replace"))
        site.rmsd, site.best, site.best_tag = parsed["rmsd"], parsed["best"], parsed["best_tag"]
        inp = out_file.parent / "findgeo.input"
        if inp.exists():
            metal, donors = parse_findgeo_input(inp.read_text(errors="replace"))
            site.atoms = donors
            if metal is not None and not site.metal:
                site.metal = metal.element
        sites.append(site)
    return sites


# ===========================================================================
# 3. Running FindGeo
# ===========================================================================

def split_command(cmd: str, windows: bool | None = None) -> list[str]:
    """Split a command string such as ``python C:\\Tools\\findgeo.py`` into arguments.

    POSIX-style splitting treats a backslash as an escape character and would
    turn ``C:\\Tools\\findgeo.py`` into ``C:Toolsfindgeo.py``.  On Windows
    (or with ``windows=True``) backslashes are therefore kept, and quotes
    around an argument (needed for paths with spaces) are removed.
    """
    if windows is None:
        windows = os.name == "nt"
    if not windows:
        return shlex.split(cmd)
    parts = shlex.split(cmd, posix=False)
    return [p[1:-1] if len(p) >= 2 and p[0] == p[-1] and p[0] in "\"'" else p for p in parts]


def _run_one(cmd: list[str], pdb: Path, workdir: Path, threshold: float | None,
             extra: list[str]) -> tuple[str, bool, str]:
    workdir.mkdir(parents=True, exist_ok=True)
    local = workdir / pdb.name
    if not local.exists():
        shutil.copy(pdb, local)
    args = cmd + ["-p", local.name, "-w", str(workdir), "-o"]
    if threshold is not None:
        args += ["-t", str(threshold)]
    args += extra
    try:
        proc = subprocess.run(args, capture_output=True, text=True, timeout=300)
    except Exception as exc:                       # noqa: BLE001
        return pdb.stem, False, str(exc)
    return pdb.stem, proc.returncode == 0, (proc.stderr or proc.stdout)[-300:]


def run_findgeo(
    findgeo_cmd: str,
    pdb_dir: str | Path,
    out_root: str | Path,
    threshold: float | None,
    workers: int = 4,
    extra_args: list[str] | None = None,
    name_upper: bool = False,
) -> dict[str, bool]:
    """Run the stand-alone FindGeo once per ``*.pdb`` file in *pdb_dir*.

    Each input is copied to ``out_root/<name>/`` and FindGeo is run there
    (``findgeo.py -p <file> -w <dir> -o [-t threshold]``), so its per-metal
    output directories end up below ``out_root/<name>/``.

    Parameters
    ----------
    findgeo_cmd : str
        Command that starts FindGeo, e.g. ``"python /path/to/findgeo.py"``.
        (Check the options of the version you installed.)
    threshold : float or None
        FindGeo ``-t`` value.  Use a value >= the largest donor distance for
        comparison A (e.g. 4.5), FindGeo's default 2.8 for comparison B, or
        None to leave FindGeo's default.
    name_upper : bool
        Upper-case the run directory name (for PDB IDs).
    """
    cmd = split_command(findgeo_cmd)
    pdb_dir, out_root = Path(pdb_dir), Path(out_root)
    files = sorted(pdb_dir.glob("*.pdb"))
    results: dict[str, bool] = {}

    def job(f: Path):
        name = f.stem.upper() if name_upper else f.stem
        return _run_one(cmd, f, out_root / name, threshold, extra_args or [])

    with ThreadPoolExecutor(max_workers=workers) as pool:
        for name, ok, msg in pool.map(job, files):
            results[name] = ok
            if not ok:
                log.warning("FindGeo failed for %s: %s", name, msg)
    log.info("FindGeo finished: %d/%d runs without error", sum(results.values()), len(results))
    return results


# ===========================================================================
# 4. Comparison
# ===========================================================================

def cohens_kappa(a: list[str], b: list[str]) -> float | None:
    """Cohen's kappa for two equal-length label lists (None if undefined)."""
    n = len(a)
    if n == 0:
        return None
    po = sum(x == y for x, y in zip(a, b)) / n
    cats = set(a) | set(b)
    pe = sum((a.count(c) / n) * (b.count(c) / n) for c in cats)
    if abs(1 - pe) < 1e-12:
        return None
    return (po - pe) / (1 - pe)


def _status(ours: str, fg_code: str) -> tuple[str, str]:
    """Classify a geometry comparison; return (status, FindGeo name in our naming)."""
    fg_name = FINDGEO_TO_OURS.get(fg_code, "")
    if not fg_code:
        return "findgeo_no_result", ""
    if ours == "Irregular":
        return "ours_irregular", fg_name
    ours_comparable = ours in OURS_TO_FINDGEO
    fg_comparable = fg_code in FINDGEO_TO_OURS
    if ours_comparable and fg_comparable:
        return ("agree" if ours == fg_name else "disagree"), fg_name
    if not ours_comparable and not fg_comparable:
        return "neither_in_other_library", fg_name
    if not ours_comparable:
        return "ours_not_in_findgeo_library", fg_name
    return "findgeo_geometry_not_in_our_library", fg_name


def _rmsd_at_3A(rmsd_a, mean_dist) -> float | None:
    try:
        if pd.isna(rmsd_a) or pd.isna(mean_dist) or float(mean_dist) <= 0:
            return None
        return round(float(rmsd_a) * 3.0 / float(mean_dist), 3)
    except (TypeError, ValueError):
        return None


_PROTEIN_LABEL = re.compile(r"^([A-Za-z0-9]{3})(-?\d+)([A-Za-z]?):(.+)$")


def parse_our_donor(label: str, protein_chain: str) -> tuple[str, str, str, str] | None:
    """Convert an ``All Donor Atoms`` entry to (resname, chain, resnum, atom).

    ``ASP8:OD1`` (protein; chain = the record's chain), ``HOH:W:101:O`` and
    ``ATP:L:501:O2G`` (hetero; chain given in the label).
    """
    parts = label.split(":")
    if len(parts) == 4:
        return parts[0], parts[1], parts[2], parts[3]
    m = _PROTEIN_LABEL.match(label)
    if m:
        return m.group(1), protein_chain, f"{m.group(2)}{m.group(3)}", m.group(4)
    return None


def _our_donor_set(row, chain_col="Chain") -> set[tuple[str, str, str, str]]:
    labels = [x for x in str(row.get("All Donor Atoms", "") or "").split(",") if x]
    out = set()
    for lab in labels:
        t = parse_our_donor(lab, str(row.get(chain_col, "")))
        if t:
            out.add(t)
    return out


def _summarise(df: pd.DataFrame) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    """Overall numbers, per-CN table and confusion matrix for a site table."""
    have = df[df["status"] != "no_findgeo_site"]
    comp = have[have["status"].isin(["agree", "disagree"])]
    n_agree = int((comp["status"] == "agree").sum())
    overall = {
        "records": int(len(df)),
        "with_findgeo_result": int(len(have)),
        "comparable": int(len(comp)),
        "agree": n_agree,
        "agree_pct_of_comparable": round(100 * n_agree / len(comp), 1) if len(comp) else None,
        "cohens_kappa_comparable": (
            None if comp.empty else
            (lambda k: None if k is None else round(k, 3))(
                cohens_kappa(list(comp["our_geometry"]), list(comp["findgeo_geometry_in_our_names"])))
        ),
        "agree_within_our_top2_pct_of_comparable": (
            round(100 * float(comp["agree_top2"].mean()), 1) if len(comp) else None),
        "not_comparable": int(len(have) - len(comp)),
    }
    for st, n in have["status"].value_counts().items():
        overall[f"status_{st}"] = int(n)

    rows = []
    for cn, g in have.groupby("our_cn"):
        c = g[g["status"].isin(["agree", "disagree"])]
        a = int((c["status"] == "agree").sum())
        rows.append({
            "CN (ours)": cn, "records": len(g), "comparable": len(c), "agree": a,
            "agree_pct_of_comparable": round(100 * a / len(c), 1) if len(c) else None,
            "not_comparable": len(g) - len(c),
        })
    per_cn = pd.DataFrame(rows)
    confusion = (
        pd.crosstab(comp["our_geometry"], comp["findgeo_geometry_in_our_names"])
        if len(comp) else pd.DataFrame()
    )
    return overall, per_cn, confusion


def compare_same_donors(motifs: pd.DataFrame, sites: list[FindGeoSite]) -> pd.DataFrame:
    """Comparison A: FindGeo run on the donor-atom PDB written for each motif.

    Rows of *motifs* are matched to FindGeo results by motif ID (the run
    directory name).  Returns one row per motif.
    """
    by_group: dict[str, list[FindGeoSite]] = {}
    for s in sites:
        by_group.setdefault(s.group, []).append(s)

    rows = []
    for _, r in motifs.iterrows():
        mid = str(r["Motif ID"])
        found = by_group.get(mid, [])
        base = {
            "motif_id": mid, "pdb_id": r.get("PDB ID"), "metal": r.get("Metal"),
            "our_cn": r.get("Coordination Number"), "our_geometry": r.get("Coordination Geometry"),
            "our_second_best": r.get("Geometry Second Best"),
            "our_rmsd_A": r.get("Geometry RMSD (A)"), "our_margin_A": r.get("Geometry Margin (A)"),
            "our_rmsd_at_3A": _rmsd_at_3A(r.get("Geometry RMSD (A)"), r.get("Mean Donor Distance (A)")),
        }
        if not found:
            rows.append({**base, "status": "no_findgeo_site"})
            continue
        s = found[0]
        status, fg_name = _status(str(base["our_geometry"]), s.best)
        rows.append({
            **base,
            "findgeo_cn": s.cn, "cn_match": s.cn == base["our_cn"],
            "findgeo_geometry_code": s.best,
            "findgeo_geometry": FINDGEO_GEOMETRIES.get(s.best, (None, ""))[1],
            "findgeo_geometry_in_our_names": fg_name,
            "findgeo_rmsd_3A": s.best_rmsd, "findgeo_tag": s.best_tag,
            "status": status,
            "agree_top2": status == "agree" or (
                status == "disagree" and fg_name == str(base["our_second_best"])),
        })
    df = pd.DataFrame(rows)
    if "agree_top2" not in df:
        df["agree_top2"] = False
    return df


def compare_own_spheres(motifs: pd.DataFrame, sites: list[FindGeoSite]) -> pd.DataFrame:
    """Comparison B: FindGeo on the full structure (its own sphere definition).

    Each motif is matched to the FindGeo result for the same ion (PDB ID, then
    atom serial, else chain + residue number + metal).  Columns report the
    coordination number, the donor-atom overlap (precision/recall of our donor
    atoms against FindGeo's) and the geometry agreement, with a flag for
    whether the two coordination spheres are identical.
    """
    by_pdb: dict[str, list[FindGeoSite]] = {}
    for s in sites:
        by_pdb.setdefault(s.group.upper(), []).append(s)

    rows = []
    for _, r in motifs.iterrows():
        pdb = str(r["PDB ID"]).upper()
        metal = str(r.get("Metal", "")).upper()
        cands = [s for s in by_pdb.get(pdb, []) if s.metal.upper() == metal]
        serial = str(r.get("Ion Serial", ""))
        match = [s for s in cands if serial and s.serial == serial]
        if not match:
            ion_chain, ion_num = str(r.get("Ion Chain", "")), str(r.get("Ion Residue Number", ""))
            match = [s for s in cands if s.chain == ion_chain and s.resnum == ion_num]
        base = {
            "motif_id": r["Motif ID"], "pdb_id": pdb, "metal": metal,
            "ion": f"{r.get('Ion Chain', '')}{r.get('Ion Residue Number', '')}",
            "our_cn": r.get("Coordination Number"), "our_geometry": r.get("Coordination Geometry"),
            "our_second_best": r.get("Geometry Second Best"),
            "our_rmsd_at_3A": _rmsd_at_3A(r.get("Geometry RMSD (A)"), r.get("Mean Donor Distance (A)")),
        }
        if not match:
            rows.append({**base, "status": "no_findgeo_site"})
            continue
        s = match[0]
        ours = _our_donor_set(r)
        theirs = {a.key for a in s.atoms}
        tp, fp, fn = len(ours & theirs), len(ours - theirs), len(theirs - ours)
        status, fg_name = _status(str(base["our_geometry"]), s.best)
        extra_fg = sorted(f"{k[0]}:{k[1]}:{k[2]}:{k[3]}" for k in theirs - ours)
        extra_ours = sorted(f"{k[0]}:{k[1]}:{k[2]}:{k[3]}" for k in ours - theirs)
        rows.append({
            **base,
            "findgeo_cn": s.cn, "cn_match": s.cn == base["our_cn"],
            "same_sphere": ours == theirs,
            "donor_precision": round(tp / (tp + fp), 3) if (tp + fp) else None,
            "donor_recall": round(tp / (tp + fn), 3) if (tp + fn) else None,
            "donors_only_in_findgeo": ";".join(extra_fg),
            "donors_only_in_ours": ";".join(extra_ours),
            "findgeo_geometry_code": s.best,
            "findgeo_geometry": FINDGEO_GEOMETRIES.get(s.best, (None, ""))[1],
            "findgeo_geometry_in_our_names": fg_name,
            "findgeo_rmsd_3A": s.best_rmsd, "findgeo_tag": s.best_tag,
            "status": status,
            "agree_top2": status == "agree" or (
                status == "disagree" and fg_name == str(base["our_second_best"])),
        })
    df = pd.DataFrame(rows)
    if "agree_top2" not in df:
        df["agree_top2"] = False
    return df


def summarise_comparison(df: pd.DataFrame, only_identical_spheres: bool = False):
    """Overall summary, per-CN table and confusion matrix.

    With ``only_identical_spheres=True`` (comparison B) only records whose
    two coordination spheres are identical are used, which isolates the
    difference between the geometry classifiers.
    """
    d = df
    if only_identical_spheres and "same_sphere" in df:
        d = df[df["same_sphere"] == True]          # noqa: E712
    return _summarise(d)
