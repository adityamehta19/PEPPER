"""
secondary_structure.py
======================
Assigns secondary structure to every residue in a motif window.

Three methods are supported in priority order:

1. DSSP (primary)
-----------------
Calls the external ``mkdssp`` binary via ``Bio.PDB.DSSP``.
Gold-standard hydrogen-bond-based algorithm (Kabsch & Sander, 1983).
Install with: sudo apt install dssp  OR  conda install -c salilab dssp

2. VMD/STRIDE (secondary — highly recommended fallback)
--------------------------------------------------------
Uses VMD in headless batch mode to run the STRIDE algorithm
(Frishman & Argos, 1995), which is bundled inside every VMD installation.
STRIDE is hydrogen-bond-based like DSSP and produces equivalent accuracy.
Because STRIDE ships with VMD, any lab that has VMD installed gets
accurate secondary structure assignment without any additional software.

VMD is invoked as:
    vmd -dispdev text -eofexit < stride_script.tcl

The Tcl script loads the PDB, calls ``mol ssrecalc`` (which invokes
STRIDE internally), reads the ``structure`` attribute of each residue,
and writes a plain-text result file.

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

import logging
import os
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
# Method 1: DSSP
# ===========================================================================

def _dssp_available() -> bool:
    for binary in ("mkdssp", "dssp"):
        try:
            subprocess.run([binary, "--version"], capture_output=True, timeout=5)
            return True
        except (FileNotFoundError, subprocess.TimeoutExpired):
            continue
    return False


def _assign_dssp(
    structure: Structure,
    pdb_file: str | Path,
    window_residues: list[Residue],
) -> list[ResidueSSRecord] | None:
    try:
        from Bio.PDB.DSSP import DSSP
        dssp = DSSP(structure[0], str(pdb_file))
        dssp_dict = dict(dssp)
        records: list[ResidueSSRecord] = []
        for pos, res in enumerate(window_residues, 1):
            chain_id = res.get_parent().id
            key = (chain_id, res.id)
            one_letter = _aa3to1(res.get_resname())
            resnum = f"{res.id[1]}{res.id[2].strip()}"
            if key in dssp_dict:
                entry = dssp_dict[key]
                ss8 = entry[2] if entry[2] != " " else "C"
                phi = round(float(entry[4]), 2) if entry[4] is not None else None
                psi = round(float(entry[5]), 2) if entry[5] is not None else None
            else:
                ss8, phi, psi = "?", None, None
            ss3    = _SS_3STATE.get(ss8, "Coil")
            ssname = _SS_FULLNAME.get(ss8, "Coil")
            records.append(ResidueSSRecord(
                position=pos, residue_name=res.get_resname().strip(),
                residue_number=resnum, one_letter=one_letter,
                ss_eight=ss8, ss_three=ss3, ss_fullname=ssname,
                phi=phi, psi=psi, method="DSSP",
            ))
        return records
    except Exception as exc:
        log.debug("DSSP failed: %s", exc)
        return None


# ===========================================================================
# Method 2: VMD/STRIDE
# ===========================================================================

def _vmd_available() -> str | None:
    """Return the VMD executable path if found, else None."""
    candidates = ["vmd", "vmd_LINUXAMD64", "vmd_WIN64.exe", "vmd_MACOSXX86_64"]
    for name in candidates:
        try:
            result = subprocess.run(
                [name, "-dispdev", "text", "-eofexit"],
                input="quit\n", capture_output=True, text=True, timeout=10,
            )
            # VMD exits 0 or 1 but will have written its banner
            if "VMD" in result.stdout or "VMD" in result.stderr:
                return name
        except (FileNotFoundError, subprocess.TimeoutExpired):
            continue
    return None


# VMD Tcl script template — runs STRIDE via mol ssrecalc and writes results
_VMD_TCL_TEMPLATE = """\
mol new {pdb_file} type pdb waitfor all
set molid [molinfo top]
mol ssrecalc $molid
set sel [atomselect $molid "protein and name CA"]
set resids  [$sel get resid]
set chains  [$sel get chain]
set resnames [$sel get resname]
set ss       [$sel get structure]
set phi_list [$sel get phi]
set psi_list [$sel get psi]
set fh [open "{out_file}" w]
foreach resid $resids chain $chains resname $resnames s $ss phi $phi_list psi $psi_list {{
    puts $fh "$chain\\t$resid\\t$resname\\t$s\\t$phi\\t$psi"
}}
close $fh
quit
"""

# VMD STRIDE single-letter to standardised 8-state code
# VMD reports: H=helix, G=3-10, I=pi, E=strand, B=bridge, T=turn, C=coil
# We map to the same alphabet used by DSSP for consistent downstream handling
_VMD_TO_SS8: dict[str, str] = {
    "H": "H", "G": "G", "I": "I",
    "E": "E", "B": "B", "b": "B",
    "T": "T", "S": "S", "C": "C",
    "": "C", "T": "T",
}


def _assign_vmd_stride(
    pdb_file: str | Path,
    window_residues: list[Residue],
    vmd_exe: str,
) -> list[ResidueSSRecord] | None:
    """Assign secondary structure using VMD's built-in STRIDE algorithm.

    VMD is run in headless batch mode (-dispdev text -eofexit) with a
    Tcl script that loads the PDB, triggers STRIDE via mol ssrecalc,
    reads per-residue structure codes, and writes them to a temp file.
    Python then reads that file and maps residue IDs to the motif window.

    Parameters
    ----------
    pdb_file : str or Path
        Full parent structure PDB file (STRIDE needs complete H-bond context).
    window_residues : list[Residue]
        Motif window residues to extract SS for.
    vmd_exe : str
        VMD executable name/path.

    Returns
    -------
    list[ResidueSSRecord] or None
        Per-residue records, or None if VMD fails.
    """
    pdb_file = str(Path(pdb_file).resolve())

    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            tcl_path = os.path.join(tmpdir, "stride.tcl")
            out_path = os.path.join(tmpdir, "stride_out.txt")

            # Write Tcl script
            tcl_script = _VMD_TCL_TEMPLATE.format(
                pdb_file=pdb_file,
                out_file=out_path,
            )
            with open(tcl_path, "w") as fh:
                fh.write(tcl_script)

            # Run VMD in batch mode
            result = subprocess.run(
                [vmd_exe, "-dispdev", "text", "-eofexit"],
                stdin=open(tcl_path),
                capture_output=True, text=True, timeout=60,
            )

            if not os.path.exists(out_path):
                log.debug(
                    "VMD/STRIDE: output file not created. VMD stderr: %s",
                    result.stderr[:500],
                )
                return None

            # Parse output: chain TAB resid TAB resname TAB ss TAB phi TAB psi
            stride_map: dict[tuple[str, int], tuple[str, float | None, float | None]] = {}
            with open(out_path) as fh:
                for line in fh:
                    parts = line.strip().split("\t")
                    if len(parts) < 6:
                        continue
                    chain, resid_str, resname, ss_raw, phi_str, psi_str = parts[:6]
                    try:
                        resid = int(resid_str)
                    except ValueError:
                        continue
                    ss8 = _VMD_TO_SS8.get(ss_raw, "C")
                    phi = float(phi_str) if phi_str not in ("", "{}") else None
                    psi = float(psi_str) if psi_str not in ("", "{}") else None
                    stride_map[(chain.strip(), resid)] = (ss8, phi, psi)

            # Build records for motif window residues
            records: list[ResidueSSRecord] = []
            for pos, res in enumerate(window_residues, 1):
                chain_id = res.get_parent().id
                resid    = res.id[1]
                resnum   = f"{resid}{res.id[2].strip()}"
                one_letter = _aa3to1(res.get_resname())

                key = (chain_id, resid)
                if key in stride_map:
                    ss8, phi, psi = stride_map[key]
                    phi = round(phi, 2) if phi is not None else None
                    psi = round(psi, 2) if psi is not None else None
                else:
                    ss8, phi, psi = "?", None, None

                ss3    = _SS_3STATE.get(ss8, "Coil")
                ssname = _SS_FULLNAME.get(ss8, "Coil")

                records.append(ResidueSSRecord(
                    position=pos,
                    residue_name=res.get_resname().strip(),
                    residue_number=resnum,
                    one_letter=one_letter,
                    ss_eight=ss8,
                    ss_three=ss3,
                    ss_fullname=ssname,
                    phi=phi,
                    psi=psi,
                    method="VMD_STRIDE",
                ))

            log.debug(
                "VMD/STRIDE assigned %d residues in motif window",
                len(records),
            )
            return records

    except Exception as exc:
        log.debug("VMD/STRIDE failed: %s", exc)
        return None


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
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                chain_residues = sorted(
                    [r for r in chain if is_aa(r)],
                    key=lambda r: (r.id[1], r.id[2]),
                )
        break
    phi_psi = _get_phi_psi(chain_residues)
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

    # --- Try DSSP first ---
    if _dssp_available():
        records = _assign_dssp(structure, pdb_file, window_residues)
        if records is not None:
            log.debug("%s: secondary structure assigned via DSSP", motif_id)
            return records
        log.warning("%s: DSSP binary found but failed — trying VMD/STRIDE", motif_id)

    # --- Try VMD/STRIDE ---
    vmd_exe = _vmd_available()
    if vmd_exe:
        records = _assign_vmd_stride(pdb_file, window_residues, vmd_exe)
        if records is not None:
            log.debug("%s: secondary structure assigned via VMD/STRIDE", motif_id)
            return records
        log.warning(
            "%s: VMD found but STRIDE failed — falling back to dihedral method",
            motif_id,
        )
    else:
        log.debug(
            "%s: VMD not found — using Bio.PDB dihedral fallback. "
            "Install VMD (https://www.ks.uiuc.edu/Research/vmd/) for "
            "accurate STRIDE-based assignment without needing mkdssp.",
            motif_id,
        )

    # --- Last resort: backbone dihedrals ---
    return _assign_biopdb(structure, window_residues)
