"""
motif.py
========
Core motif extraction logic for MetalMotifExtractor.

Algorithm Overview
------------------
For every metal ion found in a structure:

1.  Identify residues whose minimum atom-to-atom distance from the ion is
    within *distance_cutoff* Å.  Bio.PDB's :class:`NeighborSearch` is used
    so that only a KD-tree lookup is required instead of O(N·M) nested
    loops, which is critical for large proteins.

2.  Compute the *median* index of the binding residues within the chain's
    residue list.  Centering on the median (rather than the mean) is more
    robust when binding residues are non-contiguous.

3.  Extract a window of *motif_length* residues centered on that median
    index.  Windows that would extend beyond the chain termini are skipped
    (edge-of-chain artefacts).

4.  Save each motif as a PDB file including both the motif residues and the
    coordinating ion (``MotifSelect`` filter).

5.  Return a :class:`MotifRecord` dataclass for downstream CSV export.

Biological rationale
--------------------
The canonical EF-hand loop is 12 residues long with Ca²⁺ coordinated by
residues at positions 1, 3, 5, 7, 9 and 12 (X, Y, Z, -Y, -X, -Z in
crystallographic nomenclature).  A 15-residue window centred on the median
binding-residue index captures the full loop plus ≈1.5 residues of each
flanking helix, providing sequence context for motif logo construction and
machine-learning feature generation without including unrelated structural
elements.

EF-hands almost always occur in pairs; a single Ca²⁺ ion will therefore
often have binding residues drawn from two distinct regions of the sequence.
The median-centering strategy naturally positions the window between the two
coordination shells, capturing the core of the loop.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from Bio.PDB import PDBIO, PDBParser, Select, is_aa
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

from utils import AA3TO1, ensure_dirs
from coordination import analyse_coordination
from secondary_structure import assign_secondary_structure

log = logging.getLogger(__name__)


# ===========================================================================
# Data classes
# ===========================================================================

@dataclass
class MotifRecord:
    """Container for one extracted metal-binding motif.

    Attributes
    ----------
    motif_id : str
        Unique identifier composed of ``{pdb_id}_{chain}_{metal}_{index}``.
    pdb_id : str
        Four-character PDB ID.
    chain : str
        Chain ID in which the motif resides.
    metal : str
        Residue name of the metal ion (e.g. ``"CA"``).
    resolution : float or None
        Crystallographic resolution (Å).
    method : str
        Experimental method string.
    binding_residues : str
        Comma-separated three-letter names of coordinating residues.
        Only residues contributing a chemically valid donor atom within
        the cutoff are included; purely proximal residues without a donor
        atom are excluded.
    residue_numbers : str
        Comma-separated residue sequence numbers of coordinating residues.
    donor_atoms : str
        Comma-separated donor atom names for each coordinating residue,
        in the same order as ``binding_residues``.
        Examples: ``"OD1,OE1,OD1,O,OG1,OE1"`` for a typical EF-hand loop.
        Backbone oxygen donations are recorded as ``"O"``.
    coordination_type : str
        Comma-separated classification of each donor:
        ``"sidechain"`` or ``"backbone"``.
    distance_cutoff : float
        Distance cutoff used (Å).
    motif_length : int
        Number of residues in the extracted window.
    sequence : str
        One-letter amino-acid sequence of the motif window.
    output_file : str
        Path to the saved motif PDB file.
    """
    motif_id: str = ""
    pdb_id: str = ""
    chain: str = ""
    metal: str = ""
    resolution: float | None = None
    method: str = ""
    binding_residues: str = ""
    residue_numbers: str = ""
    donor_atoms: str = ""
    coordination_type: str = ""
    distance_cutoff: float = 4.0
    motif_length: int = 15       # actual extracted window length (may differ from
                                  # requested length at chain termini or when binding
                                  # residue span exceeds the requested length)
    sequence: str = ""
    output_file: str = ""
    # --- Coordination analysis (new) ---
    coordination_number: int = 0
    coordination_geometry: str = ""
    tau4: str = ""           # stored as string so "N/A" fits cleanly in CSV
    tau5: str = ""
    mean_lml_angle: float = 0.0
    min_lml_angle: float = 0.0
    max_lml_angle: float = 0.0
    angle_largest_1: float = 0.0
    angle_largest_2: float = 0.0
    # --- Secondary structure summary (new) ---
    ss_method: str = ""      # "DSSP" or "BioPDB_dihedral"
    ss_sequence: str = ""    # per-residue 8-state codes concatenated, e.g. "CCHHHHHHTCC"
    ss_helix_fraction: float = 0.0
    ss_strand_fraction: float = 0.0
    ss_coil_fraction: float = 0.0


# ===========================================================================
# PDBIO Select filter
# ===========================================================================

class MotifSelect(Select):
    """Bio.PDB ``Select`` subclass that retains only motif residues + the ion.

    Parameters
    ----------
    motif_residues : list[Residue]
        Amino-acid residues to include.
    ion : Residue
        Metal ion HETATM residue to include.
    """

    def __init__(self, motif_residues: list[Residue], ion: Residue) -> None:
        super().__init__()
        self._motif_set = set(id(r) for r in motif_residues)
        self._ion_id = id(ion)

    def accept_residue(self, residue: Residue) -> bool:
        return id(residue) in self._motif_set or id(residue) == self._ion_id


# ===========================================================================
# Metal ion identification
# ===========================================================================

def find_metal_ions(
    structure: Structure, metal: str | list[str]
) -> list[Residue]:
    """Locate all HETATM residues matching the requested metal name(s).

    Parameters
    ----------
    structure : Bio.PDB.Structure
        Parsed PDB structure.
    metal : str or list[str]
        Metal residue name(s) (case-insensitive), e.g. ``"CA"`` or
        ``["CA", "MG"]``. When a list is given, ions matching *any* of the
        listed residue names are returned.

    Returns
    -------
    list[Residue]
        All matching HETATM residues across all models and chains, sorted
        by residue name then sequence number for deterministic output.

    Notes
    -----
    PDB HETATM records have ``residue.id[0] != " "`` (the hetflag field is
    set to ``"H_<resname>"`` or ``"W"`` for waters).  We additionally check
    that the residue name matches one of *metal* to avoid confusing hetero
    amino acids or ligands with similar names.
    """
    metals = {metal.upper()} if isinstance(metal, str) else {m.upper() for m in metal}
    ions: list[Residue] = []

    for model in structure:
        for chain in model:
            for residue in chain:
                # Hetflag is non-blank for HETATM records
                if residue.id[0].strip() and residue.id[0] != "W":
                    if residue.get_resname().strip().upper() in metals:
                        ions.append(residue)

    ions.sort(key=lambda r: (r.get_resname().strip().upper(), r.id[1]))

    log.debug(
        "Found %d ion(s) matching %s in structure %s",
        len(ions), sorted(metals), structure.get_id(),
    )
    return ions


# ===========================================================================
# Donor atom definitions (coordination chemistry)
# ===========================================================================

# Side-chain donor atoms for each residue type.
# Derived from the coordination chemistry of biological metal-binding sites:
#
#   Oxygen donors  — dominate Ca²⁺, Mg²⁺, Mn²⁺ coordination
#   Nitrogen donors — common for Zn²⁺, Cu²⁺, Ni²⁺, Fe (His imidazole)
#   Sulfur donors  — Cys thiolate (Zn²⁺, Fe-S clusters, Cu²⁺)
#                    Met thioether (weaker; Cu²⁺, some Zn²⁺)
#
# References:
#   Harding (2001) Acta Cryst. D57:401 — Ca²⁺ coordination geometry
#   Andreini et al. (2008) J. Proteome Res. — metal donor census in PDB
#   Dudev & Lim (2014) Chem. Rev. — selectivity in metal-binding sites

_SIDECHAIN_DONORS: dict[str, frozenset[str]] = {
    "ASP": frozenset({"OD1", "OD2"}),   # carboxylate; bidentate possible
    "GLU": frozenset({"OE1", "OE2"}),   # carboxylate; bidentate possible
    "ASN": frozenset({"OD1"}),           # amide carbonyl oxygen
    "GLN": frozenset({"OE1"}),           # amide carbonyl oxygen
    "SER": frozenset({"OG"}),            # hydroxyl oxygen
    "THR": frozenset({"OG1"}),           # hydroxyl oxygen
    "TYR": frozenset({"OH"}),            # phenolate oxygen
    "HIS": frozenset({"ND1", "NE2"}),   # imidazole nitrogen (tautomer-dependent)
    "CYS": frozenset({"SG"}),            # thiolate sulfur
    "MET": frozenset({"SD"}),            # thioether sulfur (weaker donor)
}

# Backbone carbonyl oxygen is a legitimate donor in some EF-hand positions
# (notably position 9 of the canonical 12-residue loop, which is almost
# always a backbone O donation regardless of residue type).
# It is included but with a tighter effective cutoff (3.5 Å) applied
# internally because backbone O–metal distances are shorter on average
# than side-chain donor distances.
_BACKBONE_DONOR_ATOM = "O"
_BACKBONE_DONOR_CUTOFF_FRACTION = 0.875  # 3.5/4.0 — tighter than side-chain


# ===========================================================================
# Binding residue identification — donor-atom-aware
# ===========================================================================

def find_binding_residues(
    chain_residues: list[Residue],
    ion: Residue,
    distance_cutoff: float,
) -> list[tuple[Residue, str, str]]:
    """Identify residues coordinating the metal ion via valid donor atoms.

    Unlike a naive distance filter, this function checks that the atom
    closest to the ion on each candidate residue is a chemically valid
    donor atom (a specific side-chain O/N/S atom, or backbone carbonyl O).
    Residues that are merely geometrically proximal via their backbone Cα,
    N, or C atoms — or via non-donor side-chain atoms such as the Cβ of
    Phe or the ε-carbons of Lys — are excluded.

    Parameters
    ----------
    chain_residues : list[Residue]
        Standard amino-acid residues of the chain under examination.
    ion : Residue
        Metal ion HETATM residue.
    distance_cutoff : float
        Maximum donor-atom-to-ion distance (Å).  A tighter internal
        cutoff of ``distance_cutoff × 0.875`` is applied to backbone
        carbonyl O donations to match the shorter observed Ca²⁺–O(backbone)
        bond lengths (mean 2.4 Å vs 2.5 Å for side-chain oxygens).

    Returns
    -------
    list[tuple[Residue, str, str]]
        Each element is ``(residue, donor_atom_name, coordination_type)``
        where ``coordination_type`` is either ``"sidechain"`` or
        ``"backbone"``.  The list preserves chain order.

    Notes
    -----
    The distinction between sidechain and backbone donors is scientifically
    important: backbone donations are fixed by secondary structure and are
    insensitive to mutation at that position, whereas sidechain donations
    reflect the chemical identity of the residue and are the primary target
    of mutagenesis studies.
    """
    if not chain_residues:
        return []

    backbone_cutoff = distance_cutoff * _BACKBONE_DONOR_CUTOFF_FRACTION

    # Collect ion atom coordinates (Ca²⁺ has a single atom; some ions may
    # have multiple atoms, e.g. FE2S2 clusters, but for monatomic ions this
    # is always a single-atom residue)
    ion_coords = [
        atom.get_vector().get_array()
        for atom in ion.get_atoms()
    ]

    results: list[tuple[Residue, str, str]] = []

    for residue in chain_residues:
        resname = residue.get_resname().strip().upper()

        # Determine which atom names are valid donors for this residue type
        sidechain_donor_names = _SIDECHAIN_DONORS.get(resname, frozenset())

        best_donor_atom: str | None = None
        best_distance: float = 999.0
        best_coord_type: str = ""

        for atom in residue.get_atoms():
            atom_name = atom.get_name().strip()
            atom_coord = atom.get_vector().get_array()

            # Distance from this atom to nearest ion atom
            min_dist = min(
                float(np.linalg.norm(atom_coord - ion_coord))
                for ion_coord in ion_coords
            )

            # Check sidechain donors first (preferred over backbone)
            if atom_name in sidechain_donor_names:
                if min_dist <= distance_cutoff and min_dist < best_distance:
                    best_donor_atom = atom_name
                    best_distance = min_dist
                    best_coord_type = "sidechain"

            # Check backbone carbonyl O (applies to all residue types)
            elif atom_name == _BACKBONE_DONOR_ATOM:
                if min_dist <= backbone_cutoff and min_dist < best_distance:
                    best_donor_atom = atom_name
                    best_distance = min_dist
                    best_coord_type = "backbone"

        if best_donor_atom is not None:
            results.append((residue, best_donor_atom, best_coord_type))
            log.debug(
                "  Donor: %s %d %s via %s (%s, %.2f Å)",
                resname, residue.id[1], residue.id[2].strip(),
                best_donor_atom, best_coord_type, best_distance,
            )

    return results


# ===========================================================================
# Motif window extraction
# ===========================================================================

def extract_motif_window(
    chain_residues: list[Residue],
    binding_residues: list[Residue],
    motif_length: int,
) -> list[Residue] | None:
    """Extract a window that is guaranteed to contain all binding residues.

    Algorithm
    ---------
    The previous median-centred approach failed for asymmetrically distributed
    binding residues — specifically the canonical EF-hand pattern where
    positions 1, 3, 5, 7, 9 cluster at the N-terminal end of the loop and
    position 12 (the bidentate Glu) sits alone at the C-terminal end.  With
    binding indices {1, 3, 5, 12}, the median is 4, placing the window at
    indices [-2:10] when clamped, which excludes the critical 12th residue.

    The corrected algorithm:

    1.  Compute the span of binding residues: ``span_start = min(indices)``,
        ``span_end = max(indices)``.  This guarantees that every binding
        residue lies within [span_start, span_end].

    2.  Compute the padding available on each side:
        ``pad = max(0, motif_length - span_length) // 2``.

    3.  Expand symmetrically: ``start = span_start - pad``,
        ``end = span_end + remaining_pad + 1``.
        If motif_length is odd the extra residue goes to the C-terminal side,
        consistent with the previous convention.

    4.  Clamp to chain bounds: if ``start < 0`` shift the window right;
        if ``end > len(chain_residues)`` shift left.  After clamping, all
        binding residues are still guaranteed to be included (they are inside
        the span which is inside the window) — the only thing that changes
        is the amount of flanking context on each side.

    5.  If the binding residue span itself exceeds ``motif_length`` (unusual
        but possible for multi-site metal ions), the window is expanded to
        cover the full span without truncation, and a DEBUG log entry records
        the oversize window.

    Parameters
    ----------
    chain_residues : list[Residue]
        Ordered list of all standard amino-acid residues in the chain.
    binding_residues : list[Residue]
        Subset of *chain_residues* that coordinate the metal ion.
    motif_length : int
        Desired total number of residues in the output window.  The actual
        window may be larger if the binding residue span exceeds this value.

    Returns
    -------
    list[Residue] or None
        Motif residues in chain order, or ``None`` if the chain is too short
        to contain even the binding residues themselves.
    """
    if not binding_residues:
        return None

    id_to_index = {id(r): i for i, r in enumerate(chain_residues)}
    indices = [id_to_index[id(r)] for r in binding_residues if id(r) in id_to_index]

    if not indices:
        return None

    span_start = min(indices)
    span_end   = max(indices)
    span_len   = span_end - span_start + 1   # inclusive

    # If span already wider than requested motif_length, use span as-is
    if span_len >= motif_length:
        if span_len > motif_length:
            log.debug(
                "Binding residue span (%d) exceeds motif_length (%d) — "
                "window expanded to %d to include all binding residues",
                span_len, motif_length, span_len,
            )
        start = span_start
        end   = span_end + 1   # Python slice exclusive
    else:
        total_pad = motif_length - span_len
        pad_left  = total_pad // 2
        pad_right = total_pad - pad_left  # extra residue goes to right for odd pad
        start = span_start - pad_left
        end   = span_end   + pad_right + 1

    # Clamp to chain bounds — shift window without shrinking it
    if start < 0:
        end   -= start    # shift right by the same amount
        start  = 0
    if end > len(chain_residues):
        shift  = end - len(chain_residues)
        start  = max(0, start - shift)
        end    = len(chain_residues)

    # Sanity check: all binding residues must be inside the window
    actual_indices = set(range(start, end))
    missing = [i for i in indices if i not in actual_indices]
    if missing:
        # Should never happen after the logic above; log and continue anyway
        log.warning(
            "Window [%d:%d] still missing binding residue indices %s — "
            "chain may be too short",
            start, end, missing,
        )

    return chain_residues[start:end]


# ===========================================================================
# Sequence helper
# ===========================================================================

def motif_sequence(residues: list[Residue]) -> str:
    """Convert a list of residues to a one-letter amino-acid string.

    Parameters
    ----------
    residues : list[Residue]
        Ordered list of Bio.PDB residue objects.

    Returns
    -------
    str
        One-letter sequence (unknown residues represented as ``"X"``).
    """
    return "".join(AA3TO1.get(r.get_resname().strip().upper(), "X") for r in residues)


# ===========================================================================
# Save motif PDB
# ===========================================================================

def save_motif_pdb(
    structure: Structure,
    motif_residues: list[Residue],
    ion: Residue,
    output_path: str | Path,
) -> None:
    """Write a motif PDB file containing the window residues and the ion.

    Parameters
    ----------
    structure : Bio.PDB.Structure
        Parent structure (needed by PDBIO).
    motif_residues : list[Residue]
        Residues to include.
    ion : Residue
        Metal ion residue to include.
    output_path : str or Path
        Destination file path.
    """
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_path), MotifSelect(motif_residues, ion))
    log.debug("Saved motif → %s", output_path)


# ===========================================================================
# Main entry point — process one structure
# ===========================================================================

def extract_motifs_from_structure(
    pdb_id: str,
    pdb_file: str | Path,
    metal: str | list[str],
    distance_cutoff: float,
    motif_length: int,
    motif_dir: str | Path,
    structure_meta: dict[str, Any] | None = None,
) -> list[MotifRecord]:
    """Extract all metal-binding motifs from a single PDB file.

    Each metal ion site is processed independently, so a structure with
    three Zn²⁺ ions and two Ca²⁺ ions produces five separate motif files
    and records (assuming all five have valid, in-bounds windows).

    Parameters
    ----------
    pdb_id : str
        Four-character PDB ID (used in output filenames and records).
    pdb_file : str or Path
        Path to the PDB (or mmCIF) file.
    metal : str or list[str]
        Metal residue name(s) to search for, e.g. ``"ZN"`` or
        ``["CA", "MG"]``. Each ion is tagged in output filenames and
        records using its *own* residue name, not just the first entry
        in this list — so mixed-metal runs are never mislabeled.
    distance_cutoff : float
        Binding distance threshold (Å).
    motif_length : int
        Length of the extracted residue window.
    motif_dir : str or Path
        Output directory for motif PDB files.
    structure_meta : dict, optional
        Metadata dict from :func:`metadata.fetch_entry_metadata`.

    Returns
    -------
    list[MotifRecord]
        One record per successfully extracted motif.
    """
    pdb_id = pdb_id.upper()
    motif_dir = Path(motif_dir)
    ensure_dirs(motif_dir)

    metals_requested = [metal] if isinstance(metal, str) else list(metal)

    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure(pdb_id, str(pdb_file))
    except Exception as exc:
        log.error("Failed to parse %s (%s): %s", pdb_id, pdb_file, exc)
        return []

    ions = find_metal_ions(structure, metals_requested)
    if not ions:
        log.info(
            "%s: no ions found matching %s", pdb_id, ", ".join(metals_requested)
        )
        return []

    meta = structure_meta or {}
    resolution = meta.get("resolution")
    method = meta.get("method", "")

    records: list[MotifRecord] = []
    # Separate running motif counter per metal type, so e.g. CA and MG
    # motifs within the same structure each get their own clean numbering
    # (CA_motif1, CA_motif2, MG_motif1, ...) instead of sharing one counter.
    motif_index_by_metal: dict[str, int] = {}

    for model in structure:
        for chain in model:
            # Standard amino-acid residues in this chain, sorted by sequence number
            aa_residues: list[Residue] = sorted(
                [r for r in chain if is_aa(r)],
                key=lambda r: (r.id[1], r.id[2]),  # (seq_num, insertion_code)
            )
            if not aa_residues:
                continue

            for ion in ions:
                ion_metal = ion.get_resname().strip().upper()

                # find_binding_residues now returns list of
                # (residue, donor_atom_name, coordination_type)
                binding_tuples = find_binding_residues(
                    aa_residues, ion, distance_cutoff
                )
                if not binding_tuples:
                    continue

                # Unpack for downstream use
                binding_res   = [t[0] for t in binding_tuples]
                donor_atoms   = [t[1] for t in binding_tuples]
                coord_types   = [t[2] for t in binding_tuples]

                motif_res = extract_motif_window(aa_residues, binding_res, motif_length)
                if motif_res is None:
                    continue

                motif_index_by_metal[ion_metal] = motif_index_by_metal.get(ion_metal, 0) + 1
                idx = motif_index_by_metal[ion_metal]

                # Build output filename
                outfile = motif_dir / (
                    f"{pdb_id}_{chain.id}_{ion_metal}_motif{idx}.pdb"
                )
                save_motif_pdb(structure, motif_res, ion, outfile)

                # --- Coordination analysis -----------------------------------
                coord = analyse_coordination(binding_res, ion, distance_cutoff)

                # --- Secondary structure assignment -------------------------
                ss_records = assign_secondary_structure(
                    structure, pdb_file, motif_res,
                    motif_id=f"{pdb_id}_{chain.id}_{ion_metal}_{idx}",
                )
                ss_seq = "".join(r.ss_eight for r in ss_records)
                n_ss = len(ss_records) or 1
                helix_frac  = round(sum(1 for r in ss_records if r.ss_three == "Helix")  / n_ss, 3)
                strand_frac = round(sum(1 for r in ss_records if r.ss_three == "Strand") / n_ss, 3)
                coil_frac   = round(sum(1 for r in ss_records if r.ss_three == "Coil")   / n_ss, 3)
                ss_method_used = ss_records[0].method if ss_records else "none"

                record = MotifRecord(
                    motif_id=f"{pdb_id}_{chain.id}_{ion_metal}_{idx}",
                    pdb_id=pdb_id,
                    chain=chain.id,
                    metal=ion_metal,
                    resolution=resolution,
                    method=method,
                    binding_residues=",".join(
                        r.get_resname().strip() for r in binding_res
                    ),
                    residue_numbers=",".join(
                        f"{r.id[1]}{r.id[2].strip()}" for r in binding_res
                    ),
                    donor_atoms=",".join(donor_atoms),
                    coordination_type=",".join(coord_types),
                    distance_cutoff=distance_cutoff,
                    motif_length=len(motif_res),  # actual extracted length, not requested
                    sequence=motif_sequence(motif_res),
                    output_file=str(outfile),
                    # Coordination
                    coordination_number=coord.coordination_number,
                    coordination_geometry=coord.geometry,
                    tau4=str(coord.tau4) if coord.tau4 is not None else "N/A",
                    tau5=str(coord.tau5) if coord.tau5 is not None else "N/A",
                    mean_lml_angle=coord.mean_angle,
                    min_lml_angle=coord.min_angle,
                    max_lml_angle=coord.max_angle,
                    angle_largest_1=coord.angle_largest_1,
                    angle_largest_2=coord.angle_largest_2,
                    # Secondary structure summary
                    ss_method=ss_method_used,
                    ss_sequence=ss_seq,
                    ss_helix_fraction=helix_frac,
                    ss_strand_fraction=strand_frac,
                    ss_coil_fraction=coil_frac,
                )
                # Attach SS per-residue records for downstream CSV writing
                record._ss_records = ss_records  # type: ignore[attr-defined]

                records.append(record)
                log.info(
                    "%s chain %s: extracted %s motif %d — "
                    "%d donor(s): %s via %s",
                    pdb_id, chain.id, ion_metal, idx,
                    len(binding_res),
                    record.binding_residues,
                    record.donor_atoms,
                )

    return records
