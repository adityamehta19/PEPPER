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

5.  Coordinated water.  Every water oxygen within *water_cutoff* of the
    ion is reported (count, identifiers, distances) and its full residue
    (coordinates included) is written into the motif PDB file.  Waters are
    searched in *all* chains of the model, because deposited waters are
    frequently assigned to a different chain from the protein residues
    they contact.  The protein-donor search (steps 1-3) is unchanged.

6.  Coordinated non-protein ligands.  Every N, O or S atom of a residue that
    is neither an amino acid nor water (nucleotides, cofactors, substrates,
    anions, glycans, nucleic-acid phosphates, ...) lying within
    *ligand_cutoff* of the ion is reported, and the whole ligand residue is
    written into the motif file.  Single-atom metal ions are never counted
    as ligands.

Models
------
Each model of a multi-model file is processed on its own: an ion is only
ever matched to residues and waters of the *same* model.  By default
every model is used (``model_policy="all"``, which adds the model number
to motif IDs and file names for multi-model files); ``"first"`` uses only the
first model.

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
from scipy.spatial import cKDTree
from Bio.PDB import MMCIFIO, MMCIFParser, PDBIO, PDBParser, Select, is_aa
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

from utils import AA3TO1, ALL_METALS, COMMON_ADDITIVES, ensure_dirs, is_water_residue
from coordination import DonorAtom, analyse_coordination
from findgeo_tools import write_donor_pdb
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
        Comma-separated labels of the PROTEIN donor atoms of this chain
        (``"ASP8:OD1,ASP8:OD2,GLY10:O"``): every atom counts, so a bidentate
        carboxylate appears twice and ``binding_residues`` (one entry per
        residue) can be shorter than this list.
    coordination_type : str
        Comma-separated ``"sidechain"`` / ``"backbone"`` for each entry of
        ``donor_atoms`` (same order).
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
    coordination_number: int = 0     # TOTAL donor atoms: protein + water + non-protein ligand
    cn_protein: int = 0              # protein donor atoms only
    coordination_geometry: str = ""  # from ALL donor atoms
    geometry_protein_only: str = ""  # from the protein donor atoms alone
    geometry_rmsd: float | None = None          # RMSD (A) of the best template fit
    geometry_second: str = ""                   # runner-up template
    geometry_second_rmsd: float | None = None
    geometry_margin: float | None = None        # second_rmsd - rmsd (A)
    geometry_rmsd_protein_only: float | None = None
    all_donor_atoms: str = ""        # labels of ALL donor atoms (protein, water, ligand)
    donors_other_chains: int = 0     # protein donors of this ion in OTHER chains (not in this record)
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
    # --- Model and coordinated water (new) ---
    model_id: int = 0
    requested_length: int = 0          # requested window length (the actual one is motif_length)
    binding_span: int = 0              # residues from first to last protein ligand (list index)
    window_exceeds_requested: bool = False
    window_numbering_span: int = 0     # last - first residue number + 1 of the window
    window_missing_residues: int = 0   # residue numbers inside the window that are not modelled (gaps)
    other_chain_donor_atoms: str = ""  # protein donors of this ion in OTHER chains ("B:ASP5:OD1,...")
    ion_occupancy: float | None = None   # deposited occupancy / B-factor of the ion and the
    ion_bfactor: float | None = None     # shortest / longest metal-donor distances: provenance
    min_donor_distance: float | None = None   # for judging whether the metal assignment is plausible
    max_donor_distance: float | None = None
    ion_chain: str = ""                # chain, residue number and atom serial of the ION
    ion_resnum: str = ""               # (the other chain/ID fields describe the protein window)
    ion_serial: int = 0
    mean_donor_distance: float | None = None   # mean metal-donor distance (A)
    water_cutoff: float = 0.0          # metal...O(water) cutoff used (A); 0 = waters not searched
    water_count: int = 0               # waters with O within water_cutoff of the ion
    water_ids: str = ""                # "HOH:A:301,HOH:S:12" (resname:chain:resnum)
    water_distances: str = ""          # comma-separated Angstrom, same order as water_ids
    # --- Non-protein ligands (new) ---
    ligand_cutoff: float = 0.0         # metal...N/O/S cutoff used (A); 0 = ligands not searched
    ligand_atom_count: int = 0         # coordinated N/O/S atoms of non-protein, non-water residues
    ligand_residues: str = ""          # "ATP:A:501,GOL:A:502" (resname:chain:resnum)
    ligand_atoms: str = ""             # "ATP:A:501:O2G,..." nearest first
    ligand_distances: str = ""         # Angstrom, same order as ligand_atoms
    ligand_additives: str = ""         # ligand residues whose name is a common additive (flag only)


# ===========================================================================
# PDBIO Select filter
# ===========================================================================

class MotifSelect(Select):
    """Bio.PDB ``Select`` subclass keeping the motif residues, ion and waters.

    Only the model that contains the ion is written, so files from
    multi-model structures do not contain empty MODEL blocks.  Residues are
    matched by ``(chain id, residue id)`` rather than by Python object
    identity, because PDBIO hands ``accept_residue`` the *unpacked* children
    of disordered residues, which are different objects from the ones held
    in the motif list.

    Parameters
    ----------
    motif_residues : list[Residue]
        Amino-acid residues to include.
    ion : Residue
        Metal ion HETATM residue to include.
    water_residues : list[Residue], optional
        Coordinated water residues to include (all their atoms are kept,
        so their coordinates appear in the output file).
    model : Model, optional
        Model to write; defaults to the model that contains *ion*.
    ligand_residues : list[Residue], optional
        Coordinated non-protein ligand residues to include (all atoms).
    """

    def __init__(
        self,
        motif_residues: list[Residue],
        ion: Residue,
        water_residues: list[Residue] | None = None,
        model: Model | None = None,
        ligand_residues: list[Residue] | None = None,
    ) -> None:
        super().__init__()

        def key(res: Residue) -> tuple[str, tuple]:
            return (res.get_parent().id, res.id)

        self._keys = {key(r) for r in motif_residues}
        self._keys.add(key(ion))
        for w in water_residues or []:
            self._keys.add(key(w))
        for lig in ligand_residues or []:
            self._keys.add(key(lig))
        if model is None:
            model = ion.get_parent().get_parent()
        self._model_id = model.id
        self._chain_ids = {chain_id for chain_id, _ in self._keys}
        self._current_chain: str | None = None

    def accept_model(self, model) -> bool:
        return model.id == self._model_id

    def accept_chain(self, chain) -> bool:
        # Only chains that contain motif residues, the ion or coordinated
        # waters are written.  Returning True for every chain made PDBIO
        # raise "Chain id exceeds PDB format limit" for unrelated chains
        # with multi-character IDs (large mmCIF-only entries, e.g. "IA"),
        # even though none of their residues were being written.
        self._current_chain = chain.id
        return chain.id in self._chain_ids

    def accept_residue(self, residue: Residue) -> bool:
        return (self._current_chain, residue.id) in self._keys


# ===========================================================================
# Metal ion identification
# ===========================================================================

def find_metal_ions_in_model(
    model: Model, metal: str | list[str]
) -> list[Residue]:
    """Locate HETATM residues matching the requested metal name(s) in ONE model.

    Parameters
    ----------
    model : Bio.PDB.Model
        A single model.  Restricting the search to one model guarantees that
        an ion can only be matched to residues of its own model.
    metal : str or list[str]
        Metal residue name(s) (case-insensitive).

    Returns
    -------
    list[Residue]
        Matching HETATM residues sorted by residue name then sequence number.

    Notes
    -----
    PDB HETATM records have ``residue.id[0] != " "`` (the hetflag field is
    ``"H_<resname>"``, or ``"W"`` for waters).  The residue name must also
    match one of *metal*.
    """
    metals = {metal.upper()} if isinstance(metal, str) else {m.upper() for m in metal}
    ions: list[Residue] = []
    for chain in model:
        for residue in chain:
            if residue.id[0].strip() and residue.id[0] != "W":
                if residue.get_resname().strip().upper() in metals:
                    ions.append(residue)
    ions.sort(key=lambda r: (r.get_resname().strip().upper(), r.id[1]))
    return ions


def find_metal_ions(
    structure: Structure, metal: str | list[str]
) -> list[Residue]:
    """Locate matching HETATM residues across ALL models (compatibility helper).

    The extraction pipeline itself uses :func:`find_metal_ions_in_model`,
    one model at a time.  This function is kept for callers that only want
    to know which ions a structure contains.
    """
    ions: list[Residue] = []
    for model in structure:
        ions.extend(find_metal_ions_in_model(model, metal))
    ions.sort(key=lambda r: (r.get_resname().strip().upper(), r.id[1]))
    log.debug(
        "Found %d ion(s) matching %s in structure %s",
        len(ions), metal, structure.get_id(),
    )
    return ions


# ===========================================================================
# Coordinated water
# ===========================================================================

@dataclass
class WaterContact:
    """One water molecule whose oxygen lies within the water cutoff of an ion."""
    residue: Residue
    resname: str
    chain_id: str
    resnum: str            # residue number + insertion code
    distance: float        # metal ... O(water), Angstrom
    occupancy: float | None
    coord: np.ndarray      # oxygen coordinates

    @property
    def label(self) -> str:
        return f"{self.resname}:{self.chain_id}:{self.resnum}"

    @property
    def atom_label(self) -> str:
        return f"{self.label}:O"


@dataclass
class WaterIndex:
    """All water oxygens of one model, as arrays for fast distance queries."""
    residues: list[Residue]
    coords: np.ndarray            # shape (N, 3)


def _water_oxygen(residue: Residue):
    """Return the oxygen atom of a water residue (None if there is none)."""
    for atom in residue.get_atoms():
        element = (atom.element or "").strip().upper()
        name = atom.get_name().strip().upper()
        if element == "O" or (not element and name in ("O", "OH2", "OW")):
            return atom
    return None


def build_water_index(model: Model) -> WaterIndex:
    """Collect every water oxygen of a model (all chains).

    Deposited waters are often given their own chain ID, or the chain ID of
    a different protein chain, so waters must be gathered from the whole
    model and not from the chain being processed.
    """
    residues: list[Residue] = []
    coords: list[np.ndarray] = []
    for chain in model:
        for residue in chain:
            if not is_water_residue(residue):
                continue
            oxygen = _water_oxygen(residue)
            if oxygen is None:
                continue
            residues.append(residue)
            coords.append(oxygen.get_coord())
    arr = np.array(coords, dtype=float).reshape(-1, 3)
    return WaterIndex(residues=residues, coords=arr)


def find_coordinated_waters(
    water_index: WaterIndex | None,
    ion: Residue,
    water_cutoff: float,
) -> list[WaterContact]:
    """Return the waters whose oxygen is within *water_cutoff* of the ion.

    The distance is measured from the water oxygen to the metal (for
    multi-atom ions, to the nearest ion atom).  Results are sorted by
    distance, nearest first.
    """
    if water_index is None or len(water_index.residues) == 0:
        return []
    ion_coords = [a.get_coord() for a in ion.get_atoms()]
    if not ion_coords:
        return []

    dist = np.full(len(water_index.residues), np.inf)
    for ic in ion_coords:
        dist = np.minimum(dist, np.linalg.norm(water_index.coords - ic, axis=1))

    contacts: list[WaterContact] = []
    for idx in np.where(dist <= water_cutoff)[0]:
        res = water_index.residues[int(idx)]
        oxygen = _water_oxygen(res)
        occ = oxygen.get_occupancy() if oxygen is not None else None
        contacts.append(WaterContact(
            residue=res,
            resname=res.get_resname().strip(),
            chain_id=res.get_parent().id,
            resnum=f"{res.id[1]}{res.id[2].strip()}",
            distance=round(float(dist[idx]), 2),
            occupancy=float(occ) if occ is not None else None,
            coord=np.array(water_index.coords[int(idx)], dtype=float),
        ))
    contacts.sort(key=lambda w: w.distance)
    return contacts


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
# always a backbone O donation regardless of residue type).  It uses the
# SAME distance cutoff as the side-chain donors.
_BACKBONE_DONOR_ATOM = "O"


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
        Maximum donor-atom-to-ion distance (Å), applied equally to
        side-chain donors and backbone carbonyl O.

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

    backbone_cutoff = distance_cutoff   # same cutoff for side-chain and backbone donors

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


@dataclass
class ProteinDonor:
    """One protein donor atom bound to an ion."""
    residue: Residue
    atom_name: str
    coord_type: str        # "sidechain" or "backbone"
    distance: float
    coord: np.ndarray
    element: str = ""

    @property
    def label(self) -> str:
        r = self.residue
        return f"{r.get_resname().strip()}{r.id[1]}{r.id[2].strip()}:{self.atom_name}"


def find_protein_donors(
    chain_residues: list[Residue],
    ion: Residue,
    distance_cutoff: float,
) -> list[ProteinDonor]:
    """Return EVERY protein donor atom of one chain bound to the ion.

    This is the single definition of a protein donor atom used for the
    coordination number, the L-M-L angles, the geometry and the listed
    ``Donor Atoms``.

    * Side-chain donors (table ``_SIDECHAIN_DONORS``): within
      *distance_cutoff* of the ion.
    * Backbone carbonyl O (any residue): within the same *distance_cutoff*.

    Unlike :func:`find_binding_residues`, which keeps only the single
    closest donor atom per residue, **all** qualifying atoms are returned, so
    a bidentate Asp/Glu contributes two donors, and a residue donating both
    from its side chain and its backbone contributes both.  Donors are
    returned in chain order.
    """
    if not chain_residues:
        return []
    backbone_cutoff = distance_cutoff   # same cutoff for side-chain and backbone donors
    ion_coords = [a.get_coord() for a in ion.get_atoms()]
    if not ion_coords:
        return []

    donors: list[ProteinDonor] = []
    for residue in chain_residues:
        resname = residue.get_resname().strip().upper()
        sidechain_names = _SIDECHAIN_DONORS.get(resname, frozenset())
        for atom in residue.get_atoms():
            name = atom.get_name().strip()
            if name in sidechain_names:
                cutoff, ctype = distance_cutoff, "sidechain"
            elif name == _BACKBONE_DONOR_ATOM:
                cutoff, ctype = backbone_cutoff, "backbone"
            else:
                continue
            coord = atom.get_coord()
            dist = min(float(np.linalg.norm(coord - ic)) for ic in ion_coords)
            if dist <= cutoff:
                donors.append(ProteinDonor(
                    residue=residue, atom_name=name, coord_type=ctype,
                    distance=round(dist, 2), coord=np.array(coord, dtype=float),
                    element=_atom_element(atom),
                ))
    return donors


class ProteinDonorIndex:
    """Spatial index of every candidate protein donor atom of one model.

    ``find_protein_donors`` scans all atoms of a chain for every ion, which
    costs (chains x ions x atoms).  For a large complex with hundreds or
    thousands of ions (the slow cases reported by a reviewer) this reaches
    minutes to hours.  This class collects the candidate donor atoms
    (side-chain donor atoms from ``_SIDECHAIN_DONORS`` and backbone O) of all
    chains **once** into a k-d tree, so each ion needs one radius query.

    The result is identical to :func:`find_protein_donors` (same atoms, same
    ``<=`` cutoff, same chain / residue / atom order), only faster.
    """

    def __init__(self, aa_by_chain: dict[str, list[Residue]]) -> None:
        coords: list[np.ndarray] = []
        self._meta: list[tuple] = []     # (chain_id, residue index, atom order, residue, name, ctype, element)
        for cid, residues in aa_by_chain.items():
            for ri, residue in enumerate(residues):
                resname = residue.get_resname().strip().upper()
                sidechain_names = _SIDECHAIN_DONORS.get(resname, frozenset())
                for ai, atom in enumerate(residue.get_atoms()):
                    name = atom.get_name().strip()
                    if name in sidechain_names:
                        ctype = "sidechain"
                    elif name == _BACKBONE_DONOR_ATOM:
                        ctype = "backbone"
                    else:
                        continue
                    coords.append(atom.get_coord())
                    self._meta.append((cid, ri, ai, residue, name, ctype, _atom_element(atom)))
        self._coords = np.array(coords, dtype=float).reshape(-1, 3)
        self._tree = cKDTree(self._coords) if len(coords) else None

    def donors_for_ion(self, ion: Residue, cutoff: float) -> dict[str, list[ProteinDonor]]:
        """Protein donor atoms of the ion, grouped by chain (chain order kept)."""
        if self._tree is None:
            return {}
        ion_coords = [a.get_coord() for a in ion.get_atoms()]
        best: dict[int, float] = {}
        for ic in ion_coords:
            for i in self._tree.query_ball_point(ic, r=cutoff + 1e-6):
                d = float(np.linalg.norm(self._coords[i] - ic))
                if d <= cutoff and d < best.get(i, np.inf):
                    best[i] = d
        grouped: dict[str, list[tuple]] = {}
        for i, d in best.items():
            cid, ri, ai, residue, name, ctype, element = self._meta[i]
            grouped.setdefault(cid, []).append((ri, ai, ProteinDonor(
                residue=residue, atom_name=name, coord_type=ctype,
                distance=round(d, 2), coord=np.array(self._coords[i], dtype=float),
                element=element,
            )))
        return {cid: [t[2] for t in sorted(items, key=lambda t: (t[0], t[1]))]
                for cid, items in grouped.items()}


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


def _ion_occupancy(ion: Residue) -> float | None:
    for atom in ion.get_atoms():
        occ = atom.get_occupancy()
        return None if occ is None else round(float(occ), 2)
    return None


def _ion_bfactor(ion: Residue) -> float | None:
    for atom in ion.get_atoms():
        b = atom.get_bfactor()
        return None if b is None else round(float(b), 2)
    return None


def _ion_serial(ion: Residue) -> int:
    """PDB atom serial number of the (first atom of the) ion; 0 if unavailable."""
    for atom in ion.get_atoms():
        return int(atom.get_serial_number() or 0)
    return 0


# ===========================================================================
# Coordinated non-protein ligands
# ===========================================================================

_DONOR_ELEMENTS = frozenset({"N", "O", "S"})
_METAL_NAMES = frozenset(m.upper() for m in ALL_METALS)


def _atom_element(atom) -> str:
    """Upper-case element symbol of an atom (guessed from the name if blank)."""
    el = (atom.element or "").strip().upper()
    if el:
        return el
    letters = "".join(ch for ch in atom.get_name().strip() if ch.isalpha())
    return letters[:1].upper()


@dataclass
class LigandContact:
    """One donor atom (N, O or S) of a non-protein ligand near an ion."""
    residue: Residue
    resname: str
    chain_id: str
    resnum: str            # residue number + insertion code
    atom_name: str
    element: str
    distance: float        # metal ... donor atom, Angstrom
    coord: np.ndarray      # donor atom coordinates

    @property
    def residue_label(self) -> str:
        return f"{self.resname}:{self.chain_id}:{self.resnum}"

    @property
    def atom_label(self) -> str:
        return f"{self.residue_label}:{self.atom_name}"


@dataclass
class LigandIndex:
    """All candidate donor atoms of non-protein, non-water residues of a model."""
    residues: list[Residue]       # parent residue of each atom
    atom_names: list[str]
    elements: list[str]
    coords: np.ndarray            # shape (N, 3)


def build_ligand_index(model: Model) -> LigandIndex:
    """Collect N/O/S atoms of every non-protein, non-water residue of a model.

    "Non-protein" means: not an amino acid (``is_aa``) and not water.  This
    includes HETATM ligands (ATP, heme, sulfate, glycerol, drug-like
    molecules, ...) and nucleic-acid residues (so RNA/DNA phosphate oxygens
    binding Mg2+ are found).  Single-atom metal ions are skipped, so a
    neighbouring ion is never mistaken for a ligand.  All chains of the
    model are searched.
    """
    residues: list[Residue] = []
    names: list[str] = []
    elements: list[str] = []
    coords: list[np.ndarray] = []
    for chain in model:
        for residue in chain:
            if is_water_residue(residue) or is_aa(residue):
                continue
            atoms = list(residue.get_atoms())
            if len(atoms) == 1 and residue.get_resname().strip().upper() in _METAL_NAMES:
                continue
            for atom in atoms:
                el = _atom_element(atom)
                if el in _DONOR_ELEMENTS:
                    residues.append(residue)
                    names.append(atom.get_name().strip())
                    elements.append(el)
                    coords.append(atom.get_coord())
    arr = np.array(coords, dtype=float).reshape(-1, 3)
    return LigandIndex(residues=residues, atom_names=names, elements=elements, coords=arr)


def find_coordinated_ligands(
    ligand_index: LigandIndex | None,
    ion: Residue,
    ligand_cutoff: float,
) -> list[LigandContact]:
    """Return non-protein ligand donor atoms within *ligand_cutoff* of the ion.

    Distance is measured to the nearest ion atom.  Results are sorted by
    distance, nearest first.  Each atom is one contact, so a bidentate
    phosphate contributes two.
    """
    if ligand_index is None or len(ligand_index.residues) == 0:
        return []
    ion_coords = [a.get_coord() for a in ion.get_atoms()]
    if not ion_coords:
        return []

    dist = np.full(len(ligand_index.residues), np.inf)
    for ic in ion_coords:
        dist = np.minimum(dist, np.linalg.norm(ligand_index.coords - ic, axis=1))

    contacts: list[LigandContact] = []
    for idx in np.where(dist <= ligand_cutoff)[0]:
        i = int(idx)
        res = ligand_index.residues[i]
        contacts.append(LigandContact(
            residue=res,
            resname=res.get_resname().strip(),
            chain_id=res.get_parent().id,
            resnum=f"{res.id[1]}{res.id[2].strip()}",
            atom_name=ligand_index.atom_names[i],
            element=ligand_index.elements[i],
            distance=round(float(dist[i]), 2),
            coord=np.array(ligand_index.coords[i], dtype=float),
        ))
    contacts.sort(key=lambda c: c.distance)
    return contacts


def unique_ligand_residues(contacts: list[LigandContact]) -> list[Residue]:
    """Distinct ligand residues of a contact list, nearest contact first."""
    seen: set[tuple] = set()
    out: list[Residue] = []
    for c in contacts:
        key = (c.chain_id, c.residue.id)
        if key not in seen:
            seen.add(key)
            out.append(c.residue)
    return out


# ===========================================================================
# Structure input / motif output
# ===========================================================================

def _open_structure(pdb_id: str, pdb_file: str | Path) -> Structure:
    """Parse a structure file with the parser that matches its format.

    ``.cif`` / ``.mmcif`` files (biological assemblies, and asymmetric units
    that have no PDB-format file) go through ``MMCIFParser``; everything
    else through ``PDBParser``.
    """
    suffix = Path(pdb_file).suffix.lower()
    if suffix in (".cif", ".mmcif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure(pdb_id, str(pdb_file))


def save_motif_pdb(
    structure: Structure,
    motif_residues: list[Residue],
    ion: Residue,
    output_path: str | Path,
    water_residues: list[Residue] | None = None,
    ligand_residues: list[Residue] | None = None,
) -> Path:
    """Write the motif file: window residues, the ion, coordinated waters and
    coordinated non-protein ligands.

    Parameters
    ----------
    structure : Bio.PDB.Structure
        Parent structure (needed by PDBIO / MMCIFIO).
    motif_residues : list[Residue]
        Residues to include.
    ion : Residue
        Metal ion residue to include.
    output_path : str or Path
        Destination file path (``.pdb``).
    water_residues : list[Residue], optional
        Coordinated water residues to include, with their coordinates.
    ligand_residues : list[Residue], optional
        Coordinated non-protein ligand residues to include (all atoms).

    Returns
    -------
    Path
        The file actually written.  This is *output_path* unless a chain ID
        is longer than one character (as happens with biological assemblies,
        e.g. ``A-2``), which the PDB format cannot hold; the motif is then
        written as mmCIF with a ``.cif`` extension.
    """
    output_path = Path(output_path)
    waters = list(water_residues or [])
    ligands = list(ligand_residues or [])
    ion_model = ion.get_parent().get_parent()

    # Write a small structure that holds copies of only the needed residues.
    # (PDBIO with a Select filter walks EVERY atom of the parent structure
    # for every motif, which is slow for large complexes.)
    needed = list(motif_residues) + waters + ligands + [ion]
    by_chain: dict[str, list[Residue]] = {}
    for r in needed:
        by_chain.setdefault(r.get_parent().id, []).append(r)
    small = Structure(structure.id)
    small_model = Model(ion_model.id)
    small.add(small_model)
    seen_keys: set[tuple] = set()
    for chain in ion_model:                       # keep the original chain order
        if chain.id not in by_chain:
            continue
        new_chain = Chain(chain.id)
        group = sorted(
            by_chain[chain.id],
            key=lambda r: (r.id[0] != " ", r.id[1], r.id[2]),   # ATOM residues first, then hetero
        )
        for r in group:
            key = (chain.id, r.id)
            if key in seen_keys:
                continue
            seen_keys.add(key)
            new_chain.add(r.copy())
        small_model.add(new_chain)

    chain_ids = set(by_chain)
    if any(len(str(c)) > 1 for c in chain_ids):
        output_path = output_path.with_suffix(".cif")
        io = MMCIFIO()
    else:
        io = PDBIO()
    io.set_structure(small)
    io.save(str(output_path))
    log.debug("Saved motif -> %s", output_path)
    return output_path


# ===========================================================================
# Main entry point - process one structure
# ===========================================================================

STATUS_OK = "sites_extracted"
STATUS_NO_ION = "no_ion_of_requested_type"
STATUS_NO_LIGAND = "ion_without_protein_ligand"
STATUS_CLUSTER_ONLY = "metal_cluster_not_supported"
STATUS_PARSE_FAILURE = "parse_failure"


def extract_motifs_with_status(
    pdb_id: str,
    pdb_file: str | Path,
    metal: str | list[str],
    distance_cutoff: float,
    motif_length: int,
    motif_dir: str | Path,
    structure_meta: dict[str, Any] | None = None,
    include_water: bool = True,
    water_cutoff: float | None = None,
    model_policy: str = "all",
    include_ligands: bool = True,
    ligand_cutoff: float | None = None,
    geometry_max_rmsd: float | None = None,
    findgeo_input_dir: str | Path | None = None,
    metal_cutoffs: dict[str, float] | None = None,
) -> tuple[list[MotifRecord], str]:
    """Extract all metal-binding motifs from a single structure file.

    Each metal ion is processed independently, and each **model** is
    processed independently: an ion is matched only against the protein
    residues and waters of its own model.

    Parameters
    ----------
    pdb_id : str
        Four-character PDB ID.
    pdb_file : str or Path
        Path to the ``.pdb`` or ``.cif`` file.
    metal : str or list[str]
        Metal residue name(s) to search for.
    distance_cutoff : float
        Protein-donor distance threshold (A).
    motif_length : int
        Length of the extracted residue window.
    motif_dir : str or Path
        Output directory for motif files.
    structure_meta : dict, optional
        Metadata dict from :func:`metadata.fetch_entry_metadata`.
    include_water : bool
        Report coordinated waters and write them into the motif file.
    water_cutoff : float, optional
        Maximum metal...O(water) distance (A) for a water to be counted as
        coordinated.  Default: the same as *distance_cutoff*.
    model_policy : {"all", "first"}
        Use every model (default), or only the first model.
    include_ligands : bool
        Report coordinated non-protein ligands and write them into the
        motif file.
    ligand_cutoff : float, optional
        Maximum metal...N/O/S distance (A) for a non-protein ligand atom to
        be counted as coordinated.  Default: the same as *distance_cutoff*.
    geometry_max_rmsd : float, optional
        Label sites whose best template fit has an RMSD (A) above this value
        "Irregular".  None (default): no threshold.
    metal_cutoffs : dict, optional
        Per-metal protein-donor cutoffs (A), e.g. ``{"ZN": 2.8, "CA": 3.2}``,
        that override *distance_cutoff* for ions of that metal.  The cutoff
        actually used is stored in each record (``distance_cutoff``).
    findgeo_input_dir : str or Path, optional
        If given, write for every motif a small PDB file ``<motif_id>.pdb``
        that contains only the metal and exactly the donor atoms counted in
        the coordination number.  Running FindGeo on these files compares the
        two programs' geometry assignment on identical atoms
        (see ``findgeo_compare.py``).

    Returns
    -------
    (list[MotifRecord], str)
        The records and a status string describing the structure outcome:
        ``sites_extracted``, ``no_ion_of_requested_type``,
        ``ion_without_protein_ligand`` or ``parse_failure``.
    """
    pdb_id = pdb_id.upper()
    motif_dir = Path(motif_dir)
    ensure_dirs(motif_dir)

    metals_requested = [metal] if isinstance(metal, str) else list(metal)

    try:
        structure = _open_structure(pdb_id, pdb_file)
    except Exception as exc:
        log.error("Failed to parse %s (%s): %s", pdb_id, pdb_file, exc)
        return [], STATUS_PARSE_FAILURE

    models = structure.get_list()
    if model_policy == "first":
        models = models[:1]
    multi_model_tag = model_policy == "all" and len(models) > 1

    cutoff_by_metal = {k.upper(): float(v) for k, v in (metal_cutoffs or {}).items()}
    # water / ligand cut-offs have no default of their own: same as the protein cut-off
    if water_cutoff is None:
        water_cutoff = distance_cutoff
    if ligand_cutoff is None:
        ligand_cutoff = distance_cutoff
    meta = structure_meta or {}
    resolution = meta.get("resolution")
    method = meta.get("method", "")

    records: list[MotifRecord] = []
    any_ion = False

    for model in models:
        ions = find_metal_ions_in_model(model, metals_requested)
        if not ions:
            continue
        any_ion = True

        water_index = build_water_index(model) if include_water else None
        waters_by_ion: dict[int, list[WaterContact]] = {}
        ligand_index = build_ligand_index(model) if include_ligands else None
        ligands_by_ion: dict[int, list[LigandContact]] = {}

        # Separate running motif counter per metal type (per model)
        motif_index_by_metal: dict[str, int] = {}
        model_tag = f"m{model.id}_" if multi_model_tag else ""

        # Standard amino-acid residues of every chain, and the protein donor
        # atoms of every (chain, ion) pair.  Computed once so that the number
        # of donors in the OTHER chains of the same ion can be reported.
        aa_by_chain: dict[str, list[Residue]] = {}
        for chain in model:
            residues_sorted = sorted(
                [r for r in chain if is_aa(r)],
                key=lambda r: (r.id[1], r.id[2]),  # (seq_num, insertion_code)
            )
            if residues_sorted:
                aa_by_chain[chain.id] = residues_sorted
        donor_index = ProteinDonorIndex(aa_by_chain)
        index_in_chain = {
            cid: {id(r): i for i, r in enumerate(res_list)}
            for cid, res_list in aa_by_chain.items()
        }
        donors_by_ion: dict[int, dict[str, list[ProteinDonor]]] = {}
        n_donors_by_ion: dict[int, int] = {}
        applied_cutoff: dict[int, float] = {}
        for ion in ions:
            cut = cutoff_by_metal.get(ion.get_resname().strip().upper(), distance_cutoff)
            applied_cutoff[id(ion)] = cut
            donors_by_ion[id(ion)] = donor_index.donors_for_ion(ion, cut)
            n_donors_by_ion[id(ion)] = sum(len(v) for v in donors_by_ion[id(ion)].values())

        for chain in model:
            aa_residues = aa_by_chain.get(chain.id)
            if not aa_residues:
                continue

            for ion in ions:
                ion_metal = ion.get_resname().strip().upper()

                protein_donors = donors_by_ion[id(ion)].get(chain.id, [])
                if not protein_donors:
                    continue
                n_other_chain_donors = n_donors_by_ion[id(ion)] - len(protein_donors)
                ion_cutoff = applied_cutoff[id(ion)]

                # distinct binding residues in chain order (basis of the window)
                binding_res: list[Residue] = []
                for d in protein_donors:
                    if not binding_res or binding_res[-1] is not d.residue:
                        binding_res.append(d.residue)
                donor_labels = [d.label for d in protein_donors]
                coord_types = [d.coord_type for d in protein_donors]

                motif_res = extract_motif_window(aa_residues, binding_res, motif_length)
                if motif_res is None:
                    continue
                idxs = [index_in_chain[chain.id][id(r)] for r in binding_res]
                binding_span = max(idxs) - min(idxs) + 1
                # numbering span of the window and residue numbers that are not
                # modelled inside it (gaps; approximate when insertion codes occur)
                numbering_span = motif_res[-1].id[1] - motif_res[0].id[1] + 1
                missing_in_window = max(0, numbering_span - len(motif_res))
                other_chain_atoms = ",".join(
                    f"{cid}:{d.label}"
                    for cid in aa_by_chain if cid != chain.id
                    for d in donors_by_ion[id(ion)].get(cid, [])
                )

                # Coordinated waters (same for every chain record of this ion)
                if include_water:
                    if id(ion) not in waters_by_ion:
                        waters_by_ion[id(ion)] = find_coordinated_waters(
                            water_index, ion, water_cutoff
                        )
                    waters = waters_by_ion[id(ion)]
                else:
                    waters = []

                # Coordinated non-protein ligands (same for every chain record of this ion)
                if include_ligands:
                    if id(ion) not in ligands_by_ion:
                        ligands_by_ion[id(ion)] = find_coordinated_ligands(
                            ligand_index, ion, ligand_cutoff
                        )
                    lig_contacts = ligands_by_ion[id(ion)]
                else:
                    lig_contacts = []
                lig_residues = unique_ligand_residues(lig_contacts)

                motif_index_by_metal[ion_metal] = motif_index_by_metal.get(ion_metal, 0) + 1
                idx = motif_index_by_metal[ion_metal]

                motif_id = f"{pdb_id}_{model_tag}{chain.id}_{ion_metal}_{idx}"
                outfile = motif_dir / f"{pdb_id}_{model_tag}{chain.id}_{ion_metal}_motif{idx}.pdb"
                outfile = save_motif_pdb(
                    structure, motif_res, ion, outfile,
                    water_residues=[w.residue for w in waters],
                    ligand_residues=lig_residues,
                )

                # --- Coordination analysis ----------------------------------
                # ONE donor list: protein donors + coordinated waters (if
                # requested) + non-protein ligand atoms (if requested).  Its
                # length is the coordination number.
                donor_list: list[DonorAtom] = (
                    [DonorAtom(
                        d.label, d.coord, "protein",
                        resname=d.residue.get_resname().strip(), chain=chain.id,
                        resnum=f"{d.residue.id[1]}{d.residue.id[2].strip()}",
                        atom_name=d.atom_name, element=d.element, hetero=False,
                    ) for d in protein_donors]
                    + [DonorAtom(
                        w.atom_label, w.coord, "water",
                        resname=w.resname, chain=w.chain_id, resnum=w.resnum,
                        atom_name="O", element="O", hetero=True,
                    ) for w in waters]
                    + [DonorAtom(
                        c.atom_label, c.coord, "ligand",
                        resname=c.resname, chain=c.chain_id, resnum=c.resnum,
                        atom_name=c.atom_name, element=c.element,
                        hetero=(c.residue.id[0] != " "),
                    ) for c in lig_contacts]
                )
                coord = analyse_coordination(donor_list, ion, geometry_max_rmsd)
                if findgeo_input_dir is not None:
                    write_donor_pdb(
                        Path(findgeo_input_dir) / f"{motif_id}.pdb", ion, donor_list
                    )

                # --- Secondary structure assignment -------------------------
                ss_records = assign_secondary_structure(
                    structure, pdb_file, motif_res, motif_id=motif_id,
                )
                ss_seq = "".join(r.ss_eight for r in ss_records)
                n_ss = len(ss_records) or 1
                helix_frac  = round(sum(1 for r in ss_records if r.ss_three == "Helix")  / n_ss, 3)
                strand_frac = round(sum(1 for r in ss_records if r.ss_three == "Strand") / n_ss, 3)
                coil_frac   = round(sum(1 for r in ss_records if r.ss_three == "Coil")   / n_ss, 3)
                ss_method_used = ss_records[0].method if ss_records else "none"

                record = MotifRecord(
                    motif_id=motif_id,
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
                    donor_atoms=",".join(donor_labels),
                    coordination_type=",".join(coord_types),
                    distance_cutoff=ion_cutoff,
                    motif_length=len(motif_res),
                    requested_length=motif_length,
                    binding_span=binding_span,
                    window_exceeds_requested=len(motif_res) > motif_length,
                    window_numbering_span=numbering_span,
                    window_missing_residues=missing_in_window,
                    other_chain_donor_atoms=other_chain_atoms,
                    ion_occupancy=_ion_occupancy(ion),
                    ion_bfactor=_ion_bfactor(ion),
                    min_donor_distance=coord.min_distance,
                    max_donor_distance=coord.max_distance,
                    sequence=motif_sequence(motif_res),
                    output_file=str(outfile),
                    coordination_number=coord.coordination_number,
                    cn_protein=coord.cn_protein,
                    coordination_geometry=coord.geometry,
                    geometry_protein_only=coord.geometry_protein_only,
                    geometry_rmsd=coord.geometry_rmsd,
                    geometry_second=coord.geometry_second,
                    geometry_second_rmsd=coord.geometry_second_rmsd,
                    geometry_margin=coord.geometry_margin,
                    geometry_rmsd_protein_only=coord.geometry_rmsd_protein_only,
                    all_donor_atoms=coord.donor_atoms,
                    donors_other_chains=n_other_chain_donors,
                    mean_lml_angle=coord.mean_angle,
                    min_lml_angle=coord.min_angle,
                    max_lml_angle=coord.max_angle,
                    angle_largest_1=coord.angle_largest_1,
                    angle_largest_2=coord.angle_largest_2,
                    ss_method=ss_method_used,
                    ss_sequence=ss_seq,
                    ss_helix_fraction=helix_frac,
                    ss_strand_fraction=strand_frac,
                    ss_coil_fraction=coil_frac,
                    model_id=model.id,
                    ion_chain=ion.get_parent().id,
                    ion_resnum=f"{ion.id[1]}{ion.id[2].strip()}",
                    ion_serial=_ion_serial(ion),
                    mean_donor_distance=coord.mean_distance,
                    water_cutoff=water_cutoff if include_water else 0.0,
                    water_count=len(waters),
                    water_ids=",".join(w.label for w in waters),
                    water_distances=",".join(f"{w.distance:.2f}" for w in waters),
                    ligand_cutoff=ligand_cutoff if include_ligands else 0.0,
                    ligand_atom_count=len(lig_contacts),
                    ligand_residues=",".join(
                        dict.fromkeys(c.residue_label for c in lig_contacts)
                    ),
                    ligand_atoms=",".join(c.atom_label for c in lig_contacts),
                    ligand_distances=",".join(f"{c.distance:.2f}" for c in lig_contacts),
                    ligand_additives=",".join(
                        dict.fromkeys(
                            c.residue_label for c in lig_contacts
                            if c.resname.upper() in COMMON_ADDITIVES
                        )
                    ),
                )
                record._ss_records = ss_records  # type: ignore[attr-defined]

                records.append(record)
                log.info(
                    "%s chain %s: extracted %s motif %d - %d donor(s): %s via %s; "
                    "%d coordinated water(s); %d non-protein ligand atom(s)",
                    pdb_id, chain.id, ion_metal, idx,
                    len(binding_res), record.binding_residues,
                    record.donor_atoms, len(waters), len(lig_contacts),
                )

    if records:
        return records, STATUS_OK
    if any_ion:
        return records, STATUS_NO_LIGAND
    from metadata import UNSUPPORTED_CLUSTER_IDS
    cluster_ids = {
        r.get_resname().strip().upper()
        for model in structure for r in model.get_residues()
        if r.get_resname().strip().upper() in UNSUPPORTED_CLUSTER_IDS
    }
    if cluster_ids:
        log.info("%s: no single-atom ion matching %s, but the structure contains a metal "
                 "cluster %s, which is not supported", pdb_id, ", ".join(metals_requested),
                 sorted(cluster_ids))
        return records, STATUS_CLUSTER_ONLY
    log.info("%s: no ions found matching %s", pdb_id, ", ".join(metals_requested))
    return records, STATUS_NO_ION


def extract_motifs_from_structure(
    pdb_id: str,
    pdb_file: str | Path,
    metal: str | list[str],
    distance_cutoff: float,
    motif_length: int,
    motif_dir: str | Path,
    structure_meta: dict[str, Any] | None = None,
    include_water: bool = True,
    water_cutoff: float | None = None,
    model_policy: str = "all",
    include_ligands: bool = True,
    ligand_cutoff: float | None = None,
    geometry_max_rmsd: float | None = None,
    findgeo_input_dir: str | Path | None = None,
    metal_cutoffs: dict[str, float] | None = None,
) -> list[MotifRecord]:
    """Records-only wrapper around :func:`extract_motifs_with_status`."""
    records, _ = extract_motifs_with_status(
        pdb_id, pdb_file, metal, distance_cutoff, motif_length, motif_dir,
        structure_meta=structure_meta,
        include_water=include_water, water_cutoff=water_cutoff,
        model_policy=model_policy,
        include_ligands=include_ligands, ligand_cutoff=ligand_cutoff,
        geometry_max_rmsd=geometry_max_rmsd, findgeo_input_dir=findgeo_input_dir,
        metal_cutoffs=metal_cutoffs,
    )
    return records
