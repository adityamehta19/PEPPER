"""
coordination.py
===============
Computes coordination number, coordination geometry classification, and
key inter-ligand angles for a metal ion in a protein structure.

Coordination Number
-------------------
Defined here as the number of individual donor atoms — not donor residues —
within the distance cutoff of the metal ion.  For example, a bidentate
glutamate contributes two oxygen atoms (Oε1 and Oε2) if both are within
the cutoff, giving a coordination number of 2 from a single residue.
This is the physically meaningful definition for geometry analysis because
coordination geometry is determined by the spatial arrangement of donor
atoms, not residues.

Coordination Geometry Classification
-------------------------------------
The geometry is classified by comparing the observed arrangement of donor
atoms around the metal to ideal polyhedra using two complementary metrics:

1.  **τ₄ index** (Yang et al., 2007) for 4-coordinate centres:
    τ₄ = (360° − α − β) / 141°
    where α and β are the two largest L–M–L angles.
    τ₄ = 0 → square planar; τ₄ = 1 → tetrahedral.

2.  **τ₅ index** (Addison et al., 1984) for 5-coordinate centres:
    τ₅ = (α − β) / 60°
    where α and β are the two largest L–M–L angles.
    τ₅ = 0 → square pyramidal; τ₅ = 1 → trigonal bipyramidal.

3.  For CN=6, the mean deviation from ideal octahedral angles (90° and 180°)
    is computed and the geometry is classified as octahedral, distorted
    octahedral, or trigonal prismatic based on thresholds from the
    inorganic chemistry literature.

4.  For CN=7 (the canonical EF-hand Ca²⁺ site), pentagonal bipyramidal
    geometry is reported along with the deviation from ideal angles.

5.  For CN=2, 3, 8, 9 the geometry label is derived directly from the
    coordination number using standard IUPAC nomenclature.

All L–M–L angles (ligand–metal–ligand) between donor atoms are computed
and stored.  The mean, minimum, maximum, and the two largest angles are
reported in the motif metadata CSV.

References
----------
Addison, A. W. et al. (1984). J. Chem. Soc. Dalton Trans., 1349–1356.
Yang, L. et al. (2007). Dalton Trans., 955–964.
Alvarez, S. et al. (2005). Coord. Chem. Rev. 249, 1693–1708.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
from Bio.PDB.Residue import Residue

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Donor atom names per residue type
# Only side-chain donor atoms are included; backbone amide N and carbonyl O
# are included when present since they can coordinate in some metal sites.
# ---------------------------------------------------------------------------
_DONOR_ATOMS: dict[str, list[str]] = {
    "ASP": ["OD1", "OD2"],
    "GLU": ["OE1", "OE2"],
    "HIS": ["ND1", "NE2"],
    "CYS": ["SG"],
    "SER": ["OG"],
    "THR": ["OG1"],
    "TYR": ["OH"],
    "ASN": ["OD1", "ND2"],
    "GLN": ["OE1", "NE2"],
    "MET": ["SD"],
    "LYS": ["NZ"],
    "ARG": ["NH1", "NH2"],
    "TRP": ["NE1"],
    # Backbone atoms — relevant in some coordination environments
    "GLY": ["O"],
}
# Backbone carbonyl O is a potential donor for any residue
_BACKBONE_DONOR = "O"


# ===========================================================================
# Data class
# ===========================================================================

@dataclass
class CoordinationResult:
    """Coordination analysis result for one metal ion site.

    Attributes
    ----------
    coordination_number : int
        Number of donor atoms within the distance cutoff.
    geometry : str
        Named coordination geometry (e.g. ``"Pentagonal bipyramidal"``).
    tau4 : float or None
        τ₄ index for 4-coordinate sites (0=square planar, 1=tetrahedral).
    tau5 : float or None
        τ₅ index for 5-coordinate sites (0=sq. pyramidal, 1=trig. bipyramidal).
    mean_angle : float
        Mean L–M–L angle across all donor-atom pairs (degrees).
    min_angle : float
        Minimum L–M–L angle (degrees).
    max_angle : float
        Maximum L–M–L angle (degrees).
    angle_largest_1 : float
        Largest L–M–L angle (degrees).
    angle_largest_2 : float
        Second-largest L–M–L angle (degrees).
    donor_atoms : str
        Comma-separated list of donor atom descriptions
        (e.g. ``"ASP12-OD1,GLU17-OE1"``).
    all_angles : list[float]
        All pairwise L–M–L angles in degrees.
    """
    coordination_number: int = 0
    geometry: str = "Unknown"
    tau4: float | None = None
    tau5: float | None = None
    mean_angle: float = 0.0
    min_angle: float = 0.0
    max_angle: float = 0.0
    angle_largest_1: float = 0.0
    angle_largest_2: float = 0.0
    donor_atoms: str = ""
    all_angles: list[float] = field(default_factory=list)


# ===========================================================================
# Donor atom extraction
# ===========================================================================

def _get_donor_atoms(
    binding_residues: list[Residue],
    ion: Residue,
    distance_cutoff: float,
) -> list[tuple[str, np.ndarray]]:
    """Return all donor atoms within distance_cutoff of the ion.

    Parameters
    ----------
    binding_residues : list[Residue]
        Residues identified as coordinating (from motif extraction).
    ion : Residue
        Metal ion HETATM residue.
    distance_cutoff : float
        Maximum atom-to-atom distance (Å).

    Returns
    -------
    list[tuple[str, np.ndarray]]
        Each entry is (label, coordinate) where label is e.g. ``"ASP12-OD1"``.
    """
    ion_coords = [a.get_vector().get_array() for a in ion.get_atoms()]
    if not ion_coords:
        return []
    # Use the centroid of multi-atom ions (should be single atom for metal ions)
    ion_pos = np.mean(ion_coords, axis=0)

    donors: list[tuple[str, np.ndarray]] = []

    for res in binding_residues:
        resname = res.get_resname().strip().upper()
        seqnum = res.id[1]

        # Collect candidate donor atom names for this residue type
        candidate_names = list(_DONOR_ATOMS.get(resname, []))
        # Always include backbone carbonyl O as a candidate
        if _BACKBONE_DONOR not in candidate_names:
            candidate_names.append(_BACKBONE_DONOR)

        for atom_name in candidate_names:
            if atom_name not in res:
                continue
            atom = res[atom_name]
            coord = atom.get_vector().get_array()
            dist = float(np.linalg.norm(coord - ion_pos))
            if dist <= distance_cutoff:
                label = f"{resname}{seqnum}-{atom_name}"
                donors.append((label, coord))

    return donors


# ===========================================================================
# Angle computation
# ===========================================================================

def _lml_angles(donors: list[tuple[str, np.ndarray]],
                ion_pos: np.ndarray) -> list[float]:
    """Compute all pairwise L–M–L angles between donor atoms.

    Parameters
    ----------
    donors : list[tuple[str, np.ndarray]]
        Donor atoms (label, coordinate).
    ion_pos : np.ndarray
        Metal ion position.

    Returns
    -------
    list[float]
        Angles in degrees, one per unique pair.
    """
    angles: list[float] = []
    coords = [c for _, c in donors]
    n = len(coords)

    for i in range(n):
        for j in range(i + 1, n):
            v1 = coords[i] - ion_pos
            v2 = coords[j] - ion_pos
            n1 = np.linalg.norm(v1)
            n2 = np.linalg.norm(v2)
            if n1 < 1e-6 or n2 < 1e-6:
                continue
            cos_angle = np.clip(np.dot(v1, v2) / (n1 * n2), -1.0, 1.0)
            angle = float(np.degrees(np.arccos(cos_angle)))
            angles.append(angle)

    return sorted(angles, reverse=True)


# ===========================================================================
# Geometry classification
# ===========================================================================

def _classify_geometry(cn: int, angles: list[float]) -> tuple[str, float | None, float | None]:
    """Classify coordination geometry from coordination number and L–M–L angles.

    Parameters
    ----------
    cn : int
        Coordination number.
    angles : list[float]
        All L–M–L angles in degrees, sorted descending.

    Returns
    -------
    tuple[str, float | None, float | None]
        (geometry_name, tau4_or_None, tau5_or_None)
    """
    if cn == 0:
        return "None", None, None
    if cn == 1:
        return "Monocoordinate", None, None
    if cn == 2:
        return "Linear" if (angles and angles[0] > 150) else "Bent", None, None
    if cn == 3:
        if not angles or len(angles) < 3:
            return "Trigonal planar", None, None
        # Trigonal planar: all three angles ≈ 120°, so the largest is ≤ 135°.
        # T-shaped: one angle approaches 180° (the two trans ligands) while
        # the other two approach 90°.  The largest angle is the discriminator.
        largest = angles[0]
        return ("T-shaped" if largest > 150.0 else "Trigonal planar"), None, None

    if cn == 4:
        if len(angles) < 2:
            return "4-coordinate (insufficient angles)", None, None
        alpha, beta = angles[0], angles[1]
        tau4 = round((360.0 - alpha - beta) / 141.0, 3)
        tau4 = float(np.clip(tau4, 0.0, 1.0))
        if tau4 >= 0.85:
            geom = "Tetrahedral"
        elif tau4 <= 0.15:
            geom = "Square planar"
        elif tau4 <= 0.40:
            geom = "Square planar (distorted)"
        else:
            geom = "Seesaw / distorted tetrahedral"
        return geom, tau4, None

    if cn == 5:
        if len(angles) < 2:
            return "5-coordinate (insufficient angles)", None, None
        alpha, beta = angles[0], angles[1]
        tau5 = round((alpha - beta) / 60.0, 3)
        tau5 = float(np.clip(tau5, 0.0, 1.0))
        if tau5 >= 0.75:
            geom = "Trigonal bipyramidal"
        elif tau5 <= 0.25:
            geom = "Square pyramidal"
        else:
            geom = "Trigonal bipyramidal (distorted)"
        return geom, None, tau5

    if cn == 6:
        if len(angles) < 3:
            return "Octahedral (insufficient angles)", None, None
        # For octahedral, ideal angles are 90° (×12) and 180° (×3).
        # Compute mean absolute deviation from nearest ideal.
        deviations = [min(abs(a - 90.0), abs(a - 180.0)) for a in angles]
        mad = float(np.mean(deviations))
        if mad < 10.0:
            geom = "Octahedral"
        elif mad < 20.0:
            geom = "Octahedral (distorted)"
        else:
            geom = "Trigonal prismatic"
        return geom, None, None

    if cn == 7:
        # For seven-coordinate sites the two most common geometries in
        # protein coordination chemistry are pentagonal bipyramidal (PBP)
        # and capped octahedral (CO).  We discriminate between them by
        # comparing the observed angle distribution to ideal reference angles.
        #
        # Ideal PBP: five equatorial donors at 72° apart, two axial donors
        # at 90° from each equatorial donor and 180° from each other.
        # The $\binom{7}{2} = 21$ pairwise angles consist of:
        #   5 × 72° (equatorial–equatorial adjacent)
        #   5 × 144° (equatorial–equatorial across)
        #   10 × 90° (equatorial–axial)
        #   1 × 180° (axial–axial)
        #
        # Ideal capped octahedral: one capping donor sits above a face of
        # an octahedron.  Angles cluster around 70°, 90°, 120°, and 150°.
        #
        # We use mean absolute deviation from each ideal distribution as
        # the discrimination metric.
        if len(angles) < 6:
            return "7-coordinate (insufficient angles)", None, None

        # Reference angle sets (frequency-weighted) for each geometry
        pbp_refs  = [72.]*5 + [144.]*5 + [90.]*10 + [180.]*1  # 21 angles
        co_refs   = [70.]*3 + [90.]*6  + [120.]*6 + [150.]*3 + [180.]*3  # 21 angles

        obs = sorted(angles)
        pbp_ref_s = sorted(pbp_refs)
        co_ref_s  = sorted(co_refs)
        obs_use = (obs + [0.]*(21 - len(obs)))[:21]
        mad_pbp = float(np.mean(np.abs(np.array(sorted(obs_use)) - np.array(pbp_ref_s))))
        mad_co  = float(np.mean(np.abs(np.array(sorted(obs_use)) - np.array(co_ref_s))))

        if mad_pbp <= mad_co:
            geom = ("Pentagonal bipyramidal"
                    if mad_pbp < 15.0
                    else "Pentagonal bipyramidal (distorted)")
        else:
            geom = ("Capped octahedral"
                    if mad_co < 15.0
                    else "Capped octahedral (distorted)")
        return geom, None, None

    if cn == 8:
        # Eight-coordinate geometry in proteins most commonly falls into
        # two families: square antiprismatic (SAP) and dodecahedral (DD).
        # Full discrimination requires continuous shape measure (CSM)
        # analysis (Alvarez et al., 2005), which is beyond the scope of
        # this distance-cutoff-based tool.  We therefore report the mean
        # L–M–L angle alongside the generic label to allow manual
        # inspection.  The two geometries can be distinguished post-hoc
        # using the SHAPE software (Alvarez et al.).
        mean_a = float(np.mean(angles)) if angles else 0.0
        geom = ("Square antiprismatic"
                if mean_a > 100.0
                else "Dodecahedral (bicapped trigonal prismatic)")
        return geom, None, None

    if cn == 9:
        # Nine-coordinate sites in proteins are rare and occur almost
        # exclusively for large lanthanide ions (La³⁺, Ce³⁺, Nd³⁺).
        # The near-universal geometry is tricapped trigonal prismatic (TTP).
        # No angular discrimination is implemented; the reported L–M–L
        # angles allow post-hoc verification using SHAPE or similar tools.
        return "Tricapped trigonal prismatic", None, None

    return f"{cn}-coordinate (unclassified)", None, None


# ===========================================================================
# Main entry point
# ===========================================================================

def analyse_coordination(
    binding_residues: list[Residue],
    ion: Residue,
    distance_cutoff: float,
) -> CoordinationResult:
    """Compute full coordination analysis for a metal ion site.

    Parameters
    ----------
    binding_residues : list[Residue]
        Residues within distance_cutoff of the ion (from motif extraction).
    ion : Residue
        Metal ion HETATM residue.
    distance_cutoff : float
        Distance cutoff used for binding residue identification (Å).

    Returns
    -------
    CoordinationResult
        Full coordination analysis result.
    """
    # Get ion centroid
    ion_coords = [a.get_vector().get_array() for a in ion.get_atoms()]
    if not ion_coords:
        log.debug("Ion %s has no atoms — cannot compute coordination", ion.id)
        return CoordinationResult()

    ion_pos = np.mean(ion_coords, axis=0)

    # Get donor atoms
    donors = _get_donor_atoms(binding_residues, ion, distance_cutoff)
    cn = len(donors)

    if cn == 0:
        return CoordinationResult(
            coordination_number=0,
            geometry="No donor atoms within cutoff",
        )

    # Compute all L–M–L angles
    angles = _lml_angles(donors, ion_pos)

    # Summary statistics
    mean_angle = float(np.mean(angles)) if angles else 0.0
    min_angle  = float(np.min(angles))  if angles else 0.0
    max_angle  = float(np.max(angles))  if angles else 0.0
    a1 = angles[0] if len(angles) > 0 else 0.0
    a2 = angles[1] if len(angles) > 1 else 0.0

    # Geometry classification
    geom, tau4, tau5 = _classify_geometry(cn, angles)

    donor_str = ",".join(label for label, _ in donors)

    log.debug(
        "Ion %s: CN=%d, geometry=%s, mean_angle=%.1f°, donors=%s",
        ion.get_resname(), cn, geom, mean_angle, donor_str,
    )

    return CoordinationResult(
        coordination_number=cn,
        geometry=geom,
        tau4=tau4,
        tau5=tau5,
        mean_angle=round(mean_angle, 2),
        min_angle=round(min_angle, 2),
        max_angle=round(max_angle, 2),
        angle_largest_1=round(a1, 2),
        angle_largest_2=round(a2, 2),
        donor_atoms=donor_str,
        all_angles=[round(a, 2) for a in angles],
    )
