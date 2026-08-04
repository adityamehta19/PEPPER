#!/usr/bin/env python3
"""
metal_motif_extractor.py
=========================
MetalMotifExtractor — main orchestration script.

A publication-quality, general-purpose tool for extracting metal-ion-binding
structural motifs (e.g. EF-hand Ca2+ loops) from RCSB PDB structures.

Usage
-----
    python metal_motif_extractor.py --metal ZN --resolution 2.0 \\
        --distance 4 --motif_length 15 --mode uniprot

    python metal_motif_extractor.py --excel proteins.xlsx --metal CA

See ``cli.py`` for the full list of command-line options, and ``config.yaml``
for the equivalent persistent configuration.

Pipeline Stages
---------------
1.  Parse CLI arguments and merge with ``config.yaml``.
2.  Initialise logging (console + ``logs/log.txt``).
3.  Acquire candidate PDB IDs (Excel / text file / automatic EF-hand search).
4.  Fetch metadata and filter by resolution / experimental method.
5.  Download structures in parallel.
6.  Extract metal-binding motifs from each downloaded structure.
7.  Write structure-level and motif-level metadata CSVs.
8.  Print a summary report.
"""

from __future__ import annotations

import logging
import sys
import time
from dataclasses import asdict
from pathlib import Path

import pandas as pd

from cli import build_parser, parse_args_to_overrides
from downloader import (
    download_structures,
    ids_from_excel,
    ids_from_interpro_efhand,
    ids_from_text,
)
from metadata import (
    fetch_entry_metadata,
    has_any_metal,
    passes_filter,
)
from motif import MotifRecord, extract_motifs_from_structure
from utils import ensure_dirs, load_config, setup_logging

log = logging.getLogger(__name__)


def acquire_pdb_ids(cfg: dict) -> list[str]:
    """Dispatch to the appropriate input-mode handler to obtain PDB IDs.

    Modes
    -----
    excel   — read PDB IDs from a single-column Excel file.
    text    — read PDB IDs from a plain text file (one ID per line).
    uniprot — automatically discover protein structures by querying the
              UniProt REST API with a configurable domain query string.
              The default query targets EF-hand proteins; any protein family
              with curated Swiss-Prot domain annotations can be targeted by
              changing ``uniprot_query`` in config.yaml or via --uniprot_query.

    Parameters
    ----------
    cfg : dict
        Merged configuration dictionary.

    Returns
    -------
    list[str]
        Candidate PDB IDs.
    """
    mode = cfg["input_mode"]

    if mode == "excel":
        if not cfg.get("input_file"):
            log.error("Excel mode selected but no input_file provided.")
            sys.exit(1)
        return ids_from_excel(cfg["input_file"])

    if mode == "text":
        if not cfg.get("input_file"):
            log.error("Text mode selected but no input_file provided.")
            sys.exit(1)
        return ids_from_text(cfg["input_file"])

    if mode == "uniprot":
        uq = cfg.get("uniprot_query", '(ft_domain:"EF-hand") AND (reviewed:true)')
        log.info("UniProt discovery mode: query = %s", uq)
        non_ca_metals = [m for m in cfg["metal"] if m not in ("CA", "MG")]
        if non_ca_metals and not cfg.get("any_metal", False) and "EF-hand" in uq:
            log.warning(
                "The UniProt query targets EF-hand proteins (Ca2+/Mg2+ binding) "
                "but --metal includes %s. Use --metal any to extract all metals "
                "from EF-hand structures, or supply a family-specific "
                "--uniprot_query for the protein family of interest.",
                non_ca_metals,
            )
        return ids_from_interpro_efhand(
            max_resolution=cfg["max_resolution"],
            allowed_methods=cfg["allowed_methods"],
            include_nmr=cfg["include_nmr"],
            uniprot_query=uq,
        )

    log.error("Unknown input_mode: '%s'. Valid modes: excel, text, uniprot", mode)
    sys.exit(1)


def filter_and_fetch_metadata(
    pdb_ids: list[str], cfg: dict
) -> tuple[list[str], dict[str, dict], int]:
    """Fetch metadata for each candidate and apply all filters serially.

    Filters applied in order:
    1. Resolution and experimental method (RCSB Data API).
    2. Metal presence — confirmed via the RCSB nonpolymer entity API.
       This removes apo forms and alternative-metal crystal forms that
       UniProt cross-references include simply because the same protein
       sequence was ever crystallised with the requested metal elsewhere.

    Parameters
    ----------
    pdb_ids : list[str]
        Candidate PDB IDs from the discovery step.
    cfg : dict
        Merged configuration dictionary.

    Returns
    -------
    tuple[list[str], dict[str, dict], int]
        Accepted PDB IDs, metadata mapping, and total count filtered out.
    """
    accepted: list[str] = []
    meta_map: dict[str, dict] = {}
    filtered_resolution = 0
    filtered_no_metal = 0
    failed = 0
    metals: list[str] = cfg["metal"]
    any_metal: bool = cfg.get("any_metal", False)

    for pdb_id in pdb_ids:
        meta = fetch_entry_metadata(pdb_id)
        if not meta:
            failed += 1
            continue

        if not passes_filter(
            meta, cfg["max_resolution"], cfg["allowed_methods"], cfg["include_nmr"]
        ):
            filtered_resolution += 1
            continue

        # When any_metal=True the user wants all structures from the
        # protein family regardless of which metal (if any) is present.
        # Skip the per-structure metal presence API call entirely.
        if not any_metal:
            if not has_any_metal(pdb_id, metals):
                filtered_no_metal += 1
                continue

        accepted.append(pdb_id)
        meta_map[pdb_id] = meta

    log.info(
        "Metadata filtering: %d accepted | %d filtered by resolution/method | "
        "%d filtered (metal not present in structure) | %d fetch failures",
        len(accepted), filtered_resolution, filtered_no_metal, failed,
    )
    return accepted, meta_map, filtered_resolution + filtered_no_metal




def write_structure_metadata_csv(
    meta_map: dict[str, dict],
    downloaded: dict[str, str],
    cfg: dict,
    motif_records: list[MotifRecord] | None = None,
) -> Path:
    """Write per-structure metadata to a CSV file.

    Parameters
    ----------
    meta_map : dict
        Mapping of ``pdb_id -> metadata dict`` (from RCSB Data API).
    downloaded : dict
        Mapping of ``pdb_id -> local file path`` for successfully
        downloaded structures.
    cfg : dict
        Merged configuration dictionary.
    motif_records : list[MotifRecord], optional
        All extracted motif records, used to report which of the
        *requested* metals were actually found in each structure (rather
        than just echoing the requested list back unconditionally).

    Returns
    -------
    Path
        Path to the written CSV file.
    """
    metadata_dir = Path(cfg["metadata_dir"])
    ensure_dirs(metadata_dir)
    out_path = metadata_dir / "structure_metadata.csv"

    # Build pdb_id -> set of metal names actually found in extracted motifs
    found_metals_by_pdb: dict[str, set[str]] = {}
    for rec in motif_records or []:
        found_metals_by_pdb.setdefault(rec.pdb_id, set()).add(rec.metal)

    rows = []
    for pdb_id, meta in meta_map.items():
        if pdb_id not in downloaded:
            continue
        found = found_metals_by_pdb.get(pdb_id)
        metal_str = ",".join(sorted(found)) if found else "none found"
        rows.append({
            "PDB ID": pdb_id,
            "Resolution": meta.get("resolution"),
            "Experimental Method": meta.get("method"),
            "Release Date": meta.get("release_date"),
            "Organism": meta.get("organism"),
            "Protein Name": meta.get("protein_name"),
            "Chain IDs": meta.get("chain_ids"),
            "Metal ion found": metal_str,
        })

    df = pd.DataFrame(rows)
    df.to_csv(out_path, index=False)
    log.info("Wrote structure metadata CSV (%d rows) → %s", len(rows), out_path)
    return out_path


def write_motif_metadata_csv(records: list[MotifRecord], cfg: dict) -> Path:
    """Write per-motif metadata to a CSV file including coordination and SS.

    Parameters
    ----------
    records : list[MotifRecord]
        All extracted motif records.
    cfg : dict
        Merged configuration dictionary.

    Returns
    -------
    Path
        Path to the written CSV file.
    """
    metadata_dir = Path(cfg["metadata_dir"])
    ensure_dirs(metadata_dir)
    out_path = metadata_dir / "motif_metadata.csv"

    rows = []
    for r in records:
        d = asdict(r)
        rows.append({
            "Motif ID":                  d["motif_id"],
            "PDB ID":                    d["pdb_id"],
            "Chain":                     d["chain"],
            "Metal":                     d["metal"],
            "Resolution":                d["resolution"],
            "Experimental Method":       d["method"],
            "Binding Residues":          d["binding_residues"],
            "Residue Numbers":           d["residue_numbers"],
            "Donor Atoms":               d["donor_atoms"],
            "Coordination Type":         d["coordination_type"],
            "Distance Cutoff":           d["distance_cutoff"],
            "Motif Length":              d["motif_length"],
            "Sequence":                  d["sequence"],
            "Output Filename":           d["output_file"],
            # --- Coordination analysis ---
            "Coordination Number":       d["coordination_number"],
            "Coordination Geometry":     d["coordination_geometry"],
            "tau4":                      d["tau4"],
            "tau5":                      d["tau5"],
            "Mean LML Angle (deg)":      d["mean_lml_angle"],
            "Min LML Angle (deg)":       d["min_lml_angle"],
            "Max LML Angle (deg)":       d["max_lml_angle"],
            "Largest Angle 1 (deg)":     d["angle_largest_1"],
            "Largest Angle 2 (deg)":     d["angle_largest_2"],
            # --- Secondary structure summary ---
            "SS Method":                 d["ss_method"],
            "SS Sequence":               d["ss_sequence"],
            "Helix Fraction":            d["ss_helix_fraction"],
            "Strand Fraction":           d["ss_strand_fraction"],
            "Coil Fraction":             d["ss_coil_fraction"],
        })

    df = pd.DataFrame(rows)
    df.to_csv(out_path, index=False)
    log.info(
        "Wrote motif metadata CSV (%d rows, %d columns) → %s",
        len(rows), len(df.columns), out_path,
    )
    return out_path


def write_secondary_structure_csv(records: list[MotifRecord], cfg: dict) -> Path:
    """Write per-residue secondary structure data for all motifs.

    Produces a separate ``secondary_structure.csv`` in the metadata directory
    with one row per residue per motif — i.e. for a dataset of 400 motifs
    each with a 15-residue window, this file will have 6000 rows.

    Columns: Motif ID, PDB ID, Chain, Metal, Position (1-based within window),
    Residue Name, Residue Number, One Letter, SS 8-state, SS 3-state,
    SS Full Name, Phi (deg), Psi (deg), SS Method.

    Parameters
    ----------
    records : list[MotifRecord]
        All extracted motif records (must have ``_ss_records`` attribute
        attached by :func:`motif.extract_motifs_from_structure`).
    cfg : dict
        Merged configuration dictionary.

    Returns
    -------
    Path
        Path to the written CSV file.
    """
    metadata_dir = Path(cfg["metadata_dir"])
    ensure_dirs(metadata_dir)
    out_path = metadata_dir / "secondary_structure.csv"

    rows = []
    for rec in records:
        ss_records = getattr(rec, "_ss_records", [])
        for ss in ss_records:
            rows.append({
                "Motif ID":        rec.motif_id,
                "PDB ID":          rec.pdb_id,
                "Chain":           rec.chain,
                "Metal":           rec.metal,
                "Position":        ss.position,
                "Residue Name":    ss.residue_name,
                "Residue Number":  ss.residue_number,
                "One Letter":      ss.one_letter,
                "SS 8-state":      ss.ss_eight,
                "SS 3-state":      ss.ss_three,
                "SS Full Name":    ss.ss_fullname,
                "Phi (deg)":       ss.phi,
                "Psi (deg)":       ss.psi,
                "SS Method":       ss.method,
            })

    df = pd.DataFrame(rows)
    df.to_csv(out_path, index=False)
    log.info(
        "Wrote secondary structure CSV (%d rows) → %s", len(rows), out_path
    )
    return out_path


    return out_path


def print_summary(
    total_candidates: int,
    downloaded: dict[str, str],
    filtered_resolution: int,
    failed_downloads: int,
    metal: list[str],
    structures_with_metal: int,
    motif_records: list[MotifRecord],
    motif_dir: str,
) -> None:
    """Print a final summary report to stdout and the log.

    Parameters
    ----------
    total_candidates : int
        Total number of candidate PDB IDs considered.
    downloaded : dict
        Mapping of successfully downloaded structures.
    filtered_resolution : int
        Number of structures excluded by resolution/method filtering.
    failed_downloads : int
        Number of structures that failed to download.
    metal : list[str]
        Metal ion residue name(s) requested for this run.
    structures_with_metal : int
        Number of downloaded structures yielding at least one motif for
        any of the requested metals.
    motif_records : list[MotifRecord]
        All extracted motif records.
    motif_dir : str
        Output directory for motif PDB files.
    """
    metal_label = "/".join(metal)

    # Per-metal motif and structure breakdown
    motifs_by_metal: dict[str, int] = {}
    structures_by_metal: dict[str, set[str]] = {}
    for rec in motif_records:
        motifs_by_metal[rec.metal] = motifs_by_metal.get(rec.metal, 0) + 1
        structures_by_metal.setdefault(rec.metal, set()).add(rec.pdb_id)

    per_metal_lines = ""
    for m in metal:
        n_motifs = motifs_by_metal.get(m, 0)
        n_structs = len(structures_by_metal.get(m, set()))
        per_metal_lines += f"  {m:<4} → {n_structs} structures, {n_motifs} motifs\n"

    summary = (
        "\n=========================================\n"
        f"Total proteins processed : {total_candidates}\n"
        f"Downloaded structures    : {len(downloaded)}\n"
        f"Filtered by resolution   : {filtered_resolution}\n"
        f"Failed downloads         : {failed_downloads}\n"
        f"Metals requested         : {metal_label}\n"
        f"Structures with any motif: {structures_with_metal}\n"
        f"{per_metal_lines}"
        f"Motifs extracted (total) : {len(motif_records)}\n"
        f"Output directory         : {motif_dir}/\n"
        "=========================================\n"
    )
    print(summary)
    log.info(summary)


def main() -> None:
    """Entry point: parse arguments, run the pipeline, print summary."""
    parser = build_parser()
    args = parser.parse_args()
    overrides = parse_args_to_overrides(args)
    cfg = load_config(args.config, overrides=overrides)

    setup_logging(cfg["log_dir"])
    log.info("=" * 60)
    log.info("MetalMotifExtractor starting")
    log.info("Configuration: %s", cfg)
    log.info("=" * 60)

    ensure_dirs(
        cfg["output_dir"], cfg["pdb_dir"], cfg["motif_dir"],
        cfg["metadata_dir"], cfg["log_dir"],
    )

    start_time = time.time()

    # ------------------------------------------------------------------
    # Step 1: Acquire candidate PDB IDs
    # ------------------------------------------------------------------
    candidate_ids = acquire_pdb_ids(cfg)
    total_candidates = len(candidate_ids)
    if total_candidates == 0:
        log.error("No PDB IDs were acquired. Exiting.")
        sys.exit(1)
    log.info("Acquired %d candidate PDB IDs", total_candidates)

    # ------------------------------------------------------------------
    # Step 2: Metadata fetch + resolution/method/metal-presence filtering
    # ------------------------------------------------------------------
    filtered_ids, meta_map, filtered_resolution_count = filter_and_fetch_metadata(
        candidate_ids, cfg
    )

    if not filtered_ids:
        log.error(
            "No structures passed filtering (resolution/method/metal presence). "
            "Exiting."
        )
        sys.exit(1)

    # ------------------------------------------------------------------
    # Step 3: Parallel download
    # ------------------------------------------------------------------
    downloaded = download_structures(
        filtered_ids,
        pdb_dir=cfg["pdb_dir"],
        use_biological_assembly=cfg["use_biological_assembly"],
        max_workers=cfg["max_workers"],
    )
    failed_downloads = len(filtered_ids) - len(downloaded)

    # ------------------------------------------------------------------
    # Step 4: Motif extraction
    # ------------------------------------------------------------------
    all_records: list[MotifRecord] = []
    structures_with_metal = 0

    for pdb_id, pdb_file in downloaded.items():
        records = extract_motifs_from_structure(
            pdb_id=pdb_id,
            pdb_file=pdb_file,
            metal=cfg["metal"],
            distance_cutoff=cfg["distance_cutoff"],
            motif_length=cfg["motif_length"],
            motif_dir=cfg["motif_dir"],
            structure_meta=meta_map.get(pdb_id, {}),
        )
        if records:
            structures_with_metal += 1
            all_records.extend(records)

    # ------------------------------------------------------------------
    # Step 5: Metadata CSV export
    # ------------------------------------------------------------------
    write_structure_metadata_csv(meta_map, downloaded, cfg, motif_records=all_records)
    write_motif_metadata_csv(all_records, cfg)
    write_secondary_structure_csv(all_records, cfg)

    # ------------------------------------------------------------------
    # Step 6: Summary
    # ------------------------------------------------------------------
    elapsed = time.time() - start_time
    log.info("Pipeline completed in %.1f seconds", elapsed)

    print_summary(
        total_candidates=total_candidates,
        downloaded=downloaded,
        filtered_resolution=filtered_resolution_count,
        failed_downloads=failed_downloads,
        metal=cfg["metal"],
        structures_with_metal=structures_with_metal,
        motif_records=all_records,
        motif_dir=cfg["motif_dir"],
    )


if __name__ == "__main__":
    main()
