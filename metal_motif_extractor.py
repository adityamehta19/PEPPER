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
3.  Acquire candidate PDB IDs (Excel / text file / UniProt discovery).
    In UniProt mode, resolution/method filtering already happens here.
4.  Metadata filtering.  Metadata fetched in stage 3 is REUSED (nothing is
    fetched or filtered twice).  For Excel/text input the resolution/method
    filter is applied here.  The metal-presence filter (does the entry
    contain a requested metal?) is applied here ONLY when specific metals
    are requested; with ``--metal any`` it is skipped.  Every removal is
    counted under its own reason.
5.  Download structures in parallel (PDB format, mmCIF for assemblies).
6.  Extract metal-binding motifs (and coordinated waters) from each
    downloaded structure, model by model.
7.  Write structure-level, motif-level and per-residue CSVs, plus
    ``run_accounting.json`` and ``structure_outcomes.csv``.
8.  Print a summary report in which every count has its own label.
"""

from __future__ import annotations

import json
import logging
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict
from pathlib import Path

import pandas as pd

from cli import build_parser, parse_args_to_overrides
from downloader import (
    DiscoveryResult,
    discover_from_uniprot,
    download_structures,
    ids_from_excel,
    ids_from_text,
)
from metadata import (
    CLUSTER_ONLY_SEEN,
    enrich_entity_metadata,
    fetch_metadata_parallel,
    filter_reason,
    has_any_metal,
)
from motif import MotifRecord, extract_motifs_with_status
from utils import ConfigError, ensure_dirs, load_config, setup_logging

log = logging.getLogger(__name__)


def acquire_pdb_ids(cfg: dict) -> tuple[list[str], DiscoveryResult | None]:
    """Dispatch to the appropriate input-mode handler to obtain PDB IDs.

    Modes
    -----
    excel   - read PDB IDs from a single-column Excel file.
    text    - read PDB IDs from a plain text file (one ID per line).
    uniprot - discover protein structures by querying the UniProt REST API
              with a configurable domain query string.

    Returns
    -------
    (list[str], DiscoveryResult or None)
        Candidate PDB IDs, and (UniProt mode only) the discovery result that
        carries the already-fetched metadata and the stage counts.  For the
        UniProt mode the returned IDs are the ones that already passed
        resolution/method filtering.
    """
    mode = cfg["input_mode"]

    if mode == "excel":
        if not cfg.get("input_file"):
            log.error("Excel mode selected but no input_file provided.")
            sys.exit(1)
        return ids_from_excel(cfg["input_file"]), None

    if mode == "text":
        if not cfg.get("input_file"):
            log.error("Text mode selected but no input_file provided.")
            sys.exit(1)
        return ids_from_text(cfg["input_file"]), None

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
        disc = discover_from_uniprot(
            max_resolution=cfg["max_resolution"],
            allowed_methods=cfg["allowed_methods"],
            include_nmr=cfg["include_nmr"],
            uniprot_query=uq,
            max_workers=cfg["max_workers"],
        )
        return disc.ids, disc

    log.error("Unknown input_mode: '%s'. Valid modes: excel, text, uniprot", mode)
    sys.exit(1)


def filter_and_fetch_metadata(
    pdb_ids: list[str],
    cfg: dict,
    discovery: DiscoveryResult | None = None,
) -> tuple[list[str], dict[str, dict], dict[str, int]]:
    """Apply the remaining filters and complete the metadata.

    What is done, and what is NOT repeated
    --------------------------------------
    * **UniProt mode** (``discovery`` given): resolution and method were
      already filtered during discovery and the metadata was kept, so it is
      neither re-fetched nor re-filtered here.
    * **Excel / text mode**: metadata is fetched once (in parallel) and the
      resolution/method filter is applied here.
    * **Metal-presence filter**: applied only when specific metals are
      requested (``any_metal`` is False).  It removes structures whose entry
      does not contain any requested metal (apo forms, other-metal crystal
      forms).  With ``--metal any`` it is skipped entirely.  The count of
      structures removed by this filter is reported under its own name;
      it is never merged into the resolution count.
    * Finally, protein name / organism / chain IDs are fetched (from the
      polymer-entity endpoints) for the accepted structures only.

    Returns
    -------
    (accepted ids, metadata map for accepted ids, stats)
        ``stats`` has separate counters: ``metadata_failed``,
        ``removed_no_resolution``, ``removed_method``, ``removed_resolution``
        (only counted here for excel/text mode; for UniProt mode they come
        from discovery), ``removed_no_metal``, ``metal_filter_applied``,
        ``passed_filters``.
    """
    stats: dict[str, int] = {
        "metadata_failed": 0,
        "removed_no_resolution": 0,
        "removed_method": 0,
        "removed_resolution": 0,
        "removed_no_metal": 0,
        "removed_cluster_only": 0,
        "metal_filter_applied": 0,
        "passed_filters": 0,
    }
    metals: list[str] = cfg["metal"]
    any_metal: bool = cfg.get("any_metal", False)

    # ---- resolution / method -------------------------------------------------
    if discovery is not None:
        # Already filtered during discovery; reuse cached metadata.
        candidates = list(pdb_ids)
        meta_all = {pid: discovery.meta_map[pid] for pid in candidates if pid in discovery.meta_map}
    else:
        fetched = fetch_metadata_parallel(
            pdb_ids, max_workers=cfg["max_workers"], include_entities=False
        )
        stats["metadata_failed"] = len(pdb_ids) - len(fetched)
        candidates = []
        meta_all = {}
        for pid in pdb_ids:
            meta = fetched.get(pid)
            if not meta:
                continue
            reason = filter_reason(
                meta, cfg["max_resolution"], cfg["allowed_methods"], cfg["include_nmr"]
            )
            if reason == "":
                candidates.append(pid)
                meta_all[pid] = meta
            elif reason in ("no_resolution", "method", "resolution"):
                stats[f"removed_{reason}"] += 1

    # ---- metal presence (only when specific metals are requested) ------------
    if any_metal:
        accepted = list(candidates)
        log.info(
            "Metal-presence filter: OFF (--metal any); %d structures go to extraction",
            len(accepted),
        )
    else:
        stats["metal_filter_applied"] = 1
        CLUSTER_ONLY_SEEN.clear()
        with ThreadPoolExecutor(max_workers=cfg["max_workers"]) as pool:
            keep_flags = list(pool.map(lambda pid: has_any_metal(pid, metals), candidates))
        accepted = [pid for pid, keep in zip(candidates, keep_flags) if keep]
        stats["removed_no_metal"] = len(candidates) - len(accepted)
        removed = {pid for pid, keep in zip(candidates, keep_flags) if not keep}
        stats["removed_cluster_only"] = sum(1 for pid in removed if pid.upper() in CLUSTER_ONLY_SEEN)
        log.info(
            "Metal-presence filter: ON (%s); %d structures removed because none "
            "of the requested metals is present, %d kept",
            "/".join(metals), stats["removed_no_metal"], len(accepted),
        )

    meta_map = {pid: meta_all[pid] for pid in accepted}
    stats["passed_filters"] = len(accepted)

    # ---- protein name / organism / chains for accepted structures ------------
    if meta_map:
        enrich_entity_metadata(meta_map, max_workers=cfg["max_workers"])

    return accepted, meta_map, stats


def write_structure_metadata_csv(
    meta_map: dict[str, dict],
    downloaded: dict[str, str],
    cfg: dict,
    motif_records: list[MotifRecord] | None = None,
    outcomes: dict[str, str] | None = None,
) -> Path:
    """Write per-structure metadata to a CSV file.

    Parameters
    ----------
    meta_map : dict
        Mapping of ``pdb_id -> metadata dict`` (from the RCSB Data API).
    downloaded : dict
        Mapping of ``pdb_id -> local file path`` for downloaded structures.
    cfg : dict
        Merged configuration dictionary.
    motif_records : list[MotifRecord], optional
        All extracted motif records, used to report which of the requested
        metals were actually found in each structure.
    outcomes : dict, optional
        ``pdb_id -> extraction status`` (see ``motif.STATUS_*``).

    Returns
    -------
    Path
        Path to the written CSV file.
    """
    metadata_dir = Path(cfg["metadata_dir"])
    ensure_dirs(metadata_dir)
    out_path = metadata_dir / "structure_metadata.csv"

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
            "Extraction Status": (outcomes or {}).get(pdb_id, ""),
        })

    df = pd.DataFrame(rows)
    df.to_csv(out_path, index=False)
    log.info("Wrote structure metadata CSV (%d rows) -> %s", len(rows), out_path)
    return out_path


def write_motif_metadata_csv(records: list[MotifRecord], cfg: dict) -> Path:
    """Write per-motif metadata to a CSV file.

    Coordination columns
    --------------------
    ``Coordination Number``
        The number of DONOR ATOMS bound to the ion: protein donor atoms +
        coordinated water oxygens (unless ``--no_water``) + N/O/S atoms of
        non-protein ligands (unless ``--no_ligands``).  A bidentate
        carboxylate counts two.  With both options off it is simply the
        number of protein donor atoms of the window's chain.
    ``CN (protein donors)``, ``Water Count``, ``Ligand Atom Count``
        The three components of the coordination number.
    ``All Donor Atoms``
        Labels of every donor atom counted (``ASP8:OD1``, ``HOH:W:101:O``,
        ``ATP:L:501:O2G``, ...), so the number can always be checked by eye.
    ``Coordination Geometry``
        The ideal polyhedron that best fits ALL donor atoms (template
        fitting).  ``Geometry RMSD (A)`` is the RMSD of that fit;
        ``Geometry Second Best`` / ``Geometry Second RMSD (A)`` / ``Geometry
        Margin (A)`` give the runner-up and the gap, so marginal assignments
        are visible.  ``Coordination Geometry (protein donors only)`` and its
        RMSD show the result without water/ligands.
    L-M-L angles (mean, minimum, maximum and the two largest)
        Computed from all donor atoms.
    ``Donors in Other Chains`` / ``Other-Chain Donor Atoms``
        Number and labels (``B:ASP5:OD1``) of protein donor atoms of the same
        ion that lie in OTHER chains.  They are not part of this record's
        coordination number or geometry; a value above 0 marks an inter-chain
        site whose CN is therefore underestimated.  Every donor atom of the
        site can be recovered from the records of that ion.
    ``Window Numbering Span`` / ``Window Missing Residues``
        Last minus first residue number of the window, and how many residue
        numbers inside it are not modelled (gaps), i.e. where the window is
        not contiguous in the numbering.
    Water and ligand columns: see ``Water Count`` / ``Ligand Atom Count``
    and the ``Water ...`` / ``Ligand ...`` columns (cutoffs, identifiers,
    distances).  ``Common Additive Residues`` flags ligand residues named
    like crystallisation additives (flag only; nothing is excluded).
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
            "Model":                     d["model_id"],
            "Chain":                     d["chain"],
            "Metal":                     d["metal"],
            "Ion Chain":                 d["ion_chain"],
            "Ion Residue Number":        d["ion_resnum"],
            "Ion Serial":                d["ion_serial"],
            "Resolution":                d["resolution"],
            "Experimental Method":       d["method"],
            "Binding Residues":          d["binding_residues"],
            "Residue Numbers":           d["residue_numbers"],
            "Donor Atoms":               d["donor_atoms"],
            "Coordination Type":         d["coordination_type"],
            "Distance Cutoff":           d["distance_cutoff"],
            "Motif Length":              d["motif_length"],
            "Requested Motif Length":    d["requested_length"],
            "Binding Span (residues)":   d["binding_span"],
            "Window Exceeds Requested":  d["window_exceeds_requested"],
            "Window Numbering Span":     d["window_numbering_span"],
            "Window Missing Residues":   d["window_missing_residues"],
            "Sequence":                  d["sequence"],
            "Output Filename":           d["output_file"],
            # --- Coordination analysis (protein donors) ---
            "Coordination Number":       d["coordination_number"],
            "Expected Max CN":           d["expected_max_cn"],
            "CN Exceeds Expected Max":   d["cn_exceeds_max"],
            "Mean Donor Distance (A)":   d["mean_donor_distance"],
            "Min Donor Distance (A)":    d["min_donor_distance"],
            "Max Donor Distance (A)":    d["max_donor_distance"],
            "Ion Occupancy":             d["ion_occupancy"],
            "Ion B-factor":              d["ion_bfactor"],
            "CN (protein donors)":       d["cn_protein"],
            "Coordination Geometry":     d["coordination_geometry"],
            "Geometry RMSD (A)":         d["geometry_rmsd"],
            "Geometry Second Best":      d["geometry_second"],
            "Geometry Second RMSD (A)":  d["geometry_second_rmsd"],
            "Geometry Margin (A)":       d["geometry_margin"],
            "Coordination Geometry (protein donors only)": d["geometry_protein_only"],
            "Geometry RMSD (protein donors only) (A)": d["geometry_rmsd_protein_only"],
            "Mean LML Angle (deg)":      d["mean_lml_angle"],
            "Min LML Angle (deg)":       d["min_lml_angle"],
            "Max LML Angle (deg)":       d["max_lml_angle"],
            "Largest Angle 1 (deg)":     d["angle_largest_1"],
            "Largest Angle 2 (deg)":     d["angle_largest_2"],
            # --- Coordinated water ---
            "Water Cutoff (A)":          d["water_cutoff"],
            "Water Count":               d["water_count"],
            "Water IDs":                 d["water_ids"],
            "Water Distances (A)":       d["water_distances"],
            # --- Non-protein ligands ---
            "Ligand Cutoff (A)":         d["ligand_cutoff"],
            "Ligand Atom Count":         d["ligand_atom_count"],
            "Ligand Residues":           d["ligand_residues"],
            "Ligand Atoms":              d["ligand_atoms"],
            "Ligand Distances (A)":      d["ligand_distances"],
            "Common Additive Residues":  d["ligand_additives"],
            "All Donor Atoms":           d["all_donor_atoms"],
            "Donors in Other Chains":    d["donors_other_chains"],
            "Other-Chain Donor Atoms":   d["other_chain_donor_atoms"],
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
        "Wrote motif metadata CSV (%d rows, %d columns) -> %s",
        len(rows), len(df.columns), out_path,
    )
    return out_path


def write_secondary_structure_csv(records: list[MotifRecord], cfg: dict) -> Path:
    """Write per-residue secondary structure data for all motifs.

    One row per residue per motif.  Columns: Motif ID, PDB ID, Chain, Metal,
    Position (1-based within window), Residue Name, Residue Number, One
    Letter, SS 8-state, SS 3-state, SS Full Name, Phi (deg), Psi (deg),
    SS Method.
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
    log.info("Wrote secondary structure CSV (%d rows) -> %s", len(rows), out_path)
    return out_path


def write_accounting_files(
    accounting: dict[str, int],
    outcomes: dict[str, str],
    motif_counts: dict[str, int],
    discovery: DiscoveryResult | None,
    cfg: dict,
    timings: dict[str, float] | None = None,
) -> None:
    """Write run_accounting.json, structure_outcomes.csv and (UniProt mode)
    uniprot_pdb_xrefs.csv so every number in the summary can be audited."""
    metadata_dir = Path(cfg["metadata_dir"])
    ensure_dirs(metadata_dir)

    with (metadata_dir / "run_accounting.json").open("w", encoding="utf-8") as fh:
        json.dump(
            {"input_mode": cfg["input_mode"], "metal": cfg["metal"],
             "any_metal": cfg.get("any_metal", False),
             "max_resolution": cfg["max_resolution"],
             "water_cutoff": cfg.get("water_cutoff"),
             "include_water": cfg.get("include_water"),
             "ligand_cutoff": cfg.get("ligand_cutoff"),
             "include_ligands": cfg.get("include_ligands"),
             "geometry_max_rmsd": cfg.get("geometry_max_rmsd"),
             "metal_cutoffs": cfg.get("metal_cutoffs"),
             "use_biological_assembly": cfg.get("use_biological_assembly"),
             "asu_format": cfg.get("asu_format"),
             "write_findgeo_input": cfg.get("write_findgeo_input"),
             "model_policy": cfg.get("model_policy"),
             "counts": accounting},
            fh, indent=2,
        )

    pd.DataFrame(
        [{"PDB ID": pid, "Extraction Status": st, "Motifs": motif_counts.get(pid, 0),
          "Extraction Seconds": (timings or {}).get(pid)}
         for pid, st in sorted(outcomes.items())]
    ).to_csv(metadata_dir / "structure_outcomes.csv", index=False)

    if discovery is not None:
        pd.DataFrame(
            [{"UniProt Accession": acc, "PDB IDs": ",".join(pdbs), "N PDB IDs": len(pdbs)}
             for acc, pdbs in sorted(discovery.accession_to_pdb.items())]
        ).to_csv(metadata_dir / "uniprot_pdb_xrefs.csv", index=False)

    log.info("Wrote run_accounting.json and structure_outcomes.csv to %s", metadata_dir)


def cn_check(records: list[MotifRecord], limits: dict[str, int]) -> dict[str, dict]:
    """Compare the coordination numbers with the highest value each metal can have.

    A metal whose sites exceed its limit has been given a distance cut-off that
    is too high for the coordination that is physically possible (second-shell
    atoms are being counted).  Metals with no known limit are reported as such;
    no assumption is made for them.
    """
    out: dict[str, dict] = {}
    for metal in sorted({r.metal for r in records}):
        recs = [r for r in records if r.metal == metal]
        cns = [r.coordination_number for r in recs]
        dist: dict[int, int] = {}
        for c in cns:
            dist[c] = dist.get(c, 0) + 1
        limit = limits.get(metal)
        over = [c for c in cns if limit is not None and c > limit]
        out[metal] = {
            "sites": len(recs), "limit": limit, "highest_cn": max(cns),
            "n_over": len(over),
            "pct_over": round(100.0 * len(over) / len(recs), 1),
            "cn_distribution": dict(sorted(dist.items())),
            "status": "no_limit" if limit is None else ("too_high" if over else "ok"),
        }
    return out


def format_cn_check(check: dict[str, dict], cfg: dict) -> str:
    """Text block for the run summary: cut-offs used and the cut-off check."""
    if not check:
        return ""
    per_metal = cfg.get("metal_cutoffs") or {}
    water = f"{cfg['water_cutoff']}" if cfg.get("include_water") else "not searched"
    lig = f"{cfg['ligand_cutoff']}" if cfg.get("include_ligands") else "not searched"
    lines = [
        "Distance cut-offs used (A)   : "
        f"protein {cfg['distance_cutoff']}, water {water}, ligand {lig}"
        + (f", per metal {per_metal}" if per_metal else ""),
        "Coordination number check (is the cut-off too high for the metal?):",
    ]
    warnings = []
    for metal, c in check.items():
        dist = " ".join(f"{k}:{v}" for k, v in c["cn_distribution"].items())
        if c["status"] == "no_limit":
            lines.append(f"  {metal:<4} sites {c['sites']:>6} | highest CN {c['highest_cn']:>2} | "
                         f"no known limit, check skipped (set one with --max_cn {metal}=N)")
        else:
            flag = "OK" if c["status"] == "ok" else "TOO HIGH"
            lines.append(f"  {metal:<4} sites {c['sites']:>6} | limit {c['limit']:>2} | above limit "
                         f"{c['n_over']} ({c['pct_over']}%) | highest CN {c['highest_cn']:>2} | {flag}")
        lines.append(f"       CN distribution: {dist}")
        if c["status"] == "too_high":
            warnings.append(
                f"  *** WARNING: the distance cut-off ({cfg['distance_cutoff']} A) is too high for the "
                f"coordination physically possible for {metal}: {c['n_over']} of {c['sites']} sites "
                f"({c['pct_over']}%) have more than {c['limit']} donor atoms (highest {c['highest_cn']}). "
                f"Lower --distance (or set --metal_cutoff {metal}=<smaller value>; also check "
                f"--water_cutoff and --ligand_cutoff) and run again."
            )
    return "\n".join(lines + warnings) + "\n"


def length_summary(records: list[MotifRecord], requested: int) -> dict:
    """Summarise the lengths of the extracted motifs against the requested length.

    Returns counts and percentages of motifs that are exactly the requested
    length, longer than requested (the window was expanded to contain all
    protein donor residues) and shorter than requested (the chain is shorter
    than the window), the minimum/median/mean/maximum length, and a
    distribution table: one row per length up to the requested length, one row
    for each of the next 8 lengths, then bins of doubling width.
    """
    lengths = sorted(r.motif_length for r in records)
    n = len(lengths)
    if n == 0:
        return {"requested": requested, "total": 0}
    eq = sum(1 for x in lengths if x == requested)
    over = sum(1 for x in lengths if x > requested)
    under = sum(1 for x in lengths if x < requested)
    pct = lambda k: round(100.0 * k / n, 1)

    rows: list[tuple[str, int]] = []
    below = {}
    for x in lengths:
        if x <= requested:
            below[x] = below.get(x, 0) + 1
    for length in sorted(below):
        rows.append((str(length), below[length]))
    above = [x for x in lengths if x > requested]
    if above:
        top = max(above)
        # one row per length for the first 8 lengths above the requested one,
        # then bins that double in width (10, 20, 40, ...) so any maximum fits
        lo, width = requested + 1, 1
        while lo <= top:
            hi = lo + width - 1
            count = sum(1 for x in above if lo <= x <= hi)
            if count:
                rows.append((f"{lo}" if lo == hi else f"{lo}-{hi}", count))
            lo = hi + 1
            width = 1 if lo <= requested + 8 else (10 if lo == requested + 9 else width * 2)
    mid = n // 2
    median = lengths[mid] if n % 2 else (lengths[mid - 1] + lengths[mid]) / 2
    return {
        "requested": requested, "total": n,
        "equal_n": eq, "equal_pct": pct(eq),
        "exceeded_n": over, "exceeded_pct": pct(over),
        "shorter_n": under, "shorter_pct": pct(under),
        "min": lengths[0], "median": median,
        "mean": round(sum(lengths) / n, 2), "max": lengths[-1],
        "distribution": [{"length": lab, "motifs": cnt, "pct": pct(cnt)} for lab, cnt in rows],
    }


def format_length_summary(ls: dict) -> str:
    """Text block for the run summary."""
    if not ls.get("total"):
        return ""
    lines = [
        f"Requested motif length       : {ls['requested']}",
        f"  extracted at requested len : {ls['equal_n']} ({ls['equal_pct']}%)",
        f"  exceeded requested length  : {ls['exceeded_n']} ({ls['exceeded_pct']}%)"
        + (f"   (longest {ls['max']})" if ls["exceeded_n"] else ""),
        f"  shorter than requested     : {ls['shorter_n']} ({ls['shorter_pct']}%)"
        + ("   (chain shorter than the window)" if ls["shorter_n"] else ""),
        f"Motif length (residues)      : min {ls['min']}, median {ls['median']}, "
        f"mean {ls['mean']}, max {ls['max']}",
        "Length distribution:",
        "    length    motifs       %",
    ]
    peak = max(d["motifs"] for d in ls["distribution"])
    for d in ls["distribution"]:
        bar = "#" * max(1, round(30 * d["motifs"] / peak))
        lines.append(f"  {d['length']:>8}  {d['motifs']:>8}  {d['pct']:>6}%  {bar}")
    return "\n".join(lines) + "\n"


def write_cn_check(records: list[MotifRecord], cfg: dict) -> None:
    """Store the cut-offs used and the coordination-number check in run_accounting.json."""
    metadata_dir = Path(cfg["metadata_dir"])
    ensure_dirs(metadata_dir)
    path = metadata_dir / "run_accounting.json"
    data = json.load(path.open()) if path.exists() else {}
    data["cutoffs_A"] = {
        "protein": cfg["distance_cutoff"],
        "water": cfg["water_cutoff"] if cfg.get("include_water") else None,
        "ligand": cfg["ligand_cutoff"] if cfg.get("include_ligands") else None,
        "per_metal_protein": cfg.get("metal_cutoffs") or {},
    }
    data["cn_check"] = {m: {k: v for k, v in c.items()}
                        for m, c in cn_check(records, cfg["max_coordination_resolved"]).items()}
    data["max_coordination_used"] = cfg["max_coordination_resolved"]
    with path.open("w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2)


def write_length_distribution(records: list[MotifRecord], cfg: dict) -> None:
    """Write the motif-length summary to run_accounting.json and a CSV."""
    ls = length_summary(records, cfg["motif_length"])
    metadata_dir = Path(cfg["metadata_dir"])
    ensure_dirs(metadata_dir)
    path = metadata_dir / "run_accounting.json"
    data = json.load(path.open()) if path.exists() else {}
    data["motif_length"] = {k: v for k, v in ls.items() if k != "distribution"}
    with path.open("w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2)
    if ls.get("total"):
        pd.DataFrame(ls["distribution"]).rename(
            columns={"length": "Motif length (residues)", "motifs": "Motifs", "pct": "Percent"}
        ).to_csv(metadata_dir / "motif_length_distribution.csv", index=False)


def ss_method_summary(records: list[MotifRecord]) -> str:
    """Text block: how many motifs got their secondary structure from each method."""
    if not records:
        return ""
    counts: dict[str, int] = {}
    for r in records:
        counts[r.ss_method or "none"] = counts.get(r.ss_method or "none", 0) + 1
    total = len(records)
    parts = ", ".join(f"{m} {n} ({100.0 * n / total:.1f}%)" for m, n in sorted(counts.items()))
    text = f"Secondary structure method  : {parts}\n"
    if set(counts) <= {"BioPDB_dihedral", "none"}:
        text += ("  NOTE: DSSP and VMD/STRIDE were not used (not found or failed); only the less "
                 "accurate dihedral method was used. Run --check_ss for help.\n")
    elif counts.get("BioPDB_dihedral"):
        text += ("  NOTE: some motifs used the dihedral fallback (DSSP/STRIDE analyse only the first "
                 "model, or failed for that structure).\n")
    return text


def _model_line(model_policy: str | None) -> str:
    if model_policy == "all":
        return "Models used                  : all models of each structure (--first_model_only for the first only)\n"
    if model_policy == "first":
        return "Models used                  : first model of each structure only\n"
    return ""


def print_summary(
    accounting: dict[str, int],
    metal: list[str],
    any_metal: bool,
    motif_records: list[MotifRecord],
    motif_dir: str,
    requested_length: int | None = None,
    cn_check_text: str = "",
    model_policy: str | None = None,
) -> None:
    """Print a summary report in which every removal has its own label."""
    a = accounting
    metal_label = "any (all known metals)" if any_metal else "/".join(metal)

    motifs_by_metal: dict[str, int] = {}
    structures_by_metal: dict[str, set[str]] = {}
    for rec in motif_records:
        motifs_by_metal[rec.metal] = motifs_by_metal.get(rec.metal, 0) + 1
        structures_by_metal.setdefault(rec.metal, set()).add(rec.pdb_id)

    per_metal_lines = ""
    for m in sorted(motifs_by_metal):
        per_metal_lines += (
            f"  {m:<4} -> {len(structures_by_metal[m])} structures, "
            f"{motifs_by_metal[m]} motifs\n"
        )

    uniprot_block = ""
    if "uniprot_entries" in a:
        uniprot_block = (
            f"UniProt entries returned     : {a['uniprot_entries']}\n"
            f"  entries with PDB xrefs     : {a['entries_with_pdb']}\n"
            f"  raw PDB cross-references   : {a['raw_xrefs']}\n"
        )

    if a.get("metal_filter_applied"):
        metal_line = (
            f"Removed: requested metal absent : {a['removed_no_metal']}   "
            f"(metal filter ON: {metal_label})\n"
        )
        if a.get("removed_cluster_only"):
            metal_line += (
                f"  of these, only a metal cluster : {a['removed_cluster_only']}   "
                "(Fe-S clusters FES/SF4/F3S are multi-atom residues; NOT supported)\n"
            )
    else:
        metal_line = "Metal-presence filter        : OFF (--metal any); nothing removed\n"

    summary = (
        "\n=========================================\n"
        f"{uniprot_block}"
        f"Unique candidate PDB IDs     : {a['candidates']}\n"
        f"Removed: metadata unavailable: {a['metadata_failed']}\n"
        f"Removed: no resolution       : {a['removed_no_resolution']}\n"
        f"Removed: method not allowed  : {a['removed_method']}\n"
        f"Removed: resolution too low  : {a['removed_resolution']}\n"
        f"{metal_line}"
        f"Passed all filters           : {a['passed_filters']}\n"
        f"Failed downloads             : {a['download_failures']}\n"
        f"Downloaded structures        : {a['downloaded']}\n"
        f"  parse failures             : {a['parse_failure']}\n"
        f"  extraction errors (see log): {a['extraction_error']}\n"
        f"  no ion of requested metal  : {a['no_ion_of_requested_type']}\n"
        f"  only a metal cluster (n/s) : {a.get('metal_cluster_not_supported', 0)}\n"
        f"  ion but no protein ligand  : {a['ion_without_protein_ligand']}\n"
        f"  structures with motifs     : {a['sites_extracted']}\n"
        f"Metals requested             : {metal_label}\n"
        f"{_model_line(model_policy)}"
        f"{ss_method_summary(motif_records)}"
        f"{per_metal_lines}"
        f"Motifs extracted (total)     : {len(motif_records)}\n"
        f"Coordinated waters (total)   : {sum(r.water_count for r in motif_records)}\n"
        f"Non-protein ligand atoms     : {sum(r.ligand_atom_count for r in motif_records)}\n"
        f"{cn_check_text}"
        f"{format_length_summary(length_summary(motif_records, requested_length)) if requested_length else ''}"
        f"Output directory             : {motif_dir}/\n"
        "=========================================\n"
    )
    print(summary)
    log.info(summary)


def main() -> None:
    """Entry point: parse arguments, run the pipeline, print summary."""
    parser = build_parser()
    args = parser.parse_args()
    overrides = parse_args_to_overrides(args)
    if args.check_ss:
        import secondary_structure
        secondary_structure.configure_executables(args.dssp_exe, args.vmd_exe)
        res = secondary_structure.self_test()
        sys.exit(0 if (res["DSSP"]["ok"] or res["VMD_STRIDE"]["ok"]) else 1)
    try:
        cfg = load_config(args.config, overrides=overrides)
    except ConfigError as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(2)

    import secondary_structure
    secondary_structure.configure_executables(cfg.get("dssp_exe"), cfg.get("vmd_exe"))

    setup_logging(cfg["log_dir"])
    log.info("=" * 60)
    log.info("MetalMotifExtractor starting")
    log.info("Configuration: %s", cfg)
    log.info("Secondary structure method: %s", secondary_structure.available_methods())
    log.info("=" * 60)

    ensure_dirs(
        cfg["output_dir"], cfg["pdb_dir"], cfg["motif_dir"],
        cfg["metadata_dir"], cfg["log_dir"],
    )

    start_time = time.time()
    accounting: dict[str, int] = {}

    # ------------------------------------------------------------------
    # Step 1: Acquire candidate PDB IDs
    # ------------------------------------------------------------------
    candidate_ids, discovery = acquire_pdb_ids(cfg)
    if discovery is not None:
        accounting.update({
            k: discovery.stats[k]
            for k in ("uniprot_entries", "entries_with_pdb", "raw_xrefs")
        })
        accounting["candidates"] = discovery.stats["unique_pdb_ids"]
        accounting["metadata_failed"] = discovery.stats["metadata_failed"]
        accounting["removed_no_resolution"] = discovery.stats["removed_no_resolution"]
        accounting["removed_method"] = discovery.stats["removed_method"]
        accounting["removed_resolution"] = discovery.stats["removed_resolution"]
        if discovery.stats["unique_pdb_ids"] == 0:
            log.error("No PDB IDs were acquired. Exiting.")
            sys.exit(1)
    else:
        accounting["candidates"] = len(candidate_ids)
        if not candidate_ids:
            log.error("No PDB IDs were acquired. Exiting.")
            sys.exit(1)
    log.info("Acquired %d candidate PDB IDs", accounting["candidates"])

    # ------------------------------------------------------------------
    # Step 2: remaining filters (metadata is not fetched/filtered twice)
    # ------------------------------------------------------------------
    filtered_ids, meta_map, fstats = filter_and_fetch_metadata(
        candidate_ids, cfg, discovery
    )
    if discovery is None:
        for k in ("metadata_failed", "removed_no_resolution",
                  "removed_method", "removed_resolution"):
            accounting[k] = fstats[k]
    accounting["removed_no_metal"] = fstats["removed_no_metal"]
    accounting["removed_cluster_only"] = fstats["removed_cluster_only"]
    accounting["metal_filter_applied"] = fstats["metal_filter_applied"]
    accounting["passed_filters"] = fstats["passed_filters"]

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
        asu_format=cfg["asu_format"],
    )
    accounting["download_failures"] = len(filtered_ids) - len(downloaded)
    accounting["downloaded"] = len(downloaded)

    # ------------------------------------------------------------------
    # Step 4: Motif extraction
    # ------------------------------------------------------------------
    all_records: list[MotifRecord] = []
    outcomes: dict[str, str] = {}
    motif_counts: dict[str, int] = {}

    timings: dict[str, float] = {}
    warned_metals: set[str] = set()
    for pdb_id, pdb_file in downloaded.items():
        t_struct = time.time()
        try:
            records, status = extract_motifs_with_status(
                pdb_id=pdb_id,
                pdb_file=pdb_file,
                metal=cfg["metal"],
                distance_cutoff=cfg["distance_cutoff"],
                motif_length=cfg["motif_length"],
                motif_dir=cfg["motif_dir"],
                structure_meta=meta_map.get(pdb_id, {}),
                include_water=cfg["include_water"],
                water_cutoff=cfg["water_cutoff"],
                model_policy=cfg["model_policy"],
                include_ligands=cfg["include_ligands"],
                ligand_cutoff=cfg["ligand_cutoff"],
                geometry_max_rmsd=cfg["geometry_max_rmsd"],
                metal_cutoffs=cfg["metal_cutoffs"],
                max_cn=cfg["max_coordination_resolved"],
                findgeo_input_dir=(
                    cfg["findgeo_dir"] if cfg["write_findgeo_input"] else None
                ),
            )
        except Exception:
            # One problematic structure must not abort a multi-hour run.
            # Log the full traceback and record the structure as an error.
            log.exception("%s: extraction failed - structure skipped", pdb_id)
            records, status = [], "extraction_error"
        for rec in records:
            if rec.cn_exceeds_max and rec.metal not in warned_metals:
                warned_metals.add(rec.metal)
                log.warning(
                    "%s: %s site has %d donor atoms, above the highest coordination number "
                    "expected for %s (%d). The distance cut-off (%s A) may be too high for "
                    "this metal; see the summary for the share of sites affected.",
                    pdb_id, rec.metal, rec.coordination_number, rec.metal,
                    rec.expected_max_cn, cfg["distance_cutoff"])
        timings[pdb_id] = round(time.time() - t_struct, 2)
        if timings[pdb_id] > 60:
            log.warning("%s: extraction took %.0f s (%d motifs)", pdb_id,
                        timings[pdb_id], len(records))
        outcomes[pdb_id] = status
        motif_counts[pdb_id] = len(records)
        all_records.extend(records)

    for status in ("sites_extracted", "no_ion_of_requested_type",
                   "metal_cluster_not_supported",
                   "ion_without_protein_ligand", "parse_failure",
                   "extraction_error"):
        accounting[status] = sum(1 for v in outcomes.values() if v == status)

    # ------------------------------------------------------------------
    # Step 5: CSV export
    # ------------------------------------------------------------------
    write_structure_metadata_csv(
        meta_map, downloaded, cfg, motif_records=all_records, outcomes=outcomes
    )
    write_motif_metadata_csv(all_records, cfg)
    write_secondary_structure_csv(all_records, cfg)
    write_accounting_files(accounting, outcomes, motif_counts, discovery, cfg, timings)
    if timings:
        slow = sorted(timings.items(), key=lambda kv: -kv[1])[:3]
        log.info("Slowest structures: %s",
                 ", ".join(f"{pid} ({sec:.1f} s)" for pid, sec in slow))
        print("Slowest structures: " + ", ".join(f"{pid} ({sec:.1f} s)" for pid, sec in slow))

    # ------------------------------------------------------------------
    # Step 6: Summary
    # ------------------------------------------------------------------
    elapsed = time.time() - start_time
    log.info("Pipeline completed in %.1f seconds", elapsed)

    print_summary(
        accounting=accounting,
        metal=cfg["metal"],
        any_metal=cfg.get("any_metal", False),
        motif_records=all_records,
        motif_dir=cfg["motif_dir"],
        requested_length=cfg["motif_length"],
        cn_check_text=format_cn_check(cn_check(all_records, cfg["max_coordination_resolved"]), cfg),
        model_policy=cfg["model_policy"],
    )
    write_length_distribution(all_records, cfg)
    write_cn_check(all_records, cfg)


if __name__ == "__main__":
    main()
