"""
downloader.py
=============
Handles all three PDB ID acquisition modes and parallel structure downloads.

Input Modes
-----------
A  **excel**     — read PDB IDs from a single-column ``.xlsx`` file.
B  **text**      — read PDB IDs from a plain text file, one per line.
C  **uniprot**   — automatically discover protein structures (EF-hand by
                   default) via the UniProt workflow described below.

Recommended API Workflow for Mode C (UniProt discovery)
--------------------------------------------------------
Rather than querying InterPro directly for every metal (InterPro's REST API
returns domain entries, not individual PDB structures, making it awkward to
enumerate all associated PDB IDs), we use the following three-step pipeline
that is more reliable and faster:

**Step 1 — UniProt REST API**
Query UniProt for reviewed (Swiss-Prot) proteins annotated with an
EF-hand structural domain. UniProt curates EF-hand annotations as a
*structured DOMAIN feature* (``FT DOMAIN ... /note="EF-hand"``), not as
a free-text MOTIF feature — so the query must target ``ft_domain``, not
``ft_motif``. (An earlier version of this tool incorrectly queried
``ft_motif:"EF-hand"``, which matches almost nothing because EF-hand is
essentially never annotated as a MOTIF feature in UniProt; correcting
this to ``ft_domain:"EF-hand"`` increases the candidate set roughly two
orders of magnitude.) The REST API supports boolean queries combined
with ``reviewed:true``, returning JSON with cross-references to PDB.

Reference: https://rest.uniprot.org/uniprotkb/search

**Step 2 — Extract PDB cross-references**
Each UniProt JSON entry contains ``uniProtKBCrossReferences`` listing PDB
IDs together with chain mappings.  We collect all unique PDB IDs from
these cross-references.

**Step 3 — Resolution / method filtering (delegated to metadata.py)**
Each candidate PDB ID is checked via the RCSB Data API
(``data.rcsb.org/rest/v1/core/entry/{pdb_id}``) for resolution and
experimental method.  Only structures satisfying the user's thresholds are
retained.

Why not InterPro REST API directly?
    InterPro entries (e.g. IPR002048 — EF-hand domain) list member proteins
    via UniProt accessions, but fetching PDB IDs from InterPro requires two
    extra API hops (InterPro → UniProt → PDB) with more complex pagination.
    Querying UniProt directly with ``ft_domain:"EF-hand"`` achieves the same
    result in one step with higher-quality annotation filtering.

Why not PDBe SIFTS API?
    SIFTS is excellent for residue-level UniProt↔PDB mapping but overkill
    for bulk PDB ID discovery.  We use it only if UniProt cross-references
    are insufficient.

Biological rationale
--------------------
The EF-hand is a helix–loop–helix motif that binds Ca²⁺ (and occasionally
Mg²⁺, Sr²⁺) through a canonical 12-residue loop.  EF-hand proteins include
calmodulin, troponin C, calbindin, recoverin and many others.  Restricting
discovery to the EF-hand domain annotation (``ft_domain:"EF-hand"`` in
UniProt — the structured PROSITE-derived DOMAIN feature, e.g.
``FT DOMAIN ... /note="EF-hand 1"``) ensures that only bona-fide EF-hand
Ca²⁺-binding structures are analysed, rather than all calcium-containing
structures in the PDB.
"""

from __future__ import annotations

import logging
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd
import requests

from metadata import fetch_metadata_parallel, filter_reason
from utils import ensure_dirs, http_get, safe_json

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# API endpoints
# ---------------------------------------------------------------------------
_UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
_RCSB_PDB_URL   = "https://files.rcsb.org/download/{pdb_id}.pdb"
_RCSB_CIF_URL   = "https://files.rcsb.org/download/{pdb_id}.cif"
_RCSB_ASSEM_URL = "https://files.rcsb.org/download/{pdb_id}-assembly{n}.cif"


# ===========================================================================
# PDB ID acquisition — three input modes
# ===========================================================================

def ids_from_excel(path: str | Path) -> list[str]:
    """Read PDB IDs from the first column of an Excel file.

    Parameters
    ----------
    path : str or Path
        Path to ``.xlsx`` file.  The first column must contain PDB IDs;
        a header row is optional and is auto-detected by pandas.

    Returns
    -------
    list[str]
        Unique, upper-case PDB IDs.
    """
    path = Path(path)
    df = pd.read_excel(path, header=0, dtype=str)
    col = df.iloc[:, 0].dropna()
    ids = [str(v).strip().upper() for v in col if str(v).strip()]
    ids = list(dict.fromkeys(ids))  # preserve order, deduplicate
    log.info("Read %d PDB IDs from Excel: %s", len(ids), path)
    return ids


def ids_from_text(path: str | Path) -> list[str]:
    """Read PDB IDs from a plain-text file (one ID per line).

    Parameters
    ----------
    path : str or Path
        Path to text file.

    Returns
    -------
    list[str]
        Unique, upper-case PDB IDs.
    """
    path = Path(path)
    lines = path.read_text(encoding="utf-8").splitlines()
    ids = [ln.strip().upper() for ln in lines if ln.strip() and not ln.startswith("#")]
    ids = list(dict.fromkeys(ids))
    log.info("Read %d PDB IDs from text file: %s", len(ids), path)
    return ids


@dataclass
class DiscoveryResult:
    """Everything produced by the UniProt discovery stage.

    Attributes
    ----------
    ids : list[str]
        PDB IDs that passed resolution / method filtering (sorted).
    meta_map : dict[str, dict]
        Metadata for **every** candidate whose lookup succeeded (accepted or
        not).  Kept so that later stages do not fetch the same entries again.
    stats : dict[str, int]
        Stage counts for the run summary (see keys below).
    accession_to_pdb : dict[str, list[str]]
        UniProt accession -> PDB IDs cross-referenced by that entry.  Written
        to ``uniprot_pdb_xrefs.csv`` so the entry -> structure accounting can
        be checked independently.

    Stats keys
    ----------
    ``uniprot_entries``       entries returned by the UniProt query
    ``entries_with_pdb``      of those, entries with >= 1 PDB cross-reference
    ``raw_xrefs``             PDB cross-references before de-duplication
    ``unique_pdb_ids``        unique PDB IDs after de-duplication
    ``metadata_failed``       IDs whose RCSB metadata lookup failed
    ``removed_no_resolution`` no numeric resolution (and not accepted NMR)
    ``removed_method``        experimental method not allowed
    ``removed_resolution``    resolution worse than the cutoff
    ``accepted``              IDs passing this stage
    """
    ids: list[str] = field(default_factory=list)
    meta_map: dict[str, dict] = field(default_factory=dict)
    stats: dict[str, int] = field(default_factory=dict)
    accession_to_pdb: dict[str, list[str]] = field(default_factory=dict)


def discover_from_uniprot(
    max_resolution: float,
    allowed_methods: list[str],
    include_nmr: bool = False,
    page_size: int = 500,
    uniprot_query: str = '(ft_domain:"EF-hand") AND (reviewed:true)',
    max_workers: int = 8,
) -> DiscoveryResult:
    """Discover PDB structures for any protein family via the UniProt REST API.

    Accepts any valid UniProt query string, so it applies to any protein
    family with curated Swiss-Prot domain annotations, not just EF-hand
    proteins.

    Workflow
    --------
    1. Query UniProt with *uniprot_query*, paginating through all results
       via cursor-based navigation (``Link`` response header).
    2. Extract **every** PDB cross-reference from each entry.  An entry can
       carry many PDB IDs, many entries carry none, and one PDB ID can be
       cross-referenced from several entries, so the number of unique PDB
       IDs is not bounded by the number of entries.  IDs are upper-cased
       and de-duplicated across entries.
    3. Fetch RCSB metadata once per unique ID (in parallel) and filter by
       resolution and experimental method.  Each removal is counted under
       its own reason.

    Protein family examples and their queries
    -----------------------------------------
    EF-hand (Ca2+):      ``(ft_domain:"EF-hand") AND (reviewed:true)``
    C2H2 zinc finger:    ``(ft_domain:"Zinc finger, C2H2-type") AND (reviewed:true)``
    RING finger (Zn2+):  ``(ft_domain:"RING-type zinc finger") AND (reviewed:true)``
    LIM domain (Zn2+):   ``(ft_domain:"LIM zinc-binding") AND (reviewed:true)``
    PHD finger (Zn2+):   ``(ft_domain:"PHD-type zinc finger") AND (reviewed:true)``
    Rubredoxin (Fe2+):   ``(ft_domain:"Rubredoxin-type") AND (reviewed:true)``
    2Fe-2S ferredoxin:   ``(ft_domain:"2Fe-2S ferredoxin-type") AND (reviewed:true)``
    EGF-like Ca-binding: ``(ft_domain:"EGF-like calcium-binding") AND (reviewed:true)``
    FYVE zinc finger:    ``(ft_domain:"Zinc finger, FYVE-type") AND (reviewed:true)``

    Parameters
    ----------
    max_resolution : float
        Maximum acceptable crystallographic / cryo-EM resolution (A).
    allowed_methods : list[str]
        Whitelisted experimental methods.
    include_nmr : bool
        Whether to accept NMR structures.
    page_size : int
        Number of UniProt entries per page request.
    uniprot_query : str
        Valid UniProt REST API query string.
    max_workers : int
        Threads for the metadata lookups.

    Returns
    -------
    DiscoveryResult
    """
    log.info(
        "UniProt discovery workflow: querying with '%s' "
        "(resolution <= %.1f A) ...", uniprot_query, max_resolution
    )

    # ------------------------------------------------------------------
    # Step 1 - collect UniProt entries matching the query
    # ------------------------------------------------------------------
    pdb_ids_raw: set[str] = set()
    accession_to_pdb: dict[str, list[str]] = {}
    cursor: str | None = None
    total_entries = 0
    entries_with_pdb = 0
    raw_xrefs = 0
    page_num = 0

    params: dict[str, Any] = {
        "query": uniprot_query,
        "fields": "xref_pdb",
        "format": "json",
        "size": page_size,
    }

    while True:
        page_params = dict(params)
        if cursor:
            page_params["cursor"] = cursor

        try:
            resp = http_get(_UNIPROT_SEARCH, params=page_params, timeout=60)
        except requests.RequestException as exc:
            log.error("UniProt query failed on page %d: %s", page_num + 1, exc)
            break

        data = safe_json(resp, context=f"UniProt page {page_num + 1}")
        if data is None:
            log.error(
                "UniProt returned an unparseable response on page %d. "
                "If you are behind an institutional proxy, ensure "
                "--no_ssl_verify is set and the UniProt domain is not "
                "blocked by your network firewall.",
                page_num + 1,
            )
            break

        entries = data.get("results") or []
        total_entries += len(entries)
        page_num += 1

        # --------------------------------------------------------------
        # Step 2 - extract ALL PDB cross-references of every entry
        # --------------------------------------------------------------
        for entry in entries:
            acc = entry.get("primaryAccession", "")
            xrefs = entry.get("uniProtKBCrossReferences") or []
            entry_pdbs: list[str] = []
            for xref in xrefs:
                if xref.get("database") == "PDB":
                    pid = xref["id"].upper()
                    raw_xrefs += 1
                    entry_pdbs.append(pid)
                    pdb_ids_raw.add(pid)
            if entry_pdbs:
                entries_with_pdb += 1
                if acc:
                    accession_to_pdb[acc] = sorted(set(entry_pdbs))

        log.debug(
            "UniProt page %d: %d entries (running total %d entries, "
            "%d unique PDB IDs)",
            page_num, len(entries), total_entries, len(pdb_ids_raw),
        )

        # Pagination via Link header
        link_header = resp.headers.get("Link", "")
        if 'rel="next"' in link_header:
            match = re.search(r'cursor=([^&>]+)', link_header)
            if match:
                cursor = match.group(1)
            else:
                log.warning(
                    "Link header present but cursor could not be parsed "
                    "(page %d) - stopping pagination early: %s",
                    page_num, link_header,
                )
                break
        else:
            break

        if page_num >= 200:
            log.warning(
                "Reached pagination safety limit (%d pages, %d entries) - "
                "stopping early.", page_num, total_entries,
            )
            break

    log.info(
        "UniProt returned %d entries (%d with PDB cross-references) -> "
        "%d raw cross-references -> %d unique PDB IDs",
        total_entries, entries_with_pdb, raw_xrefs, len(pdb_ids_raw),
    )

    # ------------------------------------------------------------------
    # Step 3 - metadata (once per ID, in parallel) and filtering by reason
    # ------------------------------------------------------------------
    candidates = sorted(pdb_ids_raw)
    meta_map = fetch_metadata_parallel(
        candidates, max_workers=max_workers, include_entities=False
    )

    accepted: list[str] = []
    counts = {"no_resolution": 0, "method": 0, "resolution": 0}
    for pdb_id in candidates:
        meta = meta_map.get(pdb_id)
        if not meta:
            continue
        reason = filter_reason(meta, max_resolution, allowed_methods, include_nmr)
        if reason == "":
            accepted.append(pdb_id)
        elif reason in counts:
            counts[reason] += 1

    metadata_failed = len(candidates) - len(meta_map)

    stats = {
        "uniprot_entries": total_entries,
        "entries_with_pdb": entries_with_pdb,
        "raw_xrefs": raw_xrefs,
        "unique_pdb_ids": len(candidates),
        "metadata_failed": metadata_failed,
        "removed_no_resolution": counts["no_resolution"],
        "removed_method": counts["method"],
        "removed_resolution": counts["resolution"],
        "accepted": len(accepted),
    }
    log.info(
        "Discovery filtering: %d unique IDs | %d metadata failures | "
        "%d removed (no resolution) | %d removed (method) | "
        "%d removed (resolution) | %d accepted",
        stats["unique_pdb_ids"], metadata_failed, counts["no_resolution"],
        counts["method"], counts["resolution"], len(accepted),
    )
    return DiscoveryResult(
        ids=accepted, meta_map=meta_map, stats=stats,
        accession_to_pdb=accession_to_pdb,
    )


def ids_from_uniprot(
    max_resolution: float,
    allowed_methods: list[str],
    include_nmr: bool = False,
    page_size: int = 500,
    uniprot_query: str = '(ft_domain:"EF-hand") AND (reviewed:true)',
    max_workers: int = 8,
) -> list[str]:
    """List-returning wrapper around :func:`discover_from_uniprot`."""
    return discover_from_uniprot(
        max_resolution, allowed_methods, include_nmr, page_size,
        uniprot_query, max_workers,
    ).ids


# Backward-compatible name (the function is not EF-hand or InterPro specific)
ids_from_interpro_efhand = ids_from_uniprot


# ===========================================================================
# PDB file download
# ===========================================================================

def _download_single(
    pdb_id: str,
    pdb_dir: Path,
    use_biological_assembly: bool,
    asu_format: str = "pdb",
) -> tuple[str, str | None]:
    """Download a single PDB structure (thread-safe).

    The asymmetric unit is fetched in PDB format; if that file does not
    exist (very large structures have no PDB-format file) the mmCIF file is
    fetched instead.  Biological assembly 1 is always mmCIF.  The file
    extension of the returned path (``.pdb`` or ``.cif``) tells the
    extraction code which parser to use.

    Parameters
    ----------
    pdb_id : str
        Upper-case four-character PDB ID.
    pdb_dir : Path
        Destination directory.
    use_biological_assembly : bool
        If ``True``, download assembly 1 in mmCIF format.
    asu_format : {"pdb", "cif"}
        Format of the asymmetric unit when the assembly is not requested.
        ``"pdb"`` (default) falls back to mmCIF only if no PDB file exists.
        ``"cif"`` always downloads mmCIF, which lets the same asymmetric unit
        be analysed in both formats (the two should give identical motifs).

    Returns
    -------
    tuple[str, str | None]
        ``(pdb_id, local_path)`` - local_path is ``None`` on failure.
    """
    pdb_id = pdb_id.upper()
    low = pdb_id.lower()

    if use_biological_assembly:
        options = [(f"{low}-assembly1.cif", _RCSB_ASSEM_URL.format(pdb_id=low, n=1))]
    elif asu_format == "cif":
        options = [(f"{low}.cif", _RCSB_CIF_URL.format(pdb_id=low))]
    else:
        options = [
            (f"{low}.pdb", _RCSB_PDB_URL.format(pdb_id=low)),
            (f"{low}.cif", _RCSB_CIF_URL.format(pdb_id=low)),   # fallback
        ]

    # Cached copy?
    for filename, _ in options:
        dest = pdb_dir / filename
        if dest.exists() and dest.stat().st_size > 0:
            log.debug("%s already cached at %s", pdb_id, dest)
            return pdb_id, str(dest)

    last_exc: Exception | None = None
    for filename, url in options:
        dest = pdb_dir / filename
        try:
            resp = http_get(url, timeout=60)
            if not resp.content:
                raise requests.RequestException("empty response body")
            dest.write_bytes(resp.content)
            log.info("Downloaded %s -> %s", pdb_id, dest)
            return pdb_id, str(dest)
        except requests.RequestException as exc:
            last_exc = exc
            log.debug("Download of %s failed (%s): %s", pdb_id, url, exc)

    log.error("Failed to download %s: %s", pdb_id, last_exc)
    return pdb_id, None


def download_structures(
    pdb_ids: list[str],
    pdb_dir: str | Path,
    use_biological_assembly: bool = False,
    max_workers: int = 4,
    asu_format: str = "pdb",
) -> dict[str, str]:
    """Download multiple PDB structures in parallel.

    Parameters
    ----------
    pdb_ids : list[str]
        PDB IDs to download.
    pdb_dir : str or Path
        Directory for downloaded files.
    use_biological_assembly : bool
        Download biological assembly 1 instead of asymmetric unit.
    max_workers : int
        Number of parallel download threads.  Keep modest to respect
        RCSB's fair-use rate limits.

    Returns
    -------
    dict[str, str]
        Mapping of ``pdb_id → local_file_path`` for successful downloads.
    """
    pdb_dir = Path(pdb_dir)
    ensure_dirs(pdb_dir)

    results: dict[str, str] = {}
    failed: list[str] = []

    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = {
            pool.submit(_download_single, pid, pdb_dir, use_biological_assembly, asu_format): pid
            for pid in pdb_ids
        }
        for future in as_completed(futures):
            pdb_id, path = future.result()
            if path:
                results[pdb_id] = path
            else:
                failed.append(pdb_id)

    if failed:
        log.warning("Failed downloads (%d): %s", len(failed), ", ".join(failed))

    return results
