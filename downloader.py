"""
downloader.py
=============
Handles all three PDB ID acquisition modes and parallel structure downloads.

Input Modes
-----------
A  **excel**     — read PDB IDs from a single-column ``.xlsx`` file.
B  **text**      — read PDB IDs from a plain text file, one per line.
C  **interpro**  — automatically discover EF-hand calcium-binding proteins
                   via the recommended API workflow described below.

Recommended API Workflow for Mode C (EF-hand discovery)
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
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import pandas as pd
import requests

from metadata import fetch_entry_metadata, passes_filter
from utils import ensure_dirs, http_get, safe_json

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# API endpoints
# ---------------------------------------------------------------------------
_UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
_RCSB_PDB_URL   = "https://files.rcsb.org/download/{pdb_id}.pdb"
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


def ids_from_interpro_efhand(
    max_resolution: float,
    allowed_methods: list[str],
    include_nmr: bool = False,
    page_size: int = 500,
    uniprot_query: str = '(ft_domain:"EF-hand") AND (reviewed:true)',
) -> list[str]:
    """Discover PDB structures for any protein family via UniProt REST API.

    This function is the backbone of Mode D (``--mode uniprot``).  It
    accepts any valid UniProt query string, making it applicable to any
    protein family that has curated domain annotations in Swiss-Prot —
    not just EF-hand proteins.

    Workflow
    --------
    1. Query UniProt with the supplied *uniprot_query*, paginating through
       all results via cursor-based navigation.
    2. Extract PDB cross-references from each UniProt entry.
    3. Filter by resolution and experimental method using the RCSB Data API.

    The default query ``(ft_domain:"EF-hand") AND (reviewed:true)`` is
    retained for backward compatibility but can be overridden via the
    ``--uniprot_query`` CLI flag or ``uniprot_query`` config key.

    Protein family examples and their queries
    -----------------------------------------
    EF-hand (Ca²⁺):
        ``(ft_domain:"EF-hand") AND (reviewed:true)``
    C2H2 zinc finger:
        ``(ft_domain:"Zinc finger, C2H2-type") AND (reviewed:true)``
    RING finger (E3 ubiquitin ligase, Zn²⁺):
        ``(ft_domain:"RING-type zinc finger") AND (reviewed:true)``
    LIM domain (double Zn²⁺):
        ``(ft_domain:"LIM zinc-binding") AND (reviewed:true)``
    PHD finger (chromatin reader, Zn²⁺):
        ``(ft_domain:"PHD-type zinc finger") AND (reviewed:true)``
    Rubredoxin (Fe²⁺, all-Cys tetrahedral):
        ``(ft_domain:"Rubredoxin-type") AND (reviewed:true)``
    2Fe-2S ferredoxin:
        ``(ft_domain:"2Fe-2S ferredoxin-type") AND (reviewed:true)``
    EGF-like calcium-binding (distinct from EF-hand):
        ``(ft_domain:"EGF-like calcium-binding") AND (reviewed:true)``
    FYVE zinc finger:
        ``(ft_domain:"Zinc finger, FYVE-type") AND (reviewed:true)``

    Parameters
    ----------
    max_resolution : float
        Maximum acceptable crystallographic / cryo-EM resolution (Å).
    allowed_methods : list[str]
        Whitelisted experimental methods.
    include_nmr : bool
        Whether to accept NMR structures.
    page_size : int
        Number of UniProt entries per page request.
    uniprot_query : str
        Valid UniProt REST API query string.  The default targets EF-hand
        proteins; change this to target any other protein family.

    Returns
    -------
    list[str]
        Unique PDB IDs passing all filters.
    """
    log.info(
        "UniProt discovery workflow: querying with '%s' "
        "(resolution ≤ %.1f Å) …", uniprot_query, max_resolution
    )

    # ------------------------------------------------------------------
    # Step 1 — collect UniProt entries matching the query
    # ------------------------------------------------------------------
    pdb_ids_raw: set[str] = set()
    cursor: str | None = None
    total_fetched = 0
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
        total_fetched += len(entries)
        page_num += 1

        # ------------------------------------------------------------------
        # Step 2 — extract PDB cross-references
        # ------------------------------------------------------------------
        for entry in entries:
            xrefs = entry.get("uniProtKBCrossReferences") or []
            for xref in xrefs:
                if xref.get("database") == "PDB":
                    pdb_ids_raw.add(xref["id"].upper())

        log.debug(
            "UniProt page %d: %d entries (running total %d entries, "
            "%d unique PDB IDs)",
            page_num, len(entries), total_fetched, len(pdb_ids_raw),
        )

        # Pagination via Link header
        link_header = resp.headers.get("Link", "")
        if 'rel="next"' in link_header:
            import re
            match = re.search(r'cursor=([^&>]+)', link_header)
            if match:
                cursor = match.group(1)
            else:
                log.warning(
                    "Link header present but cursor could not be parsed "
                    "(page %d) — stopping pagination early: %s",
                    page_num, link_header,
                )
                break
        else:
            break

        if page_num >= 200:
            log.warning(
                "Reached pagination safety limit (%d pages, %d entries) — "
                "stopping early.", page_num, total_fetched,
            )
            break

    log.info(
        "UniProt returned %d entries → %d candidate PDB IDs",
        total_fetched, len(pdb_ids_raw),
    )

    # ------------------------------------------------------------------
    # Step 3 — filter by resolution and method (RCSB Data API)
    # ------------------------------------------------------------------
    accepted: list[str] = []
    filtered_resolution = 0
    filtered_method = 0
    failed = 0

    for pdb_id in sorted(pdb_ids_raw):
        meta = fetch_entry_metadata(pdb_id)
        if not meta:
            failed += 1
            continue
        if passes_filter(meta, max_resolution, allowed_methods, include_nmr):
            accepted.append(pdb_id)
        else:
            if meta.get("resolution") and float(meta["resolution"]) > max_resolution:
                filtered_resolution += 1
            else:
                filtered_method += 1

    log.info(
        "Resolution filter removed %d; method filter removed %d; "
        "fetch failures: %d; accepted: %d",
        filtered_resolution, filtered_method, failed, len(accepted),
    )
    return accepted


# ===========================================================================
# PDB file download
# ===========================================================================

def _download_single(
    pdb_id: str,
    pdb_dir: Path,
    use_biological_assembly: bool,
) -> tuple[str, str | None]:
    """Download a single PDB structure (thread-safe).

    Parameters
    ----------
    pdb_id : str
        Upper-case four-character PDB ID.
    pdb_dir : Path
        Destination directory.
    use_biological_assembly : bool
        If ``True``, download assembly 1 in mmCIF format.

    Returns
    -------
    tuple[str, str | None]
        ``(pdb_id, local_path)`` — local_path is ``None`` on failure.
    """
    pdb_id = pdb_id.upper()

    if use_biological_assembly:
        filename = f"{pdb_id.lower()}-assembly1.cif"
        url = _RCSB_ASSEM_URL.format(pdb_id=pdb_id.lower(), n=1)
    else:
        filename = f"{pdb_id.lower()}.pdb"
        url = _RCSB_PDB_URL.format(pdb_id=pdb_id.lower())

    dest = pdb_dir / filename
    if dest.exists():
        log.debug("%s already cached at %s", pdb_id, dest)
        return pdb_id, str(dest)

    try:
        resp = http_get(url, timeout=60)
        dest.write_bytes(resp.content)
        log.info("Downloaded %s → %s", pdb_id, dest)
        return pdb_id, str(dest)
    except requests.RequestException as exc:
        log.error("Failed to download %s: %s", pdb_id, exc)
        return pdb_id, None


def download_structures(
    pdb_ids: list[str],
    pdb_dir: str | Path,
    use_biological_assembly: bool = False,
    max_workers: int = 4,
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
        Number of parallel download threads.  Kept low (default 4) to
        respect RCSB's fair-use rate limits.

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
            pool.submit(_download_single, pid, pdb_dir, use_biological_assembly): pid
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
