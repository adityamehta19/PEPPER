"""
metadata.py
===========
Fetch and filter PDB structure metadata from the RCSB REST API.

Why RCSB Search + Data APIs?
-----------------------------
The RCSB provides two complementary REST endpoints:

1. **RCSB Search API** (search.rcsb.org/rcsbsearch/v2/query)
   - Full-featured JSON query DSL
   - Supports filtering by resolution, method, ligand ID, etc.
   - Returns paginated PDB ID lists — ideal for bulk retrieval

2. **RCSB Data API** (data.rcsb.org/rest/v1/core/entry/{pdb_id})
   - Returns per-structure metadata in JSON
   - Provides resolution, method, release date, organism, chains, etc.

This two-step approach (search → metadata enrichment) is more reliable and
faster than scraping PDBe or parsing raw PDB headers, and avoids the
complexity of SPARQL queries required by UniProt's RDF endpoint.

For the InterPro workflow the caller (``downloader.py``) pre-filters IDs
by querying UniProt and then delegates final resolution/method filtering to
the functions in this module.
"""

from __future__ import annotations

import logging
from typing import Any

from utils import http_get, safe_json

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# RCSB Data API base URL
# ---------------------------------------------------------------------------
_RCSB_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
_RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"


# ---------------------------------------------------------------------------
# Per-entry metadata
# ---------------------------------------------------------------------------

def fetch_entry_metadata(pdb_id: str) -> dict[str, Any]:
    """Retrieve metadata for a single PDB entry from the RCSB Data API.

    Parameters
    ----------
    pdb_id : str
        Four-character PDB identifier (case-insensitive).

    Returns
    -------
    dict
        Dictionary with keys:

        ``pdb_id``
            Normalised upper-case PDB ID.
        ``resolution``
            Crystallographic or cryo-EM resolution in Å (``None`` for NMR).
        ``method``
            Experimental method string, e.g. ``"X-RAY DIFFRACTION"``.
        ``release_date``
            ISO-formatted release date string.
        ``organism``
            Source organism scientific name (first entry).
        ``protein_name``
            Macromolecule description (first polymer entity).
        ``chain_ids``
            Comma-separated asymmetric-unit chain IDs.

    Notes
    -----
    Returns an empty dict on network or parse errors; callers should treat
    an empty dict as a failed lookup.
    """
    pdb_id = pdb_id.upper()
    url = _RCSB_ENTRY_URL.format(pdb_id=pdb_id)

    try:
        resp = http_get(url)
        data = safe_json(resp, context=f"entry metadata {pdb_id}")
        if data is None:
            return {}
    except Exception as exc:
        log.warning("Could not fetch metadata for %s: %s", pdb_id, exc)
        return {}

    # ---- resolution --------------------------------------------------------
    resolution: float | None = None
    refine = data.get("refine") or []
    if refine and isinstance(refine, list):
        resolution = refine[0].get("ls_d_res_high")  # highest-resolution shell

    # Fall back to EM resolution field
    if resolution is None:
        em_data = data.get("em_3d_reconstruction") or []
        if em_data and isinstance(em_data, list):
            resolution = em_data[0].get("resolution")

    # ---- experimental method -----------------------------------------------
    method = "UNKNOWN"
    exp_methods = data.get("exptl") or []
    if exp_methods and isinstance(exp_methods, list):
        method = exp_methods[0].get("method", "UNKNOWN").upper()

    # ---- release date ------------------------------------------------------
    release_date = (data.get("rcsb_accession_info") or {}).get(
        "initial_release_date", ""
    )

    # ---- organism ----------------------------------------------------------
    organism = ""
    entity_src = data.get("rcsb_entity_source_organism") or []
    if entity_src and isinstance(entity_src, list):
        organism = entity_src[0].get("scientific_name", "")

    # ---- protein name ------------------------------------------------------
    protein_name = ""
    polymer_entities = data.get("polymer_entities") or []
    if not polymer_entities:
        # Try the flat description field that appears in some API versions
        struct = data.get("struct") or {}
        protein_name = struct.get("title", "")
    else:
        desc = polymer_entities[0].get("rcsb_polymer_entity", {})
        protein_name = desc.get("pdbx_description", "")

    # ---- chain IDs ---------------------------------------------------------
    chains: list[str] = []
    asym = data.get("rcsb_entry_container_identifiers", {})
    chains = asym.get("polymer_entity_ids_polymer", []) or []
    # Use auth_asym_ids when available (human-readable chain letters)
    assembly_chains: list[str] = []
    for poly in data.get("polymer_entities") or []:
        for inst in poly.get("rcsb_polymer_entity_instances") or []:
            ch = inst.get("rcsb_polymer_entity_instance_container_identifiers", {}).get(
                "auth_asym_id"
            )
            if ch:
                assembly_chains.append(ch)
    if assembly_chains:
        chains = sorted(set(assembly_chains))

    return {
        "pdb_id": pdb_id,
        "resolution": resolution,
        "method": method,
        "release_date": release_date,
        "organism": organism,
        "protein_name": protein_name,
        "chain_ids": ",".join(chains),
    }


# ---------------------------------------------------------------------------
# Metal presence check
# ---------------------------------------------------------------------------

_RCSB_NONPOLY_ENTITY_URL = (
    "https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{pdb_id}/{entity_id}"
)


def fetch_ligand_ids_in_entry(pdb_id: str) -> set[str]:
    """Return the set of small-molecule / ion component IDs present in a PDB entry.

    Uses the correct two-step RCSB Data API hierarchy:

    Step 1: ``/core/entry/{pdb_id}`` → read
    ``rcsb_entry_container_identifiers.non_polymer_entity_ids``
    to get the list of non-polymer entity integer IDs.

    Step 2: For each entity ID, ``/core/nonpolymer_entity/{pdb_id}/{entity_id}``
    → read ``rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id``
    to get the chemical component ID (e.g. ``"CA"``, ``"MG"``, ``"SO4"``).

    The critical correction from a previous implementation: ``nonpolymer_entities``
    is NOT a field inside the ``/core/entry/`` response. It is a separate endpoint
    in the RCSB hierarchy, identified by entry_id + entity_id. Attempting to read
    ``data["nonpolymer_entities"]`` from the entry response always returns None,
    causing all structures to appear as having no ligands.

    Parameters
    ----------
    pdb_id : str
        Four-character PDB ID (case-insensitive).

    Returns
    -------
    set[str]
        Upper-case PDB Chemical Component IDs (e.g. ``{"MG", "CA", "SO4"}``).
        Empty set if the entry has no non-polymer entities or on API error.
    """
    pdb_id = pdb_id.upper()

    # Step 1: get list of non-polymer entity IDs from the entry endpoint
    entry_url = _RCSB_ENTRY_URL.format(pdb_id=pdb_id)
    try:
        resp = http_get(entry_url)
        entry_data = safe_json(resp, context=f"nonpolymer entity list {pdb_id}")
        if entry_data is None:
            return set()
    except Exception as exc:
        log.warning("Could not fetch entry data for %s: %s", pdb_id, exc)
        return set()

    ids_container = entry_data.get("rcsb_entry_container_identifiers") or {}
    entity_ids = ids_container.get("non_polymer_entity_ids") or []

    if not entity_ids:
        log.debug("%s: no non-polymer entities in entry", pdb_id)
        return set()

    # Step 2: fetch the chemical component ID for each non-polymer entity
    ligand_ids: set[str] = set()
    for eid in entity_ids:
        url = _RCSB_NONPOLY_ENTITY_URL.format(pdb_id=pdb_id, entity_id=eid)
        try:
            resp = http_get(url)
            entity_data = safe_json(resp, context=f"{pdb_id} entity {eid}")
            if entity_data is None:
                continue
            comp_id = (
                entity_data
                .get("rcsb_nonpolymer_entity_container_identifiers", {})
                .get("nonpolymer_comp_id", "")
            )
            if comp_id:
                ligand_ids.add(comp_id.strip().upper())
        except Exception as exc:
            log.debug(
                "%s entity %s: could not fetch nonpolymer details: %s",
                pdb_id, eid, exc,
            )

    log.debug("%s ligands present: %s", pdb_id, sorted(ligand_ids))
    return ligand_ids


def has_any_metal(pdb_id: str, metals: list[str]) -> bool:
    """Return True if the PDB entry contains at least one of the requested metals.

    Parameters
    ----------
    pdb_id : str
        Four-character PDB ID.
    metals : list[str]
        Upper-case PDB Component IDs to look for (e.g. ``["CA", "MG"]``).

    Returns
    -------
    bool
        ``True`` if any of *metals* is present as a non-polymer entity.
        Returns ``True`` (permissive) on API failure so network errors do not
        silently discard potentially valid structures.
    """
    metals_set = {m.upper() for m in metals}
    present = fetch_ligand_ids_in_entry(pdb_id)

    if not present:
        log.warning(
            "%s: could not verify metal presence via API — will attempt extraction",
            pdb_id,
        )
        return True

    found = metals_set & present
    if found:
        log.debug("%s: confirmed metal(s) present: %s", pdb_id, sorted(found))
        return True

    log.info(
        "%s: skipped — none of %s found in structure (ligands present: %s)",
        pdb_id, sorted(metals_set), sorted(present),
    )
    return False


def rcsb_batch_filter(
    pdb_ids: list[str],
    metals: list[str],
    max_resolution: float,
    allowed_methods: list[str],
    include_nmr: bool = False,
) -> list[str]:
    """Filter a list of PDB IDs using a single RCSB Search API batch query.

    Replaces the previous N × (N+1) serial API calls (one entry fetch +
    one entity fetch per non-polymer entity, per structure) with a single
    RCSB Search API POST that applies resolution, method, and metal-presence
    filters simultaneously server-side.

    The query sends the candidate PDB IDs as an entry_id set filter combined
    with resolution, method, and chemical component (ligand) filters.  RCSB
    executes the intersection server-side and returns only the IDs that pass
    all conditions — no per-structure API calls required.

    Parameters
    ----------
    pdb_ids : list[str]
        Candidate PDB IDs to filter (from UniProt cross-references).
    metals : list[str]
        Metal component IDs to require (e.g. ``["CA", "MG"]``).
    max_resolution : float
        Maximum crystallographic resolution (Å).
    allowed_methods : list[str]
        Allowed experimental methods.
    include_nmr : bool
        If True, also accept NMR structures.

    Returns
    -------
    list[str]
        PDB IDs passing all filters, upper-cased.
    """
    if not pdb_ids:
        return []

    # Build one chemical-component node per requested metal, combined with OR
    metal_nodes = [
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id",
                "operator": "exact_match",
                "value": m.upper(),
            },
        }
        for m in metals
    ]
    metal_group = (
        metal_nodes[0]
        if len(metal_nodes) == 1
        else {
            "type": "group",
            "logical_operator": "or",
            "nodes": metal_nodes,
        }
    )

    # Method filter — OR across allowed methods
    method_nodes = [
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "exptl.method",
                "operator": "exact_match",
                "value": m.upper(),
            },
        }
        for m in allowed_methods
    ]
    if include_nmr:
        method_nodes.append({
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "exptl.method",
                "operator": "contains_phrase",
                "value": "NMR",
            },
        })
    method_group = (
        method_nodes[0]
        if len(method_nodes) == 1
        else {
            "type": "group",
            "logical_operator": "or",
            "nodes": method_nodes,
        }
    )

    # Entry ID set filter — restricts to our UniProt-derived candidate set
    entry_id_node = {
        "type": "terminal",
        "service": "text",
        "parameters": {
            "attribute": "rcsb_entry_container_identifiers.entry_id",
            "operator": "in",
            "value": [pid.upper() for pid in pdb_ids],
        },
    }

    # Resolution filter (only for non-NMR entries)
    resolution_node = {
        "type": "terminal",
        "service": "text",
        "parameters": {
            "attribute": "rcsb_entry_info.resolution_combined",
            "operator": "less_or_equal",
            "value": max_resolution,
        },
    }

    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                entry_id_node,
                metal_group,
                method_group,
                resolution_node,
            ],
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": 10000},
            "results_content_type": ["experimental"],
        },
    }

    try:
        resp = http_get(_RCSB_SEARCH_URL, json=query, timeout=60)
        data = safe_json(resp, context="RCSB batch filter")
        results = (data.get("result_set") or []) if data else []
        ids = [r["identifier"].upper() for r in results]
        log.info(
            "RCSB batch filter: %d candidates → %d passed "
            "(resolution ≤ %.1f Å, metal=%s)",
            len(pdb_ids), len(ids), max_resolution, metals,
        )
        return ids
    except Exception as exc:
        log.error(
            "RCSB batch filter failed (%s) — falling back to serial filtering",
            exc,
        )
        return []


def fetch_metadata_parallel(
    pdb_ids: list[str],
    max_workers: int = 8,
) -> dict[str, dict]:
    """Fetch per-structure metadata for multiple PDB IDs in parallel.

    Uses a ThreadPoolExecutor to issue concurrent HTTP requests to the
    RCSB Data API, replacing the previous sequential loop.  For 460
    structures, parallel fetching at 8 workers reduces wall time from
    ~5 minutes to ~30–60 seconds depending on network conditions.

    Parameters
    ----------
    pdb_ids : list[str]
        PDB IDs to fetch metadata for.
    max_workers : int
        Number of concurrent threads.  8 is a reasonable default that
        respects RCSB's fair-use limits while providing substantial
        speedup.

    Returns
    -------
    dict[str, str -> dict]
        Mapping of pdb_id → metadata dict.  Failed fetches are omitted.
    """
    from concurrent.futures import ThreadPoolExecutor, as_completed

    results: dict[str, dict] = {}

    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        future_to_id = {
            pool.submit(fetch_entry_metadata, pid): pid
            for pid in pdb_ids
        }
        for future in as_completed(future_to_id):
            pid = future_to_id[future]
            try:
                meta = future.result()
                if meta:
                    results[pid] = meta
            except Exception as exc:
                log.warning("Metadata fetch failed for %s: %s", pid, exc)

    log.info(
        "Parallel metadata fetch: %d/%d succeeded",
        len(results), len(pdb_ids),
    )
    return results


def passes_filter(
    meta: dict[str, Any],
    max_resolution: float,
    allowed_methods: list[str],
    include_nmr: bool = False,
) -> bool:
    """Decide whether a structure passes the resolution and method filters.

    Parameters
    ----------
    meta : dict
        Metadata dict as returned by :func:`fetch_entry_metadata`.
    max_resolution : float
        Upper bound on resolution (Å).  Smaller is better.
    allowed_methods : list[str]
        Whitelisted experimental methods (upper-case strings).
    include_nmr : bool
        If ``True``, NMR structures (no numeric resolution) are accepted
        regardless of *max_resolution*.

    Returns
    -------
    bool
        ``True`` if the structure should be processed.
    """
    if not meta:
        return False

    method = (meta.get("method") or "").upper()
    resolution = meta.get("resolution")

    # --- NMR handling -------------------------------------------------------
    if resolution is None:
        if include_nmr and "NMR" in method:
            log.debug("%s accepted (NMR, include_nmr=True)", meta.get("pdb_id"))
            return True
        log.debug("%s skipped (no resolution, NMR not allowed)", meta.get("pdb_id"))
        return False

    # --- Method filter ------------------------------------------------------
    if allowed_methods and method not in [m.upper() for m in allowed_methods]:
        log.debug(
            "%s skipped (method=%s not in allowed list)", meta.get("pdb_id"), method
        )
        return False

    # --- Resolution filter --------------------------------------------------
    if float(resolution) > max_resolution:
        log.debug(
            "%s skipped (resolution=%.2f > max=%.2f)",
            meta.get("pdb_id"), resolution, max_resolution,
        )
        return False

    return True


# ---------------------------------------------------------------------------
# RCSB Search API — direct structure search by metal ligand
# ---------------------------------------------------------------------------

def search_rcsb_by_metal(
    metals: list[str],
    max_resolution: float,
    allowed_methods: list[str],
    include_nmr: bool = False,
) -> list[str]:
    """Query the RCSB Search API for structures containing given metal ion(s).

    This is the correct discovery method for metals that are not specific
    to the EF-hand fold (e.g. ZN, FE, MN, CU).  Unlike Mode C which
    queries UniProt for EF-hand domain annotations, this mode queries the
    RCSB structural database directly by ligand chemical component ID,
    resolution, and experimental method.  It requires only RCSB network
    access and does not contact UniProt at all, making it the reliable
    choice when UniProt is blocked by an institutional firewall.

    Parameters
    ----------
    metals : list[str]
        Metal residue names as they appear in PDB HETATM records.
    max_resolution : float
        Maximum resolution (Å) to include.
    allowed_methods : list[str]
        Allowed experimental methods.
    include_nmr : bool
        If True, NMR structures are also accepted.

    Returns
    -------
    list[str]
        Deduplicated list of PDB IDs satisfying all criteria.
    """
    all_ids: set[str] = set()

    for metal in metals:
        # Build method nodes
        method_nodes = [
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "exptl.method",
                    "operator": "exact_match",
                    "value": m.upper(),
                },
            }
            for m in allowed_methods
        ]
        if include_nmr:
            method_nodes.append({
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "exptl.method",
                    "operator": "contains_phrase",
                    "value": "NMR",
                },
            })
        method_group = (
            method_nodes[0]
            if len(method_nodes) == 1
            else {
                "type": "group",
                "logical_operator": "or",
                "nodes": method_nodes,
            }
        )

        query = {
            "query": {
                "type": "group",
                "logical_operator": "and",
                "nodes": [
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": (
                                "rcsb_nonpolymer_entity_container_identifiers"
                                ".nonpolymer_comp_id"
                            ),
                            "operator": "exact_match",
                            "value": metal.upper(),
                        },
                    },
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_entry_info.resolution_combined",
                            "operator": "less_or_equal",
                            "value": max_resolution,
                        },
                    },
                    method_group,
                ],
            },
            "return_type": "entry",
            "request_options": {"paginate": {"start": 0, "rows": 10000}},
        }

        try:
            resp = http_get(_RCSB_SEARCH_URL, json=query, timeout=60)
            data = safe_json(resp, context=f"RCSB metal search {metal}")
            results = (data.get("result_set") or []) if data else []
            ids = [r["identifier"].upper() for r in results]
            log.info(
                "RCSB search for metal=%s, resolution≤%.1f Å → %d structures",
                metal, max_resolution, len(ids),
            )
            all_ids.update(ids)
        except Exception as exc:
            log.error("RCSB search failed for metal=%s: %s", metal, exc)

    return sorted(all_ids)
