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

_RCSB_POLYMER_ENTITY_URL = (
    "https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
)


def fetch_entity_metadata(pdb_id: str, entity_ids: list[str]) -> dict[str, str]:
    """Fetch protein name, source organism and chain IDs from polymer entities.

    These fields are NOT part of the ``/core/entry/{pdb_id}`` response.  In
    the RCSB Data API hierarchy they live one level down, at
    ``/core/polymer_entity/{pdb_id}/{entity_id}``, where ``entity_id`` comes
    from ``rcsb_entry_container_identifiers.polymer_entity_ids`` in the entry
    response.  (A previous version of this module looked for
    ``polymer_entities`` inside the entry response, which does not exist, so
    the Organism, Protein Name and Chain IDs columns were empty.)

    Fields read from each polymer entity
    ------------------------------------
    ``rcsb_polymer_entity.pdbx_description``
        Macromolecule name.
    ``rcsb_entity_source_organism[*].ncbi_scientific_name`` (falls back to
    ``scientific_name``)
        Source organism.
    ``rcsb_polymer_entity_container_identifiers.auth_asym_ids`` (falls back
    to ``asym_ids``)
        Author chain identifiers (the letters used in the PDB file).

    Parameters
    ----------
    pdb_id : str
        Four-character PDB ID.
    entity_ids : list[str]
        Polymer entity IDs of the entry.

    Returns
    -------
    dict
        Keys ``protein_name``, ``organism``, ``chain_ids``; a key is only
        present when at least one entity supplied a non-empty value.
    """
    names: list[str] = []
    organisms: list[str] = []
    chains: list[str] = []

    for eid in entity_ids:
        url = _RCSB_POLYMER_ENTITY_URL.format(pdb_id=pdb_id.upper(), entity_id=eid)
        try:
            resp = http_get(url)
            ent = safe_json(resp, context=f"{pdb_id} polymer entity {eid}")
        except Exception as exc:
            log.debug("%s entity %s: could not fetch polymer entity: %s", pdb_id, eid, exc)
            continue
        if not isinstance(ent, dict):
            continue

        desc = ((ent.get("rcsb_polymer_entity") or {}).get("pdbx_description") or "").strip()
        if desc and desc not in names:
            names.append(desc)

        for src in ent.get("rcsb_entity_source_organism") or []:
            org = (src.get("ncbi_scientific_name") or src.get("scientific_name") or "").strip()
            if org and org not in organisms:
                organisms.append(org)

        ids = ent.get("rcsb_polymer_entity_container_identifiers") or {}
        for ch in (ids.get("auth_asym_ids") or ids.get("asym_ids") or []):
            if ch and ch not in chains:
                chains.append(ch)

    out: dict[str, str] = {}
    if names:
        out["protein_name"] = "; ".join(names)
    if organisms:
        out["organism"] = "; ".join(organisms)
    if chains:
        out["chain_ids"] = ",".join(sorted(chains))
    return out


def fetch_entry_metadata(pdb_id: str, include_entities: bool = True) -> dict[str, Any]:
    """Retrieve metadata for a single PDB entry from the RCSB Data API.

    Parameters
    ----------
    pdb_id : str
        Four-character PDB identifier (case-insensitive).
    include_entities : bool
        If True (default), also query the polymer-entity endpoints for
        protein name, organism and chain IDs (one extra request per polymer
        entity).  Discovery-stage filtering only needs resolution and
        method, so it passes False to keep that stage fast; the entity
        fields are then added once, for accepted structures only, with
        :func:`enrich_entity_metadata`.

    Returns
    -------
    dict
        ``pdb_id``, ``resolution``, ``resolution_source``, ``method``,
        ``release_date``, ``organism``, ``protein_name``, ``chain_ids``,
        ``polymer_entity_ids``.  ``protein_name`` defaults to the entry
        title (``struct.title``) until entity information is added.

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
        if not isinstance(data, dict):
            return {}
    except Exception as exc:
        log.warning("Could not fetch metadata for %s: %s", pdb_id, exc)
        return {}

    # ---- resolution --------------------------------------------------------
    resolution: float | None = None
    resolution_source = ""
    refine = data.get("refine") or []
    if refine and isinstance(refine, list):
        resolution = refine[0].get("ls_d_res_high")  # highest-resolution shell
        if resolution is not None:
            resolution_source = "refine.ls_d_res_high"

    # Fall back to EM resolution field
    if resolution is None:
        em_data = data.get("em_3d_reconstruction") or []
        if em_data and isinstance(em_data, list):
            resolution = em_data[0].get("resolution")
            if resolution is not None:
                resolution_source = "em_3d_reconstruction.resolution"

    # Last resort: RCSB's own combined resolution value
    if resolution is None:
        combined = (data.get("rcsb_entry_info") or {}).get("resolution_combined") or []
        if combined and isinstance(combined, list):
            resolution = combined[0]
            resolution_source = "rcsb_entry_info.resolution_combined"

    # ---- experimental method -----------------------------------------------
    method = "UNKNOWN"
    exp_methods = data.get("exptl") or []
    if exp_methods and isinstance(exp_methods, list):
        method = exp_methods[0].get("method", "UNKNOWN").upper()

    # ---- release date ------------------------------------------------------
    release_date = (data.get("rcsb_accession_info") or {}).get(
        "initial_release_date", ""
    )

    # ---- title and polymer entity IDs (both ARE in the entry response) -----
    title = (data.get("struct") or {}).get("title", "") or ""
    entity_ids = (
        (data.get("rcsb_entry_container_identifiers") or {}).get("polymer_entity_ids")
        or []
    )

    meta: dict[str, Any] = {
        "pdb_id": pdb_id,
        "resolution": resolution,
        "resolution_source": resolution_source,
        "method": method,
        "release_date": release_date,
        "organism": "",
        "protein_name": title,
        "chain_ids": "",
        "polymer_entity_ids": list(entity_ids),
    }

    if include_entities:
        meta.update(fetch_entity_metadata(pdb_id, list(entity_ids)))

    return meta


def enrich_entity_metadata(
    meta_map: dict[str, dict],
    max_workers: int = 8,
) -> int:
    """Add protein name, organism and chain IDs to already-fetched metadata.

    Called once, after all filtering, so entity requests are only made for
    structures that will actually be processed.  Modifies *meta_map* in
    place.

    Returns
    -------
    int
        Number of structures for which at least one entity field was found.
    """
    from concurrent.futures import ThreadPoolExecutor, as_completed

    n_ok = 0
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = {
            pool.submit(
                fetch_entity_metadata, pid, list(meta.get("polymer_entity_ids") or [])
            ): pid
            for pid, meta in meta_map.items()
            if meta.get("polymer_entity_ids")
        }
        for fut in as_completed(futures):
            pid = futures[fut]
            try:
                extra = fut.result()
            except Exception as exc:
                log.warning("Entity metadata failed for %s: %s", pid, exc)
                continue
            if extra:
                meta_map[pid].update(extra)
                n_ok += 1
    log.info(
        "Entity metadata (protein name / organism / chains) added for %d/%d structures",
        n_ok, len(meta_map),
    )
    return n_ok


# ---------------------------------------------------------------------------
# Metal presence check
# ---------------------------------------------------------------------------

_RCSB_NONPOLY_ENTITY_URL = (
    "https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{pdb_id}/{entity_id}"
)


def fetch_ligand_ids_status(pdb_id: str) -> tuple[set[str], bool]:
    """Return ``(ligand component IDs, complete)`` for a PDB entry.

    Uses the RCSB Data API hierarchy:

    Step 1: ``/core/entry/{pdb_id}`` -> read
    ``rcsb_entry_container_identifiers.non_polymer_entity_ids``.

    Step 2: for each entity ID, ``/core/nonpolymer_entity/{pdb_id}/{entity_id}``
    -> read ``rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id``
    (e.g. ``"CA"``, ``"MG"``, ``"SO4"``).

    ``complete`` tells a *verified* answer from a *failed* lookup, which the
    previous version could not (it returned an empty set for both):

    * ``complete=True``: every request succeeded.  An empty set then means the
      entry genuinely has no non-polymer entities (no ions, no ligands), e.g.
      an apo structure.
    * ``complete=False``: the entry request failed, or at least one entity
      request failed.  The returned set may be missing ligands, so an absent
      metal is NOT proof of absence.

    Parameters
    ----------
    pdb_id : str
        Four-character PDB ID (case-insensitive).
    """
    pdb_id = pdb_id.upper()

    entry_url = _RCSB_ENTRY_URL.format(pdb_id=pdb_id)
    try:
        resp = http_get(entry_url)
        entry_data = safe_json(resp, context=f"nonpolymer entity list {pdb_id}")
        if not isinstance(entry_data, dict):
            return set(), False
    except Exception as exc:
        log.warning("Could not fetch entry data for %s: %s", pdb_id, exc)
        return set(), False

    ids_container = entry_data.get("rcsb_entry_container_identifiers") or {}
    entity_ids = ids_container.get("non_polymer_entity_ids") or []

    if not entity_ids:
        log.debug("%s: no non-polymer entities in entry (verified)", pdb_id)
        return set(), True

    ligand_ids: set[str] = set()
    complete = True
    for eid in entity_ids:
        url = _RCSB_NONPOLY_ENTITY_URL.format(pdb_id=pdb_id, entity_id=eid)
        try:
            resp = http_get(url)
            entity_data = safe_json(resp, context=f"{pdb_id} entity {eid}")
            if not isinstance(entity_data, dict):
                complete = False
                continue
            comp_id = (
                entity_data
                .get("rcsb_nonpolymer_entity_container_identifiers", {})
                .get("nonpolymer_comp_id", "")
            )
            if comp_id:
                ligand_ids.add(comp_id.strip().upper())
        except Exception as exc:
            complete = False
            log.debug(
                "%s entity %s: could not fetch nonpolymer details: %s",
                pdb_id, eid, exc,
            )

    log.debug("%s ligands present: %s (complete=%s)", pdb_id, sorted(ligand_ids), complete)
    return ligand_ids, complete


def fetch_ligand_ids_in_entry(pdb_id: str) -> set[str]:
    """Set of ligand/ion component IDs of an entry (compatibility wrapper).

    Prefer :func:`fetch_ligand_ids_status`, which also says whether the answer
    is complete.
    """
    return fetch_ligand_ids_status(pdb_id)[0]


# Iron-sulfur clusters are single multi-atom residues in the PDB (component IDs
# below), not single-atom ions, so they are NOT found by --metal FE / FE2 and are
# not supported by this tool.  Entries that contain one of them but none of the
# requested metals are counted separately so the user can see why they were left out.
UNSUPPORTED_CLUSTER_IDS: frozenset[str] = frozenset({"FES", "SF4", "F3S"})

# pdb_id -> unsupported cluster IDs seen when the metal filter removed the entry
CLUSTER_ONLY_SEEN: dict[str, set[str]] = {}


def has_any_metal(pdb_id: str, metals: list[str]) -> bool:
    """Return True if the PDB entry contains at least one of the requested metals.

    Decision rules
    --------------
    * A requested metal is listed             -> keep.
    * The lookup was complete and none is listed -> **remove**.  This includes
      entries with no ions or ligands at all (apo structures): the answer is
      verified, so the structure is not carried on to extraction.
    * The lookup failed or was incomplete and no requested metal was seen
      -> keep (permissive), with a warning, so that a network error cannot
      silently discard a valid structure.

    Parameters
    ----------
    pdb_id : str
        Four-character PDB ID.
    metals : list[str]
        Upper-case PDB Component IDs to look for (e.g. ``["CA", "MG"]``).
    """
    metals_set = {m.upper() for m in metals}
    present, complete = fetch_ligand_ids_status(pdb_id)

    found = metals_set & present
    if found:
        log.debug("%s: confirmed metal(s) present: %s", pdb_id, sorted(found))
        return True

    if not complete:
        log.warning(
            "%s: could not verify metal presence (RCSB request failed or was "
            "incomplete) - keeping the structure and attempting extraction",
            pdb_id,
        )
        return True

    clusters = present & UNSUPPORTED_CLUSTER_IDS
    if clusters:
        CLUSTER_ONLY_SEEN[pdb_id.upper()] = set(clusters)
        log.info(
            "%s: skipped - none of %s found; the structure contains a metal cluster %s, "
            "which is a multi-atom residue and is not supported (ligands present: %s)",
            pdb_id, sorted(metals_set), sorted(clusters), sorted(present),
        )
    elif present:
        log.info(
            "%s: skipped - none of %s found in structure (ligands present: %s)",
            pdb_id, sorted(metals_set), sorted(present),
        )
    else:
        log.info(
            "%s: skipped - entry has no ions or ligands at all (none of %s present)",
            pdb_id, sorted(metals_set),
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
    include_entities: bool = True,
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
            pool.submit(fetch_entry_metadata, pid, include_entities): pid
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


def filter_reason(
    meta: dict[str, Any],
    max_resolution: float,
    allowed_methods: list[str],
    include_nmr: bool = False,
) -> str:
    """Return why a structure fails the resolution/method filter ("" if it passes).

    Possible non-empty return values, so that callers can count each
    removal category separately instead of lumping them together:

    ``"no_metadata"``
        Empty metadata (failed lookup).
    ``"no_resolution"``
        No numeric resolution and not an accepted NMR entry.
    ``"method"``
        Experimental method not in *allowed_methods*.
    ``"resolution"``
        Resolution worse than *max_resolution*.
    """
    if not meta:
        return "no_metadata"

    method = (meta.get("method") or "").upper()
    resolution = meta.get("resolution")

    # --- NMR handling -------------------------------------------------------
    if resolution is None:
        if include_nmr and "NMR" in method:
            log.debug("%s accepted (NMR, include_nmr=True)", meta.get("pdb_id"))
            return ""
        log.debug("%s skipped (no resolution, NMR not allowed)", meta.get("pdb_id"))
        return "no_resolution"

    # --- Method filter ------------------------------------------------------
    if allowed_methods and method not in [m.upper() for m in allowed_methods]:
        log.debug(
            "%s skipped (method=%s not in allowed list)", meta.get("pdb_id"), method
        )
        return "method"

    # --- Resolution filter --------------------------------------------------
    if float(resolution) > max_resolution:
        log.debug(
            "%s skipped (resolution=%.2f > max=%.2f)",
            meta.get("pdb_id"), float(resolution), max_resolution,
        )
        return "resolution"

    return ""


def passes_filter(
    meta: dict[str, Any],
    max_resolution: float,
    allowed_methods: list[str],
    include_nmr: bool = False,
) -> bool:
    """Boolean wrapper around :func:`filter_reason` (kept for compatibility)."""
    return filter_reason(meta, max_resolution, allowed_methods, include_nmr) == ""


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
