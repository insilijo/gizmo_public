"""Populate graph metabolite nodes with PubChem synonyms.

The MetaboliteMapper (gizmo/evidence/mappers.py:188) already looks for a
``synonyms`` attribute on each metabolite node when resolving a free-text
feature ID. The lookup was always intended to draw from PubChem's
synonym list, but the enrichment that populates the attribute was never
implemented — so for users with Metabolon HD4 long-form names (e.g.
"S-1-pyrroline-5-carboxylate", "homovanillate (HVA)"), the mapper had
no aliases to match against and ~95% of features missed the graph.

This module fixes that. ``enrich_pubchem_synonyms(mg)`` queries PubChem
PUG REST for every metabolite node carrying a ``pubchem_cid`` attribute,
caches the synonym lists on disk, and writes a lowercase ``synonyms``
set into each node.

Called once during ``build_human_graph()``. Cached so subsequent runs
incur zero network cost.
"""
from __future__ import annotations

import json
import logging
import time
import urllib.error
import urllib.request
from pathlib import Path

log = logging.getLogger(__name__)


PUG_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid"
DEFAULT_BATCH = 100
DEFAULT_SLEEP = 0.2  # s between batches; PubChem allows ≤5 req/sec
DEFAULT_TIMEOUT = 30
DEFAULT_CACHE = Path("data/raw/pubchem/synonyms.json")


def _query_batch(cids: list[str], timeout: float = DEFAULT_TIMEOUT) -> dict[str, list[str]]:
    """Query PubChem PUG REST for synonyms of a batch of CIDs.

    Returns ``{cid_str: [synonym, ...]}``. Missing CIDs are absent from
    the returned dict.
    """
    if not cids:
        return {}
    url = f"{PUG_BASE}/{','.join(cids)}/synonyms/JSON"
    try:
        with urllib.request.urlopen(url, timeout=timeout) as r:
            payload = json.load(r)
    except (urllib.error.URLError, json.JSONDecodeError) as exc:
        log.warning("PubChem PUG batch (%d CIDs) failed: %s", len(cids), exc)
        return {}
    out: dict[str, list[str]] = {}
    info = payload.get("InformationList", {}).get("Information", [])
    for entry in info:
        cid = str(entry.get("CID"))
        syns = entry.get("Synonym") or []
        out[cid] = [s for s in syns if s and isinstance(s, str)]
    return out


def _inchikey_to_cid(inchikey: str, timeout: float = DEFAULT_TIMEOUT) -> str | None:
    """Resolve a single InChIKey to a PubChem CID via PUG REST.

    Returns the first matching CID as a string, or None if no match.
    PubChem typically returns one CID per full 27-char InChIKey.
    """
    if not inchikey:
        return None
    url = f"{PUG_BASE.replace('/cid', '/inchikey')}/{inchikey}/cids/JSON"
    try:
        with urllib.request.urlopen(url, timeout=timeout) as r:
            payload = json.load(r)
    except (urllib.error.URLError, json.JSONDecodeError):
        return None
    cids = payload.get("IdentifierList", {}).get("CID", [])
    return str(cids[0]) if cids else None


def enrich_pubchem_synonyms(
    mg,
    *,
    cache_path: Path | str | None = None,
    batch_size: int = DEFAULT_BATCH,
    sleep_between_batches: float = DEFAULT_SLEEP,
    max_synonyms_per_node: int = 50,
) -> int:
    """Populate ``synonyms`` (lowercased set) on each metabolite node
    that has a ``pubchem_cid`` attribute.

    Parameters
    ----------
    mg
        GizmoGraph
    cache_path
        Where to persist the {cid: [synonyms]} cache. Defaults to
        ``data/raw/pubchem/synonyms.json`` (relative to cwd, like the
        Reactome cache).
    batch_size
        CIDs per PUG REST call. PubChem accepts ~200; 100 is safe.
    sleep_between_batches
        Polite throttle. PubChem rate limit is 5 req/s.
    max_synonyms_per_node
        Cap synonyms stored per node. PubChem returns hundreds for some
        common compounds (water has >500); keeping all bloats memory
        and dilutes the mapper's name index without much marginal
        recall. 50 covers the common-name + IUPAC + InChI variants.

    Returns
    -------
    int
        Number of metabolite nodes updated with synonyms.
    """
    g = mg.graph
    cache_path = Path(cache_path) if cache_path else DEFAULT_CACHE
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache: dict[str, list[str]] = {}
    if cache_path.exists():
        try:
            cache = json.loads(cache_path.read_text())
        except json.JSONDecodeError:
            log.warning("PubChem synonym cache at %s is corrupt; rebuilding",
                         cache_path)
            cache = {}

    # Collect every metabolite node + its CID. If a node has no CID but
    # has an InChIKey, look the CID up first. PubChem InChIKey→CID
    # results are cached on the inchikey side of the cache to avoid
    # re-querying.
    inchikey_cid_cache: dict[str, str] = cache.setdefault("__inchikey_cid__", {})
    cid_to_nodes: dict[str, list[str]] = {}
    n_via_inchikey = 0
    nodes_needing_inchikey_lookup: list[tuple[str, str]] = []  # (nid, inchikey)
    for nid, attrs in g.nodes(data=True):
        if attrs.get("node_type") != "metabolite":
            continue
        cid = attrs.get("pubchem_cid")
        if cid:
            try:
                cid_clean = str(int(float(cid)))
                cid_to_nodes.setdefault(cid_clean, []).append(nid)
                continue
            except (TypeError, ValueError):
                pass
        # No PubChem CID directly — try InChIKey path
        ik = attrs.get("inchikey") or ""
        if ik and len(ik) >= 14:
            if ik in inchikey_cid_cache:
                resolved = inchikey_cid_cache[ik]
                if resolved:
                    cid_to_nodes.setdefault(resolved, []).append(nid)
                    n_via_inchikey += 1
            else:
                nodes_needing_inchikey_lookup.append((nid, ik))

    # Resolve novel InChIKeys to CIDs (one network call each)
    if nodes_needing_inchikey_lookup:
        msg = (f"PubChem synonym enrichment: resolving "
               f"{len(nodes_needing_inchikey_lookup)} new InChIKeys → CID")
        log.info(msg); print(f"  {msg}", flush=True)
        for i, (nid, ik) in enumerate(nodes_needing_inchikey_lookup):
            cid = _inchikey_to_cid(ik)
            inchikey_cid_cache[ik] = cid or ""
            if cid:
                cid_to_nodes.setdefault(cid, []).append(nid)
                # Self-describe: write the resolved CID back onto the
                # node so downstream tools (and re-runs without cache)
                # can find it directly.
                if not g.nodes[nid].get("pubchem_cid"):
                    g.nodes[nid]["pubchem_cid"] = cid
                n_via_inchikey += 1
            if i % 100 == 0 and i > 0:
                log.info("  InChIKey→CID: %d/%d done",
                          i, len(nodes_needing_inchikey_lookup))
                print(f"    InChIKey→CID: {i}/{len(nodes_needing_inchikey_lookup)} done",
                      flush=True)
            time.sleep(sleep_between_batches)
        # Persist the inchikey→cid cache early in case synonym fetch fails
        cache_path.write_text(json.dumps(cache, indent=2))

    if n_via_inchikey:
        log.info("  + %d nodes routed through InChIKey→CID lookup",
                  n_via_inchikey)

    if not cid_to_nodes:
        log.info("PubChem synonym enrichment: no metabolite nodes have "
                  "pubchem_cid or resolvable InChIKey; nothing to do")
        return 0

    todo = sorted(c for c in cid_to_nodes if c not in cache)
    log.info("PubChem synonym enrichment: %d metabolite nodes with CIDs, "
              "%d need querying (%d cached)",
              len(cid_to_nodes), len(todo), len(cache))

    # Batch-query missing CIDs
    for i in range(0, len(todo), batch_size):
        batch = todo[i:i + batch_size]
        results = _query_batch(batch)
        for cid in batch:
            cache[cid] = results.get(cid, [])
        if (i // batch_size) % 10 == 0 and i > 0:
            log.info("  PubChem synonym fetch: %d/%d CIDs done", i, len(todo))
        time.sleep(sleep_between_batches)

    cache_path.write_text(json.dumps(cache, indent=2))

    # Apply synonyms to nodes. Stored as a sorted *list* (not set) because
    # gizmo.evidence.mappers.MetaboliteMapper.__init__ checks
    # ``isinstance(syns, (list, tuple))`` before iterating; a set silently
    # drops the synonyms.
    n_nodes_updated = 0
    n_synonyms_total = 0
    for cid, nids in cid_to_nodes.items():
        syns = cache.get(cid, [])[:max_synonyms_per_node]
        syn_set = {s.lower() for s in syns}
        if not syn_set:
            continue
        for nid in nids:
            existing = g.nodes[nid].get("synonyms")
            if existing is None:
                merged = syn_set
            else:
                merged = set(existing) | syn_set
            g.nodes[nid]["synonyms"] = sorted(merged)
            n_nodes_updated += 1
            n_synonyms_total += len(syn_set)

    log.info("PubChem synonym enrichment: %d nodes enriched with %d "
              "synonyms total (avg %.1f/node)",
              n_nodes_updated, n_synonyms_total,
              n_synonyms_total / max(1, n_nodes_updated))
    return n_nodes_updated
