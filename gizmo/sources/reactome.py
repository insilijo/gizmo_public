"""
Reactome REST API client and bulk species loader.

License: Reactome data is CC BY 4.0 — https://reactome.org/license
API docs: https://reactome.org/ContentService/

Key endpoints used:
  /data/pathways/top/{species}         — top-level pathways (29 for human)
  /data/pathway/{id}/containedEvents   — direct children of a pathway
  /data/query/{id}                     — generic entity / reaction lookup

API quirk: containedEvents returns a mixed list of dicts (full event objects)
and ints (bare dbIds). Ints are normalised by fetching /data/query/{dbId}.

Reaction schemaClasses included: Reaction, BlackBoxEvent, Polymerisation,
  Depolymerisation.  Sub-pathways (Pathway, TopLevelPathway) are recursed
  into but not added as nodes.
"""

from __future__ import annotations

import json
import logging
import time
from collections import deque
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Optional

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from gizmo.schema import EdgeRole, MetaboliteNode, ReactionEdge, ReactionNode

log = logging.getLogger(__name__)

_BASE = "https://reactome.org/ContentService/"
_SPECIES_DEFAULT = "Homo sapiens"

_REACTION_CLASSES = {"Reaction", "BlackBoxEvent", "Polymerisation", "Depolymerisation"}
_PATHWAY_CLASSES  = {"Pathway", "TopLevelPathway"}


def _make_session() -> requests.Session:
    """
    Build a requests.Session with:
    - Retry on connection errors and 429/500/502/503/504 with exponential backoff
    - Connection pool sized for concurrent use (max 20 connections)
    """
    session = requests.Session()
    session.headers.update({"Accept": "application/json"})
    retry = Retry(
        total=5,
        backoff_factor=1.0,          # waits: 1, 2, 4, 8, 16 s
        status_forcelist={429, 500, 502, 503, 504},
        allowed_methods={"GET"},
        raise_on_status=False,
    )
    adapter = HTTPAdapter(
        max_retries=retry,
        pool_connections=10,
        pool_maxsize=20,
    )
    session.mount("https://", adapter)
    session.mount("http://",  adapter)
    return session


# ---------------------------------------------------------------------------
# ReactomeClient — thin REST wrapper
# ---------------------------------------------------------------------------

class ReactomeClient:
    """Thin wrapper around the Reactome Content Service REST API."""

    def __init__(self, base_url: str = _BASE, timeout: int = 45) -> None:
        self.base_url = base_url
        self.session = _make_session()
        self.timeout = timeout

    # ------------------------------------------------------------------
    # Low-level
    # ------------------------------------------------------------------

    def _get(self, path: str, **params: Any) -> Any:
        url = self.base_url + path
        # Use a prepared request with manually-built URL to avoid requests
        # percent-encoding spaces in path segments (Reactome rejects %20 in
        # species names — it requires literal spaces, e.g. "Homo sapiens").
        import urllib.parse
        if params:
            qs = urllib.parse.urlencode(params)
            url = f"{url}?{qs}"
        resp = self.session.get(url, timeout=self.timeout)
        resp.raise_for_status()
        return resp.json()

    # ------------------------------------------------------------------
    # Pathways
    # ------------------------------------------------------------------

    def top_pathways(self, species: str = _SPECIES_DEFAULT) -> list[dict]:
        """Return top-level pathway stubs for a species."""
        return self._get(f"data/pathways/top/{species}")

    def pathway_events(self, pathway_stid: str) -> list[dict]:
        """
        Return events (reactions + sub-pathways) directly contained in a pathway.

        Normalises the Reactome API quirk where some items are bare int dbIds
        rather than full event objects — those are resolved via /data/query/{dbId}.
        Returns an empty list for pathways that return 404 (deleted/moved).
        """
        try:
            raw = self._get(f"data/pathway/{pathway_stid}/containedEvents")
        except requests.HTTPError as exc:
            if exc.response is not None and exc.response.status_code == 404:
                log.debug("pathway_events 404 for %s — skipping", pathway_stid)
                return []
            raise

        result: list[dict] = []
        for item in raw:
            if isinstance(item, dict):
                result.append(item)
            elif isinstance(item, int):
                # Bare dbId — resolve to a minimal stub
                try:
                    obj = self._get(f"data/query/{item}")
                    result.append(obj)
                except Exception as exc:
                    log.debug("Could not resolve dbId %d: %s", item, exc)
        return result

    # Keep old name as alias so existing code doesn't break
    def pathway_reactions(self, pathway_stid: str) -> list[dict]:
        return self.pathway_events(pathway_stid)

    # ------------------------------------------------------------------
    # Reaction detail → graph primitives
    # ------------------------------------------------------------------

    def reaction_detail(self, reaction_stid: str) -> dict:
        return self._get(f"data/query/{reaction_stid}")

    def parse_reaction(
        self,
        detail: dict,
        pathway_stids: Optional[list[str]] = None,
        catalyst_details: Optional[dict[int, dict]] = None,
        gene_symbols: Optional[list[str]] = None,
    ) -> tuple[ReactionNode, list[MetaboliteNode], list[ReactionEdge]]:
        """
        Convert a Reactome reaction detail dict into graph primitives.

        Returns (ReactionNode, [MetaboliteNode, ...], [ReactionEdge, ...]).
        Metabolite nodes are compartment-aware; node IDs are
        "CHEBI:XXXXX@compartment" when both are available, else "reactome:{stid}".

        Parameters
        ----------
        catalyst_details : dict[int, dict], optional
            Pre-fetched full CatalystActivity objects keyed by dbId.
            Used to read activity.ecNumber without inline API calls.
        gene_symbols : list[str], optional
            Pre-fetched gene symbols (from ReactomeLoader's cached
            referenceEntities call).  When absent, fetched inline via
            /data/participants/{stId}/referenceEntities.
        """
        rxn_stid = str(detail.get("stId") or detail.get("dbId") or "unknown")
        rxn_id = f"reactome:{rxn_stid}"

        # ---- EC numbers -------------------------------------------------------
        # catalystActivity stubs in reaction detail lack activity.ecNumber.
        # We need the full CatalystActivity object from /data/query/{cat_dbId}.
        ec_list: list[str] = []
        for cat in detail.get("catalystActivity", []):
            if not isinstance(cat, dict):
                continue
            cat_dbid = cat.get("dbId")
            if catalyst_details and cat_dbid in catalyst_details:
                full_cat = catalyst_details[cat_dbid]
            elif cat_dbid:
                try:
                    full_cat = self._get(f"data/query/{cat_dbid}")
                except Exception:
                    full_cat = {}
            else:
                full_cat = {}
            if ec := full_cat.get("activity", {}).get("ecNumber"):
                ec_list.append(ec)

        # ---- Gene symbols -----------------------------------------------------
        # /data/participants/{stId}/referenceEntities returns all protein
        # participants in a single call — far simpler than traversing the
        # physicalEntity → Complex → EWAS stub chain.
        if gene_symbols is not None:
            gene_symbols_final = list(gene_symbols)
        else:
            try:
                refs = self._get(f"data/participants/{rxn_stid}/referenceEntities")
                gene_symbols_final = list(dict.fromkeys(
                    e["name"][0]
                    for e in refs
                    if e.get("moleculeType") == "Protein" and e.get("name")
                ))
            except Exception:
                gene_symbols_final = []

        species_name: Optional[str] = None
        if sp := detail.get("species"):
            if isinstance(sp, list) and sp:
                species_name = sp[0].get("displayName")
            elif isinstance(sp, dict):
                species_name = sp.get("displayName")

        rxn_node = ReactionNode(
            node_id=rxn_id,
            reactome_id=rxn_stid,
            name=detail.get("displayName", rxn_stid),
            reversible=detail.get("isReversible", False),
            ec_numbers=ec_list,
            gene_symbols=gene_symbols_final,
            pathways=pathway_stids or [],
            species=species_name,
        )

        metabolites: list[MetaboliteNode] = []
        edges: list[ReactionEdge] = []

        def _parse(entries: list[dict], role: EdgeRole) -> None:
            for entry in entries:
                if not isinstance(entry, dict):
                    continue
                met_node, edge = _participant_to_node_edge(entry, rxn_id, role)
                if met_node:
                    metabolites.append(met_node)
                    edges.append(edge)

        _parse(detail.get("input", []),  EdgeRole.SUBSTRATE)
        _parse(detail.get("output", []), EdgeRole.PRODUCT)
        _parse(detail.get("catalystActivity", []), EdgeRole.MODIFIER)

        return rxn_node, metabolites, edges


# ---------------------------------------------------------------------------
# ReactomeLoader — bulk species-level loader
# ---------------------------------------------------------------------------

class ReactomeLoader:
    """
    Load all reactions for a species from Reactome via BFS pathway traversal.

    Strategy
    --------
    Phase 1 — collect:
        BFS from top-level pathways → sub-pathways → leaf reaction stIDs.
        Pathway event lists are cached locally (most are tiny JSON files).

    Phase 2 — fetch details:
        Fetch reaction detail for each unique stID, in parallel (default 5
        workers). Raw responses are cached; re-runs are near-instant.
        Workers share a session with automatic retry + exponential backoff.

    Phase 3 — parse:
        Convert each detail dict to (ReactionNode, [MetaboliteNode], [ReactionEdge])
        and add to a GizmoGraph.

    Usage
    -----
    ::

        loader = ReactomeLoader(cache_dir="data/raw/reactome")
        mg = loader.load_species("Homo sapiens")   # ~16 000 reactions, ~3 min first run
        # or
        mg = loader.load_pathways(["R-HSA-70171", "R-HSA-8964208"])
    """

    def __init__(
        self,
        client: Optional[ReactomeClient] = None,
        cache_dir: str | Path = "data/raw/reactome",
        max_workers: int = 5,
        request_delay: float = 0.05,   # seconds between requests per worker
    ) -> None:
        self.client = client or ReactomeClient()
        self.cache_dir = Path(cache_dir)
        self.request_delay = request_delay
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.max_workers = max_workers

    # ------------------------------------------------------------------
    # Caching helpers
    # ------------------------------------------------------------------

    def _cache_path(self, key: str) -> Path:
        safe = key.replace("/", "_").replace(":", "_")
        return self.cache_dir / f"{safe}.json"

    def _cached_get(self, key: str, fetch_fn) -> Any:
        path = self._cache_path(key)
        if path.exists():
            return json.loads(path.read_text())
        result = fetch_fn()
        path.write_text(json.dumps(result))
        return result

    def _fetch_events(self, pathway_stid: str) -> list[dict]:
        return self._cached_get(
            f"events_{pathway_stid}",
            lambda: self.client.pathway_events(pathway_stid),
        )

    def _fetch_gene_symbols(self, stid: str) -> list[str]:
        """
        Return protein gene symbols for a reaction via /data/participants/{stId}/referenceEntities.
        Results are disk-cached.  Returns [] on any error.
        """
        def _fetch() -> list[str]:
            try:
                refs = self.client._get(f"data/participants/{stid}/referenceEntities")
                return list(dict.fromkeys(
                    e["name"][0]
                    for e in refs
                    if e.get("moleculeType") == "Protein" and e.get("name")
                ))
            except Exception as exc:
                log.debug("referenceEntities failed %s: %s", stid, exc)
                return []

        return self._cached_get(f"genes_{stid}", _fetch)

    def _fetch_detail(self, stid: str) -> Optional[dict]:
        cache_path = self._cache_path(f"detail_{stid}")
        if cache_path.exists():
            return json.loads(cache_path.read_text())
        # Not cached — hit the network with a small inter-request delay
        if self.request_delay > 0:
            time.sleep(self.request_delay)
        try:
            result = self.client.reaction_detail(stid)
            cache_path.write_text(json.dumps(result))
            return result
        except Exception as exc:
            log.debug("Detail fetch failed %s: %s", stid, exc)
            return None

    # ------------------------------------------------------------------
    # Phase 1: collect all reaction stIDs via BFS
    # ------------------------------------------------------------------

    def _collect_reaction_stids(
        self,
        seed_events: list[dict],
        *,
        pathway_ancestry: Optional[dict[str, list[str]]] = None,
    ) -> dict[str, list[str]]:
        """
        BFS from seed events; returns {reaction_stid: [all_containing_pathway_stids]}.

        Reactome's /containedEvents endpoint returns ALL descendants of a pathway
        in a flat list (not just direct children).  This means reactions appear at
        depth-1 in every ancestor's event list.  To capture the full hierarchy we
        accumulate pathways per reaction as a set: each time a reaction is
        encountered during the processing of a particular pathway node, that
        pathway stId is added to the reaction's membership set.
        """
        if pathway_ancestry is None:
            pathway_ancestry = {}

        queue: deque[tuple[dict, list[str]]] = deque(
            (ev, []) for ev in seed_events
        )
        seen_pathways: set[str] = set()
        # Use sets for accumulation; converted to sorted lists before returning.
        reaction_pathways: dict[str, set[str]] = {}

        while queue:
            event, ancestors = queue.popleft()
            stid = str(event.get("stId") or event.get("dbId") or "")
            schema = event.get("schemaClass", "")

            if schema in _PATHWAY_CLASSES or not schema:
                if stid in seen_pathways:
                    continue
                seen_pathways.add(stid)
                children = self._fetch_events(stid)
                for child in children:
                    queue.append((child, ancestors + ([stid] if stid else [])))

            elif schema in _REACTION_CLASSES:
                if stid:
                    if stid not in reaction_pathways:
                        reaction_pathways[stid] = set(ancestors)
                    else:
                        # Accumulate: add any pathway stIds this traversal contributes.
                        reaction_pathways[stid].update(ancestors)

        return {stid: sorted(pws) for stid, pws in reaction_pathways.items()}

    # ------------------------------------------------------------------
    # Phase 2+3: fetch details + build graph
    # ------------------------------------------------------------------

    def _build_graph(
        self,
        reaction_stids: dict[str, list[str]],
        mg,
    ) -> None:
        from gizmo.graph.network import GizmoGraph

        total = len(reaction_stids)
        # Separate cached vs uncached to report accurately
        cached_count = sum(1 for stid in reaction_stids if self._cache_path(f"detail_{stid}").exists())
        uncached_count = total - cached_count
        log.info(
            "Fetching %d reactions (%d cached / %d to download, %d workers) …",
            total, cached_count, uncached_count, self.max_workers,
        )

        items = list(reaction_stids.items())

        # Each task: fetch (possibly cached) + parse, then write to graph under a lock
        import threading
        _lock = threading.Lock()

        def _task(stid: str, pathway_stids: list[str]):
            detail = self._fetch_detail(stid)
            if not detail:
                return
            try:
                # Pre-fetch full CatalystActivity objects for EC numbers (disk-cached).
                # Reaction detail stubs lack activity.ecNumber — it only appears in the
                # full CatalystActivity fetched via /data/query/{cat_dbId}.
                catalyst_details: dict[int, dict] = {}
                for cat in detail.get("catalystActivity", []):
                    if isinstance(cat, dict) and (cat_dbid := cat.get("dbId")):
                        cat_full = self._fetch_detail(str(cat_dbid))
                        if cat_full:
                            catalyst_details[cat_dbid] = cat_full

                # Pre-fetch + cache gene symbols from referenceEntities endpoint
                # (one call per reaction, avoids complex stub traversal).
                gene_syms = self._fetch_gene_symbols(stid)

                rxn_node, mets, edges = self.client.parse_reaction(
                    detail,
                    pathway_stids=pathway_stids,
                    catalyst_details=catalyst_details or None,
                    gene_symbols=gene_syms,
                )
                with _lock:
                    mg.add_reaction(rxn_node)
                    mg.add_metabolites(mets)
                    mg.add_edges(edges)
            except Exception as exc:
                log.debug("parse_reaction failed %s: %s", stid, exc)

        try:
            from rich.progress import Progress, BarColumn, TaskProgressColumn, TimeRemainingColumn, SpinnerColumn
            with Progress(SpinnerColumn(), "[progress.description]{task.description}",
                          BarColumn(), TaskProgressColumn(), TimeRemainingColumn()) as prog:
                task = prog.add_task(f"Fetching {total} reactions", total=total)
                with ThreadPoolExecutor(max_workers=self.max_workers) as pool:
                    futures = {pool.submit(_task, stid, pw): stid for stid, pw in items}
                    for fut in as_completed(futures):
                        fut.result()   # surface exceptions to log
                        prog.advance(task)
        except ImportError:
            done = 0
            with ThreadPoolExecutor(max_workers=self.max_workers) as pool:
                futures = {pool.submit(_task, stid, pw): stid for stid, pw in items}
                for fut in as_completed(futures):
                    try:
                        fut.result()
                    except Exception as exc:
                        log.debug("task failed: %s", exc)
                    done += 1
                    if done % 500 == 0:
                        log.info("  %d / %d reactions processed", done, total)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def load_species(
        self,
        species: str = _SPECIES_DEFAULT,
        mg=None,
    ):
        """
        Load all reactions for a species into a GizmoGraph.
        First run ~3 min for human (16 000 reactions); subsequent runs use cache.
        """
        from gizmo.graph.network import GizmoGraph
        if mg is None:
            mg = GizmoGraph()

        log.info("Collecting reaction stIDs for '%s' via BFS …", species)
        top = self.client.top_pathways(species)
        reaction_stids = self._collect_reaction_stids(top)
        log.info("Found %d unique reactions", len(reaction_stids))

        self._build_graph(reaction_stids, mg)
        log.info("Done. %s", mg)
        return mg

    def load_pathways(
        self,
        pathway_stids: list[str],
        mg=None,
    ):
        """
        Load reactions from a specific list of pathway stIDs.
        Recurses into sub-pathways automatically.
        """
        from gizmo.graph.network import GizmoGraph
        if mg is None:
            mg = GizmoGraph()

        seed = []
        for pid in pathway_stids:
            try:
                stub = self.client.reaction_detail(pid)   # gets the pathway object
                seed.append(stub)
            except Exception:
                seed.append({"stId": pid, "schemaClass": "Pathway"})

        reaction_stids = self._collect_reaction_stids(seed)
        log.info("Found %d reactions across %d pathways", len(reaction_stids), len(pathway_stids))
        self._build_graph(reaction_stids, mg)
        return mg


# ---------------------------------------------------------------------------
# Patch ReactomeClient with the method ReactomeLoader needs
# ---------------------------------------------------------------------------

def _client_fetch_and_parse(
    self: ReactomeClient,
    stid: str,
    pathway_stids: list[str],
) -> tuple[ReactionNode, list[MetaboliteNode], list[ReactionEdge]]:
    """Fetch reaction detail and parse — used by ReactomeLoader threadpool."""
    detail = self.reaction_detail(stid)
    return self.parse_reaction(detail, pathway_stids=pathway_stids)

ReactomeClient._fetch_and_parse = _client_fetch_and_parse  # type: ignore[attr-defined]


# ---------------------------------------------------------------------------
# Helper: gene name extraction from Reactome physicalEntity
# ---------------------------------------------------------------------------

def _extract_gene_names(pe: dict) -> list[str]:
    """
    Extract gene names from a Reactome physicalEntity dict.

    Handles:
    - EntityWithAccessionedSequence (EWAS) — referenceEntity.geneName
    - Complex / DefinedSet — iterates hasComponent / hasMember / hasCandidate
      one level deep (sub-complexes are skipped; their stubs rarely carry
      referenceEntity without an additional fetch)

    Returns a deduplicated list preserving first-occurrence order.
    """
    seen: dict[str, None] = {}

    def _from_entity(entity: dict) -> None:
        ref = entity.get("referenceEntity", {}) if isinstance(entity, dict) else {}
        for gn in ref.get("geneName", []):
            seen[gn] = None

    _from_entity(pe)

    # Flatten one level of complex members
    for key in ("hasComponent", "hasMember", "hasCandidate"):
        for member in (pe.get(key) or []):
            if isinstance(member, dict):
                _from_entity(member)

    return list(seen)


# ---------------------------------------------------------------------------
# Participant parsing helper
# ---------------------------------------------------------------------------

def _participant_to_node_edge(
    entry: dict, rxn_id: str, role: EdgeRole
) -> tuple[Optional[MetaboliteNode], Optional[ReactionEdge]]:
    """
    Map a Reactome participant dict to MetaboliteNode + ReactionEdge.
    Returns (None, None) for non-small-molecule participants.
    """
    schema_class = entry.get("schemaClass", "")
    if schema_class not in {"SimpleEntity", "ChemicalDrug", "OtherEntity"}:
        return None, None

    stid = entry.get("stId") or str(entry.get("dbId", ""))
    chebi_xrefs = [
        ref.get("identifier")
        for ref in entry.get("crossReference", [])
        if ref.get("databaseName") == "ChEBI"
    ]
    chebi_id = f"CHEBI:{chebi_xrefs[0]}" if chebi_xrefs else None

    compartment_name: Optional[str] = None
    if comp := entry.get("compartment"):
        if isinstance(comp, dict):
            compartment_name = comp.get("displayName") or comp.get("name")

    node_id = chebi_id or f"reactome:{stid}"
    if compartment_name:
        node_id = f"{node_id}@{compartment_name}"

    met_node = MetaboliteNode(
        node_id=node_id,
        chebi_id=chebi_id,
        reactome_id=str(stid),
        name=entry.get("displayName", str(stid)),
        compartment=compartment_name,
    )

    src, tgt = (node_id, rxn_id) if role == EdgeRole.SUBSTRATE else (rxn_id, node_id)

    edge = ReactionEdge(
        source=src,
        target=tgt,
        role=role,
        stoichiometry=float(entry.get("stoichiometry", 1.0)),
        compartment=compartment_name,
    )

    return met_node, edge
