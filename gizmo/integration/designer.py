"""Constrained metabolite queries for prospective panel design.

Three orthogonal constraint streams, each producing a list of candidate
metabolite node_ids with an evidence-style score:

  pathways  : metabolite ∈ reactions ∈ pathway_set
              → score = count of distinct pathways this metabolite belongs to
                (within the chosen set)

  diseases  : metabolite reachable from disease_node within max_hops
              along DiseaseEdgeType edges
              → score = count of distinct disease anchors that reach this
                metabolite

  taxa      : delegated to gemma.integration.panel_design (microbiome-only),
              not implemented here

Returned per stream: DataFrame with columns
   metabolite_node_id, score, source_ids (semicolon-joined)

The orchestrator combines streams with max-of-stream + source-tag union,
so a metabolite covered by both a pathway and a disease shows up once
with both source tags preserved.
"""
from __future__ import annotations

from collections import deque
from typing import Iterable

import pandas as pd


def metabolites_for_pathways(mg, pathway_ids: Iterable[str]) -> pd.DataFrame:
    """Metabolites participating in reactions of any given Reactome pathway.

    score = count of distinct pathways (from the input set) each metabolite
    is reached through. Higher = better-covered by the chosen set.
    """
    g = mg.graph if hasattr(mg, "graph") else mg
    pids = set(pathway_ids)
    if not pids:
        return pd.DataFrame(columns=["metabolite_node_id", "score", "source_ids"])

    reaction_to_pways: dict[str, set[str]] = {}
    for n, d in g.nodes(data=True):
        if d.get("node_type") != "reaction":
            continue
        hits = pids & set(d.get("pathways") or [])
        if hits:
            reaction_to_pways[n] = hits

    metab_hits: dict[str, set[str]] = {}
    for u, v, ed in g.edges(data=True):
        role = (ed.get("role") or ed.get("edge_type") or "").lower()
        if role not in ("substrate", "product"):
            continue
        u_t = g.nodes.get(u, {}).get("node_type")
        v_t = g.nodes.get(v, {}).get("node_type")
        if u_t == "metabolite" and v_t == "reaction":
            m, r = u, v
        elif v_t == "metabolite" and u_t == "reaction":
            m, r = v, u
        else:
            continue
        if r in reaction_to_pways:
            metab_hits.setdefault(m, set()).update(reaction_to_pways[r])

    rows = [
        {
            "metabolite_node_id": m,
            "score": float(len(hits)),
            "source_ids": ";".join(sorted(hits)),
        }
        for m, hits in metab_hits.items()
    ]
    if not rows:
        return pd.DataFrame(columns=["metabolite_node_id", "score", "source_ids"])
    return pd.DataFrame(rows).sort_values("score", ascending=False)


def metabolites_for_diseases(
    mg,
    disease_ids: Iterable[str],
    *,
    max_hops: int = 3,
    expand_via_pathways: bool = True,
) -> pd.DataFrame:
    """Metabolites reachable from any of the given disease node ids.

    Two harvesters, both optional, results unioned:

    1. **Direct graph walk** (always on): BFS up to ``max_hops`` hops from
       each disease, collect any metabolite node reached. Works cleanly
       for diseases linked to enzymes whose reactions have substrates/
       products in the graph.

    2. **Pathway expansion** (``expand_via_pathways=True``, default):
       Along the same BFS, also collect every reaction reached, then
       read each reaction's ``pathways`` attribute, then harvest all
       metabolites participating in any reaction of those pathways.
       Covers diseases linked to non-enzymatic genes whose roles are
       curated as pathway components (signaling/regulatory). Captures
       the "non-canonical metabolic function" route — a disease that
       has no direct enzymatic edge but is annotated to a Reactome
       pathway via its genes still surfaces relevant metabolites.

    Each output row carries a ``via`` column tagging the harvester(s).
    """
    g = mg.graph if hasattr(mg, "graph") else mg
    targets = [d for d in disease_ids if d in g]
    if not targets:
        return pd.DataFrame(columns=["metabolite_node_id", "score", "source_ids", "via"])

    metab_via_direct: dict[str, set[str]] = {}
    pathway_ids: set[str] = set()
    pathway_to_disease: dict[str, set[str]] = {}

    for did in targets:
        seen = {did}
        frontier: deque[tuple[str, int]] = deque([(did, 0)])
        local_metab: set[str] = set()
        local_pways: set[str] = set()
        while frontier:
            nid, depth = frontier.popleft()
            attrs = g.nodes[nid]
            ntype = attrs.get("node_type")
            if ntype == "metabolite":
                local_metab.add(nid)
            if expand_via_pathways and ntype == "reaction":
                for p in attrs.get("pathways") or []:
                    local_pways.add(str(p))
            if depth >= max_hops:
                continue
            neighbors = []
            for _, nbr in g.out_edges(nid):
                neighbors.append(nbr)
            for nbr, _ in g.in_edges(nid):
                neighbors.append(nbr)
            for nbr in neighbors:
                if nbr in seen:
                    continue
                seen.add(nbr)
                frontier.append((nbr, depth + 1))
        for m in local_metab:
            metab_via_direct.setdefault(m, set()).add(did)
        for p in local_pways:
            pathway_ids.add(p)
            pathway_to_disease.setdefault(p, set()).add(did)

    metab_via_expand: dict[str, set[str]] = {}
    if expand_via_pathways and pathway_ids:
        expand_df = metabolites_for_pathways(mg, pathway_ids)
        for _, r in expand_df.iterrows():
            mid = r["metabolite_node_id"]
            hit_pways = set(str(r["source_ids"]).split(";")) & pathway_ids
            for p in hit_pways:
                metab_via_expand.setdefault(mid, set()).update(
                    pathway_to_disease.get(p, set())
                )

    all_metab = set(metab_via_direct) | set(metab_via_expand)
    rows = []
    for m in all_metab:
        direct_diseases = metab_via_direct.get(m, set())
        expand_diseases = metab_via_expand.get(m, set())
        all_d = direct_diseases | expand_diseases
        via_tags = []
        if direct_diseases:
            via_tags.append("direct")
        if expand_diseases - direct_diseases:
            via_tags.append("pathway_expand")
        rows.append({
            "metabolite_node_id": m,
            "score": float(len(all_d)),
            "source_ids": ";".join(sorted(all_d)),
            "via": "+".join(via_tags),
        })
    if not rows:
        return pd.DataFrame(columns=["metabolite_node_id", "score", "source_ids", "via"])
    return pd.DataFrame(rows).sort_values("score", ascending=False)
