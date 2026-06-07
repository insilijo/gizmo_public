"""Branch-point scoring for metabolites in the substrate graph.

A metabolite is a high-information branch point if measuring it
simultaneously gates many pathways and many reactions. This is the
*opposite* of the hub-penalty used in reaction scoring (where high-
degree nodes are downweighted to avoid spurious signal aggregation):
for panel-design purposes, hubness *is* the signal.

Score columns
-------------
- n_reactions  : count of distinct reactions touching the metabolite
- n_pathways   : count of distinct Reactome pathways covering those
                 reactions
- pagerank     : substrate-graph PageRank (carries hub-bias)
- degree       : undirected degree in the substrate subgraph
- branch_score : log1p(n_pathways) * pagerank  (composite)

The composite leans on n_pathways because that is the load-bearing
interpretation — "this single measurement informs N pathways" — and
multiplies by PageRank so that two metabolites tied on n_pathways
break by network-topological centrality.
"""
from __future__ import annotations

from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd


def compute_metabolite_branch_points(mg) -> pd.DataFrame:
    """Walk the GizmoGraph once; emit a DataFrame indexed by metabolite node_id.

    The substrate subgraph is taken to be the bipartite (metabolite, reaction)
    induced subgraph. PageRank runs on the undirected projection so that
    forward/reverse reactions don't bias toward downstream products.
    """
    g = mg.graph if hasattr(mg, "graph") else mg

    metab_nodes = [
        n for n, d in g.nodes(data=True) if d.get("node_type") == "metabolite"
    ]
    reaction_pathways: dict[str, set[str]] = {}
    for n, d in g.nodes(data=True):
        if d.get("node_type") != "reaction":
            continue
        reaction_pathways[n] = set(d.get("pathways") or [])

    rxns_by_metab: dict[str, set[str]] = {m: set() for m in metab_nodes}
    for u, v, ed in g.edges(data=True):
        role = (ed.get("role") or ed.get("edge_type") or "").lower()
        if role not in ("substrate", "product", "modifier"):
            continue
        u_t = g.nodes.get(u, {}).get("node_type")
        v_t = g.nodes.get(v, {}).get("node_type")
        if u_t == "metabolite" and v_t == "reaction":
            rxns_by_metab[u].add(v)
        elif v_t == "metabolite" and u_t == "reaction":
            rxns_by_metab[v].add(u)

    bipartite_edges = []
    for m, rxns in rxns_by_metab.items():
        for r in rxns:
            bipartite_edges.append((m, r))
    sub = nx.Graph()
    sub.add_nodes_from(metab_nodes)
    sub.add_edges_from(bipartite_edges)

    pr = nx.pagerank(sub, alpha=0.85) if sub.number_of_edges() else {}

    rows = []
    for m in metab_nodes:
        rxns = rxns_by_metab.get(m, set())
        pways: set[str] = set()
        for r in rxns:
            pways |= reaction_pathways.get(r, set())
        attrs = g.nodes[m]
        rows.append({
            "metabolite_node_id": m,
            "name": attrs.get("name", m),
            "pubchem_title": attrs.get("pubchem_title"),
            "pubchem_cid": attrs.get("pubchem_cid"),
            "hmdb_id": attrs.get("hmdb_id"),
            "chebi_id": attrs.get("chebi_id"),
            "inchikey": attrs.get("inchikey"),
            "is_currency": bool(attrs.get("is_currency", False)),
            "n_reactions": len(rxns),
            "n_pathways": len(pways),
            "degree": sub.degree(m) if m in sub else 0,
            "pagerank": float(pr.get(m, 0.0)),
        })
    df = pd.DataFrame(rows)
    df["branch_score"] = np.log1p(df["n_pathways"]) * df["pagerank"]
    df = df.sort_values("branch_score", ascending=False).reset_index(drop=True)
    return df


def save_branch_points(df: pd.DataFrame, path: Path | str) -> None:
    df.to_parquet(Path(path), index=False)


def load_branch_points(path: Path | str) -> pd.DataFrame | None:
    p = Path(path)
    if not p.exists():
        return None
    return pd.read_parquet(p)
