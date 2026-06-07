"""
Compartment utilities for GIZMO metabolite nodes.

Node ID convention:  ``CHEBI:XXXXX@compartment``
e.g. ``CHEBI:15422@cytosol``, ``CHEBI:15422@mitochondrial matrix``

Public API
----------
parse_compartment(node_id)        → compartment string or None
strip_compartment(node_id)        → base CHEBI/pubchem/metabolon ID
compartment_summary(mg)           → {compartment: count} dict
add_compartment_attributes(mg)    → writes `compartment` attr to every node in-place
same_compartment_subgraph(mg, c)  → subgraph restricted to one compartment
"""

from __future__ import annotations

from collections import Counter
from typing import Optional


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_compartment(node_id: str) -> Optional[str]:
    """
    Extract the compartment suffix from a node ID.

    Returns None if the node ID has no ``@compartment`` suffix.

    Examples
    --------
    >>> parse_compartment("CHEBI:15422@cytosol")
    'cytosol'
    >>> parse_compartment("CHEBI:15422")
    None
    >>> parse_compartment("reactome:R-HSA-70171")
    None
    """
    if "@" in node_id:
        return node_id.split("@", 1)[1]
    return None


def strip_compartment(node_id: str) -> str:
    """
    Return the base node ID without any compartment suffix.

    Examples
    --------
    >>> strip_compartment("CHEBI:15422@cytosol")
    'CHEBI:15422'
    >>> strip_compartment("CHEBI:15422")
    'CHEBI:15422'
    """
    return node_id.split("@", 1)[0]


# ---------------------------------------------------------------------------
# Graph-level utilities
# ---------------------------------------------------------------------------

def compartment_summary(mg) -> dict[str, int]:
    """
    Count metabolite nodes per compartment.

    Returns a dict ``{compartment_label: node_count}`` sorted by count desc.
    Nodes without a compartment are grouped under ``"(none)"``.
    """
    counts: Counter = Counter()
    g = mg.graph
    for nid, attrs in g.nodes(data=True):
        if attrs.get("node_type") != "metabolite":
            continue
        comp = attrs.get("compartment") or parse_compartment(nid) or "(none)"
        counts[comp] += 1
    return dict(counts.most_common())


def add_compartment_attributes(mg) -> int:
    """
    Parse ``@compartment`` suffixes from node IDs and write the result
    into each node's ``compartment`` attribute in-place.

    Only processes nodes where ``compartment`` is currently unset.
    Returns the number of nodes updated.
    """
    g       = mg.graph
    updated = 0
    for nid, attrs in g.nodes(data=True):
        if attrs.get("node_type") != "metabolite":
            continue
        if attrs.get("compartment"):
            continue
        comp = parse_compartment(nid)
        if comp:
            g.nodes[nid]["compartment"] = comp
            updated += 1
    return updated


def same_compartment_subgraph(mg, compartment: str):
    """
    Return a GizmoGraph restricted to metabolites in ``compartment``
    plus all reactions and other node types.

    Useful for compartment-aware flux or scoring analysis.
    """
    from gizmo.graph.network import GizmoGraph

    g   = mg.graph
    keep: set[str] = set()
    for nid, attrs in g.nodes(data=True):
        ntype = attrs.get("node_type", "")
        if ntype != "metabolite":
            keep.add(nid)
            continue
        node_comp = attrs.get("compartment") or parse_compartment(nid)
        if node_comp == compartment:
            keep.add(nid)

    sub = GizmoGraph()
    sub._g = g.subgraph(keep).copy()
    return sub


def collapse_compartments(mg):
    """
    Return a copy of the graph with compartment suffixes stripped from all
    metabolite node IDs — merging compartment-specific duplicates.

    Edges are redirected to the collapsed base ID.  Duplicate edges are
    kept with the highest confidence value.

    Useful for building a compartment-agnostic reaction graph.
    """
    import networkx as nx
    from gizmo.graph.network import GizmoGraph

    g    = mg.graph
    new_g = nx.DiGraph()

    id_map: dict[str, str] = {}
    for nid, attrs in g.nodes(data=True):
        base = strip_compartment(nid)
        id_map[nid] = base
        if base not in new_g:
            new_g.add_node(base, **{**attrs, "compartment": None})
        else:
            # Keep the richer attribute set (prefer nodes with chebi_id)
            existing = new_g.nodes[base]
            if not existing.get("chebi_id") and attrs.get("chebi_id"):
                new_g.nodes[base].update(attrs)
                new_g.nodes[base]["compartment"] = None

    for u, v, eattrs in g.edges(data=True):
        bu, bv = id_map[u], id_map[v]
        if bu == bv:
            continue
        if new_g.has_edge(bu, bv):
            # Keep max confidence
            existing_conf = new_g.edges[bu, bv].get("confidence", 0)
            new_conf      = eattrs.get("confidence", 1.0)
            if new_conf > existing_conf:
                new_g.edges[bu, bv].update(eattrs)
        else:
            new_g.add_edge(bu, bv, **eattrs)

    out = GizmoGraph()
    out._g = new_g
    return out
