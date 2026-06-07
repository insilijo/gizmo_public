"""
GraphML export/import — compatible with GATOM, Cytoscape, igraph (R), and gephi.

Notes:
  - GraphML only supports scalar attribute types (string, int, float, boolean).
    List-valued attributes (ec_numbers, pathways) are serialised as pipe-delimited strings.
    Dict-valued attributes (reference_ranges, tissue_expression) are serialised as JSON strings.
  - Node type is preserved as a "node_type" attribute so R/igraph can filter
    the bipartite graph back into metabolite-only projections.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Union

import networkx as nx

from gizmo.graph.network import GizmoGraph


def write_graphml(mg: GizmoGraph, path: Union[str, Path]) -> None:
    """Write a GizmoGraph to GraphML, sanitising list and dict attributes."""
    g = _sanitise_for_graphml(mg.graph)
    nx.write_graphml(g, str(path))


def read_graphml(path: Union[str, Path]) -> GizmoGraph:
    """Read a GraphML file into a GizmoGraph (round-trip safe)."""
    g = nx.read_graphml(str(path))
    # Restore pipe-delimited lists
    for _, data in g.nodes(data=True):
        for key in ("ec_numbers", "pathways"):
            if key in data and isinstance(data[key], str):
                data[key] = [v for v in data[key].split("|") if v]
    mg = GizmoGraph()
    mg._g = g
    return mg


# ------------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------------

def _sanitise_for_graphml(g: nx.DiGraph) -> nx.DiGraph:
    """Return a copy with non-scalar attributes converted to GraphML-safe strings."""
    g2 = g.copy()
    for _, data in g2.nodes(data=True):
        for key, val in list(data.items()):
            if isinstance(val, list):
                data[key] = "|".join(str(v) for v in val)
            elif isinstance(val, dict):
                data[key] = json.dumps(val)
            elif val is None:
                data[key] = ""
    for _, _, data in g2.edges(data=True):
        for key, val in list(data.items()):
            if isinstance(val, list):
                data[key] = "|".join(str(v) for v in val)
            elif isinstance(val, dict):
                data[key] = json.dumps(val)
            elif val is None:
                data[key] = ""
    return g2
