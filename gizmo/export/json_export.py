"""
JSON export/import using NetworkX node-link format.

This is the lossless round-trip format — all Python types are preserved.
The output is compatible with D3.js force-directed graphs and
can be loaded directly into PyTorch Geometric via custom loaders.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Union

import networkx as nx
from networkx.readwrite import json_graph

from gizmo.graph.network import GizmoGraph


def write_json(mg: GizmoGraph, path: Union[str, Path], indent: int = 2) -> None:
    """Serialise a GizmoGraph to JSON node-link format."""
    import math

    def _default(obj):
        if isinstance(obj, float) and (math.isnan(obj) or math.isinf(obj)):
            return None
        return str(obj)

    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        data = json_graph.node_link_data(mg.graph)
    Path(path).write_text(json.dumps(data, indent=indent, default=_default))


def read_json(path: Union[str, Path]) -> GizmoGraph:
    """Load a GizmoGraph from JSON node-link format."""
    data = json.loads(Path(path).read_text())
    # NetworkX <3.6 uses "links"; >=3.6 uses "edges".  Normalise for compat.
    if "links" not in data and "edges" in data:
        data = dict(data, links=data["edges"])
    g = json_graph.node_link_graph(data, directed=True, multigraph=False)
    mg = GizmoGraph()
    mg._g = g
    return mg
