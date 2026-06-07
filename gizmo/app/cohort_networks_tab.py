"""GIZMO Dash app — Cohort Networks tab.

Renders the prebuilt benchmark cohort cytoscape exports (one per
cohort × design). Loads `data/processed/.../cytoscape_web/{cohort}_{design}.json`
on demand via dropdown.

Each network shows the top consensus reactions from per-patient kernel
analysis, plus their 1-hop biochemistry neighbors, color-coded by
anchor tier (anchor / confirmed / bridging / novel / neighbor).

Wire in from dash_app.py::

    from gizmo.app.cohort_networks_tab import build_tab, register_callbacks
    cohort_tab = build_tab()
    # Add cohort_tab to dcc.Tabs children
    register_callbacks(app)
"""
from __future__ import annotations

import json
from pathlib import Path

import dash_cytoscape as cyto
from dash import Input, Output, dcc, html

# Default location — relative to repo root or set via env var
import os
_RESULTS_DIR = Path(os.environ.get(
    "GIZMO_COHORT_NETWORKS_DIR",
    str(Path(__file__).resolve().parent.parent.parent / "benchmarks/results/cytoscape_web")
))


# Color scheme by anchor tier — red for canonical, blue for distant biology
TIER_COLORS = {
    "anchor":     "#27ae60",   # green — direct anchor reaction
    "confirmed":  "#52c878",   # light green — 1-hop neighbor of anchor
    "bridging":   "#f39c12",   # orange — 2-3 hop reach
    "novel":      "#e74c3c",   # red — outside curated knowledge
    "neighbor":   "#bdc3c7",   # gray — supporting neighbor (gene/metab)
}

NODE_TYPE_SHAPES = {
    "reaction":   "round-rectangle",
    "metabolite": "ellipse",
    "gene":       "diamond",
}


def _list_available_networks() -> list[dict]:
    """Scan the cytoscape_web dir for available networks."""
    out = []
    if not _RESULTS_DIR.exists():
        return out
    for f in sorted(_RESULTS_DIR.glob("*.json")):
        cohort_design = f.stem
        out.append({
            "label": cohort_design.replace("_", " / ", 1).replace("_", " ", 1),
            "value": cohort_design,
        })
    return out


def _load_network(cohort_design: str) -> tuple[list, dict]:
    """Load a cytoscape.js JSON file. Returns (elements, metadata)."""
    f = _RESULTS_DIR / f"{cohort_design}.json"
    if not f.exists():
        return [], {}
    data = json.loads(f.read_text())
    return data.get("elements", []), data.get("metadata", {})


def build_tab():
    """Build the Cohort Networks dcc.Tab."""
    options = _list_available_networks()
    default_value = options[0]["value"] if options else None

    legend = html.Div([
        html.Span("Anchor tier: "),
        *[html.Span(
            tier.upper(),
            style={
                "background": col, "color": "white",
                "padding": "2px 8px", "marginRight": "4px",
                "borderRadius": "3px", "fontSize": "11px",
            }
        ) for tier, col in TIER_COLORS.items() if tier != "neighbor"],
    ], style={"padding": "8px 0", "fontSize": "12px"})

    body = html.Div([
        html.Div([
            html.Label("Cohort / Design:",
                       style={"fontWeight": "600", "marginRight": "10px"}),
            dcc.Dropdown(
                id="cohort-network-select",
                options=options,
                value=default_value,
                clearable=False,
                style={"width": "320px", "display": "inline-block"},
            ),
            html.Span(id="cohort-network-meta",
                       style={"marginLeft": "20px", "fontSize": "13px",
                              "color": "#555"}),
        ], style={"padding": "10px"}),
        legend,
        cyto.Cytoscape(
            id="cohort-cytoscape",
            elements=[],
            layout={"name": "cose", "animate": False,
                    "nodeRepulsion": 12000, "idealEdgeLength": 90},
            style={"width": "100%", "height": "70vh",
                   "border": "1px solid #ddd", "borderRadius": "4px"},
            stylesheet=[
                {"selector": "node", "style": {
                    "label": "data(label)",
                    "font-size": "10px",
                    "text-valign": "center", "text-halign": "center",
                    "width": 40, "height": 40,
                    "background-color": (
                        "data(_color)"   # set per-node via callback
                    ),
                    "shape": "data(_shape)",
                    "color": "#1f2937",
                    "text-wrap": "wrap", "text-max-width": "80px",
                    "border-width": 1, "border-color": "#94a3b8",
                }},
                {"selector": "edge", "style": {
                    "curve-style": "bezier",
                    "target-arrow-shape": "triangle",
                    "width": 1.2,
                    "line-color": "#94a3b8",
                    "target-arrow-color": "#94a3b8",
                    "opacity": 0.7,
                }},
            ],
        ),
        html.Div(id="cohort-network-info",
                 style={"padding": "10px", "fontSize": "12px",
                         "color": "#444"}),
    ])

    return dcc.Tab(label="Cohort Networks", value="cohort_networks",
                   children=[body])


def register_callbacks(app):
    """Register the callbacks for the Cohort Networks tab."""

    @app.callback(
        Output("cohort-cytoscape", "elements"),
        Output("cohort-network-meta", "children"),
        Output("cohort-network-info", "children"),
        Input("cohort-network-select", "value"),
    )
    def update_network(cohort_design):
        if not cohort_design:
            return [], "no network selected", ""
        elements, meta = _load_network(cohort_design)
        # Annotate each node with color + shape based on tier + type
        for el in elements:
            data = el.get("data", {})
            if "source" in data:   # edge
                continue
            tier = data.get("anchor_tier", "neighbor")
            data["_color"] = TIER_COLORS.get(tier, "#bdc3c7")
            data["_shape"] = NODE_TYPE_SHAPES.get(data.get("type"), "ellipse")

        n_nodes = sum(1 for e in elements if "source" not in e.get("data", {}))
        n_edges = sum(1 for e in elements if "source" in e.get("data", {}))

        meta_text = (f"{n_nodes} nodes, {n_edges} edges  | "
                     f"{meta.get('n_consensus_reactions_used', '?')} consensus reactions  | "
                     f"{meta.get('anchor_count', '?')} curated anchor reactions for cohort")

        # Tier breakdown for info area
        from collections import Counter
        tier_counts = Counter(
            el["data"].get("anchor_tier", "neighbor")
            for el in elements if "source" not in el.get("data", {})
        )
        info = html.Div([
            html.B("Anchor tier breakdown: "),
            *[html.Span(f"{t}={c}  ", style={"marginRight": "8px"})
              for t, c in sorted(tier_counts.items(),
                                  key=lambda kv: -kv[1])],
        ])

        return elements, meta_text, info
