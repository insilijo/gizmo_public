"""
GIZMO Dash explorer — interactive graph visualisation + QC report + Metabolon curation.

Quick start (Python)::

    from gizmo.app.dash_app import create_app
    app = create_app(mg)
    app.run(debug=True, port=8050)

CLI::

    gizmo app --graph data/processed/gizmo_full.json [--port 8050]
    gizmo app --graph gizmo_full.json --metabolon-csv metabolon.csv

Production (gunicorn)::

    GIZMO_GRAPH=data/processed/gizmo_full.json \\
        gunicorn gizmo.app.wsgi:server -b 0.0.0.0:8050 -w 2 --timeout 120
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

import dash_cytoscape as cyto
from dash import (
    ALL,
    Dash,
    Input,
    Output,
    State,
    callback_context,
    dash_table,
    dcc,
    html,
    no_update,
)

log = logging.getLogger(__name__)

cyto.load_extra_layouts()

# ---------------------------------------------------------------------------
# Colour constants
# ---------------------------------------------------------------------------

NODE_COLORS: dict[str, str] = {
    "metabolite": "#4C8BE0",
    "reaction":   "#50C878",
    "disease":    "#E05252",
    "gene":       "#F5A623",
    "unknown":    "#AAAAAA",
}

EDGE_COLORS: dict[str, str] = {
    "substrate":            "#E07B4C",
    "product":              "#4CE07B",
    "modifier":             "#B04CE0",
    "protein_interaction":  "#F5A623",
    "gene_associated":      "#E05252",
    "gene_reaction":        "#F5A623",
    "biomarker":            "#9B59B6",
    "causal":               "#E74C3C",
    "pathway_associated":   "#3498DB",
    "default":              "#CCCCCC",
}

CURRENCY_COLOR  = "#E05252"
MATCHED_COLOR   = "#2196F3"
UNMATCHED_COLOR = "#9E9E9E"

_PLATFORM_PALETTE = [
    "#e41a1c", "#377eb8", "#4daf4a", "#984ea3",
    "#ff7f00", "#a65628", "#f781bf", "#999999",
]

# ---------------------------------------------------------------------------
# Cytoscape stylesheet
# ---------------------------------------------------------------------------

CYTO_STYLESHEET = [
    {
        "selector": "node",
        "style": {
            "content": "data(label)",
            "font-size": "9px",
            "text-valign": "center",
            "text-halign": "center",
            "background-color": "data(color)",
            "width": "data(size)",
            "height": "data(size)",
            "color": "#111",
            "text-outline-color": "#fff",
            "text-outline-width": "1.5px",
        },
    },
    {
        "selector": "edge",
        "style": {
            "curve-style": "bezier",
            "target-arrow-shape": "triangle",
            "arrow-scale": 1.2,
            "line-color": "data(color)",
            "target-arrow-color": "data(color)",
            "width": 1.5,
            "opacity": 0.75,
        },
    },
    {
        "selector": "node:selected",
        "style": {"border-width": 3, "border-color": "#FFD700", "border-opacity": 1},
    },
    {
        "selector": "edge:selected",
        "style": {"width": 3, "line-color": "#FFD700"},
    },
]

# ---------------------------------------------------------------------------
# Node helpers
# ---------------------------------------------------------------------------

def _node_label(nid: str, attrs: dict) -> str:
    name = (
        attrs.get("name")
        or attrs.get("metabolon_name")
        or attrs.get("symbol")
        or nid
    )
    return str(name)[:30]


_SCORE_GRADIENT = [
    "#313695", "#4575b4", "#74add1", "#abd9e9",
    "#fee090", "#fdae61", "#f46d43", "#d73027",
]


def _score_to_color(v: float | None) -> str:
    """Map a [0, 1] score to a diverging colour ramp (blue → red)."""
    if v is None:
        return "#AAAAAA"
    idx = min(int(max(v, 0.0) * len(_SCORE_GRADIENT)), len(_SCORE_GRADIENT) - 1)
    return _SCORE_GRADIENT[idx]


def _color_for_node(
    nid: str,
    attrs: dict,
    color_by: str,
    platform_index: dict,
    score_overlay: dict[str, float] | None = None,
) -> str:
    ntype = attrs.get("node_type", "unknown")
    if color_by == "node_type":
        return NODE_COLORS.get(ntype, NODE_COLORS["unknown"])
    if color_by == "is_currency":
        return CURRENCY_COLOR if attrs.get("is_currency") else NODE_COLORS.get(ntype, "#4C8BE0")
    if color_by == "has_chebi":
        if ntype != "metabolite":
            return NODE_COLORS.get(ntype, "#AAA")
        return MATCHED_COLOR if attrs.get("chebi_id") else UNMATCHED_COLOR
    if color_by == "has_metabolon":
        if ntype != "metabolite":
            return NODE_COLORS.get(ntype, "#AAA")
        return "#9C27B0" if attrs.get("metabolon_name") else NODE_COLORS["metabolite"]
    if color_by == "platform":
        plat = attrs.get("platform") or ""
        idx  = platform_index.get(plat, 0)
        return _PLATFORM_PALETTE[idx % len(_PLATFORM_PALETTE)] if plat else "#AAA"
    if color_by in ("druggability", "perturbability"):
        return _score_to_color(score_overlay.get(nid) if score_overlay else None)
    return NODE_COLORS.get(ntype, "#AAA")


def _node_size(attrs: dict) -> int:
    return {"reaction": 18, "disease": 30, "gene": 22}.get(
        attrs.get("node_type", ""), 24
    )


def build_elements(
    g,
    node_ids: list[str],
    color_by: str,
    platform_index: dict,
    edge_roles: set[str] | None = None,
    score_overlay: dict[str, float] | None = None,
) -> list[dict]:
    """Convert a node ID list into Cytoscape element dicts (nodes + edges)."""
    node_set = set(node_ids)
    elements: list[dict] = []

    for nid in node_ids:
        attrs = g.nodes[nid]
        elements.append({
            "data": {
                "id":             nid,
                "label":          _node_label(nid, attrs),
                "color":          _color_for_node(nid, attrs, color_by, platform_index, score_overlay),
                "size":           _node_size(attrs),
                "node_type":      attrs.get("node_type", ""),
                "chebi_id":       attrs.get("chebi_id") or "",
                "hmdb_id":        attrs.get("hmdb_id") or "",
                "metabolon_name": attrs.get("metabolon_name") or "",
                "is_currency":    str(attrs.get("is_currency", False)),
                "platform":       attrs.get("platform") or "",
                "mass":           str(attrs.get("mass") or ""),
                "inchikey":       attrs.get("inchikey") or "",
                "formula":        attrs.get("formula") or "",
                "ec_numbers":     ", ".join(attrs.get("ec_numbers") or []),
                "gene_symbols":   ", ".join((attrs.get("gene_symbols") or [])[:8]),
                "pathways":       ", ".join((attrs.get("pathways") or [])[:6]),
                "reversible":     str(attrs.get("reversible", "")),
                "symbol":         attrs.get("symbol") or "",
                "mondo_id":       attrs.get("mondo_id") or "",
                "is_rare":        str(attrs.get("is_rare", "")),
            }
        })

    for src, tgt, eattrs in g.edges(data=True):
        if src not in node_set or tgt not in node_set:
            continue
        # PPI edges store edge_type="protein_interaction" (no `role` key)
        role = eattrs.get("role") or eattrs.get("edge_type") or "default"
        if edge_roles is not None and role not in edge_roles:
            continue
        edata: dict = {
            "source":        src,
            "target":        tgt,
            "role":          role,
            "stoichiometry": str(eattrs.get("stoichiometry") or ""),
            "color":         EDGE_COLORS.get(role, EDGE_COLORS["default"]),
        }
        if role == "protein_interaction":
            edata["combined_score"] = str(eattrs.get("combined_score", ""))
            edata["experimental"]   = str(eattrs.get("experimental", ""))
            edata["coexpression"]   = str(eattrs.get("coexpression", ""))
            edata["database"]       = str(eattrs.get("database", ""))
            edata["textmining"]     = str(eattrs.get("textmining", ""))
        elements.append({"data": edata})

    return elements

# ---------------------------------------------------------------------------
# Graph filtering helpers
# ---------------------------------------------------------------------------

def _search_nodes(g, query: str, max_results: int = 80) -> list[str]:
    """Return node IDs whose label or raw ID contains `query` (case-insensitive)."""
    q = query.lower().strip()
    if not q:
        return []
    matches: list[str] = []
    for nid, attrs in g.nodes(data=True):
        label = (_node_label(nid, attrs) or "").lower()
        if q in label or q in nid.lower():
            matches.append(nid)
            if len(matches) >= max_results:
                break
    return matches


def _neighborhood(g, nid: str, hops: int = 1) -> set[str]:
    """All nodes within `hops` hops of `nid` (predecessors + successors)."""
    visited = {nid}
    frontier = {nid}
    for _ in range(hops):
        nxt: set[str] = set()
        for n in frontier:
            nxt.update(g.predecessors(n))
            nxt.update(g.successors(n))
        frontier = nxt - visited
        visited |= frontier
    return visited

# ---------------------------------------------------------------------------
# Pathway name index
# ---------------------------------------------------------------------------

def build_pathway_name_index(reactome_cache_dir: Path) -> dict[str, str]:
    """Scan cached `events_*.json` files to build stId → displayName."""
    names: dict[str, str] = {}
    for f in reactome_cache_dir.glob("events_*.json"):
        try:
            events = json.loads(f.read_text())
        except Exception:
            continue
        for ev in events:
            s = ev.get("stId") or ev.get("dbId")
            n = ev.get("displayName") or ev.get("name")
            if s and n:
                names[str(s)] = str(n)
    return names

# ---------------------------------------------------------------------------
# Layout micro-helpers
# ---------------------------------------------------------------------------

_BTN_BASE: dict[str, Any] = {
    "fontSize": "12px", "padding": "4px 10px",
    "border": "none", "borderRadius": "4px", "cursor": "pointer",
}


def _lbl(text: str) -> html.Div:
    return html.Div(text, style={
        "fontSize": "11px", "fontWeight": "bold", "marginBottom": "3px",
    })


def _section(title: str, *children) -> html.Div:
    return html.Div([
        html.Div(title, style={
            "fontSize": "10px", "fontWeight": "bold", "color": "#666",
            "textTransform": "uppercase", "letterSpacing": "1px",
            "marginTop": "10px", "marginBottom": "4px",
        }),
        *children,
    ])


def _prop_row(key: str, val: Any) -> html.Tr:
    return html.Tr([
        html.Td(key, style={
            "fontWeight": "bold", "paddingRight": "8px", "whiteSpace": "nowrap",
            "verticalAlign": "top", "fontSize": "11px", "color": "#555",
            "paddingBottom": "3px",
        }),
        html.Td(
            str(val) if val is not None else "—",
            style={"fontSize": "11px", "wordBreak": "break-all", "paddingBottom": "3px"},
        ),
    ])


def _stat_card(label: str, value: Any, color: str = "#1565C0") -> html.Div:
    return html.Div(style={
        "flex": "1", "minWidth": "80px", "background": "#fff",
        "border": f"2px solid {color}", "borderRadius": "6px",
        "padding": "5px 8px", "textAlign": "center",
    }, children=[
        html.Div(str(value), style={"fontSize": "18px", "fontWeight": "bold", "color": color}),
        html.Div(label, style={"fontSize": "10px", "color": "#666", "marginTop": "1px"}),
    ])


def _qc_metric(label: str, value: str, ok: bool | None = None) -> html.Div:
    indicator = ""
    indicator_color = "#555"
    if ok is True:
        indicator, indicator_color = "✓", "#2e7d32"
    elif ok is False:
        indicator, indicator_color = "⚠", "#c62828"

    return html.Div(style={
        "display": "flex", "justifyContent": "space-between", "alignItems": "center",
        "padding": "5px 0", "borderBottom": "1px solid #eee",
    }, children=[
        html.Span(label, style={"fontSize": "12px", "color": "#333"}),
        html.Span([
            html.Span(value, style={"fontSize": "12px", "fontWeight": "bold"}),
            html.Span(f" {indicator}", style={"color": indicator_color, "marginLeft": "4px"}),
        ]),
    ])

# ---------------------------------------------------------------------------
# App factory
# ---------------------------------------------------------------------------

def create_app(
    mg,
    met_loader=None,
    curator=None,
    reactome_cache_dir: str | Path = "data/raw/reactome",
    max_nodes: int = 500,
    action_scores=None,
) -> Dash:
    """
    Create and return the GIZMO Dash explorer.

    Parameters
    ----------
    mg : GizmoGraph
    met_loader : MetabolonLoader, optional   — enables curation tab
    curator : MetabolonCurator, optional     — enables curation tab
    reactome_cache_dir : directory with cached Reactome JSON (for pathway names)
    max_nodes : performance cap on simultaneously rendered nodes
    action_scores : list[ActionabilityScore], optional
        Pre-computed actionability scores; enables druggability/perturbability
        color-by modes in the graph explorer.
    """
    g = mg.graph

    # ------------------------------------------------------------------ #
    # Pre-compute indexes and static data                                  #
    # ------------------------------------------------------------------ #

    cache_dir = Path(reactome_cache_dir)
    pathway_names = build_pathway_name_index(cache_dir) if cache_dir.exists() else {}

    # Pathway → node set (reactions + their direct neighbours)
    all_pw_stids: set[str] = set()
    pathway_to_nodes: dict[str, set[str]] = {}
    for nid in g.nodes():
        pws = g.nodes[nid].get("pathways") or []
        all_pw_stids.update(pws)
        for pw in pws:
            pathway_to_nodes.setdefault(pw, set()).add(nid)
            for nbr in list(g.predecessors(nid)) + list(g.successors(nid)):
                pathway_to_nodes[pw].add(nbr)

    pw_options = sorted(
        [{"label": pathway_names.get(s, s), "value": s} for s in all_pw_stids],
        key=lambda x: x["label"],
    )

    platforms = sorted({
        g.nodes[n].get("platform") or ""
        for n in g.nodes()
        if g.nodes[n].get("platform")
    })
    platform_index = {p: i for i, p in enumerate(platforms)}

    # Graph-level stats for the stats bar
    n_mets      = sum(1 for _, d in g.nodes(data=True) if d.get("node_type") == "metabolite")
    n_rxns      = sum(1 for _, d in g.nodes(data=True) if d.get("node_type") == "reaction")
    n_disease   = sum(1 for _, d in g.nodes(data=True) if d.get("node_type") == "disease")
    n_gene      = sum(1 for _, d in g.nodes(data=True) if d.get("node_type") == "gene")
    n_chebi     = sum(1 for _, d in g.nodes(data=True)
                      if d.get("chebi_id") and d.get("node_type") == "metabolite")
    n_metabolon = sum(1 for _, d in g.nodes(data=True) if d.get("metabolon_name"))

    # Actionability score overlays: {node_id -> normalised [0,1] score}
    # gene nodes → druggability; reaction nodes → perturbability
    _drug_overlay:  dict[str, float] = {}
    _pert_overlay:  dict[str, float] = {}
    if action_scores:
        drug_vals = [a.druggability_score for a in action_scores if a.druggability_score is not None]
        pert_vals = [a.perturbability_score for a in action_scores if a.perturbability_score is not None]
        max_drug = max(drug_vals, default=1.0) or 1.0
        max_pert = max(pert_vals, default=1.0) or 1.0
        for a in action_scores:
            if a.gene_id and a.druggability_score is not None:
                _drug_overlay[a.gene_id]     = a.druggability_score / max_drug
            if a.reaction_id and a.perturbability_score is not None:
                _pert_overlay[a.reaction_id] = a.perturbability_score / max_pert

    _has_actionability = bool(_drug_overlay or _pert_overlay)

    # QC report (computed once)
    try:
        from gizmo.analysis.qc import assess_readiness
        _qc = assess_readiness(mg)
        log.info("QC report computed.")
    except Exception as exc:
        log.warning("QC report failed: %s", exc)
        _qc = None

    # Unmatched records for curation table
    def _unmatched_records() -> list[dict]:
        if curator is None:
            return []
        df = curator.unmatched
        return df[["BIOCHEMICAL", "PATHWAY", "INCHIKEY", "PUBCHEM", "MASS"]].fillna("").to_dict("records")

    # ------------------------------------------------------------------ #
    # Dash application                                                     #
    # ------------------------------------------------------------------ #

    app = Dash(
        __name__,
        title="GIZMO Explorer",
        suppress_callback_exceptions=True,
        meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}],
    )

    # -- shared style dicts --
    _SIDEBAR = {
        "width": "280px", "minWidth": "280px", "padding": "10px 14px",
        "background": "#f8f9fa", "borderRight": "1px solid #ddd",
        "overflowY": "auto", "fontSize": "13px",
    }
    _PROPS = {
        "width": "270px", "minWidth": "270px", "padding": "10px 14px",
        "background": "#f8f9fa", "borderLeft": "1px solid #ddd",
        "overflowY": "auto", "fontSize": "13px",
    }

    # ------------------------------------------------------------------ #
    # Header + stats bar                                                   #
    # ------------------------------------------------------------------ #

    header = html.Div(style={
        "background": "#1565C0", "color": "#fff",
        "padding": "6px 16px", "display": "flex", "alignItems": "center", "gap": "16px",
    }, children=[
        html.Span("GIZMO", style={"fontWeight": "bold", "fontSize": "20px"}),
        html.Span("Graph-Integrated Zone of Metabolite Operations",
                  style={"fontSize": "13px", "opacity": 0.85}),
    ])

    stats_bar = html.Div(style={
        "display": "flex", "gap": "6px", "padding": "6px 16px",
        "background": "#e3f2fd", "borderBottom": "1px solid #bbdefb", "flexWrap": "wrap",
    }, children=[
        _stat_card("Nodes",      g.number_of_nodes(), "#1565C0"),
        _stat_card("Edges",      g.number_of_edges(), "#1565C0"),
        _stat_card("Metabolites", n_mets,             "#4C8BE0"),
        _stat_card("Reactions",   n_rxns,             "#50C878"),
        _stat_card("Diseases",    n_disease,          "#E05252"),
        _stat_card("Genes",       n_gene,             "#F5A623"),
        _stat_card("w/ ChEBI",    n_chebi,            "#2196F3"),
        _stat_card("Metabolon",   n_metabolon,        "#9C27B0"),
    ])

    # ================================================================== #
    # Tab 1 — Graph Explorer                                               #
    # ================================================================== #

    graph_tab = dcc.Tab(label="Graph Explorer", value="graph", children=[
        html.Div(style={"display": "flex", "height": "calc(100vh - 120px)"}, children=[

            # --- Left sidebar ---
            html.Div(style=_SIDEBAR, children=[

                _section("Search nodes",
                    html.Div(style={"display": "flex", "gap": "4px"}, children=[
                        dcc.Input(
                            id="node-search", type="text",
                            placeholder="name or ID…", debounce=True,
                            style={
                                "flex": 1, "fontSize": "12px", "padding": "4px 6px",
                                "border": "1px solid #ccc", "borderRadius": "4px",
                            },
                        ),
                        html.Button("×", id="node-search-clear", n_clicks=0,
                                    style={**_BTN_BASE, "background": "#e0e0e0", "color": "#333"}),
                    ]),
                ),

                _section("Pathway filter",
                    dcc.Dropdown(
                        id="pw-selector",
                        options=pw_options,
                        placeholder="Select pathway(s)…",
                        multi=True,
                        style={"fontSize": "11px"},
                    ),
                ),

                html.Div(style={"marginTop": "6px"}, children=[
                    dcc.Checklist(
                        id="show-all",
                        options=[{"label": f"  Show all nodes (≤ {max_nodes})", "value": "all"}],
                        value=[],
                        style={"fontSize": "12px"},
                    ),
                ]),

                _section("Node types",
                    dcc.Checklist(
                        id="type-filter",
                        options=[
                            {"label": " Metabolites", "value": "metabolite"},
                            {"label": " Reactions",   "value": "reaction"},
                            {"label": " Diseases",    "value": "disease"},
                            {"label": " Genes",       "value": "gene"},
                        ],
                        value=["metabolite", "reaction", "disease", "gene"],
                        labelStyle={"display": "block"},
                        style={"fontSize": "12px"},
                    ),
                ),

                _section("Edge types",
                    dcc.Checklist(
                        id="edge-filter",
                        options=[
                            {"label": " Substrate",    "value": "substrate"},
                            {"label": " Product",      "value": "product"},
                            {"label": " Modifier",     "value": "modifier"},
                            {"label": " PPI (STRING)", "value": "protein_interaction"},
                            {"label": " Other",        "value": "default"},
                        ],
                        value=["substrate", "product", "modifier", "protein_interaction", "default"],
                        labelStyle={"display": "block"},
                        style={"fontSize": "12px"},
                    ),
                ),

                _section("Colour by",
                    dcc.RadioItems(
                        id="color-by",
                        options=(
                            [
                                {"label": " Node type",       "value": "node_type"},
                                {"label": " Currency status",  "value": "is_currency"},
                                {"label": " Has ChEBI ID",    "value": "has_chebi"},
                                {"label": " Has Metabolon",   "value": "has_metabolon"},
                                {"label": " Platform",        "value": "platform"},
                            ] + ([
                                {"label": " Druggability (ChEMBL)", "value": "druggability"},
                                {"label": " Perturbability",        "value": "perturbability"},
                            ] if _has_actionability else [])
                        ),
                        value="node_type",
                        labelStyle={"display": "block"},
                        style={"fontSize": "12px"},
                    ),
                ),

                _section("Layout",
                    dcc.Dropdown(
                        id="layout-selector",
                        options=[
                            {"label": "Cose-Bilkent (recommended)", "value": "cose-bilkent"},
                            {"label": "Cose (force-directed)",      "value": "cose"},
                            {"label": "Concentric",                 "value": "concentric"},
                            {"label": "Breadthfirst",               "value": "breadthfirst"},
                            {"label": "Circle",                     "value": "circle"},
                            {"label": "Grid",                       "value": "grid"},
                        ],
                        value="cose-bilkent",
                        clearable=False,
                        style={"fontSize": "12px"},
                    ),
                ),

                html.Div(id="node-count-badge", style={
                    "marginTop": "10px", "padding": "4px 8px",
                    "background": "#e9ecef", "borderRadius": "4px",
                    "fontSize": "11px", "color": "#495057",
                }),

                html.Hr(style={"margin": "8px 0"}),
                html.Div("Legend", style={"fontWeight": "bold", "fontSize": "11px", "marginBottom": "4px"}),
                html.Div(id="legend-div"),
            ]),

            # --- Centre canvas ---
            html.Div(style={"flex": "1", "position": "relative"}, children=[
                html.Div(
                    id="graph-empty-msg",
                    children=[
                        html.Div("GIZMO Explorer",
                                 style={"fontSize": "24px", "fontWeight": "bold",
                                        "color": "#1565C0", "marginBottom": "8px"}),
                        html.Div(
                            "Search for a node · select a pathway · or tick 'Show all' to begin.",
                            style={"fontSize": "14px", "color": "#888"},
                        ),
                    ],
                    style={
                        "position": "absolute", "top": "50%", "left": "50%",
                        "transform": "translate(-50%,-50%)", "textAlign": "center", "zIndex": 10,
                    },
                ),
                cyto.Cytoscape(
                    id="cytoscape",
                    elements=[],
                    layout={"name": "cose-bilkent", "animate": False, "randomize": False},
                    style={"width": "100%", "height": "100%"},
                    stylesheet=CYTO_STYLESHEET,
                    minZoom=0.05,
                    maxZoom=8,
                ),
            ]),

            # --- Right properties panel ---
            html.Div(style=_PROPS, children=[
                html.Div("Properties", style={"fontWeight": "bold", "fontSize": "13px", "marginBottom": "6px"}),
                html.Div(id="props-panel",
                         children=html.Div("Click a node or edge.", style={"color": "#999", "fontSize": "12px"})),

                html.Hr(style={"margin": "8px 0"}),
                html.Div("Explore neighbourhood", style={"fontWeight": "bold", "fontSize": "12px", "marginBottom": "4px"}),
                html.Div(style={"display": "flex", "gap": "6px", "alignItems": "center"}, children=[
                    html.Button("Show neighbors", id="btn-neighbors", n_clicks=0,
                                style={**_BTN_BASE, "background": "#1565C0", "color": "#fff"}),
                    dcc.Dropdown(
                        id="hop-count",
                        options=[{"label": f"{n} hop{'s' if n > 1 else ''}", "value": n} for n in range(1, 4)],
                        value=1,
                        clearable=False,
                        style={"width": "90px", "fontSize": "12px"},
                    ),
                ]),
                html.Div(id="explore-status",
                         style={"fontSize": "11px", "color": "#666", "marginTop": "4px"}),
            ]),
        ]),
    ])

    # ================================================================== #
    # Tab 2 — QC Report                                                    #
    # ================================================================== #

    def _qc_tab_content():
        if _qc is None:
            return html.Div("QC report unavailable.", style={"padding": "20px", "color": "#888"})

        def pct(f: float) -> str:
            return f"{f * 100:.1f}%"

        col_style = {"flex": "1", "minWidth": "280px", "padding": "0 16px"}
        section_head = lambda t: html.Div(t, style={
            "fontSize": "13px", "fontWeight": "bold", "color": "#1565C0",
            "borderBottom": "2px solid #1565C0", "marginBottom": "8px", "paddingBottom": "2px",
            "marginTop": "20px",
        })

        left = html.Div(style=col_style, children=[
            section_head("Graph Composition"),
            _qc_metric("Metabolite nodes",  str(_qc.n_metabolites), _qc.n_metabolites > 0),
            _qc_metric("Reaction nodes",    str(_qc.n_reactions),   _qc.n_reactions > 0),
            _qc_metric("Disease nodes",     str(_qc.n_diseases)),
            _qc_metric("Gene nodes",        str(_qc.n_genes)),
            _qc_metric("Edges",             str(_qc.n_edges),       _qc.n_edges > 0),

            section_head("Currency Metabolites"),
            _qc_metric("Count",            str(_qc.n_currency)),
            _qc_metric("Edge fraction",    pct(_qc.currency_edge_fraction),
                       _qc.currency_edge_fraction < 0.5),

            section_head("Structural Quality"),
            _qc_metric("Dead-end metabolites",  str(_qc.n_dead_end_metabolites),
                       _qc.n_dead_end_metabolites == 0),
            _qc_metric("Orphan reactions",      str(_qc.n_orphan_reactions),
                       _qc.n_orphan_reactions == 0),

            section_head("Connectivity"),
            _qc_metric("Weakly connected components",  str(_qc.n_weakly_connected_components),
                       _qc.n_weakly_connected_components <= 5),
            _qc_metric("Largest component fraction",   pct(_qc.largest_component_fraction),
                       _qc.largest_component_fraction > 0.8),
            _qc_metric("Isolated nodes",               str(_qc.n_isolated_nodes),
                       _qc.n_isolated_nodes == 0),
        ])

        right = html.Div(style=col_style, children=[
            section_head("Reaction Annotations"),
            _qc_metric("With EC number",   pct(_qc.reactions_with_ec_fraction),
                       _qc.reactions_with_ec_fraction > 0.5),
            _qc_metric("With gene symbol", pct(_qc.reactions_with_gene_fraction),
                       _qc.reactions_with_gene_fraction > 0.3),
            _qc_metric("With pathway",     pct(_qc.reactions_with_pathway_fraction),
                       _qc.reactions_with_pathway_fraction > 0.5),

            section_head("Metabolite Annotations"),
            _qc_metric("With molecular formula",  pct(_qc.metabolites_with_formula_fraction),
                       _qc.metabolites_with_formula_fraction > 0.7),
            _qc_metric("With InChIKey",            pct(_qc.metabolites_with_inchikey_fraction),
                       _qc.metabolites_with_inchikey_fraction > 0.7),
            _qc_metric("With ChEBI ID",            pct(_qc.metabolites_with_chebi_fraction),
                       _qc.metabolites_with_chebi_fraction > 0.8),
            _qc_metric("Compartments",             ", ".join(_qc.compartments) or "none"),

            section_head("Metabolon Coverage"),
            _qc_metric("Compounds with Metabolon name",  str(_qc.metabolon_compounds_total)),
            _qc_metric("Metabolon → ChEBI coverage",    pct(_qc.metabolon_chebi_coverage),
                       _qc.metabolon_chebi_coverage > 0.7),

            section_head("Clinical Overlay"),
            _qc_metric("Disease–reaction edges",   str(_qc.disease_reaction_edges)),
            _qc_metric("Disease–metabolite edges", str(_qc.disease_metabolite_edges)),
            _qc_metric("Disease–gene edges",       str(_qc.disease_gene_edges)),

            section_head("Overall Verdict"),
            html.Div(style={
                "textAlign": "center", "padding": "12px",
                "background": "#e8f5e9" if _qc.is_fba_ready else "#fce4ec",
                "borderRadius": "8px", "marginTop": "8px",
            }, children=[
                html.Div("FBA-READY" if _qc.is_fba_ready else "NOT FBA-READY", style={
                    "fontSize": "20px", "fontWeight": "bold",
                    "color": "#2e7d32" if _qc.is_fba_ready else "#c62828",
                }),
                html.Div(
                    "(> 100 metabolites, > 50 reactions, > 30% EC, > 50% formula, < 10 components)",
                    style={"fontSize": "10px", "color": "#666", "marginTop": "4px"},
                ),
            ]),
        ])

        return html.Div(style={"padding": "16px", "overflowY": "auto", "height": "calc(100vh - 120px)"}, children=[
            html.Div(style={"display": "flex", "flexWrap": "wrap", "gap": "0"}, children=[left, right]),
        ])

    qc_tab = dcc.Tab(label="QC Report", value="qc", children=[_qc_tab_content()])

    # ================================================================== #
    # Tab 3 — Compound Curation                                           #
    # ================================================================== #

    if curator is None:
        _cur_body = html.Div(style={"padding": "50px", "textAlign": "center", "color": "#888"}, children=[
            html.Div("Curation unavailable", style={"fontSize": "20px", "fontWeight": "bold", "marginBottom": "10px"}),
            html.Div("Pass --metabolon-csv to the CLI (or met_loader + curator to create_app()) to enable."),
        ])
    else:
        _cur_body = html.Div(style={"display": "flex", "height": "calc(100vh - 120px)"}, children=[

            # Unmatched table
            html.Div(style={"flex": "1", "padding": "14px", "overflowY": "auto"}, children=[
                html.H4("Unmatched Metabolon Compounds", style={"marginTop": 0, "fontSize": "15px"}),
                html.Div(id="cur-summary", style={"color": "#555", "marginBottom": "8px", "fontSize": "13px"}),
                dash_table.DataTable(
                    id="unmatched-table",
                    columns=[
                        {"name": "Name",     "id": "BIOCHEMICAL"},
                        {"name": "Pathway",  "id": "PATHWAY"},
                        {"name": "InChIKey", "id": "INCHIKEY"},
                        {"name": "PubChem",  "id": "PUBCHEM"},
                        {"name": "Mass",     "id": "MASS"},
                    ],
                    data=_unmatched_records(),
                    row_selectable="single",
                    selected_rows=[],
                    page_size=25,
                    style_table={"overflowX": "auto"},
                    style_cell={
                        "fontSize": "12px", "padding": "4px 8px",
                        "maxWidth": "220px", "overflow": "hidden", "textOverflow": "ellipsis",
                    },
                    style_header={"fontWeight": "bold", "background": "#e9ecef"},
                    style_data_conditional=[
                        {"if": {"row_index": "odd"}, "backgroundColor": "#f8f9fa"},
                    ],
                    filter_action="native",
                    sort_action="native",
                    tooltip_delay=0,
                    tooltip_duration=None,
                ),
            ]),

            # Curation panel
            html.Div(style={
                "width": "400px", "minWidth": "400px", "padding": "14px",
                "borderLeft": "1px solid #ddd", "background": "#f8f9fa", "overflowY": "auto",
            }, children=[
                html.H4("Curate Selected Compound", style={"marginTop": 0, "fontSize": "14px"}),
                html.Div(id="cur-compound-info",
                         children=html.Div("Select a compound from the table.", style={"color": "#999"})),
                html.Div(id="cur-pubchem-link"),

                html.Hr(),
                _lbl("Search ChEBI by name"),
                html.Div(style={"display": "flex", "gap": "6px"}, children=[
                    dcc.Input(
                        id="cur-search-input", type="text", placeholder="compound name…",
                        debounce=False,
                        style={"flex": 1, "fontSize": "13px", "padding": "4px 8px",
                               "border": "1px solid #ccc", "borderRadius": "4px"},
                    ),
                    html.Button("Search", id="cur-search-btn", n_clicks=0,
                                style={**_BTN_BASE, "fontSize": "13px", "padding": "4px 12px",
                                       "background": "#1565C0", "color": "#fff"}),
                ]),
                html.Div(id="cur-candidates", style={"marginTop": "8px"}),

                html.Hr(),
                _lbl("Search VMH/HMDB by name"),
                html.Div(style={"display": "flex", "gap": "6px"}, children=[
                    dcc.Input(
                        id="cur-hmdb-search-input", type="text", placeholder="compound name…",
                        debounce=False,
                        style={"flex": 1, "fontSize": "13px", "padding": "4px 8px",
                               "border": "1px solid #ccc", "borderRadius": "4px"},
                    ),
                    html.Button("Search", id="cur-hmdb-search-btn", n_clicks=0,
                                style={**_BTN_BASE, "fontSize": "13px", "padding": "4px 12px",
                                       "background": "#6A1B9A", "color": "#fff"}),
                ]),
                html.Div(id="cur-hmdb-candidates", style={"marginTop": "8px"}),

                html.Hr(),
                _lbl("Manual ChEBI assignment"),
                html.Div(style={"display": "flex", "gap": "6px", "marginBottom": "6px"}, children=[
                    dcc.Input(
                        id="cur-manual-chebi", type="text", placeholder="CHEBI:XXXXX",
                        style={"width": "150px", "fontSize": "13px", "padding": "4px 8px",
                               "border": "1px solid #ccc", "borderRadius": "4px"},
                    ),
                    html.Button("Assign", id="cur-assign-btn", n_clicks=0,
                                style={**_BTN_BASE, "fontSize": "13px", "padding": "4px 10px",
                                       "background": "#2e7d32", "color": "#fff"}),
                    html.Button("Exclude", id="cur-exclude-btn", n_clicks=0,
                                style={**_BTN_BASE, "fontSize": "13px", "padding": "4px 10px",
                                       "background": "#c62828", "color": "#fff"}),
                ]),
                _lbl("Manual HMDB assignment"),
                html.Div(style={"display": "flex", "gap": "6px", "marginBottom": "6px"}, children=[
                    dcc.Input(
                        id="cur-manual-hmdb", type="text", placeholder="HMDB0000001",
                        style={"width": "150px", "fontSize": "13px", "padding": "4px 8px",
                               "border": "1px solid #ccc", "borderRadius": "4px"},
                    ),
                    html.Button("Assign HMDB", id="cur-assign-hmdb-btn", n_clicks=0,
                                style={**_BTN_BASE, "fontSize": "13px", "padding": "4px 10px",
                                       "background": "#6A1B9A", "color": "#fff"}),
                ]),
                dcc.Input(
                    id="cur-notes", type="text", placeholder="Notes (optional)",
                    style={"width": "100%", "fontSize": "12px", "padding": "4px 8px",
                           "border": "1px solid #ccc", "borderRadius": "4px", "marginBottom": "6px"},
                ),
                html.Button("Save all overrides", id="cur-save-btn", n_clicks=0,
                            style={**_BTN_BASE, "fontSize": "13px", "padding": "5px 14px",
                                   "background": "#1565C0", "color": "#fff"}),
                html.Div(id="cur-status", style={"marginTop": "8px", "fontSize": "13px", "color": "#2e7d32"}),

                html.Hr(),
                html.Div(id="cur-graph-connection", style={"fontSize": "12px", "color": "#555"}),
            ]),
        ])

    curation_tab = dcc.Tab(label="Compound Curation", value="curation", children=[_cur_body])

    # ================================================================== #
    # App layout                                                           #
    # ================================================================== #

    app.layout = html.Div([
        header,
        stats_bar,
        dcc.Tabs(
            id="tabs", value="graph",
            children=[graph_tab, qc_tab, curation_tab],
            style={"height": "38px"},
        ),
        # Persistent stores
        dcc.Store(id="selected-node-store"),        # node ID for neighbour exploration
        dcc.Store(id="selected-compound-store"),    # Metabolon compound name
        dcc.Store(id="candidates-store", data=[]),  # ChEBI search candidates
        dcc.Store(id="hmdb-candidates-store", data=[]),
    ])

    # ================================================================== #
    # Callbacks — Graph Explorer                                           #
    # ================================================================== #

    @app.callback(
        Output("node-search", "value"),
        Input("node-search-clear", "n_clicks"),
        prevent_initial_call=True,
    )
    def clear_search(_: int) -> str:
        return ""

    @app.callback(
        Output("cytoscape", "elements"),
        Output("cytoscape", "layout"),
        Output("node-count-badge", "children"),
        Output("graph-empty-msg", "style"),
        Input("pw-selector", "value"),
        Input("type-filter", "value"),
        Input("color-by", "value"),
        Input("layout-selector", "value"),
        Input("show-all", "value"),
        Input("node-search", "value"),
        Input("edge-filter", "value"),
        Input("btn-neighbors", "n_clicks"),
        State("selected-node-store", "data"),
        State("hop-count", "value"),
    )
    def update_graph(
        selected_pws, type_filter, color_by, layout_name,
        show_all_val, search_query, edge_filter_val, _nbr_clicks,
        selected_node_id, hop_count,
    ):
        ctx        = callback_context
        trigger    = ctx.triggered[0]["prop_id"] if ctx.triggered else ""
        type_set   = set(type_filter or [])
        edge_roles = set(edge_filter_val or [])
        layout     = {"name": layout_name or "cose-bilkent", "animate": False, "randomize": False}
        _hidden    = {"display": "none"}
        _visible   = {
            "position": "absolute", "top": "50%", "left": "50%",
            "transform": "translate(-50%,-50%)", "textAlign": "center", "zIndex": 10,
        }

        def _type_ok(nid: str) -> bool:
            return not type_set or g.nodes[nid].get("node_type") in type_set

        def _render(node_ids: list[str], label: str):
            cap = len(node_ids) > max_nodes
            if cap:
                node_ids = node_ids[:max_nodes]
            overlay = (
                _drug_overlay if color_by == "druggability" else
                _pert_overlay if color_by == "perturbability" else
                None
            )
            elems  = build_elements(g, node_ids, color_by, platform_index, edge_roles, overlay)
            n_e    = sum(1 for e in elems if "source" in e["data"])
            badge  = f"{len(node_ids)} nodes · {n_e} edges" + (" (capped)" if cap else "")
            return elems, layout, f"{label} — {badge}", _hidden

        # Mode 1: neighbour exploration (button click)
        if "btn-neighbors" in trigger and selected_node_id and selected_node_id in g:
            nbrs = _neighborhood(g, selected_node_id, hops=hop_count or 1)
            ids  = [n for n in nbrs if _type_ok(n)]
            return _render(ids, f"{hop_count}-hop neighbourhood of {selected_node_id[:30]}")

        # Mode 2: node name / ID search
        if search_query and search_query.strip():
            hits     = _search_nodes(g, search_query)
            ctx_set: set[str] = set(hits)
            for nid in hits:
                ctx_set.update(g.predecessors(nid))
                ctx_set.update(g.successors(nid))
            ids = [n for n in ctx_set if _type_ok(n)]
            return _render(ids, f"Search '{search_query}': {len(hits)} hit(s)")

        # Mode 3: show all
        if show_all_val and "all" in show_all_val:
            ids = [n for n in g.nodes() if _type_ok(n)]
            return _render(ids, "All nodes")

        # Mode 4: pathway filter
        if selected_pws:
            node_set: set[str] = set()
            for pw in selected_pws:
                node_set.update(pathway_to_nodes.get(pw, set()))
            ids = [n for n in node_set if _type_ok(n)]
            return _render(ids, f"{len(selected_pws)} pathway(s)")

        # Empty state
        return [], layout, "", _visible

    @app.callback(
        Output("legend-div", "children"),
        Input("color-by", "value"),
    )
    def update_legend(color_by: str):
        def swatch(color: str, label: str) -> html.Div:
            return html.Div(style={"display": "flex", "alignItems": "center", "marginBottom": "3px"}, children=[
                html.Div(style={
                    "width": "12px", "height": "12px", "borderRadius": "50%",
                    "background": color, "marginRight": "6px", "flexShrink": 0,
                }),
                html.Span(label, style={"fontSize": "11px"}),
            ])

        if color_by == "node_type":
            return [swatch(c, t.capitalize()) for t, c in NODE_COLORS.items() if t != "unknown"]
        if color_by == "is_currency":
            return [swatch(CURRENCY_COLOR, "Currency"), swatch(NODE_COLORS["metabolite"], "Non-currency")]
        if color_by == "has_chebi":
            return [swatch(MATCHED_COLOR, "Has ChEBI"), swatch(UNMATCHED_COLOR, "No ChEBI")]
        if color_by == "has_metabolon":
            return [swatch("#9C27B0", "Has Metabolon"), swatch(NODE_COLORS["metabolite"], "Not Metabolon")]
        if color_by == "platform":
            return [swatch(_PLATFORM_PALETTE[i % len(_PLATFORM_PALETTE)], p)
                    for i, p in enumerate(platforms[:8])]
        if color_by in ("druggability", "perturbability"):
            labels = ["Low (0)", "·", "·", "·", "·", "·", "·", "High (1)"]
            return [swatch(_SCORE_GRADIENT[i], labels[i]) for i in range(len(_SCORE_GRADIENT))]
        return []

    @app.callback(
        Output("props-panel", "children"),
        Output("selected-node-store", "data"),
        Output("explore-status", "children"),
        Input("cytoscape", "tapNodeData"),
        Input("cytoscape", "tapEdgeData"),
    )
    def update_props(node_data, edge_data):
        ctx = callback_context
        if not ctx.triggered:
            return html.Div("Click a node or edge.", style={"color": "#999", "fontSize": "12px"}), no_update, ""
        trigger = ctx.triggered[0]["prop_id"]

        if "tapNodeData" in trigger and node_data:
            nid   = node_data.get("id", "")
            ntype = node_data.get("node_type", "")
            rows  = [
                _prop_row("ID",    nid),
                _prop_row("Type",  ntype),
                _prop_row("Label", node_data.get("label")),
            ]
            if ntype == "metabolite":
                rows += [
                    _prop_row("ChEBI",    node_data.get("chebi_id") or "—"),
                    _prop_row("HMDB",     node_data.get("hmdb_id") or "—"),
                    _prop_row("Metabolon", node_data.get("metabolon_name") or "—"),
                    _prop_row("InChIKey", node_data.get("inchikey") or "—"),
                    _prop_row("Formula",  node_data.get("formula") or "—"),
                    _prop_row("Mass",     node_data.get("mass") or "—"),
                    _prop_row("Platform", node_data.get("platform") or "—"),
                    _prop_row("Currency", node_data.get("is_currency")),
                ]
            elif ntype == "reaction":
                rows += [
                    _prop_row("EC",        node_data.get("ec_numbers") or "—"),
                    _prop_row("Genes",     node_data.get("gene_symbols") or "—"),
                    _prop_row("Pathways",  node_data.get("pathways") or "—"),
                    _prop_row("Reversible", node_data.get("reversible") or "—"),
                ]
            elif ntype == "gene":
                rows += [
                    _prop_row("Symbol",          node_data.get("symbol") or "—"),
                    _prop_row("Ensembl",         node_data.get("ensembl_id") or "—"),
                    _prop_row("HGNC",            node_data.get("hgnc_id") or "—"),
                    _prop_row("Entrez",          node_data.get("entrez_id") or "—"),
                    _prop_row("Species",         node_data.get("species") or "—"),
                    _prop_row("Expression",      node_data.get("tissue_expression") or "—"),
                ]
            elif ntype == "disease":
                rows += [
                    _prop_row("MONDO",           node_data.get("mondo_id") or "—"),
                    _prop_row("OMIM",            node_data.get("xref_omim") or "—"),
                    _prop_row("Orphanet",        node_data.get("xref_orphanet") or "—"),
                    _prop_row("ICD-10",          node_data.get("xref_icd10") or "—"),
                    _prop_row("DOID",            node_data.get("xref_doid") or "—"),
                    _prop_row("Rare",            node_data.get("is_rare") or "—"),
                    _prop_row("IEM",             node_data.get("is_inborn_error_of_metabolism") or "—"),
                    _prop_row("Definition",      node_data.get("definition") or "—"),
                ]
            elif ntype == "pathway":
                rows += [
                    _prop_row("Reactome",        node_data.get("reactome_id") or "—"),
                    _prop_row("KEGG",            node_data.get("kegg_id") or "—"),
                    _prop_row("Species",         node_data.get("species") or "—"),
                ]
            elif ntype == "phenotype":
                rows += [
                    _prop_row("HPO",             node_data.get("hpo_id") or "—"),
                    _prop_row("Metabolic",       node_data.get("is_metabolic") or "—"),
                    _prop_row("Definition",      node_data.get("definition") or "—"),
                ]
            elif ntype == "drug":
                rows += [
                    _prop_row("ChEMBL",          node_data.get("chembl_id") or "—"),
                    _prop_row("DrugBank",        node_data.get("drugbank_id") or "—"),
                    _prop_row("Max phase",       node_data.get("max_phase") or "—"),
                    _prop_row("Mechanism",       node_data.get("mechanism") or "—"),
                ]
            if nid in g:
                rows.append(_prop_row("Degree in / out",
                                      f"{g.in_degree(nid)} / {g.out_degree(nid)}"))
            table = html.Table(rows, style={"borderCollapse": "collapse", "width": "100%"})
            return table, nid, f"Selected: {nid[:40]}"

        if "tapEdgeData" in trigger and edge_data:
            rows = [
                _prop_row("From",          edge_data.get("source")),
                _prop_row("To",            edge_data.get("target")),
                _prop_row("Role",          edge_data.get("role") or "—"),
                _prop_row("Stoichiometry", edge_data.get("stoichiometry") or "—"),
            ]
            # Extra STRING-specific fields stored in edge data
            if edge_data.get("combined_score"):
                rows += [
                    _prop_row("Combined score", edge_data.get("combined_score")),
                    _prop_row("Experimental",   edge_data.get("experimental") or "—"),
                    _prop_row("Co-expression",  edge_data.get("coexpression") or "—"),
                    _prop_row("Database",       edge_data.get("database") or "—"),
                    _prop_row("Text mining",    edge_data.get("textmining") or "—"),
                ]
            return html.Table(rows, style={"borderCollapse": "collapse", "width": "100%"}), no_update, ""

        return html.Div("Click a node or edge.", style={"color": "#999", "fontSize": "12px"}), no_update, ""

    # ================================================================== #
    # Callbacks — Curation tab (only when curator is available)           #
    # ================================================================== #

    if curator is not None:

        @app.callback(
            Output("cur-summary", "children"),
            Output("unmatched-table", "data"),
            Input("cur-status", "children"),
            Input("tabs", "value"),
        )
        def refresh_curation_table(_, __):
            n_un = len(curator.unmatched)
            n_as = len(curator.assigned)
            n_ex = len(curator.excluded)
            summary = f"{n_un} unmatched · {n_as} manually assigned · {n_ex} excluded"
            return summary, _unmatched_records()

        @app.callback(
            Output("cur-compound-info", "children"),
            Output("cur-search-input", "value"),
            Output("cur-hmdb-search-input", "value"),
            Output("cur-pubchem-link", "children"),
            Output("selected-compound-store", "data"),
            Input("unmatched-table", "selected_rows"),
            State("unmatched-table", "data"),
        )
        def select_compound(selected_rows, table_data):
            if not selected_rows or not table_data:
                return html.Div("Select a compound.", style={"color": "#999"}), "", "", None, None
            row  = table_data[selected_rows[0]]
            name = row.get("BIOCHEMICAL", "")
            info = html.Table([
                _prop_row("Name",    name),
                _prop_row("Pathway", row.get("PATHWAY") or "—"),
                _prop_row("InChIKey", row.get("INCHIKEY") or "—"),
                _prop_row("Mass",    row.get("MASS") or "—"),
            ], style={"borderCollapse": "collapse", "width": "100%", "marginBottom": "6px"})
            cid  = str(row.get("PUBCHEM") or "")
            link = (
                html.A(f"PubChem CID {cid} ↗",
                       href=f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}",
                       target="_blank", style={"fontSize": "12px"})
                if cid and cid not in ("", "nan") else None
            )
            return info, name, name, link, name

        @app.callback(
            Output("cur-candidates", "children"),
            Output("candidates-store", "data"),
            Input("cur-search-btn", "n_clicks"),
            State("cur-search-input", "value"),
            prevent_initial_call=True,
        )
        def search_chebi(_, query):
            if not query:
                return html.Div("Enter a name first.", style={"color": "#999"}), []
            results = curator.search_chebi(query.strip(), n=8)
            if not results:
                return html.Div("No ChEBI results found.", style={"color": "#999"}), []
            rows = []
            for r in results:
                conn = ""
                if r.get("in_graph") is True:
                    conn = f" · {r['reaction_count']} reactions in graph"
                elif r.get("in_graph") is False:
                    conn = " · not in graph"
                rows.append(html.Div(
                    f"{r['chebi_id']} — {r['name']}{conn}",
                    id={"type": "candidate-row", "chebi": r["chebi_id"]},
                    n_clicks=0,
                    style={
                        "padding": "5px 8px", "marginBottom": "3px",
                        "background": "#e8f5e9" if r.get("in_graph") else "#fafafa",
                        "border": "1px solid #c8e6c9" if r.get("in_graph") else "1px solid #e0e0e0",
                        "borderRadius": "4px", "fontSize": "12px", "cursor": "pointer",
                    },
                ))
            return html.Div(rows), results

        @app.callback(
            Output("cur-hmdb-candidates", "children"),
            Output("hmdb-candidates-store", "data"),
            Input("cur-hmdb-search-btn", "n_clicks"),
            State("cur-hmdb-search-input", "value"),
            prevent_initial_call=True,
        )
        def search_hmdb(_, query):
            if not query:
                return html.Div("Enter a name first.", style={"color": "#999"}), []
            results = curator.search_hmdb(query.strip(), n=8)
            if not results:
                return html.Div("No VMH/HMDB results found.", style={"color": "#999"}), []
            rows = []
            for r in results:
                conn = ""
                if r.get("in_graph") is True:
                    conn = f" · {r['reaction_count']} reactions in graph"
                elif r.get("in_graph") is False:
                    conn = " · not in graph"
                chebi = f" · {r['chebi_id']}" if r.get("chebi_id") else ""
                rows.append(html.Div(
                    f"{r['hmdb_id']} — {r['name']}{chebi}{conn}",
                    id={"type": "hmdb-candidate-row", "hmdb": r["hmdb_id"]},
                    n_clicks=0,
                    style={
                        "padding": "5px 8px", "marginBottom": "3px",
                        "background": "#f3e5f5" if r.get("in_graph") else "#fafafa",
                        "border": "1px solid #e1bee7" if r.get("in_graph") else "1px solid #e0e0e0",
                        "borderRadius": "4px", "fontSize": "12px", "cursor": "pointer",
                    },
                ))
            return html.Div(rows), results

        @app.callback(
            Output("cur-manual-chebi", "value"),
            Input({"type": "candidate-row", "chebi": ALL}, "n_clicks"),
            State("candidates-store", "data"),
            prevent_initial_call=True,
        )
        def fill_chebi_from_candidate(n_clicks_list, _candidates):
            ctx = callback_context
            if not ctx.triggered or all(n == 0 for n in (n_clicks_list or [])):
                return no_update
            try:
                chebi = json.loads(ctx.triggered[0]["prop_id"].split(".")[0])["chebi"]
                return chebi
            except Exception:
                return no_update

        @app.callback(
            Output("cur-manual-hmdb", "value"),
            Input({"type": "hmdb-candidate-row", "hmdb": ALL}, "n_clicks"),
            State("hmdb-candidates-store", "data"),
            prevent_initial_call=True,
        )
        def fill_hmdb_from_candidate(n_clicks_list, _candidates):
            ctx = callback_context
            if not ctx.triggered or all(n == 0 for n in (n_clicks_list or [])):
                return no_update
            try:
                hmdb = json.loads(ctx.triggered[0]["prop_id"].split(".")[0])["hmdb"]
                return hmdb
            except Exception:
                return no_update

        @app.callback(
            Output("cur-status", "children"),
            Output("cur-graph-connection", "children"),
            Input("cur-assign-btn", "n_clicks"),
            Input("cur-assign-hmdb-btn", "n_clicks"),
            Input("cur-exclude-btn", "n_clicks"),
            Input("cur-save-btn", "n_clicks"),
            State("selected-compound-store", "data"),
            State("cur-manual-chebi", "value"),
            State("cur-manual-hmdb", "value"),
            State("cur-notes", "value"),
            prevent_initial_call=True,
        )
        def handle_curation_action(_, __, ___, ____, compound_name, chebi_input, hmdb_input, notes):
            ctx     = callback_context
            trigger = ctx.triggered[0]["prop_id"]

            if "assign-btn" in trigger:
                if not compound_name:
                    return "No compound selected.", ""
                if not chebi_input:
                    return "Enter a ChEBI ID first.", ""
                chebi = chebi_input.strip()
                curator.assign(compound_name, chebi, notes=notes or "")
                info  = curator.verify_connection(chebi)
                if info["in_graph"] is True:
                    conn = f"✓ {chebi} is in graph · {info['reaction_count']} reactions"
                    if info.get("pathways"):
                        conn += f" · Pathways: {', '.join(info['pathways'][:3])}"
                elif info["in_graph"] is False:
                    conn = f"⚠ {chebi} is not in the reaction graph."
                else:
                    conn = ""
                return f"✓ Assigned {chebi} → '{compound_name}'", conn

            elif "assign-hmdb-btn" in trigger:
                if not compound_name:
                    return "No compound selected.", ""
                if not hmdb_input:
                    return "Enter an HMDB ID first.", ""
                hmdb = hmdb_input.strip().upper()
                curator.assign(compound_name, hmdb_id=hmdb, notes=notes or "")
                info = curator.verify_hmdb_connection(hmdb)
                if info["in_graph"] is True:
                    conn = f"✓ {hmdb} is in graph · {info['reaction_count']} reactions"
                    if info.get("pathways"):
                        conn += f" · Pathways: {', '.join(info['pathways'][:3])}"
                elif info["in_graph"] is False:
                    conn = f"⚠ {hmdb} is not on any graph metabolite node."
                else:
                    conn = ""
                return f"✓ Assigned {hmdb} → '{compound_name}'", conn

            elif "exclude-btn" in trigger:
                if not compound_name:
                    return "No compound selected.", ""
                curator.exclude(compound_name, notes=notes or "")
                return f"Excluded '{compound_name}'", ""

            elif "save-btn" in trigger:
                curator.save()
                curator.apply()
                return f"✓ Saved {len(curator.overrides)} overrides.", ""

            return no_update, no_update

    return app
