"""
GIZMO WSGI entry point for production deployment (gunicorn).

Usage::

    gunicorn gizmo.app.wsgi:server -b 0.0.0.0:8050 -w 2 --timeout 120

Environment variables
---------------------
GIZMO_GRAPH             Path to the graph JSON file (required).
GIZMO_REACTOME_CACHE    Reactome cache directory (default: data/raw/reactome).
GIZMO_METABOLON_CSV     Metabolon CSV path — enables the Curation tab (optional).
GIZMO_METANETX_PROP     chem_prop.tsv path (default: data/resources/gizmo/metanetx/chem_prop.tsv).
GIZMO_METANETX_XREF     chem_xref.tsv path (default: data/resources/gizmo/metanetx/chem_xref.tsv).
GIZMO_OVERRIDES         Curation overrides JSON (default: data/resources/gizmo/curation/metabolon_overrides.json).
GIZMO_MAX_NODES         Max nodes rendered per view (default: 500).
"""

from __future__ import annotations

import logging
import os
from pathlib import Path

from gizmo.resources import (
    metabolon_csv_default,
    metanetx_prop_default,
    metanetx_xref_default,
    overrides_default,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(name)s  %(message)s",
)
log = logging.getLogger("gizmo.wsgi")

_graph_path     = os.environ.get("GIZMO_GRAPH", "data/processed/gizmo_full.json")
_reactome_cache = os.environ.get("GIZMO_REACTOME_CACHE", "data/raw/reactome")
_metabolon_csv  = os.environ.get("GIZMO_METABOLON_CSV", metabolon_csv_default())
_metanetx_prop  = os.environ.get("GIZMO_METANETX_PROP", metanetx_prop_default())
_metanetx_xref  = os.environ.get("GIZMO_METANETX_XREF", metanetx_xref_default())
_overrides_path = os.environ.get("GIZMO_OVERRIDES", overrides_default())
_max_nodes      = int(os.environ.get("GIZMO_MAX_NODES", "500"))

from gizmo.export.json_export import read_json  # noqa: E402
from gizmo.app.dash_app import create_app       # noqa: E402

log.info("Loading graph from %s …", _graph_path)
mg = read_json(_graph_path)
log.info("Graph loaded: %d nodes, %d edges", mg.graph.number_of_nodes(), mg.graph.number_of_edges())

met_loader = None
curator    = None

if _metabolon_csv and Path(_metabolon_csv).exists():
    from gizmo.sources.metabolon import MetabolonLoader              # noqa: E402
    from gizmo.curation.metabolon_curator import MetabolonCurator   # noqa: E402

    log.info("Loading Metabolon CSV from %s …", _metabolon_csv)
    met_loader = MetabolonLoader(_metabolon_csv)

    prop = Path(_metanetx_prop)
    xref = Path(_metanetx_xref)
    if prop.exists() and xref.exists():
        log.info("Building MetaNetX index (streaming ~1.4 GB) …")
        met_loader.load_metanetx_index(str(prop), str(xref))
    else:
        log.warning("MetaNetX files not found — skipping local InChIKey index.")

    curator = MetabolonCurator(
        met_loader, graph=mg, overrides_path=_overrides_path,
    )
    curator.apply()
    log.info("Curator ready — %d unmatched compounds.", len(curator.unmatched))

log.info("Creating Dash app …")
app = create_app(
    mg,
    met_loader=met_loader,
    curator=curator,
    reactome_cache_dir=_reactome_cache,
    max_nodes=_max_nodes,
)

# Gunicorn / uWSGI entry point
server = app.server

log.info("GIZMO WSGI server ready.")
