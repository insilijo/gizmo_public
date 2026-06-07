# Source modules are imported lazily in build pipeline steps to avoid hard
# dependencies on optional packages (obonet, etc.) at import time.
# Explicit imports below are for packages with stable, required dependencies only.

from gizmo.sources.reactome import ReactomeClient, ReactomeLoader
from gizmo.sources.chebi import ChebiClient
from gizmo.sources.metanetx import MetaNetXClient
from gizmo.sources.open_targets import OpenTargetsClient
from gizmo.sources.metabolon import MetabolonLoader
from gizmo.sources.vmh import VMHMapper, build_hmdb_to_vmh, build_pubchem_to_vmh, enrich_graph_vmh, fetch_vmh_metabolites

__all__ = [
    "ReactomeClient",
    "ReactomeLoader",
    "ChebiClient",
    "MetaNetXClient",
    "OpenTargetsClient",
    "MetabolonLoader",
    "VMHMapper",
    "build_hmdb_to_vmh",
    "build_pubchem_to_vmh",
    "enrich_graph_vmh",
    "fetch_vmh_metabolites",
    # Lazy-import sources (require optional deps):
    # MondoClient, OrphanetClient  — require obonet
    # HPOClient                    — require obonet
    # GTExClient                   — stdlib only
    # ClinVarClient                — stdlib only
    # ChEMBLDrugClient             — stdlib only
    # OrthologMapper               — stdlib only
]
