"""
GIZMO — Graph-Integrated Zone of Metabolite Operations

A metabolite-centered, open-licensed reaction graph for systems biology and cheminformatics.

Data lineage (all open-licensed):
  - Reactome     CC BY 4.0   https://reactome.org
  - ChEBI        CC BY 4.0   https://www.ebi.ac.uk/chebi/
  - MetaNetX     CC BY 4.0   https://www.metanetx.org
"""

from gizmo.graph.network import GizmoGraph

__all__ = ["GizmoGraph"]
__version__ = "0.1.0"
