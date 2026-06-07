"""Lipid module for GIZMO (paper 2 — lipid layer).

Foundational types and parsers for the lipid sub-graph and physiologic
state model. See ``docs/lipid_layer_design.md``.
"""

from gizmo.lipid.identity import (
    LipidCategory,
    LipidIdentity,
    ResolutionLevel,
)
from gizmo.lipid.pathways import (
    LipidTransformation,
    all_genes_referenced,
    all_pathways,
    all_subclasses_referenced,
    iter_transformations,
)
from gizmo.lipid.refmet import RefmetParser, parse_refmet
from gizmo.lipid.taxonomy import (
    SubclassEntry,
    get_subclass,
    iter_subclasses,
)

__all__ = [
    "LipidCategory",
    "LipidIdentity",
    "ResolutionLevel",
    "RefmetParser",
    "parse_refmet",
    "SubclassEntry",
    "get_subclass",
    "iter_subclasses",
    "LipidTransformation",
    "iter_transformations",
    "all_pathways",
    "all_subclasses_referenced",
    "all_genes_referenced",
]
