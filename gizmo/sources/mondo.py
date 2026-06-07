"""
MONDO (Monarch Disease Ontology) client.

License: CC BY 4.0 — https://mondo.monarchinitiative.org/
OBO file: http://purl.obolibrary.org/obo/mondo.obo

MONDO cross-maps OMIM, Orphanet, DOID, ICD-10, MeSH IDs as xrefs.
We store those IDs as strings without pulling any data from restricted sources.

Inborn errors of metabolism (IEM) are identified by MONDO subclass hierarchy:
  MONDO:0004736  inherited metabolic disorder
  MONDO:0019052  inborn errors of metabolism (Orphanet grouping)
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional
from urllib.request import urlretrieve

import obonet

from gizmo.schema import DiseaseNode

log = logging.getLogger(__name__)

_MONDO_OBO_URL = "http://purl.obolibrary.org/obo/mondo.obo"
_MONDO_SLIM_URL = "http://purl.obolibrary.org/obo/mondo/subsets/mondo-rare.obo"

# MONDO IDs for inherited metabolic disorders
_IEM_ANCESTORS = {
    "MONDO:0004736",  # inherited metabolic disorder
    "MONDO:0019052",  # inborn error of metabolism (Orphanet grouping)
    "MONDO:0005066",  # metabolic disease
}

# Xref prefix → field mapping
_XREF_MAP = {
    "OMIM":    "xref_omim",
    "Orphanet": "xref_orphanet",
    "DOID":    "xref_doid",
    "ICD10CM": "xref_icd10",
    "ICD10":   "xref_icd10",
    "MESH":    "xref_mesh",
    "MeSH":    "xref_mesh",
}


class MondoClient:
    """
    Parse MONDO OBO to produce DiseaseNode objects.

    Usage::

        client = MondoClient(cache_dir="data/raw/mondo")
        client.download()
        diseases = client.load_all()            # all diseases
        iem = client.load_iem_subset()          # IEM only
    """

    def __init__(self, cache_dir: str | Path = "data/raw/mondo") -> None:
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self._graph: Optional[object] = None   # obonet graph

    @property
    def obo_path(self) -> Path:
        return self.cache_dir / "mondo.obo"

    def download(self, force: bool = False) -> None:
        if self.obo_path.exists() and not force:
            log.info("MONDO OBO cache hit: %s", self.obo_path)
            return
        log.info("Downloading MONDO OBO → %s", self.obo_path)
        urlretrieve(_MONDO_OBO_URL, self.obo_path)

    def _load_graph(self) -> object:
        if self._graph is None:
            if not self.obo_path.exists():
                raise FileNotFoundError(f"MONDO OBO not found at {self.obo_path}. Run .download() first.")
            log.info("Parsing MONDO OBO …")
            self._graph = obonet.read_obo(str(self.obo_path))
        return self._graph

    def load_all(self, species_filter: Optional[str] = None) -> list[DiseaseNode]:
        """
        Parse all MONDO disease terms into DiseaseNode objects.
        Filters out obsolete terms.
        """
        g = self._load_graph()
        nodes: list[DiseaseNode] = []
        for node_id, data in g.nodes(data=True):
            if not node_id.startswith("MONDO:"):
                continue
            if data.get("is_obsolete"):
                continue
            dn = _term_to_disease_node(node_id, data)
            if dn:
                nodes.append(dn)
        log.info("Loaded %d MONDO disease nodes", len(nodes))
        return nodes

    def load_iem_subset(self) -> list[DiseaseNode]:
        """
        Return DiseaseNodes for inborn errors of metabolism only.
        Uses ancestor traversal in the MONDO hierarchy.
        """
        g = self._load_graph()
        iem_ids = _descendants_of(g, _IEM_ANCESTORS)

        nodes: list[DiseaseNode] = []
        for node_id in iem_ids:
            data = g.nodes.get(node_id, {})
            if data.get("is_obsolete"):
                continue
            dn = _term_to_disease_node(node_id, data, is_inborn_error=True)
            if dn:
                nodes.append(dn)
        log.info("Loaded %d IEM disease nodes from MONDO", len(nodes))
        return nodes

    def rare_diseases(self) -> list[DiseaseNode]:
        """Return DiseaseNodes for MONDO terms that have an Orphanet xref (proxy for rare)."""
        nodes = self.load_all()
        return [n for n in nodes if n.xref_orphanet]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _term_to_disease_node(
    term_id: str,
    data: dict,
    is_inborn_error: bool = False,
) -> Optional[DiseaseNode]:
    name = data.get("name", "")
    if not name:
        return None

    synonyms: list[str] = []
    for syn in data.get("synonym", []):
        # OBO synonym format: "label" SCOPE [source]
        val = syn.strip('"').split('"')[0] if '"' in syn else syn
        synonyms.append(val.strip())

    xref_omim: list[str] = []
    xref_orphanet: list[str] = []
    xref_doid: list[str] = []
    xref_icd10: list[str] = []
    xref_mesh: list[str] = []

    for xref in data.get("xref", []):
        prefix = xref.split(":")[0] if ":" in xref else ""
        field = _XREF_MAP.get(prefix)
        if field == "xref_omim":
            xref_omim.append(xref)
        elif field == "xref_orphanet":
            xref_orphanet.append(xref)
        elif field == "xref_doid":
            xref_doid.append(xref)
        elif field == "xref_icd10":
            xref_icd10.append(xref)
        elif field == "xref_mesh":
            xref_mesh.append(xref)

    return DiseaseNode(
        node_id=term_id,
        mondo_id=term_id,
        name=name,
        synonyms=synonyms,
        definition=data.get("def", "").strip('"') if data.get("def") else None,
        xref_omim=xref_omim,
        xref_orphanet=xref_orphanet,
        xref_doid=xref_doid,
        xref_icd10=xref_icd10,
        xref_mesh=xref_mesh,
        is_rare=bool(xref_orphanet),
        is_inborn_error_of_metabolism=is_inborn_error,
    )


def _descendants_of(g: object, ancestor_ids: set[str]) -> set[str]:
    """Return all descendant node IDs of a set of ancestor IDs in an obonet graph."""
    import networkx as nx
    # obonet stores is_a as directed edges parent → child? No: child→parent.
    # Descendants = nodes reachable from ancestors in *reverse* edge direction.
    descendants: set[str] = set()
    for anc in ancestor_ids:
        if anc in g:
            # nx.ancestors follows edges in reverse direction; obonet edges go child→parent
            descendants |= nx.ancestors(g, anc)  # type: ignore[attr-defined]
            descendants.add(anc)
    return descendants
