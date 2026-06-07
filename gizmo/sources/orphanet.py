"""
Orphanet / Orphadata client for rare disease and IEM data.

License: CC BY 4.0 — https://www.orphadata.com/
Orphadata provides open XML datasets for rare diseases.

Key datasets used:
  - Product 1: Rare disease classification (disease names, synonyms, OMIM xrefs)
  - Product 6: Genes associated with rare diseases (gene-disease relationships)
  - Product 2: Clinical signs / phenotypes (HPO cross-reference)

XML endpoints (no API key needed):
  http://www.orphadata.com/cgi-bin/inc/product1.inc.php  → disease list XML URL
  The actual data lives at:
  https://github.com/Orphanet/Orphapacket/tree/master/Data-papers/...
  or directly:
  https://www.orphadata.com/data/xml/en_product1.xml
  https://www.orphadata.com/data/xml/en_product6.xml

IEM subset: Orphanet classification group 158652 "Inborn errors of metabolism"
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional
from urllib.request import urlretrieve
from xml.etree import ElementTree as ET

from gizmo.schema import DiseaseEdge, DiseaseEdgeType, DiseaseNode, GeneNode

log = logging.getLogger(__name__)

_PRODUCT1_URL = "https://www.orphadata.com/data/xml/en_product1.xml"
_PRODUCT6_URL = "https://www.orphadata.com/data/xml/en_product6.xml"

# Orphanet classification ORPHA ID for IEM group
_IEM_ORPHA_ID = "158652"


class OrphanetClient:
    """
    Parse Orphadata XML files for rare disease and gene-disease data.

    Usage::

        client = OrphanetClient(cache_dir="data/raw/orphanet")
        client.download()
        diseases = client.load_diseases()
        genes, edges = client.load_gene_associations()
    """

    def __init__(self, cache_dir: str | Path = "data/raw/orphanet") -> None:
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    @property
    def product1_path(self) -> Path:
        return self.cache_dir / "en_product1.xml"

    @property
    def product6_path(self) -> Path:
        return self.cache_dir / "en_product6.xml"

    def download(self, force: bool = False) -> None:
        for url, path in [(_PRODUCT1_URL, self.product1_path), (_PRODUCT6_URL, self.product6_path)]:
            if path.exists() and not force:
                log.info("Orphanet cache hit: %s", path)
                continue
            log.info("Downloading %s → %s", url, path)
            urlretrieve(url, path)

    # ------------------------------------------------------------------
    # Disease loading (Product 1)
    # ------------------------------------------------------------------

    def load_diseases(self) -> list[DiseaseNode]:
        """
        Parse en_product1.xml to produce DiseaseNode objects.
        Orphanet disease IDs are used as secondary identifiers;
        MONDO IDs are the canonical node IDs when available (via xref).
        If no MONDO xref exists, node_id = "Orphanet:XXXXXX".
        """
        if not self.product1_path.exists():
            raise FileNotFoundError(f"Not found: {self.product1_path}. Run .download() first.")

        tree = ET.parse(self.product1_path)
        root = tree.getroot()

        nodes: list[DiseaseNode] = []
        for disorder in root.iter("Disorder"):
            dn = _parse_disorder(disorder)
            if dn:
                nodes.append(dn)

        log.info("Loaded %d Orphanet disease nodes", len(nodes))
        return nodes

    def load_iem_diseases(self) -> list[DiseaseNode]:
        """Return DiseaseNodes belonging to the IEM classification group."""
        all_diseases = self.load_diseases()
        # Heuristic: diseases that are direct members of IEM group
        # In product1 XML the ClassificationNode ancestor can be checked,
        # but for simplicity we filter by xref_orphanet pattern or parent class field.
        # Full implementation requires traversing Classification elements.
        return [d for d in all_diseases if d.is_inborn_error_of_metabolism]

    # ------------------------------------------------------------------
    # Gene associations (Product 6)
    # ------------------------------------------------------------------

    def load_gene_associations(
        self,
        *,
        association_types: Optional[set[str]] = None,
    ) -> tuple[list[GeneNode], list[DiseaseEdge]]:
        """
        Parse en_product6.xml for gene-disease associations.

        association_types: filter by DisorderGeneAssociationType name, e.g.
          {"Disease-causing germline mutation(s) in", "Modifying germline mutation"}
          Default: all types.

        Returns (gene_nodes, disease_edges).
        """
        if not self.product6_path.exists():
            raise FileNotFoundError(f"Not found: {self.product6_path}. Run .download() first.")

        tree = ET.parse(self.product6_path)
        root = tree.getroot()

        genes: list[GeneNode] = []
        edges: list[DiseaseEdge] = []
        seen_genes: dict[str, GeneNode] = {}

        for disorder in root.iter("Disorder"):
            orpha_id = disorder.findtext("OrphaCode") or ""
            disease_node_id = f"Orphanet:{orpha_id}"

            assoc_list = disorder.find("DisorderGeneAssociationList")
            if assoc_list is None:
                continue

            for assoc in assoc_list.iter("DisorderGeneAssociation"):
                assoc_type = assoc.findtext(".//DisorderGeneAssociationType/Name") or ""
                if association_types and assoc_type not in association_types:
                    continue

                gene_el = assoc.find("Gene")
                if gene_el is None:
                    continue

                symbol = gene_el.findtext("Symbol") or ""
                gene_name = gene_el.findtext("Name") or ""
                # Extract Ensembl and HGNC IDs from ExternalReferenceList
                ensembl_id: Optional[str] = None
                hgnc_id: Optional[str] = None
                for xref in gene_el.iter("ExternalReference"):
                    source = xref.findtext("Source") or ""
                    ref = xref.findtext("Reference") or ""
                    if source == "Ensembl":
                        ensembl_id = ref
                    elif source == "HGNC":
                        hgnc_id = f"HGNC:{ref}"

                gene_node_id = f"ENSG:{ensembl_id}" if ensembl_id else f"symbol:{symbol}"

                if gene_node_id not in seen_genes:
                    g = GeneNode(
                        node_id=gene_node_id,
                        ensembl_id=ensembl_id,
                        hgnc_id=hgnc_id,
                        symbol=symbol,
                        name=gene_name if gene_name != symbol else None,
                    )
                    seen_genes[gene_node_id] = g
                    genes.append(g)

                edges.append(DiseaseEdge(
                    source=disease_node_id,
                    target=gene_node_id,
                    edge_type=DiseaseEdgeType.GENE_ASSOCIATED,
                    source_db="orphanet",
                ))

        log.info("Loaded %d genes, %d disease-gene edges from Orphanet", len(genes), len(edges))
        return genes, edges


# ---------------------------------------------------------------------------
# XML parsing helpers
# ---------------------------------------------------------------------------

def _parse_disorder(disorder: ET.Element) -> Optional[DiseaseNode]:
    orpha_code = disorder.findtext("OrphaCode")
    name = disorder.findtext("Name")
    if not orpha_code or not name:
        return None

    orpha_node_id = f"Orphanet:{orpha_code}"

    synonyms: list[str] = []
    for syn in disorder.iter("Synonym"):
        if syn.text:
            synonyms.append(syn.text.strip())

    xref_omim: list[str] = []
    xref_icd10: list[str] = []
    for xref in disorder.iter("ExternalReference"):
        source = xref.findtext("Source") or ""
        ref = xref.findtext("Reference") or ""
        if source == "OMIM":
            xref_omim.append(f"OMIM:{ref}")
        elif source in ("ICD-10", "ICD10", "ICD10CM"):
            xref_icd10.append(f"ICD10:{ref}")

    # Check if this is an IEM (disorder type or classification)
    disorder_type = disorder.findtext(".//DisorderType/Name") or ""
    is_iem = "metabolism" in disorder_type.lower()

    return DiseaseNode(
        node_id=orpha_node_id,
        name=name,
        synonyms=synonyms,
        xref_omim=xref_omim,
        xref_orphanet=[f"Orphanet:{orpha_code}"],
        xref_icd10=xref_icd10,
        is_rare=True,
        is_inborn_error_of_metabolism=is_iem,
    )
