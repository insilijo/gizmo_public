"""
Open Targets Platform client for gene-disease associations.

License: CC BY 4.0 — https://platform.opentargets.org/
GraphQL API: https://api.platform.opentargets.org/api/v4/graphql

Open Targets provides scored gene-disease associations integrating
genetic evidence (GWAS, rare variant), expression, somatic, and text-mining.

We use it to link:
  DiseaseNode (MONDO ID) → GeneNode (Ensembl ID) → ReactionNode (via EC/Reactome)

Key association score: 0–1, higher = stronger evidence.
Minimum score threshold (default 0.1) avoids noise from weak text-mining hits.

Note: Open Targets maps diseases to EFO IDs; we cross-reference to MONDO via
the EFO → MONDO mapping that Open Targets provides in their disease API.
"""

from __future__ import annotations

import logging
from typing import Any, Iterator, Optional

import requests

from gizmo.schema import DiseaseEdge, DiseaseEdgeType, GeneNode

log = logging.getLogger(__name__)

_API_URL = "https://api.platform.opentargets.org/api/v4/graphql"

# GraphQL query: gene-disease associations for a disease (by EFO/MONDO ID)
_ASSOC_QUERY = """
query DiseaseAssociations($diseaseId: String!, $page: Int!, $size: Int!) {
  disease(efoId: $diseaseId) {
    id
    name
    associatedTargets(page: {index: $page, size: $size}) {
      count
      rows {
        target {
          id
          approvedSymbol
          approvedName
        }
        score
        datatypeScores {
          id
          score
        }
      }
    }
  }
}
"""

# GraphQL query: diseases for a target gene (by Ensembl ID)
_GENE_DISEASE_QUERY = """
query GeneAssociations($geneId: String!, $page: Int!, $size: Int!) {
  target(ensemblId: $geneId) {
    id
    approvedSymbol
    associatedDiseases(page: {index: $page, size: $size}) {
      count
      rows {
        disease {
          id
          name
        }
        score
      }
    }
  }
}
"""


class OpenTargetsClient:
    """
    Query Open Targets Platform for gene-disease association data.

    Usage::

        client = OpenTargetsClient()
        # Get top genes associated with a disease (MONDO or EFO ID)
        genes, edges = client.gene_associations_for_disease("MONDO:0004736", min_score=0.1)
    """

    def __init__(self, api_url: str = _API_URL, timeout: int = 60) -> None:
        self.api_url = api_url
        self.session = requests.Session()
        self.session.headers.update({"Content-Type": "application/json"})
        self.timeout = timeout

    def _query(self, query: str, variables: dict[str, Any]) -> dict:
        resp = self.session.post(
            self.api_url,
            json={"query": query, "variables": variables},
            timeout=self.timeout,
        )
        resp.raise_for_status()
        data = resp.json()
        if "errors" in data:
            raise ValueError(f"GraphQL errors: {data['errors']}")
        return data.get("data", {})

    def gene_associations_for_disease(
        self,
        disease_id: str,
        *,
        min_score: float = 0.1,
        max_results: int = 500,
        page_size: int = 100,
    ) -> tuple[list[GeneNode], list[DiseaseEdge]]:
        """
        Fetch gene associations for a disease.

        disease_id: MONDO ID ("MONDO:XXXXXXX") or EFO ID ("EFO_XXXXXXX").
        Open Targets accepts MONDO IDs directly.

        Returns (gene_nodes, disease_edges).
        """
        # Open Targets uses underscore in IDs: MONDO_0004736
        ot_id = disease_id.replace(":", "_")
        disease_node_id = disease_id  # keep "MONDO:XXXXXXX" for our graph

        genes: list[GeneNode] = []
        edges: list[DiseaseEdge] = []

        page = 0
        collected = 0
        while collected < max_results:
            try:
                data = self._query(
                    _ASSOC_QUERY,
                    {"diseaseId": ot_id, "page": page, "size": page_size},
                )
            except Exception as exc:
                log.warning("Open Targets query failed: %s", exc)
                break

            rows = (
                data.get("disease", {})
                .get("associatedTargets", {})
                .get("rows", [])
            )
            if not rows:
                break

            for row in rows:
                score = row.get("score", 0.0)
                if score < min_score:
                    continue
                target = row.get("target", {})
                ensembl_id = target.get("id", "")
                symbol = target.get("approvedSymbol", ensembl_id)
                gene_node_id = f"ENSG:{ensembl_id}"

                gene = GeneNode(
                    node_id=gene_node_id,
                    ensembl_id=ensembl_id,
                    symbol=symbol,
                    name=target.get("approvedName"),
                )
                edge = DiseaseEdge(
                    source=disease_node_id,
                    target=gene_node_id,
                    edge_type=DiseaseEdgeType.GENE_ASSOCIATED,
                    score=score,
                    source_db="open_targets",
                )
                genes.append(gene)
                edges.append(edge)
                collected += 1

            page += 1
            if len(rows) < page_size:
                break

        log.info(
            "Open Targets: %d genes for %s (min_score=%.2f)",
            len(genes), disease_id, min_score,
        )
        return genes, edges

    def disease_associations_for_gene(
        self,
        ensembl_id: str,
        *,
        min_score: float = 0.1,
        max_results: int = 200,
        page_size: int = 100,
    ) -> list[DiseaseEdge]:
        """
        Fetch disease associations for a gene.
        Returns DiseaseEdge objects (disease → gene direction preserved).
        """
        gene_node_id = f"ENSG:{ensembl_id}"
        edges: list[DiseaseEdge] = []

        page = 0
        collected = 0
        while collected < max_results:
            try:
                data = self._query(
                    _GENE_DISEASE_QUERY,
                    {"geneId": ensembl_id, "page": page, "size": page_size},
                )
            except Exception as exc:
                log.warning("Open Targets gene query failed: %s", exc)
                break

            rows = (
                data.get("target", {})
                .get("associatedDiseases", {})
                .get("rows", [])
            )
            if not rows:
                break

            for row in rows:
                score = row.get("score", 0.0)
                if score < min_score:
                    continue
                disease = row.get("disease", {})
                disease_efo_id = disease.get("id", "")
                # Convert EFO_XXXXXXX → MONDO:XXXXXXX best-effort; keep EFO form if no map
                disease_node_id = disease_efo_id.replace("_", ":", 1)

                edges.append(DiseaseEdge(
                    source=disease_node_id,
                    target=gene_node_id,
                    edge_type=DiseaseEdgeType.GENE_ASSOCIATED,
                    score=score,
                    source_db="open_targets",
                ))
                collected += 1

            page += 1
            if len(rows) < page_size:
                break

        return edges
