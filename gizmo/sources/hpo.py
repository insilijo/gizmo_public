"""
Human Phenotype Ontology (HPO) client.

License: CC BY 4.0 — https://hpo.jax.org/

Downloads:
  hp.obo                  HPO term definitions
  phenotype_to_genes.txt  HPO term → gene associations (JAX)
  phenotype.hpoa          HPO term → disease associations (OMIM/ORPHA/DECIPHER)

Nodes:   PhenotypeNode  (HP:XXXXXXX)
Edges:   PhenotypeEdge
           phenotype → disease   (phenotype_disease)
           phenotype → gene      (phenotype_gene)
"""

from __future__ import annotations

import csv
import io
import logging
from pathlib import Path
from typing import Optional
from urllib.request import urlretrieve

import obonet

from gizmo.schema import PhenotypeEdge, PhenotypeNode

log = logging.getLogger(__name__)

_HPO_OBO_URL  = "https://purl.obolibrary.org/obo/hp.obo"
_HPO_GENE_URL = (
    "https://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt"
)
_HPO_DIS_URL  = (
    "https://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa"
)

# Root term for "Abnormality of metabolism or catabolism"
_METABOLIC_ROOT = "HP:0001939"


class HPOClient:
    """
    Parse HPO OBO + annotation files into PhenotypeNode / PhenotypeEdge objects.

    Usage::

        client = HPOClient(cache_dir="data/raw/hpo")
        client.download()
        phenotypes = client.load_phenotypes()
        gene_edges = client.load_gene_edges()
        disease_edges = client.load_disease_edges(mondo_id_map)
    """

    def __init__(self, cache_dir: str | Path = "data/raw/hpo") -> None:
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self._obo_graph = None

    @property
    def obo_path(self) -> Path:
        return self.cache_dir / "hp.obo"

    @property
    def gene_path(self) -> Path:
        return self.cache_dir / "phenotype_to_genes.txt"

    @property
    def disease_path(self) -> Path:
        return self.cache_dir / "phenotype.hpoa"

    def download(self, force: bool = False) -> None:
        for url, path in [
            (_HPO_OBO_URL, self.obo_path),
            (_HPO_GENE_URL, self.gene_path),
            (_HPO_DIS_URL, self.disease_path),
        ]:
            if path.exists() and not force:
                log.info("HPO cache hit: %s", path.name)
                continue
            log.info("Downloading %s → %s", url.split("/")[-1], path)
            urlretrieve(url, path)

    def _load_obo(self):
        if self._obo_graph is None:
            if not self.obo_path.exists():
                raise FileNotFoundError(
                    f"HPO OBO not found at {self.obo_path}. Run .download() first."
                )
            log.info("Parsing HPO OBO …")
            self._obo_graph = obonet.read_obo(str(self.obo_path))
        return self._obo_graph

    def load_phenotypes(self) -> list[PhenotypeNode]:
        """
        Parse all HP terms into PhenotypeNode objects.
        Sets is_metabolic=True for terms under HP:0001939.
        """
        g = self._load_obo()

        # Find metabolic subtree
        import networkx as nx
        metabolic_ids: set[str] = set()
        if _METABOLIC_ROOT in g:
            metabolic_ids = nx.ancestors(g, _METABOLIC_ROOT)  # type: ignore
            metabolic_ids.add(_METABOLIC_ROOT)

        nodes: list[PhenotypeNode] = []
        for nid, data in g.nodes(data=True):
            if not nid.startswith("HP:"):
                continue
            if data.get("is_obsolete"):
                continue
            name = data.get("name", "")
            if not name:
                continue

            synonyms: list[str] = []
            for syn in data.get("synonym", []):
                val = syn.strip('"').split('"')[0] if '"' in syn else syn
                synonyms.append(val.strip())

            nodes.append(PhenotypeNode(
                node_id=nid,
                hpo_id=nid,
                name=name,
                definition=_strip_def(data.get("def", "")),
                synonyms=synonyms,
                is_metabolic=nid in metabolic_ids,
            ))

        log.info("Loaded %d HPO phenotype nodes.", len(nodes))
        return nodes

    def load_gene_edges(self) -> list[PhenotypeEdge]:
        """
        Parse phenotype_to_genes.txt → PhenotypeEdge(phenotype → gene).
        Gene target is stored as "symbol:{SYMBOL}" to match GeneMapper convention.
        """
        if not self.gene_path.exists():
            raise FileNotFoundError(
                f"HPO gene annotations not found at {self.gene_path}. Run .download() first."
            )

        edges: list[PhenotypeEdge] = []
        with open(self.gene_path, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 4:
                    continue
                # Format: hpo_id, hpo_name, ncbi_gene_id, gene_symbol
                hpo_id     = parts[0].strip()
                gene_sym   = parts[3].strip()
                if not hpo_id.startswith("HP:") or not gene_sym:
                    continue
                edges.append(PhenotypeEdge(
                    source=hpo_id,
                    target=f"symbol:{gene_sym}",
                    edge_type="phenotype_gene",
                    source_db="hpo",
                ))

        log.info("Loaded %d HPO–gene edges.", len(edges))
        return edges

    def load_disease_edges(
        self,
        mondo_id_map: dict[str, str] | None = None,
    ) -> list[PhenotypeEdge]:
        """
        Parse phenotype.hpoa → PhenotypeEdge(phenotype → disease).

        Parameters
        ----------
        mondo_id_map : {omim_id: mondo_id, orphanet_id: mondo_id} to resolve
                       OMIM:XXXXXX → MONDO:XXXXXXX.  Built from disease nodes
                       already in the graph if not provided.
        """
        if not self.disease_path.exists():
            raise FileNotFoundError(
                f"HPO disease annotations not found at {self.disease_path}. Run .download() first."
            )

        edges: list[PhenotypeEdge] = []
        seen: set[tuple] = set()

        with open(self.disease_path, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 4:
                    continue
                # Format: database_id, disease_name, qualifier, hpo_id, ...
                database_id = parts[0].strip()   # e.g. "OMIM:143100", "ORPHA:93"
                qualifier   = parts[2].strip()   # "NOT" means phenotype is excluded
                hpo_id      = parts[3].strip()

                if qualifier == "NOT" or not hpo_id.startswith("HP:"):
                    continue

                # Try to resolve to MONDO ID
                target = database_id
                if mondo_id_map:
                    target = mondo_id_map.get(database_id, database_id)

                key = (hpo_id, target)
                if key in seen:
                    continue
                seen.add(key)

                edges.append(PhenotypeEdge(
                    source=hpo_id,
                    target=target,
                    edge_type="phenotype_disease",
                    source_db="hpo",
                ))

        log.info("Loaded %d HPO–disease edges.", len(edges))
        return edges


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _strip_def(raw: str) -> Optional[str]:
    if not raw:
        return None
    # OBO def format: "text" [source]
    if raw.startswith('"'):
        end = raw.find('"', 1)
        if end > 0:
            return raw[1:end]
    return raw.strip()


def build_mondo_id_map(mg) -> dict[str, str]:
    """
    Build {omim_id: mondo_id, orphanet_id: mondo_id} from disease nodes in mg.
    Useful for resolving HPO disease annotations to the graph's MONDO IDs.
    """
    id_map: dict[str, str] = {}
    for nid, attrs in mg.graph.nodes(data=True):
        if attrs.get("node_type") != "disease":
            continue
        for xref in (attrs.get("xref_omim") or []):
            id_map[xref] = nid
        for xref in (attrs.get("xref_orphanet") or []):
            id_map[xref] = nid
            # Also accept "ORPHA:XXX" alias
            orpha_num = xref.replace("Orphanet:", "").replace("ORPHA:", "")
            id_map[f"ORPHA:{orpha_num}"] = nid
    return id_map
