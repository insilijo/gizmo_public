"""
Comparative Toxicogenomics Database (CTD) loader.

License: CTD data is freely available for non-commercial use.
https://ctdbase.org/about/legal.jsp

Key files (gzipped TSV, comment lines start with '#'):
  CTD_chem_gene_ixns.tsv.gz  — chemical–gene interactions
  CTD_chem_diseases.tsv.gz   — chemical–disease associations

Chemical matching strategy:
  1. CAS Registry Number → MetaboliteNode.cas_id  (exact)
  2. CAS → MetaNetX (chem_xref "cas:" prefix) → ChEBI → graph node
  3. MeSH Chemical ID → skip (no standard cross-ref in graph yet)

Only Homo sapiens interactions are loaded by default; pass
``organisms=None`` to load all species.

Usage::

    from gizmo.sources.ctd import CTDClient

    client = CTDClient(cache_dir="data/raw/ctd")
    client.download()
    client.build_cas_index(mg, chem_xref_path="data/raw/metanetx/chem_xref.tsv")
    gene_edges = client.load_chem_gene_interactions(mg)
    disease_edges = client.load_chem_disease_associations(mg)
"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path
from typing import Iterable, Optional
from urllib.request import urlretrieve

from gizmo.schema import ToxEdge, ToxEdgeType

log = logging.getLogger(__name__)

_BASE = "https://ctdbase.org/reports/"
_FILES = {
    "chem_gene":    "CTD_chem_gene_ixns.tsv.gz",
    "chem_disease": "CTD_chem_diseases.tsv.gz",
}

# CTD column indices (0-based) — stable across releases
_CGI_COLS = {
    "chem_name": 0, "chem_id": 1, "cas_rn": 2,
    "gene_symbol": 3, "gene_id": 4, "gene_forms": 5,
    "organism": 6, "organism_id": 7,
    "interaction": 8, "interaction_actions": 9, "pubmed_ids": 10,
}
_CD_COLS = {
    "chem_name": 0, "chem_id": 1, "cas_rn": 2,
    "disease_name": 3, "disease_id": 4,
    "direct_evidence": 5, "inference_gene": 6, "inference_score": 7,
    "omim_ids": 8, "pubmed_ids": 9,
}


class CTDClient:
    """
    Downloads and parses CTD chemical–gene and chemical–disease files.

    Parameters
    ----------
    cache_dir  : local directory for downloaded files
    organisms  : set of organism names to include (default: {"Homo sapiens"}).
                 Pass None to include all organisms.
    """

    def __init__(
        self,
        cache_dir: str | Path = "data/raw/ctd",
        organisms: Optional[set[str]] = None,
    ) -> None:
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.organisms: Optional[set[str]] = (
            {"Homo sapiens"} if organisms is None else organisms
        )
        # CAS → graph node_id  (populated by build_cas_index)
        self._cas_to_node: dict[str, str] = {}
        # gene symbol (upper) → graph node_id
        self._symbol_to_node: dict[str, str] = {}
        # MONDO/disease node_id → set of MeSH IDs stored as xref_mesh
        self._mesh_to_disease: dict[str, str] = {}

    # ------------------------------------------------------------------
    # Download
    # ------------------------------------------------------------------

    def download(self, force: bool = False) -> None:
        """Download CTD TSV files if not already cached."""
        for key, fname in _FILES.items():
            dest = self.cache_dir / fname
            if dest.exists() and not force:
                log.info("CTD cache hit: %s", dest)
                continue
            url = _BASE + fname
            log.info("Downloading CTD %s → %s", url, dest)
            urlretrieve(url, dest)

    # ------------------------------------------------------------------
    # Index building
    # ------------------------------------------------------------------

    def build_cas_index(
        self,
        mg,
        chem_xref_path: Optional[str | Path] = None,
    ) -> int:
        """
        Build CAS → node_id lookup from the graph and optionally MetaNetX.

        Pass ``chem_xref_path`` to MetaNetX chem_xref.tsv to extend coverage
        via CAS → MNX → ChEBI → node_id for metabolites that don't yet have
        a ``cas_id`` attribute.

        Returns the number of CAS entries indexed.
        """
        # Pass 1: graph nodes that already have cas_id
        for nid, attrs in mg.graph.nodes(data=True):
            if attrs.get("node_type") != "metabolite":
                continue
            cas = (attrs.get("cas_id") or "").strip()
            if cas:
                self._cas_to_node[cas] = nid

        # Pass 2: MetaNetX chem_xref → CAS → ChEBI → node_id
        if chem_xref_path:
            self._cas_to_node.update(
                self._cas_from_metanetx(mg, Path(chem_xref_path))
            )

        # Gene symbol index
        for nid, attrs in mg.graph.nodes(data=True):
            if attrs.get("node_type") == "gene":
                sym = (attrs.get("symbol") or "").upper()
                if sym:
                    self._symbol_to_node[sym] = nid

        # Disease MeSH xref index
        for nid, attrs in mg.graph.nodes(data=True):
            if attrs.get("node_type") != "disease":
                continue
            mesh_list = attrs.get("xref_mesh") or []
            if isinstance(mesh_list, str):
                mesh_list = [mesh_list]
            for m in mesh_list:
                # normalise "MeSH:D001234" → "D001234"
                mesh_id = m.split(":")[-1].strip()
                if mesh_id:
                    self._mesh_to_disease[mesh_id] = nid

        log.info(
            "CTD index: %d CAS, %d gene symbols, %d MeSH disease IDs",
            len(self._cas_to_node), len(self._symbol_to_node),
            len(self._mesh_to_disease),
        )
        return len(self._cas_to_node)

    def _cas_from_metanetx(self, mg, xref_path: Path) -> dict[str, str]:
        """Build CAS → node_id via MetaNetX chem_xref cas: entries."""
        from gizmo.sources.metanetx import _mnx_header_info

        if not xref_path.exists():
            log.warning("CTD: MetaNetX chem_xref not found at %s", xref_path)
            return {}

        # Build ChEBI → node_id from graph
        chebi_to_node: dict[str, str] = {}
        for nid, attrs in mg.graph.nodes(data=True):
            if attrs.get("node_type") == "metabolite":
                ch = attrs.get("chebi_id")
                if ch:
                    chebi_to_node[ch] = nid

        cols, start = _mnx_header_info(xref_path)
        mnx_idx = cols.index("ID") if "ID" in cols else 1

        # Two passes: chebi: → MNX, then cas: → MNX
        mnx_to_chebi: dict[str, str] = {}
        cas_to_mnx:   dict[str, str] = {}

        import re
        with open(xref_path, encoding="utf-8") as fh:
            for lineno, line in enumerate(fh):
                if lineno < start:
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= mnx_idx:
                    continue
                source = parts[0]
                mnx_id = parts[mnx_idx]
                if source.startswith("chebi:"):
                    chebi_id = "CHEBI:" + re.sub(r"^chebi:(?:CHEBI:)?", "", source)
                    mnx_to_chebi.setdefault(mnx_id, chebi_id)
                elif source.startswith("cas:"):
                    cas = source[4:]
                    cas_to_mnx.setdefault(cas, mnx_id)

        result: dict[str, str] = {}
        for cas, mnx_id in cas_to_mnx.items():
            chebi = mnx_to_chebi.get(mnx_id)
            if chebi:
                node_id = chebi_to_node.get(chebi)
                if node_id:
                    result[cas] = node_id
                    # Back-fill cas_id on the graph node
                    mg.graph.nodes[node_id]["cas_id"] = cas

        log.info("CTD MetaNetX CAS index: %d entries", len(result))
        return result

    # ------------------------------------------------------------------
    # Interaction loading
    # ------------------------------------------------------------------

    def load_chem_gene_interactions(
        self,
        mg,
        direct_evidence_only: bool = False,
    ) -> list[ToxEdge]:
        """
        Parse CTD_chem_gene_ixns.tsv.gz and return ToxEdge list.

        Parameters
        ----------
        direct_evidence_only : if True, only load rows with curated evidence
                               (skips inferred interactions)
        """
        path = self.cache_dir / _FILES["chem_gene"]
        if not path.exists():
            log.warning("CTD: chem_gene file not found (%s). Run .download() first.", path)
            return []

        edges: list[ToxEdge] = []
        skipped_no_node = 0
        skipped_organism = 0

        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 10:
                    continue

                organism = parts[_CGI_COLS["organism"]].strip()
                if self.organisms and organism not in self.organisms:
                    skipped_organism += 1
                    continue

                cas = parts[_CGI_COLS["cas_rn"]].strip()
                chem_node = self._cas_to_node.get(cas)
                if not chem_node:
                    skipped_no_node += 1
                    continue

                gene_sym = parts[_CGI_COLS["gene_symbol"]].strip().upper()
                gene_node = self._symbol_to_node.get(gene_sym)
                if not gene_node:
                    skipped_no_node += 1
                    continue

                interaction_actions = parts[_CGI_COLS["interaction_actions"]].strip()
                pubmed_ids = parts[_CGI_COLS["pubmed_ids"]].strip()
                ref_count = len(pubmed_ids.split("|")) if pubmed_ids else 0

                edges.append(ToxEdge(
                    source=chem_node,
                    target=gene_node,
                    edge_type=ToxEdgeType.TOX_GENE,
                    effect_type=interaction_actions or None,
                    organism=organism or None,
                    reference_count=ref_count or None,
                    source_db="ctd",
                ))

        log.info(
            "CTD chem-gene: %d edges | skipped %d (no node), %d (organism filter)",
            len(edges), skipped_no_node, skipped_organism,
        )
        return edges

    def load_chem_disease_associations(
        self,
        mg,
        direct_evidence_only: bool = True,
        min_inference_score: float = 0.0,
    ) -> list[ToxEdge]:
        """
        Parse CTD_chem_diseases.tsv.gz and return ToxEdge list.

        Parameters
        ----------
        direct_evidence_only  : if True, only return rows with DirectEvidence
                                (marker/mechanism or therapeutic), not inferred
        min_inference_score   : minimum inference score for inferred associations
                                (only used when direct_evidence_only=False)
        """
        path = self.cache_dir / _FILES["chem_disease"]
        if not path.exists():
            log.warning("CTD: chem_disease file not found (%s). Run .download() first.", path)
            return []

        edges: list[ToxEdge] = []
        skipped_no_node = 0

        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 6:
                    continue

                cas = parts[_CD_COLS["cas_rn"]].strip()
                chem_node = self._cas_to_node.get(cas)
                if not chem_node:
                    skipped_no_node += 1
                    continue

                direct = parts[_CD_COLS["direct_evidence"]].strip()
                if direct_evidence_only and not direct:
                    continue

                if not direct:
                    score_s = parts[_CD_COLS["inference_score"]].strip()
                    try:
                        score = float(score_s)
                    except (ValueError, TypeError):
                        score = 0.0
                    if score < min_inference_score:
                        continue

                disease_id_raw = parts[_CD_COLS["disease_id"]].strip()
                # CTD disease IDs look like "MESH:D012345" or "OMIM:123456"
                mesh_id = disease_id_raw.split(":")[-1].strip()
                disease_node = self._mesh_to_disease.get(mesh_id)
                if not disease_node:
                    skipped_no_node += 1
                    continue

                pubmed_ids = parts[_CD_COLS["pubmed_ids"]].strip() if len(parts) > 9 else ""
                ref_count = len(pubmed_ids.split("|")) if pubmed_ids else 0

                edges.append(ToxEdge(
                    source=chem_node,
                    target=disease_node,
                    edge_type=ToxEdgeType.TOX_DISEASE,
                    direct_evidence=direct or None,
                    reference_count=ref_count or None,
                    source_db="ctd",
                ))

        log.info(
            "CTD chem-disease: %d edges | skipped %d (no node match)",
            len(edges), skipped_no_node,
        )
        return edges
