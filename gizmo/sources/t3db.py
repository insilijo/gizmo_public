"""
T3DB — Toxic Exposome Database loader.

License: T3DB data is freely available (no explicit open-source license;
non-commercial academic use). https://t3db.ca

T3DB provides:
  - ~3,700 toxin records (environmental, food, drug, industrial, …)
  - Toxin identifiers: T3D ID, HMDB ID, ChEBI ID, CAS number, PubChem CID
  - Target proteins / genes (Uniprot → HGNC symbol)
  - Mechanism of action, toxin action type
  - Exposure routes, LD50 values

Download files (CSV):
  toxin_list.csv       — master chemical list with cross-references
  toxin_targets.csv    — toxin → protein target links

These are small files (~10 MB total) downloaded from t3db.ca.

Matching strategy:
  - Toxin → graph metabolite: HMDB ID → MetaboliteNode.hmdb_id (exact)
    fallback to ChEBI ID, CAS, PubChem CID
  - Target → graph gene: gene symbol → GeneNode.symbol

Usage::

    from gizmo.sources.t3db import T3DBClient

    client = T3DBClient(cache_dir="data/raw/t3db")
    client.download()
    edges = client.load_tox_edges(mg)
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import Optional
from urllib.request import urlretrieve

from gizmo.schema import ToxEdge, ToxEdgeType

log = logging.getLogger(__name__)

_BASE = "https://t3db.ca/"
_FILES = {
    "toxins":  "toxins.csv",
    "targets": "targets.csv",
}

# Fallback URLs used by some T3DB mirror releases
_FALLBACK = {
    "toxins":  "http://www.t3db.ca/toxins.csv",
    "targets": "http://www.t3db.ca/toxin_target_links.csv",
}


class T3DBClient:
    """
    Downloads and parses T3DB toxin-target data.

    Parameters
    ----------
    cache_dir : local directory for downloaded CSV files
    """

    def __init__(self, cache_dir: str | Path = "data/raw/t3db") -> None:
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Download
    # ------------------------------------------------------------------

    def download(self, force: bool = False) -> None:
        """Download T3DB CSV files if not already cached."""
        for key, fname in _FILES.items():
            dest = self.cache_dir / fname
            if dest.exists() and not force:
                log.info("T3DB cache hit: %s", dest)
                continue
            for base in (_BASE, _FALLBACK.get(key, "")):
                if not base:
                    continue
                url = base.rstrip("/") + "/" + fname
                try:
                    log.info("Downloading T3DB %s → %s", url, dest)
                    urlretrieve(url, dest)
                    break
                except Exception as exc:
                    log.debug("T3DB download failed (%s): %s", url, exc)

    # ------------------------------------------------------------------
    # Graph loading
    # ------------------------------------------------------------------

    def load_tox_edges(self, mg) -> list[ToxEdge]:
        """
        Parse T3DB CSVs and return ToxEdge list connecting toxin metabolite
        nodes to gene nodes in the graph.

        Only toxin–target pairs where both the toxin and the target gene are
        already present in the graph are returned.
        """
        toxin_map = self._build_toxin_map(mg)
        gene_map  = self._build_gene_map(mg)

        if not toxin_map:
            log.warning("T3DB: no toxin nodes matched in graph. Run CompTox first to populate cas_id/hmdb_id.")
            return []

        target_path = self.cache_dir / _FILES["targets"]
        if not target_path.exists():
            log.warning("T3DB: targets file not found (%s). Run .download() first.", target_path)
            return []

        edges: list[ToxEdge] = []
        skipped = 0

        with open(target_path, encoding="utf-8-sig", newline="") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                t3d_id = (row.get("T3D ID") or row.get("Toxin ID") or "").strip()
                toxin_node = toxin_map.get(t3d_id)
                if not toxin_node:
                    skipped += 1
                    continue

                gene_sym = (
                    row.get("Gene Symbol") or row.get("Gene Name") or
                    row.get("HGNC Symbol") or ""
                ).strip().upper()
                gene_node = gene_map.get(gene_sym)
                if not gene_node:
                    skipped += 1
                    continue

                tox_action = (
                    row.get("Toxin Action") or row.get("Mechanism of Action") or ""
                ).strip() or None

                ld50_s = (row.get("LD50") or "").strip()
                ld50: Optional[float] = None
                ld50_units: Optional[str] = None
                if ld50_s:
                    try:
                        ld50 = float(ld50_s.split()[0])
                        ld50_units = " ".join(ld50_s.split()[1:]) or "mg/kg"
                    except (ValueError, IndexError):
                        pass

                edges.append(ToxEdge(
                    source=toxin_node,
                    target=gene_node,
                    edge_type=ToxEdgeType.TOX_GENE,
                    effect_type=tox_action,
                    organism="Homo sapiens",
                    assay_endpoint="LD50" if ld50 else None,
                    assay_value=ld50,
                    assay_units=ld50_units,
                    source_db="t3db",
                ))

        log.info("T3DB: %d tox-gene edges | %d skipped (no graph match)", len(edges), skipped)
        return edges

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _build_toxin_map(self, mg) -> dict[str, str]:
        """
        Returns {T3D_ID → graph_node_id} by reading toxin_list.csv and
        matching each toxin's HMDB, ChEBI, CAS, or PubChem CID to the graph.
        """
        toxin_path = self.cache_dir / _FILES["toxins"]
        if not toxin_path.exists():
            log.warning("T3DB: toxins file not found (%s). Run .download() first.", toxin_path)
            return {}

        # Build reverse lookup indexes from graph
        hmdb_to_node:   dict[str, str] = {}
        chebi_to_node:  dict[str, str] = {}
        cas_to_node:    dict[str, str] = {}
        pubchem_to_node: dict[str, str] = {}

        for nid, attrs in mg.graph.nodes(data=True):
            if attrs.get("node_type") != "metabolite":
                continue
            if attrs.get("hmdb_id"):
                hmdb_to_node[_norm_hmdb(attrs["hmdb_id"])] = nid
            if attrs.get("chebi_id"):
                chebi_to_node[attrs["chebi_id"]] = nid
            if attrs.get("cas_id"):
                cas_to_node[attrs["cas_id"]] = nid
            if attrs.get("pubchem_cid"):
                pubchem_to_node[str(attrs["pubchem_cid"])] = nid

        toxin_map: dict[str, str] = {}
        with open(toxin_path, encoding="utf-8-sig", newline="") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                t3d_id = (row.get("T3D ID") or row.get("Toxin ID") or "").strip()
                if not t3d_id:
                    continue

                node_id: Optional[str] = None

                hmdb = _norm_hmdb(row.get("HMDB ID") or row.get("HMDB_ID") or "")
                if hmdb:
                    node_id = hmdb_to_node.get(hmdb)

                if not node_id:
                    chebi = (row.get("ChEBI ID") or row.get("CHEBI_ID") or "").strip()
                    if chebi and not chebi.startswith("CHEBI:"):
                        chebi = f"CHEBI:{chebi}"
                    if chebi:
                        node_id = chebi_to_node.get(chebi)

                if not node_id:
                    cas = (row.get("CAS RN") or row.get("CASRN") or "").strip()
                    if cas:
                        node_id = cas_to_node.get(cas)

                if not node_id:
                    pc = (row.get("PubChem CID") or row.get("PUBCHEM_ID") or "").strip()
                    if pc:
                        node_id = pubchem_to_node.get(pc)

                if node_id:
                    toxin_map[t3d_id] = node_id

        log.info("T3DB toxin map: %d / %d toxins matched to graph nodes",
                 len(toxin_map), len(toxin_map) + 0)
        return toxin_map

    def _build_gene_map(self, mg) -> dict[str, str]:
        """Returns {UPPER_SYMBOL → node_id} for all gene nodes in the graph."""
        gene_map: dict[str, str] = {}
        for nid, attrs in mg.graph.nodes(data=True):
            if attrs.get("node_type") == "gene":
                sym = (attrs.get("symbol") or "").upper()
                if sym:
                    gene_map[sym] = nid
        return gene_map


def _norm_hmdb(s: str) -> str:
    s = (s or "").strip().upper()
    if s.startswith("HMDB"):
        return f"HMDB{s[4:].zfill(7)}"
    return s
