"""
ChEMBL drug-node builder.

Queries ChEMBL for approved and clinical-stage drugs targeting genes already
in the graph, creating first-class DrugNode + DrugEdge records.

License: ChEMBL data is released under CC BY-SA 3.0.
         https://www.ebi.ac.uk/chembl/

Each gene node is queried by HGNC symbol via the ChEMBL REST API.
Results are cached per gene.  Only drugs at max_phase ≥ min_phase are kept.
"""

from __future__ import annotations

import json
import logging
import time
from pathlib import Path
from typing import Optional
from urllib.request import urlopen
from urllib.error import URLError

from gizmo.schema import DrugEdge, DrugNode

log = logging.getLogger(__name__)

_CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"


class ChEMBLDrugClient:
    """
    Build DrugNode + DrugEdge objects from ChEMBL compound-target data.

    Usage::

        client = ChEMBLDrugClient(cache_dir="data/raw/chembl_drugs")
        nodes, edges = client.load_for_graph(mg, min_phase=2)
    """

    def __init__(self, cache_dir: str | Path = "data/raw/chembl_drugs") -> None:
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def _cache_path(self, symbol: str) -> Path:
        safe = symbol.replace("/", "_").replace(":", "_")
        return self.cache_dir / f"{safe}.json"

    def _fetch_target(self, symbol: str) -> list[dict]:
        """Return ChEMBL target JSON for a gene symbol (cached)."""
        cache = self._cache_path(symbol)
        if cache.exists():
            return json.loads(cache.read_text())

        url = (
            f"{_CHEMBL_BASE}/target.json"
            f"?target_synonym={symbol}&organism=Homo+sapiens&limit=5"
        )
        try:
            with urlopen(url, timeout=15) as resp:
                data = json.loads(resp.read())
            targets = data.get("targets", [])
        except (URLError, Exception) as exc:
            log.debug("ChEMBL target fetch failed for %s: %s", symbol, exc)
            targets = []

        cache.write_text(json.dumps(targets))
        return targets

    def _fetch_activities(self, chembl_target_id: str) -> list[dict]:
        """Return approved/clinical compounds for a ChEMBL target ID (cached)."""
        cache = self.cache_dir / f"act_{chembl_target_id}.json"
        if cache.exists():
            return json.loads(cache.read_text())

        url = (
            f"{_CHEMBL_BASE}/mechanism.json"
            f"?target_chembl_id={chembl_target_id}&limit=200"
        )
        try:
            with urlopen(url, timeout=15) as resp:
                data = json.loads(resp.read())
            records = data.get("mechanisms", [])
        except (URLError, Exception) as exc:
            log.debug("ChEMBL mechanism fetch failed for %s: %s", chembl_target_id, exc)
            records = []

        cache.write_text(json.dumps(records))
        return records

    def _fetch_molecule(self, chembl_id: str) -> dict:
        """Return molecule metadata (name, max_phase, ATC codes)."""
        cache = self.cache_dir / f"mol_{chembl_id}.json"
        if cache.exists():
            return json.loads(cache.read_text())

        url = f"{_CHEMBL_BASE}/molecule/{chembl_id}.json"
        try:
            with urlopen(url, timeout=10) as resp:
                data = json.loads(resp.read())
        except (URLError, Exception) as exc:
            log.debug("ChEMBL molecule fetch failed for %s: %s", chembl_id, exc)
            data = {}

        cache.write_text(json.dumps(data))
        return data

    def load_for_graph(
        self,
        mg,
        min_phase: int = 2,
        sleep_s: float = 0.3,
    ) -> tuple[list[DrugNode], list[DrugEdge]]:
        """
        Query ChEMBL for drugs targeting gene nodes already in ``mg``.

        Parameters
        ----------
        mg        : GizmoGraph to scan for gene nodes
        min_phase : only include drugs at this clinical phase or higher (2 = Phase 2+)
        sleep_s   : politeness delay between API calls

        Returns
        -------
        (drug_nodes, drug_edges)
        """
        gene_symbols: list[str] = []
        for nid, attrs in mg.graph.nodes(data=True):
            if attrs.get("node_type") == "gene":
                sym = attrs.get("symbol")
                if sym:
                    gene_symbols.append((nid, sym))

        n_genes = len(gene_symbols)
        log.info("ChEMBL drugs: querying targets for %d gene symbols …", n_genes)

        nodes: list[DrugNode] = []
        edges: list[DrugEdge] = []
        seen_drugs: dict[str, str] = {}  # chembl_id → node_id

        for gi, (gene_nid, symbol) in enumerate(gene_symbols):
            if (gi + 1) % 50 == 0 or gi == 0:
                log.info("  [%d/%d] gene %s — %d drugs found so far",
                         gi + 1, n_genes, symbol, len(nodes))
            targets = self._fetch_target(symbol)
            time.sleep(sleep_s)

            for target in targets:
                target_id = target.get("target_chembl_id")
                if not target_id:
                    continue

                mechanisms = self._fetch_activities(target_id)
                time.sleep(sleep_s)

                for mech in mechanisms:
                    mol_id = mech.get("molecule_chembl_id")
                    if not mol_id:
                        continue

                    # Fetch molecule to get max_phase. ChEMBL recently began
                    # returning the field as a float-string (e.g. "4.0"),
                    # which used to parse as an int — cast via float first.
                    mol = self._fetch_molecule(mol_id)
                    try:
                        max_phase = int(float(mol.get("max_phase") or 0))
                    except (TypeError, ValueError):
                        max_phase = 0
                    if max_phase < min_phase:
                        continue
                    time.sleep(sleep_s * 0.5)

                    drug_node_id = f"CHEMBL:{mol_id}"

                    if mol_id not in seen_drugs:
                        pref_name = mol.get("pref_name") or mol_id
                        raw_syns = mol.get("molecule_synonyms") or []
                        synonyms = []
                        for s in raw_syns:
                            if isinstance(s, dict):
                                v = s.get("molecule_synonym", "")
                            else:
                                v = str(s) if s else ""
                            if v:
                                synonyms.append(v)
                        synonyms = synonyms[:5]

                        raw_atc = mol.get("atc_classifications") or []
                        atc_codes = []
                        for a in raw_atc:
                            if isinstance(a, dict):
                                v = a.get("level5", "")
                            else:
                                v = str(a) if a else ""
                            if v:
                                atc_codes.append(v)

                        nodes.append(DrugNode(
                            node_id=drug_node_id,
                            chembl_id=mol_id,
                            name=pref_name,
                            synonyms=synonyms,
                            max_phase=max_phase,
                            mechanism=mech.get("mechanism_of_action"),
                            atc_codes=atc_codes,
                        ))
                        seen_drugs[mol_id] = drug_node_id

                    edges.append(DrugEdge(
                        source=drug_node_id,
                        target=gene_nid,
                        edge_type="drug_target",
                        mechanism=mech.get("action_type"),
                        max_phase=max_phase,
                        source_db="chembl",
                    ))

        log.info(
            "ChEMBL drugs: %d drug nodes, %d drug-target edges (min_phase=%d).",
            len(nodes), len(edges), min_phase,
        )
        return nodes, edges
