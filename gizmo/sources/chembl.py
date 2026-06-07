"""
ChEMBL REST API client (CC BY-SA 4.0).

https://www.ebi.ac.uk/chembl/api/data/

Used by gizmo/actionability/ to fetch:
  - target ChEMBL IDs for gene symbols
  - approved and clinical drugs targeting each enzyme
  - mechanism of action per drug-target pair

All results are cached in memory for the lifetime of the client instance.
For persistent caching pass a ``cache_dir``.

Usage::

    from gizmo.sources.chembl import ChEMBLClient

    client = ChEMBLClient()
    targets = client.targets_for_gene("ALDH2")
    drugs   = client.drugs_for_target(targets[0]["target_chembl_id"], min_phase=2)
"""

from __future__ import annotations

import json
import logging
import time
from pathlib import Path
from typing import Any, Optional

import requests

log = logging.getLogger(__name__)

_BASE         = "https://www.ebi.ac.uk/chembl/api/data"
_PAGE_SIZE    = 100
_REQUEST_DELAY= 0.2   # seconds between requests


class ChEMBLClient:
    """
    Thin wrapper around the ChEMBL REST API.

    Parameters
    ----------
    cache_dir : if given, results are persisted to JSON files here
    """

    def __init__(self, cache_dir: Optional[str | Path] = None) -> None:
        self._session = requests.Session()
        self._session.headers.update({
            "Accept": "application/json",
            "User-Agent": "gizmo/1.0",
        })
        self._cache_dir = Path(cache_dir) if cache_dir else None
        self._mem: dict[str, Any] = {}    # in-memory cache

    # ------------------------------------------------------------------
    # Target lookup
    # ------------------------------------------------------------------

    def targets_for_gene(
        self,
        gene_symbol: str,
        organism: str = "Homo sapiens",
    ) -> list[dict]:
        """
        Return ChEMBL target records for a HGNC gene symbol.

        Each dict contains at minimum:
            target_chembl_id, pref_name, organism, target_type
        """
        cache_key = f"targets:{gene_symbol}"
        if cache_key in self._mem:
            return self._mem[cache_key]
        cached = self._load_cache(cache_key)
        if cached is not None:
            self._mem[cache_key] = cached
            return cached

        results: list[dict] = []
        # Search by gene symbol in target synonyms
        for search_field in ("target_synonym", "pref_name"):
            params = {
                "format":         "json",
                "limit":          20,
                f"{search_field}": gene_symbol,
            }
            if organism:
                params["organism"] = organism
            try:
                data = self._get(f"{_BASE}/target", params=params)
                hits = data.get("targets", [])
                for t in hits:
                    if t.get("target_type") in (
                        "SINGLE PROTEIN", "PROTEIN COMPLEX", "PROTEIN FAMILY"
                    ):
                        results.append({
                            "target_chembl_id": t["target_chembl_id"],
                            "pref_name":        t.get("pref_name", ""),
                            "organism":         t.get("organism", ""),
                            "target_type":      t.get("target_type", ""),
                        })
                if results:
                    break
            except Exception as exc:
                log.debug("ChEMBL target search failed (%s=%s): %s",
                          search_field, gene_symbol, exc)

        self._mem[cache_key] = results
        self._save_cache(cache_key, results)
        return results

    # ------------------------------------------------------------------
    # Drug / mechanism lookup
    # ------------------------------------------------------------------

    def drugs_for_target(
        self,
        target_chembl_id: str,
        min_phase: int = 0,
    ) -> list[dict]:
        """
        Return drug–mechanism records for a ChEMBL target ID.

        Each dict contains:
            molecule_chembl_id, molecule_name, max_phase,
            mechanism_of_action, action_type, direct_interaction
        """
        cache_key = f"drugs:{target_chembl_id}:{min_phase}"
        if cache_key in self._mem:
            return self._mem[cache_key]
        cached = self._load_cache(cache_key)
        if cached is not None:
            self._mem[cache_key] = cached
            return cached

        results: list[dict] = []
        offset = 0
        while True:
            try:
                data = self._get(f"{_BASE}/mechanism", params={
                    "format":             "json",
                    "target_chembl_id":   target_chembl_id,
                    "limit":              _PAGE_SIZE,
                    "offset":             offset,
                })
            except Exception as exc:
                log.debug("ChEMBL mechanisms failed (%s): %s", target_chembl_id, exc)
                break

            mechs = data.get("mechanisms", [])
            for m in mechs:
                phase = m.get("max_phase") or 0
                if phase >= min_phase:
                    results.append({
                        "molecule_chembl_id": m.get("molecule_chembl_id", ""),
                        "molecule_name":      m.get("molecule_name") or m.get("molecule_chembl_id", ""),
                        "max_phase":          phase,
                        "mechanism_of_action": m.get("mechanism_of_action", ""),
                        "action_type":        m.get("action_type", ""),
                        "direct_interaction": m.get("direct_interaction", True),
                    })

            page_meta = data.get("page_meta", {})
            if page_meta.get("next") is None:
                break
            offset += _PAGE_SIZE
            time.sleep(_REQUEST_DELAY)

        self._mem[cache_key] = results
        self._save_cache(cache_key, results)
        return results

    # ------------------------------------------------------------------
    # Toxicology / ADMET assay lookup
    # ------------------------------------------------------------------

    # Key ADMET endpoints to query — (assay name fragment, canonical label, units)
    _TOX_ENDPOINTS: list[tuple[str, str, str]] = [
        ("herg",        "hERG IC50",    "uM"),
        ("ld50",        "LD50",         "mg/kg"),
        ("ames",        "AMES",         ""),
        ("hepatotox",   "Hepatotoxicity", ""),
        ("cardiotox",   "Cardiotoxicity", ""),
        ("bbb",         "BBB Permeability", ""),
        ("cyp",         "CYP Inhibition", "uM"),
    ]

    def tox_assays_for_compound(
        self,
        molecule_chembl_id: str,
        endpoints: Optional[list[str]] = None,
    ) -> list[dict]:
        """
        Fetch ADMET / toxicology bioactivity records for a ChEMBL molecule ID.

        Parameters
        ----------
        molecule_chembl_id : e.g. "CHEMBL25"
        endpoints          : list of endpoint keywords to filter on
                             (e.g. ["herg", "ld50"]); None → all tox assays

        Returns a list of dicts with keys:
            assay_chembl_id, assay_description, standard_type,
            standard_value, standard_units, pchembl_value
        """
        cache_key = f"tox:{molecule_chembl_id}"
        if cache_key in self._mem:
            return self._mem[cache_key]
        cached = self._load_cache(cache_key)
        if cached is not None:
            self._mem[cache_key] = cached
            return cached

        results: list[dict] = []
        offset = 0
        while True:
            try:
                data = self._get(f"{_BASE}/activity", params={
                    "format":             "json",
                    "molecule_chembl_id": molecule_chembl_id,
                    "assay_type":         "A",   # A = ADMET
                    "limit":              _PAGE_SIZE,
                    "offset":             offset,
                })
            except Exception as exc:
                log.debug("ChEMBL tox assays failed (%s): %s", molecule_chembl_id, exc)
                break

            activities = data.get("activities", [])
            for a in activities:
                std_type  = (a.get("standard_type") or "").lower()
                assay_desc = (a.get("assay_description") or "").lower()

                # Filter to recognised tox endpoints
                keep = False
                if endpoints:
                    keep = any(ep.lower() in std_type or ep.lower() in assay_desc
                               for ep in endpoints)
                else:
                    keep = any(ep.lower() in std_type or ep.lower() in assay_desc
                               for ep, _, _ in self._TOX_ENDPOINTS)
                if not keep:
                    continue

                results.append({
                    "assay_chembl_id":  a.get("assay_chembl_id", ""),
                    "assay_description": a.get("assay_description", ""),
                    "standard_type":    a.get("standard_type", ""),
                    "standard_value":   a.get("standard_value"),
                    "standard_units":   a.get("standard_units", ""),
                    "pchembl_value":    a.get("pchembl_value"),
                })

            page_meta = data.get("page_meta", {})
            if page_meta.get("next") is None:
                break
            offset += _PAGE_SIZE
            time.sleep(_REQUEST_DELAY)

        self._mem[cache_key]  = results
        self._save_cache(cache_key, results)
        return results

    def tox_edges_for_graph(
        self,
        mg,
        min_phase: int = 0,
    ) -> list:
        """
        For every metabolite node with a known ChEMBL ID (dtxsid not required),
        fetch ADMET assay records and return ToxEdge objects pointing from the
        metabolite node to itself (self-loops carrying the assay annotation) or
        to a target gene node where the assay is gene-specific (e.g. hERG).

        In practice this emits a list of dicts suitable for graph annotation
        rather than formal ToxEdges, since ChEMBL ADMET data is compound-
        intrinsic rather than compound→target.

        Returns list of dicts:
            {node_id, assay_endpoint, assay_value, assay_units, chembl_id}
        """
        from gizmo.schema import ToxEdge, ToxEdgeType

        annotations = []
        herg_node = self._find_herg_gene_node(mg)

        for nid, attrs in mg.graph.nodes(data=True):
            if attrs.get("node_type") != "metabolite":
                continue
            chembl_id = attrs.get("chembl_id") or attrs.get("drug_chembl_id")
            if not chembl_id:
                continue

            records = self.tox_assays_for_compound(chembl_id)
            for r in records:
                std_type  = (r.get("standard_type") or "").lower()
                std_val   = r.get("standard_value")
                std_units = r.get("standard_units") or ""
                try:
                    val = float(std_val) if std_val is not None else None
                except (ValueError, TypeError):
                    val = None

                # hERG: add as compound→gene edge
                if "herg" in std_type and herg_node:
                    annotations.append(ToxEdge(
                        source=nid,
                        target=herg_node,
                        edge_type=ToxEdgeType.TOX_GENE,
                        effect_type="inhibition",
                        assay_endpoint="hERG IC50",
                        assay_value=val,
                        assay_units=std_units or "uM",
                        source_db="chembl_tox",
                    ))
                else:
                    # General ADMET: annotate node attribute
                    key = f"chembl_tox_{std_type.replace(' ', '_')}"
                    if val is not None:
                        mg.graph.nodes[nid][key] = val

        return annotations

    @staticmethod
    def _find_herg_gene_node(mg) -> Optional[str]:
        """Return the node_id for the hERG / KCNH2 gene if it exists in the graph."""
        for nid, attrs in mg.graph.nodes(data=True):
            if attrs.get("node_type") == "gene":
                sym = (attrs.get("symbol") or "").upper()
                if sym in ("KCNH2", "HERG", "ERG1"):
                    return nid
        return None

    def drugs_for_gene(
        self,
        gene_symbol: str,
        min_phase: int = 0,
        organism: str = "Homo sapiens",
    ) -> list[dict]:
        """
        Convenience: resolve gene symbol → target IDs → drugs.
        Returns deduplicated drug list across all matching targets.
        """
        targets  = self.targets_for_gene(gene_symbol, organism=organism)
        seen:    set[str] = set()
        drugs:   list[dict] = []
        for t in targets:
            for d in self.drugs_for_target(t["target_chembl_id"], min_phase=min_phase):
                cid = d["molecule_chembl_id"]
                if cid not in seen:
                    seen.add(cid)
                    d["target_chembl_id"] = t["target_chembl_id"]
                    d["target_name"]      = t["pref_name"]
                    drugs.append(d)
        return drugs

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _get(self, url: str, params: dict) -> dict:
        time.sleep(_REQUEST_DELAY)
        resp = self._session.get(url, params=params, timeout=30)
        resp.raise_for_status()
        return resp.json()

    def _cache_path(self, key: str) -> Optional[Path]:
        if self._cache_dir is None:
            return None
        safe = key.replace(":", "_").replace("/", "_")
        return self._cache_dir / f"chembl_{safe}.json"

    def _load_cache(self, key: str) -> Optional[Any]:
        p = self._cache_path(key)
        if p and p.exists():
            try:
                return json.loads(p.read_text())
            except Exception:
                pass
        return None

    def _save_cache(self, key: str, data: Any) -> None:
        p = self._cache_path(key)
        if p:
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(json.dumps(data, indent=2))
