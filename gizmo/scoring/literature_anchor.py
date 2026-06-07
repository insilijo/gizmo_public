"""Literature-anchored disease → metabolite harvester.

Uses NCBI PubMed E-utilities (public, license-clean) to find chemicals
co-cited with a disease in metabolomics / biomarker papers, then maps
those chemical-MeSH names to GIZMO metabolite node_ids via the
precomputed name + synonym index.

Results are cached to disk by disease node_id so the slow first-pick
penalty isn't paid twice. Subsequent picks (any user, any process) of
the same disease are O(1).

Coverage gap this fills: diseases that have **no edges** in the GIZMO
graph (degree=0 MONDO terms like generic Crohn / Alzheimer / COVID-19)
still produce a metabolite list, anchored to peer-reviewed literature
rather than the canonical-metabolic-walk that doesn't reach them.
"""
from __future__ import annotations

import json
import os
import re
import time
from pathlib import Path
from typing import Iterable

import pandas as pd


def _cache_root() -> Path:
    """User-scoped cache dir; overridable via ``SQUID_INC_CACHE_DIR``."""
    env = os.environ.get("SQUID_INC_CACHE_DIR")
    if env:
        return Path(env)
    return Path.home() / ".cache" / "squid_inc"


CACHE_DIR = _cache_root() / "pubmed_disease_metab"
NAME_INDEX_PATH = _cache_root() / "metabolite_name_index.json"

_COMPARTMENT_SUFFIX = re.compile(r"\s*\[[^\]]+\]\s*$")
_GENERIC_TERMS = {
    # Drop noise that ChemicalList contains but isn't metabolite biology
    "biomarkers", "antibodies, monoclonal", "antibodies, monoclonal, humanized",
    "cytokines", "anti-bacterial agents", "anti-inflammatory agents",
    "immunosuppressive agents", "gastrointestinal agents", "rna, ribosomal, 16s",
}


def build_metabolite_name_index(mg) -> dict:
    """Build {base_name → node_ids, synonym → node_ids} from a GizmoGraph.

    Strips ``[compartment]`` suffixes; lowercases. Synonyms are pulled
    from the metabolite node's ``synonyms`` attribute when present.
    """
    g = mg.graph if hasattr(mg, "graph") else mg
    name_to_nodes: dict[str, list[str]] = {}
    synonyms_to_nodes: dict[str, list[str]] = {}
    for n, d in g.nodes(data=True):
        if d.get("node_type") != "metabolite":
            continue
        raw_name = d.get("name", "")
        if raw_name:
            base = _COMPARTMENT_SUFFIX.sub("", str(raw_name)).strip().lower()
            if base:
                name_to_nodes.setdefault(base, []).append(n)
        for syn in (d.get("synonyms") or []):
            syn_clean = _COMPARTMENT_SUFFIX.sub("", str(syn)).strip().lower()
            if syn_clean:
                synonyms_to_nodes.setdefault(syn_clean, []).append(n)
    return {"name_to_nodes": name_to_nodes, "synonyms_to_nodes": synonyms_to_nodes}


def _load_name_index(mg=None) -> tuple[dict[str, list[str]], dict[str, list[str]]]:
    """Return (name_to_nodes, synonyms_to_nodes).

    Reads the cached JSON if present. If the cache is missing and a
    GizmoGraph ``mg`` is supplied, builds and caches the index from the
    graph (one-time ~1s pass). Returns empty dicts if neither is available.
    """
    if NAME_INDEX_PATH.exists():
        raw = json.loads(NAME_INDEX_PATH.read_text())
        return raw.get("name_to_nodes", {}), raw.get("synonyms_to_nodes", {})
    if mg is None:
        return {}, {}
    idx = build_metabolite_name_index(mg)
    NAME_INDEX_PATH.parent.mkdir(parents=True, exist_ok=True)
    NAME_INDEX_PATH.write_text(json.dumps(idx))
    return idx["name_to_nodes"], idx["synonyms_to_nodes"]


def _normalize_chem(name: str) -> str:
    return _COMPARTMENT_SUFFIX.sub("", name).strip().lower()


def _resolve_chem_to_nodes(
    chem_name: str,
    name_to_nodes: dict[str, list[str]],
    synonyms_to_nodes: dict[str, list[str]],
) -> list[str]:
    """Return GIZMO metabolite node_ids matching a chemical-MeSH name.

    Match order: exact base-name, exact synonym, exact lowercase. No
    fuzzy/substring match — too many false positives at this scale.
    """
    key = _normalize_chem(chem_name)
    if not key:
        return []
    nodes = list(dict.fromkeys(
        (name_to_nodes.get(key, []) or []) + (synonyms_to_nodes.get(key, []) or [])
    ))
    return nodes


def metabolites_via_literature(
    mg,
    disease_ids: Iterable[str],
    *,
    retmax: int = 1500,
    client=None,
) -> pd.DataFrame:
    """Harvest metabolite candidates for diseases from PubMed co-citation.

    For each disease_id:
      1. Resolve to display name via the graph's disease-node ``name`` attr.
      2. If a cached result file exists, load it.
      3. Else call PubMed: esearch (disease + metabolomics/biomarkers) →
         efetch (chemical MeSH terms) → aggregate.
      4. Map chemical names to GIZMO metabolite node_ids.

    Returns DataFrame with columns:
      metabolite_node_id, score, source_ids, via, chemical_mesh,
      citation_count, total_pmids, disease_query
    """
    from gizmo.sources.pubmed_metabolite import PubMedClient, DiseaseChemicalHits

    g = mg.graph if hasattr(mg, "graph") else mg
    name_to_nodes, synonyms_to_nodes = _load_name_index(mg)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    if client is None:
        client = PubMedClient()

    metab_hits: dict[str, dict] = {}

    for did in disease_ids:
        if did not in g:
            continue
        disease_name = g.nodes[did].get("name") or did
        cache_path = CACHE_DIR / f"{did.replace(':', '_').replace('/', '_')}.json"
        cached: DiseaseChemicalHits | None = None
        if cache_path.exists():
            try:
                blob = json.loads(cache_path.read_text())
                cached = DiseaseChemicalHits(**blob)
            except Exception:
                cached = None

        if cached is None:
            try:
                cached = client.disease_metabolite_hits(disease_name, retmax=retmax)
            except Exception as exc:
                cached = DiseaseChemicalHits(
                    disease_query=disease_name,
                    total_pmids=0,
                    fetched_pmids=[],
                    chemical_counts={"__error__": 0},
                )
                cache_path.write_text(json.dumps({
                    **cached.__dict__,
                    "_error": f"{type(exc).__name__}: {exc}",
                    "_when": time.time(),
                }))
                continue
            cache_path.write_text(json.dumps({
                **cached.__dict__,
                "_when": time.time(),
            }))

        for chem, count in cached.chemical_counts.items():
            if chem in _GENERIC_TERMS or chem.lower() in _GENERIC_TERMS:
                continue
            nodes = _resolve_chem_to_nodes(chem, name_to_nodes, synonyms_to_nodes)
            for nid in nodes:
                row = metab_hits.setdefault(nid, {
                    "metabolite_node_id": nid,
                    "citation_count": 0,
                    "source_diseases": set(),
                    "chemicals": set(),
                })
                row["citation_count"] += int(count)
                row["source_diseases"].add(did)
                row["chemicals"].add(chem)

    rows = []
    for nid, row in metab_hits.items():
        rows.append({
            "metabolite_node_id": nid,
            "score": float(row["citation_count"]),
            "source_ids": ";".join(sorted(row["source_diseases"])),
            "via": "literature·pubmed",
            "chemical_mesh": ";".join(sorted(row["chemicals"])),
            "citation_count": row["citation_count"],
        })
    if not rows:
        return pd.DataFrame(columns=[
            "metabolite_node_id", "score", "source_ids", "via",
            "chemical_mesh", "citation_count",
        ])
    return pd.DataFrame(rows).sort_values("citation_count", ascending=False)
