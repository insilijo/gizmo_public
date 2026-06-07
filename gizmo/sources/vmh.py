"""
Virtual Metabolic Human (VMH) metabolite ID mapping source.

VMH provides a table of metabolite abbreviations with cross-references to
HMDB, KEGG, BiGG, and PubChem.  This module:

1. Fetches (and caches) the VMH metabolite table.
2. Builds HMDB → VMH abbreviation and PubChem → VMH abbreviation mappings.
3. Enriches MetaboliteNode entries in a GizmoGraph with ``vmh_id`` attributes.

License: VMH data is licensed under CC BY 4.0.
Source:  https://www.vmh.life/#downloadview
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Optional

import pandas as pd

log = logging.getLogger(__name__)

_VMH_CSV_URL = "https://www.vmh.life/files/metabolites/all_metabolites.csv"
_VMH_API_URL = "https://www.vmh.life/_api/metabolites/?format=json&page=1&page_size=500"
_CACHE_DEFAULT = Path.home() / ".cache" / "gizmo" / "vmh_metabolites.csv"
# GeMMA caches the same file here — reuse it if present to avoid a redundant download
_GEMMA_CACHE = Path.home() / ".cache" / "gemma" / "vmh_metabolites.csv"


# ---------------------------------------------------------------------------
# Fetch / cache
# ---------------------------------------------------------------------------

def fetch_vmh_metabolites(
    cache_path: Optional[Path] = None,
    force_refresh: bool = False,
) -> pd.DataFrame:
    """Return the VMH metabolite table as a DataFrame.

    Columns: abbreviation, hmdb, kegg, bigg_id, pubchem_compound_id, name

    Results are cached locally.  Pass ``force_refresh=True`` to re-download.
    """
    cache = Path(cache_path) if cache_path else _CACHE_DEFAULT

    # Prefer GeMMA's already-cached copy to avoid a redundant network hit
    for candidate in (cache, _GEMMA_CACHE):
        if candidate.exists() and not force_refresh:
            df = pd.read_csv(candidate, dtype=str, low_memory=False).fillna("")
            log.debug("Loaded VMH metabolites from cache: %s", candidate)
            return df

    # Try CSV download first (single request), then paginated API fallback
    log.info("Fetching VMH metabolite table …")
    import urllib.request

    df: Optional[pd.DataFrame] = None
    try:
        import tempfile, os
        with tempfile.NamedTemporaryFile(delete=False, suffix=".csv") as tmp:
            tmp_path = tmp.name
        urllib.request.urlretrieve(_VMH_CSV_URL, tmp_path)
        df = pd.read_csv(tmp_path, dtype=str, low_memory=False).fillna("")
        os.unlink(tmp_path)
        log.info("Downloaded VMH CSV (%d rows)", len(df))
    except Exception as exc:
        log.warning("VMH CSV download failed (%s); trying paginated API …", exc)

    if df is None:
        rows = []
        next_url: Optional[str] = _VMH_API_URL
        while next_url:
            with urllib.request.urlopen(next_url, timeout=60) as resp:
                payload = json.loads(resp.read().decode())
            for item in payload.get("results", []):
                rows.append({
                    "abbreviation":        item.get("abbreviation", ""),
                    "hmdb":                item.get("hmdb", "") or "",
                    "kegg":                item.get("keggId") or item.get("kegg") or "",
                    "bigg_id":             item.get("biggId") or item.get("bigg_id") or "",
                    "pubchem_compound_id": str(item.get("pubChemId") or item.get("pubchem_compound_id") or ""),
                    "name":                item.get("fullName") or item.get("name") or "",
                })
            next_url = payload.get("next")
        df = pd.DataFrame(rows)
        log.info("Downloaded %d VMH metabolites via API", len(df))

    cache.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(cache, index=False)
    return df


# ---------------------------------------------------------------------------
# Mapping builders
# ---------------------------------------------------------------------------

def build_hmdb_to_vmh(
    vmh_df: Optional[pd.DataFrame] = None,
    *,
    cache_path: Optional[Path] = None,
) -> dict[str, str]:
    """Return ``{hmdb_id_normalised: vmh_abbreviation}`` mapping.

    HMDB IDs are normalised to 7-digit zero-padded form (HMDB0000001) so they
    match the format used by MetaboliteNode and MetaboliteMapper.
    """
    if vmh_df is None:
        vmh_df = fetch_vmh_metabolites(cache_path=cache_path)

    mapping: dict[str, str] = {}
    for _, row in vmh_df.iterrows():
        abbr = str(row.get("abbreviation", "")).strip()
        hmdb = str(row.get("hmdb", "")).strip()
        if not abbr or not hmdb:
            continue
        mapping[_normalise_hmdb(hmdb)] = abbr
    return mapping


def build_pubchem_to_vmh(
    vmh_df: Optional[pd.DataFrame] = None,
    *,
    cache_path: Optional[Path] = None,
) -> dict[str, str]:
    """Return ``{pubchem_cid: vmh_abbreviation}`` mapping."""
    if vmh_df is None:
        vmh_df = fetch_vmh_metabolites(cache_path=cache_path)

    mapping: dict[str, str] = {}
    for _, row in vmh_df.iterrows():
        abbr = str(row.get("abbreviation", "")).strip()
        cid  = str(row.get("pubchem_compound_id", "")).strip()
        if abbr and cid and cid not in ("", "nan", "0"):
            mapping[cid] = abbr
    return mapping


# ---------------------------------------------------------------------------
# Graph enrichment
# ---------------------------------------------------------------------------

def enrich_graph_vmh(mg, *, cache_path: Optional[Path] = None) -> int:
    """Back-fill ``vmh_id`` on MetaboliteNode entries in *mg*.

    Matches nodes by ``hmdb_id`` first, then by ``pubchem_cid``.  Returns the
    number of nodes enriched.

    Parameters
    ----------
    mg:
        A :class:`~gizmo.graph.network.GizmoGraph` instance.
    cache_path:
        Override the default VMH cache location.
    """
    vmh_df = fetch_vmh_metabolites(cache_path=cache_path)
    hmdb_map   = build_hmdb_to_vmh(vmh_df)
    pubchem_map = build_pubchem_to_vmh(vmh_df)

    enriched = 0
    updates: list[tuple[str, str]] = []   # (node_id, vmh_abbr)

    for nid, attrs in mg.graph.nodes(data=True):
        if attrs.get("node_type") != "metabolite":
            continue
        if attrs.get("vmh_id"):
            continue  # already set

        vmh_abbr: Optional[str] = None

        hmdb = attrs.get("hmdb_id", "")
        if hmdb:
            vmh_abbr = hmdb_map.get(_normalise_hmdb(str(hmdb)))

        if not vmh_abbr:
            cid = attrs.get("pubchem_cid", "")
            if cid:
                vmh_abbr = pubchem_map.get(str(cid).split(".")[0])

        if vmh_abbr:
            updates.append((nid, vmh_abbr))

    for nid, vmh_abbr in updates:
        mg.graph.nodes[nid]["vmh_id"] = vmh_abbr
        enriched += 1

    log.info("enrich_graph_vmh: %d/%d metabolite nodes enriched with VMH IDs",
             enriched, sum(1 for _, a in mg.graph.nodes(data=True) if a.get("node_type") == "metabolite"))
    return enriched


# ---------------------------------------------------------------------------
# Convenience: direct HMDB → VMH lookup without a graph
# ---------------------------------------------------------------------------

class VMHMapper:
    """Lightweight HMDB/PubChem → VMH abbreviation mapper (no graph required).

    Useful when you have a list of HMDB IDs (e.g. from Metabolon ``HMDB`` column)
    and just need the corresponding VMH abbreviations for GeMMA pathway enrichment.

    Example
    -------
    >>> mapper = VMHMapper()
    >>> mapper.hmdb_to_vmh("HMDB0000122")   # glucose
    'glc_D'
    """

    def __init__(
        self,
        vmh_df: Optional[pd.DataFrame] = None,
        *,
        cache_path: Optional[Path] = None,
    ) -> None:
        df = vmh_df if vmh_df is not None else fetch_vmh_metabolites(cache_path=cache_path)
        self._hmdb    = build_hmdb_to_vmh(df)
        self._pubchem = build_pubchem_to_vmh(df)

    def hmdb_to_vmh(self, hmdb_id: str) -> Optional[str]:
        """Return VMH abbreviation for *hmdb_id*, or None."""
        return self._hmdb.get(_normalise_hmdb(hmdb_id))

    def pubchem_to_vmh(self, cid: str) -> Optional[str]:
        """Return VMH abbreviation for PubChem CID, or None."""
        return self._pubchem.get(str(cid).split(".")[0])

    def map_hmdb_list(self, hmdb_ids: list[str]) -> dict[str, Optional[str]]:
        """Map a list of HMDB IDs → VMH abbreviations (None for unmapped)."""
        return {h: self.hmdb_to_vmh(h) for h in hmdb_ids}

    def coverage(self, hmdb_ids: list[str]) -> float:
        """Fraction of *hmdb_ids* that resolve to a VMH abbreviation."""
        if not hmdb_ids:
            return 0.0
        mapped = sum(1 for h in hmdb_ids if self.hmdb_to_vmh(h) is not None)
        return mapped / len(hmdb_ids)


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _normalise_hmdb(s: str) -> str:
    """Normalise HMDB IDs to 7-digit zero-padded form (HMDB0000001)."""
    s = s.strip().upper()
    if s.startswith("HMDB"):
        return f"HMDB{s[4:].zfill(7)}"
    return s
