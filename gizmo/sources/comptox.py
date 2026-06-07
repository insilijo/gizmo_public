"""
EPA CompTox / DSSTox chemical data (CC BY 4.0).

https://comptox.epa.gov/dashboard
https://api.epa.gov/CCTE/

CompTox provides:
  - DTXSID identifiers (EPA's canonical chemical ID)
  - CAS Registry Numbers
  - Hazard flags: carcinogenicity, reproductive, developmental, neurotoxicity, etc.
  - Physicochemical properties
  - Bioactivity summaries

This module enriches existing MetaboliteNodes in the graph with:
  - cas_id        (CAS Registry Number)
  - dtxsid        (stored as a node attribute, not in schema — added dynamically)
  - EPA hazard flags (stored as node attribute "epa_hazard_flags": list[str])

Matching: InChIKey → CompTox API → DTXSID + CAS.
Bulk CSV mode is also supported for offline/batch workflows.

Usage::

    from gizmo.sources.comptox import CompToxClient

    client = CompToxClient(cache_dir="data/raw/comptox")
    n = client.enrich_graph(mg)
    print(f"Enriched {n} nodes")
"""

from __future__ import annotations

import json
import logging
import time
from pathlib import Path
from typing import Optional

import requests

log = logging.getLogger(__name__)

_API_BASE = "https://api.epa.gov/CCTE"
_BY_INCHIKEY = _API_BASE + "/chemical/detail/search/by-inchikey/{inchikey}"
_RATE_LIMIT_S = 0.25   # conservative: ~240 req/min


class CompToxClient:
    """
    Enriches MetaboliteNodes with EPA CompTox/DSSTox data via the CCTE REST API.

    Parameters
    ----------
    cache_dir  : directory for caching API responses (JSON per InChIKey)
    rate_limit : seconds between API calls
    """

    def __init__(
        self,
        cache_dir: str | Path = "data/raw/comptox",
        rate_limit: float = _RATE_LIMIT_S,
    ) -> None:
        self.cache_dir  = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.rate_limit = rate_limit
        self._session   = requests.Session()
        self._session.headers.update({
            "Accept": "application/json",
            "User-Agent": "gizmo/1.0",
        })

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def enrich_graph(self, mg, force: bool = False) -> int:
        """
        Query CompTox for every metabolite node that has an InChIKey but no
        cas_id or dtxsid, and annotate those nodes in-place.

        Returns the number of nodes enriched.
        """
        enriched = 0
        for nid, attrs in mg.graph.nodes(data=True):
            if attrs.get("node_type") != "metabolite":
                continue
            inchikey = (attrs.get("inchikey") or "").strip()
            if not inchikey:
                continue
            if not force and attrs.get("cas_id"):
                continue   # already have CAS, skip

            data = self._lookup_inchikey(inchikey)
            if not data:
                continue

            changed = False
            cas = data.get("casrn") or data.get("casNumber") or ""
            if cas and not attrs.get("cas_id"):
                mg.graph.nodes[nid]["cas_id"] = cas.strip()
                changed = True

            dtxsid = data.get("dtxsid") or data.get("dsstoxSubstanceId") or ""
            if dtxsid and not attrs.get("dtxsid"):
                mg.graph.nodes[nid]["dtxsid"] = dtxsid.strip()
                changed = True

            hazard_flags = self._extract_hazard_flags(data)
            if hazard_flags:
                existing = set(attrs.get("epa_hazard_flags") or [])
                merged = sorted(existing | set(hazard_flags))
                mg.graph.nodes[nid]["epa_hazard_flags"] = merged
                changed = True

            if changed:
                enriched += 1

        log.info("CompTox: enriched %d metabolite nodes.", enriched)
        return enriched

    def lookup_inchikey(self, inchikey: str) -> Optional[dict]:
        """Public wrapper — look up a single InChIKey, returns raw API dict or None."""
        return self._lookup_inchikey(inchikey)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _lookup_inchikey(self, inchikey: str) -> Optional[dict]:
        cache_path = self.cache_dir / f"{inchikey}.json"
        if cache_path.exists():
            try:
                return json.loads(cache_path.read_text())
            except Exception:
                pass

        url = _BY_INCHIKEY.format(inchikey=inchikey)
        try:
            time.sleep(self.rate_limit)
            resp = self._session.get(url, timeout=20)
            if resp.status_code == 404:
                cache_path.write_text("null")
                return None
            resp.raise_for_status()
            data = resp.json()
            # API may return a list; take first element
            if isinstance(data, list):
                data = data[0] if data else None
            cache_path.write_text(json.dumps(data, indent=2))
            return data
        except Exception as exc:
            log.debug("CompTox lookup failed for %s: %s", inchikey, exc)
            return None

    @staticmethod
    def _extract_hazard_flags(data: dict) -> list[str]:
        """
        Extract EPA hazard category flags from a CompTox API response.

        The CCTE API nests hazard data under various keys depending on the
        endpoint version.  We inspect known fields and return a list of
        short flag strings (e.g. "carcinogen", "reproductive_toxicant").
        """
        flags: list[str] = []
        if not data:
            return flags

        # hazardSummary dict (older API)
        hazard = data.get("hazardSummary") or {}
        _flag_map = {
            "carcinogenicity":        "carcinogen",
            "reproductiveToxicity":   "reproductive_toxicant",
            "developmentalToxicity":  "developmental_toxicant",
            "neurotoxicity":          "neurotoxicant",
            "acuteMammalianToxicity": "acute_toxicant",
            "skinEyeIrritation":      "irritant",
            "mutagenicity":           "mutagen",
        }
        for api_key, flag in _flag_map.items():
            val = hazard.get(api_key)
            if val and str(val).lower() not in ("", "none", "null", "false", "0"):
                flags.append(flag)

        # hazardData list (newer API)
        for entry in data.get("hazardData") or []:
            endpoint = (entry.get("hazardEndpoint") or "").lower()
            value    = (entry.get("hazardValue")    or "").lower()
            if "cancer" in endpoint and value not in ("", "none", "not classified"):
                flags.append("carcinogen")
            elif "repro" in endpoint and value not in ("", "none", "not classified"):
                flags.append("reproductive_toxicant")
            elif "devel" in endpoint and value not in ("", "none", "not classified"):
                flags.append("developmental_toxicant")
            elif "neuro" in endpoint and value not in ("", "none", "not classified"):
                flags.append("neurotoxicant")

        return sorted(set(flags))
