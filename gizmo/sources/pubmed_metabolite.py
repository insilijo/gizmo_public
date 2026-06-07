"""NCBI PubMed E-utilities client for disease ↔ chemical-MeSH co-citation.

Free public API (no auth required). Honors ``NCBI_API_KEY`` env var if
set (raises throughput from 3 to 10 req/sec). User-Agent identifies
Insilijo Science as the caller per NCBI policy.
"""
from __future__ import annotations

import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter
from dataclasses import dataclass

_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

_USER_AGENT = (
    "SQuID-INC-PanelDesigner/1.0 (joseph.j.gardner@gmail.com; "
    "https://insilijo.github.io)"
)
_TIMEOUT = 30


@dataclass
class DiseaseChemicalHits:
    disease_query: str
    total_pmids: int
    fetched_pmids: list[str]
    chemical_counts: dict[str, int]


class PubMedClient:
    """Thin client over E-utilities. Polite rate limiting + retries."""

    def __init__(self, api_key: str | None = None, *, min_interval: float | None = None):
        self.api_key = api_key or os.environ.get("NCBI_API_KEY") or None
        # 10 req/s with key, 3 req/s without → 0.105 / 0.34 s between calls
        if min_interval is None:
            min_interval = 0.11 if self.api_key else 0.34
        self.min_interval = min_interval
        self._last = 0.0

    def _wait(self) -> None:
        delta = time.time() - self._last
        if delta < self.min_interval:
            time.sleep(self.min_interval - delta)
        self._last = time.time()

    def _get(self, url: str, params: dict, *, retries: int = 3) -> bytes:
        params = dict(params)
        if self.api_key:
            params["api_key"] = self.api_key
        req = urllib.request.Request(
            url + "?" + urllib.parse.urlencode(params),
            headers={"User-Agent": _USER_AGENT, "Accept": "application/json"},
        )
        last_err: Exception | None = None
        for attempt in range(retries):
            self._wait()
            try:
                with urllib.request.urlopen(req, timeout=_TIMEOUT) as r:
                    return r.read()
            except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as exc:
                last_err = exc
                time.sleep(0.6 * (2 ** attempt))
        assert last_err is not None
        raise last_err

    def search_disease_metabolomics(
        self, disease_term: str, *, retmax: int = 300
    ) -> tuple[int, list[str]]:
        """Find PubMed papers that link a disease to metabolomics / biomarkers /
        small-molecule biology. Returns (total_count, fetched_pmids).

        Query design tries the disease as a MeSH heading first (cleanest
        match) and falls back to title/abstract free-text inside the OR:
            ("<disease>"[MeSH] OR "<disease>"[Title/Abstract])
            AND ("metabolomics"[MeSH] OR "biomarkers"[MeSH]
                 OR "metabolite"[Title/Abstract])
        """
        term = (
            f'("{disease_term}"[MeSH] OR "{disease_term}"[Title/Abstract])'
            ' AND ("metabolomics"[MeSH] OR "biomarkers"[MeSH]'
            ' OR "metabolite"[Title/Abstract])'
        )
        params = {"db": "pubmed", "term": term, "retmax": retmax, "retmode": "json"}
        raw = self._get(_ESEARCH, params)
        body = json.loads(raw)
        result = body.get("esearchresult", {})
        total = int(result.get("count", "0") or 0)
        pmids = result.get("idlist", []) or []
        return total, pmids

    def fetch_chemicals(self, pmids: list[str]) -> dict[str, int]:
        """For a batch of PMIDs, parse ChemicalList entries and return a
        ``{chemical_name: article_count}`` aggregation. Empty if no PMIDs.

        Batches PMIDs to 200 per efetch call (NCBI limit).
        """
        if not pmids:
            return {}
        counts: Counter[str] = Counter()
        for i in range(0, len(pmids), 200):
            batch = pmids[i:i + 200]
            params = {
                "db": "pubmed",
                "id": ",".join(batch),
                "rettype": "xml",
                "retmode": "xml",
            }
            xml_bytes = self._get(_EFETCH, params)
            root = ET.fromstring(xml_bytes)
            for article in root.iter("PubmedArticle"):
                # Each article contributes 1 to each distinct ChemicalList entry.
                seen_in_article: set[str] = set()
                for chem in article.iter("Chemical"):
                    name_el = chem.find("NameOfSubstance")
                    if name_el is None or not name_el.text:
                        continue
                    nm = name_el.text.strip()
                    if nm and nm not in seen_in_article:
                        seen_in_article.add(nm)
                        counts[nm] += 1
        return dict(counts)

    def disease_metabolite_hits(
        self, disease_term: str, *, retmax: int = 300
    ) -> DiseaseChemicalHits:
        total, pmids = self.search_disease_metabolomics(disease_term, retmax=retmax)
        chem_counts = self.fetch_chemicals(pmids)
        return DiseaseChemicalHits(
            disease_query=disease_term,
            total_pmids=total,
            fetched_pmids=pmids,
            chemical_counts=chem_counts,
        )
