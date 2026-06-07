"""
Metabolon data dictionary loader and ChEBI mapper.

Input: Metabolon global data dictionary CSV (PMC Open Access subset, 2024-04-14)
Columns: BIOCHEMICAL, PATHWAY, ChemRICHClass, PUBCHEM, INCHIKEY, PLATFORM, MASS, RI

ChEBI mapping strategy (in order of preference):
  1. InChIKey → MetaNetX chem_prop.tsv local lookup (fast, no rate limiting)
  2. InChIKey → ChEBI OLS4 REST API
  3. PubChem CID → PubChem PUG REST → InChIKey → step 1/2
  4. Unmatched: node_id = "metabolon:{BIOCHEMICAL}" (name-based fallback)

Match confidence levels:
  "exact_inchikey"   — InChIKey matched in MetaNetX / OLS4
  "pubchem_inchikey" — resolved via PubChem CID → InChIKey → ChEBI
  "unmatched"        — no ChEBI found

Usage::

    from gizmo.sources.metabolon import MetabolonLoader
    loader = MetabolonLoader("data/resources/gizmo/sources/metabolon_data_dictionary_PMC_OA_subset_4.14.2024.csv")
    # Fast path: use shared MetaNetX resources for bulk InChIKey matching
    loader.load_metanetx_index(
        "data/resources/gizmo/metanetx/chem_prop.tsv",
        "data/resources/gizmo/metanetx/chem_xref.tsv",
    )
    nodes, report = loader.to_metabolite_nodes()
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pandas as pd
import requests

from gizmo.schema import MetaboliteNode
from gizmo.sources.metanetx import _mnx_header_info

log = logging.getLogger(__name__)

_OLS4_SEARCH    = "https://www.ebi.ac.uk/ols4/api/search"
_PUBCHEM_REST   = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/InChIKey/JSON"
_PUBCHEM_NAME   = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/InChIKey,CID/JSON"


@dataclass
class MatchReport:
    total: int = 0
    exact_inchikey: int = 0        # InChIKey exact match via MetaNetX (local)
    connectivity_inchikey: int = 0 # InChIKey connectivity-layer match via MetaNetX (local)
    pubchem_mnx: int = 0           # PubChem CID match via MetaNetX chem_xref (local)
    api_inchikey: int = 0          # InChIKey match via ChEBI OLS4 API
    pubchem_inchikey: int = 0      # PubChem CID → InChIKey → ChEBI via API
    name_pubchem: int = 0          # biochemical name → PubChem → InChIKey → ChEBI
    unmatched: int = 0

    hmdb_matched: int = 0

    @property
    def chebi_coverage(self) -> float:
        if self.total == 0:
            return 0.0
        matched = (self.exact_inchikey + self.connectivity_inchikey
                   + self.pubchem_mnx + self.api_inchikey + self.pubchem_inchikey
                   + self.name_pubchem)
        return matched / self.total

    @property
    def hmdb_coverage(self) -> float:
        return self.hmdb_matched / self.total if self.total else 0.0

    def __str__(self) -> str:
        return (
            f"Metabolon mapping: {self.total} compounds | "
            f"InChIKey(exact): {self.exact_inchikey} | "
            f"InChIKey(connectivity): {self.connectivity_inchikey} | "
            f"PubChem(local): {self.pubchem_mnx} | "
            f"InChIKey(API): {self.api_inchikey} | "
            f"PubChem(API): {self.pubchem_inchikey} | "
            f"Name→PubChem: {self.name_pubchem} | "
            f"Unmatched: {self.unmatched} | "
            f"ChEBI coverage: {self.chebi_coverage*100:.1f}% | "
            f"HMDB coverage: {self.hmdb_coverage*100:.1f}%"
        )


class MetabolonLoader:
    """
    Load Metabolon data dictionary and map compounds to MetaboliteNode objects.

    The core matching is done via InChIKey against the MetaNetX chem_prop table
    (which is cached locally). This avoids any API rate-limiting for bulk loads.
    """

    def __init__(self, csv_path: str | Path) -> None:
        self.csv_path = Path(csv_path)
        self._df: Optional[pd.DataFrame] = None
        # InChIKey → ChEBI ID (from MetaNetX chem_prop + chem_xref, exact match)
        self._inchikey_to_chebi: dict[str, str] = {}
        # InChIKey connectivity prefix (14 chars) → ChEBI ID (fallback for stereoisomers)
        self._connectivity_to_chebi: dict[str, str] = {}
        # PubChem CID string → ChEBI ID (from MetaNetX chem_xref pubchem: entries)
        self._pubchem_to_chebi: dict[str, str] = {}
        # MNX compound ID → ChEBI ID (from chem_xref)
        self._mnx_to_chebi: dict[str, str] = {}
        # InChIKey → HMDB ID (7-digit VMH-compatible, from MetaNetX chem_xref hmdb: entries)
        self._inchikey_to_hmdb: dict[str, str] = {}
        # InChIKey connectivity prefix → HMDB ID
        self._connectivity_to_hmdb: dict[str, str] = {}
        # PubChem CID → HMDB ID
        self._pubchem_to_hmdb: dict[str, str] = {}
        # InChIKey → MNX compound ID
        self._inchikey_to_mnx: dict[str, str] = {}
        # InChIKey → SMILES (from MetaNetX chem_prop)
        self._inchikey_to_smiles: dict[str, str] = {}
        # InChIKey → InChI
        self._inchikey_to_inchi: dict[str, str] = {}
        # InChIKey → molecular formula
        self._inchikey_to_formula: dict[str, str] = {}
        # Name → (InChIKey, PubChem CID) from PubChem name search API
        # Populated by run_auto_matching() or to_metabolite_nodes(name_fallback=True)
        self._name_to_inchikey: dict[str, str] = {}   # biochemical name → InChIKey
        self._name_to_cid:      dict[str, str] = {}   # biochemical name → PubChem CID
        # Manual overrides populated by MetabolonCurator.apply().
        self._name_overrides: dict[str, str] = {}
        self._name_hmdb_overrides: dict[str, str] = {}

    # ------------------------------------------------------------------
    # Index loading
    # ------------------------------------------------------------------

    def load_metanetx_index(
        self,
        chem_prop_path: str | Path,
        chem_xref_path: str | Path,
    ) -> int:
        """
        Build InChIKey → ChEBI lookup from MetaNetX flat files.

        Streams both files line-by-line to avoid loading 700MB+ DataFrames.
        Returns number of indexed InChIKey → ChEBI entries.

        chem_prop.tsv columns (v4.x): ID, name, reference, formula, charge, mass, InChI, InChIKey, SMILES
        chem_xref.tsv columns: source, ID, description
          (source looks like "chebi:CHEBI:15422" or "chebi:15422")
        """
        import re

        chem_prop_path = Path(chem_prop_path)
        chem_xref_path = Path(chem_xref_path)

        if not chem_prop_path.exists() or not chem_xref_path.exists():
            log.warning(
                "MetaNetX files not found (%s, %s). Run MetaNetXClient().download() first.",
                chem_prop_path, chem_xref_path,
            )
            return 0

        # --- Pass 1: stream chem_xref ---
        #   chebi: rows  → {MNX_ID: ChEBI_ID}
        #   pubchem: rows → {pubchem_cid: MNX_ID}  (joined with chebi map after pass)
        xref_cols, xref_start = _mnx_header_info(chem_xref_path)
        # Typical columns: source (0), ID (1), description (2)
        mnx_idx = xref_cols.index("ID") if "ID" in xref_cols else 1

        mnx_to_chebi: dict[str, str] = {}
        mnx_to_hmdb: dict[str, str] = {}         # MNX_ID → HMDB ID (first seen, normalised)
        mnx_from_pubchem: dict[str, str] = {}   # pubchem_cid → mnx_id (temporary)
        with open(chem_xref_path, encoding="utf-8") as fh:
            for lineno, line in enumerate(fh):
                if lineno < xref_start:
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= mnx_idx:
                    continue
                source = parts[0]
                mnx_id = parts[mnx_idx]
                if source.startswith("chebi:"):
                    # "chebi:CHEBI:15422" → "CHEBI:15422";  "chebi:15422" → "CHEBI:15422"
                    chebi_id = "CHEBI:" + re.sub(r"^chebi:(?:CHEBI:)?", "", source)
                    mnx_to_chebi[mnx_id] = chebi_id
                elif source.startswith("hmdb:"):
                    if mnx_id not in mnx_to_hmdb:   # keep first (canonical) entry
                        mnx_to_hmdb[mnx_id] = _normalise_hmdb(source[5:])
                elif source.startswith("pubchem:"):
                    cid = source[8:]   # strip "pubchem:"
                    mnx_from_pubchem[cid] = mnx_id

        if not mnx_to_chebi:
            log.warning("No ChEBI entries found in chem_xref.tsv")
            return 0

        self._mnx_to_chebi = mnx_to_chebi
        # Join pubchem → mnx with mnx → chebi to get pubchem → chebi
        self._pubchem_to_chebi = {
            cid: mnx_to_chebi[mid]
            for cid, mid in mnx_from_pubchem.items()
            if mid in mnx_to_chebi
        }
        # Join pubchem → mnx with mnx → hmdb to get pubchem → hmdb
        self._pubchem_to_hmdb = {
            cid: mnx_to_hmdb[mid]
            for cid, mid in mnx_from_pubchem.items()
            if mid in mnx_to_hmdb
        }
        log.info(
            "chem_xref: %d MNX→ChEBI, %d MNX→HMDB, %d PubChem→ChEBI entries",
            len(mnx_to_chebi), len(mnx_to_hmdb), len(self._pubchem_to_chebi),
        )

        # --- Pass 2: stream chem_prop → build {InChIKey: ChEBI_ID/SMILES/InChI/formula/MNX} ---
        prop_cols, prop_start = _mnx_header_info(chem_prop_path)

        def _col(name: str, fallback: int = -1) -> int:
            try:
                return next(i for i, c in enumerate(prop_cols) if c.lower() == name.lower())
            except StopIteration:
                return fallback

        id_idx = _col("ID", 0)
        ik_idx = _col("InChIKey")
        smiles_idx = _col("SMILES")
        inchi_idx = _col("InChI")
        formula_idx = _col("formula")

        if ik_idx < 0:
            log.warning("No InChIKey column found in chem_prop.tsv")
            return 0

        self._inchikey_to_chebi = {}
        self._connectivity_to_chebi = {}
        self._inchikey_to_hmdb = {}
        self._connectivity_to_hmdb = {}
        self._inchikey_to_mnx = {}
        self._inchikey_to_smiles = {}
        self._inchikey_to_inchi = {}
        self._inchikey_to_formula = {}
        count = 0
        with open(chem_prop_path, encoding="utf-8") as fh:
            for lineno, line in enumerate(fh):
                if lineno < prop_start:
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= max(id_idx, ik_idx):
                    continue
                mnx_id = parts[id_idx]
                inchikey = parts[ik_idx].strip()
                if not inchikey:
                    continue
                pfx = inchikey[:14]

                if mnx_id in mnx_to_chebi:
                    chebi = mnx_to_chebi[mnx_id]
                    self._inchikey_to_chebi[inchikey] = chebi
                    self._inchikey_to_mnx.setdefault(inchikey, mnx_id)
                    # Store structural data only for ChEBI-mapped compounds (keeps memory bounded)
                    if smiles_idx >= 0 and smiles_idx < len(parts):
                        smiles = parts[smiles_idx].strip()
                        if smiles:
                            self._inchikey_to_smiles.setdefault(inchikey, smiles)
                    if inchi_idx >= 0 and inchi_idx < len(parts):
                        inchi = parts[inchi_idx].strip()
                        if inchi:
                            self._inchikey_to_inchi.setdefault(inchikey, inchi)
                    if formula_idx >= 0 and formula_idx < len(parts):
                        formula = parts[formula_idx].strip()
                        if formula:
                            self._inchikey_to_formula.setdefault(inchikey, formula)
                    # Connectivity prefix index (first 14 chars = connectivity layer).
                    # Only store the first entry seen — gives a canonical ChEBI for
                    # any stereoisomer / charge variant with the same skeleton.
                    if pfx not in self._connectivity_to_chebi:
                        self._connectivity_to_chebi[pfx] = chebi
                    count += 1
                if mnx_id in mnx_to_hmdb:
                    hmdb = mnx_to_hmdb[mnx_id]
                    self._inchikey_to_hmdb[inchikey] = hmdb
                    if pfx not in self._connectivity_to_hmdb:
                        self._connectivity_to_hmdb[pfx] = hmdb

        log.info(
            "MetaNetX index: %d exact InChIKey→ChEBI, %d InChIKey→HMDB, %d connectivity entries",
            count, len(self._inchikey_to_hmdb), len(self._connectivity_to_chebi),
        )
        return count

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    def _load_df(self) -> pd.DataFrame:
        if self._df is None:
            self._df = pd.read_csv(self.csv_path, encoding="utf-8-sig")
        return self._df

    def raw_df(self) -> pd.DataFrame:
        """Return the raw Metabolon CSV as a DataFrame."""
        return self._load_df().copy()

    # ------------------------------------------------------------------
    # Automated matching
    # ------------------------------------------------------------------

    def run_auto_matching(
        self,
        cache_path: Optional[str | Path] = None,
        rate_limit_s: float = 0.15,
        max_api_calls: int = 5000,
    ) -> int:
        """
        Run all automated matching tiers on still-unmatched compounds.

        Tier 1 (local) — already done by load_metanetx_index().
        Tier 2 (API) — InChIKey → ChEBI OLS4 + PubChem CID → InChIKey → ChEBI.
        Tier 3 (name) — biochemical name → PubChem name search → InChIKey → ChEBI.

        Results are cached to ``cache_path`` (JSON) so subsequent runs skip
        already-resolved names.  Pass the same path in subsequent sessions to
        accumulate coverage incrementally.

        Parameters
        ----------
        cache_path     : JSON file to persist name→InChIKey results across runs
        rate_limit_s   : seconds between API calls (default 0.15 ≈ 400 req/min)
        max_api_calls  : hard cap on API calls per run (safety valve)

        Returns number of new name→InChIKey mappings resolved this run.
        """
        import json as _json
        import time

        # Load persisted cache
        cache: dict[str, dict] = {}
        if cache_path:
            p = Path(cache_path)
            if p.exists():
                try:
                    cache = _json.loads(p.read_text())
                except Exception:
                    pass

        df = self._load_df()
        n_calls = 0
        n_new   = 0

        for _, row in df.iterrows():
            if n_calls >= max_api_calls:
                log.warning(
                    "run_auto_matching: hit max_api_calls=%d cap. "
                    "Run again to continue.", max_api_calls,
                )
                break

            biochemical = str(row.get("BIOCHEMICAL", "")).strip()
            if not biochemical:
                continue

            inchikey    = _clean(row.get("INCHIKEY"))
            pubchem_cid = _clean(row.get("PUBCHEM"))

            # Skip if already matched locally
            already_matched = (
                (inchikey and inchikey in self._inchikey_to_chebi)
                or (inchikey and len(inchikey) >= 14 and inchikey[:14] in self._connectivity_to_chebi)
                or (pubchem_cid and pubchem_cid in self._pubchem_to_chebi)
                or (biochemical in self._name_overrides)
            )
            if already_matched:
                continue

            # Already resolved in cache from a previous run
            if biochemical in cache:
                entry = cache[biochemical]
                ik  = entry.get("inchikey")
                cid = entry.get("cid")
                if ik:
                    self._name_to_inchikey[biochemical] = ik
                if cid:
                    self._name_to_cid[biochemical] = cid
                continue

            # Tier 2: OLS4 InChIKey lookup (if we have an InChIKey but it wasn't in MetaNetX)
            if inchikey and inchikey not in self._inchikey_to_chebi:
                chebi = _chebi_from_inchikey_api(inchikey)
                time.sleep(rate_limit_s)
                n_calls += 1
                if chebi:
                    # Backfill the local index so to_metabolite_nodes() picks it up
                    self._inchikey_to_chebi[inchikey] = chebi
                    if len(inchikey) >= 14:
                        self._connectivity_to_chebi.setdefault(inchikey[:14], chebi)
                    cache[biochemical] = {"inchikey": inchikey, "chebi": chebi, "tier": 2}
                    n_new += 1
                    continue

            # Tier 2b: PubChem CID → InChIKey → OLS4
            if pubchem_cid and pubchem_cid not in self._pubchem_to_chebi:
                resolved_ik = _inchikey_from_pubchem(pubchem_cid)
                time.sleep(rate_limit_s)
                n_calls += 1
                if resolved_ik:
                    chebi = _chebi_from_inchikey_api(resolved_ik)
                    time.sleep(rate_limit_s)
                    n_calls += 1
                    if chebi:
                        self._pubchem_to_chebi[pubchem_cid] = chebi
                        self._inchikey_to_chebi[resolved_ik]  = chebi
                        if len(resolved_ik) >= 14:
                            self._connectivity_to_chebi.setdefault(resolved_ik[:14], chebi)
                        cache[biochemical] = {
                            "inchikey": resolved_ik, "cid": pubchem_cid,
                            "chebi": chebi, "tier": "2b",
                        }
                        n_new += 1
                        continue

            # Tier 3: name → PubChem
            result = _inchikey_from_name_pubchem(biochemical)
            time.sleep(rate_limit_s)
            n_calls += 1
            if result:
                ik, cid = result
                self._name_to_inchikey[biochemical] = ik
                if cid:
                    self._name_to_cid[biochemical] = cid

                chebi = (
                    self._inchikey_to_chebi.get(ik)
                    or (self._connectivity_to_chebi.get(ik[:14]) if len(ik) >= 14 else None)
                    or _chebi_from_inchikey_api(ik)
                )
                if chebi:
                    time.sleep(rate_limit_s)
                    n_calls += 1
                    self._inchikey_to_chebi[ik] = chebi
                    if len(ik) >= 14:
                        self._connectivity_to_chebi.setdefault(ik[:14], chebi)
                    if cid:
                        self._pubchem_to_chebi[cid] = chebi
                    cache[biochemical] = {
                        "inchikey": ik, "cid": cid,
                        "chebi": chebi, "tier": 3,
                    }
                    n_new += 1
                else:
                    cache[biochemical] = {"inchikey": ik, "cid": cid, "tier": 3}
            else:
                # Mark as tried-and-failed so we don't re-query next run
                cache[biochemical] = {"tier": 3, "result": None}

        # Persist cache
        if cache_path:
            p = Path(cache_path)
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(_json.dumps(cache, indent=2))

        log.info(
            "run_auto_matching: %d API calls, %d new resolutions. "
            "Total name→IK entries: %d",
            n_calls, n_new, len(self._name_to_inchikey),
        )
        return n_new

    # ------------------------------------------------------------------
    # Node generation
    # ------------------------------------------------------------------

    def to_metabolite_nodes(
        self,
        *,
        api_fallback: bool = False,
        name_fallback: bool = True,
        api_rate_limit_s: float = 0.2,
    ) -> tuple[list[MetaboliteNode], MatchReport]:
        """
        Convert Metabolon CSV rows to MetaboliteNode objects with ChEBI IDs where possible.

        Recommended workflow for maximum automated coverage:
            1. loader.load_metanetx_index(...)         # local, fast
            2. loader.run_auto_matching(cache_path=...) # API tiers 2+3, one-time
            3. nodes, report = loader.to_metabolite_nodes()
            4. curator.widget()                        # manual for remaining ~10-15%

        Parameters
        ----------
        api_fallback : bool
            If True, call ChEBI OLS4 / PubChem REST inline for rows not yet matched.
            Prefer run_auto_matching() instead — it caches results and is resumable.
        name_fallback : bool
            If True, use any name→InChIKey results cached by run_auto_matching().
            Has no effect if run_auto_matching() has not been called.
        api_rate_limit_s : float
            Seconds to sleep between API calls (only when api_fallback=True).

        Returns
        -------
        (nodes, report)
        """
        import time

        df = self._load_df()
        report = MatchReport(total=len(df))
        nodes: list[MetaboliteNode] = []

        for _, row in df.iterrows():
            biochemical = str(row.get("BIOCHEMICAL", "")).strip()
            inchikey = _clean(row.get("INCHIKEY"))
            pubchem_cid = _clean(row.get("PUBCHEM"))
            pathway = _clean(row.get("PATHWAY"))
            platform = _clean(row.get("PLATFORM"))
            mass_raw = row.get("MASS")
            mass = float(mass_raw) if _is_numeric(mass_raw) else None
            ri_raw = row.get("RI")
            ri = float(ri_raw) if _is_numeric(ri_raw) else None

            chebi_id: Optional[str] = None
            hmdb_id: Optional[str] = None
            matched_inchikey: Optional[str] = inchikey
            confidence = "unmatched"
            manually_reviewed = False

            # --- Match 1: MetaNetX InChIKey exact (local) ---
            if inchikey and inchikey in self._inchikey_to_chebi:
                chebi_id = self._inchikey_to_chebi[inchikey]
                hmdb_id = self._inchikey_to_hmdb.get(inchikey)
                confidence = "exact_inchikey"

            # --- Match 1b: MetaNetX InChIKey connectivity layer (local) ---
            # Catches stereoisomers / charge variants with the same carbon skeleton.
            elif inchikey and len(inchikey) >= 14 and inchikey[:14] in self._connectivity_to_chebi:
                pfx = inchikey[:14]
                chebi_id = self._connectivity_to_chebi[pfx]
                hmdb_id = self._connectivity_to_hmdb.get(pfx)
                confidence = "connectivity_inchikey"

            # --- Match 2: MetaNetX PubChem CID index (local) ---
            elif pubchem_cid and pubchem_cid in self._pubchem_to_chebi:
                chebi_id = self._pubchem_to_chebi[pubchem_cid]
                hmdb_id = self._pubchem_to_hmdb.get(pubchem_cid)
                confidence = "pubchem_mnx"

            # --- Match 3: Manual override (from MetabolonCurator) ---
            elif biochemical in self._name_overrides or biochemical in self._name_hmdb_overrides:
                chebi_id = self._name_overrides.get(biochemical)
                hmdb_id = self._name_hmdb_overrides.get(biochemical)
                confidence = "exact_inchikey"  # treat as exact for reporting
                manually_reviewed = True

            # --- Match 3b: Name → PubChem cache (from run_auto_matching) ---
            elif name_fallback and biochemical in self._name_to_inchikey:
                resolved_ik = self._name_to_inchikey[biochemical]
                matched_inchikey = resolved_ik
                chebi_id = (
                    self._inchikey_to_chebi.get(resolved_ik)
                    or (self._connectivity_to_chebi.get(resolved_ik[:14])
                        if len(resolved_ik) >= 14 else None)
                )
                hmdb_id = (
                    self._inchikey_to_hmdb.get(resolved_ik)
                    or (self._connectivity_to_hmdb.get(resolved_ik[:14])
                        if len(resolved_ik) >= 14 else None)
                )
                if not pubchem_cid and biochemical in self._name_to_cid:
                    pubchem_cid = self._name_to_cid[biochemical]
                if chebi_id:
                    confidence = "name_pubchem"

            # --- Match 4: API fallback ---
            elif api_fallback:
                if inchikey:
                    chebi_id = _chebi_from_inchikey_api(inchikey)
                    if chebi_id:
                        confidence = "api_inchikey"

                if not chebi_id and pubchem_cid:
                    try:
                        resolved_ik = _inchikey_from_pubchem(pubchem_cid)
                        if resolved_ik:
                            matched_inchikey = resolved_ik
                            chebi_id = _chebi_from_inchikey_api(resolved_ik)
                            if chebi_id:
                                confidence = "pubchem_inchikey"
                    except Exception as exc:
                        log.debug("PubChem lookup failed for %s: %s", pubchem_cid, exc)

                if api_rate_limit_s > 0:
                    time.sleep(api_rate_limit_s)

            # --- Build node_id ---
            if chebi_id:
                node_id = chebi_id
            elif pubchem_cid:
                node_id = f"pubchem:{pubchem_cid}"
            else:
                # Sanitise biochemical name for use as ID
                safe_name = biochemical.replace(" ", "_").replace("/", "-")[:80]
                node_id = f"metabolon:{safe_name}"

            ik = matched_inchikey or ""
            node = MetaboliteNode(
                node_id=node_id,
                chebi_id=chebi_id,
                hmdb_id=hmdb_id,
                metanetx_id=self._inchikey_to_mnx.get(ik) if ik else None,
                pubchem_cid=pubchem_cid,
                metabolon_name=biochemical,
                name=biochemical,
                formula=self._inchikey_to_formula.get(ik) if ik else None,
                inchi=self._inchikey_to_inchi.get(ik) if ik else None,
                inchikey=matched_inchikey,
                smiles=self._inchikey_to_smiles.get(ik) if ik else None,
                mass=mass,
                platform=platform,
                retention_index=ri,
                manually_reviewed=manually_reviewed,
            )
            nodes.append(node)

            if confidence == "exact_inchikey":
                report.exact_inchikey += 1
            elif confidence == "connectivity_inchikey":
                report.connectivity_inchikey += 1
            elif confidence == "pubchem_mnx":
                report.pubchem_mnx += 1
            elif confidence == "api_inchikey":
                report.api_inchikey += 1
            elif confidence == "pubchem_inchikey":
                report.pubchem_inchikey += 1
            elif confidence == "name_pubchem":
                report.name_pubchem += 1
            else:
                report.unmatched += 1
            if hmdb_id:
                report.hmdb_matched += 1

        log.info("%s", report)
        return nodes, report


# ---------------------------------------------------------------------------
# API helpers
# ---------------------------------------------------------------------------

def _normalise_hmdb(hmdb_id: str) -> str:
    """Normalise HMDB IDs to 7-digit zero-padded format (HMDB0000001) for VMH compatibility."""
    hmdb_id = str(hmdb_id).strip().upper()
    if hmdb_id.startswith("HMDB"):
        digits = hmdb_id[4:]
        return f"HMDB{digits.zfill(7)}"
    return hmdb_id


def _chebi_from_inchikey_api(inchikey: str) -> Optional[str]:
    """Query ChEBI OLS4 REST for a ChEBI ID by InChIKey. Returns first match or None."""
    try:
        resp = requests.get(
            _OLS4_SEARCH,
            params={
                "q": inchikey,
                "ontology": "chebi",
                "fieldList": "id,obo_id",
                "exact": "true",
                "rows": 1,
            },
            timeout=15,
        )
        resp.raise_for_status()
        docs = resp.json().get("response", {}).get("docs", [])
        if docs:
            obo_id = docs[0].get("obo_id", "")
            if obo_id.startswith("CHEBI:"):
                return obo_id
    except Exception as exc:
        log.debug("ChEBI OLS4 lookup failed for %s: %s", inchikey, exc)
    return None


def _inchikey_from_name_pubchem(name: str) -> Optional[tuple[str, str]]:
    """
    Query PubChem PUG REST by compound name.

    Returns (InChIKey, CID) for the top hit, or None if not found.
    Takes the first (best-ranked) result only.
    """
    import urllib.parse
    try:
        encoded = urllib.parse.quote(name, safe="")
        url = _PUBCHEM_NAME.format(name=encoded)
        resp = requests.get(url, timeout=15)
        if resp.status_code == 404:
            return None
        resp.raise_for_status()
        props = resp.json().get("PropertyTable", {}).get("Properties", [])
        if props:
            ik  = props[0].get("InChIKey", "")
            cid = str(props[0].get("CID", ""))
            if ik:
                return ik, cid
    except Exception as exc:
        log.debug("PubChem name lookup failed for %r: %s", name, exc)
    return None


def _inchikey_from_pubchem(cid: str) -> Optional[str]:
    """Query PubChem PUG REST to get InChIKey for a CID."""
    try:
        url = _PUBCHEM_REST.format(cid=cid.strip())
        resp = requests.get(url, timeout=15)
        resp.raise_for_status()
        props = resp.json().get("PropertyTable", {}).get("Properties", [])
        if props:
            return props[0].get("InChIKey")
    except Exception as exc:
        log.debug("PubChem lookup failed for CID %s: %s", cid, exc)
    return None


def _clean(val: object) -> Optional[str]:
    if val is None or (isinstance(val, float) and pd.isna(val)):
        return None
    # Convert float-encoded integers (e.g. 5793.0) to int strings
    if isinstance(val, float) and val == int(val):
        val = int(val)
    s = str(val).strip()
    return s if s else None


def _is_numeric(val: object) -> bool:
    try:
        float(val)  # type: ignore[arg-type]
        return True
    except (TypeError, ValueError):
        return False
