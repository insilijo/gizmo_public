"""
MetaNetX FTP/REST client for cross-referencing and ID mapping.

License: MetaNetX data is CC BY 4.0 — https://www.metanetx.org/mnxdoc/mnxref.html

MetaNetX provides flat TSV files via FTP and a REST API.
Key files used:
  chem_xref.tsv   — maps external IDs (ChEBI, KEGG, HMDB, …) → MNX IDs
  reac_xref.tsv   — maps external reaction IDs → MNXR IDs
  chem_prop.tsv   — name, formula, charge, InChIKey for each MNX compound
  reac_prop.tsv   — stoichiometry, direction, EC for each MNXR reaction

We download and cache these files locally; they are released under CC BY 4.0
and do NOT carry KEGG or HMDB downstream restrictions.
"""

from __future__ import annotations

import gzip
import logging
import os
from pathlib import Path
from typing import Iterator, Optional
from urllib.request import urlretrieve

import pandas as pd

log = logging.getLogger(__name__)

_FTP_BASE = "https://www.metanetx.org/cgi-bin/mnxget/mnxref/"
_FILES = {
    "chem_xref": "chem_xref.tsv",
    "reac_xref": "reac_xref.tsv",
    "chem_prop": "chem_prop.tsv",
    "reac_prop": "reac_prop.tsv",
}


def _mnx_header_info(path: Path) -> tuple[list[str], int]:
    """
    Scan a MetaNetX TSV and return (column_names, data_start_line).

    MetaNetX files prefix their column-header row with '#' (e.g. '#ID\\t...')
    alongside hundreds of other '#'-comment lines.  This function finds the
    last '#'-prefixed line with ≥3 non-empty tab-separated fields whose first
    field is a bare word — that is the column header.
    """
    header_cols: list[str] = []
    data_start: int = 0
    with open(path, encoding="utf-8") as fh:
        for lineno, line in enumerate(fh):
            if line.startswith("#"):
                parts = line[1:].rstrip("\n").split("\t")
                non_empty = [p for p in parts if p.strip()]
                if len(non_empty) >= 3 and parts[0] and " " not in parts[0]:
                    header_cols = parts
                    data_start = lineno + 1
            else:
                break
    return header_cols, data_start


def _read_mnx_tsv(path: Path) -> pd.DataFrame:
    """
    Read a MetaNetX flat-file TSV correctly.

    MetaNetX files have hundreds of '#'-prefixed comment lines followed by a
    column header that is *also* '#'-prefixed (e.g. '#ID\\tname\\t...' or
    '#source\\tID\\t...').  Using pd.read_csv(comment='#') drops that header,
    so the first data row becomes the column names — producing garbage.
    """
    header_cols, data_start = _mnx_header_info(path)
    if not header_cols:
        return pd.read_csv(path, sep="\t", comment="#", header=0, low_memory=False)
    return pd.read_csv(
        path, sep="\t", skiprows=data_start, header=None, names=header_cols, low_memory=False
    )


class MetaNetXClient:
    """
    Downloads and parses MetaNetX reference files for ID mapping and
    stoichiometry enrichment.
    """

    def __init__(self, cache_dir: str | Path = "data/raw/metanetx") -> None:
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self._chem_xref: Optional[pd.DataFrame] = None
        self._chem_prop: Optional[pd.DataFrame] = None
        self._reac_xref: Optional[pd.DataFrame] = None
        self._reac_prop: Optional[pd.DataFrame] = None

    # ------------------------------------------------------------------
    # Download helpers
    # ------------------------------------------------------------------

    def download(self, force: bool = False) -> None:
        """Download all reference TSVs if not already cached."""
        for key, fname in _FILES.items():
            dest = self.cache_dir / fname
            if dest.exists() and not force:
                log.info("Cache hit: %s", dest)
                continue
            url = _FTP_BASE + fname
            log.info("Downloading %s → %s", url, dest)
            urlretrieve(url, dest)

    def _load(self, key: str) -> pd.DataFrame:
        path = self.cache_dir / _FILES[key]
        if not path.exists():
            raise FileNotFoundError(f"MetaNetX file not found: {path}. Run .download() first.")
        return _read_mnx_tsv(path)

    # ------------------------------------------------------------------
    # Compound cross-references
    # ------------------------------------------------------------------

    @property
    def chem_xref(self) -> pd.DataFrame:
        if self._chem_xref is None:
            df = self._load("chem_xref")
            # columns: source, ID, mnx_id, description
            self._chem_xref = df
        return self._chem_xref

    @property
    def chem_prop(self) -> pd.DataFrame:
        if self._chem_prop is None:
            self._chem_prop = self._load("chem_prop")
        return self._chem_prop

    def chebi_to_mnx(self, chebi_ids: list[str]) -> dict[str, str]:
        """Map ChEBI IDs to MetaNetX MNX compound IDs."""
        xref = self.chem_xref
        # source column looks like "chebi:CHEBI:15422" or "chebi:15422"
        chebi_col = xref[xref["source"].str.startswith("chebi:", na=False)].copy()
        chebi_col["chebi_id"] = chebi_col["source"].str.replace(
            r"^chebi:(?:CHEBI:)?", "CHEBI:", regex=True
        )
        mapping = chebi_col.set_index("chebi_id")["ID"].to_dict()
        return {cid: mapping[cid] for cid in chebi_ids if cid in mapping}

    def mnx_compound_properties(self, mnx_ids: list[str]) -> pd.DataFrame:
        """Return chem_prop rows for a list of MNX compound IDs."""
        prop = self.chem_prop
        id_col = prop.columns[0]  # "ID" or "mnx_id" depending on version
        return prop[prop[id_col].isin(mnx_ids)].copy()

    # ------------------------------------------------------------------
    # Reaction cross-references
    # ------------------------------------------------------------------

    @property
    def reac_xref(self) -> pd.DataFrame:
        if self._reac_xref is None:
            self._reac_xref = self._load("reac_xref")
        return self._reac_xref

    @property
    def reac_prop(self) -> pd.DataFrame:
        if self._reac_prop is None:
            self._reac_prop = self._load("reac_prop")
        return self._reac_prop

    def chebi_to_hmdb(self, chebi_ids: list[str] | None = None) -> dict[str, str]:
        """
        Map ChEBI IDs to VMH-compatible HMDB IDs (7-digit zero-padded).

        Uses chem_xref to build ChEBI → MNX → HMDB in one DataFrame pass.
        If chebi_ids is None, returns the full mapping for all ChEBI entries.
        """
        xref = self.chem_xref
        chebi_filter = set(chebi_ids) if chebi_ids is not None else None

        # MNX → ChEBI
        chebi_rows = xref[xref["source"].str.startswith("chebi:", na=False)].copy()
        chebi_rows["chebi_id"] = chebi_rows["source"].str.replace(
            r"^chebi:(?:CHEBI:)?", "CHEBI:", regex=True
        )
        if chebi_filter:
            chebi_rows = chebi_rows[chebi_rows["chebi_id"].isin(chebi_filter)]
        mnx_to_chebi: dict[str, str] = chebi_rows.set_index("ID")["chebi_id"].to_dict()

        if not mnx_to_chebi:
            return {}

        # MNX → HMDB (keep first/canonical per MNX)
        hmdb_rows = xref[
            xref["source"].str.startswith("hmdb:", na=False)
            & xref["ID"].isin(mnx_to_chebi)
        ].copy()

        def _norm(s: str) -> str:
            s = s[5:].strip().upper()   # strip "hmdb:"
            return f"HMDB{s[4:].zfill(7)}" if s.startswith("HMDB") else s

        hmdb_rows["hmdb_id"] = hmdb_rows["source"].apply(_norm)
        mnx_to_hmdb: dict[str, str] = (
            hmdb_rows.groupby("ID")["hmdb_id"].first().to_dict()
        )

        # Join: ChEBI → MNX → HMDB (first hit per ChEBI wins)
        result: dict[str, str] = {}
        for mnx_id, chebi in mnx_to_chebi.items():
            if mnx_id in mnx_to_hmdb and chebi not in result:
                result[chebi] = mnx_to_hmdb[mnx_id]
        return result

    def reactome_to_hmdb(self, reactome_ids: list[str] | None = None) -> dict[str, str]:
        """Map Reactome stIDs (e.g. ``R-ALL-29438``) to HMDB IDs via MetaNetX.

        Uses ``chem_xref.tsv`` entries with ``reactome:R-ALL-XXXXX`` source
        values to route: Reactome stId → MNX compound ID → HMDB ID.

        Parameters
        ----------
        reactome_ids:
            List of bare Reactome stIds (without ``reactome:`` prefix).
            If None, builds the full mapping for all Reactome entries.

        Returns
        -------
        ``{reactome_stid: hmdb_id_normalised}``
        """
        xref = self.chem_xref
        reactome_rows = xref[xref["source"].str.startswith("reactome:", na=False)].copy()
        reactome_rows["reactome_stid"] = reactome_rows["source"].str.replace(
            "reactome:", "", regex=False
        )
        if reactome_ids is not None:
            reactome_rows = reactome_rows[reactome_rows["reactome_stid"].isin(reactome_ids)]

        mnx_to_reactome: dict[str, str] = (
            reactome_rows.groupby("ID")["reactome_stid"].first().to_dict()
        )
        # Invert: reactome stid → MNX
        reactome_to_mnx: dict[str, str] = {v: k for k, v in mnx_to_reactome.items()}

        # MNX → HMDB
        hmdb_rows = xref[
            xref["source"].str.startswith("hmdb:", na=False)
            & xref["ID"].isin(reactome_to_mnx.values())
        ].copy()

        def _norm(s: str) -> str:
            s = s[5:].strip().upper()  # strip "hmdb:"
            return f"HMDB{s[4:].zfill(7)}" if s.startswith("HMDB") else s

        hmdb_rows["hmdb_id"] = hmdb_rows["source"].apply(_norm)
        mnx_to_hmdb: dict[str, str] = (
            hmdb_rows.groupby("ID")["hmdb_id"].first().to_dict()
        )

        result: dict[str, str] = {}
        for stid, mnx_id in reactome_to_mnx.items():
            if mnx_id in mnx_to_hmdb:
                result[stid] = mnx_to_hmdb[mnx_id]
        return result

    def reactome_to_chebi(self, reactome_ids: list[str] | None = None) -> dict[str, str]:
        """Map Reactome stIds to ChEBI IDs via MetaNetX.

        ReactomeLoader doesn't extract ChEBI from participant JSON because
        SimpleEntity entries omit ``crossReference``; this fills the gap by
        joining Reactome stId → MNX compound ID → ChEBI through chem_xref.

        Parameters
        ----------
        reactome_ids:
            Bare Reactome stIds (without ``reactome:`` prefix). If None,
            returns the full mapping for all Reactome entries with a ChEBI
            cross-reference.

        Returns
        -------
        ``{reactome_stid: "CHEBI:NNNNN"}``  (first ChEBI per stId wins)
        """
        xref = self.chem_xref

        # Reactome stId → MNX
        reactome_rows = xref[xref["source"].str.startswith("reactome:", na=False)].copy()
        reactome_rows["reactome_stid"] = reactome_rows["source"].str.replace(
            "reactome:", "", regex=False
        )
        if reactome_ids is not None:
            reactome_rows = reactome_rows[reactome_rows["reactome_stid"].isin(reactome_ids)]

        reactome_to_mnx: dict[str, str] = (
            reactome_rows.groupby("reactome_stid")["ID"].first().to_dict()
        )

        # MNX → ChEBI (chem_xref source like "chebi:NNNNN" or "chebi:CHEBI:NNNNN")
        chebi_rows = xref[
            xref["source"].str.startswith("chebi:", na=False)
            & xref["ID"].isin(reactome_to_mnx.values())
        ].copy()
        chebi_rows["chebi_id"] = chebi_rows["source"].str.replace(
            r"^chebi:(?:CHEBI:)?", "CHEBI:", regex=True
        )
        mnx_to_chebi: dict[str, str] = (
            chebi_rows.groupby("ID")["chebi_id"].first().to_dict()
        )

        return {
            stid: mnx_to_chebi[mnx_id]
            for stid, mnx_id in reactome_to_mnx.items()
            if mnx_id in mnx_to_chebi
        }

    def enrich_graph_chebi(self, mg) -> int:
        """Back-fill ``chebi_id`` on metabolite nodes that lack it.

        Reactome's loader omits ChEBI for ``SimpleEntity`` participants
        because their JSON doesn't carry ``crossReference``. We restore
        ChEBI IDs via two routes:

        1. ``reactome_id`` → MetaNetX → ChEBI  (primary)
        2. ``inchikey`` → MetaNetX InChIKey[:14] → MNX → ChEBI  (structural
           backup; requires :meth:`enrich_graph_inchikey` to have run)

        .. warning::
           MetaNetX's ``chem_xref.tsv`` cross-references many entries to
           **obsolete ChEBI IDs** (e.g. CHEBI:10789 for ATP rather than
           the current CHEBI:15422/30616). The IDs populated here will
           NOT match canonical ChEBI lists like
           :data:`gizmo.analysis.currency.KNOWN_CURRENCY_CHEBI` until an
           obsolete→primary remap is layered on top. Use InChIKey-based
           matching where structural identity matters.

        Mutates NetworkX node attributes in place. Returns total enriched.
        """
        enriched = 0

        # ── Route 1: reactome_id → ChEBI ──────────────────────────────────
        needs_reactome = {
            nid: attrs["reactome_id"]
            for nid, attrs in mg.graph.nodes(data=True)
            if attrs.get("node_type") == "metabolite"
            and attrs.get("reactome_id")
            and not attrs.get("chebi_id")
        }
        if needs_reactome:
            reactome_map = self.reactome_to_chebi(list(set(needs_reactome.values())))
            for nid, stid in needs_reactome.items():
                if stid in reactome_map:
                    mg.graph.nodes[nid]["chebi_id"] = reactome_map[stid]
                    enriched += 1

        # ── Route 2: inchikey[:14] → MNX → ChEBI ──────────────────────────
        prop = self.chem_prop
        xref = self.chem_xref
        prop_clean = prop[prop["InChIKey"].notna() & (prop["InChIKey"] != "")]
        mnx_to_ik: dict[str, str] = prop_clean.set_index("ID")["InChIKey"].to_dict()

        chebi_rows = xref[xref["source"].str.startswith("chebi:", na=False)].copy()
        chebi_rows["chebi_id"] = chebi_rows["source"].str.replace(
            r"^chebi:(?:CHEBI:)?", "CHEBI:", regex=True
        )
        mnx_to_chebi: dict[str, str] = (
            chebi_rows.groupby("ID")["chebi_id"].first().to_dict()
        )

        ik14_to_chebi: dict[str, str] = {}
        for mnx_id, ik in mnx_to_ik.items():
            if mnx_id in mnx_to_chebi:
                ik14 = ik[:14]
                ik14_to_chebi.setdefault(ik14, mnx_to_chebi[mnx_id])

        needs_inchikey = {
            nid: attrs["inchikey"]
            for nid, attrs in mg.graph.nodes(data=True)
            if attrs.get("node_type") == "metabolite"
            and attrs.get("inchikey")
            and not attrs.get("chebi_id")
        }
        ik_enriched = 0
        for nid, ik in needs_inchikey.items():
            chebi = ik14_to_chebi.get(str(ik)[:14])
            if chebi:
                mg.graph.nodes[nid]["chebi_id"] = chebi
                enriched += 1
                ik_enriched += 1

        log.info(
            "enrich_graph_chebi: enriched %d nodes total "
            "(reactome route: %d candidates; inchikey route: %d candidates → %d hits)",
            enriched, len(needs_reactome), len(needs_inchikey), ik_enriched,
        )
        return enriched

    def enrich_graph_inchikey(self, mg) -> int:
        """Populate ``inchikey`` on metabolite nodes that lack it.

        Route: ``reactome_id`` → MetaNetX ``chem_xref`` → MNX compound ID
        → ``chem_prop`` InChIKey.

        Mutates NetworkX node attribute dicts directly.  Returns the number
        of nodes enriched.
        """
        xref = self.chem_xref
        prop = self.chem_prop

        # Build reactome_id → MNX compound ID
        react_rows = xref[xref["source"].str.startswith("reactome:", na=False)].copy()
        react_rows["reactome_id"] = react_rows["source"].str.replace(
            "reactome:", "", regex=False
        )
        reactome_to_mnx: dict[str, str] = (
            react_rows.groupby("reactome_id")["ID"].first().to_dict()
        )

        # Build MNX compound ID → InChIKey from chem_prop
        prop_clean = prop[prop["InChIKey"].notna() & (prop["InChIKey"] != "")]
        mnx_to_ik: dict[str, str] = prop_clean.set_index("ID")["InChIKey"].to_dict()

        enriched = 0
        for nid, attrs in mg.graph.nodes(data=True):
            if attrs.get("node_type") != "metabolite":
                continue
            if attrs.get("inchikey"):
                continue
            rid = attrs.get("reactome_id")
            if not rid:
                continue
            mnx_id = reactome_to_mnx.get(rid)
            if not mnx_id:
                continue
            ik = mnx_to_ik.get(mnx_id)
            if ik:
                mg.graph.nodes[nid]["inchikey"] = ik
                enriched += 1

        log.info("enrich_graph_inchikey: populated inchikey on %d nodes", enriched)
        return enriched

    def enrich_graph_hmdb(self, mg) -> int:
        """Back-fill ``hmdb_id`` on metabolite nodes that have ``chebi_id``,
        ``reactome_id``, or ``inchikey`` but no ``hmdb_id``.

        Tries three routes in order:
        1. ``chebi_id`` → MetaNetX → HMDB
        2. ``reactome_id`` → MetaNetX → HMDB  (for Reactome-sourced graphs
           where ChEBI wasn't populated by the loader)
        3. ``inchikey`` → MetaNetX InChIKey[0:14] → MNX → HMDB  (structural
           matching; requires :meth:`enrich_graph_inchikey` to have been run
           first, but also works with any pre-existing ``inchikey`` values)

        Mutates NetworkX node attribute dicts directly.  Returns total
        number of nodes enriched.
        """
        enriched = 0

        # ── Route 1: chebi_id → HMDB ──────────────────────────────────────
        needs_chebi = {
            nid: attrs["chebi_id"]
            for nid, attrs in mg.graph.nodes(data=True)
            if attrs.get("node_type") == "metabolite"
            and attrs.get("chebi_id")
            and not attrs.get("hmdb_id")
        }
        if needs_chebi:
            chebi_map = self.chebi_to_hmdb(list(set(needs_chebi.values())))
            for nid, chebi in needs_chebi.items():
                if chebi in chebi_map:
                    mg.graph.nodes[nid]["hmdb_id"] = chebi_map[chebi]
                    enriched += 1

        # ── Route 2: reactome_id → HMDB ───────────────────────────────────
        needs_reactome = {
            nid: attrs["reactome_id"]
            for nid, attrs in mg.graph.nodes(data=True)
            if attrs.get("node_type") == "metabolite"
            and attrs.get("reactome_id")
            and not attrs.get("hmdb_id")
        }
        if needs_reactome:
            reactome_map = self.reactome_to_hmdb(list(set(needs_reactome.values())))
            for nid, stid in needs_reactome.items():
                if stid in reactome_map:
                    mg.graph.nodes[nid]["hmdb_id"] = reactome_map[stid]
                    enriched += 1

        # ── Route 3: inchikey → MNX → HMDB ───────────────────────────────
        # Build InChIKey[0:14] → HMDB from chem_prop + chem_xref
        prop = self.chem_prop
        xref = self.chem_xref
        prop_clean = prop[prop["InChIKey"].notna() & (prop["InChIKey"] != "")]
        mnx_to_ik: dict[str, str] = prop_clean.set_index("ID")["InChIKey"].to_dict()

        hmdb_rows = xref[xref["source"].str.startswith("hmdb:", na=False)].copy()

        def _norm_hmdb(s: str) -> str:
            s = s[5:].strip().upper()
            return f"HMDB{s[4:].zfill(7)}" if s.startswith("HMDB") else s

        hmdb_rows["hmdb_id"] = hmdb_rows["source"].apply(_norm_hmdb)
        mnx_to_hmdb: dict[str, str] = (
            hmdb_rows.groupby("ID")["hmdb_id"].first().to_dict()
        )

        ik14_to_hmdb: dict[str, str] = {}
        for mnx_id, ik in mnx_to_ik.items():
            if mnx_id in mnx_to_hmdb:
                ik14 = ik[:14]
                if ik14 not in ik14_to_hmdb:
                    ik14_to_hmdb[ik14] = mnx_to_hmdb[mnx_id]

        needs_inchikey = {
            nid: attrs["inchikey"]
            for nid, attrs in mg.graph.nodes(data=True)
            if attrs.get("node_type") == "metabolite"
            and attrs.get("inchikey")
            and not attrs.get("hmdb_id")
        }
        inchikey_enriched = 0
        for nid, ik in needs_inchikey.items():
            hmdb = ik14_to_hmdb.get(str(ik)[:14])
            if hmdb:
                mg.graph.nodes[nid]["hmdb_id"] = hmdb
                enriched += 1
                inchikey_enriched += 1

        log.info(
            "enrich_graph_hmdb: enriched %d nodes "
            "(chebi: %d candidates, reactome: %d candidates, inchikey: %d candidates → %d hits)",
            enriched, len(needs_chebi), len(needs_reactome),
            len(needs_inchikey), inchikey_enriched,
        )
        return enriched

    def reactome_to_mnxr(self, reactome_stids: list[str]) -> dict[str, str]:
        """Map Reactome stIDs to MetaNetX MNXR reaction IDs."""
        xref = self.reac_xref
        reactome_rows = xref[xref["source"].str.startswith("reactome:", na=False)].copy()
        reactome_rows["reactome_stid"] = reactome_rows["source"].str.replace(
            "reactome:", "", regex=False
        )
        mapping = reactome_rows.set_index("reactome_stid")["ID"].to_dict()
        return {sid: mapping[sid] for sid in reactome_stids if sid in mapping}
