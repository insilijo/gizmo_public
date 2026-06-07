"""
Metabolon compound curation interface.

Resolves unmatched Metabolon compounds (those that didn't hit MetaNetX InChIKey
lookup) by providing ChEBI and VMH-compatible HMDB search, manual
assignment, graph connection verification, and a Jupyter widget UI.

Usage in a notebook::

    from gizmo.curation.metabolon_curator import MetabolonCurator

    curator = MetabolonCurator(
        loader=met_loader,
        graph=mg,
        overrides_path='data/resources/gizmo/curation/metabolon_overrides.json',
    )
    curator.widget()          # interactive Jupyter widget

    # Batch workflow:
    curator.export_csv('data/curation/unmatched.csv')
    # ... fill in CHEBI_CURATED / HMDB_CURATED columns ...
    curator.import_csv('data/curation/unmatched.csv')
    curator.save()
    curator.apply()           # push assignments into loader
"""

from __future__ import annotations

import json
import logging
import time
from dataclasses import asdict, dataclass, field
from datetime import date
from pathlib import Path
from typing import Optional

import pandas as pd
import requests

from gizmo.resources import overrides_default
from gizmo.sources.metabolon import MetabolonLoader

log = logging.getLogger(__name__)

_OLS4_SEARCH = "https://www.ebi.ac.uk/ols4/api/search"


# ---------------------------------------------------------------------------
# ChEBI name search helper
# ---------------------------------------------------------------------------

def _search_chebi_name(query: str, n: int = 5) -> list[dict]:
    """Search ChEBI OLS4 by compound name. Returns [{chebi_id, name, score}]."""
    try:
        resp = requests.get(
            _OLS4_SEARCH,
            params={
                "q": query,
                "ontology": "chebi",
                "fieldList": "obo_id,label,score",
                "rows": n,
            },
            timeout=15,
        )
        resp.raise_for_status()
        results = []
        for doc in resp.json().get("response", {}).get("docs", []):
            obo_id = doc.get("obo_id", "")
            if obo_id.startswith("CHEBI:"):
                results.append({
                    "chebi_id": obo_id,
                    "name": doc.get("label", ""),
                    "score": round(float(doc.get("score", 0)), 3),
                })
        return results
    except Exception as exc:
        log.debug("ChEBI name search failed for %r: %s", query, exc)
        return []


# ---------------------------------------------------------------------------
# Override data model
# ---------------------------------------------------------------------------

@dataclass
class CurationEntry:
    chebi_id: Optional[str] = None
    hmdb_id: Optional[str] = None
    inchikey: Optional[str] = None
    confidence: str = "manual"  # manual | api_name | api_pubchem
    notes: str = ""
    excluded: bool = False
    assigned_at: str = field(default_factory=lambda: str(date.today()))
    # Edge-level curation: reaction stIDs that are explicitly excluded or approved.
    # Empty excluded_reactions → all inherited edges are kept.
    # Non-empty → those reaction node IDs are removed from the graph on apply().
    excluded_reactions: list[str] = field(default_factory=list)
    approved_reactions: list[str] = field(default_factory=list)  # informational


# ---------------------------------------------------------------------------
# Curator
# ---------------------------------------------------------------------------

class MetabolonCurator:
    """
    Interactive curation tool for unmatched Metabolon compounds.

    Wraps a MetabolonLoader and adds:
    - Persistent manual ChEBI / HMDB assignments (JSON overrides file)
    - ChEBI name search via OLS4
    - HMDB/VMH search via graph metabolite nodes
    - Graph connection verification
    - Jupyter widget UI
    - CSV export / import for batch review

    After curation, call .apply() to push assignments into the loader so that
    loader.to_metabolite_nodes() reflects the manual assignments.
    """

    def __init__(
        self,
        loader: MetabolonLoader,
        graph=None,                  # Optional[GizmoGraph]
        overrides_path: str | Path = overrides_default(),
    ) -> None:
        self.loader = loader
        self.graph = graph
        self.overrides_path = Path(overrides_path)
        self.overrides: dict[str, CurationEntry] = {}
        self._load()

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------

    def _load(self) -> None:
        if self.overrides_path.exists():
            raw = json.loads(self.overrides_path.read_text())
            for name, data in raw.get("assignments", {}).items():
                self.overrides[name] = CurationEntry(**data)
            log.info(
                "Loaded %d curation overrides from %s",
                len(self.overrides), self.overrides_path,
            )

    def save(self) -> None:
        """Persist all overrides to the JSON file."""
        self.overrides_path.parent.mkdir(parents=True, exist_ok=True)
        data = {"assignments": {k: asdict(v) for k, v in self.overrides.items()}}
        self.overrides_path.write_text(json.dumps(data, indent=2))
        log.info("Saved %d curation overrides → %s", len(self.overrides), self.overrides_path)

    # ------------------------------------------------------------------
    # Status / reporting
    # ------------------------------------------------------------------

    def _is_auto_matched(self, row: pd.Series) -> bool:
        """Return True if the MetabolonLoader auto-matches this row without overrides."""
        ik = str(row.get("INCHIKEY", "")).strip()
        if not ik or ik == "nan":
            return False
        if ik in self.loader._inchikey_to_chebi:
            return True
        if len(ik) >= 14 and ik[:14] in self.loader._connectivity_to_chebi:
            return True
        return False

    @property
    def unmatched(self) -> pd.DataFrame:
        """Compounds with no auto-match and no accepted manual override."""
        df = self.loader._load_df().copy()
        accepted = {k for k, v in self.overrides.items() if not v.excluded}
        rows = []
        for _, row in df.iterrows():
            name = str(row.get("BIOCHEMICAL", "")).strip()
            if name in accepted:
                continue
            if name in self.overrides and self.overrides[name].excluded:
                continue
            if self._is_auto_matched(row):
                continue
            rows.append(row)
        return pd.DataFrame(rows).reset_index(drop=True)

    @property
    def assigned(self) -> pd.DataFrame:
        """All manually accepted assignments."""
        entries = [
            {"metabolon_name": k, **asdict(v)}
            for k, v in self.overrides.items()
            if not v.excluded
        ]
        return pd.DataFrame(entries)

    @property
    def excluded(self) -> list[str]:
        return [k for k, v in self.overrides.items() if v.excluded]

    def summary(self) -> None:
        df = self.loader._load_df()
        total = len(df)
        auto = sum(1 for _, row in df.iterrows() if self._is_auto_matched(row))
        manual = sum(
            1 for v in self.overrides.values()
            if not v.excluded and (v.chebi_id or v.hmdb_id)
        )
        excl = len(self.excluded)
        remaining = total - auto - manual - excl
        coverage = (auto + manual) / total * 100
        print(f"Total compounds:          {total}")
        print(f"  Auto-matched (MetaNetX): {auto}")
        print(f"  Manually curated:        {manual}")
        print(f"  Excluded:                {excl}")
        print(f"  Still unmatched:         {remaining}")
        print(f"  ChEBI/HMDB coverage:     {coverage:.1f}%")

    # ------------------------------------------------------------------
    # Assignment
    # ------------------------------------------------------------------

    def assign(
        self,
        metabolon_name: str,
        chebi_id: str | None = None,
        hmdb_id: str | None = None,
        notes: str = "",
        inchikey: str = "",
        confidence: str = "manual",
    ) -> None:
        """Manually assign a ChEBI and/or HMDB ID to a Metabolon compound."""
        chebi = (chebi_id or "").strip() or None
        hmdb = (hmdb_id or "").strip().upper() or None
        if hmdb and hmdb.startswith("HMDB"):
            hmdb = f"HMDB{hmdb[4:].zfill(7)}"
        if not (chebi or hmdb):
            raise ValueError("assign() requires a chebi_id and/or hmdb_id")
        self.overrides[metabolon_name] = CurationEntry(
            chebi_id=chebi,
            hmdb_id=hmdb,
            inchikey=inchikey.strip() or None,
            notes=notes,
            confidence=confidence,
        )

    def exclude(self, metabolon_name: str, notes: str = "") -> None:
        """Mark a compound as intentionally excluded (no ChEBI mapping expected)."""
        self.overrides[metabolon_name] = CurationEntry(excluded=True, notes=notes)

    def unapply(self, metabolon_name: str) -> None:
        """Remove a manual override."""
        self.overrides.pop(metabolon_name, None)

    # ------------------------------------------------------------------
    # ChEBI search
    # ------------------------------------------------------------------

    def search_chebi(self, query: str, n: int = 5) -> list[dict]:
        """
        Search ChEBI by name and annotate results with graph connection info.
        Returns [{chebi_id, name, score, in_graph, reaction_count}].
        """
        results = _search_chebi_name(query, n)
        for r in results:
            info = self.verify_connection(r["chebi_id"])
            r["in_graph"] = info["in_graph"]
            r["reaction_count"] = info["reaction_count"]
        return results

    def search_hmdb(self, query: str, n: int = 5) -> list[dict]:
        """Search graph metabolite nodes for VMH-compatible HMDB IDs by name."""
        if self.graph is None:
            return []
        q = query.strip().lower()
        if not q:
            return []
        hits: list[dict] = []
        seen: set[tuple[str, str]] = set()
        for node_id, attrs in self.graph.graph.nodes(data=True):
            if attrs.get("node_type") != "metabolite":
                continue
            hmdb_id = str(attrs.get("hmdb_id") or "").strip().upper()
            if not hmdb_id:
                continue
            name = str(attrs.get("name") or attrs.get("metabolon_name") or node_id)
            haystacks = [name.lower(), hmdb_id.lower(), str(attrs.get("chebi_id") or "").lower()]
            if not any(q in h for h in haystacks):
                continue
            key = (hmdb_id, name)
            if key in seen:
                continue
            seen.add(key)
            info = self.verify_hmdb_connection(hmdb_id)
            hits.append({
                "hmdb_id": hmdb_id,
                "name": name,
                "chebi_id": attrs.get("chebi_id") or "",
                "node_id": node_id,
                "in_graph": info["in_graph"],
                "reaction_count": info["reaction_count"],
            })
            if len(hits) >= n:
                break
        return hits

    # ------------------------------------------------------------------
    # Graph verification
    # ------------------------------------------------------------------

    def verify_connection(self, chebi_id: str) -> dict:
        """
        Check whether a ChEBI ID is present in the reaction graph.
        Returns {in_graph, reaction_count, pathways}.
        """
        if self.graph is None:
            return {"in_graph": None, "reaction_count": None, "pathways": []}
        g = self.graph.graph
        if chebi_id not in g:
            return {"in_graph": False, "reaction_count": 0, "pathways": []}
        nbrs = list(g.predecessors(chebi_id)) + list(g.successors(chebi_id))
        rxns = [n for n in nbrs if g.nodes[n].get("node_type") == "reaction"]
        pathways = []
        for rxn in rxns[:10]:
            pathways.extend(g.nodes[rxn].get("pathways") or [])
        return {
            "in_graph": True,
            "reaction_count": len(rxns),
            "pathways": list(dict.fromkeys(pathways)),  # deduplicated, ordered
        }

    def verify_hmdb_connection(self, hmdb_id: str) -> dict:
        """Check whether an HMDB ID is present on any metabolite node in the graph."""
        if self.graph is None:
            return {"in_graph": None, "reaction_count": None, "pathways": [], "node_ids": []}
        g = self.graph.graph
        target = hmdb_id.strip().upper()
        matches = [
            nid for nid, attrs in g.nodes(data=True)
            if attrs.get("node_type") == "metabolite"
            and str(attrs.get("hmdb_id") or "").strip().upper() == target
        ]
        if not matches:
            return {"in_graph": False, "reaction_count": 0, "pathways": [], "node_ids": []}
        rxns: list[str] = []
        pathways: list[str] = []
        for nid in matches:
            nbrs = list(g.predecessors(nid)) + list(g.successors(nid))
            local_rxns = [n for n in nbrs if g.nodes[n].get("node_type") == "reaction"]
            rxns.extend(local_rxns)
            for rxn in local_rxns[:10]:
                pathways.extend(g.nodes[rxn].get("pathways") or [])
        return {
            "in_graph": True,
            "reaction_count": len(set(rxns)),
            "pathways": list(dict.fromkeys(pathways)),
            "node_ids": matches,
        }

    # ------------------------------------------------------------------
    # Edge review
    # ------------------------------------------------------------------

    def reaction_connections(self, metabolon_name: str) -> list[dict]:
        """
        Return the reaction nodes the curated compound would connect to if its
        assigned ChEBI ID is already in the graph.

        Returns a list of dicts with keys:
            reaction_id, name, role, pathways, ec_numbers, gene_symbols
        """
        if self.graph is None:
            return []
        entry = self.overrides.get(metabolon_name)
        chebi = (entry.chebi_id if entry else None) or self.loader._name_overrides.get(metabolon_name)
        if not chebi or chebi not in self.graph.graph:
            return []
        g = self.graph.graph
        results = []
        for nbr in list(g.predecessors(chebi)) + list(g.successors(chebi)):
            if g.nodes[nbr].get("node_type") != "reaction":
                continue
            edge_data = g.get_edge_data(nbr, chebi) or g.get_edge_data(chebi, nbr) or {}
            results.append({
                "reaction_id":   nbr,
                "name":          g.nodes[nbr].get("name", nbr),
                "role":          edge_data.get("role", ""),
                "pathways":      g.nodes[nbr].get("pathways") or [],
                "ec_numbers":    g.nodes[nbr].get("ec_numbers") or [],
                "gene_symbols":  g.nodes[nbr].get("gene_symbols") or [],
                "excluded":      nbr in (entry.excluded_reactions if entry else []),
                "approved":      nbr in (entry.approved_reactions if entry else []),
            })
        return results

    def exclude_reaction(self, metabolon_name: str, reaction_id: str) -> None:
        """Mark a reaction as excluded for a curated compound."""
        entry = self.overrides.get(metabolon_name)
        if entry is None:
            raise KeyError(f"No curation entry for {metabolon_name!r} — assign() first")
        if reaction_id not in entry.excluded_reactions:
            entry.excluded_reactions.append(reaction_id)
        if reaction_id in entry.approved_reactions:
            entry.approved_reactions.remove(reaction_id)

    def approve_reaction(self, metabolon_name: str, reaction_id: str) -> None:
        """Explicitly approve a reaction for a curated compound (informational)."""
        entry = self.overrides.get(metabolon_name)
        if entry is None:
            raise KeyError(f"No curation entry for {metabolon_name!r} — assign() first")
        if reaction_id in entry.excluded_reactions:
            entry.excluded_reactions.remove(reaction_id)
        if reaction_id not in entry.approved_reactions:
            entry.approved_reactions.append(reaction_id)

    def restore_reaction(self, metabolon_name: str, reaction_id: str) -> None:
        """Remove a reaction from both excluded and approved lists (reset to default)."""
        entry = self.overrides.get(metabolon_name)
        if entry is None:
            return
        entry.excluded_reactions = [r for r in entry.excluded_reactions if r != reaction_id]
        entry.approved_reactions = [r for r in entry.approved_reactions if r != reaction_id]

    # ------------------------------------------------------------------
    # Apply overrides to loader
    # ------------------------------------------------------------------

    def apply(self) -> int:
        """
        Push manual assignments into loader._name_overrides so that
        loader.to_metabolite_nodes() picks them up.  Also removes any
        explicitly excluded reaction edges from the graph in-place.
        Returns number of node assignments applied.
        """
        count = 0
        for name, entry in self.overrides.items():
            if entry.excluded or not (entry.chebi_id or entry.hmdb_id):
                continue
            if entry.chebi_id:
                self.loader._name_overrides[name] = entry.chebi_id
            else:
                self.loader._name_overrides.pop(name, None)
            if entry.hmdb_id:
                self.loader._name_hmdb_overrides[name] = entry.hmdb_id
            else:
                self.loader._name_hmdb_overrides.pop(name, None)
            count += 1

            # Remove explicitly excluded reaction edges
            if entry.excluded_reactions and self.graph is not None and entry.chebi_id:
                chebi = entry.chebi_id
                g = self.graph.graph
                if chebi in g:
                    for rxn_id in entry.excluded_reactions:
                        if g.has_edge(chebi, rxn_id):
                            g.remove_edge(chebi, rxn_id)
                            log.debug("Removed edge %s → %s (excluded by curator)", chebi, rxn_id)
                        if g.has_edge(rxn_id, chebi):
                            g.remove_edge(rxn_id, chebi)
                            log.debug("Removed edge %s → %s (excluded by curator)", rxn_id, chebi)

        log.info("Applied %d manual overrides to loader", count)
        return count

    # ------------------------------------------------------------------
    # CSV batch workflow
    # ------------------------------------------------------------------

    def export_csv(self, path: str | Path, include_edge_review: bool = True) -> None:
        """
        Export unmatched compounds to CSV for offline review.

        Columns: BIOCHEMICAL, PATHWAY, INCHIKEY, PUBCHEM, MASS,
                 CHEBI_CURATED, HMDB_CURATED, NOTES,
                 REACTION_CONNECTIONS (pipe-delimited, read-only reference),
                 EXCLUDED_REACTIONS (pipe-delimited, editable)
        """
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        df = self.unmatched[["BIOCHEMICAL", "PATHWAY", "INCHIKEY", "PUBCHEM", "MASS"]].copy()
        df["CHEBI_CURATED"] = ""
        df["HMDB_CURATED"] = ""
        df["NOTES"] = ""
        if include_edge_review and self.graph is not None:
            rxn_cols = []
            excl_cols = []
            for _, row in df.iterrows():
                name = str(row["BIOCHEMICAL"]).strip()
                conns = self.reaction_connections(name)
                rxn_cols.append("|".join(c["reaction_id"] for c in conns))
                entry = self.overrides.get(name)
                excl_cols.append("|".join(entry.excluded_reactions if entry else []))
            df["REACTION_CONNECTIONS"] = rxn_cols
            df["EXCLUDED_REACTIONS"] = excl_cols
        df.to_csv(path, index=False)
        print(f"Exported {len(df)} unmatched compounds → {path}")
        print("Fill in CHEBI_CURATED and/or HMDB_CURATED (or write 'excluded').")
        if include_edge_review:
            print("To exclude reactions: add pipe-delimited reaction IDs to EXCLUDED_REACTIONS.")

    def import_csv(self, path: str | Path) -> int:
        """
        Import curated CSV back.

        Reads CHEBI_CURATED, HMDB_CURATED (write 'excluded' to mark a compound as
        intentionally unmatched) and EXCLUDED_REACTIONS (pipe-delimited reaction IDs).
        Returns number of assignments imported.
        """
        df = pd.read_csv(path)
        count = 0
        for _, row in df.iterrows():
            name = str(row.get("BIOCHEMICAL", "")).strip()
            chebi = str(row.get("CHEBI_CURATED", "")).strip()
            hmdb = str(row.get("HMDB_CURATED", "")).strip()
            notes = str(row.get("NOTES", "")).strip()
            excl_rxns_raw = str(row.get("EXCLUDED_REACTIONS", "")).strip()
            if not name:
                continue
            chebi_norm = chebi.lower()
            hmdb_norm = hmdb.lower()
            if chebi_norm == "excluded" or hmdb_norm == "excluded":
                self.exclude(name, notes=notes)
            elif chebi_norm in ("", "nan") and hmdb_norm in ("", "nan"):
                # Only edge exclusions, no node assignment change
                pass
            else:
                self.assign(
                    name,
                    chebi_id=None if chebi_norm in ("", "nan") else chebi,
                    hmdb_id=None if hmdb_norm in ("", "nan") else hmdb,
                    notes=notes,
                )
                count += 1

            # Apply reaction exclusions from the EXCLUDED_REACTIONS column
            if excl_rxns_raw and excl_rxns_raw not in ("", "nan") and name in self.overrides:
                for rxn_id in excl_rxns_raw.split("|"):
                    rxn_id = rxn_id.strip()
                    if rxn_id:
                        self.exclude_reaction(name, rxn_id)

        print(f"Imported {count} curation entries from {path}")
        return count

    # ------------------------------------------------------------------
    # Jupyter widget
    # ------------------------------------------------------------------

    def widget(self):
        """Launch the interactive curation widget in a Jupyter notebook."""
        try:
            import ipywidgets as widgets
            from IPython.display import display
        except ImportError:
            print("ipywidgets not installed. Run:  pip install ipywidgets")
            print("Falling back to text preview. Use .export_csv() for batch review.")
            self._text_preview()
            return None
        return self._build_widget(widgets, display)

    def _text_preview(self, n: int = 5) -> None:
        um = self.unmatched.head(n)
        print(f"=== Unmatched compounds (showing {len(um)} of {len(self.unmatched)}) ===\n")
        for _, row in um.iterrows():
            name = row["BIOCHEMICAL"]
            print(f"Name:     {name}")
            print(f"InChIKey: {row.get('INCHIKEY', '—')}")
            print(f"PubChem:  {row.get('PUBCHEM', '—')}")
            cands = self.search_chebi(name, 3)
            for i, c in enumerate(cands):
                conn = f"{c['reaction_count']} reactions" if c.get("in_graph") else "not in graph"
                print(f"  [{i+1}] {c['chebi_id']} — {c['name']}  ({conn})")
            print()
        print("Use .assign(name, chebi_id=..., hmdb_id=...) / .exclude(name) then .save() and .apply().")

    def _build_widget(self, widgets, display):
        # ---- state ----
        df_all = self.loader._load_df().set_index("BIOCHEMICAL")
        state = {
            "rows": self.unmatched["BIOCHEMICAL"].tolist(),
            "idx": 0,
            "candidates": [],
        }

        # ---- layout helpers ----
        W = widgets.Layout(width="100%")
        BTN = widgets.Layout(width="auto", min_width="90px")

        # ---- widgets ----
        w_progress   = widgets.HTML()
        w_info       = widgets.HTML()
        w_candidates = widgets.VBox()
        w_hmdb_candidates = widgets.VBox()
        w_custom     = widgets.Text(
            placeholder="CHEBI:XXXXX", description="ChEBI ID:",
            layout=widgets.Layout(width="220px"),
        )
        w_custom_hmdb = widgets.Text(
            placeholder="HMDB0000001", description="HMDB ID:",
            layout=widgets.Layout(width="220px"),
        )
        w_notes      = widgets.Text(
            placeholder="optional notes", description="Notes:",
            layout=widgets.Layout(width="380px"),
        )
        w_status     = widgets.Output(layout=widgets.Layout(height="60px", overflow="auto"))
        w_edge_review = widgets.VBox(layout=widgets.Layout(margin="8px 0 0 0"))

        btn_search   = widgets.Button(description="Search ChEBI", button_style="info",   icon="search", layout=BTN)
        btn_search_hmdb = widgets.Button(description="Search VMH/HMDB", button_style="info", icon="search", layout=BTN)
        btn_edges    = widgets.Button(description="Review Edges",  button_style="info",  icon="link",   layout=BTN)
        btn_skip     = widgets.Button(description="Skip",         button_style="warning",               layout=BTN)
        btn_exclude  = widgets.Button(description="Exclude",      button_style="danger",                layout=BTN)
        btn_assign   = widgets.Button(description="Assign",       button_style="success",               layout=BTN)
        btn_assign_hmdb = widgets.Button(description="Assign HMDB", button_style="success",             layout=BTN)
        btn_prev     = widgets.Button(description="◀ Prev",                                             layout=BTN)
        btn_save     = widgets.Button(description="💾 Save",      button_style="primary",               layout=BTN)

        # ---- helpers ----
        def current_name() -> Optional[str]:
            rows = state["rows"]
            i = state["idx"]
            return rows[i] if i < len(rows) else None

        def refresh():
            state["rows"] = self.unmatched["BIOCHEMICAL"].tolist()
            total = len(state["rows"])
            i = state["idx"]
            name = current_name()

            if name is None:
                w_progress.value = "<b style='color:green'>✓ All compounds reviewed!</b>"
                w_info.value = ""
                w_candidates.children = []
                w_hmdb_candidates.children = []
                return

            # clamp index
            if i >= total:
                state["idx"] = total - 1
                name = current_name()

            row = df_all.loc[name] if name in df_all.index else pd.Series(dtype=object)
            ik   = str(row.get("INCHIKEY", "") or "").strip()
            cid  = str(row.get("PUBCHEM", "") or "").strip()
            pw   = str(row.get("PATHWAY", "") or "").strip()
            mass = str(row.get("MASS", "") or "").strip()

            cid_str = (
                f'<a href="https://pubchem.ncbi.nlm.nih.gov/compound/{cid}" target="_blank">{cid}</a>'
                if cid and cid != "nan" else "—"
            )
            ik_str = (
                f'<a href="https://www.ebi.ac.uk/chembl/compound_report_card/{ik}" target="_blank">{ik}</a>'
                if ik and ik != "nan" else "—"
            )

            curated = len(self.overrides)
            w_progress.value = (
                f"<b>{i+1} / {total}</b> unmatched &nbsp;·&nbsp; "
                f"<span style='color:green'>{curated} curated</span> &nbsp;·&nbsp; "
                f"{len(self.excluded)} excluded"
            )
            w_info.value = (
                "<table style='font-size:13px;border-collapse:collapse'>"
                f"<tr><td style='padding:2px 12px 2px 0'><b>Name</b></td><td>{name}</td></tr>"
                f"<tr><td><b>Pathway</b></td><td>{pw or '—'}</td></tr>"
                f"<tr><td><b>InChIKey</b></td><td>{ik_str}</td></tr>"
                f"<tr><td><b>PubChem CID</b></td><td>{cid_str}</td></tr>"
                f"<tr><td><b>Mass</b></td><td>{mass or '—'}</td></tr>"
                "</table>"
            )
            w_candidates.children = [
                widgets.Label("Click 'Search ChEBI' to fetch name-based candidates.")
            ]
            w_hmdb_candidates.children = [
                widgets.Label("Click 'Search VMH/HMDB' to fetch graph-backed HMDB candidates.")
            ]

        def advance():
            if state["idx"] < len(state["rows"]) - 1:
                state["idx"] += 1
            refresh()

        # ---- event handlers ----
        def on_search(_):
            name = current_name()
            if not name:
                return
            w_candidates.children = [widgets.Label("Searching…")]
            cands = self.search_chebi(name, 6)
            state["candidates"] = cands
            if not cands:
                w_candidates.children = [widgets.Label("No ChEBI candidates found.")]
                return
            btns = []
            for c in cands:
                conn = (
                    f"{c['reaction_count']} reactions in graph"
                    if c.get("in_graph") else "not in graph"
                )
                label = f"{c['chebi_id']}  —  {c['name']}  [{conn}]"
                b = widgets.Button(
                    description=label[:100],
                    button_style="success",
                    layout=widgets.Layout(width="750px"),
                )
                chebi_id = c["chebi_id"]
                def _make_assign_fn(cid_):
                    def fn(_):
                        self.assign(current_name(), cid_, notes=w_notes.value)
                        advance()
                    return fn
                b.on_click(_make_assign_fn(chebi_id))
                btns.append(b)
            w_candidates.children = btns

        def on_skip(_):
            advance()

        def on_search_hmdb(_):
            name = current_name()
            if not name:
                return
            w_hmdb_candidates.children = [widgets.Label("Searching…")]
            cands = self.search_hmdb(name, 6)
            if not cands:
                w_hmdb_candidates.children = [widgets.Label("No VMH/HMDB candidates found.")]
                return
            btns = []
            for c in cands:
                conn = (
                    f"{c['reaction_count']} reactions in graph"
                    if c.get("in_graph") else "not in graph"
                )
                chebi = f" · {c['chebi_id']}" if c.get("chebi_id") else ""
                label = f"{c['hmdb_id']}  —  {c['name']}{chebi}  [{conn}]"
                b = widgets.Button(
                    description=label[:100],
                    button_style="success",
                    layout=widgets.Layout(width="750px"),
                )
                hmdb_id = c["hmdb_id"]
                def _make_assign_hmdb_fn(hid_):
                    def fn(_):
                        self.assign(current_name(), hmdb_id=hid_, notes=w_notes.value)
                        advance()
                    return fn
                b.on_click(_make_assign_hmdb_fn(hmdb_id))
                btns.append(b)
            w_hmdb_candidates.children = btns

        def on_exclude(_):
            name = current_name()
            if name:
                self.exclude(name, notes=w_notes.value)
            advance()

        def on_assign(_):
            name = current_name()
            chebi = w_custom.value.strip()
            if not (name and chebi):
                with w_status:
                    print("Enter a CHEBI ID first.")
                return
            info = self.verify_connection(chebi)
            self.assign(name, chebi_id=chebi, notes=w_notes.value)
            with w_status:
                conn = f"{info['reaction_count']} reactions" if info.get("in_graph") else "not in graph"
                print(f"✓ Assigned {chebi} to '{name}'  ({conn})")
            w_custom.value = ""
            advance()

        def on_assign_hmdb(_):
            name = current_name()
            hmdb = w_custom_hmdb.value.strip()
            if not (name and hmdb):
                with w_status:
                    print("Enter an HMDB ID first.")
                return
            info = self.verify_hmdb_connection(hmdb)
            self.assign(name, hmdb_id=hmdb, notes=w_notes.value)
            with w_status:
                conn = f"{info['reaction_count']} reactions" if info.get("in_graph") else "not in graph"
                print(f"✓ Assigned {hmdb} to '{name}'  ({conn})")
            w_custom_hmdb.value = ""
            advance()

        def on_review_edges(_):
            name = current_name()
            if not name:
                return
            connections = self.reaction_connections(name)
            if not connections:
                w_edge_review.children = [
                    widgets.HTML("<i>No reaction connections found — assign a ChEBI ID first, "
                                 "or the matched ChEBI is not yet in the graph.</i>")
                ]
                return

            rows = []
            header = widgets.HTML(
                f"<b>{len(connections)} reaction connections</b> for <i>{name}</i> — "
                "click to exclude (red) or approve (green):"
            )
            rows.append(header)

            for rxn in connections:
                rid = rxn["reaction_id"]
                rname = rxn["name"][:60]
                role = rxn["role"] or "?"
                ec = ", ".join(rxn["ec_numbers"][:3]) or "—"
                genes = ", ".join(rxn["gene_symbols"][:4]) or "—"
                is_excl = rxn["excluded"]
                is_appr = rxn["approved"]

                style = "background:#fdd" if is_excl else ("background:#dfd" if is_appr else "")
                label = f"[{'✗' if is_excl else ('✓' if is_appr else '·')}]  {rname}  [{role}]  EC:{ec}  genes:{genes}"

                btn_excl = widgets.Button(
                    description="Exclude",
                    button_style="danger" if not is_excl else "",
                    layout=widgets.Layout(width="80px"),
                    disabled=is_excl,
                )
                btn_appr = widgets.Button(
                    description="Approve",
                    button_style="success" if not is_appr else "",
                    layout=widgets.Layout(width="80px"),
                    disabled=is_appr,
                )
                btn_reset = widgets.Button(
                    description="Reset",
                    layout=widgets.Layout(width="70px"),
                    disabled=not (is_excl or is_appr),
                )
                lbl = widgets.HTML(
                    f"<span style='font-size:12px;{style}'>{label}</span>",
                    layout=widgets.Layout(flex="1"),
                )

                def _make_excl(n_, r_):
                    def fn(_):
                        self.exclude_reaction(n_, r_)
                        on_review_edges(None)
                    return fn
                def _make_appr(n_, r_):
                    def fn(_):
                        self.approve_reaction(n_, r_)
                        on_review_edges(None)
                    return fn
                def _make_reset(n_, r_):
                    def fn(_):
                        self.restore_reaction(n_, r_)
                        on_review_edges(None)
                    return fn

                btn_excl.on_click(_make_excl(name, rid))
                btn_appr.on_click(_make_appr(name, rid))
                btn_reset.on_click(_make_reset(name, rid))

                rows.append(widgets.HBox(
                    [lbl, btn_appr, btn_excl, btn_reset],
                    layout=widgets.Layout(align_items="center", margin="1px 0"),
                ))

            w_edge_review.children = rows

        def on_prev(_):
            if state["idx"] > 0:
                state["idx"] -= 1
            refresh()

        def on_save(_):
            self.save()
            with w_status:
                print(f"✓ Saved {len(self.overrides)} overrides → {self.overrides_path}")

        btn_search.on_click(on_search)
        btn_search_hmdb.on_click(on_search_hmdb)
        btn_edges.on_click(on_review_edges)
        btn_skip.on_click(on_skip)
        btn_exclude.on_click(on_exclude)
        btn_assign.on_click(on_assign)
        btn_assign_hmdb.on_click(on_assign_hmdb)
        btn_prev.on_click(on_prev)
        btn_save.on_click(on_save)

        refresh()

        ui = widgets.VBox([
            w_progress,
            w_info,
            widgets.HBox([btn_prev, btn_search, btn_search_hmdb, btn_edges, btn_skip, btn_exclude, btn_save]),
            widgets.HBox([w_custom, btn_assign, w_notes]),
            widgets.HBox([w_custom_hmdb, btn_assign_hmdb]),
            w_candidates,
            w_hmdb_candidates,
            w_edge_review,
            w_status,
        ], layout=widgets.Layout(padding="10px"))

        display(ui)
        return ui
