"""
STRING protein–protein interaction (PPI) database source.

Data:    https://string-db.org   (CC BY 4.0)
API:     https://string-db.org/help/api/

STRING assigns each interaction a combined confidence score (0–1000).
Common thresholds:
  400  medium confidence   (default)
  700  high confidence
  900  highest confidence

The loader adds bidirectional PPI edges between *existing* GeneNode entries.
No new nodes are created — gene coverage in the graph depends on how the
graph was built (Reactome, Open Targets, Orphanet, etc.).

Usage::

    from gizmo.sources.stringdb import StringDBLoader

    loader = StringDBLoader(min_score=0.4)   # medium confidence
    n_added = loader.enrich(mg)              # adds edges in-place, returns edge count

    # Or restrict to high-confidence only:
    n_added = loader.enrich(mg, min_score=0.7)
"""

from __future__ import annotations

import logging
import time
from pathlib import Path
from typing import Optional

import requests

log = logging.getLogger(__name__)

_API_BASE = "https://string-db.org/api"
_SPECIES_HUMAN = 9606
_REQUEST_DELAY = 0.5   # seconds between multi-chunk requests (polite use)
_CHUNK_SIZE    = 1000  # STRING API practical limit per request


# ---------------------------------------------------------------------------
# REST client
# ---------------------------------------------------------------------------

class StringDBClient:
    """
    Thin wrapper around the STRING REST API v12.
    No authentication required.
    """

    def __init__(self) -> None:
        self._session = requests.Session()
        self._session.headers["User-Agent"] = "gizmo/1.0"

    def get_interactions(
        self,
        identifiers: list[str],
        species: int = _SPECIES_HUMAN,
        required_score: int = 400,
    ) -> list[dict]:
        """
        Return all within-set interactions among a list of gene symbols or Ensembl IDs.

        Parameters
        ----------
        identifiers   : gene symbols or ENSG / ENSP IDs
        species       : NCBI taxonomy (9606 = human)
        required_score: minimum combined score 0–1000

        Returns list of dicts with keys:
            stringId_A, stringId_B, preferredName_A, preferredName_B,
            score, nscore, fscore, pscore, ascore, escore, dscore, tscore
        """
        results: list[dict] = []
        chunks = [
            identifiers[i : i + _CHUNK_SIZE]
            for i in range(0, len(identifiers), _CHUNK_SIZE)
        ]
        for idx, chunk in enumerate(chunks):
            if idx > 0:
                time.sleep(_REQUEST_DELAY)
            resp = self._session.post(
                f"{_API_BASE}/json/network",
                data={
                    "identifiers":     "\r".join(chunk),
                    "species":         species,
                    "required_score":  required_score,
                    "caller_identity": "gizmo",
                },
                timeout=90,
            )
            resp.raise_for_status()
            results.extend(resp.json())
            log.debug("Chunk %d/%d → %d interactions so far", idx + 1, len(chunks), len(results))

        return results

    def resolve_ids(
        self,
        identifiers: list[str],
        species: int = _SPECIES_HUMAN,
    ) -> list[dict]:
        """
        Map gene symbols / Ensembl IDs to canonical STRING IDs.

        Returns list of dicts with:
            queryIndex, stringId, ncbiTaxonId, taxonName, preferredName, annotation
        """
        resp = self._session.post(
            f"{_API_BASE}/json/get_string_ids",
            data={
                "identifiers":     "\r".join(identifiers),
                "species":         species,
                "caller_identity": "gizmo",
            },
            timeout=60,
        )
        resp.raise_for_status()
        return resp.json()


# ---------------------------------------------------------------------------
# Graph enrichment loader
# ---------------------------------------------------------------------------

class StringDBLoader:
    """
    Enriches a GizmoGraph with PPI edges from STRING.

    Only adds edges between gene nodes that already exist in the graph.
    Gene nodes are matched by their `symbol` attribute (HGNC symbol).

    Parameters
    ----------
    min_score : minimum combined confidence score [0, 1]  (default 0.4)
    species   : NCBI taxonomy ID  (default 9606 = Homo sapiens)
    """

    def __init__(
        self,
        min_score: float = 0.4,
        species: int = _SPECIES_HUMAN,
    ) -> None:
        self.min_score = min_score
        self.species   = species
        self._client   = StringDBClient()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def enrich(self, mg, min_score: float | None = None) -> int:
        """
        Query STRING and add PPI edges between existing gene nodes.

        Parameters
        ----------
        mg        : GizmoGraph
        min_score : override instance default [0, 1]

        Returns
        -------
        int : number of unique undirected interactions added
              (each stored as two directed edges in the DiGraph)
        """
        threshold   = min_score if min_score is not None else self.min_score
        g           = mg.graph
        int_thresh  = int(threshold * 1000)

        # Collect gene symbols from gene nodes already in the graph
        symbol_to_nid: dict[str, str] = {}
        for nid, attrs in g.nodes(data=True):
            if attrs.get("node_type") == "gene":
                sym = attrs.get("symbol")
                if sym:
                    symbol_to_nid[sym] = nid

        if not symbol_to_nid:
            log.warning("No gene nodes with symbols found — nothing to enrich.")
            return 0

        symbols = list(symbol_to_nid)
        log.info(
            "Querying STRING for %d gene symbols (min_score=%.2f / %d) …",
            len(symbols), threshold, int_thresh,
        )

        try:
            raw_interactions = self._client.get_interactions(
                symbols, species=self.species, required_score=int_thresh,
            )
        except requests.HTTPError as exc:
            log.error("STRING HTTP error: %s", exc)
            return 0
        except requests.RequestException as exc:
            log.error("STRING request failed: %s", exc)
            return 0

        log.info("STRING returned %d raw interactions.", len(raw_interactions))

        added  = 0
        seen:   set[frozenset] = set()

        for row in raw_interactions:
            sym_a = row.get("preferredName_A", "")
            sym_b = row.get("preferredName_B", "")
            nid_a = symbol_to_nid.get(sym_a)
            nid_b = symbol_to_nid.get(sym_b)

            if not nid_a or not nid_b or nid_a == nid_b:
                continue

            pair = frozenset({nid_a, nid_b})
            if pair in seen:
                continue
            seen.add(pair)

            # Normalise STRING integer scores (0–1000) to floats (0–1)
            attrs = {
                "edge_type":      "protein_interaction",
                "combined_score": round(row.get("score",   0) / 1000, 4),
                "experimental":   round(row.get("escore",  0) / 1000, 4),
                "coexpression":   round(row.get("ascore",  0) / 1000, 4),
                "database":       round(row.get("dscore",  0) / 1000, 4),
                "textmining":     round(row.get("tscore",  0) / 1000, 4),
                "source_db":      "stringdb",
            }

            # Store as bidirectional in the DiGraph (PPI is undirected)
            g.add_edge(nid_a, nid_b, **attrs)
            g.add_edge(nid_b, nid_a, **attrs)
            added += 1

        log.info(
            "Added %d PPI interactions (%d directed edges). "
            "%d pairs skipped (gene not in graph).",
            added, added * 2,
            len(raw_interactions) - added - len([r for r in raw_interactions
                                                  if r.get("preferredName_A") == r.get("preferredName_B")]),
        )
        return added

    # ------------------------------------------------------------------
    # Convenience: score distribution summary
    # ------------------------------------------------------------------

    def summarise(self, mg) -> None:
        """Print a quick summary of PPI edges already in the graph."""
        g    = mg.graph
        ppi  = [(u, v, d) for u, v, d in g.edges(data=True)
                if d.get("edge_type") == "protein_interaction"]
        if not ppi:
            print("No PPI edges in graph.")
            return

        n_unique = len(ppi) // 2
        scores   = [d["combined_score"] for _, _, d in ppi]
        print(f"PPI edges (STRING):  {n_unique} interactions ({len(ppi)} directed)")
        print(f"  combined_score — min {min(scores):.3f}  mean {sum(scores)/len(scores):.3f}  max {max(scores):.3f}")
        genes_with_ppi = {u for u, _, d in ppi if d.get("edge_type") == "protein_interaction"}
        print(f"  Gene nodes with ≥1 PPI partner: {len(genes_with_ppi)}")
