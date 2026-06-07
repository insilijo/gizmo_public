"""
Enzyme / gene target druggability scoring.

Queries ChEMBL for known drugs and tool compounds targeting each gene
node in the graph, and optionally enriches with Open Targets tractability
annotations.

Usage::

    from gizmo.actionability.druggability import score_druggability

    scores = score_druggability(mg, cache_dir="data/raw/chembl")
    scores.sort(key=lambda d: -d.score)
    for d in scores[:10]:
        print(d)
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

from gizmo.actionability.model import DruggabilityScore
from gizmo.sources.chembl import ChEMBLClient

log = logging.getLogger(__name__)

# Max clinical phase → base druggability score
_PHASE_SCORE = {4: 1.0, 3: 0.75, 2: 0.55, 1: 0.30, 0: 0.10}


def score_druggability(
    mg,
    gene_ids: Optional[list[str]] = None,
    cache_dir: Optional[str | Path] = None,
    min_phase: int = 0,
    enrich_open_targets: bool = False,
) -> list[DruggabilityScore]:
    """
    Score druggability for gene nodes in ``mg``.

    Parameters
    ----------
    mg                   : GizmoGraph
    gene_ids             : restrict to these node IDs (default: all gene nodes)
    cache_dir            : directory for ChEMBL response caching
    min_phase            : minimum ChEMBL max_phase to include (0 = all)
    enrich_open_targets  : also fetch tractability buckets from Open Targets

    Returns
    -------
    List of DruggabilityScore, one per gene with a resolved ChEMBL target.
    Genes with no ChEMBL hits are omitted.
    """
    g       = mg.graph
    client  = ChEMBLClient(cache_dir=cache_dir)

    targets = gene_ids or mg.gene_nodes()
    results: list[DruggabilityScore] = []

    for gene_id in targets:
        attrs  = g.nodes.get(gene_id, {})
        symbol = attrs.get("symbol") or ""
        if not symbol:
            continue

        log.debug("ChEMBL lookup: %s", symbol)
        drugs = client.drugs_for_gene(symbol, min_phase=min_phase)
        if not drugs:
            continue

        approved  = [d for d in drugs if d.get("max_phase", 0) == 4]
        clinical  = [d for d in drugs if d.get("max_phase", 0) >= 2]
        max_phase = max((d.get("max_phase", 0) for d in drugs), default=0)
        best      = max(drugs, key=lambda d: d.get("max_phase", 0))

        # Unique action types
        action_types = sorted({
            d.get("action_type", "")
            for d in drugs if d.get("action_type")
        })

        score = _PHASE_SCORE.get(max_phase, 0.05)

        ds = DruggabilityScore(
            gene_id          = gene_id,
            symbol           = symbol,
            target_chembl_id = best.get("target_chembl_id"),
            n_approved_drugs = len(approved),
            n_clinical_drugs = len(clinical),
            max_phase        = max_phase,
            best_drug_name   = best.get("molecule_name"),
            best_action_type = best.get("action_type"),
            score            = round(score, 3),
            drugs            = drugs,
        )

        # Optional Open Targets tractability enrichment
        if enrich_open_targets:
            _enrich_ot_tractability(ds, attrs.get("ensembl_id") or "")

        results.append(ds)

    results.sort(key=lambda d: -d.score)
    log.info(
        "Druggability scored %d / %d gene nodes with ChEMBL data.",
        len(results), len(targets),
    )
    return results


def _enrich_ot_tractability(ds: DruggabilityScore, ensembl_id: str) -> None:
    """Fetch Open Targets tractability and write into ds in-place."""
    if not ensembl_id:
        return
    try:
        from gizmo.sources.open_targets import OpenTargetsClient
        client = OpenTargetsClient()
        result = client._query(
            """
            query Tractability($id: String!) {
              target(ensemblId: $id) {
                tractability { label modality value }
              }
            }
            """,
            {"id": ensembl_id},
        )
        rows = (result.get("data", {})
                      .get("target", {})
                      .get("tractability") or [])
        for row in rows:
            if row.get("modality") == "SM" and row.get("value"):
                ds.tractability_sm = row.get("label")
            elif row.get("modality") == "AB" and row.get("value"):
                ds.tractability_ab = row.get("label")
    except Exception as exc:
        log.debug("OT tractability failed for %s: %s", ensembl_id, exc)
