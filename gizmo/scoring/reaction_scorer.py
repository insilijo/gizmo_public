"""
Transparent reaction evidence scorer.

Algorithm
---------
For each reaction R the scorer aggregates signed, weighted evidence
from three sources:

  1. Substrate metabolites
       depleted  (effect_size < 0) → reaction consuming substrate  → +
       accumulated (effect_size > 0) → reaction may be blocked      → −

  2. Product metabolites
       accumulated (effect_size > 0) → reaction producing product   → +
       depleted    (effect_size < 0) → reaction not producing        → −

  3. Catalysing genes / enzymes
       upregulated   → more enzyme activity                          → +
       downregulated → less enzyme activity                          → −

Reversibility handling
----------------------
For reactions marked ``reversible=True`` the direction of flux is
ambiguous from metabolomics alone.  In this case:

  • Score magnitude uses ``|contrib|`` — any metabolite shift is
    evidence of the reaction being active regardless of direction.
  • Direction is still computed from the signed metabolite pattern:
    substrate depletion + product accumulation → net forward (+);
    substrate accumulation + product depletion → net reverse (−).
  • The ``reversible`` flag is propagated to ``ReactionScore`` so
    callers can suppress or weight the direction signal appropriately.

Each contribution is weighted by:
  • confidence   — mapping quality [0, 1]
  • hub penalty  — 1 / sqrt(reaction_degree) to down-weight promiscuous metabolites
  • currency flag — currency metabolites contribute 0 (excluded)

Raw score = Σ(weighted contributions) / sqrt(max(n_evidence, 1))
            (√n normalisation prevents pure-count inflation)

Direction = mean signed direction contribution, independent of magnitude.

Usage::

    from gizmo.scoring.reaction_scorer import score_reactions

    scores = score_reactions(mg, ctx)
    scores.sort(key=lambda r: abs(r.score), reverse=True)
    for r in scores[:20]:
        print(r)
"""

from __future__ import annotations

import math
import logging
from dataclasses import dataclass, field
from typing import Optional

from gizmo.evidence.model import EvidenceRecord, SampleContext

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Result dataclass
# ---------------------------------------------------------------------------

@dataclass
class ReactionScore:
    """
    Evidence-weighted score for one reaction node.

    Attributes
    ----------
    reaction_id          : graph node ID  (e.g. "reactome:R-HSA-70171")
    score                : signed, √n-normalised aggregate score
                           positive → reaction likely more active
                           negative → reaction likely suppressed
    direction            : net signed direction [-1, +1], independent of magnitude
    evidence_count       : total number of evidence features contributing
    supporting_metabolites : list of {node_id, role, effect_size, contribution, confidence}
    supporting_genes       : list of {node_id, symbol, effect_size, contribution}
    contradictory_metabolites : features whose signal conflicts with the majority direction
    pathway_ids          : Reactome stIds the reaction belongs to
    ec_numbers           : EC numbers for the reaction (if annotated)
    confidence           : mean mapping confidence across all contributing features
    notes                : human-readable score explanation
    """

    reaction_id:               str
    score:                     float = 0.0
    direction:                 float = 0.0        # net direction [-1, 1]
    reversible:                bool  = False       # from ReactionNode.reversible
    evidence_count:            int   = 0
    supporting_metabolites:    list[dict] = field(default_factory=list)
    supporting_genes:          list[dict] = field(default_factory=list)
    contradictory_metabolites: list[dict] = field(default_factory=list)
    pathway_ids:               list[str]  = field(default_factory=list)
    ec_numbers:                list[str]  = field(default_factory=list)
    confidence:                float = 0.0
    notes:                     str   = ""

    def __repr__(self) -> str:
        direction_arrow = "▲" if self.direction > 0.1 else ("▼" if self.direction < -0.1 else "~")
        rev_tag = "⇌" if self.reversible else ""
        return (
            f"ReactionScore({self.reaction_id!r}, "
            f"score={self.score:+.3f}, {direction_arrow}{rev_tag}, "
            f"n={self.evidence_count})"
        )


# ---------------------------------------------------------------------------
# Scorer
# ---------------------------------------------------------------------------

def score_reactions(
    mg,
    ctx: SampleContext,
    min_evidence: int = 1,
    currency_weight: float = 0.0,
    hub_penalty: bool = True,
    gene_weight: float = 1.0,
    substrate_weight: float = 1.0,
    product_weight: float = 1.0,
) -> list[ReactionScore]:
    """
    Score all reactions in ``mg`` against multiomic evidence in ``ctx``.

    Parameters
    ----------
    mg               : GizmoGraph
    ctx              : SampleContext with mapped EvidenceRecords
    min_evidence     : skip reactions with fewer than this many evidence features
    currency_weight  : weight for currency metabolites (default 0 = exclude)
    hub_penalty      : apply 1/√degree down-weighting for high-degree metabolites
    gene_weight      : scaling factor for gene/enzyme contributions
    substrate_weight : scaling factor for substrate contributions
    product_weight   : scaling factor for product contributions

    Returns
    -------
    List of ReactionScore, one per reaction with at least min_evidence features.
    Reactions with no mapped evidence are omitted.
    """
    g = mg.graph

    # Pre-compute metabolite reaction-degrees for hub penalty
    met_degree: dict[str, int] = {}
    if hub_penalty:
        for nid, attrs in g.nodes(data=True):
            if attrs.get("node_type") == "metabolite":
                rxn_nbrs = sum(
                    1 for nbr in list(g.predecessors(nid)) + list(g.successors(nid))
                    if g.nodes[nbr].get("node_type") == "reaction"
                )
                met_degree[nid] = max(rxn_nbrs, 1)

    def _hub_w(nid: str) -> float:
        if not hub_penalty:
            return 1.0
        return 1.0 / math.sqrt(met_degree.get(nid, 1))

    results: list[ReactionScore] = []

    for rxn_id, rxn_attrs in g.nodes(data=True):
        if rxn_attrs.get("node_type") != "reaction":
            continue

        reversible: bool = bool(rxn_attrs.get("reversible", False))

        contributions:    list[float] = []   # magnitude scores for activity
        dir_contributions: list[float] = []  # signed direction signals
        conf_sum:          float       = 0.0
        sup_mets:          list[dict]  = []
        sup_genes:         list[dict]  = []
        contra_mets:       list[dict]  = []

        # ---- Substrate metabolites (predecessors with role=substrate) ----
        for met_id in g.predecessors(rxn_id):
            edata = g.edges[met_id, rxn_id]
            role  = edata.get("role", "")
            if role != "substrate":
                continue
            attrs = g.nodes[met_id]
            if attrs.get("node_type") != "metabolite":
                continue

            is_currency = bool(attrs.get("is_currency", False))
            w_currency  = currency_weight if is_currency else 1.0
            w_hub       = _hub_w(met_id)
            w_base      = substrate_weight * w_currency * w_hub

            for rec in ctx.by_node(met_id):
                if w_base == 0.0:
                    continue
                # Depleted substrate → positive signal (forward direction)
                signed_contrib = -rec.effect_size * rec.confidence * w_base
                if reversible:
                    # For reversible reactions: magnitude = activity level;
                    # direction signal: depletion = forward (+), accumulation = reverse (−)
                    contrib = abs(signed_contrib)
                    dir_contrib = signed_contrib
                else:
                    contrib = signed_contrib
                    dir_contrib = signed_contrib
                contributions.append(contrib)
                dir_contributions.append(dir_contrib)
                conf_sum += rec.confidence
                entry = {
                    "node_id":      met_id,
                    "role":         "substrate",
                    "effect_size":  rec.effect_size,
                    "contribution": round(contrib, 4),
                    "signed_contribution": round(signed_contrib, 4),
                    "confidence":   rec.confidence,
                    "is_currency":  is_currency,
                    "hub_weight":   round(w_hub, 3),
                }
                # For reversible reactions all metabolite shifts are supporting;
                # for irreversible, contradictory = opposing the expected direction
                if reversible or signed_contrib >= 0:
                    sup_mets.append(entry)
                else:
                    contra_mets.append(entry)

        # ---- Product metabolites (successors with role=product) ----
        for met_id in g.successors(rxn_id):
            edata = g.edges[rxn_id, met_id]
            role  = edata.get("role", "")
            if role != "product":
                continue
            attrs = g.nodes[met_id]
            if attrs.get("node_type") != "metabolite":
                continue

            is_currency = bool(attrs.get("is_currency", False))
            w_currency  = currency_weight if is_currency else 1.0
            w_hub       = _hub_w(met_id)
            w_base      = product_weight * w_currency * w_hub

            for rec in ctx.by_node(met_id):
                if w_base == 0.0:
                    continue
                # Accumulated product → positive signal (forward direction)
                signed_contrib = rec.effect_size * rec.confidence * w_base
                if reversible:
                    # For reversible reactions: magnitude = activity level;
                    # direction signal: accumulation = forward (+), depletion = reverse (−)
                    contrib = abs(signed_contrib)
                    dir_contrib = signed_contrib
                else:
                    contrib = signed_contrib
                    dir_contrib = signed_contrib
                contributions.append(contrib)
                dir_contributions.append(dir_contrib)
                conf_sum += rec.confidence
                entry = {
                    "node_id":      met_id,
                    "role":         "product",
                    "effect_size":  rec.effect_size,
                    "contribution": round(contrib, 4),
                    "signed_contribution": round(signed_contrib, 4),
                    "confidence":   rec.confidence,
                    "is_currency":  is_currency,
                    "hub_weight":   round(w_hub, 3),
                }
                if reversible or signed_contrib >= 0:
                    sup_mets.append(entry)
                else:
                    contra_mets.append(entry)

        # ---- Catalysing genes (predecessors with node_type=gene) ----
        for gene_id in g.predecessors(rxn_id):
            if g.nodes[gene_id].get("node_type") != "gene":
                continue
            for rec in ctx.by_node(gene_id):
                contrib = rec.effect_size * rec.confidence * gene_weight
                contributions.append(abs(contrib) if reversible else contrib)
                dir_contributions.append(contrib)
                conf_sum += rec.confidence
                sup_genes.append({
                    "node_id":      gene_id,
                    "symbol":       g.nodes[gene_id].get("symbol") or gene_id,
                    "effect_size":  rec.effect_size,
                    "contribution": round(contrib, 4),
                    "confidence":   rec.confidence,
                    "assay_type":   rec.assay_type,
                })

        n = len(contributions)
        if n < min_evidence:
            continue

        raw_score = sum(contributions) / math.sqrt(n)
        # Direction from signed signals (independent of reversibility magnitude treatment)
        direction = sum(1 if c > 0 else -1 for c in dir_contributions) / len(dir_contributions)
        mean_conf = conf_sum / n

        # Build a brief notes string
        n_s = len([e for e in sup_mets if e["role"] == "substrate"])
        n_p = len([e for e in sup_mets if e["role"] == "product"])
        n_c = len(contra_mets)
        n_g = len(sup_genes)
        rev_note = " [reversible — magnitude only]" if reversible else ""
        notes = (
            f"{n_s} substrate + {n_p} product signals, "
            f"{n_g} gene signals, {n_c} contradictory{rev_note}"
        )

        results.append(ReactionScore(
            reaction_id               = rxn_id,
            score                     = round(raw_score, 4),
            direction                 = round(direction, 3),
            reversible                = reversible,
            evidence_count            = n,
            supporting_metabolites    = sorted(sup_mets,   key=lambda x: -abs(x["contribution"])),
            supporting_genes          = sorted(sup_genes,  key=lambda x: -abs(x["contribution"])),
            contradictory_metabolites = sorted(contra_mets, key=lambda x: abs(x["contribution"])),
            pathway_ids               = list(rxn_attrs.get("pathways") or []),
            ec_numbers                = list(rxn_attrs.get("ec_numbers") or []),
            confidence                = round(mean_conf, 3),
            notes                     = notes,
        ))

    log.info(
        "Scored %d reactions from %d evidence records (%d nodes mapped).",
        len(results), len(ctx), len(ctx.mapped_nodes()),
    )
    return results
