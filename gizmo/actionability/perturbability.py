"""
Reaction perturbability scoring.

Combines enzyme druggability with reaction-level structural features
(cofactor dependence, metabolic redundancy) to estimate how actionable
each reaction is with known chemistry.

Usage::

    from gizmo.actionability.perturbability import score_perturbability
    from gizmo.actionability.druggability import score_druggability

    drug_scores = score_druggability(mg, cache_dir="data/raw/chembl")
    perturb     = score_perturbability(mg, drug_scores)

    # Combine with reaction evidence:
    from gizmo.scoring import score_reactions
    rxn_scores = score_reactions(mg, ctx)
    actions    = combine_actionability(perturb, rxn_scores)
    actions.sort(key=lambda a: -a.priority_score)
"""

from __future__ import annotations

import logging
from typing import Optional

from gizmo.actionability.model import (
    ActionabilityScore,
    DruggabilityScore,
    PerturbabilityScore,
)

log = logging.getLogger(__name__)


def score_perturbability(
    mg,
    drug_scores: list[DruggabilityScore],
    redundancy_penalty: float = 0.15,
) -> list[PerturbabilityScore]:
    """
    Score reaction perturbability from enzyme druggability.

    A reaction scores higher when:
      - at least one catalysing enzyme has approved / clinical drugs
      - there are few alternative routes (low metabolic redundancy)
      - the reaction has a cofactor requirement (additional lever point)

    Parameters
    ----------
    mg                : GizmoGraph
    drug_scores       : output of score_druggability()
    redundancy_penalty: score discount per alternative reaction route

    Returns
    -------
    List of PerturbabilityScore, one per reaction that has at least one
    catalysing enzyme with known chemistry.  Reactions with no druggable
    enzymes are omitted.
    """
    g = mg.graph

    # Index druggability by gene_id
    drug_index: dict[str, DruggabilityScore] = {d.gene_id: d for d in drug_scores}

    # Count alternative routes: for each substrate, how many reactions consume it?
    substrate_reaction_count: dict[str, int] = {}
    for nid, attrs in g.nodes(data=True):
        if attrs.get("node_type") != "reaction":
            continue
        for met_id in g.predecessors(nid):
            if (g.nodes[met_id].get("node_type") == "metabolite"
                    and g.edges[met_id, nid].get("role") == "substrate"):
                substrate_reaction_count[met_id] = (
                    substrate_reaction_count.get(met_id, 0) + 1
                )

    results: list[PerturbabilityScore] = []

    for rxn_id, rxn_attrs in g.nodes(data=True):
        if rxn_attrs.get("node_type") != "reaction":
            continue

        # Gene nodes catalysing this reaction
        gene_ids = [
            n for n in g.predecessors(rxn_id)
            if g.nodes[n].get("node_type") == "gene"
        ]

        # Druggable enzymes
        druggable = [drug_index[gid] for gid in gene_ids if gid in drug_index]
        if not druggable:
            continue

        # Best drug across all enzymes
        best_enz   = max(druggable, key=lambda d: d.score)
        max_phase  = max(d.max_phase for d in druggable)
        all_actions = sorted({
            at
            for d in druggable
            for at in (d.best_action_type,) if at
        })

        # Cofactor dependence: reaction has modifier metabolites?
        has_modifier = any(
            g.edges[n, rxn_id].get("role") == "modifier"
            for n in g.predecessors(rxn_id)
            if g.nodes[n].get("node_type") == "metabolite"
        )

        # Metabolic redundancy: mean number of reactions sharing each substrate
        substrates = [
            n for n in g.predecessors(rxn_id)
            if (g.nodes[n].get("node_type") == "metabolite"
                and not g.nodes[n].get("is_currency", False)
                and g.edges[n, rxn_id].get("role") == "substrate")
        ]
        if substrates:
            mean_alt = sum(
                substrate_reaction_count.get(s, 1) - 1
                for s in substrates
            ) / len(substrates)
        else:
            mean_alt = 0.0

        # Perturbability score
        base_score    = best_enz.score
        redundancy    = min(mean_alt * redundancy_penalty, 0.5)
        cofactor_bonus= 0.05 if has_modifier else 0.0
        score         = max(0.0, min(1.0, base_score - redundancy + cofactor_bonus))

        notes = (
            f"Best compound: {best_enz.best_drug_name or '—'} "
            f"(phase {max_phase}, {best_enz.symbol})"
        )
        if mean_alt > 0:
            notes += f"; {mean_alt:.1f} alternative routes"

        results.append(PerturbabilityScore(
            reaction_id          = rxn_id,
            n_druggable_enzymes  = len(druggable),
            max_phase            = max_phase,
            best_drug_name       = best_enz.best_drug_name,
            best_gene_symbol     = best_enz.symbol,
            action_types         = all_actions,
            cofactor_dependent   = has_modifier,
            n_alternative_routes = round(mean_alt),
            score                = round(score, 3),
            notes                = notes,
        ))

    results.sort(key=lambda p: -p.score)
    log.info("Perturbability scored %d reactions.", len(results))
    return results


def combine_actionability(
    perturbability_scores: list[PerturbabilityScore],
    reaction_scores: Optional[list] = None,       # list[ReactionScore]
    drug_scores: Optional[list[DruggabilityScore]] = None,
) -> list[ActionabilityScore]:
    """
    Combine perturbability with reaction evidence scores.

    Parameters
    ----------
    perturbability_scores : output of score_perturbability()
    reaction_scores       : output of score_reactions() — optional
    drug_scores           : output of score_druggability() — optional, for enzyme detail

    Returns
    -------
    List of ActionabilityScore sorted by priority_score descending.
    """
    rxn_score_map: dict[str, float] = {}
    if reaction_scores:
        rxn_score_map = {rs.reaction_id: rs.score for rs in reaction_scores}

    drug_index: dict[str, DruggabilityScore] = {}
    if drug_scores:
        drug_index = {d.gene_id: d for d in drug_scores}

    actions: list[ActionabilityScore] = []
    for ps in perturbability_scores:
        rxn_score   = rxn_score_map.get(ps.reaction_id, 0.0)
        priority    = ps.score * (1.0 + abs(rxn_score))

        explanation = ps.notes
        if rxn_score != 0.0:
            dir_word = "elevated" if rxn_score > 0 else "suppressed"
            explanation += f"; reaction evidence score {rxn_score:+.3f} ({dir_word})"

        actions.append(ActionabilityScore(
            reaction_id    = ps.reaction_id,
            perturbability = ps,
            combined_score = round(ps.score, 3),
            reaction_score = round(rxn_score, 4),
            priority_score = round(priority, 4),
            explanation    = explanation,
        ))

    actions.sort(key=lambda a: -a.priority_score)
    return actions


def print_actionability_report(
    actions: list[ActionabilityScore],
    top_n: int = 20,
) -> None:
    """Print a ranked actionability report."""
    print(f"\n{'Rank':<5} {'Priority':>9}  {'Perturb':>8}  {'Phase':>5}  "
          f"{'RxnScore':>9}  Reaction")
    print("─" * 85)
    for i, a in enumerate(actions[:top_n], 1):
        ps = a.perturbability
        print(f"{i:<5} {a.priority_score:>9.4f}  {ps.score:>8.3f}  "
              f"{ps.max_phase:>5}  {a.reaction_score:>+9.3f}  {a.reaction_id}")
        if ps.best_drug_name:
            print(f"      {'':>9}  {'':>8}  {'':>5}  {'':>9}  "
                  f"    ↳ {ps.best_drug_name} ({ps.best_gene_symbol})")
