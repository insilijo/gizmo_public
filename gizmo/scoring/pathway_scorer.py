"""
Pathway-level summarisation from per-reaction scores.

Groups ReactionScores by Reactome pathway stId and computes aggregate
statistics: mean/median score, directional consistency, leading-edge
reactions, and total evidence count.

Usage::

    from gizmo.scoring.reaction_scorer import score_reactions
    from gizmo.scoring.pathway_scorer import summarise_pathways

    rxn_scores = score_reactions(mg, ctx)
    pw_summaries = summarise_pathways(rxn_scores)
    pw_summaries.sort(key=lambda p: abs(p.mean_score), reverse=True)
    for p in pw_summaries[:10]:
        print(p)
"""

from __future__ import annotations

import math
import statistics
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class PathwaySummary:
    """
    Aggregated evidence score for one Reactome pathway.

    Attributes
    ----------
    pathway_id         : Reactome stId
    pathway_name       : display name (populated if pathway_names dict provided)
    n_reactions        : number of scored reactions in this pathway
    n_reactions_total  : total reactions in this pathway in the graph (if mg provided)
    mean_score         : mean reaction score (signed)
    median_score       : median reaction score (signed)
    std_score          : standard deviation of reaction scores
    direction          : mean direction across reactions [-1, 1]
    consistency        : fraction of reactions with the same sign as mean_score [0, 1]
    evidence_count     : total evidence features across all reactions
    top_reactions      : top 5 reactions by |score| — list of ReactionScore
    """

    pathway_id:        str
    pathway_name:      str   = ""
    n_reactions:       int   = 0
    n_reactions_total: int   = 0
    mean_score:        float = 0.0
    median_score:      float = 0.0
    std_score:         float = 0.0
    direction:         float = 0.0
    consistency:       float = 0.0
    evidence_count:    int   = 0
    top_reactions:     list  = field(default_factory=list)

    def __repr__(self) -> str:
        name  = self.pathway_name or self.pathway_id
        arrow = "▲" if self.direction > 0.1 else ("▼" if self.direction < -0.1 else "~")
        return (
            f"PathwaySummary({name!r}, "
            f"mean={self.mean_score:+.3f}, {arrow}, "
            f"n={self.n_reactions}, consistency={self.consistency:.0%})"
        )


def summarise_pathways(
    reaction_scores: list,           # list[ReactionScore]
    pathway_names: Optional[dict[str, str]] = None,
    mg=None,
    top_n_reactions: int = 5,
) -> list[PathwaySummary]:
    """
    Aggregate per-reaction scores into pathway-level summaries.

    Parameters
    ----------
    reaction_scores : output of score_reactions()
    pathway_names   : optional dict {stId: displayName} for readable labels
                      (can be built with build.pipeline or app.dash_app helpers)
    mg              : if provided, total reaction counts per pathway are included
    top_n_reactions : how many top reactions to include in PathwaySummary.top_reactions

    Returns
    -------
    List of PathwaySummary, one per pathway, sorted by |mean_score| descending.
    """
    names = pathway_names or {}

    # Count total reactions per pathway in the graph
    total_per_pw: dict[str, int] = {}
    if mg is not None:
        g = mg.graph
        for nid, attrs in g.nodes(data=True):
            if attrs.get("node_type") != "reaction":
                continue
            for pw in (attrs.get("pathways") or []):
                total_per_pw[pw] = total_per_pw.get(pw, 0) + 1

    # Group reaction scores by pathway
    pw_groups: dict[str, list] = {}
    for rs in reaction_scores:
        for pw in rs.pathway_ids:
            pw_groups.setdefault(pw, []).append(rs)

    summaries: list[PathwaySummary] = []
    for pw_id, group in pw_groups.items():
        scores = [rs.score for rs in group]
        dirs   = [rs.direction for rs in group]
        n      = len(scores)
        mean_s = statistics.mean(scores)
        med_s  = statistics.median(scores)
        std_s  = statistics.stdev(scores) if n > 1 else 0.0
        mean_d = statistics.mean(dirs)
        # Consistency: fraction of reactions with same sign as mean_score
        dominant_sign = 1 if mean_s >= 0 else -1
        consistent    = sum(1 for s in scores if (s >= 0) == (dominant_sign >= 0))
        consistency   = consistent / n

        top = sorted(group, key=lambda r: abs(r.score), reverse=True)[:top_n_reactions]

        summaries.append(PathwaySummary(
            pathway_id        = pw_id,
            pathway_name      = names.get(pw_id, ""),
            n_reactions       = n,
            n_reactions_total = total_per_pw.get(pw_id, 0),
            mean_score        = round(mean_s, 4),
            median_score      = round(med_s, 4),
            std_score         = round(std_s, 4),
            direction         = round(mean_d, 3),
            consistency       = round(consistency, 3),
            evidence_count    = sum(rs.evidence_count for rs in group),
            top_reactions     = top,
        ))

    summaries.sort(key=lambda p: abs(p.mean_score), reverse=True)
    return summaries


def print_pathway_report(
    summaries: list[PathwaySummary],
    top_n: int = 20,
    show_reactions: int = 3,
) -> None:
    """Print a human-readable ranked pathway report."""
    print(f"\n{'Rank':<5} {'Score':>8}  {'Dir':>4}  {'Con%':>5}  {'Rxns':>5}  Pathway")
    print("-" * 75)
    for i, p in enumerate(summaries[:top_n], 1):
        arrow = "▲" if p.direction > 0.1 else ("▼" if p.direction < -0.1 else "~")
        label = (p.pathway_name or p.pathway_id)[:45]
        print(f"{i:<5} {p.mean_score:>8.3f}  {arrow:>4}  {p.consistency:>4.0%}  "
              f"{p.n_reactions:>5}  {label}")
        if show_reactions:
            for rs in p.top_reactions[:show_reactions]:
                arrow2 = "▲" if rs.direction > 0.1 else "▼"
                print(f"       {'':8}  {arrow2:>4}  {'':>5}  {'':>5}    "
                      f"  {rs.reaction_id}  ({rs.score:+.3f})")
