"""
Causal chain and mechanism ranking (Phase 5).

Builds and scores interpretable causal paths through the GIZMO graph.
Supported path shapes:

    disease → gene → reaction → metabolite
    gene → reaction → metabolite
    metabolite → reaction → metabolite   (propagation)

Each path becomes a ``PathHypothesis`` with a composite confidence score
based on:

  - reaction evidence score (from Phase 4 scorer)
  - edge association confidence (disease-gene score, gene-reaction evidence)
  - path length decay  (longer paths are penalised)
  - hub penalty        (intermediate high-degree nodes are down-weighted)
  - terminal metabolite evidence (bonus when the endpoint has omics support)

Usage::

    from gizmo.scoring.reaction_scorer import score_reactions
    from gizmo.scoring.chain_ranker import rank_chains, print_chain_report

    rxn_scores = score_reactions(mg, ctx)
    chains = rank_chains(mg, rxn_scores, ctx=ctx)
    print_chain_report(chains[:20])
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional

from gizmo.evidence.model import SampleContext


# ---------------------------------------------------------------------------
# PathHypothesis
# ---------------------------------------------------------------------------

@dataclass
class PathHypothesis:
    """
    A ranked causal path through the GIZMO graph.

    Attributes
    ----------
    nodes           : node IDs in causal order
    node_types      : type label for each node
    edge_types      : type label for each edge
    score           : composite ranked score (higher = more confident)
    reaction_id     : key reaction node in this path
    reaction_score  : signed evidence score of that reaction
    disease_id      : originating disease node (if applicable)
    evidence_support: sum of |effect_size × confidence| at terminal metabolite
    confidence      : mean edge confidence along path (0–1)
    explanation     : human-readable one-liner
    """

    nodes:            list[str]
    node_types:       list[str]
    edge_types:       list[str]
    score:            float
    reaction_id:      str
    reaction_score:   float
    disease_id:       Optional[str] = None
    evidence_support: float = 0.0
    confidence:       float = 1.0
    explanation:      str   = ""

    def __repr__(self) -> str:
        arrow = " → ".join(
            f"{ntype}({nid.split(':')[-1][:20]})"
            for nid, ntype in zip(self.nodes, self.node_types)
        )
        return f"PathHypothesis(score={self.score:+.3f}, {arrow})"


# ---------------------------------------------------------------------------
# Main ranker
# ---------------------------------------------------------------------------

def rank_chains(
    mg,
    reaction_scores: list,               # list[ReactionScore]
    ctx: Optional[SampleContext] = None,
    max_chains: int = 1000,
    length_decay: float = 0.80,
    min_reaction_score: float = 0.0,
    include_metabolite_propagation: bool = False,
) -> list[PathHypothesis]:
    """
    Build and rank causal chain hypotheses from per-reaction evidence scores.

    Parameters
    ----------
    mg                             : GizmoGraph
    reaction_scores                : output of score_reactions()
    ctx                            : SampleContext for terminal evidence lookup
    max_chains                     : cap on returned hypotheses
    length_decay                   : per-hop score decay factor [0, 1]
    min_reaction_score             : skip reactions below this |score|
    include_metabolite_propagation : also build met→reaction→met paths

    Returns
    -------
    List of PathHypothesis sorted by score descending.
    """
    g = mg.graph

    rxn_map = {rs.reaction_id: rs for rs in reaction_scores
               if abs(rs.score) >= min_reaction_score}

    def _hub_w(nid: str) -> float:
        deg = g.in_degree(nid) + g.out_degree(nid)
        return 1.0 / math.sqrt(max(deg, 1))

    def _terminal_evidence(node_id: str) -> float:
        if ctx is None:
            return 0.0
        return sum(abs(r.effect_size * r.confidence) for r in ctx.by_node(node_id))

    chains: list[PathHypothesis] = []

    for rxn_id, rs in rxn_map.items():
        if rxn_id not in g:
            continue

        rxn_attrs = g.nodes[rxn_id]
        abs_score = abs(rs.score)

        # Non-currency product metabolites
        products = [
            n for n in g.successors(rxn_id)
            if (g.nodes[n].get("node_type") == "metabolite"
                and not g.nodes[n].get("is_currency", False)
                and g.edges[rxn_id, n].get("role") == "product")
        ][:6]   # cap fan-out

        # Non-currency substrate metabolites (for propagation)
        substrates = [
            n for n in g.predecessors(rxn_id)
            if (g.nodes[n].get("node_type") == "metabolite"
                and not g.nodes[n].get("is_currency", False)
                and g.edges[n, rxn_id].get("role") == "substrate")
        ][:6]

        # Catalysing genes
        genes = [
            n for n in g.predecessors(rxn_id)
            if g.nodes[n].get("node_type") == "gene"
        ]

        ec = rxn_attrs.get("ec_numbers") or []
        rxn_label = rxn_id.split(":")[-1]

        # ------------------------------------------------------------------
        # Shape 1: gene → reaction → metabolite
        # ------------------------------------------------------------------
        for gene_id in genes:
            gene_sym = g.nodes[gene_id].get("symbol") or gene_id
            w_gene   = _hub_w(gene_id)
            base     = abs_score * length_decay * w_gene

            for prod_id in products:
                ev = _terminal_evidence(prod_id)
                w_prod = _hub_w(prod_id)
                score  = round(base * length_decay * w_prod + ev * 0.05, 4)

                chains.append(PathHypothesis(
                    nodes           = [gene_id, rxn_id, prod_id],
                    node_types      = ["gene", "reaction", "metabolite"],
                    edge_types      = ["gene_reaction", "product"],
                    score           = score,
                    reaction_id     = rxn_id,
                    reaction_score  = rs.score,
                    evidence_support= round(ev, 3),
                    confidence      = round(w_gene, 3),
                    explanation     = (
                        f"Gene {gene_sym} catalyses {rxn_label}"
                        + (f" [EC {ec[0]}]" if ec else "")
                        + f" → {prod_id.split(':')[-1]}"
                    ),
                ))

            # ------------------------------------------------------------------
            # Shape 2: disease → gene → reaction → metabolite
            # ------------------------------------------------------------------
            for dis_id in g.predecessors(gene_id):
                if g.nodes[dis_id].get("node_type") != "disease":
                    continue
                dis_edge  = g.edges.get((dis_id, gene_id), {})
                dis_conf  = float(dis_edge.get("score") or 0.4)
                w_dis     = _hub_w(dis_id)
                dis_base  = abs_score * (length_decay ** 2) * dis_conf * w_gene * w_dis

                for prod_id in products[:3]:
                    ev = _terminal_evidence(prod_id)
                    w_prod = _hub_w(prod_id)
                    score  = round(dis_base * length_decay * w_prod + ev * 0.05, 4)

                    dis_name = (
                        g.nodes[dis_id].get("name")
                        or g.nodes[dis_id].get("label")
                        or dis_id.split(":")[-1]
                    )
                    chains.append(PathHypothesis(
                        nodes           = [dis_id, gene_id, rxn_id, prod_id],
                        node_types      = ["disease", "gene", "reaction", "metabolite"],
                        edge_types      = ["disease_gene", "gene_reaction", "product"],
                        score           = score,
                        reaction_id     = rxn_id,
                        reaction_score  = rs.score,
                        disease_id      = dis_id,
                        evidence_support= round(ev, 3),
                        confidence      = round(dis_conf, 3),
                        explanation     = (
                            f"{dis_name} → {gene_sym} → {rxn_label}"
                            + f" → {prod_id.split(':')[-1]}"
                        ),
                    ))

        # ------------------------------------------------------------------
        # Shape 3: metabolite → reaction → metabolite  (propagation)
        # ------------------------------------------------------------------
        if include_metabolite_propagation:
            for sub_id in substrates:
                ev_sub = _terminal_evidence(sub_id)
                if ev_sub == 0.0:
                    continue   # only trace from metabolites with evidence
                w_sub  = _hub_w(sub_id)
                base_p = abs_score * length_decay * w_sub * (1.0 + ev_sub * 0.1)

                for prod_id in products:
                    ev_prod = _terminal_evidence(prod_id)
                    w_prod  = _hub_w(prod_id)
                    score   = round(base_p * length_decay * w_prod + ev_prod * 0.05, 4)

                    chains.append(PathHypothesis(
                        nodes           = [sub_id, rxn_id, prod_id],
                        node_types      = ["metabolite", "reaction", "metabolite"],
                        edge_types      = ["substrate", "product"],
                        score           = score,
                        reaction_id     = rxn_id,
                        reaction_score  = rs.score,
                        evidence_support= round(ev_sub + ev_prod, 3),
                        confidence      = round(w_sub, 3),
                        explanation     = (
                            f"{sub_id.split(':')[-1]} → {rxn_label}"
                            f" → {prod_id.split(':')[-1]}"
                        ),
                    ))

    # Deduplicate (same node sequence) and sort
    seen:   set[tuple] = set()
    unique: list[PathHypothesis] = []
    for ch in sorted(chains, key=lambda x: -x.score):
        key = tuple(ch.nodes)
        if key not in seen:
            seen.add(key)
            unique.append(ch)
        if len(unique) >= max_chains:
            break

    return unique


# ---------------------------------------------------------------------------
# Report helper
# ---------------------------------------------------------------------------

def print_chain_report(
    chains: list[PathHypothesis],
    top_n: int = 20,
    group_by_disease: bool = False,
) -> None:
    """Print a ranked causal chain report."""
    if group_by_disease:
        from collections import defaultdict
        by_dis: dict = defaultdict(list)
        for ch in chains:
            by_dis[ch.disease_id or "(no disease)"].append(ch)
        for dis, group in sorted(by_dis.items(), key=lambda x: -max(c.score for c in x[1])):
            print(f"\n── {dis} ──")
            for ch in sorted(group, key=lambda c: -c.score)[:5]:
                print(f"  {ch.score:+.3f}  {ch.explanation}")
        return

    print(f"\n{'Rank':<5} {'Score':>8}  {'RxnScore':>9}  {'Evid':>6}  Path")
    print("─" * 80)
    for i, ch in enumerate(chains[:top_n], 1):
        print(f"{i:<5} {ch.score:>8.4f}  {ch.reaction_score:>+9.3f}  "
              f"{ch.evidence_support:>6.2f}  {ch.explanation}")
