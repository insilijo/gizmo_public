"""
Enzyme-availability gate for reaction-level scores.

GIZMO BP propagates substrate-side metabolite evidence through the
reaction graph to estimate which reactions are perturbed. Transcriptomic
evidence (gene expression of catalyzing enzymes) is asymmetric and
direction-agnostic with respect to reaction state:

  - Gene UP / NORMAL  → enzyme is available; reaction state is determined
                        by substrate evidence.
  - Gene DOWN         → enzyme is capacity-limited; the reaction is less
                        likely to be running regardless of what the
                        substrate side suggests.

Pushing reaction direction with gene evidence (theta!=0 catalysis edges
in BP) creates sign conflicts with substrate evidence — IDH1 is DOWN in
IDH1-mut tumors but its mutant reaction (R-HSA-880053) is UP via
2-HG production. Coupling the gene's DOWN to the reaction's UP through a
+1 sign edge cancels the substrate-side signal.

Solution: keep BP gene-edge-agnostic (theta=0), but apply a post-BP
multiplicative gate on reaction perturbation magnitudes:

    final_score = bp_perturbation * enzyme_gate(reaction)

where ``enzyme_gate ∈ [0, 1]`` is < 1 only when the catalyzing gene is
strongly DOWN. Up-regulated or normal genes pass through unchanged.
"""
from __future__ import annotations

import logging
from typing import Iterable, Optional

log = logging.getLogger(__name__)


def compute_enzyme_gate(
    mg,
    gene_evidence: dict[str, float],
    *,
    down_threshold: float = -0.5,
    floor: float = 0.3,
    combine: str = "min",
) -> dict[str, float]:
    """
    Compute a per-reaction enzyme-availability gate from gene evidence.

    For each reaction R catalyzed by gene(s) {g1, g2, …}:
        per-gene gate(g) = 1.0                              if e_g >= 0
                         = 1.0 + (e_g - down_threshold) /
                                  abs(down_threshold)        if down_threshold < e_g < 0
                         = floor                             if e_g <= down_threshold

    Multiple catalyzing genes are combined by ``combine``:
        - ``"min"``   — most-restrictive (any one DOWN dampens the reaction)
        - ``"mean"``  — geometric mean (compensatory enzyme pools)
        - ``"max"``   — least-restrictive (any one available is sufficient)

    Parameters
    ----------
    mg
        GizmoGraph with gene→reaction (``catalysis``) edges.
    gene_evidence
        ``{gene_node_id: signed_effect_size}`` from the trans assay
        (e.g. each gene's log2 fold change).
    down_threshold
        Effect size at which gene gate hits ``floor`` (default -0.5,
        i.e. ~30% lower expression saturates the gate).
    floor
        Minimum gate value for a strongly-down gene (default 0.3 — even
        with no enzyme, leave 30% of the perturbation signal in case
        the gate is wrong about which gene catalyzes the reaction).
    combine
        How to combine multiple catalyzing genes per reaction.

    Returns
    -------
    {reaction_id: gate_value} for every reaction that has at least one
    catalyzing gene with evidence. Reactions without gene evidence are
    not in the dict (caller should treat them as gate=1.0).
    """
    if not gene_evidence:
        return {}
    if down_threshold >= 0:
        raise ValueError("down_threshold must be negative")
    g = mg.graph

    # Gather per-reaction list of (gene_id, effect_size) for genes with evidence
    rxn_genes: dict[str, list[float]] = {}
    for gene_id, effect in gene_evidence.items():
        if not g.has_node(gene_id):
            continue
        for _, rxn_id, ed in g.out_edges(gene_id, data=True):
            if g.nodes[rxn_id].get("node_type") != "reaction":
                continue
            role = (ed.get("edge_type") or ed.get("role") or "").lower()
            # Only count catalysis / gene_reaction / gene_associated edges
            if role not in ("catalysis", "gene_reaction", "gene_associated"):
                continue
            rxn_genes.setdefault(rxn_id, []).append(float(effect))

    # Convert per-gene effect → per-gene gate ∈ [floor, 1.0]
    def _gene_gate(e: float) -> float:
        if e >= 0.0:
            return 1.0
        if e <= down_threshold:
            return floor
        # Linear interpolation between (down_threshold, floor) and (0, 1)
        frac = (e - down_threshold) / (-down_threshold)
        return floor + (1.0 - floor) * frac

    out: dict[str, float] = {}
    for rxn_id, effects in rxn_genes.items():
        gates = [_gene_gate(e) for e in effects]
        if combine == "min":
            g_val = min(gates)
        elif combine == "max":
            g_val = max(gates)
        elif combine == "mean":
            # Geometric mean
            from math import prod
            g_val = prod(gates) ** (1.0 / len(gates))
        else:
            raise ValueError(f"unknown combine={combine!r}")
        out[rxn_id] = float(g_val)

    log.info(
        "compute_enzyme_gate: %d reactions gated; %d at floor (%.2f)",
        len(out), sum(1 for v in out.values() if v <= floor + 1e-9), floor,
    )
    return out


def apply_enzyme_gate(
    reaction_scores: dict[str, float],
    gate: dict[str, float],
    default: float = 1.0,
) -> dict[str, float]:
    """
    Apply ``gate`` (per-reaction multiplier) to ``reaction_scores``.

    Reactions absent from ``gate`` retain their original score (multiplied
    by ``default=1.0`` — pass-through).

    Returns a new dict; does not mutate inputs.
    """
    return {
        rxn_id: score * gate.get(rxn_id, default)
        for rxn_id, score in reaction_scores.items()
    }
