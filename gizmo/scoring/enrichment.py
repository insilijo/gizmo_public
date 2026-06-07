"""
Statistical pathway enrichment tests.

Two methods are provided:

ORA (Over-Representation Analysis)
    Classic hypergeometric test: given a set of *hit* features and a
    background population, tests whether each pathway is more populated
    by hits than expected by chance.  BH FDR is applied across all
    pathways tested.

Pre-ranked enrichment
    A lightweight Kolmogorov-Smirnov style test over a ranked feature
    list (e.g. sorted by |fold-change| or -log10 p-value).  Reports the
    enrichment score (ES), normalised ES (NES), and a permutation-free
    p-value estimated from the leading-edge fraction and list depth.

Usage::

    from gizmo.scoring.enrichment import run_ora, run_preranked

    ora  = run_ora(hit_ids, background_ids, mg, pathway_names)
    prer = run_preranked(ranked_ids, ranked_scores, mg, pathway_names)
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _bh_fdr(pvalues: list[float]) -> list[float]:
    """Benjamini-Hochberg FDR correction.  Returns adjusted p-values (list of floats)."""
    n = len(pvalues)
    if n == 0:
        return []
    indexed = sorted(enumerate(pvalues), key=lambda x: x[1])
    adjusted = [0.0] * n
    prev = 1.0
    for rank, (orig_idx, pv) in enumerate(reversed(indexed), 1):
        adj = min(prev, pv * n / (n - rank))
        adjusted[orig_idx] = adj
        prev = adj
    return adjusted


def _pathway_to_nodes(mg, min_size: int = 3, max_size: int = 500) -> dict[str, set[str]]:
    """
    Build a {pathway_id: {node_id, ...}} map from the graph.

    Pathway membership is derived from reaction nodes whose ``pathways``
    attribute lists Reactome stIds.  For each such reaction, all
    neighbour metabolite and gene nodes are collected as pathway members.

    If the graph has first-class PathwayNode nodes (from add_pathway_nodes),
    those are used to harvest ``n_reactions`` but membership is still
    reaction-neighbour based for accuracy.
    """
    g = mg.graph

    # reaction → pathway stIds
    rxn_to_pathways: dict[str, list[str]] = {}
    for nid, attrs in g.nodes(data=True):
        if attrs.get("node_type") != "reaction":
            continue
        pws = attrs.get("pathways") or []
        if pws:
            rxn_to_pathways[nid] = pws

    # pathway → {member node_ids}
    pw_members: dict[str, set[str]] = {}
    for rxn_id, pws in rxn_to_pathways.items():
        for nbr in g.neighbors(rxn_id):
            nt = g.nodes[nbr].get("node_type", "")
            if nt in ("metabolite", "gene"):
                for pw in pws:
                    pw_members.setdefault(pw, set()).add(nbr)
        for pred in g.predecessors(rxn_id):
            nt = g.nodes[pred].get("node_type", "")
            if nt in ("metabolite", "gene"):
                for pw in pws:
                    pw_members.setdefault(pw, set()).add(pred)

    return {
        pw: nodes
        for pw, nodes in pw_members.items()
        if min_size <= len(nodes) <= max_size
    }


def _pathway_names_from_graph(mg) -> dict[str, str]:
    """Extract {pathway_id: name} from PathwayNode nodes if present."""
    names: dict[str, str] = {}
    g = mg.graph
    for nid, attrs in g.nodes(data=True):
        if attrs.get("node_type") == "pathway":
            stid = attrs.get("stId") or attrs.get("pathway_id") or nid
            name = attrs.get("name") or attrs.get("displayName") or stid
            names[stid] = name
    return names


# ---------------------------------------------------------------------------
# ORA
# ---------------------------------------------------------------------------

@dataclass
class ORAResult:
    """Result for one pathway from over-representation analysis."""
    pathway_id:     str
    pathway_name:   str   = ""
    n_hit:          int   = 0   # hits mapping to this pathway
    n_pathway:      int   = 0   # total pathway members in background
    n_background:   int   = 0   # total background features
    n_total_hits:   int   = 0   # total hit features
    p_value:        float = 1.0
    p_adj:          float = 1.0
    fold_enrichment: float = 1.0
    hit_features:   list  = field(default_factory=list)


def run_ora(
    hit_node_ids: set[str],
    background_node_ids: set[str],
    mg,
    pathway_names: Optional[dict[str, str]] = None,
    min_pathway_size: int = 3,
    max_pathway_size: int = 500,
    min_hits: int = 1,
) -> list[ORAResult]:
    """
    Over-Representation Analysis using the hypergeometric distribution.

    Parameters
    ----------
    hit_node_ids        : graph node IDs for the significant feature set
    background_node_ids : graph node IDs for all measured features
    mg                  : GizmoGraph instance
    pathway_names       : optional {stId: displayName} override; auto-derived
                          from the graph's PathwayNode nodes if None
    min_pathway_size    : ignore pathways with fewer members in background
    max_pathway_size    : ignore pathways with more members in background
    min_hits            : ignore pathways with fewer than this many hit members

    Returns
    -------
    list[ORAResult] sorted by p_adj ascending, then fold_enrichment descending
    """
    try:
        from scipy.stats import hypergeom
    except ImportError:
        raise ImportError("scipy is required for ORA — pip install scipy")

    names = pathway_names or _pathway_names_from_graph(mg)
    pw_members = _pathway_to_nodes(mg, min_size=min_pathway_size, max_size=max_pathway_size)

    N = len(background_node_ids)     # population size
    K_total = len(hit_node_ids)      # total successes in population

    results: list[ORAResult] = []
    for pw_id, members in pw_members.items():
        # Restrict pathway to background
        pw_in_bg  = members & background_node_ids
        pw_hits   = members & hit_node_ids
        K_pathway = len(pw_in_bg)
        k_hit     = len(pw_hits)

        if K_pathway < min_pathway_size or k_hit < min_hits:
            continue

        # P(X >= k_hit) = hypergeom survival function
        pv = hypergeom.sf(k_hit - 1, N, K_pathway, K_total)

        expected = (K_total * K_pathway) / max(N, 1)
        fe = k_hit / expected if expected > 0 else 0.0

        results.append(ORAResult(
            pathway_id      = pw_id,
            pathway_name    = names.get(pw_id, ""),
            n_hit           = k_hit,
            n_pathway       = K_pathway,
            n_background    = N,
            n_total_hits    = K_total,
            p_value         = round(float(pv), 6),
            fold_enrichment = round(fe, 3),
            hit_features    = sorted(pw_hits),
        ))

    if not results:
        return []

    # BH FDR
    padjs = _bh_fdr([r.p_value for r in results])
    for r, adj in zip(results, padjs):
        r.p_adj = round(adj, 6)

    results.sort(key=lambda r: (r.p_adj, -r.fold_enrichment))
    return results


# ---------------------------------------------------------------------------
# Pre-ranked enrichment (lightweight KS-style)
# ---------------------------------------------------------------------------

@dataclass
class PrerankedResult:
    """Result for one pathway from pre-ranked enrichment."""
    pathway_id:      str
    pathway_name:    str   = ""
    enrichment_score: float = 0.0   # leading-edge ES (signed)
    nes:             float = 0.0    # normalised ES (ES / mean(|ES|) across permutations)
    p_value:         float = 1.0    # asymptotic p-value from KS distribution
    p_adj:           float = 1.0    # BH-corrected
    n_hit:           int   = 0
    n_pathway:       int   = 0
    leading_edge:    list  = field(default_factory=list)


def run_preranked(
    ranked_node_ids: list[str],
    ranked_scores: list[float],
    mg,
    pathway_names: Optional[dict[str, str]] = None,
    min_pathway_size: int = 3,
    max_pathway_size: int = 500,
    p_exponent: float = 1.0,
) -> list[PrerankedResult]:
    """
    Pre-ranked enrichment analysis (KS-style enrichment score).

    Parameters
    ----------
    ranked_node_ids : node IDs sorted by descending score (most positive first)
    ranked_scores   : corresponding scores (same order); used to weight hits
    mg              : GizmoGraph instance
    p_exponent      : weight exponent for hit scores (1.0 = classic weighted KS)

    Returns
    -------
    list[PrerankedResult] sorted by |enrichment_score| descending
    """
    try:
        from scipy.stats import ks_1samp
    except ImportError:
        raise ImportError("scipy is required for pre-ranked enrichment — pip install scipy")

    names = pathway_names or _pathway_names_from_graph(mg)
    pw_members = _pathway_to_nodes(mg, min_size=min_pathway_size, max_size=max_pathway_size)

    N = len(ranked_node_ids)
    score_arr = list(ranked_scores)

    results: list[PrerankedResult] = []

    for pw_id, members in pw_members.items():
        # Positions of pathway members in ranked list
        hit_idx = [i for i, nid in enumerate(ranked_node_ids) if nid in members]
        if len(hit_idx) < min_pathway_size:
            continue

        n_hit = len(hit_idx)
        n_miss = N - n_hit

        # Weighted running sum
        hit_weights = [abs(score_arr[i]) ** p_exponent for i in hit_idx]
        total_weight = sum(hit_weights) or 1.0
        miss_penalty = 1.0 / max(n_miss, 1)

        running = 0.0
        es = 0.0
        peak_idx = 0
        member_set = set(hit_idx)
        hw_idx = {i: w / total_weight for i, w in zip(hit_idx, hit_weights)}

        for pos in range(N):
            if pos in member_set:
                running += hw_idx[pos]
            else:
                running -= miss_penalty
            if abs(running) > abs(es):
                es = running
                peak_idx = pos

        # Leading edge: hits before the peak
        leading = [ranked_node_ids[i] for i in hit_idx if i <= peak_idx]

        # Asymptotic p-value from KS distribution (|ES| scaled by sqrt(N))
        # KS statistic D = |ES| — use scipy KS two-sided cdf
        ks_stat = abs(es)
        # Approximate p from Kolmogorov distribution
        lam = ks_stat * math.sqrt(n_hit)
        if lam <= 0:
            pv = 1.0
        else:
            # P-value: 2 * sum_{k=1}^{inf} (-1)^{k-1} exp(-2 k^2 lam^2)  (Kolmogorov distribution)
            pv = 0.0
            for k in range(1, 20):
                term = ((-1) ** (k - 1)) * math.exp(-2 * k * k * lam * lam)
                pv += term
                if abs(term) < 1e-10:
                    break
            pv = max(0.0, min(1.0, 2 * pv))

        results.append(PrerankedResult(
            pathway_id       = pw_id,
            pathway_name     = names.get(pw_id, ""),
            enrichment_score = round(es, 4),
            p_value          = round(pv, 6),
            n_hit            = n_hit,
            n_pathway        = len(members),
            leading_edge     = leading[:20],
        ))

    if not results:
        return []

    # Normalise ES: divide by mean |ES| across all pathways
    mean_abs_es = sum(abs(r.enrichment_score) for r in results) / len(results)
    for r in results:
        r.nes = round(r.enrichment_score / mean_abs_es, 3) if mean_abs_es > 0 else 0.0

    # BH FDR
    padjs = _bh_fdr([r.p_value for r in results])
    for r, adj in zip(results, padjs):
        r.p_adj = round(adj, 6)

    results.sort(key=lambda r: -abs(r.enrichment_score))
    return results
