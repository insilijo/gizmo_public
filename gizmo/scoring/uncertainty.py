"""
Bootstrap-based uncertainty estimation for causal chain ranking.

The deterministic chain ranker produces a single score per chain. In a small
or noisy cohort, that single number is misleading: a chain may have ranked
high because one or two features happened to push its reaction hard. By
resampling the evidence set with replacement we can tell whether a chain is
"stable" (recurs across bootstraps) or "fragile" (appears only when specific
features are drawn).

Each bootstrap redraws the ``EvidenceRecord`` list, reruns
``score_reactions`` + ``rank_chains``, and records the chain's score + rank
per iteration. We then emit:

  - score_mean / score_lo / score_hi : 95% percentile CI
  - stability                         : fraction of bootstraps the chain
                                        appeared in the top-N
  - rank_mean                         : average position when it appeared

Typical use: n_boot=50 and top_n=30 — a few seconds of runtime on a small
sample, much faster than a permutation null.
"""

from __future__ import annotations

import random
from typing import Optional

from gizmo.evidence.model import SampleContext


def bootstrap_chain_uncertainty(
    mg,
    ctx: SampleContext,
    *,
    n_boot: int = 50,
    top_n: int = 30,
    length_decay: float = 0.80,
    min_reaction_score: float = 0.0,
    seed: Optional[int] = 17,
) -> dict[tuple, dict]:
    """
    Bootstrap-resample ``ctx``'s evidence records and report per-chain
    score distribution + rank stability.

    Parameters
    ----------
    mg              : GizmoGraph
    ctx             : SampleContext containing evidence records
    n_boot          : number of bootstrap iterations (default 50)
    top_n           : keep top-N chains per bootstrap (default 30)
    length_decay    : passed to rank_chains
    min_reaction_score : passed to score_reactions filter
    seed            : RNG seed for reproducibility; set None for non-deterministic

    Returns
    -------
    dict keyed on ``tuple(chain.nodes)`` mapping to stats:
      - n_appearances, n_boot, stability
      - score_mean, score_lo, score_hi   (95% percentile CI)
      - rank_mean
    """
    from gizmo.scoring.reaction_scorer import score_reactions
    from gizmo.scoring.chain_ranker import rank_chains

    rng = random.Random(seed) if seed is not None else random

    base_records = ctx.records(mapped_only=True)
    n = len(base_records)
    if n == 0 or n_boot <= 0:
        return {}

    accum: dict[tuple, dict] = {}

    for _ in range(n_boot):
        sampled = [rng.choice(base_records) for _ in range(n)]
        boot_ctx = SampleContext(sample_id=ctx.sample_id, cohort_id=ctx.cohort_id)
        boot_ctx.add_all(sampled)

        try:
            boot_rxns = score_reactions(mg, boot_ctx, min_evidence=1)
        except Exception:
            continue
        try:
            boot_chains = rank_chains(
                mg, boot_rxns, ctx=boot_ctx,
                length_decay=length_decay,
                min_reaction_score=min_reaction_score,
            )
        except Exception:
            continue

        for rank_idx, ch in enumerate(boot_chains[:top_n]):
            key = tuple(ch.nodes)
            entry = accum.setdefault(key, {"scores": [], "ranks": []})
            entry["scores"].append(float(ch.score))
            entry["ranks"].append(rank_idx)

    out: dict[tuple, dict] = {}
    for key, e in accum.items():
        scores = sorted(e["scores"])
        ranks = e["ranks"]
        n_app = len(scores)
        if n_app >= 5:
            lo_idx = max(0, int(n_app * 0.025))
            hi_idx = min(n_app - 1, int(n_app * 0.975))
            lo, hi = scores[lo_idx], scores[hi_idx]
        else:
            lo, hi = scores[0], scores[-1]
        out[key] = {
            "n_appearances": n_app,
            "n_boot":        n_boot,
            "stability":     round(n_app / n_boot, 2),
            "score_mean":    round(sum(scores) / n_app, 4),
            "score_lo":      round(lo, 4),
            "score_hi":      round(hi, 4),
            "rank_mean":     round(sum(ranks) / n_app, 1),
        }
    return out
