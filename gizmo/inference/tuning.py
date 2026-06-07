"""
Data-driven tuning of edge-type coupling constants for the BP model.

The default couplings shipped in :mod:`gizmo.inference.model` are informed
priors, not learned parameters. This module implements a simple gradient
descent procedure that refines them against held-out labelled pairs.

Problem setup
-------------
Given a dataset of *(observations, labels)* pairs, where:

  observations : dict[node_id → effect_size]     — what we gave BP
  labels       : dict[node_id → one_of({DOWN, NORMAL, UP})]  — what we expect

we want to pick coupling constants θ_t per edge type that maximise the
mean log-posterior at the labelled nodes after BP convergence:

    L(θ) = Σ_{pairs} Σ_{labelled_i} log P_BP(X_i = label_i | obs; θ)

Gradients of BP marginals w.r.t. θ don't have a closed form, so we use a
simple coordinate-wise finite-difference step. With only ~10–15 edge types
this converges in a few dozen BP runs per epoch.

This is a v1: we do grid-free coordinate descent with early stopping.
Pseudo-likelihood (per-node local evidence) would be faster but sacrifices
the very graph-propagation effect we care about.
"""

from __future__ import annotations

import logging
import math
from dataclasses import replace
from typing import Iterable, Optional

import numpy as np

from gizmo.inference.bp import run_bayesian_inference, BPConfig
from gizmo.inference.model import (
    DOWN, NORMAL, UP, Couplings, DEFAULT_COUPLINGS, DEFAULT_ASSAY_SIGMAS,
)
from gizmo.evidence.model import SampleContext, EvidenceRecord

log = logging.getLogger(__name__)


def _coerce_obs(nid: str, value) -> tuple[float, str]:
    """Normalise an obs entry into (effect_size, assay_type).

    Accepts a bare float (legacy, assumed metabolomics) or a
    ``(effect, assay_type)`` tuple/list.
    """
    if isinstance(value, (tuple, list)) and len(value) == 2:
        return float(value[0]), str(value[1])
    return float(value), "metabolomics"


def log_posterior_score(
    mg, pairs: Iterable[tuple[dict, dict]], cfg: BPConfig,
) -> float:
    """Mean log-posterior at labelled nodes across the pair list.

    ``pairs`` is an iterable of ``(observations, labels)`` where each
    ``observations`` dict maps ``node_id → effect_size`` *or*
    ``node_id → (effect_size, assay_type)``. The tuple form is required for
    tuning per-assay σ; the scalar form defaults to metabolomics.
    """
    total = 0.0
    count = 0
    for obs, labels in pairs:
        ctx = SampleContext(sample_id="tune")
        for nid, value in obs.items():
            eff, atype = _coerce_obs(nid, value)
            ctx.add(EvidenceRecord(
                feature_id=nid, node_id=nid,
                node_type=mg.graph.nodes.get(nid, {}).get("node_type") or "metabolite",
                effect_size=eff, confidence=1.0,
                assay_type=atype, sample_id="tune",
            ))
        res = run_bayesian_inference(mg, ctx, cfg)
        for nid, lab in labels.items():
            post = res.posteriors.get(nid)
            if post is None:
                continue
            p = max(1e-6, float(post[lab]))
            total += math.log(p)
            count += 1
    return total / max(count, 1)


def tune_couplings(
    mg,
    pairs: list[tuple[dict, dict]],
    *,
    base: Optional[Couplings] = None,
    step: float = 0.5,
    min_step: float = 0.05,
    edge_types: Optional[list[str]] = None,
    bp_iters: int = 40,
    max_epochs: int = 4,
) -> Couplings:
    """
    Coordinate-descent tune of coupling θ per edge type.

    Parameters
    ----------
    mg         : GizmoGraph
    pairs      : list of (observations, labels) per training example
    base       : starting coupling values (default DEFAULT_COUPLINGS)
    step       : initial coordinate step size
    min_step   : terminate when step size falls below this
    edge_types : restrict tuning to these edge types (default all in base)
    bp_iters   : max BP iterations per evaluation
    max_epochs : number of coordinate-descent passes

    Returns
    -------
    Refined Couplings instance.
    """
    current = Couplings(values=dict((base or DEFAULT_COUPLINGS).values))
    cfg = BPConfig(max_iter=bp_iters, log_every=0, couplings=current)

    targets = edge_types or list(current.values.keys())
    best_score = log_posterior_score(mg, pairs, cfg)
    log.info("Initial log-posterior = %.4f", best_score)

    for epoch in range(max_epochs):
        improved = False
        for t in targets:
            theta0, sign = current.values[t]
            for delta in (+step, -step):
                trial = theta0 + delta
                if trial < 0:
                    continue
                current.values[t] = (trial, sign)
                cfg_try = replace(cfg, couplings=current)
                score = log_posterior_score(mg, pairs, cfg_try)
                if score > best_score + 1e-4:
                    best_score = score
                    improved = True
                    log.info("  θ(%s) %.3f → %.3f  score=%.4f", t, theta0, trial, score)
                    break
            else:
                # Neither direction improved; revert
                current.values[t] = (theta0, sign)
        if not improved:
            step *= 0.5
            log.info("Epoch %d: no improvement; shrinking step to %.3f", epoch + 1, step)
            if step < min_step:
                break

    log.info("Final log-posterior = %.4f", best_score)
    return current


def tune_sigmas(
    mg,
    pairs: list[tuple[dict, dict]],
    *,
    couplings: Optional[Couplings] = None,
    base_sigmas: Optional[dict[str, float]] = None,
    assay_types: Optional[list[str]] = None,
    step: float = 0.5,
    min_step: float = 0.05,
    bp_iters: int = 40,
    max_epochs: int = 4,
) -> dict[str, float]:
    """
    Coordinate-descent tune of per-assay observation σ.

    Holds couplings fixed; varies the σ used in ``gaussian_obs_likelihood``
    per ``assay_type``. Observations in ``pairs`` must carry their assay type
    (pass ``obs`` as ``{node_id: (effect_size, assay_type)}``) or the
    per-assay signal is indistinguishable.

    Parameters
    ----------
    mg          : GizmoGraph
    pairs       : list of (observations, labels); see :func:`log_posterior_score`
    couplings   : edge-type θ to hold fixed (default DEFAULT_COUPLINGS)
    base_sigmas : starting σ dict (default DEFAULT_ASSAY_SIGMAS)
    assay_types : restrict tuning to these assays (default: all in base_sigmas)
    step        : initial coordinate step size on σ
    min_step    : terminate when step size falls below this
    bp_iters    : max BP iterations per evaluation
    max_epochs  : number of coordinate-descent passes

    Returns
    -------
    Refined ``assay_sigmas`` dict.
    """
    current = dict(base_sigmas or DEFAULT_ASSAY_SIGMAS)
    cpl = couplings or Couplings(values=dict(DEFAULT_COUPLINGS.values))
    cfg = BPConfig(
        max_iter=bp_iters, log_every=0,
        couplings=cpl, assay_sigmas=current,
    )

    targets = assay_types or list(current.keys())
    best_score = log_posterior_score(mg, pairs, cfg)
    log.info("Initial log-posterior (σ tune) = %.4f", best_score)

    for epoch in range(max_epochs):
        improved = False
        for t in targets:
            sigma0 = current[t]
            for delta in (+step, -step):
                trial = sigma0 + delta
                if trial <= 0.1:   # match the floor in gaussian_obs_likelihood
                    continue
                current[t] = trial
                cfg_try = replace(cfg, assay_sigmas=current)
                score = log_posterior_score(mg, pairs, cfg_try)
                if score > best_score + 1e-4:
                    best_score = score
                    improved = True
                    log.info("  σ(%s) %.3f → %.3f  score=%.4f", t, sigma0, trial, score)
                    break
            else:
                current[t] = sigma0
        if not improved:
            step *= 0.5
            log.info("Epoch %d (σ): no improvement; shrinking step to %.3f", epoch + 1, step)
            if step < min_step:
                break

    log.info("Final log-posterior (σ tune) = %.4f", best_score)
    return current


__all__ = ["tune_couplings", "tune_sigmas", "log_posterior_score"]
