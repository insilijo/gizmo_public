"""Limma-style moderated t-test (Smyth 2004 empirical Bayes).

Vectorized implementation that takes per-feature case/control values
and returns a moderated t-statistic per feature. The moderated
statistic shrinks per-feature variance estimates toward a global
log-variance prior estimated from the full panel — much more stable
for small N than Welch's t-test.

Why this matters for the GIZMO Laplacian framework: when we feed
``effect_size = log2_fc`` directly into the Laplacian, low-variance
features with tiny fold-changes contribute as much as high-variance
features with large fold-changes. Using ``effect_size = moderated_t``
weighs evidence by statistical confidence — features with consistent
direction across replicates count proportionally more.

Adapted from gold-standard limma (R/Bioconductor) and mirrors the
Python port already in `grandma_metabolomics.stats.run_ttest`. Kept
self-contained so GIZMO benchmarks don't import GrAndMA.
"""
from __future__ import annotations

import numpy as np
from scipy import stats
from scipy.special import digamma, polygamma


def moderated_t_panel(
    case_matrix: np.ndarray,
    control_matrix: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute limma-moderated t-statistic per feature.

    Parameters
    ----------
    case_matrix
        2D array, samples (case) × features
    control_matrix
        2D array, samples (control) × features. Same column order.

    Returns
    -------
    log2_fc, mod_t, p_value
        Each shape (n_features,). NaN where insufficient data.
    """
    if case_matrix.shape[1] != control_matrix.shape[1]:
        raise ValueError("case and control must have same n_features")

    n_features = case_matrix.shape[1]
    log2_fc = np.full(n_features, np.nan)
    mod_t = np.full(n_features, np.nan)
    p_value = np.full(n_features, np.nan)

    n1_f = np.sum(~np.isnan(case_matrix), axis=0).astype(float)
    n2_f = np.sum(~np.isnan(control_matrix), axis=0).astype(float)
    d1 = np.maximum(n1_f - 1, 0)
    d2 = np.maximum(n2_f - 1, 0)
    d_pool = d1 + d2

    mean1 = np.nanmean(case_matrix, axis=0)
    mean2 = np.nanmean(control_matrix, axis=0)
    log2_fc[:] = mean1 - mean2  # already-log values: difference IS log fold-change

    var1 = np.nanvar(case_matrix, axis=0, ddof=1)
    var2 = np.nanvar(control_matrix, axis=0, ddof=1)
    s2_pool = np.where(
        d_pool > 0,
        (d1 * var1 + d2 * var2) / np.maximum(d_pool, 1),
        np.nan,
    )

    valid_v = np.isfinite(s2_pool) & (s2_pool > 0) & (d_pool > 0)
    if valid_v.sum() < 3:
        # Not enough features for empirical Bayes; fall back to plain t
        with np.errstate(divide="ignore", invalid="ignore"):
            se = np.sqrt(s2_pool * (1.0 / np.maximum(n1_f, 1)
                                     + 1.0 / np.maximum(n2_f, 1)))
            mod_t = np.where(se > 0, log2_fc / se, np.nan)
            p_value = np.where(np.isfinite(mod_t),
                                2.0 * stats.t.sf(np.abs(mod_t),
                                                  df=np.maximum(d_pool, 1)),
                                np.nan)
        return log2_fc, mod_t, p_value

    # Estimate prior hyperparameters d0 and s0² from log(s²) distribution
    # using Smyth's method-of-moments on the log-variance distribution
    log_s2 = np.log(s2_pool[valid_v])
    d_i = d_pool[valid_v]
    mean_correction = np.mean(digamma(d_i / 2) - np.log(d_i / 2))
    log_s02 = np.mean(log_s2) - mean_correction
    s02 = float(np.exp(log_s02))
    var_log_s2 = float(np.var(log_s2, ddof=1))
    var_correct = float(np.mean(polygamma(1, d_i / 2)))
    excess_var = var_log_s2 - var_correct
    if excess_var > 0:
        try:
            from scipy.optimize import brentq
            d0 = float(2 * brentq(
                lambda x: float(polygamma(1, x)) - excess_var,
                1e-6, 1e6,
            ))
        except Exception:
            d0 = 4.0
    else:
        d0 = 1e6  # essentially no shrinkage needed (low between-feature variability)

    # Posterior moderated variance and t-statistic
    s2_mod = np.where(
        d_pool > 0,
        (d_pool * s2_pool + d0 * s02) / (d_pool + d0),
        s02,
    )
    se_mod = np.sqrt(
        s2_mod * (1.0 / np.maximum(n1_f, 1) + 1.0 / np.maximum(n2_f, 1))
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        mod_t = np.where(se_mod > 0, log2_fc / se_mod, np.nan)
    df_mod = d_pool + d0
    p_value = np.where(
        np.isfinite(mod_t),
        2.0 * stats.t.sf(np.abs(mod_t), df=df_mod),
        np.nan,
    )

    return log2_fc, mod_t, p_value


def moderated_t_per_feature(
    case_values: list[list[float]],
    control_values: list[list[float]],
    feature_names: list[str],
    min_n: int = 3,
) -> list[dict]:
    """Convenience wrapper for the loader pattern: take per-feature
    case/control value lists (jagged — different features may have
    different sample coverage), pad with NaN, run moderated t, return
    a panel-of-records suitable for the Laplacian.

    Each record has: feature_name, log2_fc, mod_t, p_value. Records
    with too few samples (<min_n per group) are dropped.
    """
    # Determine max group sizes
    n1_max = max((len(v) for v in case_values), default=0)
    n2_max = max((len(v) for v in control_values), default=0)
    n_feats = len(feature_names)

    # Pad with NaN
    case_mat = np.full((n1_max, n_feats), np.nan)
    ctrl_mat = np.full((n2_max, n_feats), np.nan)
    for j, v in enumerate(case_values):
        case_mat[:len(v), j] = v
    for j, v in enumerate(control_values):
        ctrl_mat[:len(v), j] = v

    log2_fc, mod_t, p_value = moderated_t_panel(case_mat, ctrl_mat)

    out = []
    for j, name in enumerate(feature_names):
        n_case = len(case_values[j])
        n_ctrl = len(control_values[j])
        if n_case < min_n or n_ctrl < min_n:
            continue
        if not np.isfinite(mod_t[j]):
            continue
        out.append({
            "feature_name": name,
            "log2_fc": float(log2_fc[j]),
            "mod_t": float(mod_t[j]),
            "p_value": float(p_value[j]) if np.isfinite(p_value[j]) else 1.0,
        })
    return out
