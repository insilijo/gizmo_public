"""Signed PLS-DA loadings as effect_size for the Laplacian framework.

GrAndMA's `run_plsda` returns VIP scores (always positive). For the
Laplacian we need *signed* per-feature evidence — the sign carries
direction-of-change relative to the case/control contrast.

This module wraps sklearn's PLSRegression to expose:
- ``signed_loading``: x_weights_[:, 0] reoriented so positive = "up in
  case" relative to the encoded class label
- ``vip``: standard Variable Importance in Projection (magnitude only)
- ``effect_size``: signed_loading × vip (signed magnitude — usable as
  drop-in effect_size in the same scale as moderated-t)

The first PLS component captures the direction of maximum class
separation. Sign of x_weights_[:, 0] tells us whether each feature
goes up in case or control. We disambiguate orientation by checking
the mean component-1 score per class and flipping if needed so that
"case" maps to the positive direction.
"""
from __future__ import annotations
from typing import Sequence
import numpy as np


def moderated_t_panel_via_plsda(
    case_matrix: np.ndarray,
    control_matrix: np.ndarray,
    n_components: int = 2,
    use_vip_weighting: bool = True,
) -> dict:
    """Run PLS-DA and return signed loadings + VIP per feature.

    Parameters
    ----------
    case_matrix : array (n_case, n_features)
    control_matrix : array (n_control, n_features)
    n_components : PLS components to fit (default 2; first is signed-load source)
    use_vip_weighting : if True, return effect_size = signed_loading * VIP;
        otherwise return signed_loading directly. VIP weighting amplifies
        high-importance features and dampens low-importance ones — closer
        to limma in spirit.

    Returns
    -------
    dict with arrays of shape (n_features,):
      ``signed_loading``  — first-component x_weights_, oriented so case > 0
      ``vip``             — Variable Importance in Projection (magnitude)
      ``effect_size``     — signed_loading * vip if use_vip_weighting else signed_loading
      ``log2_fc``         — case_mean - control_mean (in input units)
      ``p_value``         — placeholder 0.05; PLS-DA has no per-feature p
    """
    from sklearn.cross_decomposition import PLSRegression
    from sklearn.preprocessing import StandardScaler

    case = np.asarray(case_matrix, dtype=float)
    ctrl = np.asarray(control_matrix, dtype=float)
    if case.shape[1] != ctrl.shape[1]:
        raise ValueError(
            f"feature dim mismatch: case {case.shape[1]} vs ctrl {ctrl.shape[1]}"
        )
    n_feat = case.shape[1]

    X = np.vstack([case, ctrl])
    # case = +1, control = -1 — sign convention makes "up in case" positive
    y = np.concatenate([np.ones(len(case)), -np.ones(len(ctrl))])

    # NaN imputation to feature mean (PLS doesn't tolerate NaN)
    col_mean = np.nanmean(X, axis=0)
    nan_mask = np.isnan(X)
    if nan_mask.any():
        X = X.copy()
        X[nan_mask] = np.take(col_mean, np.where(nan_mask)[1])

    X_scaled = StandardScaler().fit_transform(X)

    n_comp = min(n_components, X_scaled.shape[0] - 1, n_feat)
    if n_comp < 1:
        return {
            "signed_loading": np.zeros(n_feat),
            "vip": np.zeros(n_feat),
            "effect_size": np.zeros(n_feat),
            "log2_fc": np.zeros(n_feat),
            "p_value": np.full(n_feat, 0.5),
        }

    pls = PLSRegression(n_components=n_comp, scale=False)
    pls.fit(X_scaled, y)

    w1 = pls.x_weights_[:, 0]  # (n_features,) — signed first-component loading

    # Orient so case has higher mean component-1 score than control.
    scores = pls.transform(X_scaled)[:, 0]
    case_mean = scores[: len(case)].mean()
    ctrl_mean = scores[len(case):].mean()
    if case_mean < ctrl_mean:
        w1 = -w1
        scores = -scores

    vip = _vip_scores(pls)

    if use_vip_weighting:
        effect_size = w1 * vip
    else:
        effect_size = w1.copy()

    case_raw_mean = np.nanmean(case, axis=0)
    ctrl_raw_mean = np.nanmean(ctrl, axis=0)
    log2_fc = case_raw_mean - ctrl_raw_mean

    return {
        "signed_loading": w1,
        "vip": vip,
        "effect_size": effect_size,
        "log2_fc": log2_fc,
        "p_value": np.full(n_feat, 0.05),
    }


def plsda_per_feature(
    case_values: Sequence[Sequence[float]],
    control_values: Sequence[Sequence[float]],
    feature_names: Sequence[str],
    min_n: int = 3,
    use_vip_weighting: bool = True,
) -> list[dict]:
    """Per-feature PLS-DA evidence given ragged per-feature value lists.

    Aligns inputs (drops samples with any-feature missingness for the PLS
    fit, but reports per-feature log2_fc on full populations).

    Parameters
    ----------
    case_values : list of length n_features, each a list of case values
    control_values : list of length n_features, each a list of ctrl values
    feature_names : list of length n_features
    min_n : require ≥min_n in each arm; otherwise feature is dropped
    use_vip_weighting : pass-through to moderated_t_panel_via_plsda

    Returns
    -------
    list of dicts: feature_name, log2_fc, signed_loading, vip, effect_size,
    p_value, mod_t (= effect_size, for evidence-model uniformity)
    """
    keep_idx = []
    for i, (c, h) in enumerate(zip(case_values, control_values)):
        if len(c) >= min_n and len(h) >= min_n:
            keep_idx.append(i)

    if not keep_idx:
        return []

    n_case = max(len(case_values[i]) for i in keep_idx)
    n_ctrl = max(len(control_values[i]) for i in keep_idx)
    n_feat = len(keep_idx)

    case_mat = np.full((n_case, n_feat), np.nan)
    ctrl_mat = np.full((n_ctrl, n_feat), np.nan)
    for j, i in enumerate(keep_idx):
        c = case_values[i]
        h = control_values[i]
        case_mat[: len(c), j] = c
        ctrl_mat[: len(h), j] = h

    res = moderated_t_panel_via_plsda(
        case_mat, ctrl_mat,
        n_components=2,
        use_vip_weighting=use_vip_weighting,
    )

    out = []
    for j, i in enumerate(keep_idx):
        out.append({
            "feature_name": feature_names[i],
            "log2_fc": float(res["log2_fc"][j]),
            "signed_loading": float(res["signed_loading"][j]),
            "vip": float(res["vip"][j]),
            "effect_size": float(res["effect_size"][j]),
            "mod_t": float(res["effect_size"][j]),
            "p_value": float(res["p_value"][j]),
        })
    return out


def _vip_scores(model) -> np.ndarray:
    """Variable Importance in Projection for fitted PLSRegression."""
    t = model.x_scores_
    w = model.x_weights_
    q = model.y_loadings_

    p, h = w.shape
    s = np.diag(t.T @ t @ q.T @ q)
    total = s.sum()
    if total == 0:
        return np.ones(p)

    vips = np.zeros(p)
    for i in range(p):
        w_norm = np.array([
            (w[i, j] / np.linalg.norm(w[:, j])) ** 2 for j in range(h)
        ])
        vips[i] = float(np.sqrt(p * np.dot(s, w_norm) / total))
    return vips
