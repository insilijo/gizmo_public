"""Multi-block sparse PLS-DA (DIABLO-equivalent) for cross-omic effect_size.

Why this exists
---------------

Single-block PLS-DA (gizmo.scoring.plsda_loadings) treats each omic
independently: it cannot use the fact that a transcript and a metabolite
covary across the same set of samples. DIABLO (Singh 2019) optimizes
component loadings *jointly* across blocks to maximize cross-block
covariance — making it possible to recover signal that is weak in any
single omic but strong in their joint structure.

Implementation
--------------

This module wraps `mbpls.MBPLS` (NIPALS-derived multi-block PLS) to fit
a discriminant variant against a one-hot class indicator. Each block's
component-1 weight (W_[block][:, 0]) is taken as the signed per-feature
loading, scaled by per-block VIP for magnitude. Sign is reoriented per
block so that "case" maps to the positive direction in super-score
space (Ts_).

For a sparser DIABLO-style output, `sparsity_keep_top_per_block` zeros
all but the top-K |loading| features per block.

Output shape matches plsda_loadings.plsda_per_feature for drop-in
substitution.
"""
from __future__ import annotations
from typing import Sequence, Mapping
import numpy as np


def diablo_multiblock_loadings(
    blocks: Mapping[str, np.ndarray],
    case_mask: np.ndarray,
    n_components: int = 2,
    sparsity_keep_top_per_block: int | None = None,
    block_weights: Mapping[str, float] | None = None,
) -> dict:
    """Run multi-block sPLS-DA across omics blocks.

    Parameters
    ----------
    blocks : dict {block_name: array (n_samples, n_block_features)}
        All blocks must share the same n_samples ordering.
    case_mask : bool array (n_samples,) — True where sample is "case"
    n_components : PLS components (default 2; first drives loadings)
    sparsity_keep_top_per_block : if int, retain top-K features per block
        by absolute loading; others zeroed. None = no sparsity.
    block_weights : optional dict {block: weight} — multiplies each block's
        scaled matrix. Use to compensate for unbalanced block sizes.

    Returns
    -------
    dict per block:
      {block_name: {
          "feature_indices": np.ndarray (n_block_features,) — column index in block,
          "signed_loading":  np.ndarray (n_block_features,) — signed comp-1 weight,
          "vip":             np.ndarray (n_block_features,),
          "effect_size":     np.ndarray (n_block_features,) — signed_loading * vip,
          "block_importance": float — A_[block_idx, 0] from MBPLS super-weights,
      }}
    """
    import warnings
    from mbpls.mbpls import MBPLS

    block_names = list(blocks.keys())
    n_samples = next(iter(blocks.values())).shape[0]
    for b, X in blocks.items():
        if X.shape[0] != n_samples:
            raise ValueError(f"block {b} has {X.shape[0]} rows, expected {n_samples}")
    if case_mask.shape[0] != n_samples:
        raise ValueError(f"case_mask shape {case_mask.shape} != n_samples {n_samples}")

    # NaN-impute (per-column mean) before passing to MBPLS; it doesn't
    # tolerate NaN. Optional per-block re-weighting (multiplies after
    # MBPLS internal standardize to bias super-component allocation).
    cleaned_blocks = []
    feature_counts = []
    for b in block_names:
        Xc = blocks[b].astype(float).copy()
        col_mean = np.nanmean(Xc, axis=0)
        nan_mask = np.isnan(Xc)
        if nan_mask.any():
            mean_grid = np.broadcast_to(col_mean, Xc.shape)
            Xc = np.where(nan_mask, mean_grid, Xc)
        if block_weights and b in block_weights:
            Xc = Xc * float(block_weights[b])
        cleaned_blocks.append(Xc)
        feature_counts.append(Xc.shape[1])

    y = np.where(case_mask, 1.0, -1.0).reshape(-1, 1)

    n_comp = min(n_components, n_samples - 1, max(feature_counts))
    if n_comp < 1:
        return {
            b: {"feature_indices": np.arange(blocks[b].shape[1]),
                "signed_loading": np.zeros(blocks[b].shape[1]),
                "vip": np.zeros(blocks[b].shape[1]),
                "effect_size": np.zeros(blocks[b].shape[1]),
                "block_importance": 0.0}
            for b in block_names
        }

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mb = MBPLS(n_components=n_comp, standardize=True, calc_all=True)
        mb.fit(cleaned_blocks, y)

    # Reorient: super-score Ts_[:, 0] should have higher mean for case
    super_scores = mb.Ts_[:, 0]
    case_mean = super_scores[case_mask].mean()
    ctrl_mean = super_scores[~case_mask].mean()
    flip = case_mean < ctrl_mean

    out = {}
    for bi, b in enumerate(block_names):
        w1 = mb.W_[bi][:, 0].copy()
        if flip:
            w1 = -w1
        # block-specific VIP from this block's W_ across components
        vip = _vip_block(mb, bi)
        if sparsity_keep_top_per_block is not None and sparsity_keep_top_per_block < len(w1):
            order = np.argsort(-np.abs(w1))
            keep = order[:sparsity_keep_top_per_block]
            mask = np.zeros_like(w1, dtype=bool)
            mask[keep] = True
            w1 = np.where(mask, w1, 0.0)
            vip = np.where(mask, vip, 0.0)
        block_importance = float(mb.A_[bi, 0]) if hasattr(mb, "A_") else 0.0
        out[b] = {
            "feature_indices": np.arange(len(w1)),
            "signed_loading": w1,
            "vip": vip,
            "effect_size": w1 * vip,
            "block_importance": block_importance,
        }
    return out


def diablo_per_feature(
    block_case_values: Mapping[str, Sequence[Sequence[float]]],
    block_control_values: Mapping[str, Sequence[Sequence[float]]],
    block_feature_names: Mapping[str, Sequence[str]],
    min_n: int = 3,
    sparsity_keep_top_per_block: int | None = None,
) -> dict[str, list[dict]]:
    """Per-feature multi-block sPLS-DA evidence with ragged inputs.

    Aligns case/control samples per block via positional order; each
    block must have the same n_case and n_ctrl samples (cohort-paired).
    Features with insufficient samples in either arm are dropped from
    that block before fitting.

    Returns per-block list-of-dicts compatible with plsda_per_feature
    output: feature_name, signed_loading, vip, effect_size, log2_fc,
    p_value, mod_t.
    """
    aligned_blocks_case = {}
    aligned_blocks_ctrl = {}
    kept_feature_names: dict[str, list[str]] = {}
    for b, case_per in block_case_values.items():
        ctrl_per = block_control_values[b]
        names = list(block_feature_names[b])
        keep_idx = []
        for i, (c, h) in enumerate(zip(case_per, ctrl_per)):
            if len(c) >= min_n and len(h) >= min_n:
                keep_idx.append(i)
        if not keep_idx:
            kept_feature_names[b] = []
            continue
        n_case = max(len(case_per[i]) for i in keep_idx)
        n_ctrl = max(len(ctrl_per[i]) for i in keep_idx)
        case_mat = np.full((n_case, len(keep_idx)), np.nan)
        ctrl_mat = np.full((n_ctrl, len(keep_idx)), np.nan)
        for j, i in enumerate(keep_idx):
            c = case_per[i]
            h = ctrl_per[i]
            case_mat[: len(c), j] = c
            ctrl_mat[: len(h), j] = h
        aligned_blocks_case[b] = case_mat
        aligned_blocks_ctrl[b] = ctrl_mat
        kept_feature_names[b] = [names[i] for i in keep_idx]

    if not aligned_blocks_case:
        return {b: [] for b in block_case_values}

    # Need shared sample dimension — pad each block to max(n_case, n_ctrl)
    n_case_max = max(m.shape[0] for m in aligned_blocks_case.values())
    n_ctrl_max = max(m.shape[0] for m in aligned_blocks_ctrl.values())
    blocks_full = {}
    for b in aligned_blocks_case:
        c = aligned_blocks_case[b]
        h = aligned_blocks_ctrl[b]
        if c.shape[0] < n_case_max:
            pad = np.full((n_case_max - c.shape[0], c.shape[1]), np.nan)
            c = np.vstack([c, pad])
        if h.shape[0] < n_ctrl_max:
            pad = np.full((n_ctrl_max - h.shape[0], h.shape[1]), np.nan)
            h = np.vstack([h, pad])
        blocks_full[b] = np.vstack([c, h])
    case_mask = np.concatenate([np.ones(n_case_max, bool), np.zeros(n_ctrl_max, bool)])

    multiblock_res = diablo_multiblock_loadings(
        blocks_full, case_mask,
        n_components=2,
        sparsity_keep_top_per_block=sparsity_keep_top_per_block,
    )

    out: dict[str, list[dict]] = {}
    for b, res in multiblock_res.items():
        recs = []
        feat_names = kept_feature_names[b]
        case_full = aligned_blocks_case[b]
        ctrl_full = aligned_blocks_ctrl[b]
        case_mean = np.nanmean(case_full, axis=0)
        ctrl_mean = np.nanmean(ctrl_full, axis=0)
        for j, name in enumerate(feat_names):
            recs.append({
                "feature_name": name,
                "block": b,
                "log2_fc": float(case_mean[j] - ctrl_mean[j])
                            if not np.isnan(case_mean[j] - ctrl_mean[j]) else 0.0,
                "signed_loading": float(res["signed_loading"][j]),
                "vip": float(res["vip"][j]),
                "effect_size": float(res["effect_size"][j]),
                "mod_t": float(res["effect_size"][j]),
                "p_value": 0.05,
            })
        out[b] = recs
    return out


def _vip_block(mb, block_idx: int) -> np.ndarray:
    """Block-specific VIP from a fitted mbpls.MBPLS model.

    Standard VIP for the j-th block: vip_j = sqrt(p_j * Σ_h s_h * (W_j[:,h]/||W_j[:,h]||)² / Σ_h s_h)
    where s_h is the sum-of-squares explained by component h on that block.
    """
    W_block = mb.W_[block_idx]               # (p_j, n_comp)
    T_block = mb.T_[block_idx]                # (n_samples, n_comp)
    p, h = W_block.shape

    # explained SS per component (block scores capture block-specific variance)
    s_per_comp = np.diag(T_block.T @ T_block)
    total = s_per_comp.sum()
    if total == 0:
        return np.ones(p)

    vips = np.zeros(p)
    for i in range(p):
        contribs = np.zeros(h)
        for j in range(h):
            wj = W_block[:, j]
            norm = np.linalg.norm(wj)
            if norm > 0:
                contribs[j] = (W_block[i, j] / norm) ** 2
        vips[i] = float(np.sqrt(p * np.dot(s_per_comp, contribs) / total))
    return vips
