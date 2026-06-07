"""baseline_mofa_streaming.py — memory-constrained MOFA+ runs.

Two strategies for the 5 panel cohorts that don't fit standard MOFA+ on a
WSL2 / 16 GB host:

  (A) Single-block cohorts (GSE65391_SLE, GSE65682_sepsis): under a
      Gaussian likelihood and no ARD, single-view MOFA+ collapses to
      probabilistic PCA. We substitute `sklearn.decomposition.IncrementalPCA`,
      which streams in batches of N samples × full features. The output
      JSON is shaped identically to mofapy2's so downstream tooling
      (multi_pc_vs_mofa_factors.py) picks it up without modification.

  (B) Multi-block cohorts (CPTAC trio): patient-subsample MOFA+ — fit
      mofapy2 on a random N_sub patient subset (peak memory bounded by
      N_sub × Σ block features), then project the remaining patients onto
      the frozen weights via OLS posterior mean. Standard out-of-sample
      projection from the MOFA+ literature.

Outputs land in `benchmarks/results/unsupervised/mofa_weights/` under the
same `mofa_weights_<cohort>.json` filename so the augmenter pipeline is
unchanged.

Invocation:
  python3 benchmarks/baseline_mofa_streaming.py [COHORT_NAME]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "benchmarks"))

OUT_DIR = REPO / "benchmarks" / "results" / "unsupervised" / "mofa_weights"
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Shared utilities
# ---------------------------------------------------------------------------

def _block_to_matrix(data, common, substrate_mappable=None):
    """{sid: {feat: val}} + ordered sample list → (samples × features) ndarray
    + feature name list.

    If `substrate_mappable` is provided (set of feature names that map to
    substrate nodes), restrict the feature universe to that intersection —
    this is the apples-to-apples version that gives MOFA+ the same input
    universe GIZMO sees, eliminating the "MOFA+ got extra non-substrate
    features" advantage.
    """
    all_feats = sorted({f for s in common if s in data for f in data[s]})
    if substrate_mappable is not None:
        all_feats = [f for f in all_feats if f in substrate_mappable]
    X = np.zeros((len(common), len(all_feats)), dtype=np.float32)
    for i, sid in enumerate(common):
        if sid not in data: continue
        row = data[sid]
        for j, f in enumerate(all_feats):
            X[i, j] = row.get(f, 0.0)
    return X, all_feats


_SUBSTRATE_FEATURE_CACHE = None


def _load_substrate_mappable_features():
    """Return the set of gene-symbol strings that map to substrate nodes
    (~16k gene-substrate nodes after the human_full graph filter). Loaded
    lazily and cached because this can be a slow import."""
    global _SUBSTRATE_FEATURE_CACHE
    if _SUBSTRATE_FEATURE_CACHE is not None:
        return _SUBSTRATE_FEATURE_CACHE
    from gizmo.export.json_export import read_json
    print("    [substrate-matched] loading human_full substrate to identify "
          "mappable gene-symbol features…", flush=True)
    mg = read_json(REPO / "data/processed/human_full/graph.json")
    feats = set()
    for nid, attrs in mg.graph.nodes(data=True):
        if attrs.get("node_type") != "gene":
            continue
        # Substrate gene-node identifier is the gene symbol or attribute
        sym = attrs.get("symbol")
        if sym: feats.add(str(sym))
        nm = attrs.get("name")
        if nm: feats.add(str(nm))
        # Also include the bare nid stripped of any prefix
        if ":" in nid:
            feats.add(nid.split(":", 1)[1])
        else:
            feats.add(nid)
    _SUBSTRATE_FEATURE_CACHE = feats
    print(f"    [substrate-matched] {len(feats)} gene-symbol candidates "
          "for MOFA+ feature universe", flush=True)
    return feats


def _variance_filter(X, feat_names, max_features):
    if X.shape[1] <= max_features:
        return X, feat_names
    var = np.nanvar(X, axis=0)
    keep = np.argsort(var)[-max_features:]
    keep.sort()
    return X[:, keep], [feat_names[j] for j in keep]


# ---------------------------------------------------------------------------
# Strategy A — IncrementalPCA for single-block cohorts
# ---------------------------------------------------------------------------

def fit_incremental_pca_single_block(data, common, n_factors=10,
                                        max_features=10000, batch_size=64,
                                        block_name="rna",
                                        substrate_matched=False):
    """sklearn IncrementalPCA on a single-block input. Returns the mofa_weights
    JSON schema so the downstream augmenter doesn't care which method ran."""
    from sklearn.decomposition import IncrementalPCA

    substrate = _load_substrate_mappable_features() if substrate_matched else None
    X, feats = _block_to_matrix(data, common, substrate_mappable=substrate)
    X, feats = _variance_filter(X, feats, max_features)
    print(f"  IncrementalPCA: {X.shape[0]} samples × {X.shape[1]} features, "
          f"n_factors={n_factors}, batch={batch_size}", flush=True)

    # IncrementalPCA streams: partial_fit on chunks, then transform.
    n_components = min(n_factors, X.shape[0] - 1, X.shape[1])
    ipca = IncrementalPCA(n_components=n_components, batch_size=batch_size)
    n_samples = X.shape[0]
    for start in range(0, n_samples, batch_size):
        chunk = X[start:start + batch_size]
        if chunk.shape[0] < n_components:
            # Last chunk may be too small to partial_fit on its own; merge
            # with the prior chunk for the final pass.
            if start == 0:
                continue
            chunk = np.vstack([X[max(0, start - batch_size):start], chunk])
        ipca.partial_fit(chunk)
    Z = ipca.transform(X)                       # (n_samples × n_factors)
    W = ipca.components_.T                       # (n_features × n_factors)

    weights = {block_name: {feats[j]: W[j, :].tolist()
                              for j in range(W.shape[0])}}
    top_features = {block_name: {}}
    for f in range(W.shape[1]):
        abs_w = np.abs(W[:, f])
        order = np.argsort(abs_w)[::-1][:20]
        top_features[block_name][int(f)] = [
            (feats[j], float(abs_w[j])) for j in order
        ]

    return {
        "n_factors": int(Z.shape[1]),
        "samples": list(common),
        "factor_scores": Z.tolist(),
        "factor_variance_explained": {
            block_name: ipca.explained_variance_ratio_.tolist()
        },
        "weights": weights,
        "top_features_per_factor": top_features,
        "_method": "IncrementalPCA (single-view MOFA+ ≡ probabilistic PCA)",
    }


# ---------------------------------------------------------------------------
# Strategy B — patient-subsample MOFA+ for multi-block cohorts
# ---------------------------------------------------------------------------

def fit_subsample_mofa_multiblock(block_dicts, common, n_factors=10,
                                     n_subsample=80, max_features_per_block=8000,
                                     seed=0, substrate_matched=False):
    """Multi-block MOFA+ with patient-subsample fitting + OLS out-of-sample
    projection. Memory budget bounded by N_subsample × Σ feature counts."""
    from mofapy2.run.entry_point import entry_point

    substrate = _load_substrate_mappable_features() if substrate_matched else None
    rng = np.random.default_rng(seed)
    n_total = len(common)
    n_sub = min(n_subsample, n_total)
    sub_idx = rng.choice(n_total, size=n_sub, replace=False)
    sub_idx.sort()
    sub_samples = [common[i] for i in sub_idx]
    rest_samples = [s for i, s in enumerate(common) if i not in set(sub_idx)]
    print(f"  Subsample fit: {n_sub}/{n_total} patients; projecting "
          f"{len(rest_samples)} held-out", flush=True)

    # Build per-block matrices (variance-filtered) on the FULL cohort —
    # both training subset and projection set use the same feature set.
    block_names = []
    block_feats = {}
    block_full = {}
    for name, data in block_dicts.items():
        if not data: continue
        X, feats = _block_to_matrix(data, common, substrate_mappable=substrate)
        X, feats = _variance_filter(X, feats, max_features_per_block)
        block_full[name] = X
        block_feats[name] = feats
        block_names.append(name)
        print(f"    {name}: {X.shape[0]} × {X.shape[1]} features", flush=True)

    # Fit MOFA+ on the subsample
    train_data = [[block_full[n][sub_idx]] for n in block_names]
    ent = entry_point()
    ent.set_data_options(scale_views=True, scale_groups=False)
    ent.set_data_matrix(train_data, views_names=block_names,
                          groups_names=["all"], samples_names=[sub_samples])
    n_factors_eff = max(2, min(n_factors, n_sub - 2))
    ent.set_model_options(factors=n_factors_eff, spikeslab_weights=True,
                            ard_factors=True, ard_weights=True)
    ent.set_train_options(iter=300, convergence_mode="medium",
                            dropR2=0.001, gpu_mode=False, seed=seed,
                            verbose=False)
    ent.build()
    ent.run()
    expect = ent.model.getExpectations()
    Z_sub_raw = expect["Z"]["E"]
    Z_sub = np.asarray(Z_sub_raw[0]) if isinstance(Z_sub_raw, list) \
            else np.asarray(Z_sub_raw)
    if Z_sub.ndim == 1: Z_sub = Z_sub.reshape(-1, 1)
    Ws = [np.asarray(expect["W"][m]["E"]) for m in range(len(block_names))]
    # W shape per block = (n_features × n_factors)

    # Out-of-sample projection. For each held-out sample:
    #   Z_test = (Σ_b X_b W_b) @ inv(Σ_b W_b^T W_b)
    # where each block has been scaled the same way MOFA+ scaled it
    # internally — we re-apply per-feature mean centering (skipping scale
    # since MOFA+'s `scale_views` already normalizes block variance and we
    # don't have access to its internal scaler post-fit; this gives the
    # OLS-projection approximation, the standard out-of-sample formula in
    # the MOFA+ tutorial).
    WtW = sum(W.T @ W for W in Ws)
    WtW_inv = np.linalg.pinv(WtW)

    # Center each full block by the SUBSAMPLE mean (same centering MOFA+
    # would have applied internally).
    Z_full = np.zeros((n_total, Z_sub.shape[1]))
    Z_full[sub_idx] = Z_sub
    for i_glob in range(n_total):
        if i_glob in set(sub_idx): continue
        # Build X W contribution across blocks for this held-out sample
        XW = np.zeros(Z_sub.shape[1])
        for b, name in enumerate(block_names):
            X_full_b = block_full[name]
            mu_b = X_full_b[sub_idx].mean(axis=0)
            x_centered = X_full_b[i_glob] - mu_b
            XW += x_centered @ Ws[b]
        Z_full[i_glob] = XW @ WtW_inv

    weights = {}
    top_features = {}
    for b, name in enumerate(block_names):
        W = Ws[b]
        feats = block_feats[name]
        weights[name] = {feats[j]: W[j, :].tolist() for j in range(W.shape[0])}
        tpf = {}
        for f in range(W.shape[1]):
            abs_w = np.abs(W[:, f])
            order = np.argsort(abs_w)[::-1][:20]
            tpf[int(f)] = [(feats[j], float(abs_w[j])) for j in order]
        top_features[name] = tpf

    return {
        "n_factors": int(Z_full.shape[1]),
        "samples": list(common),
        "factor_scores": Z_full.tolist(),
        "factor_variance_explained": {
            name: np.zeros(Z_full.shape[1]).tolist() for name in block_names
        },
        "weights": weights,
        "top_features_per_factor": top_features,
        "_method": f"Subsample MOFA+ (N_sub={n_sub}, OLS projection for the rest)",
    }


# ---------------------------------------------------------------------------
# Cohort loaders (single-block + multi-block)
# ---------------------------------------------------------------------------

def _load_cptac(name):
    cohort_dir = REPO / "data" / "cohorts" / name
    sys.path.insert(0, str(cohort_dir))
    from importlib import import_module, reload
    mod_name = "load_and_map" if (cohort_dir / "load_and_map.py").exists() \
                else "run_gizmo_full"
    mod = import_module(mod_name)
    reload(mod)
    sys.path.pop(0)
    if name == "CPTAC_CCRCC":
        rna, prot, phos, _yl, samples = mod.load_ccrcc()
        return {"rna": rna, "prot": prot, "phospho": phos}, samples
    if name == "CPTAC_COAD":
        rna, prot, _yl, samples = mod.load_coad()
        return {"rna": rna, "prot": prot}, samples
    if name == "CPTAC_OV":
        rna, prot, _yl, samples = mod.load_ov()
        return {"rna": rna, "prot": prot}, samples
    raise ValueError(f"unknown CPTAC cohort {name}")


def _load_gse(cohort):
    from drug_sim_multi_cohort import load_sle_unified, load_sepsis_unified
    if cohort == "GSE65391_SLE":
        data, _y, common, _kind, _active = load_sle_unified()
    elif cohort == "GSE65682_sepsis":
        data, _y, common, _kind, _active = load_sepsis_unified()
    else:
        raise ValueError(f"unknown GSE cohort {cohort}")
    return {"rna": data}, common


COHORT_PLAN = {
    # cohort_name → (strategy_label, loader_fn, blocks_kind)
    "GSE65391_SLE":    ("incremental_pca",   lambda: _load_gse("GSE65391_SLE")),
    "GSE65682_sepsis": ("incremental_pca",   lambda: _load_gse("GSE65682_sepsis")),
    "CPTAC_CCRCC":     ("subsample_mofa",    lambda: _load_cptac("CPTAC_CCRCC")),
    "CPTAC_COAD":      ("subsample_mofa",    lambda: _load_cptac("CPTAC_COAD")),
    "CPTAC_OV":        ("subsample_mofa",    lambda: _load_cptac("CPTAC_OV")),
}


def main():
    # CLI: `python3 baseline_mofa_streaming.py [COHORT] [--substrate-matched]`
    argv = list(sys.argv[1:])
    substrate_matched = False
    if "--substrate-matched" in argv:
        argv.remove("--substrate-matched")
        substrate_matched = True
    target = argv[0] if argv else None
    plan = ({target: COHORT_PLAN[target]} if target in COHORT_PLAN
             else COHORT_PLAN)
    suffix = "_sm" if substrate_matched else ""
    for cohort, (strategy, loader) in plan.items():
        out_path = OUT_DIR / f"mofa_weights_{cohort}{suffix}.json"
        if out_path.exists():
            print(f"\n=== {cohort} ({strategy}{', substrate-matched' if substrate_matched else ''}) "
                  f"=== already cached; skipping", flush=True)
            continue
        print(f"\n=== {cohort} ({strategy}{', substrate-matched' if substrate_matched else ''}) ===",
              flush=True)
        block_dicts, common = loader()
        if not common:
            print("  no samples", flush=True); continue
        if strategy == "incremental_pca":
            (block_name, data) = next(iter(block_dicts.items()))
            result = fit_incremental_pca_single_block(
                data, common, n_factors=10, block_name=block_name,
                substrate_matched=substrate_matched)
        else:  # subsample_mofa
            result = fit_subsample_mofa_multiblock(
                block_dicts, common, n_factors=10, n_subsample=80,
                max_features_per_block=8000,
                substrate_matched=substrate_matched)
        if result is None:
            print("  no usable data", flush=True); continue
        out_path.write_text(json.dumps(result, indent=2))
        print(f"  n_factors: {result['n_factors']}, "
              f"samples: {len(result['samples'])}, "
              f"blocks: {list(result['weights'].keys())}", flush=True)
        print(f"  wrote {out_path}", flush=True)


if __name__ == "__main__":
    main()
