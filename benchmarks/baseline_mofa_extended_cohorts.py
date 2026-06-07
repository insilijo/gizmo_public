"""baseline_mofa_extended_cohorts.py — MOFA+ runs for the 5 panel cohorts
not covered by `baseline_mofa_with_weights.py`:

  CPTAC_CCRCC (3-block: RNA + prot + phospho)
  CPTAC_COAD  (2-block: RNA + prot)
  CPTAC_OV    (2-block: RNA + prot)
  GSE65391_SLE (1-block: microarray RNA)
  GSE65682_sepsis (1-block: microarray RNA)

Saves `mofa_weights_<cohort>.json` in the same schema as the original wrapper
so downstream tooling (`multi_pc_vs_mofa_factors.py`, the Bonferroni augmenter)
picks them up without modification.
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


def fit_mofa_arbitrary_blocks(block_dicts, common, n_factors: int = 10,
                                 max_features_per_block: int = 10000):
    """Fit MOFA+ on an arbitrary collection of {block_name: {sample_id:
    {feature: value}}} dicts. Returns the same JSON-serializable dict shape
    `fit_mofa_with_weights` returns.

    To fit in memory on cohorts with very wide feature universes
    (GSE65391_SLE: 29,737 microarray probes; GSE65682_sepsis: 19,000+
    probes), each block is variance-filtered to the top-K features by
    variance across the cohort. K=10000 is the MOFA+ paper's recommended
    default upper bound for a single block on a workstation; this still
    leaves the substrate-mappable features intact because they're
    consistently among the highest-variance.
    """
    from mofapy2.run.entry_point import entry_point

    block_names = []
    blocks = []
    block_feature_names = {}
    for name, data in block_dicts.items():
        if not data:
            continue
        all_feats = sorted({f for s in common if s in data
                            for f in data[s]})
        if not all_feats:
            continue
        X = np.zeros((len(common), len(all_feats)), dtype=np.float32)
        for i, sid in enumerate(common):
            if sid not in data: continue
            row = data[sid]
            for j, f in enumerate(all_feats):
                X[i, j] = row.get(f, 0.0)
        # Variance-filter to top-K features per block
        if X.shape[1] > max_features_per_block:
            var = np.nanvar(X, axis=0)
            keep_idx = np.argsort(var)[-max_features_per_block:]
            keep_idx.sort()
            X = X[:, keep_idx]
            all_feats = [all_feats[j] for j in keep_idx]
            print(f"    {name}: variance-filtered "
                  f"{len(var)} → {len(all_feats)} features (top-K by var)",
                  flush=True)
        blocks.append(X)
        block_names.append(name)
        block_feature_names[name] = all_feats

    if not blocks:
        return None

    ent = entry_point()
    ent.set_data_options(scale_views=True, scale_groups=False)
    ent.set_data_matrix(
        [[b] for b in blocks],
        views_names=block_names,
        groups_names=["all"],
        samples_names=[list(common)],
    )
    n_factors_eff = max(2, min(n_factors, len(common) - 2))
    ent.set_model_options(factors=n_factors_eff,
                           spikeslab_weights=True, ard_factors=True,
                           ard_weights=True)
    ent.set_train_options(iter=300, convergence_mode="medium",
                           dropR2=0.001, gpu_mode=False, seed=0,
                           verbose=False)
    ent.build()
    ent.run()

    expect = ent.model.getExpectations()
    Z_raw = expect["Z"]["E"]
    Z = np.asarray(Z_raw[0]) if isinstance(Z_raw, list) else np.asarray(Z_raw)
    if Z.ndim == 1:
        Z = Z.reshape(-1, 1)
    Ws = [np.asarray(expect["W"][m]["E"]) for m in range(len(block_names))]

    weights = {}
    top_features = {}
    for bi, name in enumerate(block_names):
        W = Ws[bi]
        feats = block_feature_names[name]
        weights[name] = {feats[j]: W[j, :].tolist() for j in range(W.shape[0])}
        tpf = {}
        for f in range(W.shape[1]):
            abs_w = np.abs(W[:, f])
            order = np.argsort(abs_w)[::-1][:20]
            tpf[int(f)] = [(feats[j], float(abs_w[j])) for j in order]
        top_features[name] = tpf

    ve_per_block_factor = np.zeros((Z.shape[1], len(blocks)))
    return {
        "n_factors": int(Z.shape[1]),
        "samples": list(common),
        "factor_scores": Z.tolist(),
        "factor_variance_explained": {
            name: ve_per_block_factor[:, b].tolist()
            for b, name in enumerate(block_names)
        },
        "weights": weights,
        "top_features_per_factor": top_features,
    }


# -------------------- cohort loaders -----------------------------------------

def _load_cptac(name):
    """CPTAC cohorts have their loaders in `run_gizmo_full.py` / `load_and_map.py`
    next to the per-cohort raw data. Path-insert that dir + import."""
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
    """RNA-only microarray cohorts (GSE65391 SLE, GSE65682 sepsis).
    Returns ({block_name: data}, list_of_common_samples)."""
    from drug_sim_multi_cohort import load_sle_unified, load_sepsis_unified
    if cohort == "GSE65391_SLE":
        data, y, common, _kind, _active = load_sle_unified()
    elif cohort == "GSE65682_sepsis":
        data, y, common, _kind, _active = load_sepsis_unified()
    else:
        raise ValueError(f"unknown GSE cohort {cohort}")
    return {"rna": data}, common


COHORT_LOADERS = {
    "CPTAC_CCRCC":   lambda: _load_cptac("CPTAC_CCRCC"),
    "CPTAC_COAD":    lambda: _load_cptac("CPTAC_COAD"),
    "CPTAC_OV":      lambda: _load_cptac("CPTAC_OV"),
    "GSE65391_SLE":  lambda: _load_gse("GSE65391_SLE"),
    "GSE65682_sepsis": lambda: _load_gse("GSE65682_sepsis"),
}


def main():
    target = sys.argv[1] if len(sys.argv) > 1 else None
    loaders = ({target: COHORT_LOADERS[target]}
               if target in COHORT_LOADERS else COHORT_LOADERS)
    for cohort, loader in loaders.items():
        out_path = OUT_DIR / f"mofa_weights_{cohort}.json"
        if out_path.exists():
            print(f"=== {cohort} === already cached at {out_path}; skipping",
                  flush=True)
            continue
        print(f"\n=== {cohort} ===", flush=True)
        block_dicts, common = loader()
        print(f"  {len(common)} samples; blocks: "
              f"{[(n, len(d)) for n, d in block_dicts.items()]}", flush=True)
        result = fit_mofa_arbitrary_blocks(block_dicts, common, n_factors=10)
        if result is None:
            print("  no usable data", flush=True); continue
        out_path.write_text(json.dumps(result, indent=2))
        print(f"  factors: {result['n_factors']}, "
              f"samples: {len(result['samples'])}, "
              f"blocks: {list(result['weights'].keys())}", flush=True)
        print(f"  wrote {out_path}", flush=True)


if __name__ == "__main__":
    main()
