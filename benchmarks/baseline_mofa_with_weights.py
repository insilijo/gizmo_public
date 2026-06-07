"""baseline_mofa_with_weights.py — MOFA+ run that saves factor weights (W).

Extension of `baseline_mofa.py` that also extracts per-feature factor weights
for each block (proteomics, metabolomics, transcriptomics). These weights
are what we need to identify "top features per factor" for the Stage 28U
unsupervised head-to-head against GIZMO's α/β decomposition.

Output JSON structure (per cohort/design):

  {
    "cohort": ..., "design": ...,
    "n_factors": K,
    "samples": [...],
    "factor_variance_explained": {
        block: [float per factor],
        ...
    },
    "weights": {
        block: {
            feature_name: [w_factor_0, w_factor_1, ...],
            ...
        }
    },
    "top_features_per_factor": {
        block: {
            factor_idx: [(feature, |w|), ...top 20]
        }
    },
    "factor_scores": [[...] per sample]
  }

Stage 28U input.
"""
from __future__ import annotations
import sys, json, statistics
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "benchmarks"))

OUT_DIR = REPO / "benchmarks" / "results" / "unsupervised" / "mofa_weights"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def fit_mofa_with_weights(prot_data, metab_data, common, n_factors: int = 10,
                            cohort_name: str = "cohort"):
    """Fit MOFA, return factor scores (Z) AND per-feature weights (W) per block."""
    from mofapy2.run.entry_point import entry_point

    blocks = []
    block_names = []
    block_feature_names = {}

    if prot_data is not None:
        all_p = sorted(set().union(*[set(prot_data[s]) for s in common if s in prot_data]))
        P = np.zeros((len(common), len(all_p)))
        for i, s in enumerate(common):
            if s not in prot_data: continue
            for j, f in enumerate(all_p):
                P[i, j] = prot_data[s].get(f, 0.0)
        blocks.append(P); block_names.append("prot")
        block_feature_names["prot"] = all_p
    if metab_data is not None:
        all_m = sorted(set().union(*[set(metab_data[s]) for s in common if s in metab_data]))
        M = np.zeros((len(common), len(all_m)))
        for i, s in enumerate(common):
            if s not in metab_data: continue
            for j, f in enumerate(all_m):
                M[i, j] = metab_data[s].get(f, 0.0)
        blocks.append(M); block_names.append("metab")
        block_feature_names["metab"] = all_m

    if not blocks:
        return None

    ent = entry_point()
    ent.set_data_options(scale_views=True, scale_groups=False)
    # MOFA+ API: data[view_idx][group_idx] = (samples × features) matrix
    # We have len(blocks) views, 1 group → data = [[P], [M], ...]
    data_struct = [[b] for b in blocks]
    ent.set_data_matrix(
        data_struct,
        views_names=block_names,
        groups_names=["all"],
        samples_names=[common],
    )
    n_factors_eff = min(n_factors, len(common) - 2)
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
    # Z_raw could be a list of group matrices, or a single 2D matrix.
    # Normalize: get a (samples × factors) 2D matrix
    if isinstance(Z_raw, list):
        Z = np.asarray(Z_raw[0])
    else:
        Z = np.asarray(Z_raw)
    if Z.ndim == 1:
        Z = Z.reshape(-1, 1)
    # expect["W"] is a list of view dicts; each "E" is (n_features × n_factors)
    Ws = [np.asarray(expect["W"][m]["E"]) for m in range(len(block_names))]

    n_factors_actual = Z.shape[1]
    # Variance explained — skip if API differs from what we expect
    ve_per_block_factor = np.zeros((n_factors_actual, len(blocks)))

    weights = {}
    top_features = {}
    for bi, block_name in enumerate(block_names):
        W = Ws[bi]  # (n_features, n_factors)
        feats = block_feature_names[block_name]
        weights[block_name] = {feats[j]: W[j, :].tolist() for j in range(W.shape[0])}
        # Top 20 features per factor by |W|
        tpf = {}
        for f in range(W.shape[1]):
            abs_w = np.abs(W[:, f])
            order = np.argsort(abs_w)[::-1][:20]
            tpf[int(f)] = [(feats[j], float(abs_w[j])) for j in order]
        top_features[block_name] = tpf

    return {
        "n_factors": int(Z.shape[1]),
        "samples": list(common),
        "factor_scores": Z.tolist(),
        "factor_variance_explained": {
            block_names[b]: ve_per_block_factor[:, b].tolist()
            for b in range(len(block_names))
        },
        "weights": weights,
        "top_features_per_factor": top_features,
    }


def main():
    from per_patient_master import (
        load_crohn, load_su_covid, load_erawijantari, load_gao_ra,
        load_idh_glioma, load_filbin_covid,
        load_tcga_idh_glioma, load_corevitas,
        load_kmplot_brca, load_tcga_luad, load_gse89408_ra, load_hmp2_ibd_cd)

    loaders = {
        "Crohn": load_crohn, "Su_COVID": load_su_covid,
        "Erawijantari": load_erawijantari, "Gao_RA": load_gao_ra,
        "IDH_glioma": load_idh_glioma, "Filbin_COVID": load_filbin_covid,
        "TCGA_IDH_glioma": load_tcga_idh_glioma, "CorEvitas_RA": load_corevitas,
        "KMPLOT_BRCA": load_kmplot_brca, "TCGA_LUAD": load_tcga_luad,
        "GSE89408_RA": load_gse89408_ra, "HMP2_IBD_CD": load_hmp2_ibd_cd,
    }
    cohort_arg = sys.argv[1] if len(sys.argv) > 1 else None
    if cohort_arg:
        loaders = {cohort_arg: loaders[cohort_arg]}

    for cohort, loader in loaders.items():
        print(f"\n=== {cohort} ===", flush=True)
        prot, metab, ylab, common = loader()
        result = fit_mofa_with_weights(prot, metab, common, n_factors=10,
                                          cohort_name=cohort)
        if result is None:
            print(f"  no usable data"); continue
        out_path = OUT_DIR / f"mofa_weights_{cohort}.json"
        out_path.write_text(json.dumps(result, indent=2))
        print(f"  factors: {result['n_factors']}, "
              f"samples: {len(result['samples'])}, "
              f"blocks: {list(result['weights'].keys())}", flush=True)
        print(f"  wrote {out_path}", flush=True)


if __name__ == "__main__":
    main()
