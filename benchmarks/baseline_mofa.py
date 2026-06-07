"""Baseline: MOFA+ multi-block factor decomposition → SVC classifier.

For each cohort with paired multi-omics: train MOFA model on (samples ×
features) per omic, extract per-sample factor loadings, classify with
class-balanced SVC + 5-seed CV/LOO/hold-out.

Compares against per_patient_master.py — isolates "graph-aware"
contribution by replacing graph propagation with factor decomposition.

Output: benchmarks/results/baseline_mofa.tsv
"""
from __future__ import annotations
import sys, math, statistics
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "benchmarks" / "results"
sys.path.insert(0, str(REPO / "benchmarks"))


def fit_mofa(prot_data, metab_data, common, n_factors: int = 10):
    """Fit MOFA on paired multi-omic data; return sample × factor matrix.

    If only one omic provided, falls back to PCA-equivalent decomposition
    via the same factor model.
    """
    from mofapy2.run.entry_point import entry_point

    # Build per-omic matrices: features × samples
    blocks = []
    block_names = []
    if prot_data is not None:
        all_p = sorted(set().union(*[set(prot_data[s]) for s in common
                                       if s in prot_data]))
        P = np.zeros((len(common), len(all_p)))
        for i, s in enumerate(common):
            if s not in prot_data: continue
            for j, f in enumerate(all_p):
                P[i, j] = prot_data[s].get(f, 0.0)
        blocks.append(P); block_names.append("prot")
    if metab_data is not None:
        all_m = sorted(set().union(*[set(metab_data[s]) for s in common
                                       if s in metab_data]))
        M = np.zeros((len(common), len(all_m)))
        for i, s in enumerate(common):
            if s not in metab_data: continue
            for j, f in enumerate(all_m):
                M[i, j] = metab_data[s].get(f, 0.0)
        blocks.append(M); block_names.append("metab")

    if not blocks:
        return None

    # MOFA expects views as list of (group × samples × features) or similar
    # Use single-group setup
    ent = entry_point()
    ent.set_data_options(scale_views=True, scale_groups=False)
    ent.set_data_matrix(
        [[b for b in blocks]],  # one group, multiple views
        views_names=block_names,
        groups_names=["all"],
        samples_names=[common],
    )
    ent.set_model_options(factors=min(n_factors, len(common) - 2),
                           spikeslab_weights=True, ard_factors=True,
                           ard_weights=True)
    ent.set_train_options(iter=300, convergence_mode="medium",
                           dropR2=0.001, gpu_mode=False, seed=0,
                           verbose=False)
    ent.build()
    ent.run()

    # Extract per-sample factor scores: Z is samples × factors
    Z = ent.model.getExpectations()["Z"]["E"][0]
    return Z


def cv_classifier_svc(X: np.ndarray, y: np.ndarray, n_seeds: int = 5):
    from sklearn.svm import SVC
    from sklearn.model_selection import StratifiedKFold, train_test_split
    from sklearn.preprocessing import StandardScaler
    from sklearn.metrics import roc_auc_score
    cv_means, cv_sds, loo_aucs, holdout_aucs = [], [], [], []
    for seed in range(n_seeds):
        Xs = StandardScaler().fit_transform(X)
        idx_tr, idx_te = train_test_split(
            np.arange(len(y)), test_size=0.3, stratify=y, random_state=seed)
        clf = SVC(kernel="rbf", probability=True, C=1.0,
                   random_state=seed, class_weight="balanced")
        clf.fit(Xs[idx_tr], y[idx_tr])
        p = clf.predict_proba(Xs[idx_te])[:, 1]
        holdout_aucs.append(roc_auc_score(y[idx_te], p))
        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
        fold_aucs = []
        for tr, te in skf.split(Xs, y):
            c = SVC(kernel="rbf", probability=True, C=1.0,
                     random_state=seed, class_weight="balanced")
            try:
                c.fit(Xs[tr], y[tr])
                fold_aucs.append(roc_auc_score(y[te], c.predict_proba(Xs[te])[:, 1]))
            except Exception:
                pass
        cv_means.append(statistics.mean(fold_aucs) if fold_aucs else float("nan"))
        cv_sds.append(statistics.stdev(fold_aucs) if len(fold_aucs) > 1 else 0.0)
        loo_preds = np.zeros(len(y))
        for i in range(len(y)):
            tr = np.array([j for j in range(len(y)) if j != i])
            c = SVC(kernel="rbf", probability=True, C=1.0,
                     random_state=seed, class_weight="balanced")
            try:
                c.fit(Xs[tr], y[tr])
                loo_preds[i] = c.predict_proba(Xs[i:i+1])[0, 1]
            except Exception:
                pass
        try:
            loo_aucs.append(roc_auc_score(y, loo_preds))
        except Exception:
            loo_aucs.append(float("nan"))
    return {
        "cv_auc": statistics.mean(cv_means),
        "cv_sd": statistics.mean(cv_sds),
        "cv_seed_sd": statistics.stdev(cv_means) if len(cv_means) > 1 else 0.0,
        "loo_auc": statistics.mean(loo_aucs),
        "loo_seed_sd": statistics.stdev(loo_aucs) if len(loo_aucs) > 1 else 0.0,
        "holdout_auc": statistics.mean(holdout_aucs),
        "holdout_seed_sd": statistics.stdev(holdout_aucs) if len(holdout_aucs) > 1 else 0.0,
    }


def main():
    print("=" * 80)
    print("Baseline: MOFA+ multi-block factor decomposition → SVC")
    print("=" * 80)

    from per_patient_master import (
        load_crohn, load_su_covid, load_erawijantari, load_gao_ra,
    )
    cohorts = [
        ("Crohn",        load_crohn),
        ("Su_COVID",     load_su_covid),
        ("Erawijantari", load_erawijantari),
        ("Gao_RA",       load_gao_ra),
    ]

    rows = []
    for cohort_name, loader in cohorts:
        print(f"\n{'#' * 60}\n# {cohort_name} (MOFA)\n{'#' * 60}")
        try:
            prot_data, metab_data, y_label, common_all = loader()
        except Exception as exc:
            print(f"  load failed: {exc}")
            continue

        # Use all paired samples for MOFA (best signal)
        common = [s for s in common_all
                   if (prot_data and s in prot_data) or (metab_data and s in metab_data)]
        # For multi-block: need samples present in ALL omics
        if prot_data and metab_data:
            common_paired = [s for s in common if s in prot_data and s in metab_data]
            if common_paired:
                common = common_paired
        y = np.array([1 if y_label[s] == "active" else 0 for s in common])
        if len(set(y)) < 2 or len(common) < 10:
            print(f"  insufficient samples / class diversity; skipping")
            continue
        print(f"  n={len(common)} {Counter(y_label[s] for s in common)}", flush=True)

        try:
            Z = fit_mofa(prot_data, metab_data, common, n_factors=10)
            if Z is None:
                continue
            print(f"    MOFA factor scores: {Z.shape}", flush=True)
        except Exception as exc:
            print(f"    MOFA fit failed: {exc}", flush=True)
            continue

        res = cv_classifier_svc(Z, y, n_seeds=5)
        rows.append({"cohort": cohort_name, "method": "MOFA+SVC",
                       "n": len(common), **res})
        print(f"    CV={res['cv_auc']:.3f}±{res['cv_seed_sd']:.3f}  "
              f"LOO={res['loo_auc']:.3f}±{res['loo_seed_sd']:.3f}  "
              f"hold-out={res['holdout_auc']:.3f}±{res['holdout_seed_sd']:.3f}",
              flush=True)

    out_path = RESULTS / "baseline_mofa.tsv"
    cols = ["cohort", "method", "n",
            "cv_auc", "cv_sd", "cv_seed_sd",
            "loo_auc", "loo_seed_sd",
            "holdout_auc", "holdout_seed_sd"]
    out_path.write_text(
        "\t".join(cols) + "\n" +
        "\n".join("\t".join(str(r.get(c, "")) for c in cols) for r in rows) + "\n"
    )
    print(f"\nWrote: {out_path}")


if __name__ == "__main__":
    main()
