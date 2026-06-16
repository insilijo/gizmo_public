# v7 manuscript reproducibility — script index

19 scripts produce all numerical claims and figures in MANUSCRIPT_v7.docx (8,597 words, final v7 docx after R1/R2/R3/R4 hostile-review + v8-substantive corrections).

## Headline result scripts (§3 — interpretability)

| Script | Manuscript section | Output |
|---|---|---|
| `check_v7_interpretability_eval.py` | §3.3 — 24/110 cells, 10/11 cohorts pass degree-null | `v7_interpretability_eval_v4.json` |
| `check_v7_pathway_matched_null.py` | §3.3 — pathway-matched null (rigor fix) | `v7_pathway_matched_null.tsv` |
| `check_v7_hmp2_pathway_null_patch.py` | §3.3 — HMP2 fix (was silently skipped); 7/11 cohorts | merged into `v7_pathway_matched_null.tsv` |
| `check_v7_anchor_recovery_degree_null.py` | §3.4 utility | `v7_anchor_recovery_degree_null.json` |
| `check_v7_mito_2hg_degree_null.py` | §3.4 + §5.2 — 2HG anchor rank under z-score | `v7_mito_2hg_degree_null.json` |

## Method comparison scripts (§4.1 — apples-to-apples)

| Script | Manuscript section | Output |
|---|---|---|
| `check_v7_pca_per_pc_ablation.py` | §4.2 PCA per-PC baseline (for Fig 6) | `v7_pca_per_pc_ablation.json` |
| `check_v7_v5_F_per_pc_ablation.py` | §4.2 F-α-PC per-PC (for Fig 6) | `v7_v5_F_per_pc_ablation.json` |
| `check_v7_v5_canonical_F_all3.py` | §4.2 F-features under v5-canonical preprocessing | `v7_v5_canonical_F_all3.json` |
| `check_v7_survival_substrate_matched_pca.py` | §4.2 substrate-matched PCA baseline | `v7_survival_substrate_matched.json` |
| `check_v7_toppr_comparison.py` | §4.1 TopPR vs F-α-PC — overturns v6 'dominates' claim | `v7_toppr_comparison.tsv` |

## Per-basin biology scripts (§4.3, 4.4, 4.5)

| Script | Manuscript section | Output |
|---|---|---|
| `check_v7_basin_conservation.py` | §4.5 149 basins, 36 conserved pairs | `v7_basin_conservation.json` + `v7_basin_conservation_breakdown.tsv` |
| `check_v7_basin_threshold_sensitivity.py` | §2.3 5% threshold sits in conservation plateau | inline reporting |
| `check_v7_survival_filbin_covid.py` | §4.3 Filbin AUC 0.776 case study | `v7_survival_filbin_28d.json` |
| `check_v7_f_as_feature_selection.py` | §4.4 F-basin vs PCA orthogonal panels | `v7_f_as_feature_selection.json` |

## β/α decomposition + scope (§2.3, §5.1)

| Script | Manuscript section | Output |
|---|---|---|
| `check_v7_beta_alpha_ablation.py` | §2.3 β captures < 2% F variance; α-PCs ≈ PCA-on-F | inline reporting |
| `check_v7_beta_discrimination_test.py` | §2.3 F-features vs PCA-on-F head-to-head | inline reporting |
| `check_v7_beta_severity_panel.py` | §2.3 β-vs-severity exploratory (5 cohorts) | inline reporting |
| `check_v7_smoothness_statistic.py` | §5.1 s_cohort table; honest 1-dim limitation | `v7_smoothness_diagnostic.tsv` |

## WGCNA legacy script (§4.1)

| Script | Manuscript section | Output |
|---|---|---|
| `unsupervised_stratification/legacy/stage31_v14c_wgcna_same_universe.py` | §4.1 corrected same-gene-universe (n = 7 cohorts filtered): median tpf ratio 1.02 — at parity | `stage31_v14c_wgcna_same_universe.json` |

## Figure generation

| Script | Output | Coverage |
|---|---|---|
| `check_v7_make_figures.py` | `figures_v7/figure{1..11}.png` | 11 main + 1 supp (modality confound test) |

## v8 still-deferred

Items requiring data we don't have locally:
- Held-out z-score validation cohort (CGGA / GSE16011 / CPTAC_GBM + cBioPortal) — needs network
- Su_COVID β-severity ID mapping — F pids are numeric, clinical INCOV001-1
- HMP2/Crohn disease activity scores — zenodo metadata is minimal (patient_id + sex only)

## v6 → v7 ablations preserved

The pre-existing v14b WGCNA comparison (`stage31_v14b_wgcna_coherence.json`) reported the v6 narrative that this manuscript now corrects (per-module coherence is at parity, not GIZMO-favored).

## Running full reproduction

```bash
cd benchmarks
# In dependency order
for s in check_v7_anchor_recovery_degree_null.py \
         check_v7_interpretability_eval.py \
         check_v7_pathway_matched_null.py \
         check_v7_hmp2_pathway_null_patch.py \
         check_v7_mito_2hg_degree_null.py \
         check_v7_basin_conservation.py \
         check_v7_basin_threshold_sensitivity.py \
         check_v7_pca_per_pc_ablation.py \
         check_v7_v5_F_per_pc_ablation.py \
         check_v7_v5_canonical_F_all3.py \
         check_v7_survival_substrate_matched_pca.py \
         check_v7_survival_filbin_covid.py \
         check_v7_f_as_feature_selection.py \
         check_v7_beta_alpha_ablation.py \
         check_v7_beta_discrimination_test.py \
         check_v7_beta_severity_panel.py \
         check_v7_smoothness_statistic.py \
         check_v7_toppr_comparison.py \
         check_v7_make_figures.py; do
  PYTHONPATH=/home/jgardner/GIZMO /usr/bin/python3 "$s"
done

# WGCNA legacy (separate path)
PYTHONPATH=/home/jgardner/GIZMO /usr/bin/python3 \
  unsupervised_stratification/legacy/stage31_v14c_wgcna_same_universe.py
```

Estimated total compute: 5-6 hours (dominated by interpretability eval, pathway-matched null, and WGCNA v14c).
