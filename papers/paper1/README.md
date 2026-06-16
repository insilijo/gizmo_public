# Paper 1 — A biochemistry substrate as a fixed coordinate system for per-patient multi-omic state

This folder is self-contained. Everything cited by the manuscript lives here: text, figures, scripts, results, and the curated ground-truth file. Substrate (`../../substrate/graph.json`) and framework code (`../../gizmo/`) are shared across all papers and live at the repo root.

## Contents

- **[MANUSCRIPT.md](MANUSCRIPT.md)** / **[MANUSCRIPT.docx](MANUSCRIPT.docx)** — the manuscript (v7; 8,955 words, 11 embedded figures, accepted after five rounds of hostile peer review).
- **[SUPPLEMENT.md](SUPPLEMENT.md)** / **[SUPPLEMENT.docx](SUPPLEMENT.docx)** — extended methods: signed-Laplacian PSD proof, per-modality Σ specification, full substrate construction details, λ-sensitivity, 13 reproducibility-table pointers.
- **[FIGURES.md](FIGURES.md)** / **[FIGURES.docx](FIGURES.docx)** — figure captions + per-figure reproducibility table.
- **`figures/`** — 11 main figures (`figure1.png` … `figure11.png`).
- **`scripts/`** — 19 v7 manuscript-load-bearing analysis scripts + WGCNA v14b/v14c (via symlinks back to `../../benchmarks/unsupervised_stratification/legacy/`). See [`scripts/README.md`](scripts/README.md) for script-to-section mapping.
- **`results/`** — 20 result deposits (JSONs + TSVs) backing every manuscript number.
- **`data/curation/v7_cohort_key_genes.tsv`** — 75-entry curated source-paper key-gene file (blinded ground truth for the interpretability evaluation in §3).

## Headline results

- **§3.3 interpretability:** 10 of 11 cohorts pass at p<0.05 under degree-PageRank-preserving null; **7 of 11 survive the stricter Reactome-pathway-membership-preserving null** that controls for annotation density (Crohn thiopurine AUROC 0.911, HMP2 autophagy 0.834, IDH-mut catalysis 0.839, Filbin COVID 0.814, Su COVID 0.868, TCGA_LUAD 0.823, TCGA_IDH 0.730).
- **§4.3 per-basin mechanism-named survival:** Filbin α-PC1- proteoglycan-degradation basin reaches AUC 0.776 (vs PCA-on-Olink 0.717) as a named-mechanism case study. Ribosome basin C 0.480 = chance, correctly identifying housekeeping.
- **§4.4 F-vs-PCA orthogonality:** F's best prognostic basin shares 0–2/68 genes with PCA's best prognostic component across three cohorts — comparable discrimination, different named axes.
- **§4.5 cross-cohort conservation:** 24 conserved cross-disease basin pairs at Jaccard ≥ 0.30 (100% cross-disease after excluding interpretively-unvalidated Gao_RA basins).

## Three v6 claims overturned in v7

Re-running the cited ablations during hostile-review rounds R4/v8 corrected three v6 claims:

1. **TopPR doesn't dominate at every K** on v7 source-paper-curated truth (vs v6's CTD-anchor truth). F-α-PC median fold = 8.17 at K=200/500 vs TopPR 0/5.45 (`results/v7_toppr_comparison.tsv`).
2. **WGCNA per-module pathway coherence is at parity with GIZMO** under matched gene universe (median ratio 1.02 at proper filter, n=7 cohorts; `results/stage31_v14c_wgcna_same_universe.json`). The v6 "more biochemistry-coherent" claim is retracted.
3. **"GIZMO modules span 8–40× more Reactome pathways per module"** was a granularity-vs-coherence misframing. The actual difference is granularity (GIZMO 164 modules per cohort vs WGCNA 6, at parity per-module coherence).

## Reproducing

```bash
# From the repo root
cd /path/to/gizmo_public

# Each script can be run with PYTHONPATH at repo root
PYTHONPATH=. python3 papers/paper1/scripts/check_v7_interpretability_eval.py
PYTHONPATH=. python3 papers/paper1/scripts/check_v7_pathway_matched_null.py
PYTHONPATH=. python3 papers/paper1/scripts/check_v7_basin_conservation.py
# ... full mapping in scripts/README.md
```

Most scripts read from the deposited F matrices on Zenodo (DOI TBD on acceptance) plus the substrate at `../../substrate/graph.json`. Result deposits in `results/` are pre-computed; the scripts regenerate them.

Estimated total compute for full reproduction: 5–6 hours (dominated by `check_v7_interpretability_eval.py`, `check_v7_pathway_matched_null.py`, and the WGCNA v14c symlink).
