---
title: "v7 Figure organization and captions"
date: "2026-06-11"
status: "Phase 6 — figure deliverables for v7 manuscript (11 figures)"
---

# v7 Figure index

All figure source code at `/tmp/v7_figures.py`, `/tmp/v7_more_figures.py`, `/tmp/v7_network_figures.py`, `/tmp/v7_chord_figures.py`; assets at `/home/jgardner/GIZMO/figures_v7/`.

| # | File | Type | Cited in | Source data |
|---|---|---|---|---|
| 1 | figure1_pipeline_schematic.png | Diagram (4-panel) | §2 | substrate.json + conceptual |
| 2 | figure2_interpretability.png | Heatmap | §3.3 | `v7_interpretability_eval_v4.json` |
| 3 | figure3_activation_matrix.png | Heatmap | §4.5 | `v7_basin_activation_matrix.tsv` |
| 4 | figure4_basin_survival.png | Bar charts (4-panel) | §4.3 | `v7_basin_survival_broad.tsv` |
| 5 | figure5_f_pca_orthogonal.png | Bar charts (3-panel) | §4.4 | inline gene-overlap data |
| 6 | figure6_f_vs_pca_per_pc.png | Line plot (4-panel) | §4.2 | `v7_v5_F_per_pc_ablation.json` + `v7_pca_per_pc_ablation.json` |
| 7 | figure7_preprocessing_scope.png | Bar charts (2-panel) | §5.1 | `v7_mito_2hg_degree_null.json` + survival results |
| 8 | figure8_oxphos_basin_network.png | Network graph (2-panel) | §4.5 + §3.5 | `v7_basin_conservation.json` + substrate |
| 9 | figure9_filbin_ecm_basin_network.png | Network graph | §4.3 | `v7_basin_conservation.json` + substrate |
| 10 | figure10_chord_conservation.png | Chord diagram | §4.5 | `v7_basin_conservation.json` |
| 11 | figure11_patient_f_space.png | Scatter plot | §2.4 | 11 cohort F matrices |

**Visual variety**: 2 heatmaps, 5 bar-chart panels, 1 line plot, 3 network sub-graphs, 1 chord diagram, 1 patient-space scatter, 1 schematic.

---

# Figure captions

## Figure 1 — GIZMO pipeline: substrate construction, MAP projection, β/α decomposition, signed-basin extraction

`figure1_pipeline_schematic.png` (18 × 5 inches, 150 dpi)

Four-panel schematic.

- **(A)** Substrate sub-graph showing gene (blue), reaction (orange), and metabolite (green) nodes from a representative connected component of the 38,148-node Reactome + StringDB + HMDB + KEGG merged graph (CC-BY 4.0).
- **(B)** MAP solve equation: data-fidelity + signed-Laplacian smoothness + ridge regularizer; strictly convex with unique F per patient.
- **(C)** β/α decomposition: F_p = β_p × log_PR + α_p, where log_PR is the substrate-fixed log-PageRank direction and α is the orthogonal residual.
- **(D)** Signed-basin extraction: largest connected sub-graph of same-sign top-5% loading nodes; red = positive sign, blue = negative.

## Figure 2 — Quantitative interpretability evaluation across 11 cohorts

`figure2_interpretability.png` (11 × 6 inches, 150 dpi)

AUROC heatmap for 11 cohorts × 5 α-PCs × 2 signs = 110 cells. AUROC of (curated source-paper key-gene membership) vs (signed |loading|), tested against degree- and PageRank-preserving 1000-resample null. **24 cells significant at p < 0.05 (\*) including 17 at p < 0.01 (\*\*).** **10 of 11 cohorts pass at least one cell.** Best: Crohn α-PC1+ (AUROC 0.980, p < 0.001 — thiopurine metabolism MPG/APEX1/TPMT/NUDT15); Su_COVID α-PC3+ (0.919 \*\*); TCGA_IDH_glioma α-PC4- (0.874 \*\* — IDH-mut catalysis axis IDH1/2/D2HGDH/L2HGDH). Color scale RdBu_r with 0.50 as midpoint. Gao_RA fails the panel (AUROC 0.406 below chance — serum proteomic panel doesn't form coherent substrate sub-graph).

## Figure 3 — Disease × biochemical-neighborhood activation matrix

`figure3_activation_matrix.png` (14 × 9 inches, 150 dpi)

Heatmap of basin gene-count per disease cohort × biochemical neighborhood category. 11 cohorts (3 RA, 2 IBD, 2 COVID, 1 BRCA, 1 LUAD, 2 Glioma) × ~14 neighborhood categories. **149 total basins; 36 conserved cross-cohort at Jaccard ≥ 0.30.** Cell annotations: best α-PC × sign carrying each basin per cohort. Conservation graded by biology: housekeeping (ribosome 0.83–0.94 cross-cohort) > canonical disease programs (T cell/MHC 0.28–0.41) > disease-specific subspecialty (Filbin ECM-degradation vs LUAD ECM-structural share 0/19 genes despite both being "Collagen/ECM").

## Figure 4 — Per-basin patient activation scores → survival/mortality discrimination

`figure4_basin_survival.png` (13 × 10 inches, 150 dpi)

Four-panel horizontal bar charts: per-basin Cox C-index (TCGA_LUAD, TCGA_IDH, KMPLOT_BRCA) or AUC (Filbin_COVID 28-day mortality), with F-features 7-ensemble and PCA-on-input reference baselines (dashed lines). **Filbin_COVID PC1-** (proteoglycan + cathepsin + MMP; 7 observed genes) → **AUC 0.776 exceeds PCA-on-Olink 0.717 by Δ +0.06** — only cohort where F-basin BEATS PCA discrimination AND names the mechanism (lung interstitial ECM degradation by neutrophil enzymes). TCGA_LUAD T cell/MHC + complement basin → C 0.608 (beats F-7 ensemble 0.580, ties PCA-on-input 0.599). Ribosome basin (TCGA_LUAD PC2+) → C 0.480 (chance, correctly identifies housekeeping).

## Figure 5 — F-basins and PCA-best-PC select orthogonal gene panels

`figure5_f_pca_orthogonal.png` (15 × 5.5 inches, 150 dpi)

Three-panel bar chart: count of genes in (F-basin only / Shared / PCA-PC top-50 only) for each of TCGA_LUAD, TCGA_IDH, Filbin_COVID. **All three cohorts: ≤ 2 shared genes between F's best prognostic basin and PCA's best prognostic component despite comparable discrimination.** TCGA_LUAD F-basin = complement immune infiltrate (68g); PCA-PC1 = lung surfactant SFTPC/SCGB1A1 (50g); 2 shared (C6, C7). TCGA_IDH F-basin = ciliary dynein DNAAF (7g); PCA-PC1 = cellular activation CLIC1/S100A11; 0 shared. Filbin F-basin = ECM degradation; PCA-PC1 = organ-damage acute-phase NTproBNP/FGF23/FABP1; 0 shared. **F and PCA read orthogonal prognostic biology.**

## Figure 6 — Per-α-PC F vs PCA-on-input ablation under matched preprocessing

`figure6_f_vs_pca_per_pc.png` (15 × 4.5 inches, 150 dpi)

Four-panel line plot: Cox C-index (or AUC for Filbin) per PC index for F-α-PCs (red) vs PCA-on-input PCs (blue) across 4 survival cohorts. F-α-PCs match or exceed PCA-PCs at the best-single-PC level in all 4: KMPLOT Δ +0.001, TCGA_LUAD Δ +0.033, **TCGA_IDH +0.138** (substrate-projection on IDH-mut catalysis axis), Filbin -0.002. Each panel annotates the best F-PC and best PCA-PC. Substrate-projection's largest discrimination gain appears where the cohort's training label aligns with graph-coherent biology.

## Figure 7 — Scope conditions: preprocessing choice is task-specific

`figure7_preprocessing_scope.png` (14 × 5 inches, 150 dpi)

Two-panel demonstration. **(A)** Mitochondrial-2HG anchor rank on TCGA_IDH glioma under canonical preprocessing (rank 6,355, out of top 5%) vs within-patient z-score (rank 26, well within top 5%). **(B)** Cox C-index for F-features under z-score vs v5-canonical preprocessing vs PCA-on-substrate-matched-input across KMPLOT_BRCA, TCGA_LUAD, TCGA_IDH. **Canonical preprocessing preserves inter-patient magnitudes required for magnitude-driven phenotypes (overall survival); z-score destroys magnitudes but recovers within-patient relative-rank signal for subtype/anchor recovery (§3.4).**

## Figure 8 — Mitochondrial OXPHOS basin: same biology surfaces as substrate sub-graph in unrelated cohorts

`figure8_oxphos_basin_network.png` (15 × 8 inches, 150 dpi)

Two-panel network sub-graph visualization. TCGA_IDH α-PC3+ (16 gene + 28 reaction nodes; e.g., MT-ND1-6, MT-CO1-3, MT-CYB, COA3, NDUFAF5, PET100, TIMMDC1) and KMPLOT_BRCA α-PC1- (19 gene + 17 reaction nodes; NDUFA10/11, NDUFAF2/5, MT-CO/ND core, MT-CYB) are different cohorts with completely different diseases (high-grade glioma vs ER+ breast cancer) but the same biochemical neighborhood emerges as the dominant α-PC basin. Genes (blue), reactions (orange), metabolites (green). Edges = substrate edges. Cross-cohort Jaccard 0.46.

## Figure 9 — Filbin_COVID ECM-degradation basin: substrate sub-graph naming a mechanism

`figure9_filbin_ecm_basin_network.png` (11 × 9 inches, 150 dpi)

Network visualization of Filbin_COVID α-PC1- basin (19 nodes; 7 genes mappable to Olink). Proteoglycan core (BGN, DCN, ACAN, VCAN, LUM, FMOD, KERA, CSPG4/5) connected to cathepsins (CTSL, CTSK) and matrix metalloproteinases (MMP13, MMP20) via substrate edges. **The basin is the SUBSTRATE-COHERENT mechanism — neutrophil-derived enzymes degrading lung interstitial proteoglycan matrix — that achieves AUC 0.776 for 28-day mortality (exceeds PCA-on-Olink 0.717 by Δ +0.06).** The framework recovers ARDS-pathology biology unsupervised, as a connected substrate sub-graph, named.

## Figure 10 — Cross-cohort basin conservation as a substrate-fixed activation vocabulary

`figure10_chord_conservation.png` (11 × 11 inches, 150 dpi)

Circular chord diagram. 11 cohorts as nodes on a circle, colored by disease group (RA = pink, IBD = orange, COVID = blue, breast = green, lung = red, glioma = purple). **36 conserved basin pairs as chord edges (curves), colored by biology category** (ribosome blue, OXPHOS red, T cell/MHC orange, DNA damage purple, collagen green, complement pink, cilia brown). Line width ∝ Jaccard similarity. **Shows: same biochemical neighborhoods (ribosome, OXPHOS, T cell/MHC, DNA damage) activate across cohorts of unrelated diseases.** The substrate's coordinate system is a shared vocabulary.

## Figure 11 — All 2,732 patients from 11 cohorts in shared F-space coordinates

`figure11_patient_f_space.png` (11 × 9 inches, 150 dpi)

PCA-2 of L2-normalized patient F-vectors across all 11 deposited cohorts (n = 2,732 patients). Each cohort is a distinct cluster colored by disease group. Top-2 PCs hold 38% of cross-cohort variance. **Substrate-coordinate disease separation is visible in just 2 dimensions** — F-space coordinates carry per-cohort identity that is comparable across cohorts because the basis is shared. This is the patient-side counterpart of Figure 10's basin-side cross-cohort comparability.

---

# Supplementary figures (deferred for v8)

- **Supp Fig 1**: Phase 2 degree-preserving null sensitivity — null AUROC histograms per cohort × α-PC vs observed AUROC at p = 0.05, 0.01 thresholds
- **Supp Fig 2**: Per-cohort source-paper key-gene curation table (75 entries; cohort, gene, rationale, confidence)
- **Supp Fig 3**: Phase 4 pre-registration timestamps + locked prediction
- **Supp Fig 4**: Reproducibility table — file path → JSON deposit → figure → claim

---

# Reproducibility table (compressed for v7 docx)

| Claim | Source code | Result deposit | Figure |
|---|---|---|---|
| 10/11 cohorts pass degree-null p<0.05; 7/11 pass pathway-null | `check_v7_interpretability_eval.py` + `check_v7_pathway_matched_null.py` + `check_v7_hmp2_pathway_null_patch.py` | `v7_interpretability_eval_v4.json`, `v7_pathway_matched_null.tsv` | Fig 2 |
| 149 basins, 36 conserved at Jaccard ≥ 0.30 (33 cross-disease, 3 within) | `check_v7_basin_conservation.py` | `v7_basin_conservation.json` + `v7_basin_conservation_breakdown.tsv` | Fig 3, 10 |
| Per-basin survival (Filbin AUC 0.776 > PCA 0.717) | `check_v7_survival_filbin_covid.py` + inline | `v7_survival_filbin_28d.json` + `v7_basin_survival_broad.tsv` | Fig 4, 9 |
| F-basin vs PCA gene overlap 0–2/68 | `check_v7_f_as_feature_selection.py` + inline | `v7_f_as_feature_selection.json` | Fig 5 |
| TCGA_IDH α-PC4 best-single-PC C 0.762 > PCA-PC1 0.624 | `check_v7_v5_F_per_pc_ablation.py` | `v7_v5_F_per_pc_ablation.json` | Fig 6 |
| Z-score 2HG anchor rank 26 (p_degree 0.011) | `check_v7_mito_2hg_degree_null.py` | `v7_mito_2hg_degree_null.json` | Fig 7 §A |
| MAP solver signed-Laplacian PSD proof | `MANUSCRIPT_v7_METHODS_RIGOR.md` | — | §7 Methods |
| 2,732 patients across 11 cohorts in shared F-space | `check_v7_make_figures.py` from 11 F matrices | `_pre_zscore_snapshot_20260607/stage3_F_*.npz` | Fig 11 |
| **β captures <2% F variance; α-PC ≈ PCA-on-F** | `check_v7_beta_alpha_ablation.py` | inline | §2.3 |
| **F-features vs PCA-on-F discrimination head-to-head** | `check_v7_beta_discrimination_test.py` | inline | §2.3 |
| **β-severity panel across 5 cohorts (Filbin/IDH/BRCA/LUAD)** | `check_v7_beta_severity_panel.py` | inline | §2.3 |
| **5% basin threshold sensitivity (conservation plateau)** | `check_v7_basin_threshold_sensitivity.py` | inline | §2.3 |
| **s_cohort smoothness statistic across 11 cohorts** | `check_v7_smoothness_statistic.py` | `v7_smoothness_diagnostic.tsv` | §5.1 |
| **TopPR vs F-α-PC on v7 substrate (overturns v6 'dominates everywhere')** | `check_v7_toppr_comparison.py` | `v7_toppr_comparison.tsv` | §4.1 |
| **WGCNA v14c (corrected same-universe; median ratio 1.02 at proper filter)** | `unsupervised_stratification/legacy/stage31_v14c_wgcna_same_universe.py` | `stage31_v14c_wgcna_same_universe.json` | §4.1 |
