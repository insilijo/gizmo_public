---
title: "v7 Supplement — extended methods, proofs, and reproducibility tables"
date: "2026-06-13"
status: "v7 supplement accompanying MANUSCRIPT_v7.docx"
---

# S1 Extended Methods

## S1.1 Substrate construction in detail

The substrate is the merge of four curated biological knowledge sources, deposited at `data/processed/human_full/graph.json` under CC-BY 4.0.

**Reactome** (release 87, June 2024): reactions filtered to `species == "Homo sapiens"`; protein participants resolved to UniProt accessions then mapped to HGNC gene symbols via the Reactome `Ensembl2Reactome_All_Levels.txt` table; small-molecule participants resolved to ChEBI IDs and matched to HMDB IDs via the ChEBI-to-HMDB cross-reference; reaction nodes preserve catalysis direction (substrate → reaction → product); regulator/modifier annotations are encoded as gene → reaction MODIFIER edges.

**StringDB** (v12.0): protein-protein interaction edges at confidence ≥ 700; multi-edges between gene-symbol pairs are reduced to single highest-confidence edge; STRING-only edges that already appear in Reactome (catalytic chain) are de-duplicated.

**HMDB** (v5.0): metabolite identifiers; manual de-duplication by InChIKey14; HMDB metabolites that don't map to any Reactome reaction are retained as isolated metabolite nodes; cross-database identifier collisions (HMDB → KEGG → ChEBI) resolved by precedence Reactome > HMDB > KEGG.

**KEGG** (latest snapshot, 2024-Q4): pathway annotations cross-referenced to gene + metabolite nodes; KEGG pathway hierarchy is stored as a node-level attribute but is *not* an edge type in the propagation-eligible substrate (KEGG pathway co-membership would introduce dense same-pathway smoothing that is biology-redundant with Reactome catalysis edges).

**Hub-capping:** nodes with degree > 200 in the merged graph are downgraded to peripheral status during MAP smoothing (their Laplacian rows are scaled down). This prevents the smoothing from being dominated by super-hub nodes (TP53, EGFR, top-50 inflammatory cytokines) whose connections to the rest of the graph would otherwise wash out all per-patient signal through over-smoothing.

**Final graph size:** 38,148 propagation-eligible nodes — 12,872 gene nodes (HGNC), 14,032 reaction nodes (Reactome reaction IDs), 11,244 metabolite nodes (HMDB IDs). 168,335 edges.

## S1.2 Per-modality observation covariance Σ specification

The MAP solve `||x - A_obs F||²_Σ⁻¹ + λ F^T L_signed F + ρ ||F||²` requires a per-modality observation covariance to weight data fidelity. We use a diagonal Σ with per-modality scalar variance σ²_m:

$$ \boldsymbol{\Sigma}^{-1} = \mathrm{diag}\left(\frac{1}{\sigma^2_m} \cdot \mathbb{1}_{[\text{node } j \in \text{modality } m]}\right) $$

The per-modality scalar σ²_m is set to the cohort-empirical residual variance after a one-pass MAP solve:

1. Initialize σ²_m = 1 for all modalities; solve MAP for all patients
2. Compute residual r_m = ||x_m - (A_obs F)_m||² / N_m per modality
3. Set σ²_m ← r_m / global_residual_geomean, normalizing so the geometric mean of all σ²_m is 1
4. Re-solve MAP with the calibrated Σ

This calibration handles the multi-modality scale problem (RNA-seq features at log-RSEM ~6, Olink NPX ~3, metabolomics intensity ~7, etc.) by reweighting each modality's data-fidelity contribution to its empirical noise floor rather than the raw input magnitude. A patient with a noisy single proteomic measurement contributes proportionally less to F than a patient with a clean RNA-seq panel.

For single-modality cohorts (KMPLOT-BRCA, TCGA-LUAD, TCGA-IDH-glioma, Filbin-Olink), Σ degenerates to σ² · I and the calibration step is trivial. For multi-modality cohorts (Su_COVID, Crohn, Filbin with metabolomics, IDH-glioma), the calibration matters and is documented per-cohort in `benchmarks/results/unsupervised/_pre_zscore_snapshot_20260607/per_cohort_sigma.json`.

## S1.3 Signed-Laplacian definiteness proof

The substrate's edges carry directed roles from Reactome: SUBSTRATE, PRODUCT, MODIFIER. The signed adjacency matrix A assigns nonzero entries by role:

$$
A_{ij} = \begin{cases} +w_{\text{prod}} \cdot s_{ij} & \text{if } (i,j) \text{ is PRODUCT (catalytic)} \\ -w_{\text{sub}} \cdot s_{ij} & \text{if } (i,j) \text{ is SUBSTRATE (substrate-of-reaction)} \\ +w_{\text{mod}} & \text{if } (i,j) \text{ is MODIFIER (gene → reaction)} \\ 0 & \text{otherwise} \end{cases}
$$

with stoichiometry s_ij and edge weights w_prod, w_sub, w_mod > 0. The matrix is symmetrized A ← (A + Aᵀ)/2.

Following Kunegis (2010), we use the **absolute-degree convention**: D_ii = Σ_j |A_ij|. The normalized signed Laplacian is:

$$ \mathbf{L}_{\text{signed}} = \mathbf{I} - \mathbf{D}^{-1/2} \mathbf{A} \mathbf{D}^{-1/2} $$

**Claim:** L_signed is positive semi-definite for any signed graph under this convention, regardless of structural balance.

**Proof.** Let x̃_i = x_i / √D_ii. Direct expansion:

$$ \mathbf{x}^T \mathbf{L}_{\text{signed}} \mathbf{x} = \mathbf{x}^T \mathbf{x} - \mathbf{x}^T \mathbf{D}^{-1/2} \mathbf{A} \mathbf{D}^{-1/2} \mathbf{x} $$

The cross-term unpacks as Σ_ij A_ij x̃_i x̃_j. Splitting A_ij = sign(A_ij) · |A_ij|:

$$ \sum_{ij} A_{ij} \tilde{x}_i \tilde{x}_j = \sum_{ij} \mathrm{sign}(A_{ij}) |A_{ij}| \tilde{x}_i \tilde{x}_j $$

Apply the identity 2 sign(a) ab = -((a - sign(a)b)² - a² - b²) when |sign|=1:

$$ \mathbf{x}^T \mathbf{L}_{\text{signed}} \mathbf{x} = \frac{1}{2} \sum_{ij} |A_{ij}| \left( \tilde{x}_i - \mathrm{sign}(A_{ij}) \tilde{x}_j \right)^2 $$

Each term is non-negative. Therefore xᵀ L_signed x ≥ 0 for all x.

Equality holds iff x̃_i = sign(A_ij) x̃_j on every nonzero edge, i.e., x is in the null space of L_signed. The null space is trivial except when the signed graph is structurally balanced (Harary 1953). For our substrate, balance fails because substrate-product edge pairs introduce sign-discordant cycles (a substrate is connected to its reaction by a negative edge and the product by a positive edge; the substrate-reaction-product chain has sign product (-1)(+1) = -1, while the substrate-PPI-product chain has sign (+1)). So L_signed has a non-zero spectral gap and trivial null space. □

**MAP system Hessian positive-definiteness.** The MAP minimizes:

$$ \mathcal{L}(\mathbf{F}) = (\mathbf{x} - \mathbf{A}_{\text{obs}} \mathbf{F})^T \boldsymbol{\Sigma}^{-1} (\mathbf{x} - \mathbf{A}_{\text{obs}} \mathbf{F}) + \lambda \mathbf{F}^T \mathbf{L}_{\text{signed}} \mathbf{F} + \rho \|\mathbf{F}\|^2 $$

Hessian H = Aₒbsᵀ Σ⁻¹ Aₒbs + λ L_signed + ρI. Each term is PSD: Aₒbsᵀ Σ⁻¹ Aₒbs because Σ⁻¹ is positive-diagonal; L_signed by the claim; ρI is positive-definite for any ρ > 0. The sum is **strictly positive definite**, the objective is **strictly convex**, and the MAP solution is **unique**.

Implementation: `gizmo/inference/laplacian.py` builds L_signed at lines 279–283; `gizmo/inference/projection.py` adds ρI = 10⁻³ to the system before CG.

## S1.4 λ-sensitivity (researcher degree of freedom)

The smoothing strength λ trades data fidelity against substrate-coherence. The default per-cohort selection is λ_diag = 0.1 × (mean data weight) / (median Laplacian eigenvalue × n_nodes).

| Cohort | λ_diag | β AUC range (λ ∈ [λ_diag/4, 4λ_diag]) | α-PC1 AUC range |
|---|---|---|---|
| Filbin_COVID | 8.2 × 10⁻⁴ | 0.638-0.651 (Δ 0.013) | 0.770-0.793 (Δ 0.023) |
| TCGA_IDH | 6.1 × 10⁻⁴ | (β not OS-predictive) | 0.812-0.834 (Δ 0.022) |
| TCGA_LUAD | 9.4 × 10⁻⁴ | 0.555-0.583 (Δ 0.028) | 0.582-0.612 (Δ 0.030) |
| Crohn | 4.7 × 10⁻⁴ | (n/a, no continuous outcome) | (interpretability stable) |

Median β AUC variation ≤ 0.03, median α-PC1 AUC variation ≤ 0.04 across the 16-fold λ range. **The choice of λ is not a major researcher degree of freedom.**

---

# S2 Reproducibility Tables

| Table | File | What it documents |
|---|---|---|
| S1 | `data/curation/v7_cohort_key_genes.tsv` | 75 curated source-paper key genes across 12 cohorts |
| S2 | `v7_interpretability_eval_v4.json` | 110 cell × cohort AUROC under degree-PR null |
| S3 | `v7_pathway_matched_null.tsv` | 110 cell × cohort AUROC under Reactome-pathway-degree null |
| S4 | `v7_basin_conservation.json` + `v7_basin_conservation_breakdown.tsv` | 149 basins + 36 cross-cohort pairs |
| S5 | `v7_basin_activation_matrix.tsv` | Disease × neighborhood activation matrix |
| S6 | `v7_smoothness_diagnostic.tsv` | s_cohort across 11 cohorts |
| S7 | `v7_toppr_comparison.tsv` | TopPR vs F-α-PC fold-enrichment per cohort × K |
| S8 | `v7_v5_F_per_pc_ablation.json` + `v7_pca_per_pc_ablation.json` + `v7_v5_canonical_F_all3.json` + `v7_survival_substrate_matched.json` | Per-PC ablation F-α-PC vs PCA-on-input; F-features ensemble + PCA-substrate-matched baselines per cohort (§4.2) |
| S9 | `v7_survival_filbin_28d.json` + `v7_basin_survival_broad.tsv` | Per-basin survival per cohort |
| S10 | `v7_f_as_feature_selection.json` | F-basin vs PCA-best-PC gene overlap |
| S11 | `stage31_v14c_wgcna_same_universe.json` + `stage31_v14b_wgcna_coherence.json` | WGCNA module comparison (corrected same-universe v14c + legacy v14b that the §4.1 self-correction refers to) |
| S12 | `v7_mito_2hg_degree_null.json` | 2HG anchor rank under canonical vs z-score |
| S13 | per-cohort λ sensitivity | Default + 4× scan |

---

# S3 Cited Software and Versions

- `gizmo` (the deposited package) v7.0; commit pinned per release
- `numpy` 1.26, `scipy` 1.11, `pandas` 2.1, `scikit-learn` 1.3, `networkx` 3.2
- `lifelines` 0.30.0 (Cox proportional-hazards, Harrell's C-index)
- `pandoc` 2.9.2.1 (markdown → docx conversion)
- `Python` 3.10

---

# S4 References

Kunegis et al. 2010. *Spectral analysis of signed graphs for clustering, prediction and visualization.* SIAM SDM.
Harary 1953. *On the notion of balance of a signed graph.* Michigan Math J.
Hofree et al. 2013. *Network-based stratification of tumor mutations.* Nature Methods.
Cowen et al. 2017. *Network propagation: a universal amplifier of genetic associations.* Nature Reviews Genetics.
Vanunu et al. 2010. *Associating genes and protein complexes with disease via network propagation.* PLoS Comput Biol.
Cantini et al. 2021. *Benchmarking joint multi-omics dimensionality reduction approaches.* Nature Communications.
Argelaguet et al. 2018. *MOFA: factor analysis for multi-omic data.* Mol Syst Biol.
Singh et al. 2019. *DIABLO: an integrative multi-omics approach for cancer subtyping.* Bioinformatics.
Wang et al. 2014. *Similarity Network Fusion.* Nature Methods.
Molenaar et al. 2014. *Radioprotection of IDH1-mutated cancer cells by the IDH1-mutant inhibitor AGI-5198.* Cancer Research.
Tateishi et al. 2015. *The Alkylating Chemotherapeutic Temozolomide Induces Metabolic Stress in IDH1-Mutant Cancers.* Cancer Research.
Cabrera-Benitez et al. 2014. *Mechanical Ventilation-associated Lung Fibrosis in ARDS: A Significant Contributor to Poor Outcome.* Anesthesiology.
Kratzer et al. 2008. *Bovine glomerular basement membrane.* Journal of Biological Chemistry.
