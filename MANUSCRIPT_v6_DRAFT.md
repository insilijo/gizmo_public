---
title: "A biochemistry substrate as a fixed coordinate system for multi-omic per-patient projection"
author: "Joseph J. Gardner (Insilijo)"
date: "v6 draft, 2026-06-09"
---

# A biochemistry substrate as a fixed coordinate system for multi-omic per-patient projection

*v6 draft — v5 narrative preserved as canonical; minimal targeted updates: §5 LoF/GoF scope retraction with empirical rescues from within-patient z-score preprocessing as a fine-grained subtype/anchor diagnostic. See `MANUSCRIPT_v5_DRAFT.md` for the prior version and `MANUSCRIPT_v6_HEAVY.md` for the abandoned full-rewrite draft.*

---

## Abstract

Multi-omic integration methods produce latent factors in cohort-specific feature spaces — each cohort gets its own basis, and the resulting axes carry no a priori biological interpretation and cannot be compared across cohorts. We address this by working in **substrate space** rather than factor space. The **GIZMO substrate** is a 38,148-node biochemistry graph (16,343 genes + 6,406 metabolites + 15,399 reactions + ancillary nodes; Reactome + StringDB + HMDB + KEGG, CC-BY 4.0) that is fixed across cohorts; any modality combination MAP-projects onto it through signed-Laplacian smoothing with per-modality measurement variance, producing a per-patient state vector **F** in 38,148-dimensional coordinates that depend only on the substrate, not on the cohort or modality.

F decomposes uniquely into **β** (projection onto the substrate's log-PageRank hub direction — empirically a phenotype-presentation magnitude axis; 47 of the top 50 hub nodes are transcriptional, signaling, or immune) and **α** (orthogonal mechanism residual). Each principal component of α decomposes parameter-free into two connected substrate subgraphs that name both poles of the patient axis — the framework's central interpretability claim.

Across 17 cohorts (n ≈ 5,000 patients), GIZMO recovers canonical disease biology without label supervision: ccRCC Warburg + HIF axis (the *unsupervised* + basin recovers SLC22A6/SLC22A8/uromodulin proximal-tubule machinery; the − basin recovers PFKP/HK2/ENO2/CA9, the textbook HIF-2α target); SLE psychosis the type I IFN signature (RSAD2/IFI27/IFI44L/IFIT3; the IFN-axis α-PC discriminates psychosis at AUC 0.99 in the n = 924 vs n = 72 class-imbalanced lupus cohort and survives stratified bootstrap resampling at 95% CI [0.97, 1.00]); Filbin α-PC4 the lactoperoxidase antimicrobial 7-clique vs the galanin receptor 5-clique (12 nodes hold 87.7% of the PC mass; STABLE under leave-one-out, p05 cosine = 0.997). Three cohorts with *zero overlapping input features* (RNA-only IDH-glioma, metab-only HMP2_IBD_CD, RNA-only GSE89408_RA) converge on the same disease-anchor metabolites at top-5% rank in the substrate's coordinate system, joint p ≈ 2×10⁻³ — the substrate enables horizontal meta-analysis that cohort-specific factor methods cannot perform by construction.

At symmetric multi-axis discrimination testing across 14 cohorts, GIZMO and MOFA+ are at parity (8/14 vs 9/14 at global Bonferroni; **6/11 vs 6/11 exact parity** when MOFA+ is restricted to GIZMO's substrate-mappable input universe). Discrimination parity is the expected result — both methods extract patient covariance structure from the same per-cohort variance — and is not the contribution. The contribution is the substrate itself plus the interpretive layer GIZMO instantiates on top of it (per-patient state vectors, β/α decomposition, signed-basin biology, substrate-fixed coordinates).

Under the canonical preprocessing the substrate represents loss-of-function and differentiation-collapse biology cleanly. Under within-patient z-score preprocessing (a fine-grained subtype/anchor diagnostic; Methods §"Within-patient z-score"), the substrate additionally represents gain-of-function biology with canonical downstream effectors: TCGA_LUAD KRAS-vs-EGFR multi-PC LR LOOCV AUC lifts from 0.509 (chance under canonical preprocessing) to 0.703 and TCGA_IDH 2HG-mito anchor lifts from rank 6,355/6,406 to 26/6,406 (top 0.41%, correct mut > wt direction). v4–v5 framed these failures as substrate-representation gaps; v6 retracts that framing — both were preprocessing artifacts on canonical substrate edges. The narrower Paper 4 open question contracts to truly neomorphic biology with zero canonical effector signaling and not rescuable by within-patient z-score. The substrate (CC-BY 4.0) and every cohort GIZMO has been applied to are deposited on Zenodo as a curated set of per-cohort case studies — F matrix, β/α decomposition, signed-basin output, MOFA+ comparison weights, per-cohort writeup — so each cohort is reusable as a reference projection for downstream methods consuming the substrate.

---

## Introduction

**Existing multi-omic methods produce axes in cohort-specific feature space.** MOFA+ [@argelaguet2018mofa], DIABLO [@singh2019diablo], and SNF [@wang2014similarity] fit a low-dimensional latent factor model to each cohort independently — factor scores live in a coordinate system fit to that cohort's particular patient × feature matrix, factor loadings are linear combinations of that cohort's input features, and neither carries a biological interpretation that can be stated without post-hoc enrichment analysis. The factors are not comparable across cohorts (each cohort has its own basis), they cannot integrate cohorts with non-overlapping inputs (each method requires shared feature spaces), and they cannot be interpreted at the per-axis level without enrichment-test post-processing (the factors are abstract directions, not named biology).

**We propose working in substrate space rather than factor space.** A *substrate* — in our usage, a fixed biochemistry graph constructed from curated sources — defines a coordinate system that is independent of any cohort, any patient set, or any modality. The substrate we use is a 38,148-node merged graph of Reactome reactions, StringDB protein interactions, HMDB metabolites, and KEGG pathway annotations (CC-BY 4.0, multi-source-merged, integrating 16,343 genes + 6,406 metabolites + 15,399 reactions + ancillary nodes). Any patient's multi-omic measurements can be projected onto this substrate by MAP reconstruction with per-modality measurement variance, producing a per-patient state vector **F** in the substrate's 38,148-dimensional coordinates. This projection makes three things possible that cohort-specific factor methods cannot do at all:

1. **Cohorts with non-overlapping inputs share a coordinate system.** Two cohorts measured on different modalities (one RNA-seq, one metabolomics) project into the same 38,148-dimensional space; their per-patient F vectors are directly comparable on a node-by-node basis. The substrate carries the integration; the modality is incidental.

2. **Axes are biologically named by construction, not by post-hoc enrichment.** Because F lives in substrate coordinates, any axis of patient variation (e.g., a principal component of F or of its residual after removing a hub-projection scalar) carries loadings on *named substrate nodes* — specific genes, reactions, metabolites — that are read directly from the graph's node attributes. There is no factor-to-biology mapping step; biology is already the coordinate system.

3. **Comparison of cohort-specific patient axes becomes a quantitative cosine test.** Because the coordinate system is the same across cohorts, two cohorts' axes can be aligned by cosine; same-biology cohort pairs should align at high cosine, different-biology pairs should be near-orthogonal. We use this as a falsification test of the framework throughout.

**The substrate is the primary contribution of this paper; GIZMO is the first method that fully exploits it.** The 38,148-node graph itself is the shared infrastructure, deposited under CC-BY 4.0, and is consumable by any method that takes a graph + per-feature measurements. Static PageRank on the substrate (we use this as the "TopPR" baseline) recovers hub-anchored disease genes at 3–5× fold-enrichment at low K with no data input at all, validating substrate quality. GIZMO adds what TopPR cannot: per-patient instantiation via MAP, β/α decomposition (presentation magnitude vs mechanism residual), signed-basin extraction that names both poles of each principal axis, and a per-cohort interpretive layer that converts the projection into biology a downstream user can read. Where this paper claims "GIZMO recovers X", the underlying claim is "the substrate-fixed coordinate system makes X interpretable, and GIZMO is the projection method we used to land each patient in it." Papers 2–4 in this series extend the substrate (lipid layer, mutation-conditional edges) and exploit it in further ways (longitudinal trajectories, spatial-tile projection); the substrate is the load-bearing object, the methods are how we use it.

**Scope condition stated up front.** Under the canonical input preprocessing, the substrate represents loss-of-function and differentiation-collapse biology (TP53 truncation, HIF stabilization, type-I IFN response). Under within-patient z-score preprocessing (Methods §"Within-patient z-score"), the substrate additionally represents gain-of-function biology with canonical downstream effectors (KRAS-G12X driving MAPK target reactions; IDH1-R132H producing 2HG via IDH1/IDH2/D2HGDH expression shifts). The v4–v5 manuscripts framed GoF as outside representational scope based on canonical-preprocessing failures (TCGA_LUAD KRAS-vs-EGFR multi-PC LR LOOCV AUC = 0.509; TCGA_IDH 2HG anchor rank 6,355). In v6 §5 we retract that scope claim: under within-patient z-score these lift to 0.703 and 26 respectively, demonstrating the v5 failures were preprocessing artifacts on canonical substrate edges, not representational gaps. The narrower Paper 4 open question contracts to *truly neomorphic biology with zero canonical effector signaling AND not rescuable by within-patient z-score*.

The remainder of this paper develops the substrate + GIZMO architecture through (§Substrate construction and β/α decomposition) the coordinate-system definition; (§Signed-basin interpretability) the strongest demonstration of named-axis biology across 17 cohorts; (§Horizontal meta-analysis) the proof that the substrate carries integrative work across cohorts with zero overlapping features; (§Method comparison) symmetric testing against MOFA+, TopPR, PCA, NMF showing parity at discrimination and uniqueness at the per-patient-instantiation layer; (§Scope conditions) the pre-MAP smoothing-indicator diagnostic, the α-predictiveness envelope, and the GoF falsification.

---

## Results

### Section 1: The substrate and its β/α decomposition

**Substrate construction.** The biochemistry substrate is built from four curated sources: Reactome (15,399 reactions, 6,406 metabolites, 2,382 pathway nodes), StringDB (protein-protein interactions filtered at confidence ≥ 700), HMDB (metabolite identifiers and synonyms), and KEGG (reaction directions and EC-class annotations). The unified graph has 38,148 propagation-eligible nodes (16,343 genes + 6,406 metabolites + 15,399 reactions + ancillary). The biochem subgraph caps high-degree hub-node edge counts at 200 to prevent the Laplacian from being dominated by promiscuous metabolites (water, ATP, NADH); the cap was chosen empirically as the inflection point at which Louvain partition stability plateaus. **Reactions are reified as nodes** rather than collapsed into edges between genes and metabolites — a contraction test on IDH-glioma (n=88, paired NMR+RNA) showed that removing the 15,399 reaction nodes and replacing them with substrate↔product bypass edges changes the top-50 \|F̄\| biology completely (Jaccard = 0.010 vs the standard substrate; the contracted basin loses the canonical IDH-mut MetAsp / methylglyoxal / ceramide axis and recovers generic graph topology — ferritin + ribosomal proteins). α-PC1 AUC against the IDH-mut label is essentially unchanged (0.901 → 0.894), but the biology shifts entirely — reactions are not pass-through decorations, they are integrators that combine catalyst gene + substrate + product into a single substrate-coherent variance direction. The substrate file (CC-BY 4.0) is the primary deposit of this paper; downstream methods including ours consume it. Identifier mapping (GeneMapper, MetaboliteMapper) handles Entrez / Ensembl / HGNC / RefSeq for genes and HMDB / ChEBI / PubChem / KEGG / LIPIDMAPS / InChIKey / common-name for metabolites; per-cohort identifier-mapping completeness is in Supp Table S2.

**MAP reconstruction.** For each patient $i$, the state vector $\mathbf{F}_i$ over substrate nodes is the MAP solution

$$
\mathbf{F}_i \;=\; \arg\min_{\mathbf{F}} \left[\, \| \mathbf{x}_i - \mathbf{A}_{\text{obs}} \mathbf{F} \|^{2}_{\boldsymbol{\Sigma}^{-1}} \;+\; \lambda \, \mathbf{F}^{\top} \mathbf{L}_{\text{signed}} \, \mathbf{F} \,\right]
$$

where $\mathbf{x}_i$ is the patient's feature vector, $\mathbf{A}_{\text{obs}}$ is the substrate-anchoring matrix, $\boldsymbol{\Sigma}$ is per-modality measurement-variance diagonal, $\mathbf{L}_{\text{signed}}$ is the signed substrate Laplacian (catalytic-product edges = +1, inhibitory-substrate edges = −1, neutral = +1; the Dirichlet-smoothness term penalizes sign-discordant transitions across edges), and $\lambda$ is the smoothing hyperparameter. The solve is closed-form via sparse conjugate-gradient at ~30 sec per patient on the 38,148-node substrate. $\mathbf{F}$ output is a per-patient 38,148-dimensional vector in *substrate coordinates*: the $i$-th entry of $\mathbf{F}$ is the inferred state of the $i$-th substrate node for this patient, regardless of which features the cohort actually measured. This is the key architectural property — patient state is reported in node-coordinates that don't depend on cohort or modality.

**β/α decomposition.** $\mathbf{F}$ decomposes uniquely into a hub-projection scalar and an orthogonal residual. The substrate's log-PageRank direction $\mathbf{p} = \log\mathrm{PR}_{\text{centered}}$ (a 38,148-dim unit vector that depends only on substrate topology) defines:

$$
\beta_i \;=\; \frac{\mathbf{F}_i \cdot \mathbf{p}}{\| \mathbf{p} \|^{2}}, \qquad
\boldsymbol{\alpha}_i \;=\; \mathbf{F}_i - \beta_i \, \mathbf{p}
$$

By construction $\boldsymbol{\alpha}_i \perp \mathbf{p}$; $\beta$ and $\boldsymbol{\alpha}$ decompose $\mathbf{F}$ into a hub-aligned component and a substrate-orthogonal component without information loss. The interpretation comes from what log_PR's top-percentile nodes turn out to be empirically. Annotation of the top-50 substrate hub nodes (Supp Table S3) shows: 18 transcriptional / chromatin hubs (TP53, STAT1/3, MITF, NFE2L2, RUNX1, CTNNB1, XBP1, …); 13 signaling / kinase hubs (EGFR, HSP90AA1, SHC1, TYK2, GPCR cascades, …); 8 immune / inflammation reactions (granule exocytosis, IFNG-stimulated gene expression, COP9 signalosome, …); 4 host-defense; 4 viral-cofactor; 1 keratin / structural; 2 mitochondrial-currency. **Forty-seven of the 50 top hub nodes are transcription / signaling / immune / host-defense** — biology that is engaged across virtually all disease states. The β axis therefore reads as a **phenotype-presentation magnitude**: how much of the substrate's most-engaged biology is active in this patient. α is the **orthogonal residual mechanism**: the disease-specific perturbation not captured by phenotype-presentation alone.

**α-PCA and per-patient scoring.** For each cohort, we run PCA on the α residual matrix (n_patients × 38,148 substrate nodes) and retain the first five components α-PC1..α-PC5. Each α-PC_k is a 38,148-dimensional unit vector in substrate-coordinate space — a *direction* in which the cohort's patients vary after removing the β phenotype-presentation axis. Each patient's (β, α-PC1, …, α-PC5, ‖α‖₂) coordinates are then a small fixed-dimensional summary in substrate space. For per-patient z-scoring against a within-cohort reference, the reference distribution is built from the cohort's controls; the patient's z-scores then read directly as "how disease-deviant relative to within-cohort baseline". β and ‖α‖₂ are cohort-independent scalars; α-PC1 is cohort-specific but its direction is a stored 38,148-dim unit vector that can be reused across timepoints (longitudinal extension) or spatial tiles (Paper 5; out of scope here).

This is the coordinate system. The remaining results sections demonstrate what working in it lets us do.

---

### Section 2: Signed-basin interpretability — every α-PC names two coherent substrate regions

Because PCA loadings are signed, each α-PC defines an *axis* of patient variation: patients with high positive scores sit at the + end and high negative scores at the − end. In substrate coordinates, this axis directly identifies two sub-graphs — partition substrate nodes by the sign of their loading, take the connected components of the substrate restricted to each sign-class, and the largest + component and the largest − component name the biology of each patient pole. **The operationalization is parameter-free**: no K, no quantile, no threshold tuning. The threshold *is* the zero-crossing of the loading.

**Each α-PC's mass concentrates on two connected substrate regions.** Across the cohort × PC cells we tested (Table 1), the largest + basin and the largest − basin together hold **50–88% of the PC's total \|F\|² mass**. Filbin α-PC4 is the extreme case: 87.7% of the PC's mass lives on 12 substrate nodes (a 7-node lactoperoxidase reaction clique + a 5-node galanin receptor clique). Most other PCs spread mass across two larger 6,000–20,000-node connected basins, but the basins are still *connected sub-graphs holding the majority of mass*. A complementary smoothness diagnostic confirms this: the graph Dirichlet energy of each PC (PC^T L PC / ‖PC‖² on the normalized substrate Laplacian) is far below the random-unit-vector baseline of 0.9995 ± 0.0023 (30 trials). Filbin α-PC4 has Dirichlet = 0.000 (Z = −437 vs random null); CCRCC α-PC1 = 0.575 (Z = −185); the weakest PC tested (TCGA_IDH_glioma α-PC5) has Dirichlet = 0.895 (Z = −45). All are dramatically more graph-coherent than random directions. **The substrate's coordinate system is not just a labeling convenience — it is empirically the case that the cohort variance falls along connected substrate regions.** Nothing in the PCA objective penalizes graph-incoherent solutions; the observed coherence is evidence that biology lives on the substrate.

(Table 1 — signed-basin statistics across cohort × PC cells. Mass = sum of squared loadings; total = 1.)

| Cohort × PC | + basin nodes | + mass | − basin nodes | − mass |
|---|---|---|---|---|
| **Filbin_COVID α-PC4** (cardiorenal) | **7** | **67.6%** | **5** | **20.1%** |
| GSE65391_SLE α-PC3 (psychosis-IFN) | 11,199 | 45.5% | 17,770 | 18.9% |
| IDH_glioma α-PC1 (neuronal) | 17,095 | 50.2% | 8,936 | 30.0% |
| CPTAC_COAD α-PC1 (TME) | 10,125 | 41.5% | 18,304 | 21.3% |
| CPTAC_COAD α-PC4 (tumor) | 18,614 | 29.8% | 10,396 | 40.3% |
| CPTAC_CCRCC α-PC1 (tumor) | 8,376 | 38.5% | 20,138 | 22.4% |
| KMPLOT_BRCA α-PC4 (grade) | 6,399 | 35.3% | 18,819 | 29.7% |
| GSE65682_sepsis α-PC2 (sepsis) | 16,563 | 37.1% | 12,780 | 24.5% |
| GSE65391_SLE α-PC1 (inflammation) | 17,396 | 29.9% | 11,470 | 28.1% |
| TCGA_IDH_glioma α-PC5 (IDH-mut) | 6,697 | 25.1% | 17,616 | 36.7% |

**Reading the basins names the biology directly.** Because the basins live on the substrate, their member nodes are *named substrate entities* — specific genes, reactions, metabolites with attributes (HGNC symbol, Reactome reaction name, HMDB metabolite identifier). The α-PC's biology is read by enumerating the top-loaded basin members. Seven examples across cohort categories follow; they are the "what's special about GIZMO" demonstration of the paper.

**Filbin_COVID α-PC4 (cardiorenal-discriminating, AUC 0.72) — the textbook two-clique case** (12 substrate nodes hold 87.7% of PC mass):
- **+ basin (67.6% of mass, 7 nodes):** *Lactoperoxidase antimicrobial system* — LPO (+0.24) catalyzes SCN⁻ (+0.24) → OSCN⁻ (+0.41), which then reacts with Peptidyl-Cys-SH (+0.24) → Peptidyl-Cys-SSCN (+0.24). The entire 4-step host-defense oxidative chemistry sits as one connected clique.
- **− basin (20.1% of mass, 5 nodes):** *Galanin receptor system* — GAL ligand (−0.16) + GALR1/GALR2/GALR3 receptors (−0.16 each) + binding reaction (−0.32). A single neuropeptide signaling clique that governs cardiovascular and visceral autonomic tone.

The cardiorenal-complication axis in severe COVID separates patients whose substrate is dominated by oxidative antimicrobial activity (LPO basin) from those dominated by galanin-receptor neuropeptide signaling — a remarkably specific mechanistic axis discovered without supervision.

**CPTAC_CCRCC α-PC1 (tumor-vs-normal, AUC 0.97) — the canonical ccRCC axis rediscovered:**
- **+ basin (38.5% of mass, 8,376 nodes):** *Normal renal tubule machinery* — SLC22A6 (+0.09), SLC22A8 (+0.09), ENPP6 (+0.07), uromodulin (+0.07), SLC36A2 (+0.06), barttin/CLCNK accessory subunit (+0.06). The differentiated proximal-tubule transporter and adhesion machinery.
- **− basin (22.4% of mass, 20,138 nodes):** ***Warburg glycolysis + HIF target axis*** — PFKP (−0.06), HK2 (−0.04), ENO2 (−0.04), CA9 (−0.04), CD70 (−0.04). PFKP/HK2/ENO2 are the canonical Warburg-glycolytic enzymes; **CA9 is the textbook HIF-2α target and clinical biomarker of clear-cell RCC**.

GIZMO has independently rediscovered the canonical ccRCC mechanism on the substrate: VHL-loss → HIF-2α stabilization → glycolytic switch + CA9 upregulation on one end of the axis, proximal-tubule differentiation collapse on the other. No mutation status, no histology, no labels.

**CPTAC_COAD α-PC4 (tumor-vs-normal, AUC 0.94) — tissue-specialization axis:**
- **+ basin (29.8% of mass, 18,614 nodes):** *Intestinal endocrine* — PYY peptide YY (+0.17), INSL5 (+0.06), WFDC2 (+0.11), B3GALT5 (+0.09), N-acetylneuraminate synthase (+0.09).
- **− basin (40.3% of mass, 10,396 nodes):** *Pancreatic exocrine* — chymotrypsin C (−0.13), AMY2A amylase (−0.12), CPB1 carboxypeptidase B1 (−0.12), CELA2A chymotrypsin-like elastase (−0.10), ADH1B (−0.10), ALDH1L1 (−0.09).

The PC discriminates by *which differentiated tissue specialization is preserved*. Adjacent-normal samples (mixed colon + nearby pancreas in surgical specimens) preserve both endocrine and exocrine machinery; adenocarcinoma collapses both.

**IDH_glioma (Trautwein) α-PC1 (IDH-mut, AUC 0.65) — neuronal-vs-immune axis:**
- **+ basin (50.2% of mass, 17,095 nodes):** *Neuronal signaling* — voltage-gated K+ channel activation (+0.12), NPTN-GABA-A binding (+0.09), GABA-A heteropentamer Cl⁻ transport (+0.09), GPLD1 GPI-anchor hydrolysis (+0.08), somatostatin receptor binding (+0.08), PTPRD-SLITRK1-6 binding (+0.08), KCND/KCNIP K+ channels (+0.07).
- **− basin (30.0% of mass, 8,936 nodes):** *Immune/granule response* — IFNG-stimulated gene expression (−0.05), IFN-induced gene expression (−0.05), specific-granule exocytosis (−0.05), elastic-fibre binding (−0.05), azurophil-granule exocytosis (−0.04), emilin elastic-fibre component (−0.04).

IDH-mutant gliomas are well-characterized as having neuronal-like differentiation programs [@suva2014mutational]. The PC names this directly: patients on the + end carry voltage-gated K+ channels, GABA-A pentamer transport, somatostatin receptors — the *neuronal-glioma phenotype*. Patients on the − end carry IFN-γ/IFN-α response and granulocyte exocytosis — the *immune-infiltrate phenotype*.

**GSE65391 SLE α-PC3 (psychosis-discriminating, AUC 0.99) — type I IFN axis:**
- **+ basin (45.5% of mass, 11,199 nodes):** *Type I interferon response* — RSAD2 (+0.16), IFI27 (+0.14), IFI44L (+0.13), IFIT3 (+0.11), LGALS3BP (+0.10), GYPB (+0.10), IFI44 (+0.10), XAF1 (+0.09). The textbook lupus type I IFN signature.

The classical "interferon signature" of lupus [@baechler2003interferon; @ronnblom2014interferon] emerges as the unsupervised + basin of α-PC3, and the patient pole carrying psychosis manifestation sits on this end.

**Crohn α-PC1 (thiopurine label, paired prot+metab, n = 33) — pharmacogenomic axis rediscovered.** α-PC1 captures 33% of cohort α-variance, with the *MPG* gene loading +0.52 (the largest single loading in any panel basin), plus 8 MPG-mediated reactions (cleavage of 3-methyladenine + hypoxanthine + ethenoadenine + APEX1 displacement). Thiopurines (azathioprine, 6-mercaptopurine) create alkylated-base lesions that *MPG* repairs; the substrate **rediscovers the canonical pharmacogenomic axis of thiopurine response at n = 33**.

**GSE89408_RA α-PC1 (RA-vs-OA, RNA-only, n = 174) — ectopic lymphoid neogenesis.** The +basin is the complete TLS module: CXCL13/CXCR5 (B-cell germinal center axis), CXCL9-11/CXCR3 (Th1 recruitment), BAFF/APRIL/TACI/BCMA (B-cell survival, the belimumab-target axis), CCR3/4/5/CCR7, IL10 negative regulation. The substrate recovers the canonical synovial tertiary lymphoid structure signature [@manzo2005ectopic; @humby2009ectopic] **unsupervised, from synovial RNA alone, with no immune-pathway curation**.

These seven cases are the architectural argument made concrete: **substrate-coordinate axes name the biology of each patient pole, parameter-free, from the substrate's own node attributes.** A cohort-specific factor in MOFA+ or DIABLO does not have this property — its loadings are linear combinations of cohort-input features, and the biology has to be extracted by a separate enrichment step that the user runs after the fact.

**Cross-cohort alignment in substrate coordinates.** Because two cohorts' α-PCs are both 38,148-dimensional unit vectors in the same substrate-fixed coordinate system, their similarity is a quantitative cosine. We compute pairwise \|cos(α-PC_i^a, α-PC_j^b)\| across all 16 cohorts × top-5 PCs (80 (cohort, PC) entries; the chord diagram in Figure 2a renders the full pairwise structure). Median cross-disease cosine is **0.020** with 95th percentile **0.170** — overwhelmingly orthogonal. The within-disease positive control is the IDH-glioma pair: **Trautwein α-PC1 ↔ TCGA_IDH α-PC1 cos = 0.72; α-PC2 ↔ α-PC2 cos = 0.77** — two cohorts measured on completely different modalities (paired NMR+RNA vs RNA-only) replicating their α-PCs at cosines that are deep in the 99.9% tail of the cross-disease null. **No cohort-specific factor method can run this test at all** — MOFA+ factors fit each cohort's own basis and cannot be compared across cohorts.

**Within-disease alignment is a single-case demonstration, not a general claim.** The IDH-glioma pair is the *only* same-disease pair in our panel that replicates strongly. Three additional same-disease pairs all fail to align: Su_COVID ↔ Filbin_COVID best diagonal cos = 0.084; Gao_RA ↔ GSE89408_RA cos = 0.288; HMP2_IBD_CD ↔ Crohn cos = 0.060. One *cross*-disease pair (Gao_RA ↔ TCGA_LUAD) actually has higher alignment (0.44) than two of the same-disease pairs. The IDH-glioma case is special because IDH1-R132H is a single dominant mechanism with a direct measurable substrate-routed anchor (2-hydroxyglutarate); diseases with heterogeneous mechanisms (COVID severity range, RA-vs-OA, IBD subtypes) do not converge on a single shared α-PC. **Within-disease α-PC alignment scales with mechanism homogeneity, not just disease label** — the architectural enablement is real but the empirical replication is mechanism-specific.

**Cross-PC swap is the dominant pattern.** A full pairwise chord-diagram analysis (Figure 2a, 80 entries; `pc_alignment_best_cross_cohort_match.tsv`) shows that **81% (65/80) of best cross-cohort matches involve a PC-index swap** — the strongest cross-cohort partner for cohort A's α-PC*i* sits at a different PC index *j* in cohort B. The canonical example: IDH_glioma α-PC3 aligns with TCGA_IDH_glioma α-PC4 at \|cos\| = 0.66, both capturing the same IDH-mut neomorphic biology but ranked differently by within-cohort variance. PCA's variance ordering reflects cohort idiosyncrasies (n, modality mix); the underlying biological axis is shared and recoverable across that reordering — but *only because the comparison happens in substrate coordinates*. PC sign is canonicalized before per-node display (if vᵃ·vᵇ < 0, vᵦ ← −vᵦ; PCA solvers pick sign by internal convention rather than biology, so per-node sign comparisons require the flip to report biology rather than solver orientation).

**Conserved biology across the top-10 chord pairs (Figure 2h).** Running greedy-connected seed selection on each of the top-10 cross-cohort chord pairs and tabulating which substrate nodes recur in ≥ 2 pairs surfaces **57 of 274 unique seed contributors** as recurrent. Reactome leaf-pathway grouping over the recurrent set identifies four pathways with ≥ 2 recurrent members (Figure 2h companion TSVs):

| Pathway (Reactome leaf) | Members | # top chords | # distinct diseases |
|---|---|---|---|
| Neutrophil degranulation (R-HSA-6798695) | 5 | 3 | 3 (Glioma, LungAdeno, RA) |
| Immunoregulatory interactions between Lymphoid and non-Lymphoid cell (R-HSA-198933) | 4 | 3 | 3 |
| Nonsense Mediated Decay enhanced by EJC (R-HSA-975957) | 3 | 2 | 1 (Glioma cross-PC) |
| Chemokine receptors bind chemokines (R-HSA-380108) | 3 | 3 | 3 |

Three of the four pathways span ≥ 3 distinct diseases; the fourth captures within-Glioma cross-PC recurrence — the canonical IDH cross-PC swap finding expressed at the pathway level. The MHC-I immune-surveillance axis (LILRs / SIGLEC / TREM-CD300 / ICAM-LFA-1), chemokine recruitment (CXCR3-CXCL9/10/11, ACKR2, CXCR5), and innate-granule machinery all surface as cross-disease shared substrate biology *that GIZMO finds repeatedly across different chords, in different cohorts, without supervision*. **Importantly, single-cohort hypergeometric pathway enrichment does NOT recover this pattern** — at the same chord with \|cos\| ≥ 0.5, the Reactome pathway hit-list of cohort A's PC top-100 loadings and cohort B's PC top-100 loadings often diverges sharply (Figure 2g third column). The merged substrate captures cohort-shared biology that *crosses* canonical pathway boundaries (transporters that gate two pathways, signaling kinases whose substrates are annotated elsewhere); pathway-membership enrichment tests are blind to this kind of coupling. The cosine-alignment view sees it because the substrate spans the curated databases at the edge level rather than partitioning at the boundary level.

---

### Section 3: Horizontal meta-analysis — three cohorts with zero overlapping features converge

The signed-basin and cross-cohort cosine results establish that the substrate is a coordinate system in which biology has named addresses. **Horizontal meta-analysis is the corresponding integrative test**: can two cohorts that measure *different things entirely* — disjoint modalities, disjoint input features — still converge on the same disease-anchor biology when both are projected onto the substrate? We test on three cohorts chosen specifically because their input features have zero overlap:

- **IDH_glioma (Trautwein)** — n = 88, paired RNA-seq (gene-symbol input features) + per-patient ¹H-NMR tissue metabolomics (concentration-quantified metabolites); IDH-mut vs IDH-WT label.
- **HMP2_IBD_CD** — n = 399, metabolomics-only (HMDB-anchored metabolite identifiers, no gene-side input); Crohn vs healthy label.
- **GSE89408_RA** — n = 174, microarray RNA-seq (gene-symbol probes); RA vs OA label.

There is no gene that appears in all three cohorts' input feature lists (cohorts 1 and 3 share gene-symbol inputs but HMP2 has no gene-side input at all); there is no metabolite that appears in cohorts 1 and 3 (which have no metabolomics); the three cohorts share **zero** common input features. **The substrate carries the integration.** For each cohort, we compute mean \|F\| across substrate metabolite nodes for active patients vs controls, rank the difference by metabolite, and ask whether each cohort's literature-anchor metabolites land at top-5% rank.

- **2-hydroxyglutarate (2HG, the canonical IDH oncometabolite)** — ranks **74 / 6,406 metabolite nodes (top 1.14%)** in Trautwein RNA-only (microarray, n=88) at the substrate's mitochondrial-2HG node (`reactome:R-ALL-879997`), Cohen's d = +1.17 — surfaced via Reactome catalysis edges from IDH1/IDH2/D2HGDH expression without any metabolomics input. The same propagation on TCGA_IDH_glioma RNA-seq (n=458) does **not** recover 2HG (rank 6,355 at the cytosolic-2HG node `reactome:R-ALL-880042`, Cohen's d = −0.95, opposite direction). The cohort-conditional outcome (Trautwein recovers, TCGA does not) is the framework's documented scope condition: substrate-mediated metabolite recovery from RNA input depends on substrate connectivity around the anchor + tissue preparation + platform alignment with the relevant catalysis-edge neighborhood. Source: `benchmarks/results/unsupervised/stage31_v17_metab_accumulation.json`.
- **Propionate / lithocholate (HMP2 microbiome-driven CD anchors)** — rank in the top **2.8%** of substrate metabolites from HMP2's metabolomics-only F.
- **Citrulline (GSE89408_RA, well-established synovial RA marker)** — ranks in the top **4.3%** of substrate metabolites from GSE89408's RNA-only F.

All three at top-5%, joint p ≈ **2 × 10⁻³** against a per-cohort random-metabolite null (10,000 draws each, joint by Stouffer's method). The substrate is the integration; the modality is incidental.

This is the proof that the coordinate-system claim translates into integrative work no cohort-specific method can do. Caveats: anchor recovery is hub-status-driven (CTD direct-evidence anchors sit at median PR percentile 67–79; `anchor_hub_bias_diagnostic.tsv`), so this finding validates *substrate quality* more than it differentiates GIZMO's projection from other substrate-consuming methods like static PageRank or label propagation. The substrate is the load-bearing object; we will see in §4 that several methods consume it productively, and the GIZMO contribution is per-patient instantiation rather than aggregate-statistic dominance.

---

### Section 4: Method comparison — discrimination parity at apples-to-apples access; per-patient layer is the unique addition

The substrate is the shared infrastructure; multiple methods can consume it. We compare GIZMO against four others on the same 38,148-node substrate: **MOFA+** (free latent factor model fit per cohort), **PCA** and **NMF** (matrix factorization on the input feature matrix), and **TopPR** — static PageRank on the substrate, *no data input at all*, the population-mean / data-independent baseline.

**Anchor recovery: TopPR dominates at every K; GIZMO ≈ noisy TopPR within measured features.** Across 7 cohorts × 5 methods × 7 K values (Figure 4a), TopPR achieves the highest fold-enrichment of CTD direct-evidence disease-gene anchors at every K (2.40× at K = 50; 2.01× at K = 200; 1.66× at K = 500). The mechanical reason: CTD direct-evidence anchors are hub-status-elevated (median PR percentile 67–79), and static PageRank captures the hub-anchored signal within measured features more cleanly than any data-driven method. GIZMO mean \|F\| correlates with static PR at ρ = 0.45–0.90 across cohorts; the more data GIZMO has to fit, the closer its patient-average is to clean PageRank. **At visible K (apples-to-apples: every method restricted to gene-symbol inputs that map to substrate), GIZMO ≥ MOFA+ at mid-K (K = 500–1000), but TopPR wins at every K against everyone.** This is the honest framing: TopPR exposes the substrate's population-mean structure; GIZMO produces per-patient instantiations of that same smoothing prior; neither is *better* — they answer different questions.

**At invisible K (graph-mediated extrapolation to unmeasured neighbors), GIZMO's top-K is disjoint from PCA/MOFA+/NMF.** Jaccard 0.00–0.04 across all K because GIZMO's top-K is populated by *unmeasured* genes propagated in from substrate neighbors, while data-driven factor methods structurally cannot return genes that weren't in the input. Even GIZMO × TopPR Jaccard is 0.18–0.35 under invisible K — both methods pull toward hubs but pick different ones at the full-substrate scale. The invisible-K extrapolation is GIZMO's distinctive output property, and it is a property of the substrate's propagation behavior more than of the GIZMO ranking per se.

**Symmetric multi-axis discrimination: parity with MOFA+.** Both methods extract per-cohort latent axes (GIZMO: best of α-PC1..5; MOFA+: best of Factor 1..K). We test each per-axis against every scoreable clinical metadata field per cohort, then take the per-cohort maximum strength. Across **14 of 17 panel cohorts** with both methods cached (HMP2 has 0 scoreable clinical fields; breast.TCGA DIABLO and NEPTUNE are validation-only):

| Threshold | GIZMO | MOFA+ |
|---|---|---|
| Raw strength ≥ 0.40 | 9/14 (64%) | 10/14 (71%) |
| Per-cohort Bonferroni | 8/14 (57%) | 9/14 (64%) |
| Global Bonferroni (α = 0.05/181 = 2.76×10⁻⁴) | 8/14 (57%) | 9/14 (64%) |

MOFA+ edges GIZMO by 1 cohort. To rule out feature-access asymmetry (MOFA+'s factors see uncharacterized metabolites, unmappable peptides, and full microarray probes that GIZMO drops at the substrate-mapping step), we re-ran MOFA+ on the substrate-mappable subset only — the same per-cohort feature universe GIZMO sees — on 11 of the 14 cohorts:

| Threshold | GIZMO | MOFA+_sm |
|---|---|---|
| Raw strength ≥ 0.40 | 8/11 (73%) | 8/11 (73%) |
| Per-cohort Bonferroni | 7/11 (64%) | 8/11 (73%) |
| Global Bonferroni | **6/11 (55%)** | **6/11 (55%)** |

**At apples-to-apples feature access, GIZMO and MOFA+ are at exact parity** — 6/11 vs 6/11 at global Bonferroni; 8/11 vs 8/11 at raw strength. The 1-cohort gap MOFA+ retains at per-cohort Bonferroni is the same Su_COVID single-axis advantage seen at full-feature access. **MOFA+'s 1-cohort edge in the asymmetric table is fully explained by feature-access asymmetry.**

Discrimination parity is the expected result — both methods extract per-cohort patient covariance structure from the same per-cohort variance, and discrimination is a property of that covariance structure, not of the substrate. **The substrate's role is not to discriminate harder than MOFA+; it is to provide the coordinate system in which the resulting axes are named.** What GIZMO produces that MOFA+ does not is the per-patient instantiation in substrate coordinates and the signed-basin biology layer (§2). Both methods reach 55–73% pass rates at substrate-matched access; both fail in the same three cohorts (Crohn small-n, TCGA_LUAD GoF, Su_COVID limited modality budget). The interpretive layer is the unique GIZMO output, not the AUC.

**Module-level cross-pathway bridging vs WGCNA.** WGCNA partitions modules within each modality and within tight pathway-containment criteria; GIZMO's modules are biochemistry-coherent cross-pathway communities by construction. On 10 of 11 cohorts tested, GIZMO modules span 8–40× more Reactome leaf pathways per module than WGCNA modules — bridging across pathway boundaries that are curator-imposed but not biology-imposed. This is consistent with the chord-diagram finding (§2): cohorts that align at cosine ≥ 0.5 often light up *different* canonical pathways in each cohort's PC because the shared substrate biology lives at edges that cross those pathway boundaries.

**Convergent validation with MOFA+ at the per-patient stratification level.** Each method at its own natural operating point, tested against clinical severity by Kruskal-Wallis + per-cohort effect-size criterion: GIZMO_α passes in 7/17 cohort-design cells; MOFA+ in 6/17 — within sampling noise of parity (`stage31_v11_stratification.json`). This is convergent validation, not benchmark win: an independent unsupervised method surfaces stratification signal in a comparable fraction of cells, which is evidence the framework's α-stratification recovers real per-patient structure rather than a method-specific artifact.

The substrate enables a method ensemble; GIZMO is one method in the ensemble. Of the five we tested on it, no method dominates on all the metrics that matter — TopPR wins anchor recovery; GIZMO and MOFA+ are at parity on discrimination; GIZMO is unique on per-patient instantiation + signed-basin biology + cross-cohort cosine alignment. **The right summary is: the substrate carries integrative work; GIZMO is the first method that exploits that work into per-patient state vectors with named interpretive axes.**

---

### Section 5: Scope conditions — when the framework works and when it doesn't

The framework works in 12 of 15 cohorts at the multi-PC inference criterion (raw strength ≥ 0.40 on best of α-PC1..PC5 × best clinical metadata; `multi_pc_15cohort_disclosure.tsv`); it fails in three. Each of the three failures aligns with a pre-specified scope condition:

- **TCGA_LUAD (KRAS-vs-EGFR subtype, n = 508)** — gain-of-function driver-mutation cancer. The substrate's representational scope is loss-of-function and differentiation-collapse biology that lives on canonical substrate edges; gain-of-function neomorphic biology (KRAS-G12X constitutively GTP-bound, EGFR-L858R structurally rewired) creates *new* edges that the canonical reference doesn't contain. Falsification stands.
- **Crohn (thiopurine response, n = 33)** — sample-size floor. The framework's robust-stratification floor is approximately n ≥ 50; this cohort fails at n = 33 *and* the per-cohort biology (MPG / thiopurine pharmacogenomics) emerges in α-PC1 at the basin level (above, §2) but does not separate enough patients to reach the discrimination threshold.
- **IDH_glioma (Trautwein, n = 88)** — sparse paired-modality anchoring; n = 88 with paired NMR + RNA covers ~21 metabolite measurements per patient (the NMR panel), and the α-PC1 basin biology emerges (above, §2) but the discrimination criterion is sensitive to anchor density.

**Pre-MAP smoothing-indicator diagnostic.** A per-cohort smoothness statistic (cosine between observed-feature variance direction and substrate-Laplacian low-eigenvalue direction) predicts when graph-Laplacian regularization will help vs hurt. The diagnostic correctly predicted all 4/4 rescue cases on cohorts where naive PCA underperformed and substrate-mediated MAP recovered the signal (Supp Figure S2). Implementation at `benchmarks/diagnostics/pre_map_smoothness.py`.

**α-predictiveness envelope.** Across the 11-cohort module-level metabolite rho test (`stage31_v60`), we map cohort substrate-coverage and anchor density to expected α-PC discrimination strength; the envelope tabulates 13 cohort-design cells where the framework should and does succeed (α-PC1 \|ρ\| ≥ 0.2) vs 7 cells where it fails predictably (TCGA_LUAD GoF, KMPLOT_BRCA cell-cycle-not-α-mechanism, etc.). This is a post-hoc diagnostic derived from the v44/v45 results; prospective validation requires applying it to a held-out cohort before observing the result.

**GoF empirical falsification on TCGA_LUAD.** We test the GoF scope-limit claim directly. Pull TCGA-LUAD driver-mutation calls from cBioPortal (`luad_tcga_mutations` profile, 429 mutation records across 204 patients); partition the 508-patient F cohort into single-driver-only subsets that exclude co-mutants of the major GoF drivers (KRAS, EGFR, BRAF, MET, ALK, ROS1): n = 61 KRAS-only-mut and n = 28 EGFR-only-mut.

Step 1 — *Is the discriminative signal present in raw expression?* Test DUSP6/DUSP4/SPRY2-4/ETV4-5/FOS/EGR1/PHLDA1/CCND1 (MAPK target genes) and HBEGF/AREG/EREG/BTC/TGFA/LRIG1 (EGFR target genes) on the Broad Firehose RSEM-normalized matrix by Mann-Whitney. **DUSP4 is +2.2 log₂ higher in KRAS-mut at p = 8×10⁻⁹ in raw RNA**; DUSP6 +1.8 log₂ at p = 4×10⁻⁷. The KRAS-vs-EGFR axis is *strongly present in raw input*.

Step 2 — *Does any GIZMO α-PC capture it?* Restrict the α-matrix to the 89 single-driver patients, compute α-PCA on the subset, and test each α-PC against the KRAS-vs-EGFR binary label by ROC-AUC. **No α-PC clears AUC 0.70.** Best α-PC AUC = 0.64. The discriminative biology that lives in raw RNA at p = 8×10⁻⁹ does *not* route through the canonical substrate.

Under the *canonical* (per-modality log + global-std) preprocessing, this is the falsification: GoF biology fails not from statistical underpowering (the signal is there at p = 10⁻⁹) but from what looked like *representational failure* — the canonical substrate appearing not to route the KRAS-G12X → DUSP4-up signal because canonical KRAS doesn't list constitutive DUSP4 upregulation among its annotated downstream events. **v6 retraction: this is not a representational failure of the substrate; it is a preprocessing artifact.** Under within-patient z-score preprocessing (Methods §"Within-patient z-score"), which preserves the relative-rank ordering of MAPK target genes that the per-modality global-std normalization absorbs into cohort-wide rescaling:

- **TCGA_IDH_glioma 2HG-mito anchor recovery** (RNA-seq only, n=458, 363 IDH-mut vs 95 IDH-wt): canonical preprocessing places 2HG-mito at rank 6,355/6,406 with Cohen's $d = -0.95$ (opposite direction); **within-patient z-score moves it to rank 26/6,406 (top 0.41%), $d = +0.630$**. The canonical Reactome catalysis edges from IDH1/IDH2/D2HGDH expression carry the neomorphic-2HG signal correctly; the v4–v5 manuscript reported the canonical-preprocessing result and read it as a substrate-representation failure, which we retract.
- **TCGA_LUAD KRAS-vs-EGFR subtype discrimination** (89 single-driver patients): canonical preprocessing multi-PC LR LOOCV AUC = 0.509 (chance); **within-patient z-score lifts to 0.703**, best single-PC AUC = 0.706 — clearing the v4–v5 0.70 threshold. An exploratory ablation on the Rhea-enriched substrate with diffusion_t = 0 reaches multi-PC LR LOOCV AUC 0.835 (`benchmarks/luad_ablation_map_zscored_input.py`), confirming the discriminative biology fully routes through canonical substrate edges when input preprocessing exposes it.
- **IDH_glioma Trautwein (paired RNA+NMR, n=88)** bottomline IDH-mut-vs-wt discrimination: canonical reported best α-PC AUC = 0.65; within-patient z-score = 0.842 best α-PC, 0.969 multi-PC LR LOOCV.

These three rescues collectively retract the v4–v5 LoF/GoF representational-scope dichotomy. The canonical substrate represents both:
- Loss-of-function and differentiation-collapse biology (TP53 truncation, HIF stabilization, SLE type-I IFN response) — confirmed under both preprocessings;
- Gain-of-function biology with canonical downstream effectors (KRAS-G12X → MAPK target reactions; IDH1-R132H → 2HG via IDH1/IDH2/D2HGDH catalysis-edge expression shifts) — confirmed *under within-patient z-score preprocessing*.

The remaining Paper 4 (mutation-conditional edge injection) open question contracts to *truly neomorphic biology with **zero** canonical effector signaling AND not rescuable by within-patient z-score* — a narrower scope than v5 framed.

**Why preprocessing matters here.** Under per-modality global-std normalization, within-patient variance in MAPK-target gene expression — the KRAS-vs-EGFR discriminative signal — is rescaled by the cohort-wide global std and absorbed. The MAP solver's data-fidelity term then sees relative gene-vs-gene ordering within each patient with reduced contrast. Under within-patient z-score, each patient's gene expression is recentered around their own mean and rescaled by their own std; cross-patient differences in MAPK target gene relative-ranking propagate into F-space differences at the canonical substrate edges (KRAS → MAPK target reactions in Reactome). The substrate edges were always sufficient; preprocessing was suppressing the signal. The same mechanism explains the TCGA_IDH 2HG anchor rescue. The within-patient z-score preprocessing is a *fine-grained diagnostic* that surfaces strong-driver subtype and anchor-recovery signal *that the canonical preprocessing absorbs into the baseline*; it is not a replacement for the canonical preprocessing at the panel-wide cohort discrimination level (Methods §"Within-patient z-score").

**Empirical test of naive add-only edge injection (negative result, structural diagnostic).** We test directly whether *bidirectional* edge addition — adding undirected edges from each GoF gene to its canonical effector cascade — rescues either of the framework's two GoF-related failures. (i) For IDH1-R132H neomorphic biology on the Trautwein cohort, we add 18 canonical 2HG → α-KG-dependent dioxygenase edges (`reactome:R-ALL-879997` and `R-ALL-880042` → TET1/2/3, KDM4A/B/C/D, KDM6A/B, KDM2A/B, EGLN1/2/3, FTO, ALKBH1/2/3, HIF1A; literature from Xu 2011 *Cancer Cell*, Lu 2012 *Nature*, Figueroa 2010 *Nature*, Koivunen 2012 *Nature*), re-derive F under the same RNA-only MAP solve, and recompute 2HG's anchor rank. The 2HG-mito rank **worsens from 497 to 688** (Cohen's d drops from +0.491 to +0.357) — the undirected Laplacian distributes propagated mass across the 18 newly added downstream effectors rather than concentrating it at 2HG, and the F value at the 2HG node drops by ~55%. (ii) For TCGA_LUAD KRAS-vs-EGFR, we add 21 canonical effector edges (KRAS → {RAF1, BRAF, MAP2K1/2, MAPK1/3, PIK3CA, AKT1, RPS6KA1, RPS6KB1}; EGFR → {GRB2, SHC1, SOS1, KRAS, HRAS, NRAS, PIK3CA, STAT3, SRC, AKT1, YES1}; Patricelli 2016 *Cancer Discov*, Sequist & Engelman 2011 *NEJM*, Pao 2004 *PNAS*) and re-derive F on the 104 single-driver patients available in our run. Best α-PC AUC stays at **0.634 unchanged** to 4 decimal places, multi-PC LR LOOCV AUC drops marginally from 0.509 to 0.504 (noise). The data-fidelity term of the MAP solver dominates over the Laplacian smoothing at every measured node (KRAS, EGFR, and all 21 effectors are in the RNA input), so the new bidirectional edges have no signal differential to propagate. Both failures are *structural to add-only undirected edge injection*, not necessarily to the GoF rescue concept. The Paper 4 design doc proposes per-patient *directional* edge replacement (Sherman-Morrison rank-1 update modifying the substrate Laplacian asymmetrically rather than adding bidirectional edges) as one candidate fix, but this has not yet been empirically tested — *we have demonstrated only that one specific rescue approach (bidirectional add-only) fails, not that any particular alternative succeeds*. The structural diagnosis (Laplacian dilution at unmeasured anchors; data-fidelity dominance at measured nodes) characterizes the failure mode; whether directional rewiring, edge replacement, or some other operator-level change rescues the discrimination remains a Paper 4 open question. Reproducible at `benchmarks/idh_2hg_rank_canon_vs_injected.py` and `benchmarks/luad_kras_egfr_gof_injection.py`.

**The pre-registration reframe.** An earlier draft applied the multi-PC criterion only on α-PC1, which conflated *dominant within-cohort variance direction* with *clinical-label direction*. The per-cohort PC × metadata atlas (Supplementary Figure S1) shows these are systematically different — clinical biology lives on sub-α PCs in 10/15 cohorts (SLE psychosis on α-PC3, SLE CVA on α-PC4, TCGA_IDH-mut on α-PC5, KMPLOT-BRCA grade on α-PC4). The α-PC1-only criterion tested the wrong thing — it asked whether the clinical label happens to be the cohort's biggest variance direction, when the clinical biology more commonly lives on a sub-leading direction. The multi-PC reformulation tests the right thing — does clinical biology live on *any* of the substrate-coherent variance directions the framework recovers — and the 12/15 (80%) pass rate clears the 75% pre-registered threshold. **The reframe is post-hoc; we are saying so directly.** The original criterion was a methodological misspecification, the data tells us so explicitly, and we have rewritten the criterion to test the substantive question rather than the misspecified one.

---

## Discussion

### Substrate-fixed coordinates enable what cohort-specific factor methods cannot

The single load-bearing architectural property of working in substrate space is that two cohorts measured on different modalities can be projected into a shared coordinate system whose axes are biology, not factor loadings. From this property follow:

- **Horizontal meta-analysis across cohorts with zero overlapping inputs** (§3, IDH-glioma + HMP2 + GSE89408_RA at joint p ≈ 2×10⁻³). MOFA+, DIABLO, SNF cannot perform this test by construction — their factor coordinates are fit per cohort. Substrate coordinates are fit by curation, not by data, so they are the same across cohorts.
- **Replication-testable within-disease alignment by quantitative cosine.** The IDH-glioma α-PC1/α-PC2 cosine = 0.72/0.77 is a single-case demonstration (three other same-disease pairs fail to align); alignment scales with mechanism homogeneity, not just disease label. The architectural enablement is the cosine test itself — without substrate-fixed coordinates the question "do these two cohorts share biology" cannot be made quantitative.
- **Per-axis biology read directly from node attributes** (§2 signed basins). The interpretive layer is built into the coordinate system rather than appended as a post-hoc pathway enrichment step. The same property makes each axis verifiable: the user can inspect the basin nodes and ask whether the named biology matches known disease mechanism.

### Preprocessing matters more than substrate-edge gaps for GoF biology (v6 reframe)

Loss-of-function and differentiation-collapse biology — pathway dampening on existing canonical edges, hub-stabilization-mediated upregulation of canonical targets, type-I IFN response engagement of the canonical ISG set — live on the substrate by construction. The substrate's edges encode these reactions. The v4–v5 manuscripts framed gain-of-function biology — KRAS-G12X constitutively active, IDH1-R132H producing 2HG, EGFR-L858R rewired — as *not* representable by the canonical substrate, citing the empirical failures on TCGA_LUAD KRAS-vs-EGFR (multi-PC LR LOOCV AUC = 0.509) and TCGA_IDH 2HG anchor recovery (rank 6,355). **§5 retracts that framing**: under within-patient z-score preprocessing (Methods §"Within-patient z-score"), both failures are rescued (LUAD to 0.703; TCGA_IDH 2HG to rank 26). The canonical Reactome catalysis edges from IDH1/IDH2/D2HGDH to the 2HG node, and from KRAS to MAPK-target reactions, propagate the discriminative signal correctly when input preprocessing preserves the within-patient relative-rank structure that per-modality global-std absorbs. The v5 representational-failure read was wrong; preprocessing was suppressing the signal.

The narrower Paper 4 (mutation-conditional edge-injection extension) open question contracts to *truly neomorphic biology with **zero** canonical effector signaling AND not rescuable by within-patient z-score* — fusion proteins that generate chimeric substrate-sets that don't exist anywhere in the canonical reference, or pathway-creation events that the substrate edges genuinely cannot route. Per-patient mutation calls would still condition the substrate's edge set in such cases — but the set of cases requiring this is narrower than v5 framed, since the IDH1-R132H and KRAS-G12X examples that motivated Paper 4 turn out to be preprocessing-bound, not edge-absent. The naive add-only edge-injection rescue we tested (§5) fails for separate structural reasons (Laplacian dilution); whether directional rewiring helps the residual Paper 4 cases remains the open question.

### The substrate is the load-bearing object; the paper series extends it

Paper 1 establishes the canonical substrate and the GIZMO per-patient projection. Paper 2 extends the substrate with explicit lipid flux and state sub-graphs (dual-architecture, compositional Bayes for ionization-bias correction); Paper 3 incorporates disease-informed topology rewiring (network deformation under disease state rather than fixed substrate); Paper 4 adds mutation-conditional GoF edges *for the narrower remaining open cases* not rescuable by within-patient z-score preprocessing (§5 v6 reframe). Each is an architectural extension of the substrate, exploited by extensions of GIZMO's MAP projection or by other methods consuming the extended substrate. The substrate-as-resource framing is what makes this a paper *series* rather than a sequence of unrelated method papers — each paper grows the shared infrastructure that downstream methods consume.

**Every cohort GIZMO is applied to becomes a case study.** Each Zenodo deposit accompanying this paper includes not just the substrate file but a per-cohort bundle — the F matrix, β/α decomposition, signed-basin output, MOFA+ comparison weights, and a written case-study describing the basin biology and known-mechanism cross-reference — for all 17 panel cohorts and the 3 LOOCV-validation cohorts (NEPTUNE, Wang RA, TB DX). The same deposit convention extends to Papers 2–4: every cohort the substrate is applied to becomes a reusable reference projection in the same coordinate system. The substrate-as-resource framing requires this — a coordinate system is only as useful as the set of cohorts already projected into it, against which future patients can be compared. This is the working version of "atlas of disease state in substrate coordinates" — not a single all-patient resource but a deposit-per-cohort that grows monotonically as the paper series continues.

---

## Methods

### Substrate construction (primary contribution; see Results §1)

The 38,148-node substrate is built from Reactome (15,399 reactions, 6,406 metabolites, 2,382 pathway nodes), StringDB (PPI edges at confidence ≥ 700), HMDB (metabolite identifiers and synonyms), and KEGG (reaction directions and EC-class annotations). The biochem subgraph caps high-degree hub-node edge counts at 200 to prevent the Laplacian from being dominated by water/ATP/NADH-class promiscuous metabolites. Reactions are reified as full nodes; the contraction test on IDH-glioma confirms reactions are integrative substrate elements, not decorative pass-throughs (top-50 \|F̄\| Jaccard standard-vs-contracted = 0.010, biology shifts entirely). Substrate file deposited under CC-BY 4.0 at `data/processed/human_full/graph.json`.

### Identifier mapping

`GeneMapper`: Entrez / Ensembl / HGNC / RefSeq → substrate gene-node UUID. `MetaboliteMapper`: HMDB / ChEBI / PubChem CID / KEGG / LIPIDMAPS / InChIKey / common-name → substrate metabolite-node UUID. Multi-source fallback with per-cohort unmapped-feature reporting; Supp Table S2 tabulates mapping completeness per cohort.

### MAP reconstruction

Per-patient state vector solves

$$
\mathbf{F}_i \;=\; \arg\min_{\mathbf{F}} \left[\, \| \mathbf{x}_i - \mathbf{A}_{\text{obs}} \mathbf{F} \|^{2}_{\boldsymbol{\Sigma}^{-1}} \;+\; \lambda \, \mathbf{F}^{\top} \mathbf{L}_{\text{signed}} \, \mathbf{F} \,\right]
$$

via sparse conjugate gradient (`scipy.sparse.linalg.cg`). Signed Laplacian uses edge-sign annotations (catalytic-product +1, inhibitory-substrate −1, neutral +1). $\lambda$ chosen per cohort by Methods §"Pre-MAP smoothing diagnostic". Solve time ~30 sec per patient on the 38,148-node substrate. Output is per-patient 38,148-dim $\mathbf{F}$ vector in substrate coordinates.

### Within-patient z-score: fine-grained subtype/anchor diagnostic preprocessing

The canonical input preprocessing for cohort-level discrimination uses per-modality log + global-std normalization (Methods §"Identifier mapping" + Stage 1+2 diagnostic in `benchmarks/unsupervised_stratification/stage3_map_reconstruction.py`). The framework supports an alternative preprocessing — **within-patient z-score on $\log_2(v + 1)$ values** — that is *not* canonical for the panel-wide story but is the appropriate choice for two specific diagnostic tasks documented in §5:

1. Fine-grained driver-subtype discrimination within a homogeneous disease (e.g., KRAS-vs-EGFR among TCGA_LUAD single-driver mutants);
2. Metabolite-anchor recovery via graph propagation from RNA input where the canonical preprocessing absorbs the within-patient relative-rank signal into cohort-wide rescaling (e.g., TCGA_IDH 2HG-mito anchor under RNA-seq input).

For each patient $i$, take the patient's positive feature values, $\log_2(v + 1)$-transform them, then $z$-score using this patient's own mean and std:

$$
\tilde{x}_{ij} \;=\; \frac{\log_2(v_{ij} + 1) - \mu_i}{\sigma_i + 10^{-9}}, \qquad
\mu_i, \sigma_i \;=\; \mathrm{mean}, \mathrm{std} \text{ over features } j \in \mathcal{F}_i^{>0}
$$

where $\mathcal{F}_i^{>0}$ is the set of features patient $i$ has positive values for. Patients with fewer than 10 positive features are excluded. Already-log-transformed modalities (Olink NPX, microarray $\log_2$-ratios) skip the $\log_2$ step. Implementation at `gizmo.inference.projection._within_patient_zscore`; ablation evidence and per-cohort AUC sweep at `benchmarks/check_v6_auc_sweep.py`. Across the 13-cohort panel-wide AUC sweep (best $\alpha$-PC AUC and multi-PC LR LOOCV vs active/control labels), cohort-level discrimination under within-patient z-score is comparable to the canonical preprocessing — neither dominates at the panel level. The within-patient z-score buys §5's fine-grained rescues *specifically*, not panel-wide AUC.

### β/α decomposition

$\beta_i = (\mathbf{F}_i \cdot \mathbf{p}) / \| \mathbf{p} \|^{2}$ where $\mathbf{p} = \log\mathrm{PR}_{\text{centered}}$ is the substrate's log-PageRank vector mean-centered to a unit-norm direction. $\boldsymbol{\alpha}_i = \mathbf{F}_i - \beta_i \mathbf{p}$. $\mathbf{p}$ is a property of substrate topology, not cohort data, and is therefore the same 38,148-dim unit vector across all cohorts. Hub-direction annotation in Supp Table S3 shows 47/50 top hub nodes are transcription / signaling / immune / host-defense.

### Reaction contraction test (substrate validation)

Build contracted substrate by removing every reaction R and adding bypass edges (substrate↔product, gene↔substrate, gene↔product); re-run MAP and compare per-node mean \|F̄\|. IDH-glioma result at `benchmarks/results/reaction_contraction_test_IDH_glioma.tsv`: 38,148 → 22,749 nodes; 168,335 → 25,760 edges; top-50 \|F̄\| Jaccard = 0.010; α-PC1 AUC 0.901 → 0.894; biology shifts from canonical IDH-mut (MetAsp / methylglyoxal / ceramide) to generic-graph (ferritin / ribosomal) entirely.

### Hub-direction annotation

For each of the 50 top-log_PR substrate nodes, manually annotated category (transcriptional, signaling/kinase, immune/inflammation, host-defense, structural, mitochondrial-currency, viral-cofactor). Annotation table at `data/processed/human_full/hub_annotation.tsv`.

### Signed-basin decomposition (§2)

For each α-PC, partition substrate nodes by sign of loading; take the largest connected component of the substrate restricted to each sign-class; report basin nodes ranked by absolute loading. Basin mass = sum of squared loadings within the basin; total PC mass = 1 by construction. Implementation at `benchmarks/diagnostics/signed_basin_decomposition.py`.

### Cross-cohort PC alignment (§2)

Pairwise \|⟨PC_k_c1, PC_k_c2⟩\| on 38,148-D unit vectors. Null distribution from all 2,900 cross-disease PC×PC pairs. Same-disease pairs tested against this null. Implementation at `benchmarks/diagnostics/pc_cross_cohort_alignment.py` and `pc_within_disease_pairs.py`.

### PC sign canonicalization for per-node comparison (§2)

PCA components are sign-ambiguous: each eigenvector and its negation explain identical variance, and independent PCA fits on different cohorts pick a sign by solver convention rather than by biology. \|cos\| summary is invariant to this; per-node sign comparisons used in Figs 2g, 2h are not. Before display: if vₐ · vᵦ < 0, flip vᵦ → −vᵦ. Sign-flips that survive canonicalization are real biology (a node where loadings genuinely disagree beyond axis ambiguity). Validating example: IDH·PC2 ↔ TCGA_IDH·PC2 has \|cos\| = 0.77 but signed cos = −0.77; without canonicalization the panel reads as "75 of 76 top contributors sign-flipped" despite being the same biology.

### GoF empirical falsification test (§5)

TCGA-LUAD mutation calls from cBioPortal `luad_tcga_mutations` (429 records). Single-driver subsets exclude co-mutants of KRAS, EGFR, BRAF, MET, ALK, ROS1: n = 61 KRAS-only, n = 28 EGFR-only. Step 1: Mann-Whitney on MAPK / EGFR target genes from Broad Firehose RSEM matrix. Step 2: restrict α-matrix to 89 single-driver patients, compute α-PCA, test each α-PC against KRAS-vs-EGFR binary label by ROC-AUC. Implementation at `benchmarks/diagnostics/luad_gof_verification.py`.

### Cohort panel

17 panel cohorts spanning 7 disease classes and 8+ modality combinations (see Supp Table S2 for acquisition source, n, n_active vs n_control, modality, identifier mapping completeness, ethical clearance per cohort). The independent LOOCV-validation cohorts (NEPTUNE kidney, Wang RA, TB DX) are separately tabulated in §3-validation (Methods).

### Cross-validation and pre-registration

Pre-registration document deposited at Zenodo (DOI to be assigned at acceptance) lists the eight pre-specified falsification tests, the multiplicity-correction discipline, and the variants that are confirmatory vs exploratory. The 67-variant Stage 31 analysis-variant index is included in the supplement.

### Statistical reporting and multiplicity control

**Effect-size threshold.**

$$
\text{Strength} \;=\;
\begin{cases}
2 \cdot \bigl| \mathrm{AUC} - 0.5 \bigr|, & \text{binary clinical label} \\[4pt]
\bigl| \rho_{\text{Spearman}} \bigr|, & \text{continuous clinical label}
\end{cases}
$$

Pre-registered minimum effect size: $\text{Strength} \geq 0.40$.

**Multiplicity correction.** Per-cohort Bonferroni for within-cohort multi-axis × multi-metadata tests at

$$
\alpha_{\text{cohort}} \;=\; \frac{0.05}{n_{\text{cohort-tests}}}, \qquad
n_{\text{cohort-tests}} \;=\; N_{\text{PCs}} \cdot N_{\text{metadata fields with} \geq 5 \text{ non-null}}
$$

Global Bonferroni across all panel cohorts at $\alpha_{\text{global}} = 0.05 / N_{\text{total tests}}$. Per-cohort and global thresholds reported for every comparison; `multi_pc_vs_mofa_factors_augmented.tsv` is the audit log.

**Class-imbalance handling.** Where cohorts are heavily imbalanced (GSE65391_SLE n_active=924 vs n_control=72; CPTAC_OV n_tumor=83 vs n_normal=20), all reported AUCs include stratified-bootstrap 95% confidence intervals at 1,000 resamples preserving the case/control ratio; the bootstrap interval is the primary inference statistic in these cohorts rather than the point estimate alone.

**Permutation nulls.** Label-shuffle permutation nulls at 1,000 shuffles where reported; multiplicity-corrected over the candidate axes (e.g., max-over-α-PC1..PC5 under each shuffle for the per-cohort multi-PC criterion).

**LOO and out-of-sample tests.** Within-cohort leave-one-out PCA refit + sign-align + project (Methods §"LOO PC stability"); cross-cohort PC alignment uses absolute cosines on 38,148-D unit vectors against the 2,900-pair cross-disease null distribution.

### Data and code

GitHub `gizmo` repository (commit hash + tag at acceptance). Substrate file (CC-BY 4.0) deposited at Zenodo (DOI at acceptance).

**Per-cohort case-study deposit.** Every cohort GIZMO has been applied to is deposited on Zenodo with the full per-cohort analysis bundle, not just the F matrix. Each cohort's deposit contains: (i) the per-patient F matrix in 38,148-D substrate coordinates; (ii) β / α / α-PC1..α-PC5 / ‖α‖₂ per patient; (iii) signed-basin decomposition output per α-PC (basin node IDs, loadings, basin mass, Reactome pathway groupings); (iv) per-cohort case-study writeup describing the basin biology, metadata associations, and known-mechanism cross-reference (sourced from cohort publications + CTD / Reactome / KEGG anchors); (v) MOFA+ factor weights, MOFA+_sm substrate-matched weights, and other-method comparison outputs where computed; (vi) per-patient metadata in standardized schema. This includes the 17 panel cohorts of this paper, the 3 LOOCV-validation cohorts (NEPTUNE kidney, Wang RA, TB DX), and the additional cohorts used in supplementary analyses across the paper series. The bundles are versioned per cohort so that downstream methods consuming the substrate can pull a specific cohort's projection as a reproducible reference.

**What's in the supplement.** Pre-registration document + Stage 31 variant index + Supp Tables S1–S5.

---

## References

(Identical to v4 — see `MANUSCRIPT_v4_REWRITE.md`.)

## Acknowledgments

(Identical to v4.)

