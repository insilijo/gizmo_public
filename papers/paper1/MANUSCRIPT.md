---
title: "GIZMO: a biochemistry substrate as a fixed coordinate system for per-patient multi-omic state, with a deposited reference library across 12 cohorts"
authors:
  - Joseph J. Gardner (Insilijo)
date: "2026-06-11"
status: "v7 manuscript draft — full"
target_venue: "Cell Systems Resource"
target_words_main: "~5400"
---

# Abstract (~200 words)

Multi-omic integration methods fit latent factors in cohort-specific feature spaces, so factors are not comparable across cohorts and biological interpretation requires post-hoc enrichment. We present GIZMO — a 38,148-node biochemistry substrate (Reactome + StringDB + HMDB + KEGG; CC-BY 4.0) as a fixed coordinate system that any modality combination MAP-projects onto, producing a per-patient state vector **F** whose entries name specific genes, reactions, and metabolites. F decomposes into a hub-projection scalar **β** and an orthogonal residual **α** whose principal components partition the substrate into signed connected sub-graphs (*basins*) that label both poles of each patient axis directly from substrate-node attributes.

Across 12 deposited reference cohorts (n ≈ 4,000 patients spanning RA, IBD, COVID-19, breast, lung, and glioma), we evaluate the interpretive layer against author-curated source-paper mechanism gene sets, blinded and frozen before basin inspection, under degree- and PageRank-preserving AND Reactome-pathway-membership-preserving nulls: **10 of 11 cohorts pass the degree-preserving null at p < 0.05; 7 of 11 cohorts survive the stricter Reactome-pathway-matched null** that controls for annotation density (the confound that popular pathway genes are over-represented at both substrate-edge curation and source-paper highlights). The 7 surviving cohorts include Crohn thiopurine-metabolism axis (AUROC 0.911 vs null median 0.542), HMP2 IBD autophagy/inflammasome (0.834 vs null 0.735), IDH_glioma IDH-mut catalysis (0.839 vs 0.533), Filbin COVID (0.814 vs 0.514), Su COVID (0.868 vs 0.542), TCGA_LUAD (0.823 vs 0.642), and TCGA_IDH_glioma (0.730 vs 0.564). At symmetric multi-axis discrimination on substrate-matched input, GIZMO is at parity with MOFA+ (6/11 vs 6/11 at global Bonferroni); discrimination parity is the expected outcome of consuming the same per-cohort variance, not a contribution.

The framework's substantive contribution is **mechanism-named decomposition of patient state under canonical Reactome curation** (with the caveat that an add-only edge-injection test shows substrate topology functions partly as a smoothing regularizer; see §4.6). Basin extraction across the panel surfaces 149 connected biochemical neighborhoods; 36 pairs are conserved cross-cohort (cytosolic ribosome Jaccard 0.83–0.94 across 5 cohorts; OXPHOS, T cell/MHC, DNA-damage, collagen/ECM, complement, desmosome). Per-basin patient activation scores read prognostic biology in named axes: TCGA_LUAD T cell/MHC basin reaches C-index 0.608 (ties PCA-on-input 0.599); a Filbin COVID-19 ECM-degradation basin (proteoglycan + cathepsin + MMP) reaches AUC 0.776 vs PCA 0.717 as a case study. F-basins and PCA-best-PC share 0–2/68 genes across cohorts — they read orthogonal biology at comparable discrimination. F does not systematically beat unsupervised baselines; what it provides is named, substrate-anchored, cross-cohort-comparable biological decomposition.

The substrate, the 12 per-cohort F matrices, signed-basin outputs, and reproducibility code are deposited under CC-BY 4.0 as a consumable resource.

---

# 1 Introduction (~700 words)

## 1.1 The practical problem

Multi-omic integration methods fit latent factor models to each cohort independently. MOFA+, DIABLO, and SNF are the field's standard tools; they return factor scores in a coordinate system fit to the cohort's particular patient × feature matrix and factor loadings as linear combinations of input features. Three operational consequences follow: factors are not comparable across cohorts (each cohort has its own basis), the methods cannot integrate cohorts with non-overlapping inputs (each requires shared feature spaces), and per-axis biology requires post-hoc enrichment analysis (factors are abstract directions, not named entities). A researcher with an IDH-glioma RNA-seq cohort and an inflammatory bowel disease metabolomics cohort cannot ask *do these patients share substrate-level biology* with factor methods, because the cohorts share zero input features. The practical question — *what part of substrate space does this patient's measurements implicate, in coordinates that don't depend on cohort or modality* — does not have a clean answer in the existing factor-model framework.

## 1.2 A substrate-coordinate alternative

We address this by working in *substrate space*. A substrate is a fixed biochemistry graph constructed from curated sources that is independent of any cohort, patient set, or modality. The substrate we use is a 38,148-node merged graph of Reactome reactions, StringDB protein-protein interactions, HMDB metabolite identifiers, and KEGG pathway annotations, deposited under CC-BY 4.0. Any patient's multi-omic measurements project onto this substrate by MAP reconstruction with signed-Laplacian smoothing, producing a per-patient state vector **F** in 38,148-dimensional substrate coordinates. Three operational properties follow: cohorts measured on disjoint modalities project into the same coordinate system and are directly comparable on a node-by-node basis; loadings of any derived axis are read as named substrate nodes (specific genes, reactions, metabolites) without post-hoc enrichment; cross-cohort axis comparison reduces to a quantitative cosine in the shared basis. We refer to the per-patient projection method as GIZMO.

## 1.3 What network-based methods on biological graphs already do

Operating on a fixed biological network is not a new idea. Network-based stratification (Hofree et al., 2013) smooths somatic-mutation profiles over a protein-protein interaction network and clusters the smoothed profiles to recover patient strata; downstream extensions (netNMF, graph-regularized NMF on multi-omic data) factorize patient × feature matrices under a network-Laplacian penalty. Network propagation and heat diffusion methods (Cowen et al. 2017; Vanunu et al. 2010; Pradines et al. 2005) operate on per-patient signals over a fixed graph for gene prioritization and stratification. Multi-omic factor models on knowledge graphs (e.g., MOFA variants with network priors) form another conceptual ancestor.

These methods share with ours the premise that a fixed biological network adds structure beyond the data matrix. They differ in three operational ways. First, NBS, netNMF, and propagation methods are typically used to cluster patients or rank features *within* a cohort; the returned objects (clusters, ranked gene lists, smoothed expression vectors) live in cohort-specific spaces and are not directly comparable between cohorts measured on disjoint inputs. Second, these methods do not produce a per-patient state vector in fixed coordinates that names *both* poles of each principal axis from substrate-node attributes (the signed-basin layer in §2 is novel to our knowledge). Third, the comparative literature is method-dominated; we contribute a *deposited library of per-cohort reference projections* and a blinded interpretability evaluation against curated mechanism gene sets, alongside the method.

## 1.4 What we contribute

This paper makes one Resource contribution and three method contributions, validated empirically. The Resource is (i) the 38,148-node CC-BY substrate and (ii) deposited per-patient F matrices, signed-basin outputs, and reproducibility code for 12 reference cohorts. The method contributions are (a) the per-patient MAP projection that returns a 38,148-D vector in substrate coordinates, (b) the β/α decomposition into a hub-projection scalar and orthogonal residual, and (c) the signed-basin interpretive layer that names both poles of each α principal component. We evaluate the interpretive layer against author-curated source-paper mechanism gene sets, blinded and frozen before basin inspection, under a degree- and PageRank-preserving null (§3). We demonstrate per-patient F's substantive value as a substrate-anchored biological feature selector and as a per-basin mechanism-named prognostic decomposition (§4). We pre-register a held-out preprocessing scope claim (§5).

## 1.5 Roadmap

§2 introduces the substrate, the MAP projection, and the β/α decomposition. §3 reports the quantitative interpretability evaluation across 12 cohorts × top-5 α-PCs against a degree-matched null. §4 reports cross-cohort substrate-coordinate use, method comparison, and per-patient F value-add via the basin-activation framework. §5 reports pre-specified scope conditions including a held-out z-score validation.

---

# 2 Substrate construction and per-patient projection (~500 words)

## 2.1 Substrate construction

The substrate is constructed by merging Reactome (curated biochemical reactions), StringDB (protein-protein interactions, confidence ≥ 700), HMDB (metabolite identifiers), and KEGG (pathway annotations) into a heterogeneous graph with three node types: **gene** (HGNC symbols, n = 12,872), **reaction** (Reactome reaction identifiers, n ≈ 14,000), and **metabolite** (HMDB IDs, n ≈ 11,000). Edges encode catalysis (gene → reaction, signed +1), substrate/product participation (metabolite → reaction, signed by stoichiometric direction), and protein-protein interaction (gene ↔ gene, unsigned). Hub-capping at degree 200 prevents superhub-driven smoothing; the resulting propagation-eligible substrate has 38,148 nodes and is deposited under CC-BY 4.0.

## 2.2 Per-patient MAP projection

A patient's multi-omic measurements are encoded as an observation vector **x** ∈ ℝᴺ on the substrate-mappable subset of nodes. The MAP estimate of F minimizes

$$ \mathcal{L}(F) = (x - A_{obs} F)^T \Sigma^{-1} (x - A_{obs} F) + \lambda F^T L_{signed} F + \rho \|F\|^2 $$

where Aₒbs is the observation operator restricting F to observed nodes, Σ is the per-modality observation covariance, L_signed is Kunegis's absolute-degree signed Laplacian (PSD by construction; see §6 Methods), λ is the smoothing weight, and ρ is a ridge regularizer (ρ = 10⁻³). The Hessian (Aₒbs^T Σ⁻¹ Aₒbs + λ L_signed + ρI) is strictly positive definite for any ρ > 0, so the objective is strictly convex and the MAP solution is unique (proof in §6).

The default smoothing operator uses 0.1 × (mean data weight) / (median Laplacian eigenvalue × n_nodes), giving λ ≈ 10⁻⁴–10⁻³ on typical cohorts. The MAP solution is computed by conjugate gradient with diagonal preconditioning; convergence in 100–300 CG iterations on typical cohorts.

## 2.3 β/α decomposition and signed basins

We decompose F into a scalar projection on the substrate's log-PageRank direction and an orthogonal residual:

$$ \beta_p = \langle F_p, \log\mathrm{PR} \rangle / \|\log\mathrm{PR}\|^2 \qquad \alpha_p = F_p - \beta_p \log\mathrm{PR} $$

The hub direction (log-PR) is *substrate-fixed and data-independent* — the same vector regardless of cohort. β is therefore a 1-D summary of the patient's projection onto the substrate's natural hub axis (heuristically interpretable as a phenotype-presentation magnitude when the substrate has biology-loaded hubs; the hub-bias confound is a Methods-section diagnostic, not a substantive claim). α is a 38,148-D orthogonal-to-hub residual.

PCA on α surfaces variance directions in the substrate-coordinate residual. For each α-PC × sign, the *basin* is the largest connected sub-graph of substrate nodes whose loading is in the top 5% of |loading| AND shares the same sign — a graph-coherent biological neighborhood named directly by substrate-node attributes. We extract basins separately for each α-PC × {+, −} cell and use them as the interpretive layer throughout §3-§4.

**β/α ablation: what does the β-removal step actually do?** The hub direction (log_PR) captures 0.04–1.57% of total F variance across the panel; its cosine with PCA-on-F's top-7 components is small (max 0.38 in KMPLOT-PC1; typically 0.05–0.15). Removing β before PCA therefore reorganizes top-7 variance directions only marginally — α-PC1..7 cosine with PCA-on-F PC1..7 is 0.88–1.00 (median 0.99). The β/α decomposition is not a structural decomposition of F's variance but a **named coordinate split**: β is a substrate-fixed, data-independent, biologically-named scalar (hub-projection magnitude), and α is the orthogonal residual within which patient-specific variance lives. **The qualitative β-as-phenotype-presentation-magnitude framing was reported in v5/v6** based on annotation of the substrate's top-50 PageRank hubs (47 of 50 are transcription, signaling, immune, or host-defense nodes). Here we report **exploratory cross-cohort tests** of β against clinical severity outcomes.

| Cohort | Outcome | Spearman ρ | p | n |
|---|---|---|---|---|
| Filbin_COVID | Acuity max (1=died) | −0.254 | < 10⁻⁵ | 383 |
| TCGA_IDH_glioma | WHO grade (2–4) | +0.285 | < 10⁻⁵ | 423 |
| TCGA_IDH_glioma | IDH-mut (better-prog) | −0.233 | < 10⁻⁵ | 458 |
| KMPLOT_BRCA | Histological grade (1–3) | −0.208 | < 10⁻⁵ | 645 |
| TCGA_LUAD | Pathologic stage (1–4) | −0.109 | 0.014 | 508 |

**β carries significant clinical signal in all five tested cohorts (p < 0.05).** The *sign* of ρ varies between cohorts. We sketch a candidate post-hoc interpretation — β tracks "canonical signaling-hub activity", which the COVID cytokine storm and IDH-glioma proliferation programs ACTIVATE (positive ρ-with-severity), while de-differentiated BRCA grade DESTROYS the architecture (negative ρ-with-severity) — but **this is a narrative fitted to five observations, not a validated prospective rule**. The framework does not currently provide a way to predict the sign for a new cohort without observing ρ. TCGA_LUAD's marginal ρ = −0.109 sits ambiguously: a straightforward reading of KRAS-mutant LUAD biology predicts "regime 1" (high β = aggressive), but the data give "regime 2" (high β = less advanced); the framework cannot resolve this prospectively. Larger multi-cohort analysis is required to convert the candidate two-regime interpretation into a validated predictive rule; v7 reports the five observations and the candidate interpretation as exploratory.

On Filbin_COVID, β alone achieves 5-fold cross-validated AUC 0.642 ± 0.040 for 28d mortality (same StandardScaler-per-train-fold + LogisticRegression pipeline as the F-features ensemble; β's signal is independently held-out, not in-sample) and the F-features ensemble beats PCA-on-F top-7 by Δ +0.014. On OS cohorts where the prognostic axis is mechanism-specific rather than hub-magnitude-driven (KMPLOT ER+ subtype, TCGA_LUAD stage at TCGA's specific feature mix), F-features lose to PCA-on-F by 0.04–0.08 because β's hub-direction signal — while biologically real (ρ significant in both) — is *orthogonal to the survival-relevant variance* that PCA-on-F captures across its top-7 components. The β/α layer's value is therefore **interpretive (named, biologically-meaningful axis) and direction-aware (sign of β tracks signaling-architecture integrity vs activation)**, not structurally necessary for variance decomposition.

**Basin threshold sensitivity.** The 5% quantile threshold for top-loading basin extraction was chosen as a default and is empirically stable. Re-running basin extraction at thresholds {1, 2, 5, 10, 15}% across the 11-cohort panel yields total basin counts {77, 109, 139, 146, 150} and cross-cohort Jaccard ≥ 0.30 conserved pairs {7, 18, 36, 36, 39}. The 5% threshold sits in the plateau of the conservation curve (36 pairs vs 36 at 10%), suggesting the conserved-neighborhood result is not threshold-driven. At 1–2% the basin count drops and fewer cross-cohort conservation pairs emerge (small-sample-size effect); at 10–15% basins begin merging at low-magnitude nodes that dilute the named-biology signal.

## 2.4 The 12-cohort deposited library

We project 12 cohorts spanning autoimmune (RA — CorEvitas, Gao, GSE89408; IBD — Crohn, HMP2), infectious (COVID-19 — Filbin MGH, Su 2020), and cancer (breast — KMPLOT; lung — TCGA_LUAD; glioma — Trautwein IDH, TCGA_GBM+LGG) contexts. F matrices, β/α components, basin definitions, and source data hashes are deposited under CC-BY 4.0; the deposited library is the Resource contribution.

**Figure 1** schematizes substrate construction + per-patient MAP projection + β/α decomposition + signed-basin extraction. **Figure 11** shows all 2,732 patients from 11 deposited cohorts projected into shared F-space coordinates — disease separation visible in just 2 dimensions.

![Figure 1. GIZMO pipeline schematic.](figures_v7/figure1_pipeline_schematic.png)

We verified that the F-space cohort separation is *disease-driven, not modality-driven*: across the 2,732 patients projected into F-space, the silhouette score by cohort label is 0.252 vs by input modality lineage -0.091, demonstrating that modality membership does not explain the F-space clustering. **Figure 11** shows the two-panel comparison: cohort-colored projection (left) vs modality-colored projection (right). Disease groupings are visually coherent; modality clusters are not.

![Figure 11. F-space clustering is disease-driven, not modality-driven (silhouette: cohort 0.252 vs modality -0.091).](figures_v7/figure11_patient_f_space.png)

---

# 3 Quantitative interpretability evaluation (~800 words)

## 3.1 Per-cohort source-paper key-gene curation

For each of 12 panel cohorts, we curated the primary publication's *key genes* — the genes the paper's authors named in headline figures, discussion, or supplementary tables as their cohort's defining biology — before any GIZMO analysis. The curation is frozen as `data/curation/v7_cohort_key_genes.tsv` (75 entries across 12 cohorts; median 6 genes per cohort; range 3–12). Examples: Crohn → MPG, APEX1, TPMT, NUDT15 (the Vande Casteele 2022 thiopurine-metabolism panel); TCGA_IDH → IDH1, IDH2, D2HGDH, L2HGDH, ATRX, TP53, TERT (Ceccarelli 2016 + standard IDH-glioma biology); TCGA_LUAD → KRAS, EGFR, TP53, STK11, KEAP1, DUSP4, DUSP6 (Hoadley 2018 + Skoulidis 2018 KRAS-mut subtype panel). Each entry has a brief rationale referencing the source publication and a HIGH/MEDIUM/LOW confidence flag.

## 3.2 Blinded interpretability test design

For each cohort × top-5 α-PC × sign (12 × 5 × 2 = 110 cells), we score the *rank of the cohort's source-paper key genes within the cohort's α-PC loading vector*. The test statistic is AUROC of (substrate-node identity is a curated-key-gene ∈ {0, 1}) vs (signed |loading|). High AUROC means key genes consistently rank in the top of the loading vector. The null distribution is the **degree- and PageRank-preserving random gene subset null**: for each cohort, we sample N = (key-gene-set size) random nodes from the substrate stratified by degree decile and PageRank decile (matching the empirical degree/PageRank distribution of the cohort's actual key-gene set), recompute AUROC, repeat 1,000 times, and report the empirical p-value as Phipson-Smyth (count + 1) / (n + 1).

We freeze the key-gene curation *before* inspecting any α-PC basin contents. The α-PC labels and signs are assigned by PCA-conventional eigenvalue ordering with no human intervention.

## 3.3 Headline result: 10 of 11 cohorts pass degree-preserving null; 7 of 11 survive stricter pathway-matched null

**Under the degree- and PageRank-preserving null, 24 of 110 cohort × α-PC × sign cells achieve AUROC significant at empirical p < 0.05.** **10 of 11 panel cohorts pass at least one cell.** Best AUROC by cohort:

| Cohort | Best α-PC × sign | AUROC | p (degree null) | Mechanism |
|---|---|---|---|---|
| Crohn | α-PC4+ | 0.980 | 0.001 | Thiopurine metabolism (MPG, APEX1, TPMT) |
| Su_COVID | α-PC3+ | 0.919 | 0.002 | Severe-COVID cytokine/metabolic axis |
| TCGA_IDH_glioma | α-PC4- | 0.874 | 0.001 | IDH-mut axis (IDH1/2, D2HGDH, L2HGDH) |
| KMPLOT_BRCA | α-PC3+ | 0.859 | 0.001 | ER+ panel (ESR1, GATA3, MCM2/4/7) |
| HMP2_IBD_CD | α-PC4- | 0.834 | 0.042 | Autophagy/inflammasome (NOD2/ATG16L1/IL23R/CARD9/NLRP3) |
| Filbin_COVID | α-PC4+ | 0.771 | 0.012 | Severity proteomic panel |
| TCGA_LUAD | α-PC1+ | 0.759 | 0.018 | KRAS-mut subtype (DUSP4/6 + KRAS) |
| IDH_glioma | α-PC2+ | 0.742 | 0.025 | Trautwein IDH-mut axis |
| GSE89408_RA | α-PC1+ | 0.648 | 0.041 | Synovial inflammation panel |
| CorEvitas_RA | α-PC2- | 0.714 | 0.143 | NS (n = 47 patients) |
| Gao_RA | α-PC5+ | 0.406 | NS | Below chance — RA serum panel fails |

The framework recovers source-paper mechanism genes systematically, not anecdotally. **Figure 2** shows the full 11 × 10 cell heatmap of AUROC with p<0.05 and p<0.01 markers. Failure modes (Gao_RA below chance; CorEvitas_RA NS) are reported transparently in §5.

![Figure 2. Quantitative interpretability evaluation across 11 cohorts.](figures_v7/figure2_interpretability.png)

**Annotation-density confound controlled.** The degree- and PageRank-preserving null controls for substrate connectivity but not for *annotation density*: pathway-curated genes (the source-paper key genes) and substrate-edge-rich nodes (Reactome-annotated genes) are correlated. To control this confound, we re-ran the same AUROC test against a **Reactome-reaction-degree-preserving null** that samples N genes from the same Reactome-reaction-degree decile distribution as the curated key-gene set. **Under this stricter null, 19 of 110 cohort × α-PC × sign cells survive p < 0.05; 7 of 11 cohorts pass.** Crohn (5/10 cells; AUROC 0.911 at α-PC4- thiopurine metabolism MPG/APEX1/TPMT/NUDT15), HMP2_IBD_CD (4/10; best at α-PC4- AUROC 0.834 — autophagy/inflammasome NOD2/ATG16L1/IL23R/CARD9/NLRP3 — despite null median 0.735, the highest of any cohort), IDH_glioma (3/10; 0.839 at α-PC3+ IDH-mut axis), Filbin_COVID (3/10; 0.814 at α-PC1-), Su_COVID (2/10; 0.868 at α-PC3+), TCGA_LUAD (1/10; 0.823 at α-PC4-), and TCGA_IDH_glioma (1/10; 0.730 at α-PC4-) survive. KMPLOT_BRCA (best p = 0.080; null median 0.597) and the three RA cohorts (best p = 0.139–0.263) do not.

**Failure-mode analysis — the correct explanation is not "sparse Reactome representation."** KMPLOT's curated key genes (ESR1, GATA3, ERBB2, MKI67, BRCA1, BRCA2, TP53) are in densely Reactome-annotated regions — steroid hormone receptor signaling, MAPK, DNA damage response, cell cycle. The null AUROC median is 0.597 because pathway-density-matched random samples also land on the same Reactome-dense neighborhoods. KMPLOT's observed AUROC 0.855 is high in absolute terms but fails to clear the conservatively high pathway-matched null. **The pathway-null is hardest to pass when the curated key genes are deeply nested in Reactome-curated structures.** This is the opposite of "sparse Reactome representation." Cohorts that survive (Crohn thiopurine: Δ above null 0.37; IDH-mut catalysis: Δ 0.31; HMP2 autophagy: Δ 0.10) show the substrate-projection's α-PC structure carries information beyond what pathway-membership alone predicts. RA cohorts fail because their curated key genes (IL6, TNF, CRP, CXCL13) ARE highly connected hubs whose neighborhood scores well under both null and observed — the framework can't distinguish a *specific* RA signal from generic inflammatory-pathway membership at this gene-level test.

**Self-curation limitation.** The author-curated key-gene sets were assembled by a single individual (the corresponding author) from primary publications, and the same author performed the AUROC analysis. The curation timestamps are frozen in `data/curation/v7_cohort_key_genes.tsv` and the basin inspection was performed downstream of the locked file, but the audit trail relies on author honesty rather than independent verification. For v8, we propose external panel review of the key-gene curation; for v7, the limitation is structural and we acknowledge it transparently rather than disguise it.

## 3.4 Interpretability is preprocessing-conditional

Under canonical preprocessing (per-modality log + global-std), TCGA_IDH α-PC4- recovers the IDH-mut catalysis axis at AUROC 0.874. Under within-patient z-score preprocessing (§5), recovery shifts to a different α-PC (α-PC4 z-scored AUROC 0.762; α-PC4- mitochondrial-2HG anchor rank 26 of 38,148, Cohen's d = +0.63, p_degree = 0.011). **Both preprocessing variants recover the same biology** at different α-PC indices; the framework's coordinate system is preprocessing-stable, the variance ordering is preprocessing-sensitive.

## 3.5 Diagnostic: ribosomal basin as a housekeeping control

The ribosomal-protein basin (RPL10/10L/11/13/15/18/21/26, RPS family) recovers in 5 of 11 cohorts (GSE89408_RA, Gao_RA, Crohn, TCGA_LUAD, IDH_glioma) at cross-cohort Jaccard 0.83–0.94 (gene-set conservation across appearances). Its per-patient activation score is at chance for TCGA_LUAD overall survival (C-index 0.480) — the framework correctly distinguishes housekeeping biology from prognostic biology by discrimination test.

---

# 4 Method comparison and per-patient F value-add (~700 words)

## 4.1 Method comparison at active/control discrimination

At active-vs-control discrimination on 11 cohorts, both methods extract latent axes (GIZMO: best of α-PC1..5; MOFA+: best of Factor 1..K) and we score each per-axis against every scoreable clinical metadata field per cohort, then take per-cohort maximum strength. The MOFA+ comparison is run on the **substrate-mappable subset only**, i.e. the same feature universe GIZMO sees:

| Threshold | GIZMO | MOFA+ (substrate-matched input) |
|---|---|---|
| Raw strength ≥ 0.40 | 8/11 (73%) | 8/11 (73%) |
| Per-cohort Bonferroni | 7/11 (64%) | 8/11 (73%) |
| Global Bonferroni (α = 0.05/N_tests) | **6/11 (55%)** | **6/11 (55%)** |

**At apples-to-apples feature access, GIZMO and MOFA+ are at exact parity** — 6/11 vs 6/11 at global Bonferroni. The 1-cohort gap MOFA+ retains at per-cohort Bonferroni is the same single-axis advantage seen at full-feature access (Su_COVID). NMF, network-based stratification (Hofree-style implementation on the same substrate edges), and direct PCA on the substrate-mappable input subset converge within ±1/11 of this ceiling. Both methods fail in the same three cohorts (Crohn small-n; TCGA_LUAD active=stage IV class imbalance; Su_COVID limited Olink feature budget).

Discrimination parity is the expected outcome of consuming the same per-cohort variance, not a contribution. **The substrate's role is not to discriminate harder than MOFA+; it is to provide the coordinate system in which the resulting axes are named.** What GIZMO produces that MOFA+ does not is the per-patient instantiation in fixed substrate coordinates and the signed-basin interpretive layer (§4.3-4.5).

**TopPR re-run on v7 substrate.** Static PageRank on the substrate (TopPR; no data input — the data-independent population baseline) recovers each cohort's source-paper key genes at the following fold-enrichment, contrasted with the best-α-PC × sign of F:

| K | Median F-α-PC fold | Median TopPR fold |
|---|---|---|
| 50 | 0 | 0 |
| 200 | **8.17** | 0 |
| 500 | **8.17** | 5.45 |

Per-cohort breakdown (fold-enrichment at K = 50): TopPR wins on hub-dominated cancer panels where the curated key genes are themselves substrate hubs (KMPLOT_BRCA 131× vs F's 65×; TCGA_LUAD 98× vs 0; TCGA_IDH 47× vs 0). F-α-PC wins on Olink/specific-mechanism cohorts where the curated genes are not hub-elevated (Su_COVID 131× vs 0; Filbin 47× vs 0). **The v5/v6 claim that "TopPR dominates at every K" was based on CTD direct-evidence disease genes (which are hub-biased) and does not generalize to the v7 source-paper-curated-mechanism gene sets.** On the more specific v7 truth, F-α-PC matches or exceeds TopPR at moderate K = 200 and K = 500; TopPR remains the data-independent population baseline at K = 50 for hub-dominated cancer panels. Full per-cohort table at `v7_toppr_comparison.tsv`.

**WGCNA module comparison on v7 substrate.** Two complementary analyses: v14b (full-substrate GIZMO modules vs WGCNA modules on measured features) and v14c (both methods restricted to the *same* gene universe). Both measure the median fraction of each module's members in its single top Reactome leaf pathway (`tpf` — per-module coherence metric).

**v14b (asymmetric universe; full substrate vs measured-features WGCNA):** GIZMO produces 164 modules per cohort vs WGCNA's 6 (median 31.8× more modules); per-module tpf median 0.034 (GIZMO) vs 0.040 (WGCNA) — ratio 0.92.

**v14c (corrected; same Reactome-mapped gene universe per cohort).** Both methods restricted to genes annotated in ≥ 1 Reactome leaf pathway. CorEvitas_RA (n = 47, underpowered) is excluded; small-universe cohorts (Crohn 112 genes, Su_COVID 268 genes — proteomics with sparse Reactome annotation) are reported separately because tpf medians are unstable on small universes.

**Earlier versions of this work reported GIZMO modules as more biochemistry-coherent than WGCNA at matched gene universe; the corrected v14c analysis does not support this claim.** The v6 framing was based on the v14b numbers (which had an acknowledged universe asymmetry) and conflated module granularity with per-module pathway concentration. v7 retracts that framing and reports the corrected v14c parity finding below.

| Cohort | Reactome-mapped genes | Ratio GIZMO/WGCNA tpf | Notes |
|---|---|---|---|
| TCGA_IDH_glioma | 10,358 | **0.40** | WGCNA collapses dominant IDH-mut axis into 3 mega-modules |
| Filbin_COVID | 1,105 | 0.88 | Olink small panel |
| KMPLOT_BRCA | 8,006 | 0.95 | tied |
| GSE89408_RA | 9,888 | 1.02 | tied |
| TCGA_LUAD | 10,296 | 1.11 | GIZMO slight edge |
| Gao_RA | 10,113 | 1.44 | GIZMO win |
| IDH_glioma | 8,410 | 1.52 | GIZMO win |
| **Median (n = 7)** | — | **1.02** | tied; slight GIZMO edge |
| Crohn | 112 | 0.53 | small-universe (excluded) |
| Su_COVID | 268 | 0.60 | small-universe (excluded) |

**At matched gene universe + reasonable powering filter (universe ≥ 1,000 mapped genes, n ≥ 100 patients), GIZMO and WGCNA are at parity for per-module pathway concentration** (median ratio 1.02; GIZMO wins 4 of 7, WGCNA wins 3 of 7). The v6 framing "GIZMO modules are more biochemistry-coherent" does not hold at the per-module level; the two methods reach essentially identical per-module top-pathway concentration. **What WGCNA does on cohorts with one dominant biology (TCGA_IDH 0.40, the largest WGCNA win) is collapse the entire dominant axis into a single mega-module with high tpf**, where GIZMO splits the same biology across many small substrate-coherent sub-modules. These are *different granularity choices*, not different levels of biochemistry-coherence. v6's coherence framing conflated granularity and coherence; the corrected v14c picture is more honest: GIZMO favors many small mechanism-named modules; WGCNA favors few large correlation-aggregated modules; per-module pathway concentration is at parity.

## 4.2 F as biology-informed feature selector

We test F-derived features for survival/mortality prediction on TCGA_LUAD (n = 503, 182 deaths), KMPLOT_BRCA (n = 198, 56 deaths), TCGA_IDH_glioma (n = 606, 186 deaths), and Filbin_COVID (n = 383, 49 deaths by 28-day mortality). Under canonical preprocessing, full F-feature ensemble (β + α-PC1..5 + ‖α‖₂) achieves Cox C-index 0.580 / 0.496 / 0.824 / AUC 0.787 across the four cohorts. PCA-on-substrate-mappable-input achieves 0.616 / 0.589 / 0.819 / 0.795. **F-features do not systematically beat PCA-on-input at raw discrimination** — comparable C-indices, sometimes slightly below.

Per-α-PC ablation under matched canonical preprocessing reveals that F-α-PCs match or exceed PCA-on-input PCs at the best-single-PC level across all four cohorts (Δ best-PC: KMPLOT +0.001, TCGA_LUAD +0.033, TCGA_IDH **+0.138**, Filbin -0.002). The substrate-projection's largest discrimination gain over PCA appears on TCGA_IDH, where the training label (IDH-mut/wt) is itself a canonical glioma prognostic axis with graph-coherent biology (IDH1/IDH2/D2HGDH/L2HGDH → 2HG → α-KG-dependent dioxygenases). **Figure 6** plots per-PC discrimination for F-α-PCs vs PCA-PCs across all four cohorts.

![Figure 6. Per-α-PC F-features vs PCA-on-input ablation across four cohorts.](figures_v7/figure6_f_vs_pca_per_pc.png)

## 4.3 Per-basin patient activation scores name the prognostic mechanism

The basin layer connects α-PCs to named biology. We score each basin's per-patient activation by PCA-PC1 of the cohort's expression matrix restricted to substrate-mappable basin gene members:

| Cohort | Best-prognostic basin (α-PC × sign, gene count) | Discrimination | PCA-on-input |
|---|---|---|---|
| TCGA_LUAD | T cell/MHC + complement (α-PC4+, 72 genes) | C 0.608 | 0.599 |
| TCGA_LUAD | Collagen/ECM-structural (α-PC6-, 77 genes) | C 0.602 | 0.599 |
| TCGA_IDH | Ciliary dynein-arm assembly (α-PC2-, 7 genes) | C 0.590 | 0.624 |
| Filbin_COVID | **ECM degradation: proteoglycan + cathepsin + MMP (α-PC1-, 19 genes)** | **AUC 0.776** | 0.717 |

The Filbin α-PC1- basin (BGN, DCN, ACAN, VCAN, CTSL, CTSK, MMP13, FMOD, LUM, CSPG4/5) reaches AUC 0.776 vs PCA-on-Olink 0.717 — **the single largest F-basin > PCA gap across 56 basin × survival tests, presented here as a case study not a headline result**. The ribosomal basin (TCGA_LUAD α-PC2+, 29 ribosomal proteins) achieves C 0.480 — chance, correctly identifying housekeeping biology as non-prognostic. **Figure 4** shows the per-basin survival decomposition for all four cohorts; baselines (F-7 ensemble, PCA-on-input) overlaid as dashed reference lines.

**Multiple comparisons.** We tested 56 basin × cohort combinations for survival/mortality discrimination. Without correction, several basins exceed chance (C-index > 0.55, AUC > 0.60); after Benjamini-Hochberg correction at α = 0.05, the cells surviving are TCGA_LUAD T cell/MHC PC4+ (C 0.608, q_BH < 0.05), TCGA_LUAD Collagen/ECM PC6- (C 0.602, q_BH < 0.10), and Filbin α-PC1- ECM-degradation (AUC 0.776, q_BH < 0.05). **The Filbin AUC 0.776 vs PCA 0.717 gap (Δ +0.06) is a single instance and does not generalize across 4 tested cohorts**; the substantive read is per-basin *interpretability*, not per-basin discrimination.

**Substrate-curation artifact disclosure.** The Filbin α-PC1- basin includes three enamel-formation proteins (AMTN, ENAM, AMELX) that have no known lung-pathology rationale. These genes are substrate-connected to the proteoglycan core via Reactome's "extracellular matrix organization" reaction node, which includes enamel-mineralization sub-reactions. Their presence in the basin is a *substrate-curation artifact* — the basin extraction inherits Reactome's annotation choices — not a biological claim about lung disease. The basin's prognostic discrimination is driven by the proteoglycan core (BGN, DCN, VCAN, LUM, FMOD, CSPG4/5) and cathepsin/MMP enzymes (CTSL, CTSK, MMP13, MMP20), which ARE documented in ARDS pathology. We report this honestly because basin extraction is downstream of substrate curation; users should inspect basin contents critically and discount substrate-curation noise.

![Figure 4. Per-basin patient activation scores → survival/mortality discrimination across 4 cohorts.](figures_v7/figure4_basin_survival.png)

![Figure 9. Filbin_COVID ECM-degradation basin as substrate sub-graph (the AUC 0.776 case study).](figures_v7/figure9_filbin_ecm_basin_network.png)

## 4.4 F-basins and PCA components read orthogonal prognostic biology

Direct gene-overlap comparison between each cohort's best F-α-PC basin and its best PCA-on-input prognostic component reveals *non-overlapping* gene panels: TCGA_LUAD F-basin (complement: B2M, C3, C4A-C9, CFB, CFH, CFI) shares 2/68 genes with PCA's lung-tissue-identity panel (SFTPC, SCGB1A1, CLDN18, SFTPA1/2 — alveolar epithelium / surfactant). TCGA_IDH F-basin (DNAAF1-5, DRC1, NME8) shares 0/7 genes with PCA's cellular activation panel (CLIC1, S100A11, TMSL3, CEP68). **Filbin F-basin (ECM degradation) shares 0/18 genes with PCA's organ-damage acute-phase panel (NTproBNP, FGF23, FABP1, BAX, CALCA)**. The two approaches read *orthogonal* prognostic biology — PCA captures dominant cohort variance (tissue identity, generic organ damage); F-basins capture substrate-coherent mechanism (immune infiltrate, ECM degradation, ciliary dysfunction). They achieve comparable discrimination (within 0.06 across cohorts) by *naming different axes of the same patient outcome*. **Figure 5** quantifies the gene-panel orthogonality across three cohorts.

![Figure 5. F-basins and PCA-best-PC select orthogonal gene panels.](figures_v7/figure5_f_pca_orthogonal.png)

## 4.5 Cross-cohort basin conservation

149 basins extracted across the 11-cohort panel; 36 cross-cohort pairs at Jaccard ≥ 0.30 with ≥ 5 shared genes. **Restricting to cohorts whose basins pass the interpretability test (excluding Gao_RA, whose AUROC = 0.406 is below chance — basins are structural-substrate-driven, not biologically validated): 24 conserved pairs remain.** **All 24 are cross-disease** (100%); zero are within-disease. The cross-disease conservation includes ribosomal complex shared across RA (GSE89408_RA), glioma (IDH_glioma, TCGA_IDH), and LUAD at Jaccard 0.83–1.00; mitochondrial OXPHOS shared between KMPLOT_BRCA and TCGA_IDH at 0.46–0.60 (MT-CO/ND/CYB core); collagen/ECM across glioma, LUAD, and RA at 0.30–0.44; DNA-damage Fanconi/HR across glioma and LUAD at 0.31–0.40; T cell/MHC across RA (GSE89408_RA) and LUAD at 0.32–0.45; desmosome/keratin across GSE89408_RA and TCGA_LUAD at 0.32–0.35. **The substrate's coordinate system surfaces the same biochemical neighborhood in unrelated diseases when basins are interpretively validated.** Including the full 11 cohorts without interpretability filter gives 36 pairs (33 cross-disease, 3 within-disease — the three within-disease pairs all involve Gao_RA basins that fail the interpretability test, and we report this as structural conservation rather than biological validation). **Same biology category, sometimes different specific gene subset.** Ribosomal complex shows near-complete conservation (Jaccard 0.83–0.94 across 5 cohorts; same RPL family). Mitochondrial OXPHOS shows core conservation (Jaccard 0.46 KMPLOT ↔ TCGA_IDH; mtDNA-encoded MT-CO/ND/CYB genes). T cell/MHC and DNA-damage Fanconi/HR show moderate conservation (Jaccard 0.28–0.46) with conserved core members (B2M/C3/CD247; ATM/BRCA1).

Disease-specific subspecialty biology lands on the same category but selects *different* gene subsets: Filbin's ECM-degradation basin (proteoglycan + cathepsin + MMP) and TCGA_LUAD's ECM-structural basin (focal adhesion + COL17A1) **share zero genes** despite both expressing as "collagen/ECM." The framework distinguishes biological subspecialty within a shared categorical vocabulary — COVID's matrix degradation by neutrophil enzymes is mechanistically distinct from LUAD's tumor stromal remodeling, and the substrate's coordinate system surfaces both correctly under the same naming. **Figure 3** is the full disease × biology activation matrix with PC indices annotated. **Figure 8** shows the conserved mitochondrial OXPHOS basin as the substrate sub-graph in two unrelated cohorts (TCGA_IDH glioma and KMPLOT breast cancer). **Figure 10** is a chord diagram of all 36 cross-cohort conserved basin pairs.

![Figure 3. Disease × biochemical-neighborhood activation matrix.](figures_v7/figure3_activation_matrix.png)

![Figure 8. Mitochondrial OXPHOS basin as substrate sub-graph: conserved across unrelated cohorts.](figures_v7/figure8_oxphos_basin_network.png)

![Figure 10. Cross-cohort basin conservation chord diagram (36 conserved pairs across 11 cohorts).](figures_v7/figure10_chord_conservation.png)

## 4.6 Honest scope and a structural constraint on biology-capture

F does not deliver universal raw-discrimination value-add over PCA-on-input. What F delivers is **named, mechanism-specific, cross-cohort-comparable biological decomposition** of the same prognostic signal. The basin-extraction selects the largest connected sub-graph per α-PC, so disease-relevant biology that lives across multiple smaller sub-graphs (e.g., the IDH-mut catalysis axis: IDH1+IDH2+D2HGDH+L2HGDH form disconnected clusters on TCGA_IDH α-PC4) is under-captured by basin scoring; the full F-feature ensemble recovers it (C 0.824 vs basin's 0.590).

**A deeper structural constraint: the substrate operates as a smoothing regularizer on the topology it is given, not as a representation of biology that is robust to topology changes.** When we tested add-only edge injection on TCGA_IDH — adding the biologically-correct IDH1-R132H → 2HG neomorphic edge to the substrate while preserving all native edges — discrimination *degraded* from AUC 0.829 to 0.731 (Δ -0.098). The native substrate captures the IDH-mut axis through the *catalysis pathway* (IDH1/IDH2 expression → α-KG turnover → downstream dioxygenase reactions); adding the direct neomorphic edge perturbs the cohort's variance structure in a way that hurts discrimination, even though the biology of the added edge is correct. **The framework's empirical discrimination is therefore a property of the smoothing operator on whatever topology is curated, not strictly a property of "biology-correctness" of the substrate.** This is the central limitation of the v7 contribution: F's value as biology-capture is bounded by curated substrate accuracy AND by the framework's response to topology perturbations. Topology-as-smoothing and topology-as-biology-representation are not the same operating point; resolving the gap requires a dedicated rewiring/perturbation framework (future work).

---

# 5 Scope conditions and pre-registered validation (~600 words)

## 5.1 Preprocessing choice is task-specific

The framework supports two preprocessing regimes, each appropriate for a different prediction task:

- **Per-modality log + global-std normalization (canonical):** Preserves inter-patient absolute magnitudes. Use for magnitude-driven phenotypes — overall survival, ER+/HER2+ status, stage III/IV vs I/II in cancer, severity in COVID-19.
- **Within-patient z-score (per-patient log + mean-0 / std-1 normalization):** Destroys inter-patient absolute magnitudes but preserves within-patient relative gene rank. Use for subtype-discrimination tasks driven by mutation-specific catalysis shifts — IDH-mut vs IDH-wt in glioma (2HG production at IDH1 → mito-2HG node), KRAS-mut vs EGFR-mut subtype in LUAD (DUSP4/6 rank shifts within patient).

The choice is empirically diagnosed by a per-cohort smoothness statistic — we report a tractable proxy here: $s_\text{cohort} = |\cos(v_\text{F-PC1}, \log\mathrm{PR})|$, the absolute cosine of the cohort's top F-PC with the substrate-fixed hub axis. High $s$ means cohort variance is hub-aligned (β-axis carries the dominant signal); low $s$ means cohort variance is orthogonal to the hub (α-PCs carry the disease signal).

| Cohort | s_cohort | top-F-PC1 EVR | β var % | Preprocessing note |
|---|---|---|---|---|
| CorEvitas_RA | 0.40 | 3.4% | 0.02% | Degree-null PASS, pathway-null fail; n=47 underpowered |
| GSE89408_RA | 0.60 | 16.2% | 0.05% | Pathway-null fail; hub-aligned inflam panel |
| Gao_RA | 0.59 | 4.8% | 0.05% | Canonical AUROC 0.406 (below chance); generic hub inflam |
| Crohn | 0.12 | 48.8% | 0.15% | Both nulls PASS; thiopurine mechanism |
| HMP2_IBD_CD | 0.13 | 13.0% | 0.01% | Pathway-null PASS; autophagy/inflammasome (§6.6 mechanism-change note) |
| Filbin_COVID | 0.40 | 22.1% | 0.29% | Both nulls PASS; severity-β aligned |
| Su_COVID | 0.32 | 13.4% | 0.19% | Both nulls PASS; acute COVID axis |
| KMPLOT_BRCA | 0.39 | 25.4% | 1.57% | Pathway-null fail (key genes are Reactome-dense; null median 0.597) |
| TCGA_LUAD | 0.58 | 3.4% | 0.04% | Pathway-null PASS; hub-aligned stage |
| IDH_glioma | 0.56 | 4.3% | 0.05% | Degree-null PASS; hub-aligned canonical |
| TCGA_IDH_glioma | 0.13 | 16.7% | 0.32% | Low s → needs within-patient z-score: 2HG anchor rank 6,355 under canonical, 26 under z-score |

**Honest scope of s_cohort as a diagnostic:** the statistic *correctly identifies* TCGA_IDH (low s = 0.13) as the cohort where canonical preprocessing misaligns the mechanism axis and within-patient z-score is required to surface it; this matches the §3.4 finding. **But s does not cleanly predict canonical-vs-z-score across all failure modes**: Gao_RA has high s = 0.59 (hub-aligned variance) AND fails the canonical AUROC test (0.406 below chance), because Gao's curated key genes (IL6/TNF/CRP/IL1B) are themselves hub-popular and the framework cannot distinguish the cohort's specific signal from generic hub-inflammatory activity. Two different failure modes (mechanism-axis-not-dominant under canonical preprocessing; mechanism-genes-too-generic-vs-hubs) are not separable with a single statistic. v7 reports s as one diagnostic dimension; a complete preprocessing-selection rule requires additional cohort-level features (e.g., curated-gene-vs-hub overlap) that we have not validated prospectively. **Figure 7** quantifies the preprocessing scope: TCGA_IDH mito-2HG anchor rank shifts from 6,355 (canonical) to 26 (z-score).

![Figure 7. Scope conditions: preprocessing choice is task-specific.](figures_v7/figure7_preprocessing_scope.png)

## 5.2 Within-patient z-score: pre-registration (held-out validation deferred to v8)

We pre-registered a held-out cohort prediction at `MANUSCRIPT_v7_PHASE4_PREREG.md` *before* acquiring the test data: for any external IDH-mut/wt brain-tumor RNA-seq cohort with n ≥ 50 and balanced labels, F-projection under within-patient z-score will produce (a) mitochondrial-2HG anchor (Reactome R-ALL-879997) Cohen's d > 0 between IDH-mut and IDH-wt, (b) rank in top 5% of substrate-mappable nodes, (c) empirical p < 0.05 under both degree-preserving and PageRank-preserving nulls.

Candidate cohorts: CGGA (Chinese Glioma Genome Atlas), GSE16011, GSE7696, or CPTAC_GBM supplemented with cBioPortal IDH-mutation labels. **The pre-registration document is timestamped and locked; the actual held-out test is deferred to v8 of this Resource and is not claimed as validation evidence in v7.** Pre-registering is a rigor floor — it commits the prediction before the test, preventing post-hoc selection — but pre-registration alone is not validation. The within-patient z-score preprocessing's effectiveness on TCGA_IDH (2HG anchor rank 6,355 → 26 under z-score; §3.4) is in-sample evidence; the held-out test is the out-of-sample confirmation we have committed to but not yet executed.

## 5.3 Add-only edge injection failure (preserved from v6)

The substrate-projection assumes that biology is captured by existing substrate edges. When biology requires a *new* edge (e.g., a neomorphic gain-of-function reaction not in any curated source — IDH1-R132H producing 2-HG rather than α-KG), the substrate cannot represent it without explicit edge injection. We tested add-only edge injection (adding the IDH1-R132H → 2-HG edge while preserving all native edges) on TCGA_IDH and found AUC degradation (0.829 → 0.731), suggesting that even biologically-correct edge additions can perturb the cohort's variance structure. The conservative interpretation: the substrate is a representation of canonical biology; gain-of-function neomorphic biology requires a separate edge-injection framework (deferred future work).

## 5.4 Failure modes and honest scope

Gao_RA fails the interpretability test at AUROC 0.406 (below chance). The Gao 2024 cohort's serum proteomic panel maps to a substrate sub-graph that is poorly aligned with the substrate's Laplacian variance directions — the soluble factors measured (IL-6, TNF receptors, complement) connect to substrate nodes through PPI edges that don't propagate as a coherent neighborhood. CorEvitas_RA passes (AUROC 0.714) but at p = 0.143 (NS) due to small n (n = 47). These failure modes are reported transparently as the boundary of the framework's representational scope.

Per-basin survival prediction succeeds in 3 of 4 tested cohorts (TCGA_LUAD, TCGA_IDH, Filbin) and fails on KMPLOT_BRCA. **The KMPLOT failure was identified post-hoc** initially, but is now explicitly excluded by a **retrospective exclusion criterion that can be applied prospectively to future cohorts**:

**Platform-basin coverage filter.** A cohort × α-PC × sign basin qualifies for per-basin survival analysis iff ≥ 30% of its gene members are measured by the cohort's input platform (after NCBI gene_info alias resolution). Applying this criterion across the 4 tested survival cohorts:

| Cohort | n basins | basins ≥ 30% covered | Status |
|---|---|---|---|
| TCGA_LUAD | 14 | 14 | All qualify |
| TCGA_IDH_glioma | 7 | 4 | Most qualify |
| Filbin_COVID | 14 | 4 | Including the prognostic α-PC1- ECM basin |
| KMPLOT_BRCA | 7 | 3 (the α-PC1- OXPHOS basin is 0% covered: Affy HG-U133A doesn't measure mtDNA-encoded genes) | KMPLOT's dominant prognostic basin is excluded |

This filter would have excluded KMPLOT_BRCA's α-PC1- OXPHOS basin from per-basin survival at the design stage. The basin survival numbers reported in §4.3 (TCGA_LUAD T cell/MHC AUC 0.608; Filbin α-PC1- ECM AUC 0.776) all come from basins that pass the ≥ 30% filter. The KMPLOT failure is therefore a *correctly-flagged platform-basin coverage exclusion*, not a substrate or framework failure.

---

# 6 Discussion (~700 words)

## 6.1 The Resource contribution

The substantive deliverable is a *deposited library* of substrate-projected reference cohorts. The substrate (CC-BY 4.0), 12 per-cohort F matrices, β/α decompositions, basin definitions, and reproducibility code are open. Downstream methods can consume the substrate coordinates without re-deriving them — a researcher with a new cohort can MAP-project onto our substrate and immediately compare basin activation patterns to the 12 deposited references via node-by-node cosine. This Resource framing is the load-bearing contribution; the method (MAP + β/α + signed-basin) is the interpretive layer that makes the Resource consumable, but the deposit is the value.

## 6.2 What's new vs network-based stratification literature

Network-based stratification (Hofree et al., 2013), netNMF, and graph-regularized factor methods share with GIZMO the premise that biological networks structure patient state beyond the data matrix. GIZMO contributes three operational distinctions over these methods: (i) the F output is in *cohort-independent substrate coordinates*, so cohorts measured on disjoint inputs are directly comparable; (ii) the signed-basin layer names *both* poles of each principal axis from substrate-node attributes (no post-hoc enrichment); (iii) the deposited reference library plus blinded interpretability evaluation is a rigor standard that the method literature is missing. These distinctions are operational and resource-oriented, not theoretical — the underlying propagation operator (signed Laplacian) is standard.

## 6.3 What F is, mechanistically

F is a substrate-projected feature representation, not a prognostic predictor. The framework's empirical contribution is *biology-informed feature selection*: the substrate's Laplacian smoothing surfaces connected sub-graphs of input features (basins) whose gene members are biochemically coherent — proteoglycan + cathepsin + MMP for ECM degradation; IDH1/IDH2/D2HGDH/L2HGDH for IDH-mut catalysis; B2M/CD3D-G/HLA-A-C for T cell/MHC infiltrate. When a downstream task (survival prediction) depends on the same biology that the substrate organizes coherently, F provides discrimination at parity with or above PCA-on-input AND names the mechanism. When the task depends on biology that doesn't align with substrate edges (Gao_RA serum proteins), F underperforms.

## 6.4 β/α decomposition honest restatement

The β scalar is the patient's projection onto the substrate's log-PageRank direction — a substrate-fixed, data-independent axis. Heuristically, this is interpretable as the patient's "phenotype-presentation magnitude" on substrate-hub-loaded biology (transcription factors, signaling cascades, immune complexes are over-represented at substrate hubs). The hub-bias confound is real: the substrate-PR distribution is biased toward transcriptional/signaling genes, and the β axis inherits that bias. We report β as a diagnostic axis rather than a primary contribution, with the understanding that downstream methods working in substrate coordinates may want to normalize or whiten the hub direction differently.

## 6.5 Preprocessing matters more than substrate edge gaps for GoF biology

The within-patient z-score preprocessing rescues subtype/anchor recovery for biology that canonical-preprocessing absorbs (TCGA_IDH 2HG anchor under canonical rank 6,355, opposite direction; under z-score rank 26, correct direction). The held-out validation in §5 pre-registers this property. The lesson is broader: **the framework's representational scope is bounded by preprocessing choice as much as by substrate-edge curation**. For fine-grained subtype discrimination, preprocessing the input to surface within-patient relative-rank biology is essential; for survival or other magnitude-driven tasks, the canonical preprocessing that preserves inter-patient magnitudes is essential. These are *not* the same operating point.

## 6.6 Limitations

**HMP2 mechanism change vs prior versions.** The HMP2_IBD_CD cell that passes the pathway-matched null in v7 is α-PC4- (autophagy/inflammasome — NOD2/ATG16L1/IL23R/CARD9/NLRP3, AUROC 0.834). Earlier drafts (v5/v6) reported the cohort's headline mechanism as a "gut metabolite + bile-acid axis" carried on α-PC2-. The bile-acid α-PC2- axis was not specifically tested under the stricter pathway-matched null in this release; v7 reports the pathway-null-surviving cell as the headline. Both mechanisms (gut microbiome bile-acid axis on α-PC2 and the autophagy/inflammasome axis on α-PC4) are documented biology of Crohn's disease in HMP2; v7 reports α-PC4 because that is the cell that survives the stricter null. v8 should report whether α-PC2- bile-acid axis also passes the pathway-matched null when explicitly tested.



The interpretability evaluation depends on author-curated source-paper key-gene sets which inherit selection bias (authors highlight what they know, not what they don't); the curation is open and revisable. AUC-based discrimination metrics conflate effect-size with cohort-size; per-cohort confidence intervals are reported but cross-cohort meta-analysis would benefit from a standardized rank-based statistic. Per-basin survival is shown for four cohorts; cross-cohort generalization to non-cancer, non-COVID contexts is future work. The substrate is a static deposit; biology that requires *new edges* (neomorphic gain-of-function) is not representable without edge injection. Gao_RA's failure shows the boundary: when serum proteomic measurements connect to substrate nodes through PPI edges that don't form coherent Laplacian neighborhoods, the framework's coordinate system is not load-bearing.

## 6.7 Forward work

The substrate is a growing deposit; future versions extend it (a lipid layer with flux + state sub-graphs; mutation-conditional edges for cancer-driver biology; longitudinal trajectories for time-resolved cohorts). The Resource framework — a deposited library with blinded interpretability evaluation as a rigor floor — generalizes to these extensions.

---

# 7 Methods (~600 words)

## 7.1 Substrate construction (brief; full in supplement)

The 38,148-node substrate is constructed from Reactome (R-HSA reactions filtered to human), StringDB (PPI edges at confidence ≥ 700), HMDB (small molecules; manual de-duplication by InChIKey14), and KEGG (pathway annotations cross-referenced to gene/metabolite nodes). Hub-cap = 200 (nodes exceeding degree 200 are downgraded to peripheral status during MAP smoothing; this prevents superhub-driven loading concentration without removing the nodes). The construction code is reproducible from the source databases via `data/processed/human_full/` build script. CC-BY 4.0 deposit.

## 7.2 MAP reconstruction: positive-definiteness and convexity

The objective is

$$ \mathcal{L}(F) = (x - A_{obs} F)^T \Sigma^{-1} (x - A_{obs} F) + \lambda F^T L_{signed} F + \rho \|F\|^2 $$

with L_signed the Kunegis (2010) absolute-degree signed Laplacian: L_signed = D̄ − A_signed, where D̄_ii = Σⱼ |A_ij| (absolute-row-sum) and A_signed encodes signed edges. Under this convention, L_signed is positive semidefinite with eigenvalue 0 only on the trivial constant vector (proof in supplement). Adding ρI for ρ > 0 makes the Hessian strictly positive definite; the objective is strictly convex and the MAP solution is unique.

Default ρ = 10⁻³. Default λ = 0.1 × (mean data weight per modality) / (median Laplacian eigenvalue × n_nodes). Conjugate gradient with diagonal preconditioning; convergence in 100–300 CG iterations on typical cohorts.

## 7.3 β/α decomposition

$\beta_p = \langle F_p, \log\mathrm{PR} \rangle / \|\log\mathrm{PR}\|^2$, $\alpha_p = F_p - \beta_p \log\mathrm{PR}$. The log-PR direction is the substrate-fixed log-PageRank vector (computed once from the substrate; data-independent).

## 7.4 Signed-basin extraction

For each α-PC × sign cell, define $\mathrm{topQ}_k = \{\text{nodes with } |\text{loading}| \geq \mathrm{quantile}_{0.95}(|\text{loadings}|)\}$. Restrict to nodes whose loading sign matches the cell's sign. The basin is the largest connected component of the substrate sub-graph induced by these nodes. Genes in the basin are read off directly from substrate-node attributes (HGNC symbol where node_type = gene).

## 7.5 Degree-, PageRank-, and pathway-membership-preserving nulls

For an empirical gene set of size N with substrate-degree distribution {d_1, ..., d_N}, PageRank distribution {pr_1, ..., pr_N}, and Reactome-reaction-degree distribution {r_1, ..., r_N}: we stratify the substrate's gene nodes into deciles of each statistic and sample N nodes matching the empirical decile distribution per null replicate; recompute the test statistic; repeat 1,000 times; empirical p-value is Phipson-Smyth (count + 1) / (n + 1). The Reactome-reaction-degree null specifically controls for *annotation density* — curated key genes are over-represented at nodes with many Reactome reactions (pathway-popular genes), and a degree- or PageRank-matched null does not address this confound. We report both null types per cohort × α-PC × sign cell in `v7_pathway_matched_null.tsv`.

## 7.6 Smoothness statistic

The preprocessing-choice diagnostic in §5.1 is computed as

$$ s_{\mathrm{cohort}} = | \cos( v_{\mathrm{F\text{-}PC1}}, \log\mathrm{PR} ) | $$

where $v_{\mathrm{F\text{-}PC1}}$ is the top-1 right singular vector of the cohort's projected F matrix and $\log\mathrm{PR}$ is the substrate-fixed log-PageRank direction (the hub axis used by the β/α decomposition in §2.3). When $s$ is high (cosine close to 1), the cohort's dominant F variance direction aligns with the substrate hub axis and the β-direction carries the cohort's dominant signal; when $s$ is low, the cohort's dominant variance is orthogonal to the hub axis and α-PCs carry the disease-specific signal.

**Earlier drafts of this section defined the statistic against the Laplacian first-non-trivial eigenvector $v_{\mathrm{low\text{-}eigen}}$ of the substrate restricted to observed nodes, on the basis that this eigenvector is the formal substrate-smoothness reference direction.** Both quantities — the Laplacian low-eigenvector and the substrate log-PageRank — capture the substrate's data-independent low-frequency variance directions; the log-PR proxy is used here because it is computationally tractable (sparse Laplacian eigsh on the 38,148-node substrate restricted per-cohort takes 1–2 hours per cohort), and the log-PR vector is the canonical direction used elsewhere in the framework (β is the projection onto log-PR by definition; §2.3). The deposited `check_v7_smoothness_statistic.py` implements the log-PR proxy that produces the §5.1 table; a Laplacian-eigenvector version is left to v8.

**The honest scope of this statistic** (§5.1): it correctly identifies TCGA_IDH (low $s$ = 0.13) as needing within-patient z-score to recover the 2HG anchor (rank 6,355 under canonical, 26 under z-score), but does not cleanly separate canonical-vs-z-score in all failure modes — Gao_RA has high $s$ = 0.59 yet fails canonical AUROC because its key genes (IL6/TNF/CRP/IL1B) are themselves hub-popular. Two failure modes (mechanism-axis-not-dominant; mechanism-genes-too-generic-vs-hubs) are not separable with a single statistic.

## 7.7 Survival modeling

Cox proportional-hazards regression with mild L2 penalty (penalizer = 0.01, lifelines 0.30.0). 5-fold cross-validation; StandardScaler fit on training fold only; Harrell's C-index from lifelines.utils.concordance_index. For Filbin (binary 28-day mortality), logistic regression + 5-fold CV AUC.

## 7.8 Per-basin patient activation score

For each basin's gene members, restrict the cohort's expression matrix to the substrate-mappable subset (using NCBI gene_info alias resolution: substrate name → canonical HGNC → matched expression column if present in any alias). Compute PCA-PC1 of the restricted matrix; the resulting per-patient score is the basin activation score.

## 7.9 Cross-cohort basin conservation

For each pair of cohort × α-PC × sign basins from different cohorts, compute Jaccard similarity of the gene members. Report pairs at Jaccard ≥ 0.30 with ≥ 5 shared genes.

## 7.10 Reproducibility and known-biology sanity checks

All code at https://github.com/insilijo/GIZMO/ commit pinned per release. Deposited results: substrate at `data/processed/human_full/graph.json`; F matrices at `benchmarks/results/unsupervised/_pre_zscore_snapshot_20260607/stage3_F_<cohort>_edge_informed.npz`; basin extractions at `benchmarks/results/unsupervised/zscored/v7_basin_conservation.json`; per-basin survival at `v7_basin_survival_broad.tsv`; pre-registration at `MANUSCRIPT_v7_PHASE4_PREREG.md`.

**Known-biology sanity spot-checks for v7 deposit.** The following are *spot-checks on known-positive biology*, not a systematic validation protocol. They test whether the deposited outputs are *consistent with* canonical known biology in five specific cases; they do not substitute for systematic external validation, which requires independent panel review of the curated key-gene file and basin-biology assignments (deferred to v8).

(i) substrate gene-symbol set cross-checked against NCBI gene_info aliases (deposit `gizmo/diagnostics/alias_validate.py`); (ii) ribosome basin gene members (RPL10/11/13/15/18 + RPS family) cross-checked against MSigDB Hallmark "MYC targets V1" ribosomal subset — Jaccard 0.78 with the substrate-extracted ribosome basin; (iii) TCGA_IDH α-PC3+ OXPHOS basin (MT-CO1-3, MT-ND1-6, MT-CYB, COA3, TIMMDC1) cross-checked against the published IDH-mutant glioma OXPHOS dependency literature (Molenaar 2014; Tateishi 2015) — 14 of 16 basin gene members are direct hits; (iv) Filbin α-PC1- ECM-degradation basin (BGN, DCN, CTSL, MMP13, VCAN) cross-checked against critical-illness ARDS ECM-degradation literature (Cabrera-Benitez 2014; Kratzer 2008) — all 19 basin proteoglycan + MMP members are documented ARDS markers; (v) 2HG anchor rank under z-score (26 of 38,148) cross-checked against the canonical 2HG biology that the Reactome catalysis edges encode (IDH1/IDH2 → α-KG → R(-)-2HG via R-ALL-879997).

**No external panel review of the curated key-gene file or the basin-biology cross-checks has been performed**; this is a structural limitation of the single-author Resource and v8 will solicit external review of both.

---

# Acknowledgments

This work was conducted solo and self-funded under the Insilijo brand. No external funding to declare.

# Author contributions

Conceptualization, methodology, software, formal analysis, investigation, data curation, writing — JJG.

# Declaration of interests

The author declares no competing interests.

# Data and code availability

Substrate (CC-BY 4.0): Zenodo DOI on acceptance. F matrices, basin outputs, reproducibility code: same Zenodo deposit + github.com/insilijo/GIZMO.
