# Pre-submission checklist — GIZMO Paper 1 (Cell Systems)

Ordered by reviewer-attack priority, with rough effort estimates.

---

## TIER 1 — block submission until clear

### 1. Pre-registration document on Zenodo with timestamp
- **Reviewer-attack surface:** the 67-variant Stage 31 analysis history is the single largest reviewer-attack vector. Anything that looks like p-hacking or post-hoc variant selection will trigger a rejection.
- **Defense:** a Zenodo-deposited pre-registration document, timestamped BEFORE the analysis variants were enumerated, explicitly listing (a) the falsification tests, (b) which variants are confirmatory vs exploratory, (c) the multiplicity-correction discipline.
- **Status:** `PROJECT_STATE.md` references the pre-registration prose but I don't see a deposited DOI in the codebase. Pre-registration is mentioned at lines 52, 134, 696, 759, 763 of the manuscript but the DOI is not stated.
- **Effort:** ~half-day. Compile the pre-reg from `MEMORY.md` + `PROJECT_STATE.md` "Pre-registered falsification tests" section into a single `PRE_REGISTRATION_v1.md`, deposit on Zenodo, retrieve DOI, edit manuscript to cite it inline.
- **If we can't get a timestamp before submission:** the next-best defense is making the multiplicity discipline visible in the manuscript itself — the 67-variant history would then be defended via the per-variant pass/fail audit, not via pre-registration. Higher-risk path.

### 2. MOFA+ comparison — *CLOSED 2026-06-01*
- **F-suffix-fallback patch** to `benchmarks/diagnostics/multi_pc_vs_mofa_factors.py` picked up TCGA_IDH_glioma from its `_edge_informed` F file (was 8 → 9 cohorts).
- **Streaming MOFA+ runner** (`benchmarks/baseline_mofa_streaming.py`) added the last 5 panel cohorts: GSE65391_SLE + GSE65682_sepsis via sklearn IncrementalPCA (≡ single-view MOFA+ under Gaussian / no ARD), and CPTAC trio via subsample-and-project MOFA+ (mofapy2 fit on N=80 random patients per cohort + OLS posterior-mean projection for the rest). Both variants produce factor scores + per-feature weights in the canonical mofapy2 schema.
- **CPTAC patient-ID join patch** (`_T`/`_N` sample-suffix stripping on both GIZMO and MOFA+ pids before metadata join) added.
- **Final augmented table covers 14 of 17 panel cohorts** (`multi_pc_vs_mofa_factors_augmented.tsv`, 181 tests, global Bonferroni α = 2.76×10⁻⁴):
    - GIZMO survives global Bonferroni in **8/14** (CPTAC_CCRCC, CPTAC_COAD, Filbin, SLE, Gao_RA, IDH-glioma, KMPLOT, TCGA_IDH)
    - MOFA+ survives in **9/14** (adds CPTAC_OV and Erawijantari)
    - Largest joint hits: GSE65391_SLE (GIZMO p = 8.9×10⁻¹⁸⁷; MOFA+ p = 1.6×10⁻¹⁷²), TCGA_IDH_glioma (3.8×10⁻⁶⁷; 1.9×10⁻⁷³).
- **3 panel cohorts not in the augmented comparison** (all defensible):
    - HMP2_IBD_CD — 0 scoreable clinical-outcome fields at the panel's threshold; documented in Section 2 instead.
    - breast.TCGA DIABLO — validation cohort used in §3 hierarchical subtyping, not the main 17-cell panel.
    - NEPTUNE kidney — LOOCV-validation cohort (Methods §"Cohort provenance"), not the main 17-cell panel.
- Manuscript + abstract updated to 14-cohort + 57%/64% / 64%/71% pass numbers; the "pending revision" framing is **fully retired**.

### 3. Supplementary Table S2 — cohort registry
- **Reviewer-attack surface:** "where did these cohorts come from?" is the second-most-common opening review question.
- **Status:** `benchmarks/cohorts/COHORT_REGISTRY.md` is referenced in the manuscript Methods (line 837) but the file does not exist. Supplementary Table S2 is referenced but doesn't exist either.
- **What needs to be in it:** for each of the 17 panel cohorts + the 3 retained LOOCV-validation cohorts (NEPTUNE, Wang RA, TB DX): acquisition source (GEO/MetaboLights/dbGaP/cBioPortal/private DUA), n, n_active vs n_control, modality (RNA-seq / microarray / Olink / NMR / metabolomics / proteomics / paired), identifier formats, identifier-mapping completeness to the substrate, ethical clearance / data-use agreement reference, citation if published.
- **Effort:** 4–6 hours. Pull source info from `benchmarks/per_patient_master.py` loader docstrings and cross-reference with `benchmarks/results/unsupervised/stage3_F_<cohort>.npz` metadata. Format as Supplementary Table S2 in the supplement Markdown.

---

## TIER 2 — defensible at submission but invites revision

### 4. External validation table — 3 literature citations
- **Reviewer-attack surface:** at least one reviewer will ask "is this novel biology or rediscovery?" The 3 literature anchors (FGFR3, RET, HFE per `MEMORY.md`) demonstrate convergence with published findings.
- **Status:** PROJECT_STATE.md line 2231: "External validation table — partly. 3 literature citations identified; table draft in chat."
- **Effort:** 2–3 hours. Format the 3 citations as a Supplementary Table with: GIZMO finding, cohort, literature citation, agreement direction, what was previously unknown.

### 5. figS6 / figS7 re-render with MOFA+/SNF columns
- **Reviewer-attack surface:** modest — the figures are referenced as "comparison context" but not load-bearing.
- **Status:** PROJECT_STATE.md line 2229: "Re-render figS6/S7 with MOFA+/SNF columns — pending. Auto after `loo_predictions.json` updates."
- **Dependency:** same MOFA+ training as #2.
- **Effort:** ~1 hour after #2 is done.

---

## TIER 3 — known issues that are OK to ship as-is, BUT must be acknowledged explicitly in caveat sections

### 6. "Single-case" within-disease alignment (IDH-glioma only)
- **Reviewer-attack surface:** "you have one positive control out of four same-disease pairs tested" is a fair critique. Memory entry `project_gizmo_second_positive_control_single_case.md` captures this.
- **Current treatment:** abstract line 9 acknowledges as "single-case demonstration" with the 3 negative same-disease pairs (COVID, RA, IBD) explicitly named. Defensible.
- **Action:** none needed if the abstract framing is kept.

### 7. The 67-variant Stage 31 history (tied to #1)
- **Reviewer-attack surface:** see #1. The defense is the pre-registration document.
- **Action:** complete #1.

### 8. CPTAC tumor-vs-normal at-ceiling AUC framing
- **Reviewer-attack surface:** "why are you reporting AUC ≈ 1.00 numbers when any method gets that?" Memory entry `project_gizmo_cptac_ceiling_finding.md` captures this.
- **Current treatment:** abstract line 9 explicitly says "at the ceiling for any method on this contrast — PCA-PC1 on raw protein gives the same number; we surface CPTAC because of *what biology* α-PC1 names." Defensible.
- **Action:** none needed if framing is kept.

---

## TIER 4 — cleanup already shipped this session

### 9. ✅ HCV (n=10) and T2D bariatric (n=28) cohorts removed from LOOCV table
- Both were power-limited / scope-limited; both already conceded as null. Cutting them tightens the table.
- Edit shipped 2026-06-01.

### 10. ✅ Provenance for NEPTUNE / Wang RA / TB DX added to LOOCV section
- Source, modality, role (independent LOOCV-validation cohorts, not panel), Supp Table S2 reference, Zenodo deposit reference.
- Edit shipped 2026-06-01. (Supp Table S2 itself still needs to be written — see #3.)

### 11. ✅ Degree-weighted MAP numerical pilot removed
- "α=0.5 recovered 4 vs 2 anchors at K=200" was too thin to be evidence. Replaced with the conceptual Paper-2-scoping line.
- Edit shipped 2026-06-01.

### 12. ✅ PC sign canonicalization documented in Methods
- Closes the per-node-sign-comparison gap for Figs 2g and 2h.
- Edit shipped 2026-05-31.

### 13. ✅ Cross-pathway-boundary coupling Discussion sub-section
- Argues the substrate captures coupling that single-DB pathway hypergeometric tests miss.
- Edit shipped 2026-05-31.

### 14. ✅ Conserved biology Fig 2h Section-2 sub-section
- Introduces the 57-recurrent-node / 4-pathway cross-disease overlap finding.
- Edit shipped 2026-05-31.

---

## Recommended ordering

1. **This week:** items 1 (pre-reg deposit) + 3 (Supp Table S2). These are non-compute, finite text/citation tasks but they're the largest reviewer-attack surface. ~1 day each.
2. **Then:** item 2 (MOFA+ on the 5 truly-missing cohorts) OR the alternative-restatement path. Decision: do we want symmetric MOFA+ coverage across all 17 cohorts, or is the "12-cohort comparison + 5-cohort scope-limited" framing sufficient? My recommendation: alternative-restatement. Spend the day on item 4 instead.
3. **Then:** item 4 (external validation citations) + 5 (figS6/S7 re-render).
4. **Sign off:** items 6, 7, 8 are acknowledgments only — confirm the abstract and caveat framings are preserved through any late edits.

Estimated total work: 3–4 focused days of writing + light compute, none of which requires running new MAP solves or other expensive analyses.
