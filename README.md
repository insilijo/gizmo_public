# GIZMO — public repository for Paper 1

A biochemistry substrate as a fixed coordinate system for multi-omic per-patient projection.

This repository contains the public-facing artifacts of GIZMO Paper 1: the substrate, the framework code, the manuscript, the figures, and the per-cohort case-study deposit machinery. The 38,148-node biochemistry substrate (`substrate/graph.json`) is the primary contribution and is licensed CC-BY 4.0.

**Current revision: v7 (Accept after five rounds of hostile peer review).**

---

## Contents

- **[MANUSCRIPT.md](MANUSCRIPT.md)** — the Paper 1 manuscript (v7, current; 8,955 words, 11 embedded figures).
- **[MANUSCRIPT_v7.docx](MANUSCRIPT_v7.docx)** — Word version of the v7 manuscript with embedded figures.
- **[MANUSCRIPT_v7_SUPPLEMENT.md](MANUSCRIPT_v7_SUPPLEMENT.md)** / **.docx** — Extended methods: PSD proof for L_signed, per-modality Σ specification, full substrate construction details, λ-sensitivity table, and 13 reproducibility-table pointers.
- **[MANUSCRIPT_v7_FIGURES.md](MANUSCRIPT_v7_FIGURES.md)** / **.docx** — Figure captions + per-figure reproducibility table.
- **[MANUSCRIPT_v6_DRAFT.md](MANUSCRIPT_v6_DRAFT.md)** / **.docx** — Prior v6 draft, retained for revision history.
- **[PRE_SUBMISSION_CHECKLIST.md](PRE_SUBMISSION_CHECKLIST.md)** — work remaining before submission.
- **`substrate/graph.json`** — the 38,148-node GIZMO biochemistry substrate (Reactome + StringDB + HMDB + KEGG, CC-BY 4.0). Load via `gizmo.export.json_export.read_json`. This is the primary deposit of the paper; downstream methods consume it.
- **`figures/v7/`** — v7 main + supplementary figures (11 PNGs covering substrate schematic, interpretability AUROC heatmap, activation vocabulary matrix, per-basin survival, F-vs-PCA orthogonal panels, per-α-PC ablation, preprocessing scope, OXPHOS basin sub-graph, Filbin ECM basin sub-graph, cross-cohort chord diagram, patient F-space embedding).
- **`figures/`** — legacy v6 manuscript figures + `FIGURE_MANIFEST.md`.
- **`gizmo/`** — Python package implementing MAP reconstruction, β/α decomposition, signed-basin extraction, and the supporting graph/scoring/integration machinery.
- **`benchmarks/v7_paper/`** — the 19 scripts that produce every v7 manuscript number, indexed in `benchmarks/v7_paper/README.md`.
- **`benchmarks/unsupervised_stratification/legacy/stage31_v14{b,c}_wgcna_*.py`** — WGCNA comparison scripts cited in §4.1 (v14b uncorrected; v14c same-universe correction).
- **`benchmarks/`** — additional analysis code (MOFA+ baselines, per-patient master loader, diagnostics, figure builders).
- **`results/v7/`** — v7 manuscript-source data tables: interpretability eval (degree + pathway-matched nulls), basin conservation, per-basin survival, F-vs-PCA orthogonal-panel data, per-α-PC ablation, smoothness statistic, TopPR comparison, WGCNA v14b coherence.
- **`results/cohorts/`**, **`results/zscored_v6/`** — legacy v6 deposits.
- **`data/curation/v7_cohort_key_genes.tsv`** — the 75-entry curated source-paper key-gene set (frozen blinded ground truth for the interpretability evaluation).

---

## What's new in v7 (vs v6)

v7 is a Resource paper recast that addresses five rounds of hostile peer review. Major changes:

1. **Quantitative interpretability evaluation** (§3) — 110-cell AUROC test of curated source-paper key genes vs F-α-PC loadings, under degree-PageRank-preserving AND Reactome-pathway-membership-preserving nulls. **10/11 cohorts pass at p<0.05 under the degree null; 7/11 survive the stricter pathway-matched null** (Crohn thiopurine 0.911, HMP2 autophagy 0.834, IDH-mut catalysis 0.839, Filbin/Su COVID, TCGA_LUAD, TCGA_IDH).

2. **Per-basin patient activation scores** (§4.3) — each named biochemical neighborhood produces a per-patient activation score with mechanism-named prognostic value. TCGA_LUAD T cell/MHC basin C-index 0.608; Filbin α-PC1- proteoglycan-degradation basin AUC 0.776 (vs PCA-on-Olink 0.717, named-mechanism case study). Ribosome basin C 0.480 = chance, correctly identifying housekeeping.

3. **F-basins and PCA components read orthogonal prognostic biology** (§4.4) — F's best prognostic basin shares 0–2/68 genes with PCA's best prognostic component across three cohorts; comparable discrimination, different named axes.

4. **Cross-cohort basin conservation** (§4.5) — 149 basins across 11 cohorts; **24 conserved cross-disease pairs at Jaccard ≥ 0.30, 100% cross-disease after excluding interpretively-unvalidated Gao_RA basins**.

5. **Three v6 claims overturned by re-running their cited ablations** (§4.1):
   - TopPR doesn't "dominate at every K" on v7 curated truth — F-α-PC wins at K=200/500
   - WGCNA per-module pathway coherence is at parity with GIZMO under matched gene universe (v14c median ratio 1.02), not GIZMO-favored
   - "GIZMO modules span 8–40× more Reactome pathways per module" was a granularity-vs-coherence misframing

6. **Quantitative β-vs-severity panel** (§2.3) — β captures < 2% of F variance but tracks clinical severity at p < 0.05 in all 5 cohorts tested (Filbin Acuity ρ −0.254, IDH grade ρ +0.285, BRCA grade ρ −0.208, LUAD stage ρ −0.109; IDH-mut ρ −0.233). Two-regime interpretation flagged as exploratory.

7. **Honest scope statements** — KMPLOT post-hoc coverage filter, self-curation structural limitation, edge-injection biology-capture constraint, held-out validation deferred (pre-registered + timestamped only).

---

## Quick start

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Load the substrate
python3 -c "
from gizmo.export.json_export import read_json
mg = read_json('substrate/graph.json')
print(f'{mg.graph.number_of_nodes()} nodes, {mg.graph.number_of_edges()} edges')
"

# 3. End-to-end pipeline on one cohort
python3 run_paper1.py Crohn

# 4. Replay β/α + basins from a Zenodo-deposited F matrix
python3 run_paper1.py --from-F /path/to/zenodo/F.npz --out results/Crohn_replay

# 5. Reproduce a v7 manuscript number (see benchmarks/v7_paper/README.md)
PYTHONPATH=. python3 benchmarks/v7_paper/check_v7_interpretability_eval.py
PYTHONPATH=. python3 benchmarks/v7_paper/check_v7_pathway_matched_null.py
PYTHONPATH=. python3 benchmarks/v7_paper/check_v7_basin_conservation.py
# ... full list in benchmarks/v7_paper/README.md
```

Full install instructions: [`docs/INSTALL.md`](docs/INSTALL.md). Step-by-step walkthrough: [`docs/USAGE.md`](docs/USAGE.md). All 19 v7 manuscript scripts: [`benchmarks/v7_paper/README.md`](benchmarks/v7_paper/README.md).

---

## What the substrate is

A fixed 38,148-node biochemistry graph (12,872 gene nodes + 11,244 metabolite nodes + 14,032 reaction nodes), merged from four curated sources:

- **Reactome** (release 87, June 2024) — 14,032 reactions
- **StringDB** (v12.0) — protein-protein interactions, confidence ≥ 700
- **HMDB** (v5.0) — metabolite identifiers; de-duplicated by InChIKey14
- **KEGG** (2024-Q4) — pathway annotations cross-referenced

Any patient's multi-omic measurements MAP-project onto this substrate, producing a per-patient state vector **F** in 38,148-dim substrate coordinates that don't depend on cohort or modality. The substrate is the coordinate system; downstream method choices (GIZMO's MAP + β/α decomposition + signed basins, or alternatives like static PageRank / MOFA+ / WGCNA / PCA) operate within it.

See **[MANUSCRIPT.md](MANUSCRIPT.md)** §2 for full substrate construction, **[MANUSCRIPT_v7_SUPPLEMENT.md](MANUSCRIPT_v7_SUPPLEMENT.md)** §S1.1 for full curation details, and `substrate/README.md` for a focused substrate-only intro.

---

## Per-cohort case-study deposit (Zenodo)

Every cohort GIZMO has been applied to is deposited at Zenodo with the full analytical bundle: F matrix, β/α decomposition, signed-basin output, per-cohort metadata, and an auto-populated case-study writeup. The bundle is regenerated by `benchmarks/make_zenodo_bundle.py`. See **[PRE_SUBMISSION_CHECKLIST.md](PRE_SUBMISSION_CHECKLIST.md)** for upload status.

Zenodo DOI: *TBD on acceptance.*

---

## Paper series

This repository is the Paper 1 deposit. Subsequent papers extend the substrate and exploit it in further ways:

- **Paper 2** — lipid layer dual-architecture (flux + state sub-graphs, compositional Bayes for ionization-bias correction), degree-weighted MAP.
- **Paper 3** — disease-informed topology rewiring (network deformation under disease state rather than fixed substrate).
- **Paper 4** — mutation-conditional GoF edge injection (the explicit fix for the structural constraint identified in v7 §4.6).
- **Paper 5** — longitudinal + spatial-tile projection (architectural enablement from substrate-fixed coordinates).

---

## License

- **Substrate + analytical outputs**: CC-BY 4.0.
- **Code**: Apache 2.0.
- **Upstream source data** (per-cohort raw inputs): governed by each source's original license; not redistributed here. See **[LICENSE](LICENSE)**.

If you use the substrate or code, please cite:

> Gardner JJ (2026). *GIZMO: a biochemistry substrate as a fixed coordinate system for per-patient multi-omic state, with a deposited reference library across 12 cohorts.* [submitted]; substrate v1.0 at https://doi.org/[Zenodo DOI TBD].
