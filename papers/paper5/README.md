# Paper 5 — Longitudinal patient-state tracking and cross-cohort substrate-axis projection

**Version: v0.1** (scaffolded; outline-stage; substantial pre-existing analyses staged but manuscript prose not yet written)

This folder accumulates the work for a paper covering (i) longitudinal patient-state tracking under intervention (Filbin acute COVID D0→D3, TB antibiotic cure, HCV DAA cure, MTX RA, T2D RYGB bariatric — five paired cohorts), (ii) cross-cohort substrate-axis projection methodology (baseline-α PC basis defined in cohort A, projected onto cohort B's trajectory), and (iii) future architectural extension to spatial-tile projection.

Substrate (`../../substrate/graph.json`) and framework code (`../../gizmo/`) are shared. β/α decomposition is from [Paper 1](../paper1/MANUSCRIPT.md).

## Status snapshot

| Element | Status | Source |
|---|---|---|
| **OUTLINE.md** | Drafted (2026-06-01); abstract + section structure for cross-cohort axis projection sub-paper | `OUTLINE.md` |
| **LONGITUDINAL_SUMMARY.md** | Drafted (2026-05-26); six findings across 5 cohorts | `LONGITUDINAL_SUMMARY.md` |
| Longitudinal calibration on 5 cohorts (n=363 paired observations) | **Run + audited (LOOCV + DeLong CI + permutation null)** | `scripts/loocv_longitudinal_audit*.py`, `results/loocv_longitudinal_audit*.json` |
| Filbin COVID D0→D3 (n=214) — α-PC3 sign-flip + α-PC1 magnitude tracking | **Computed** | `scripts/filbin_d0_d3_calibration.py`, `results/filbin_d0_d3_calibration.json` |
| TB cure calibration GSE89403 (n=70) — β 0.826→0.686 (57% normalization) | **Computed** | `scripts/gse89403_tb_cure_calibration.py`, `results/gse89403_tb_cure_calibration.json` |
| HCV DAA cure (n=4) — β 0.833→0.625 (62%) | **Computed** | embedded in cure_bundle scripts/results |
| Cross-cohort axis projection (Filbin-PC5 TNF → TB cure outcome) | **Stress-validated** (jackknife p<.014, random-axis null 0/300, β-orthogonal, baseline-balanced) | `scripts/paper4_validation_battery.py`, `paper4_loo_stability.py`, `paper4_followups.py` |
| TB progressor cohorts (GSE94438 Zak GC6-74, GSE107994 Leicester) | **Downloaded + verified** (downloaded 2026-06-02); replication analysis run | `scripts/tb_progressor_replication.py`, `results/tb_progressor_replication.json` |
| Reachability × elasticity 4-category typology | **Drafted** in outline | (see OUTLINE.md §3.4) |
| β persistence vs α relaxation rate (~2× difference) | **Quantified** | `results/cure_bundle_*.json` |
| Manuscript prose (Abstract, Intro, Results §3, Discussion, Methods) | **NOT YET WRITTEN** | — |
| Spatial-tile projection extension | **Scoped only**; no work yet | flagged in OUTLINE.md Discussion |
| Protein-state substrate extension | **Scoped** (~1 week kinase phospho-substrate variant; ~3 weeks full Reactome protein-state ingestion) | flagged in `project_gizmo_paper5_longitudinal_spatial.md` memory |

## Contents

```
papers/paper5/
├── VERSION              # v0.1
├── README.md            # this file
├── OUTLINE.md           # 2026-06-01 outline: thesis, abstract, section structure
├── LONGITUDINAL_SUMMARY.md  # 2026-05-26 findings summary across 5 cohorts
├── figures/             # 5 existing longitudinal figures
├── scripts/             # 21 analysis scripts (longitudinal calibration + cross-cohort projection + paper4 battery)
├── results/             # 10 result deposits (JSONs + 1 CSV)
├── data/                # (placeholder for cohort registry, curated axes, etc.)
└── docs/                # (placeholder for design docs, decision logs)
```

## Headline findings to-date

1. **α-PCs are clinically-decomposable longitudinal axes — emerging without supervision.** Filbin COVID α-PC3 sign-flips Improved (−0.120) vs Worsened (+0.239); α-PC1 magnitude is 16× larger in Worsened. Wang RA α-PC1 is 2.2× larger in EULAR Responders; α-PC2 8× larger.

2. **β achieves substantial but incomplete normalization under curative interventions.** TB antibiotic cure β 0.826 → 0.686 (57% relaxation, residual 43% disease metabolic set-point persists); HCV DAA cure 0.833 → 0.625 (62% relaxation). Month-scale α relaxes ~2× faster than β.

3. **Cross-cohort axis projection correctly partitions drug mechanism of action.** A TNF/integrin axis derived from acute COVID (Filbin α-PC5) is NOT engaged by methotrexate (Wang RA) or DAA (HCV), but IS strongly engaged by TB antibiotics (GSE89403, p_grp=0.003). Conjunction p=0.014 (2000-axis empirical null). Mechanism aligns with TNF/granuloma TB immunology — including correct directionality (higher baseline TNF activity predicts treatment failure, matching TB pharmacology rather than naive immunodeficiency intuition).

4. **Multi-axis coordinated cure response.** TNF + IFN-γ + mitochondrial NDUF/TIMM form three orthogonal axes with bootstrap-validated pairwise ρ 0.54–0.75 — recapitulating granuloma-resolution biology.

5. **Reachability × elasticity 4-category typology** (drafted): Category 3 inelastic axes = disease-state monitors; Category 2 semi-elastic = stratification biomarkers; pharmacodynamic-elastic axes = drug-on confirmation. β = prototype inelastic disease-endemic axis.

## What's needed to reach v1.0 (manuscript-ready)

- [ ] Write Abstract, Introduction, Results §3 (cross-cohort projection), Discussion, Methods §2 prose
- [ ] Build figure set (currently 5 longitudinal figures; need Figure 1 conceptual + chord/grid for axis projection)
- [ ] Resolve scope: this paper is currently "cross-cohort axis projection" with longitudinal as substrate. Either keep that scope OR widen to include the full longitudinal calibration as a co-headline. Decision pending.
- [ ] External replication: TB progressor cohorts (GSE94438, GSE107994) — replication analysis already run; integrate findings into Results §3 or supplement.
- [ ] Spatial-tile section: either (a) drop from this paper and defer to Paper 6, or (b) include as architectural-enablement Discussion paragraph with no empirical claim.

## Reproducing the staged analyses

```bash
cd /path/to/gizmo_public

# Filbin D0→D3 longitudinal calibration
PYTHONPATH=. python3 papers/paper5/scripts/filbin_d0_d3_calibration.py

# TB cure calibration
PYTHONPATH=. python3 papers/paper5/scripts/gse89403_tb_cure_calibration.py

# Cross-cohort axis projection (Filbin → TB) — the stress-validated proof-of-concept
PYTHONPATH=. python3 papers/paper5/scripts/paper4_validation_battery.py
PYTHONPATH=. python3 papers/paper5/scripts/paper4_loo_stability.py
PYTHONPATH=. python3 papers/paper5/scripts/paper4_followups.py

# TB progressor replication
PYTHONPATH=. python3 papers/paper5/scripts/tb_progressor_replication.py
```

Most scripts depend on cohort F matrices from the Paper 1 Zenodo deposit (DOI TBD on Paper 1 acceptance); see Paper 1 [`papers/paper1/README.md`](../paper1/README.md) for substrate + F-matrix paths.

## Related memory

- `project_gizmo_paper5_longitudinal_spatial.md` — original scope (longitudinal + spatial + protein-state extension)
- `project_gizmo_alpha_pc_longitudinal.md` — Filbin D0→D3 findings
- `project_gizmo_beta_persistence_under_cure.md` — β-persistence audit + normalization rates
- `project_gizmo_cross_cohort_axis_projection.md` — stress-validation framework
- `project_gizmo_axis_elasticity_typology.md` — 4-category reachability × elasticity grid
- `project_gizmo_tb_progressor_cohorts.md` — GSE94438 + GSE107994 verified cohorts
