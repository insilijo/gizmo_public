# v6 session notes (2026-06-09)

v5 narrative preserved as canonical; v6 adds a targeted §5 retraction of the LoF/GoF representational-scope claim, with within-patient z-score introduced as a **fine-grained subtype/anchor diagnostic preprocessing** (Methods §"Within-patient z-score") — *not* as a canonical replacement. Manuscript shipped as docx with native OMML equations.

## v6 reframe

After running the within-patient z-score preprocessing on 13 of 17 panel cohorts and sweeping discrimination AUC + v5↔v6 PC alignment, the empirical pattern that emerged:

- **At active-vs-control discrimination, both preprocessings are comparable** (panel-wide AUC sweep, `benchmarks/check_v6_auc_sweep.py`). Mean best-PC AUC: heterogeneous 0.773, strong-driver 0.739. Neither dominates.
- **Within-patient z-score specifically rescues fine-grained subtype + anchor recovery** within homogeneous-driver diseases:
  - TCGA_IDH 2HG-mito anchor rank 6,355 → **26** ($d$ sign reverses)
  - TCGA_LUAD KRAS-vs-EGFR subtype multi-PC LR LOOCV 0.509 → **0.703** (best 1-PC 0.634 → 0.706)
  - IDH_glioma Trautwein IDH-mut-vs-wt best α-PC 0.65 → **0.842**
- **v5↔v6 within-cohort PC alignment is low** (\|cos\| 0.24–0.44, Jaccard@50 0.05–0.22 across 6 testable cohorts) — basin biology restructures at the node level even when discrimination is preserved.

Conclusion: z-score is a *diagnostic preprocessing* for fine-grained tasks where v5's per-modality global-std preprocessing absorbed within-patient relative-rank signal. It is *not* the canonical preprocessing, and the v5 panel-wide story (basin examples, MOFA+ comparison, horizontal meta-analysis, figures) survives unchanged.

## Manuscript v6 changes from v5

- **Abstract scope paragraph** — rewritten to retract LoF/GoF representational-scope claim with the three verified rescue numbers as evidence.
- **§5** — LoF/GoF falsification narrative replaced with z-score rescue narrative (preprocessing-artifact interpretation), with v5-canonical numbers preserved for the alternative-preprocessing comparison.
- **Discussion §"LoF/GoF scope boundary…"** retitled to *"Preprocessing matters more than substrate-edge gaps for GoF biology (v6 reframe)"* and rewritten.
- **Discussion paper-series paragraph** — Paper 4 scope narrowed to *"truly neomorphic biology with zero canonical effector signaling AND not rescuable by within-patient z-score"*.
- **Methods §"Within-patient z-score"** added — introduces the preprocessing with formula and reference implementation pointer.
- **All v5 plaintext equations** converted to LaTeX so pandoc renders them as native Word OMML equations. v6 docx ships with **6 display + 42 inline OMML equations**.

§1, §2 (basin examples), §3 (horizontal meta-analysis with 3 cohorts), §4 (MOFA+ comparison) are preserved verbatim from v5.

## Z-scored F matrices generated (13 of 17 cohorts)

For users who want to reproduce the §5 z-score rescues, F matrices under within-patient z-score preprocessing are at `benchmarks/results/unsupervised/stage3_F_<cohort>_zscored.npz`; downstream artifacts at `benchmarks/results/unsupervised/zscored/<cohort>/`.

Cohorts with z-scored F: Crohn, Su_COVID, Gao_RA, IDH_glioma, TCGA_IDH_glioma, Filbin_COVID, Erawijantari, KMPLOT_BRCA, TCGA_LUAD, GSE89408_RA, HMP2_IBD_CD, CorEvitas_RA, Wang_RA.

Cohorts deferred (4 of 17, not blocking v6 — manuscript uses v5 canonical numbers for these):
- GSE65391_SLE (n=996), GSE65682_sepsis (n≈800) — WSL2 7.6 GB memory bound on current load; need streaming-loader patch.
- CPTAC_CCRCC / COAD / OV / GBM — 3-modality custom runner at `data/cohorts/CPTAC_*/run_gizmo_full.py` needs z-score preprocessing patch.
- breast_TCGA_DIABLO — separate runner.

## Load-bearing scripts cited in v6 manuscript

| Script | Cited at | Purpose |
|---|---|---|
| `benchmarks/batch_rerun_zscored.py` | Methods §"Within-patient z-score" | Per-cohort z-score MAP solve + β/α + signed-basin pipeline |
| `benchmarks/check_2hg_rank_fast.py` | §5 | TCGA_IDH 2HG anchor recovery test |
| `benchmarks/check_luad_multipc_loocv.py` | §5 | LUAD KRAS-vs-EGFR multi-PC LR LOOCV test |
| `benchmarks/check_v6_auc_sweep.py` | Methods §"Within-patient z-score" | 13-cohort panel-wide AUC comparison |
| `benchmarks/idh_ablation_map_zscored_heat.py` | reference impl | Original within-patient z-score reference function |
| `benchmarks/luad_ablation_map_zscored_input.py` | §5 | LUAD ablation showing AUC 0.835 under z-score + diffusion_t=0 |
| `gizmo/inference/projection.py` — `_within_patient_zscore` | Methods §"Within-patient z-score" | Canonical implementation |

## Iterative / exploratory scripts (archived)

Moved to `benchmarks/_iterative/2026_06_09_session/`:
- `check_2hg_rank_zscored.py` — slow version superseded by `check_2hg_rank_fast.py`
- `check_section3_anchors.py` — exploratory; v5 anchor node IDs don't resolve in substrate as deposited (informed future re-resolution work but not cited)
- `check_v5_v6_pc_alignment.py` — exploratory diagnostic; produced cosine/Jaccard alignment numbers that informed the v6 reframe decision but isn't a manuscript citation
- `interrogate_glyoxalase_axis.py` — exploratory; confirmed glyoxalase axis on Wang_RA/Filbin_COVID is real biology (not hub attractor), informed discussion but not cited

## Manuscript drafts on disk

- `MANUSCRIPT_v5_DRAFT.md` — prior version, unchanged
- `MANUSCRIPT_v6_DRAFT.md` / `.docx` — v6 canonical (v5 + targeted §5 + Methods + abstract patches)
- `MANUSCRIPT_v6_HEAVY.md` / `.docx` — abandoned full-rewrite draft (reverted because empirical evidence didn't support reshaping the whole paper around z-score)
- `MANUSCRIPT_v5_EQUATIONS.md` / `.docx` — earlier standalone equations file from §"can you go through the manuscript and write the equations in docx" (June session)

## Backups

- Pre-z-score F matrix snapshot: `benchmarks/results/unsupervised/_pre_zscore_snapshot_20260607/stage3_F_*.npz` (49 files, 1.1 GB)
- Heavy v6 draft (preserved): `MANUSCRIPT_v6_HEAVY.{md,docx}`

## Verified results files

- `benchmarks/results/unsupervised/zscored/twohg_anchor_recovery.json` — TCGA_IDH 2HG rank 26 / d=+0.630
- `benchmarks/results/unsupervised/zscored/luad_multipc_loocv.json` — LUAD multi-PC AUC 0.7026, best 1-PC 0.7055
- `benchmarks/results/unsupervised/zscored/v6_auc_sweep.json` — 13-cohort panel-wide AUC sweep
- `benchmarks/results/unsupervised/zscored/v5_v6_pc_alignment.json` — 6-cohort v5↔v6 cosine + Jaccard alignment
- `benchmarks/results/unsupervised/zscored/glyoxalase_interrogation.json` — glyoxalase axis diagnostic
- `benchmarks/results/unsupervised/zscored/batch_summary.json` — batch 1 summary
- `benchmarks/results/unsupervised/luad_ablation_map_zscored.json` — prior LUAD ablation (AUC 0.835)
