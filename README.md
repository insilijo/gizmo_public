# GIZMO — public repository

A biochemistry substrate as a fixed coordinate system for per-patient multi-omic state, with extensions across multiple papers.

This repository hosts (i) the shared substrate + framework code, and (ii) self-contained per-paper deposits. Each `papers/paperN/` folder contains the manuscript, supplement, figures, scripts, and result deposits for one paper. The 38,148-node biochemistry substrate (`substrate/graph.json`) is the Paper 1 deposit and is licensed CC-BY 4.0.

---

## Repository layout

```
gizmo_public/
├── README.md              # This file — portal + framework intro
├── LICENSE
├── requirements.txt
├── run_paper1.py          # End-to-end pipeline runner
├── substrate/             # SHARED — the 38,148-node biochemistry substrate
│   └── graph.json
├── gizmo/                 # SHARED — Python package (MAP, β/α, signed-basin, scoring)
├── benchmarks/            # SHARED — analysis machinery (MOFA baselines, per_patient_master, diagnostics, figure builders)
├── docs/                  # SHARED — INSTALL, USAGE, design docs
└── papers/
    └── paper1/            # Paper 1: A biochemistry substrate as a fixed coordinate system
        ├── README.md
        ├── MANUSCRIPT.{md,docx}
        ├── SUPPLEMENT.{md,docx}
        ├── FIGURES.{md,docx}
        ├── figures/       # 11 main figures
        ├── scripts/       # 19 manuscript-load-bearing analysis scripts + 2 WGCNA legacy
        ├── results/       # 20 result deposits (JSONs + TSVs)
        └── data/curation/ # Curated key-gene ground truth
```

Future papers (Paper 2, 3, 4, 5; see below) will each get their own self-contained `papers/paperN/` folder following this same structure.

---

## Papers in this repository

| Paper | Folder | Status | Topic |
|---|---|---|---|
| **Paper 1** | [`papers/paper1/`](papers/paper1/) | **v7, accepted** (5 review rounds) | Biochemistry substrate as a fixed coordinate system for per-patient multi-omic state; quantitative interpretability evaluation; substrate-anchored cross-disease activation vocabulary |
| Paper 2 | *future* | scoped | Lipid layer dual-architecture (flux + state sub-graphs, compositional Bayes for ionization bias, degree-weighted MAP) |
| Paper 3 | *future* | scoped | Disease-informed topology rewiring (network deformation under disease state) |
| Paper 4 | *future* | scoped | Mutation-conditional GoF edge injection (fix for the structural constraint in Paper 1 §4.6) |
| **Paper 5** | [`papers/paper5/`](papers/paper5/) | **v0.1, scaffolded** (5-cohort longitudinal calibration audited; cross-cohort axis projection stress-validated; manuscript prose not yet written) | Longitudinal patient-state tracking under intervention + cross-cohort substrate-axis projection methodology; spatial-tile extension scoped for follow-up |

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

# 3. End-to-end Paper 1 pipeline on one cohort
python3 run_paper1.py Crohn

# 4. Replay β/α + basins from a Zenodo-deposited F matrix
python3 run_paper1.py --from-F /path/to/zenodo/F.npz --out results/Crohn_replay

# 5. Reproduce a v7 manuscript number
PYTHONPATH=. python3 papers/paper1/scripts/check_v7_interpretability_eval.py
PYTHONPATH=. python3 papers/paper1/scripts/check_v7_basin_conservation.py
# ... full mapping in papers/paper1/scripts/README.md
```

Full install instructions: [`docs/INSTALL.md`](docs/INSTALL.md). Step-by-step walkthrough: [`docs/USAGE.md`](docs/USAGE.md). Paper 1 manuscript: [`papers/paper1/MANUSCRIPT.md`](papers/paper1/MANUSCRIPT.md).

---

## What the substrate is

A fixed 38,148-node biochemistry graph (12,872 gene nodes + 11,244 metabolite nodes + 14,032 reaction nodes), merged from four curated sources:

- **Reactome** (release 87, June 2024) — 14,032 reactions
- **StringDB** (v12.0) — protein-protein interactions, confidence ≥ 700
- **HMDB** (v5.0) — metabolite identifiers; de-duplicated by InChIKey14
- **KEGG** (2024-Q4) — pathway annotations

Any patient's multi-omic measurements MAP-project onto this substrate, producing a per-patient state vector **F** in 38,148-dim substrate coordinates that don't depend on cohort or modality. The substrate is the coordinate system; downstream method choices (GIZMO's MAP + β/α decomposition + signed basins, or alternatives) operate within it.

See [`papers/paper1/MANUSCRIPT.md`](papers/paper1/MANUSCRIPT.md) §2 for substrate construction, [`papers/paper1/SUPPLEMENT.md`](papers/paper1/SUPPLEMENT.md) §S1.1 for full curation details, and [`substrate/README.md`](substrate/README.md) for a focused substrate-only intro.

---

## Per-cohort case-study deposit (Zenodo)

Every cohort GIZMO has been applied to is deposited at Zenodo with the full analytical bundle: F matrix, β/α decomposition, signed-basin output, per-cohort metadata, and an auto-populated case-study writeup. The bundle is regenerated by `benchmarks/make_zenodo_bundle.py`. See [`PRE_SUBMISSION_CHECKLIST.md`](PRE_SUBMISSION_CHECKLIST.md) for upload status.

Zenodo DOI: *TBD on acceptance.*

---

## License

- **Substrate + analytical outputs**: CC-BY 4.0.
- **Code**: Apache 2.0.
- **Upstream source data** (per-cohort raw inputs): governed by each source's original license; not redistributed here. See [`LICENSE`](LICENSE).

If you use the substrate or code, please cite the relevant paper. For Paper 1:

> Gardner JJ (2026). *GIZMO: a biochemistry substrate as a fixed coordinate system for per-patient multi-omic state, with a deposited reference library across 12 cohorts.* [submitted]; substrate v1.0 at https://doi.org/[Zenodo DOI TBD].
