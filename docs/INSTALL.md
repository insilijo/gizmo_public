# Install

GIZMO Paper 1 codebase. Python 3.10+ required.

## 1. Clone the repository

```bash
git clone https://github.com/insilijo/gizmo_public.git
cd gizmo_public
```

If the substrate file is missing or is a Git LFS pointer (size <1 KB), pull the LFS object:

```bash
git lfs install
git lfs pull
```

The substrate file is `substrate/graph.json` and should be ~89 MB when fully present.

## 2. Create a Python environment

We tested on Python 3.10. A clean venv is recommended.

```bash
python3 -m venv .venv
source .venv/bin/activate   # bash/zsh
pip install --upgrade pip
```

## 3. Install runtime dependencies

```bash
pip install -r requirements.txt
```

This installs numpy, scipy, networkx, pandas, scikit-learn, mofapy2, matplotlib, Pillow, statsmodels at pinned major-version bounds.

## 4. Install the `gizmo` Python package in editable mode

The `gizmo/` directory is the source of the package. From the repo root:

```bash
pip install -e .
```

If `pyproject.toml`/`setup.py` is missing (the repo is intentionally minimal), use the workaround:

```bash
export PYTHONPATH="$(pwd):$PYTHONPATH"
```

Add to your shell rc to persist.

## 5. Smoke test — load the substrate

```bash
python3 -c "
from gizmo.export.json_export import read_json
mg = read_json('substrate/graph.json')
print(f'substrate loaded: {mg.graph.number_of_nodes()} nodes, '
      f'{mg.graph.number_of_edges()} edges')
"
```

Expected output:

```
substrate loaded: 38148 nodes, ~170000 edges
```

If you see this, the install is good.

## 6. (Optional) Verify benchmark imports

The benchmark scripts in `benchmarks/` use relative imports from the repo root. Smoke test:

```bash
python3 -c "
import sys; sys.path.insert(0, 'benchmarks')
import baseline_mofa_streaming
print('benchmarks import OK')
"
```

## 7. Cohort data

This repository does **not** include the per-cohort raw input data — those are governed by each cohort's source license (TCGA, GEO, CPTAC, NEPTUNE, etc.). The Zenodo deposit (DOI TBD on paper acceptance) ships the per-cohort processed F matrices + β/α decomposition + signed-basin outputs for direct use. To reproduce from scratch:

1. Download each cohort's raw input data from its original source (see `MANUSCRIPT.md` §Methods §"Cohort panel" + `benchmarks/per_patient_master.py` loader docstrings for source URLs).
2. Run the per-cohort GIZMO MAP solve to produce the cohort's F matrix (see `docs/USAGE.md` walkthrough).

## Troubleshooting

- **`mofapy2` install fails**: known issue with some pip versions. Try `pip install mofapy2==0.7.1 --no-deps`, then install its deps manually (`pandas`, `numpy`, `scipy`).
- **`scipy.sparse.linalg.cg` reports non-convergence on MAP**: increase `maxiter` in `gizmo/inference/map_solve.py` from default 500 to 2000, or check that the substrate file isn't truncated (size ~89 MB).
- **`networkx` 3.6 warns about `edges` kwarg**: cosmetic; the code uses the legacy convention. Suppress via `import warnings; warnings.filterwarnings('ignore', category=FutureWarning)`.
- **Out of memory on MOFA+ for large cohorts**: use the streaming runner (`benchmarks/baseline_mofa_streaming.py`) which implements IncrementalPCA for single-block + patient-subsample MOFA+ for multi-block. See `docs/USAGE.md` §"MOFA+ comparison" for invocation.

For all other issues, open a GitHub issue with the error trace + Python version + OS.
