# Usage — reproducing the paper's claims end-to-end

This walkthrough takes one cohort (Crohn, n = 33, smallest in the panel) and reproduces the manuscript's Crohn α-PC1 thiopurine pharmacogenomic finding (Manuscript §2). The same recipe applies to every other cohort with the cohort name substituted.

Repeat for any other cohort by changing `--cohort Crohn` to the target cohort. Cohort names: `IDH_glioma`, `TCGA_IDH_glioma`, `TCGA_LUAD`, `KMPLOT_BRCA`, `CPTAC_CCRCC`, `CPTAC_COAD`, `CPTAC_OV`, `Su_COVID`, `Filbin_COVID`, `GSE65391_SLE`, `GSE65682_sepsis`, `Gao_RA`, `GSE89408_RA`, `Crohn`, `HMP2_IBD_CD`, `Erawijantari`, `NEPTUNE_kidney`, `Wang_RA`, `TB_DX`.

## Step 0 — install

Follow [INSTALL.md](INSTALL.md). Skip ahead if `python3 -c "from gizmo.export.json_export import read_json; mg = read_json('substrate/graph.json'); print(mg.graph.number_of_nodes())"` prints `38148`.

## Step 1 — load the cohort

Each cohort has a loader in `benchmarks/per_patient_master.py` that returns proteomics dict + metabolomics dict + label dict + sample list. For Crohn:

```python
import sys; sys.path.insert(0, "benchmarks")
from per_patient_master import load_crohn
prot, metab, ylab, samples = load_crohn()
print(f"Crohn loaded: {len(samples)} samples, "
      f"{len(prot[samples[0]])} prot features, "
      f"{len(metab[samples[0]])} metab features")
```

If the cohort data file is missing, the loader will raise; see `benchmarks/per_patient_master.py` source for the expected file path and provenance.

## Step 2 — MAP-solve to produce F matrix

The MAP solve (`F = argmin_F ‖x − A_obs F‖²_Σ⁻¹ + λ·F^T L_signed F`) takes the per-patient input and projects onto the substrate's 38,148-dim coordinate system. The solver is in the `gizmo` package; in the published per-cohort pipeline scripts (`data/cohorts/<cohort>/run_gizmo_full.py` in the GIZMO research repo) this is wrapped per cohort. For an interactive solve:

```python
import numpy as np
from gizmo.export.json_export import read_json
from gizmo.evidence.mappers import GeneMapper, MetaboliteMapper
from gizmo.inference.map_solve import solve_per_patient_F

mg = read_json("substrate/graph.json")
gmap, mmap = GeneMapper(mg), MetaboliteMapper(mg)

F_per_patient = {}
for sid in samples:
    F_per_patient[sid] = solve_per_patient_F(
        prot.get(sid, {}), metab.get(sid, {}),
        mg=mg, gene_mapper=gmap, metab_mapper=mmap,
        lam=1.0)
```

For 38,148-node MAP at ~30 sec/patient, n = 33 takes ~17 min. Output is a per-patient 38,148-dim vector in substrate coordinates.

Save the F matrix:

```python
F = np.stack([F_per_patient[s] for s in samples])
np.savez_compressed("results/F_Crohn.npz",
                     F=F, patient_ids=np.array(samples, dtype=object))
```

The published Zenodo deposit ships every cohort's F matrix at this format; rebuild only if you've changed substrate or solver settings.

## Step 3 — β/α decomposition

Given F, compute β (hub-projection scalar) and α (orthogonal residual):

```python
import networkx as nx

# log-PageRank direction — depends only on substrate, not cohort
sub = mg.graph.to_undirected() if mg.graph.is_directed() else mg.graph
nodes = list(sub.nodes())
pr = nx.pagerank(sub)
log_pr = np.log10(np.array([pr.get(n, 0.0) for n in nodes]) + 1e-15)

# β = scalar projection; α = orthogonal residual
F_unit = F / (np.linalg.norm(F, axis=1, keepdims=True) + 1e-12)
x = log_pr; x_mean = x.mean(); x_var = x.var() + 1e-12
F_mean = F_unit.mean(axis=1, keepdims=True)
cov = ((F_unit - F_mean) * (x - x_mean)).mean(axis=1, keepdims=True)
beta = (cov / x_var).ravel()                # (n_patients,)
alpha = F_unit - F_mean - beta[:, None] * (x - x_mean)[None, :]
```

## Step 4 — α-PCA

```python
from sklearn.decomposition import PCA
pca = PCA(n_components=5, random_state=0)
scores = pca.fit_transform(alpha)            # (n_patients, 5)
loadings = pca.components_                    # (5, 38148)
```

Each row of `loadings` is an α-PC — a unit vector in substrate coordinates. The k-th α-PC's loading on the i-th substrate node is `loadings[k-1, i]`.

## Step 5 — signed-basin decomposition (the headline output)

For each α-PC, partition substrate nodes by loading sign, find the largest connected component in each sign-class. This produces the named-biology layer:

```python
import networkx as nx

def signed_basins(pc, sub, nodes):
    pos_sub = sub.subgraph([nodes[i] for i in range(len(nodes)) if pc[i] > 0])
    neg_sub = sub.subgraph([nodes[i] for i in range(len(nodes)) if pc[i] < 0])
    pos_cc = max(nx.connected_components(pos_sub), key=len, default=set())
    neg_cc = max(nx.connected_components(neg_sub), key=len, default=set())
    return pos_cc, neg_cc

pos_cc, neg_cc = signed_basins(loadings[0], sub, nodes)   # α-PC1
# Top + basin members by |loading|:
idx_of = {n: i for i, n in enumerate(nodes)}
top_pos = sorted(pos_cc, key=lambda n: -abs(loadings[0, idx_of[n]]))[:10]
for n in top_pos:
    attrs = mg.graph.nodes.get(n, {})
    sym = attrs.get("symbol") or attrs.get("name") or n
    print(f"  {sym:30}  loading={loadings[0, idx_of[n]]:+.3f}")
```

For Crohn α-PC1, this should print MPG at the top with loading ~+0.52 (the canonical thiopurine pharmacogenomic gene), followed by MPG-mediated reactions (3-methyladenine cleavage, hypoxanthine cleavage, ethenoadenine cleavage, APEX1 displacement). This is the Manuscript §2 Crohn finding — rediscovered without supervision.

## Step 6 — MOFA+ comparison (§4)

For symmetric multi-axis discrimination testing, run MOFA+ on the same cohort and compare per-axis discrimination against clinical metadata.

```bash
# Standard MOFA+ (memory-permitting on modest cohorts)
python3 benchmarks/baseline_mofa_with_weights.py Crohn

# Memory-bounded streaming variant (large microarray / multi-block cohorts)
python3 benchmarks/baseline_mofa_streaming.py Crohn

# Substrate-matched feature universe (apples-to-apples vs GIZMO)
python3 benchmarks/baseline_mofa_streaming.py Crohn --substrate-matched
```

Output lands in `benchmarks/results/unsupervised/mofa_weights/mofa_weights_<cohort>{_sm}.json`.

After running for the panel cohorts you care about, run the symmetric comparison pipeline:

```bash
python3 benchmarks/diagnostics/multi_pc_vs_mofa_factors.py
python3 benchmarks/diagnostics/multi_pc_variance_bonferroni.py
python3 benchmarks/diagnostics/multi_pc_vs_mofa_substrate_matched.py
```

This produces `multi_pc_vs_mofa_factors_augmented.tsv` (asymmetric: MOFA+ at full input) and `multi_pc_vs_mofa_substrate_matched.tsv` (apples-to-apples). Manuscript §4 reports 8/14 vs 9/14 global Bonferroni asymmetric and 6/11 vs 6/11 substrate-matched.

## Step 7 — figures

```bash
# Figure 2a — cross-cohort chord diagram
python3 benchmarks/figures/build_pc_alignment_chord.py

# Figure 2b–2f — individual basin diagrams per cohort
python3 benchmarks/figures/build_basin_signed_v2.py

# Figure 2g — chord co-structure (network + dot plot + hypergeometric pathway)
python3 benchmarks/figures/build_chord_costructure.py

# Figure 2h — cross-disease conserved biology + companion TSVs
# (produced by build_chord_costructure.py as a side effect)
```

Output PNGs go to `benchmarks/results/figures/`; copy to `figures/main/` for manuscript inclusion.

## Step 8 — per-cohort Zenodo case-study bundle

Once F + β/α + MOFA+ outputs are in place, generate the per-cohort case-study deposit:

```bash
python3 benchmarks/make_zenodo_bundle.py --cohorts Crohn        # one cohort
python3 benchmarks/make_zenodo_bundle.py                         # all 19
```

Output lands in `zenodo_deposit/cohorts/<cohort>/` with the F matrix, β/α per patient, signed-basin TSVs per PC, MOFA+ comparison weights, metadata, and an auto-populated case-study README ready for hand-edit before final Zenodo upload.

## Manuscript claims by step

| Claim | Section | Step(s) |
|---|---|---|
| 38,148-node substrate, four-source merge, CC-BY 4.0 | §1 | Step 0 (substrate file) |
| Per-patient F in substrate coordinates | §1 | Steps 1–2 |
| β = phenotype-presentation magnitude (47/50 top hubs T/S/I) | §1 | Step 3 |
| Signed-basin decomposition per α-PC | §2 | Steps 4–5 |
| Crohn α-PC1 MPG / thiopurine | §2 | Step 5 (cohort = Crohn) |
| ccRCC Warburg/HIF | §2 | Step 5 (cohort = CPTAC_CCRCC) |
| SLE α-PC3 type-I IFN | §2 | Step 5 (cohort = GSE65391_SLE) |
| Filbin α-PC4 LPO/galanin | §2 | Step 5 (cohort = Filbin_COVID) |
| Cross-cohort α-PC cosine alignment | §2 | Step 7 (chord) |
| Horizontal meta-analysis joint p ≈ 2×10⁻³ | §3 | Steps 1–3 across IDH+HMP2+GSE89408 |
| TopPR baseline + GIZMO ≈ noisy TopPR | §4 | (separate benchmark, see `benchmarks/operating_curve_full.py`) |
| MOFA+ parity: 6/11 vs 6/11 substrate-matched | §4 | Step 6 |
| GoF falsification on TCGA_LUAD | §5 | `benchmarks/diagnostics/luad_gof_verification.py` |

## Common gotchas

- **Cohort metadata loaders may need local data files** that aren't in this repo (e.g., `data/cohorts/GSE65391_SLE/GSE65391_series_matrix.txt.gz`). Each loader's docstring lists what it expects.
- **PC sign canonicalization** is needed for any cross-cohort per-node sign comparison; see Methods §"PC sign canonicalization." If `vᵃ · vᵇ < 0`, flip `vᵦ ← −vᵦ` before per-node display.
- **Cosine direction** at high \|cos\| is the *biology* — basin assignment can flip if cohort B's PC was solver-flipped; canonicalize before reading basins side-by-side.
- **F matrix file naming** — some cohorts saved their F with suffix (`stage3_F_TCGA_IDH_glioma_edge_informed.npz` vs `stage3_F_TCGA_IDH_glioma.npz`). The script `benchmarks/diagnostics/multi_pc_vs_mofa_factors.py` already handles fallback; check `find_F_path()` if writing new tooling.
