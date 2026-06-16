"""Test β vs severity-like continuous outcomes across multiple cohorts."""
import sys, json, openpyxl
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr
sys.path.insert(0,'/home/jgardner/GIZMO')
from gizmo.export.json_export import read_json
from gizmo.inference.projection import build_biochem_subgraph
from sklearn.decomposition import PCA

REPO = Path('/home/jgardner/GIZMO')
SNAPSHOT = REPO / 'benchmarks/results/unsupervised/_pre_zscore_snapshot_20260607'

mg = read_json(REPO / 'data/processed/human_full/graph.json')
geom = build_biochem_subgraph(mg, hub_cap=200)
log_pr = geom.log_pr


def compute_beta(F):
    F_unit = F / (np.linalg.norm(F, axis=1, keepdims=True) + 1e-12)
    x = log_pr; xm=x.mean(); xv=x.var()+1e-12
    Fm = F_unit.mean(axis=1, keepdims=True)
    cov = ((F_unit - Fm)*(x-xm)).mean(axis=1, keepdims=True)
    return (cov/xv).ravel()


def load_F_pids(cohort, lower=False):
    fp = SNAPSHOT / f'stage3_F_{cohort}_edge_informed.npz'
    if not fp.exists():
        fp = SNAPSHOT / f'stage3_F_{cohort}.npz'
    if not fp.exists(): return None, None
    npz = np.load(fp, allow_pickle=True)
    F = npz['F']
    if F.shape[1] != len(geom.nodes): return None, None
    pids = [str(p).lower() if lower else str(p) for p in npz['patient_ids']]
    return F, pids


# === Filbin (control test - already known) ===
print('=== Filbin_COVID: β vs Acuity max ===')
F, pids = load_F_pids('Filbin_COVID')
beta = compute_beta(F)
wb = openpyxl.load_workbook('/home/jgardner/.cache/filbin_mgh_covid/Clinical_Metadata.xlsx', read_only=True)
rows = list(wb['Subject-level metadata'].iter_rows(values_only=True))
clin = pd.DataFrame(rows[1:], columns=rows[0])
clin['patient_id'] = clin['Public ID'].astype(str)
for sev_col in ['Acuity 0', 'Acuity 3', 'Acuity 7', 'Acuity 28', 'Acuity max']:
    clin[sev_col + '_num'] = pd.to_numeric(clin[sev_col], errors='coerce')
clin = clin.set_index('patient_id')
common = sorted(set(pids) & set(clin.index))
pid_idx = {p:i for i,p in enumerate(pids)}
idx = np.array([pid_idx[p] for p in common])
b_vals = beta[idx]
for sev_col in ['Acuity 0', 'Acuity 3', 'Acuity 7', 'Acuity 28', 'Acuity max']:
    sev = clin.loc[common][sev_col + '_num'].values
    mask = ~np.isnan(sev)
    if mask.sum() < 30: continue
    rho, p = spearmanr(b_vals[mask], sev[mask])
    print(f'  β vs {sev_col}:  ρ = {rho:+.3f}  p = {p:.4f}  (n = {mask.sum()})')
print('  [Acuity: 1=died, 5=discharged, so negative ρ = high β tracks sicker]')

# === Su_COVID: WHO Ordinal Scale ===
print('\n=== Su_COVID: β vs WHO Ordinal Scale ===')
F, pids = load_F_pids('Su_COVID')
beta = compute_beta(F)
wb = openpyxl.load_workbook('/home/jgardner/GIZMO/benchmarks/Table S1. Human subject details, plasma proteomic and metabolomic datasets and analysis, and CITE-seq antibodies. Related to Figures 1 and S1.xlsx', read_only=True)
ws = wb['S1.1 Patient Clinical Data']
rows = list(ws.iter_rows(values_only=True))
clin = pd.DataFrame(rows[1:], columns=rows[0])
clin['who'] = pd.to_numeric(clin['Who Ordinal Scale'], errors='coerce')
# Su pids look like INCOV001-1 (with timepoint suffix); F pids?
print(f'  Su F pids sample (5): {pids[:5]}')
print(f'  Su clin Sample ID sample (5): {clin["Sample ID"].iloc[:5].tolist()}')

# Try Sample ID matching
clin['Sample ID'] = clin['Sample ID'].astype(str)
clin = clin.set_index('Sample ID')
common = sorted(set(pids) & set(clin.index))
if len(common) > 30:
    pid_idx = {p:i for i,p in enumerate(pids)}
    idx = np.array([pid_idx[p] for p in common])
    b_vals = beta[idx]
    who = clin.loc[common]['who'].values
    mask = ~np.isnan(who)
    rho, p = spearmanr(b_vals[mask], who[mask])
    print(f'  β vs WHO ordinal:  ρ = {rho:+.3f}  p = {p:.4f}  (n = {mask.sum()})')
    print(f'  [WHO higher = sicker, so positive ρ = high β tracks sicker]')
else:
    print(f'  ID format mismatch — common = {len(common)} (need 30+)')

# === TCGA_LUAD: pathologic stage ===
print('\n=== TCGA_LUAD: β vs pathologic stage ===')
F, pids = load_F_pids('TCGA_LUAD', lower=True)
beta = compute_beta(F)
CLIN = Path.home()/'.cache'/'tcga_luad'/'gdac.broadinstitute.org_LUAD.Clinical_Pick_Tier1.Level_4.2016012800.0.0'/'LUAD.clin.merged.picked.txt'
cdf = pd.read_csv(CLIN, sep='\t', header=None, low_memory=False)
attrs = cdf.iloc[:,0].astype(str).tolist()
clin_pids = [str(p).strip().lower() for p in cdf.iloc[0,1:].tolist()]
stage_str = [str(v).strip() for v in cdf.iloc[attrs.index('pathologic_stage'),1:].tolist()]
def parse_stage(s):
    if not isinstance(s,str): return None
    s = s.lower()
    if 'iv' in s: return 4
    if 'iii' in s: return 3
    if 'ii' in s and 'iii' not in s: return 2
    if 'i' in s and 'ii' not in s: return 1
    return None
stage_map = {p: parse_stage(s) for p, s in zip(clin_pids, stage_str)}
common = sorted(set(pids) & set(stage_map))
pid_idx = {p:i for i,p in enumerate(pids)}
idx = np.array([pid_idx[p] for p in common])
b_vals = beta[idx]
stage = np.array([stage_map[p] for p in common], dtype=float)
mask = ~np.isnan(stage)
rho, p = spearmanr(b_vals[mask], stage[mask])
print(f'  β vs stage (1-4):  ρ = {rho:+.3f}  p = {p:.4f}  (n = {mask.sum()})')
print('  [stage higher = more advanced]')

# === TCGA_IDH: WHO grade ===
print('\n=== TCGA_IDH_glioma: β vs WHO grade ===')
F, pids = load_F_pids('TCGA_IDH_glioma', lower=True)
beta = compute_beta(F)
CLIN = Path.home()/'.cache'/'tcga_idh'/'lgggbm_tcga_pub_clinical.tsv'
c = pd.read_csv(CLIN, sep='\t')
c['patient_id'] = c['patient_id'].astype(str).str.lower()
# Grade column
def parse_grade(s):
    if not isinstance(s, str): return None
    s = s.upper().strip()
    if 'IV' in s: return 4
    if 'III' in s: return 3
    if 'II' in s: return 2
    return None
c['grade_num'] = c['GRADE'].apply(parse_grade)
c = c.set_index('patient_id')
common = sorted(set(pids) & set(c.index))
pid_idx = {p:i for i,p in enumerate(pids)}
idx = np.array([pid_idx[p] for p in common])
b_vals = beta[idx]
grade = c.loc[common]['grade_num'].values.astype(float)
mask = ~np.isnan(grade)
rho, p = spearmanr(b_vals[mask], grade[mask])
print(f'  β vs WHO grade:  ρ = {rho:+.3f}  p = {p:.4f}  (n = {mask.sum()})')
print('  [grade higher = more aggressive]')

# === KMPLOT_BRCA: histological grade ===
print('\n=== KMPLOT_BRCA: β vs grade + ER status ===')
F, pids = load_F_pids('KMPLOT_BRCA')
beta = compute_beta(F)
KS = Path('/home/jgardner/gitlab-old/d2f483672c0239f6d7dd3c9ecee6deacbcd59185855625902a8b1c1a3bd67440/KMPLOT_BRCA_EXPRESSION/DATA/KMPLOT_BRCA_SURVIVAL.txt')
s = pd.read_csv(KS, sep='\t')
s['AffyID'] = s['AffyID'].astype(str)
s = s.set_index('AffyID')
common = sorted(set(pids) & set(s.index))
pid_idx = {p:i for i,p in enumerate(pids)}
idx = np.array([pid_idx[p] for p in common])
b_vals = beta[idx]
# Show columns to find grade
print('  KMPLOT survival columns:', list(s.columns))
grade = pd.to_numeric(s.loc[common]['Grade'], errors='coerce').values
mask = ~np.isnan(grade)
if mask.sum() > 30:
    rho, p = spearmanr(b_vals[mask], grade[mask])
    print(f'  β vs Grade (1-3): ρ = {rho:+.3f}  p = {p:.4f}  (n = {mask.sum()})')

