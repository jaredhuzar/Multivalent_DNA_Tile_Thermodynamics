#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Regenerate the CSVs for SI Appendix Tables S1 and S4, and verify them against the
published values in PNAS_SI_Appendix_clean.docx.

  Table S1  Model performances   (cross-validated RMSD and Spearman rho of 25 models)
  Table S4  Fitted parameters     (MV3 multivalent model + two-state A2, with s.e. and 95% CI)

Every NUMBER is read from the pipeline outputs that run_all.sh regenerates, so the tables
refresh automatically.  The descriptive text columns (model name, formula, parameter count,
meaning) are the fixed definitions from the SI Appendix and are carried verbatim.

Number sources (all produced earlier in run_all.sh):
  MV1-5, MV1d-3d, A1-A7, B1-B4 : Supplementary/TableS1_Model_Comparison_All/TableS1_all_models.csv
                                   (make_tableS1_combined.py <- unified_models.py, Figure4 fit)
  M0 (additive NN two-state)   : Supplementary/PredictedVsMeasured_NN_Multivalent/data.csv
  ML rows (5)                  : Figure6_MachineLearning/Fig6d_ModelComparison_SpearmanRho/
                                   all_models_incl_summedNN.csv  (make_figure6.py)
  Table S4 parameters          : Supplementary/TableS4_Parameter_CIs/TableS4_parameters.csv
                                   (make_tableS4.py)

Outputs (written to <repo>/Pseudo_Figures/SI/, the SI deliverables folder):
  TableS1_ModelPerformance.csv
  TableS4_FittedParameters.csv

Verification ("check if the scripts produce the same numbers"): each regenerated number is
compared with the value printed in PNAS_SI_Appendix_clean.docx (hard-coded in REF_S1 / REF_S4
at the document's displayed precision).  A per-row PASS/FAIL report is printed.  The script
exits 0 (report only) unless --strict is given, in which case any mismatch exits 1 (for CI).

Run from anywhere:  python3 Supplementary/SI_Tables/make_si_tables.py [--strict]
"""
import os, re, sys, math
import numpy as np
import pandas as pd

# ----------------------------------------------------------------------------- repo root
def _find_repo():
    d = os.path.dirname(os.path.abspath(__file__))
    while d != '/' and not os.path.exists(os.path.join(d, 'Filtered_DF_Tm.csv')):
        d = os.path.dirname(d)
    return d
REPO = _find_repo()
OUTDIR = os.path.join(REPO, 'Pseudo_Figures', 'SI')
os.makedirs(OUTDIR, exist_ok=True)
STRICT = '--strict' in sys.argv[1:]

def _req(path):
    p = os.path.join(REPO, path)
    if not os.path.exists(p):
        sys.exit(f"ERROR: required source not found: {path}\n  (run run_all.sh first so the "
                 f"pipeline writes it)")
    return p

# ----------------------------------------------------------------------------- load sources
TS1 = pd.read_csv(_req('Supplementary/TableS1_Model_Comparison_All/TableS1_all_models.csv'))
TS1 = TS1.set_index('id')                                    # cvRMS, cvrho, n_params keyed by id
ML  = pd.read_csv(_req('Figure6_MachineLearning/Fig6d_ModelComparison_SpearmanRho/'
                       'all_models_incl_summedNN.csv')).set_index('Model')   # OOF_RMS_C, OOF_Rho
NN  = pd.read_csv(_req('Supplementary/PredictedVsMeasured_NN_Multivalent/data.csv'))

# M0: additive nearest-neighbour two-state Tm vs measured (sequence-only, no fitted parameters)
from scipy.stats import spearmanr
_a, _m = NN['Tm_NN_additive_C'].values, NN['Tm_meas_C'].values
M0_RMS = float(np.sqrt(np.mean((_a - _m) ** 2)))
M0_RHO = float(spearmanr(_a, _m)[0])

def ts1(idx):  return float(TS1.loc[idx, 'cvRMS']), float(TS1.loc[idx, 'cvrho'])
def ml(name):  return float(ML.loc[name, 'OOF_RMS_C']), float(ML.loc[name, 'OOF_Rho'])

# ----------------------------------------------------------------------------- Table S1 spec
# (Row, Family, Model, Weighting and tethering, N_params, getter)  -- order follows the SI doc.
MV = 'Multivalent binding models (this work)'
TS = 'Two-state models (Tm = dH_eff / dS_eff)'
MLF = 'Machine learning (11-site nearest-neighbor dG encoding)'
S1_SPEC = [
    ('MV1',  MV, 'One uniform C_eff',                 'C_eff,i = C_ref',                              '1', lambda: ts1('M1')),
    ('MV2',  MV, '+ distance',                        'C_eff,i = C_ref*exp(-d_i/lambda)',             '2', lambda: ts1('M2')),
    ('MV3',  MV, '+ valency (model of Fig. 4)',       'C_eff,i = C_ref*exp(-d_i/lambda)*exp(-kappa(N-4))', '3', lambda: ts1('M3')),
    ('MV4',  MV, 'Bond-type C_eff instead of kappa',  'C_eff,i = C_ref*phi_type*exp(-d_i/lambda)',    '4', lambda: ts1('M4')),
    ('MV5',  MV, 'Bond-type C_eff and kappa',         'C_eff,i = C_ref*phi_type*exp(-d_i/lambda)*exp(-kappa(N-4))', '5', lambda: ts1('M5')),
    ('MV1d', MV, 'MV1, docking readout',              'theta = f (w = 1)',                            '1', lambda: ts1('M1d')),
    ('MV2d', MV, 'MV2, docking readout',              'theta = f (w = 1)',                            '2', lambda: ts1('M2d')),
    ('MV3d', MV, 'MV3, docking readout',              'theta = f (w = 1)',                            '3', lambda: ts1('M3d')),
    ('M0',   TS, 'Sequence-only sum',                 'none (no free parameters)',                    '0', lambda: (M0_RMS, M0_RHO)),
    ('B1',   TS, 'M1, independent sites',             'g_k = 1; one dS_bond per bond',                '2', lambda: ts1('B1')),
    ('A2',   TS, 'M3, occupancy weighting',           'g_k = exp(-d_k/lambda) on dH_NN and dS_NN',    '3', lambda: ts1('A2')),
    ('A1',   TS, 'dS_bond inside the sum',            'g_k applied to dS_NN,k + dS_bond',             '3', lambda: ts1('A1')),
    ('A3',   TS, 'Enthalpy weighting only',           'g_k on dH_NN,k only',                          '3', lambda: ts1('A3')),
    ('A4',   TS, 'Nearest-bond weighting',            'g_k = exp(-g_k^NN/lambda)',                    '3', lambda: ts1('A4')),
    ('A5',   TS, 'Per-length tethering',              'M3 with one dS_bond per bond length',          '4', lambda: ts1('A5')),
    ('A6',   TS, 'Split decay lengths',               'separate lambda for dH_NN and dS_NN',          '4', lambda: ts1('A6')),
    ('A7',   TS, 'M3 + rotational penalty',           'M3 plus a distance-weighted rotational term',  '4', lambda: ts1('A7')),
    ('B2',   TS, 'M1 + rotational penalty',           'g_k = 1 plus the rotational term',             '3', lambda: ts1('B2')),
    ('B3',   TS, 'Weighted + WLC tether',             'g_k = exp(-d_k/lambda); WLC tethering entropy','2', lambda: ts1('B3')),
    ('B4',   TS, 'WLC tether only',                   'g_k = 1; WLC tethering entropy',               '2', lambda: ts1('B4')),
    ('-',    MLF,'Summed NN dG (linear in the sum)',  'sum of the 11 site dG values',                 '2', lambda: ml('Summed NN ΔG')),
    ('-',    MLF,'Linear regression',                 '11-element encoding',                          '-', lambda: ml('Linear regression')),
    ('-',    MLF,'Quadratic regression',              'encoding plus pairwise products',              '-', lambda: ml('Quadratic regression')),
    ('-',    MLF,'Random forest',                     '11-element encoding',                          '-', lambda: ml('Random forest')),
    ('-',    MLF,'XGBoost',                           '11-element encoding',                          '-', lambda: ml('XGBoost')),
]

# Published values in PNAS_SI_Appendix_clean.docx: key -> (RMSD, rmsd_dp, rho, rho_dp)
REF_S1 = {
    'MV1':(2.40,2,0.824,3),'MV2':(1.57,2,0.933,3),'MV3':(1.34,2,0.934,3),'MV4':(1.34,2,0.935,3),
    'MV5':(1.34,2,0.935,3),'MV1d':(2.22,2,0.831,3),'MV2d':(1.51,2,0.932,3),'MV3d':(1.41,2,0.931,3),
    'M0':(23.9,1,0.36,2),'B1':(2.26,2,0.783,3),'A2':(1.30,2,0.933,3),'A1':(1.33,2,0.930,3),
    'A3':(1.31,2,0.934,3),'A4':(1.98,2,0.840,3),'A5':(1.33,2,0.932,3),'A6':(1.29,2,0.932,3),
    'A7':(1.27,2,0.937,3),'B2':(1.37,2,0.933,3),'B3':(3.18,2,0.602,3),'B4':(3.17,2,0.613,3),
    'Summed NN dG':(4.10,2,0.09,2),'Linear regression':(3.25,2,0.66,2),
    'Quadratic regression':(2.98,2,0.81,2),'Random forest':(2.10,2,0.84,2),'XGBoost':(1.72,2,0.90,2),
}
def _refkey(row, model):           # ML rows share Row='-', so key them by model name
    return row if row != '-' else model

# ----------------------------------------------------------------------------- build + check S1
s1_rows, checks = [], []
for row, fam, model, weight, npar, get in S1_SPEC:
    rms, rho = get()
    s1_rows.append(dict(Family=fam, Row=row, Model=model, Weighting_and_tethering=weight,
                        N_params=npar, RMSD_C=round(rms, 4), rho=round(rho, 4)))
    key = _refkey(row, model)
    if key in REF_S1:
        r_rms, dr, r_rho, dp = REF_S1[key]
        checks.append((f'S1 {key:<22s} RMSD', rms, r_rms, dr))
        checks.append((f'S1 {key:<22s} rho ', rho, r_rho, dp))

pd.DataFrame(s1_rows).to_csv(os.path.join(OUTDIR, 'TableS1_ModelPerformance.csv'), index=False)

# ----------------------------------------------------------------------------- Table S4
S4 = pd.read_csv(_req('Supplementary/TableS4_Parameter_CIs/TableS4_parameters.csv'))

def _first_float(s):
    s = str(s)                                             # prefer a decimal (avoid the "10" in "log10")
    m = re.search(r'-?\d*\.\d+', s) or re.search(r'-?\d+', s)
    return float(m.group()) if m else float('nan')
def _se_type(s):
    s = str(s).lower()
    for k in ('fixed', 'log10', 'jackknife', 'gn'):
        if k in s: return {'gn': 'Gauss-Newton'}.get(k, k)
    return ''
def _ci(s):
    s = str(s)
    if s.strip() in ('-', '—', '', 'nan'): return (np.nan, np.nan)
    parts = re.split(r'\s[–—-]\s', s)                 # split on spaced dash
    if len(parts) != 2: return (np.nan, np.nan)
    return _first_float(parts[0]), _first_float(parts[1].replace('mM', '').replace('nm', ''))

SECT = {'Multivalent (Fig. 4)': 'Multivalent model MV3 (Fig. 4)',
        'Two-state A2 (SI)':     'Two-state model M3 (row A2)'}
s4_rows, s4_check = [], []
for _, r in S4.iterrows():
    lo, hi = _ci(r['CI'])
    se = _first_float(r['SE_a']); se = np.nan if str(r['SE_a']).strip().lower() == 'fixed' else se
    s4_rows.append(dict(Section=SECT.get(r['model'], r['model']), Parameter=r['parameter'],
                        Units=r['units'] if str(r['units']) != 'nan' else '',
                        Estimate=round(float(r['estimate']), 4),
                        SE=('' if math.isnan(se) else round(se, 5)),
                        SE_type=_se_type(r['SE_a']),
                        CI95_low=('' if math.isnan(lo) else lo),
                        CI95_high=('' if math.isnan(hi) else hi),
                        Meaning=r['meaning']))
    s4_check.append((r['model'], r['parameter'], float(r['estimate']), se, lo, hi))
pd.DataFrame(s4_rows).to_csv(os.path.join(OUTDIR, 'TableS4_FittedParameters.csv'), index=False)

# Published S4 values: (model, parameter) -> (estimate, est_dp, se, se_dp, ci_lo, ci_hi, ci_dp)
REF_S4 = {
    ('Multivalent (Fig. 4)','C_ref'):   (26.7,1, 0.042,3, 22.0,31.3,1),
    ('Multivalent (Fig. 4)','λ'):  (12.8,1, 0.69,2,  11.8,14.3,1),
    ('Multivalent (Fig. 4)','κ'):  (0.079,3,0.0085,4,0.065,0.092,3),
    ('Multivalent (Fig. 4)','w'):       (0.57,2, None,0,  None,None,0),
    ('Two-state A2 (SI)','λ'):     (90.1,1, 4.3,1,   82.6,99.1,1),
    ('Two-state A2 (SI)','ΔS_Bond'):(-8.38,2,0.12,2, -8.63,-8.14,2),
    ('Two-state A2 (SI)','ΔS_0'):  (7.31,2, 0.84,2,  5.69,8.87,2),
}
for model, param, est, se, lo, hi in s4_check:
    if (model, param) not in REF_S4: continue
    r_est, de, r_se, dse, r_lo, r_hi, dci = REF_S4[(model, param)]
    tag = f'S4 {param:<9s} ({"MV3" if "Multi" in model else "A2"})'
    checks.append((f'{tag} value', est, r_est, de))
    if r_se is not None and not math.isnan(se): checks.append((f'{tag} s.e. ', se, r_se, dse))
    if r_lo is not None and not math.isnan(lo): checks.append((f'{tag} CI_lo', lo, r_lo, dci))
    if r_hi is not None and not math.isnan(hi): checks.append((f'{tag} CI_hi', hi, r_hi, dci))

# ----------------------------------------------------------------------------- cross-source check
# MV3 (row M3 of Table S1) must equal the "Multivalent model (Fig. 4)" row of the ML comparison.
mv3 = ts1('M3'); mvx = ml('Multivalent model (Fig. 4)')
checks.append(('X  MV3 RMSD == ML-table multivalent', mv3[0], mvx[0], 4))
checks.append(('X  MV3 rho  == ML-table multivalent', mv3[1], mvx[1], 4))

# ----------------------------------------------------------------------------- report
print(f"SI tables regenerated from the pipeline and checked against PNAS_SI_Appendix_clean.docx")
print(f"  {os.path.relpath(os.path.join(OUTDIR,'TableS1_ModelPerformance.csv'), REPO)}  "
      f"({len(s1_rows)} models)")
print(f"  {os.path.relpath(os.path.join(OUTDIR,'TableS4_FittedParameters.csv'), REPO)}  "
      f"({len(s4_rows)} parameters)")
print("-" * 72)
n_fail = 0
for label, got, ref, dp in checks:
    tol = 0.5 * 10 ** (-dp) + 1e-6
    ok = abs(got - ref) <= tol
    if not ok:
        n_fail += 1
        print(f"  MISMATCH  {label}: computed {got:.4f}  vs doc {ref:.{dp}f}  (|d|={abs(got-ref):.4f} > {tol:.4f})")
print("-" * 72)
if n_fail == 0:
    print(f"All {len(checks)} checked values match the clean SI Appendix (within displayed precision).")
else:
    print(f"{n_fail} of {len(checks)} values DIFFER from the clean SI Appendix (see MISMATCH lines above).")
    if STRICT:
        sys.exit(1)
