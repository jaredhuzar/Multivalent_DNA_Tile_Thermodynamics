#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Regenerate Dataset S1 / Table S3 -- the per-interface summary of the analysis set.

One row per interface of the 236-interface set, built by aggregating the authoritative
per-well fits in Filtered_DF_Tm.csv (Tm, transition width, plate, confidence, predicted
quantities) and the robust van't Hoff fits in Filtered_DF_HS.csv (dH, dS, dG37, exactly as
Experimental_Trends/construct_medians.csv).  Geometry (sites, gaps, span, distances) is a
function of the 11-character construct code at the unified pitch of 6.0 nm per site.

Columns (order follows the SI Appendix "Table S3 (separate file)" description):
  code                            construct code
  N_bonds                         bond count N
  bond_length_nt                  bond length (nt; mean over sites if mixed)
  occupied_sites_0indexed         0-indexed occupied sites (';'-joined)
  mean_gap_nm                     mean gap g-bar
  bond_span_nm                    bond span g*
  mean_nearest_bond_distance_nm   mean over bonds of the distance to the nearest other bond
  mean_pairwise_distance_nm       mean distance over all bond pairs
  dG_seq_additiveNN_kcal          additive nearest-neighbor dG_seq
  Tm_predicted_from_seq_C         Tm predicted from dG_seq
  n_wells_QC                      number of quality-control-passing wells
  Tm_median_C                     median Tm over those wells
  Tm_sd_C                         standard deviation of Tm over those wells (ddof=1)
  transition_width_median_C       median transition width (Boltzmann dx)
  dH_vantHoff_kcal                median apparent van't Hoff dH     (robust HS set; NaN if absent)
  dS_vantHoff_cal_per_K           median apparent van't Hoff dS     (robust HS set; NaN if absent)
  dG37_vantHoff_kcal              median apparent van't Hoff dG37   (robust HS set; NaN if absent)
  confidence_class                interface confidence class (mode over wells; ties -> 'standard')
  plate_of_origin                 plate(s) of origin (';'-joined if more than one)

Output (-> <repo>/Pseudo_Figures/SI/, the SI deliverables folder):
  Dataset_S1_InterfaceSummary.csv

Consistency check ("same numbers as the pipeline"): the shared columns
(Tm, N_bonds, bond_length, mean_gap, bond_span, dH, dS, dG37) are compared row-by-row with
Experimental_Trends/construct_medians.csv; mean_pairwise_distance is cross-checked against the
Distances_nm column carried in Filtered_DF_Tm.  A PASS/FAIL report is printed.  Exit 0
(report only) unless --strict (exit 1 on any mismatch).

Run from anywhere:  python3 Supplementary/SI_Tables/make_dataset_s1.py [--strict]
"""
import os, re, sys
import numpy as np
import pandas as pd

def _find_repo():
    d = os.path.dirname(os.path.abspath(__file__))
    while d != '/' and not os.path.exists(os.path.join(d, 'Filtered_DF_Tm.csv')):
        d = os.path.dirname(d)
    return d
REPO = _find_repo()
OUTDIR = os.path.join(REPO, 'Pseudo_Figures', 'SI')
os.makedirs(OUTDIR, exist_ok=True)
STRICT = '--strict' in sys.argv[1:]
PITCH_NM = 6.0                                  # unified nm per site (Fig 3a/3d, OLS, Fig 4D)
LEN = {'5': 5, 'X': 5, '6': 6, '7': 7, '9': 9, 'A': 7, 'B': 7}   # nt per bond by code char

def _req(path):
    p = os.path.join(REPO, path)
    if not os.path.exists(p):
        sys.exit(f"ERROR: required source not found: {path} (run run_all.sh first)")
    return p

# ----------------------------------------------------------------------------- geometry from code
def geom(code):
    pos = [i for i, ch in enumerate(code) if ch != '0']
    n = len(pos)
    bl = float(np.mean([LEN.get(ch, np.nan) for ch in code if ch != '0'])) if n else np.nan
    if n > 1:
        sp = sorted(pos)
        mean_gap = float(np.mean(np.diff(sp))) * PITCH_NM
        span = (sp[-1] - sp[0]) * PITCH_NM
        nearest = float(np.mean([min(abs(p - q) for q in pos if q != p) for p in pos])) * PITCH_NM
        pairs = [abs(a - b) for i, a in enumerate(pos) for b in pos[i + 1:]]
        pairwise = float(np.mean(pairs)) * PITCH_NM
    else:
        mean_gap = np.nan; span = 0.0; nearest = np.nan; pairwise = np.nan
    return dict(N_bonds=n, bond_length_nt=bl, occupied_sites_0indexed=';'.join(map(str, pos)),
                mean_gap_nm=mean_gap, bond_span_nm=span,
                mean_nearest_bond_distance_nm=nearest, mean_pairwise_distance_nm=pairwise)

# ----------------------------------------------------------------------------- per-well aggregation
tm = pd.read_csv(_req('Filtered_DF_Tm.csv'))
tm['code'] = tm['Code'].astype(str).str.partition('.')[0]

def _conf(s):                                    # interface class = mode over wells, tie -> 'standard'
    vc = s.value_counts()
    top = list(vc[vc == vc.max()].index)
    return 'standard' if 'standard' in top else top[0]

def _plate(s):                                   # cleaned plate stem(s), joined if >1
    labs = sorted({re.sub(r'(_Formatted)?\.(csv|txt)$', '', str(x)) for x in s})
    return ';'.join(labs)

agg = tm.groupby('code').agg(
    n_wells_QC=('n_reps', 'first'),
    Tm_median_C=('Tm', 'median'),
    Tm_sd_C=('Tm', 'std'),
    transition_width_median_C=('dx', 'median'),
    dG_seq_additiveNN_kcal=('Predicted_deltaG', 'first'),
    Tm_predicted_from_seq_C=('Predicted_tm', 'first'),
    pairwise_nm_src=('Distances_nm', 'first'),           # for cross-check only
).reset_index()
agg['confidence_class'] = tm.groupby('code')['Confidence'].apply(_conf).values
agg['plate_of_origin'] = tm.groupby('code')['Experiment'].apply(_plate).values

# ----------------------------------------------------------------------------- robust van't Hoff (HS)
hs = pd.read_csv(_req('Filtered_DF_HS.csv'))
hs['code'] = hs['Code'].astype(str).str.partition('.')[0]
hs = hs[(hs['Reliable'] == True) & (hs['Confidence'] == 'high')]
hsG = hs.groupby('code').agg(dG37_vantHoff_kcal=('Exp_DeltaG', 'median'),
                             dS_vantHoff_cal_per_K=('Entropy', 'median'),
                             dH_vantHoff_kcal=('Enthalpy', 'median')).reset_index()

# ----------------------------------------------------------------------------- assemble
G = pd.DataFrame([dict(code=c, **geom(c)) for c in agg['code']])
df = (G.merge(agg, on='code', how='left').merge(hsG, on='code', how='left'))
COLS = ['code', 'N_bonds', 'bond_length_nt', 'occupied_sites_0indexed', 'mean_gap_nm',
        'bond_span_nm', 'mean_nearest_bond_distance_nm', 'mean_pairwise_distance_nm',
        'dG_seq_additiveNN_kcal', 'Tm_predicted_from_seq_C', 'n_wells_QC', 'Tm_median_C',
        'Tm_sd_C', 'transition_width_median_C', 'dH_vantHoff_kcal', 'dS_vantHoff_cal_per_K',
        'dG37_vantHoff_kcal', 'confidence_class', 'plate_of_origin']
out = df.sort_values('code').reset_index(drop=True)
rounded = out.copy()
for c in COLS:
    if c not in ('code', 'occupied_sites_0indexed', 'confidence_class', 'plate_of_origin', 'N_bonds', 'n_wells_QC'):
        rounded[c] = rounded[c].astype(float).round(4)
rounded[COLS].to_csv(os.path.join(OUTDIR, 'Dataset_S1_InterfaceSummary.csv'), index=False)

# ----------------------------------------------------------------------------- consistency check
cm = pd.read_csv(_req('Experimental_Trends/construct_medians.csv'))
cm['code'] = cm['code'].astype(str)
chk = out.merge(cm, on='code', how='outer', suffixes=('', '_cm'))
pairs = [('Tm_median_C', 'Tm'), ('N_bonds', 'Num_Bonds'), ('bond_length_nt', 'bond_len'),
         ('mean_gap_nm', 'MeanGap_nm'), ('bond_span_nm', 'span_nm'),
         ('dH_vantHoff_kcal', 'Enthalpy'), ('dS_vantHoff_cal_per_K', 'Entropy'),
         ('dG37_vantHoff_kcal', 'dG37')]
fails = []
if len(chk) != 236 or out['code'].nunique() != 236:
    fails.append(f"row count: Dataset S1 has {out['code'].nunique()} interfaces (expected 236)")
for mine, theirs in pairs:
    a = pd.to_numeric(chk[mine], errors='coerce').values
    b = pd.to_numeric(chk[theirs], errors='coerce').values
    bad = ~(np.isclose(a, b, atol=1e-4, equal_nan=True))
    if bad.any():
        i = np.where(bad)[0][0]
        fails.append(f"{mine} != construct_medians.{theirs} "
                     f"(e.g. {chk['code'].values[i]}: {a[i]} vs {b[i]}; {int(bad.sum())} rows)")
# cross-check mean pairwise distance against the Distances_nm carried in Filtered_DF_Tm
multi = out['N_bonds'] > 1
bad = ~np.isclose(out.loc[multi, 'mean_pairwise_distance_nm'], out.loc[multi, 'pairwise_nm_src'], atol=1e-3)
if bad.any():
    fails.append(f"mean_pairwise_distance_nm != Filtered_DF_Tm.Distances_nm ({int(bad.sum())} rows)")

# ----------------------------------------------------------------------------- report
rel = os.path.relpath(os.path.join(OUTDIR, 'Dataset_S1_InterfaceSummary.csv'), REPO)
print(f"Dataset S1 / Table S3 regenerated: {rel}")
print(f"  {out['code'].nunique()} interfaces x {len(COLS)} columns  "
      f"({int((out['dG37_vantHoff_kcal'].isna()).sum())} without robust van't Hoff -> NaN, "
      f"matching construct_medians)")
print("-" * 72)
if not fails:
    print(f"Consistency: all shared columns match construct_medians.csv and Filtered_DF_Tm (236/236).")
else:
    print(f"{len(fails)} consistency problem(s):")
    for f in fails:
        print("  MISMATCH ", f)
    if STRICT:
        sys.exit(1)
