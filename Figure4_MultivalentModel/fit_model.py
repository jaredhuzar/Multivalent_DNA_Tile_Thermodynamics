"""
Fit the multivalent model (mv_model.py) to the 236 canonical interfaces and evaluate it exactly like Table S1:
10-fold cross-validation over interfaces, KFold(10, shuffle=True, random_state=0).

Nested variants x readouts:   1 uniform C_eff  |  2 + distance (lambda)  |  3 + valency (kappa)
                               intermediate readout (w = 0.57, main)  and  docking (w = 1)
Main model = 3 parameters, intermediate readout.  Its 95% CIs: Gauss-Newton (t, df = n - 3) and a construct
bootstrap (--boot N, default 200; --boot 0 skips it).

Writes model_fit/{cv_comparison.csv, params.csv, oof_predictions.csv, bootstrap_params.csv}.
Run from anywhere (~5 min, plus ~3 s per bootstrap replicate).
"""
import os, sys, numpy as np, pandas as pd
from scipy.stats import spearmanr, t as tdist
from sklearn.model_selection import KFold
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mv_model as mv

OUT = os.path.join(mv.HERE, 'model_fit'); os.makedirs(OUT, exist_ok=True)
NBOOT = int(sys.argv[sys.argv.index('--boot') + 1]) if '--boot' in sys.argv else 200
A = mv.arrays([mv.UM['Bl'][i] for i in range(mv.N_ALL)]); YK = mv.Y_C + 273.15; n = len(YK)
FOLDS = list(KFold(10, shuffle=True, random_state=0).split(np.arange(n)))
NAMES = {1: 'one uniform C_eff', 2: '+ distance (lambda)', 3: '+ valency (kappa)'}

rows, preds = [], {}
for rname, w in (('intermediate', mv.W_INT), ('docking', 1.0)):
    for npar in (1, 2, 3):
        r = mv.fit(A, YK, npar, w); fp = mv.tm(r.x, A, w) - 273.15
        oof = np.full(n, np.nan)
        for tr, te in FOLDS:
            q = mv.fit(A, YK, npar, w, idx=tr).x; oof[te] = mv.tm(q, mv._sub(A, te), w) - 273.15
        rows.append(dict(model=f'multivalent {NAMES[npar]}', readout=rname, w=w, n_params=npar,
                         fitRMS=np.sqrt(np.mean((fp - mv.Y_C) ** 2)), cvRMS=np.sqrt(np.mean((oof - mv.Y_C) ** 2)),
                         cvrho=spearmanr(mv.Y_C, oof)[0], log10_Cref_M=r.x[0],
                         lambda_nm=r.x[1] if npar > 1 else np.nan, kappa_per_site=r.x[2] if npar > 2 else np.nan))
        preds[(rname, npar)] = (fp, oof, r)
        print(rows[-1], flush=True)
T1 = pd.read_csv(os.path.join(mv.REPO, 'Supplementary/TableS1_Model_Comparison_All/model_comparison_all_235.csv'))
for key, lab in (('A2', 'phenomenological A2 (SI; former Fig 4)'), ('B1', 'phenomenological B1, no distance weighting (SI)')):
    r = T1[T1.model.str.startswith(key)].iloc[0]
    rows.append(dict(model=lab, readout='two-state', w=np.nan, n_params=3 if key == 'A2' else 2,
                     fitRMS=r.fitRMS, cvRMS=r.cvRMS, cvrho=r.cvrho))
CV = pd.DataFrame(rows); CV.round(4).to_csv(os.path.join(OUT, 'cv_comparison.csv'), index=False)

# ---- main model: parameters + CIs ----
fp, oof, r = preds[('intermediate', 3)]
J = r.jac; dof = n - J.shape[1]; cov = (r.fun @ r.fun) / dof * np.linalg.inv(J.T @ J); se = np.sqrt(np.diag(cov))
tc = tdist.ppf(.975, dof)
PAR = ['log10_Cref_M', 'lambda_nm', 'kappa_per_site']
P = pd.DataFrame(dict(parameter=PAR, point_estimate=r.x, GN_SE=se, GN_CI95_low=r.x - tc * se, GN_CI95_high=r.x + tc * se))
P['n_constructs'] = n; P['readout_w'] = mv.W_INT
P.round(5).to_csv(os.path.join(OUT, 'params.csv'), index=False)

O = pd.DataFrame(dict(code=mv.CODES, group=[mv.group(c) for c in mv.CODES], n_sites=[mv.n_sites(c) for c in mv.CODES],
                      Tm_meas_C=mv.Y_C, Tm_pred_oof_C=oof, Tm_pred_fit_C=fp,
                      width_pred_C=mv.width(r.x, A)))
for (rname, npar), (fp_, oof_, _) in preds.items():
    if (rname, npar) != ('intermediate', 3): O[f'oof_{rname}_{npar}par'] = oof_
O.round(4).to_csv(os.path.join(OUT, 'oof_predictions.csv'), index=False)
print('main model:', dict(zip(PAR, np.round(r.x, 4))), 'GN SE', np.round(se, 4))

if NBOOT > 0:
    rng = np.random.default_rng(1); B = []
    for b in range(NBOOT):
        idx = rng.integers(0, n, n); B.append(mv.fit(A, YK, 3, mv.W_INT, idx=idx, p0=list(r.x)).x)
        if (b + 1) % 20 == 0: print(f'  bootstrap {b + 1}/{NBOOT}', flush=True)
    B = np.array(B); pd.DataFrame(B, columns=PAR).to_csv(os.path.join(OUT, 'bootstrap_params.csv'), index=False)
    P['boot_SE'] = B.std(0, ddof=1); P['boot_CI95_low'] = np.percentile(B, 2.5, 0); P['boot_CI95_high'] = np.percentile(B, 97.5, 0)
    P['n_bootstrap'] = NBOOT
    P.round(5).to_csv(os.path.join(OUT, 'params.csv'), index=False)
    print(P.round(4).to_string())
