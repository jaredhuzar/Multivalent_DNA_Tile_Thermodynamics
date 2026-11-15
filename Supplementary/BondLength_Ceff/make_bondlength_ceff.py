"""
SI figure (replaces S11_BondLength_Entropy, which plotted the noise-biased two-state van't Hoff entropy by bond
length): does the effective concentration -- i.e. the per-bond tethering penalty RT ln(1 M / C_eff) -- depend on
bond length?

Model: the Figure 4 multivalent model with one extra factor per bond type,
    C_eff,i = C_ref * phi_type(i) * exp(-d_i/lambda) * exp(-kappa (N - 4)),   phi_7nt = 1 (7-nt and A/B bonds)
fitted (log10 phi_5, log10 phi_6 free) with and without the valency term, intermediate readout, canonical 236
interfaces, 10-fold CV with the Table S1 folds.  Compared with a worm-like-chain estimate of the tether effective
concentration (Table S1 model B4, C_WLC = exp(dS_teth / R)), which predicts phi_5 = C_WLC(5)/C_WLC(7) ~ 10.
Caveat: every interface carries a 5-nt reporter site, and 5-nt bonds dominate the 9-10-site interfaces, so phi_5
and kappa are partly confounded (their correlation is reported).
Writes figure.svg, figure.png, fits.csv here.
"""
import os, sys, numpy as np, pandas as pd
from scipy.optimize import least_squares
from scipy.stats import spearmanr, t as tdist
from sklearn.model_selection import KFold
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(REPO, 'Figure4_MultivalentModel')); import mv_model as mv

Bl = mv.UM['Bl']; A = mv.arrays(Bl); YK = mv.Y_C + 273.15; n = len(YK)
TYPE = np.zeros_like(A['H'], int)                          # 0: 7-nt / A / B, 1: 5-nt, 2: 6-nt
for i, b in enumerate(Bl):
    for k, x in enumerate(b): TYPE[i, k] = {'5': 1, '6': 2}.get(x[4], 0)
A['T'] = TYPE

def theta(T, p, A, w=mv.W_INT, kappa_on=True):
    lc, lam, l5, l6 = p[:4]; kap = p[4] if kappa_on else 0.0
    phi = np.choose(A['T'], [1.0, 10 ** l5, 10 ** l6])
    C = 10 ** lc * phi * np.exp(-A['D'] / lam) * np.exp(-kap * (A['N'] - 4))[:, None]
    x = np.where(A['mask'], np.exp(np.clip(-(A['H'] - T[:, None] * A['S']) / (mv.R * T[:, None]), -700, 700)) * C, 0.0)
    a = np.expm1(np.sum(np.log1p(x), 1)) / 10 ** lc * mv.C_TILE
    f = np.where(a > 1e-12, ((2 * a + 1) - np.sqrt(4 * a + 1)) / (2 * np.maximum(a, 1e-300)), a)
    prep = np.where(A['hasrep'], np.sum(np.where(A['rep'], x / (1 + x), 0), 1) / A['nrep'], 1.0)
    return f * (w + (1 - w) * prep)

def tm(p, A, kappa_on):
    m = len(A['N']); lo = np.full(m, 230.0); hi = np.full(m, 420.0)
    for _ in range(48):
        mid = (lo + hi) / 2; up = theta(mid, p, A, kappa_on=kappa_on) > .5
        lo = np.where(up, mid, lo); hi = np.where(up, hi, mid)
    return (lo + hi) / 2

def fit(idx, kappa_on):
    As = mv._sub(A, idx)
    p0 = [-1.6, 13.0, 0.0, 0.0] + ([0.08] if kappa_on else [])
    lb = [-8, 1, -3, -3] + ([-2] if kappa_on else []); ub = [1, 1e4, 3, 3] + ([2] if kappa_on else [])
    return least_squares(lambda p: tm(p, As, kappa_on) - YK[idx], p0, bounds=(lb, ub), diff_step=1e-4)

FOLDS = list(KFold(10, shuffle=True, random_state=0).split(np.arange(n)))
rows = []
for kappa_on in (True, False):
    r = fit(np.arange(n), kappa_on); oof = np.full(n, np.nan)
    for tr, te in FOLDS: oof[te] = tm(fit(tr, kappa_on).x, mv._sub(A, te), kappa_on)
    J = r.jac; dof = n - J.shape[1]; cov = (r.fun @ r.fun) / dof * np.linalg.inv(J.T @ J); se = np.sqrt(np.diag(cov))
    tc = tdist.ppf(.975, dof); corr = cov / np.outer(se, se)
    rows.append(dict(model='+ bond-type C_eff' + (' (with kappa)' if kappa_on else ' (no kappa)'), n_params=len(r.x),
                     cvRMS=np.sqrt(np.mean((oof - YK) ** 2)), cvrho=spearmanr(YK, oof)[0],
                     log10_Cref_M=r.x[0], lambda_nm=r.x[1], kappa_per_site=r.x[4] if kappa_on else np.nan,
                     phi5=10 ** r.x[2], phi5_lo=10 ** (r.x[2] - tc * se[2]), phi5_hi=10 ** (r.x[2] + tc * se[2]),
                     phi6=10 ** r.x[3], phi6_lo=10 ** (r.x[3] - tc * se[3]), phi6_hi=10 ** (r.x[3] + tc * se[3]),
                     corr_phi5_kappa=corr[2, 4] if kappa_on else np.nan))
    print(rows[-1], flush=True)
base = pd.read_csv(os.path.join(REPO, 'Figure4_MultivalentModel/model_fit/cv_comparison.csv'))
b3 = base[(base.readout == 'intermediate') & (base.n_params == 3)].iloc[0]
R_CAL = 1.98720425864; dS_teth = mv.UM['dS_teth']
cw = {L: np.exp(dS_teth(L) / R_CAL) for L in ('5', '6', '7')}
F = pd.DataFrame(rows); F['WLC_phi5'] = cw['5'] / cw['7']; F['WLC_phi6'] = cw['6'] / cw['7']
F['WLC_Ceff_5nt_mM'] = cw['5'] * 1e3; F['WLC_Ceff_6nt_mM'] = cw['6'] * 1e3; F['WLC_Ceff_7nt_mM'] = cw['7'] * 1e3
F['base_cvRMS'] = b3.cvRMS
F.round(5).to_csv(os.path.join(HERE, 'fits.csv'), index=False)

# ---- figure ----
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e6e5e0'
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'svg.fonttype': 'none'})
fig, ax = plt.subplots(1, 2, figsize=(10.2, 4.0), gridspec_kw={'width_ratios': [1.25, 1]})
a = ax[0]; xs = np.arange(3); lab = ['7-nt (reference)', '5-nt', '6-nt']
a.axhline(1, color=GRID, lw=1, zorder=0)
for k, (r, mk, col, dx) in enumerate(zip(rows, ('o', 's'), ('#1baf7a', '#8a8984'), (-.12, .12))):
    v = [1, r['phi5'], r['phi6']]; lo = [1, r['phi5_lo'], r['phi6_lo']]; hi = [1, r['phi5_hi'], r['phi6_hi']]
    a.errorbar(xs + dx, v, yerr=[np.array(v) - lo, np.array(hi) - v], fmt=mk, color=col, ms=7, capsize=3, lw=1.3,
               label='fitted, ' + ('with valency term' if 'with' in r['model'] else 'no valency term'))
a.scatter(xs, [1, cw['5'] / cw['7'], cw['6'] / cw['7']], marker='_', s=500, color='#e34948', lw=2.5, zorder=4,
          label='worm-like-chain tether estimate')
a.set_yscale('log'); a.set_xticks(xs); a.set_xticklabels(lab)
a.set(ylabel='effective concentration relative to 7-nt bonds', ylim=(.05, 40))
a.legend(fontsize=7.5, frameon=False, loc='upper left'); a.grid(True, axis='y', color=GRID, lw=.6)
a.set_title('A   bond length barely changes the effective concentration', loc='left', fontsize=9.5, weight='bold')
a = ax[1]; names = ['Fig. 4 model\n(one $C_{ref}$)', '+ bond-type\n(with κ)', '+ bond-type\n(no κ)']
vals = [b3.cvRMS, rows[0]['cvRMS'], rows[1]['cvRMS']]
a.bar(range(3), vals, color=['#1baf7a', '#7fd0b0', '#bdbdbd'], width=.6, zorder=2)
for i, v in enumerate(vals): a.text(i, v + .03, f'{v:.2f} °C', ha='center', fontsize=8)
a.set_xticks(range(3)); a.set_xticklabels(names, fontsize=8); a.set(ylabel='10-fold CV RMS error in $T_m$ (°C)', ylim=(0, max(vals) * 1.25))
a.grid(True, axis='y', color=GRID, lw=.6, zorder=0)
a.set_title('B   cross-validated error', loc='left', fontsize=9.5, weight='bold')
fig.tight_layout(); fig.savefig(os.path.join(HERE, 'figure.svg')); fig.savefig(os.path.join(HERE, 'figure.png'), dpi=160)
print('WLC C_eff (mM):', {k: round(v * 1e3, 3) for k, v in cw.items()})
