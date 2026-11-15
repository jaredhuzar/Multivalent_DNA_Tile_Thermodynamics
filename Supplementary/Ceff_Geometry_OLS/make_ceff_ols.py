"""
SI figure (replaces the avidity-alpha OLS, Alpha_Geometry_Controlled / S2_OLS_alpha_geometry): geometry-controlled
regression of each interface's effective concentration.

Response: log10 C_eff (M) of every interface, fitted on its own (one uniform C_eff per construct, no distance or
valency terms) to its measured Tm with the multivalent model of Figure 4 (intermediate readout, w = 0.57).
alpha is not used: it is referenced to the 1 M standard state, so alpha < 1 and its fall with N are built in.
Predictors: number of occupied sites N, bond span (nm), mean bond length (nt, over occupied sites), GC fraction
(mean over the interface's duplex sequences).  OLS on standardized predictors (beta, SE, t, p, VIF, R^2),
partial Spearman on residuals, added-variable panels.  Canonical 236-interface set (Table S1 / Figure 4).
Writes figure.svg, figure.png, data.csv, ols_coefficients.csv, partial_correlations.csv here.
"""
import os, sys, numpy as np, pandas as pd
from scipy.stats import spearmanr, t as tdist
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(REPO, 'Figure4_MultivalentModel')); import mv_model as mv

A = mv.arrays(mv.UM['Bl']); lc = mv.fit_ceff_per_construct(A, mv.Y_C + 273.15)
SEQ = mv.UM['SEQ']; LEN = {'5': 5, '6': 6, '7': 7, '9': 9, 'A': 7, 'B': 7}
def gc(code):
    s = [SEQ[ch][2 * i + k] for i, ch in enumerate(code) if ch != '0' for k in (1, 2)]
    return np.mean([sum(b in 'GC' for b in q) / len(q) for q in s])
rows = []
for code, l, tm in zip(mv.CODES, lc, mv.Y_C):
    st = [i for i, ch in enumerate(code) if ch != '0']
    rows.append(dict(code=code, group=mv.group(code), Tm_C=tm, log10_Ceff_M=l, n_sites=len(st),
                     span_nm=(st[-1] - st[0]) * 6.0, bond_len=np.mean([LEN[code[i]] for i in st]), GC=gc(code)))
d = pd.DataFrame(rows)
P = ['n_sites', 'span_nm', 'bond_len', 'GC']; LAB = {'n_sites': 'number of sites', 'span_nm': 'bond span (nm)',
                                                     'bond_len': 'mean bond length (nt)', 'GC': 'GC fraction'}
z = lambda v: (v - v.mean()) / v.std()
X = np.column_stack([np.ones(len(d))] + [z(d[c]).values for c in P]); y = d.log10_Ceff_M.values
b, *_ = np.linalg.lstsq(X, y, rcond=None); dof = len(d) - X.shape[1]
s2 = np.sum((y - X @ b) ** 2) / dof; se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X))); tv = b / se
pv = 2 * (1 - tdist.cdf(np.abs(tv), dof)); R2 = 1 - np.sum((y - X @ b) ** 2) / np.sum((y - y.mean()) ** 2)
vif = {}
for c in P:
    o = [k for k in P if k != c]; Xo = np.column_stack([np.ones(len(d))] + [d[k].values for k in o])
    bo, *_ = np.linalg.lstsq(Xo, d[c].values, rcond=None)
    vif[c] = 1 / (np.sum((d[c].values - Xo @ bo) ** 2) / np.sum((d[c].values - d[c].mean()) ** 2))
pd.DataFrame(dict(predictor=['intercept'] + P, beta_std=b, SE=se, t=tv, p=pv,
                  VIF=[np.nan] + [vif[c] for c in P])).assign(R2=R2, n=len(d)).round(5) \
  .to_csv(os.path.join(HERE, 'ols_coefficients.csv'), index=False)

INK2, GRID = '#52514e', '#e6e5e0'
GC_ = {'7-nt': '#3B6AA0', '5-nt': '#E8A33D', '6-nt': '#4DA167', 'A/B': '#8a6fbf'}
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'svg.fonttype': 'none'})
fig, ax = plt.subplots(1, 4, figsize=(15, 3.9)); pc = []
for a, c in zip(ax, P):
    ctrl = [k for k in P if k != c]; Xc = np.column_stack([np.ones(len(d))] + [d[k].values for k in ctrl])
    res = lambda v: v - Xc @ np.linalg.lstsq(Xc, v, rcond=None)[0]
    rx, ry = res(d[c].values.astype(float)), res(y); r_, p_ = spearmanr(rx, ry); mr = spearmanr(d[c], y)[0]
    d[f'xres_{c}'] = rx; d[f'yres_vs_{c}'] = ry
    pc.append(dict(predictor=c, marginal_rho=mr, partial_rho=r_, partial_p=p_, n=len(d)))
    for gname, col in GC_.items():
        m = (d.group == gname).values
        a.scatter(rx[m], ry[m], s=13, color=col, alpha=.6, edgecolor='none', label=gname)
    bb = np.polyfit(rx, ry, 1); xs = np.linspace(rx.min(), rx.max(), 30); a.plot(xs, np.polyval(bb, xs), color='#0b0b0b', lw=1.5)
    a.set(xlabel=f'{LAB[c]} (residual)', ylabel='log$_{10}$ $C_{eff}$ (residual)')
    a.set_title(f'{LAB[c]}\npartial ρ = {r_:+.2f}, P = {p_:.1g}', fontsize=9); a.grid(True, color=GRID, lw=.6)
ax[0].legend(fontsize=7, frameon=False, loc='lower left')
fig.suptitle(f'Effective concentration of each interface vs geometry, each factor controlling the other three '
             f'(OLS R² = {R2:.2f}, n = {len(d)})', fontsize=10, y=1.02)
fig.tight_layout(); fig.savefig(os.path.join(HERE, 'figure.svg'), bbox_inches='tight'); fig.savefig(os.path.join(HERE, 'figure.png'), dpi=160, bbox_inches='tight')
pd.DataFrame(pc).round(5).to_csv(os.path.join(HERE, 'partial_correlations.csv'), index=False)
d.round(5).to_csv(os.path.join(HERE, 'data.csv'), index=False)
print(f'C_eff OLS: n={len(d)}, R2={R2:.3f}, median C_eff {10 ** np.median(y) * 1e3:.2f} mM; ' +
      '; '.join(f"{r['predictor']}: partial rho {r['partial_rho']:+.2f} (P={r['partial_p']:.1g})" for r in pc))
