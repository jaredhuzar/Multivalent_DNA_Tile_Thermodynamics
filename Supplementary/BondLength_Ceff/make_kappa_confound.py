"""
SI Fig. S11: the valency term kappa cannot be separated from a bond-length effect in this library.

  (composition table kept in kappa_confound_composition.csv; not plotted)
  A  10-fold CV RMS error (Table S1 folds) of the Figure 4 model with kappa, with a separate 5-nt (and 6-nt)
     effective concentration instead of kappa, and with both          <- fits.csv + Figure 4 cv_comparison.csv
  B  profile along kappa: for each fixed kappa, the other parameters (C_ref, lambda, phi5, phi6) are refitted to all
     236 interfaces; the fitted 5-nt factor phi5 and the RMS error are shown -- the error stays flat while phi5
     trades off against kappa (the two are interchangeable)
Model as in make_bondlength_ceff.py: C_eff,i = C_ref phi_type(i) e^{-d_i/lambda} e^{-kappa (N-4)}, intermediate readout.
Writes kappa_confound.svg and kappa_confound_{composition,cv,profile}.csv here.
"""
import os, sys, numpy as np, pandas as pd
from scipy.optimize import least_squares
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(REPO, 'Figure4_MultivalentModel')); import mv_model as mv

Bl = mv.UM['Bl']; A = mv.arrays(Bl); YK = mv.Y_C + 273.15; n = len(YK)
TYPE = np.zeros_like(A['H'], int)
for i, b in enumerate(Bl):
    for k, x in enumerate(b): TYPE[i, k] = {'5': 1, '6': 2}.get(x[4], 0)

def tm(lc, lam, l5, l6, kap):
    phi = np.choose(TYPE, [1.0, 10 ** l5, 10 ** l6])
    C = 10 ** lc * phi * np.exp(-A['D'] / lam) * np.exp(-kap * (A['N'] - 4))[:, None]
    lo = np.full(n, 230.0); hi = np.full(n, 420.0)
    for _ in range(48):
        T = (lo + hi) / 2
        x = np.where(A['mask'], np.exp(np.clip(-(A['H'] - T[:, None] * A['S']) / (mv.R * T[:, None]), -700, 700)) * C, 0.0)
        a = np.expm1(np.sum(np.log1p(x), 1)) / 10 ** lc * mv.C_TILE
        f = np.where(a > 1e-12, ((2 * a + 1) - np.sqrt(4 * a + 1)) / (2 * np.maximum(a, 1e-300)), a)
        prep = np.where(A['hasrep'], np.sum(np.where(A['rep'], x / (1 + x), 0), 1) / A['nrep'], 1.0)
        up = f * (mv.W_INT + (1 - mv.W_INT) * prep) > .5
        lo = np.where(up, T, lo); hi = np.where(up, hi, T)
    return (lo + hi) / 2

# ---- C: profile along kappa ----
prof, p0 = [], [-1.5, 12.7, -0.2, -0.1]
for kap in np.round(np.arange(-0.02, 0.141, 0.01), 3):
    r = least_squares(lambda q: tm(q[0], q[1], q[2], q[3], kap) - YK, p0, bounds=([-8, 1, -3, -3], [1, 1e4, 3, 3]), diff_step=1e-4)
    p0 = list(r.x)
    prof.append(dict(kappa_per_site=kap, phi5=10 ** r.x[2], phi6=10 ** r.x[3], lambda_nm=r.x[1], log10_Cref_M=r.x[0],
                     fit_RMS_C=float(np.sqrt(np.mean(r.fun ** 2)))))
P = pd.DataFrame(prof); P.round(5).to_csv(os.path.join(HERE, 'kappa_confound_profile.csv'), index=False)

# ---- A, B data ----
comp = pd.DataFrame(dict(group=[mv.group(c) for c in mv.CODES], n_sites=[mv.n_sites(c) for c in mv.CODES]))
CT = comp.groupby(['group', 'n_sites']).size().rename('n_interfaces').reset_index(); CT.to_csv(os.path.join(HERE, 'kappa_confound_composition.csv'), index=False)
F = pd.read_csv(os.path.join(HERE, 'fits.csv')); CV = pd.read_csv(os.path.join(REPO, 'Figure4_MultivalentModel/model_fit/cv_comparison.csv'))
k3 = CV[(CV.readout == 'intermediate') & (CV.n_params == 3)].iloc[0]
fn, fk = F[F.model.str.contains('no kappa')].iloc[0], F[F.model.str.contains('with kappa')].iloc[0]
B = pd.DataFrame([dict(model='valency term κ (Fig. 4 model)', n_params=3, cvRMS=k3.cvRMS, cvrho=k3.cvrho, kappa=k3.kappa_per_site, phi5=1.0),
                  dict(model='separate 5-/6-nt C_eff, no κ', n_params=4, cvRMS=fn.cvRMS, cvrho=fn.cvrho, kappa=0.0, phi5=fn.phi5),
                  dict(model='both', n_params=5, cvRMS=fk.cvRMS, cvrho=fk.cvrho, kappa=fk.kappa_per_site, phi5=fk.phi5)])
B.round(4).to_csv(os.path.join(HERE, 'kappa_confound_cv.csv'), index=False)

# ---- figure ----
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e6e5e0'; GC = {'7-nt': '#3B6AA0', '5-nt': '#E8A33D', '6-nt': '#4DA167', 'A/B': '#8a6fbf'}
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'svg.fonttype': 'none'})
fig, ax = plt.subplots(1, 2, figsize=(10.2, 4.1), gridspec_kw={'width_ratios': [.85, 1.15]})
b = ax[0]
b.bar(range(3), B.cvRMS, color=['#1baf7a', '#E8A33D', '#bdbdbd'], width=.6, zorder=2)
for i, r in B.iterrows(): b.text(i, r.cvRMS + .03, f'{r.cvRMS:.2f} °C\nρ {r.cvrho:.3f}', ha='center', fontsize=7.5)
b.set_xticks(range(3)); b.set_xticklabels(['valency κ\n(Fig. 4)', '5-/6-nt $C_{eff}$\n(no κ)', 'both'], fontsize=8)
b.set(ylabel='10-fold CV RMS error in $T_m$ (°C)', ylim=(0, B.cvRMS.max() * 1.3)); b.grid(True, axis='y', color=GRID, lw=.6, zorder=0)
b.set_title('A   equally good fits', loc='left', fontsize=9.5, weight='bold')
c = ax[1]; c2 = c.twinx(); c2.spines['right'].set_visible(True)
c.plot(P.kappa_per_site, P.phi5, color=GC['5-nt'], lw=2, marker='o', ms=4, label='fitted 5-nt factor φ5')
c2.plot(P.kappa_per_site, P.fit_RMS_C, color=INK2, lw=1.6, ls=(0, (3, 2)), label='RMS error (all data)')
c.axvline(k3.kappa_per_site, color='#1baf7a', lw=1, ls=':'); c.axhline(1, color=GRID, lw=1)
c.text(k3.kappa_per_site, c.get_ylim()[1] if False else P.phi5.max() * 1.02, ' Fig. 4 κ', fontsize=7, color='#1baf7a', va='bottom')
c.set(xlabel='valency term κ held fixed (per site)', ylabel='5-nt effective concentration relative to 7-nt (φ5)')
c2.set_ylabel('RMS error in $T_m$ (°C)', color=INK2); c2.set_ylim(P.fit_RMS_C.min() - .3, P.fit_RMS_C.min() + .6)
h1, l1 = c.get_legend_handles_labels(); h2, l2 = c2.get_legend_handles_labels(); c.legend(h1 + h2, l1 + l2, fontsize=7.5, frameon=False, loc='lower left')
c.grid(True, axis='y', color=GRID, lw=.6)
c.set_title('B   κ and φ5 trade off at constant error', loc='left', fontsize=9.5, weight='bold')
fig.tight_layout(); fig.savefig(os.path.join(HERE, 'kappa_confound.svg'))
if '--png' in sys.argv: fig.savefig(sys.argv[sys.argv.index('--png') + 1], dpi=150)
print(B.round(3).to_string(index=False)); print(P[['kappa_per_site', 'phi5', 'fit_RMS_C']].round(3).to_string(index=False))
