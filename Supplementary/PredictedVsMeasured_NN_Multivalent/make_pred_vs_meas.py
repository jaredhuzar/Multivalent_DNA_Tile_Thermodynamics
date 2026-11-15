"""
SI figure (replaces S12_PredictedVsMeasured_Tm): sequence-predicted vs measured Tm.
  A  additive nearest-neighbor two-state Tm: all duplexes summed, each referenced to 1 M,
     Tm = sum dH / (sum dS + R ln(C/2)), C = 4.5 nM per tile (salt-corrected NN, SantaLucia 1998)
  B  multivalent model of Figure 4, out-of-fold (10-fold CV) -- the same sequences, with bonds after the first
     at the fitted effective concentration
Canonical 236 interfaces.  Writes figure.svg, figure.png, data.csv here.
"""
import os, sys, numpy as np, pandas as pd
from scipy.stats import spearmanr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(REPO, 'Figure4_MultivalentModel')); import mv_model as mv

A = mv.arrays(mv.UM['Bl']); nn = mv.nn_twostate_tm(A) - 273.15
O = pd.read_csv(os.path.join(REPO, 'Figure4_MultivalentModel/model_fit/oof_predictions.csv'), dtype={'code': str})
assert list(O.code) == list(mv.CODES)
D = pd.DataFrame(dict(code=mv.CODES, group=O.group, n_sites=O.n_sites, Tm_meas_C=mv.Y_C, Tm_NN_additive_C=nn,
                      Tm_multivalent_oof_C=O.Tm_pred_oof_C))
D.round(4).to_csv(os.path.join(HERE, 'data.csv'), index=False)
GC = {'7-nt': '#3B6AA0', '5-nt': '#E8A33D', '6-nt': '#4DA167', 'A/B': '#8a6fbf'}; INK2, GRID = '#52514e', '#e6e5e0'
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'svg.fonttype': 'none'})
fig, ax = plt.subplots(1, 2, figsize=(10, 4.4))
for a, col, ttl in ((ax[0], 'Tm_NN_additive_C', 'A   additive nearest-neighbor sum (1 M reference)'),
                    (ax[1], 'Tm_multivalent_oof_C', 'B   multivalent model (Fig. 4), out-of-fold')):
    x = D[col].values; y = D.Tm_meas_C.values
    for g, c in GC.items():
        m = (D.group == g).values; a.scatter(x[m], y[m], s=15, color=c, alpha=.75, edgecolor='white', lw=.3, label=f'{g} (n={m.sum()})')
    lo, hi = min(x.min(), y.min()) - 1, max(x.max(), y.max()) + 1
    a.plot([lo, hi], [lo, hi], color=INK2, lw=1); a.set(xlim=(lo, hi), ylim=(lo, hi) if col != 'Tm_NN_additive_C' else (28, 56),
                                                         xlabel='predicted $T_m$ (°C)', ylabel='measured $T_m$ (°C)')
    rms = np.sqrt(np.mean((x - y) ** 2)); rho = spearmanr(x, y)[0]
    a.text(.04, .96, f'RMS {rms:.1f} °C\nSpearman ρ = {rho:.2f}', transform=a.transAxes, va='top', fontsize=8.5)
    a.set_title(ttl, loc='left', fontsize=9.5, weight='bold'); a.grid(True, color=GRID, lw=.6)
ax[1].legend(fontsize=7, frameon=False, loc='lower right')
fig.tight_layout(); fig.savefig(os.path.join(HERE, 'figure.svg')); fig.savefig(os.path.join(HERE, 'figure.png'), dpi=160)
print('NN additive: median by group', D.groupby('group').Tm_NN_additive_C.median().round(1).to_dict(),
      f'rho {spearmanr(D.Tm_NN_additive_C, D.Tm_meas_C)[0]:.2f}; multivalent rho {spearmanr(D.Tm_multivalent_oof_C, D.Tm_meas_C)[0]:.3f}')
