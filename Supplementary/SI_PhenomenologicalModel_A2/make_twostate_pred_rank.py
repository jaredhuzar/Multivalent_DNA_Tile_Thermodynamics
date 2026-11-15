"""
SI: two-state model A2 (distance weight on each bond's NN enthalpy and entropy), out-of-fold (10-fold CV, 236 interfaces):
  A predicted vs measured Tm,  B predicted vs measured rank.
Reads Phenomenological_Model/pred_vs_exp.csv and rank_rank.csv (written by TmPredictor_NewTm.ipynb).
Writes Phenomenological_Model/pred_and_rank.svg + pred_and_rank.csv.
"""
import os, numpy as np, pandas as pd
from scipy.stats import spearmanr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Phenomenological_Model')
P = pd.read_csv(os.path.join(D, 'pred_vs_exp.csv'), dtype={'code': str}); R = pd.read_csv(os.path.join(D, 'rank_rank.csv'), dtype={'code': str})
M = P.merge(R[['code', 'true_rank', 'pred_rank']], on='code')
rms = float(np.sqrt(np.mean((M.pred_Tm_oof - M.exp_Tm) ** 2))); rho = spearmanr(M.exp_Tm, M.pred_Tm_oof)[0]
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'svg.fonttype': 'none'})
fig, ax = plt.subplots(1, 2, figsize=(9.4, 4.2)); C, INK2 = '#e34948', '#52514e'
a = ax[0]; lim = [min(M.pred_Tm_oof.min(), M.exp_Tm.min()) - 1, max(M.pred_Tm_oof.max(), M.exp_Tm.max()) + 1]
a.scatter(M.pred_Tm_oof, M.exp_Tm, s=14, color=C, alpha=.55, edgecolor='none'); a.plot(lim, lim, color=INK2, lw=1)
a.set(xlim=lim, ylim=lim, xlabel='predicted $T_m$ (°C, held-out fold)', ylabel='measured $T_m$ (°C)')
a.text(.04, .96, f'RMSD {rms:.2f} °C\nSpearman ρ = {rho:.3f}\nn = {len(M)}', transform=a.transAxes, va='top', fontsize=8.5)
a.set_title('A   two-state model, out-of-fold', loc='left', fontsize=9.5, weight='bold')
b = ax[1]; nmax = len(M)
b.scatter(M.true_rank, M.pred_rank, s=14, color=C, alpha=.55, edgecolor='none'); b.plot([1, nmax], [1, nmax], color=INK2, lw=1)
b.set(xlabel='measured rank', ylabel='predicted rank', xlim=(0, nmax + 1), ylim=(0, nmax + 1))
b.text(.04, .96, f'Spearman ρ = {rho:.3f}', transform=b.transAxes, va='top', fontsize=8.5)
b.set_title('B   rank vs rank', loc='left', fontsize=9.5, weight='bold')
for x in ax: x.grid(True, color='#e6e5e0', lw=.6)
fig.tight_layout(); fig.savefig(os.path.join(D, 'pred_and_rank.svg')); plt.close(fig)
M.round(4).to_csv(os.path.join(D, 'pred_and_rank.csv'), index=False)
print(f'two-state pred+rank: n={len(M)}, RMSD {rms:.2f}, rho {rho:.3f}')
