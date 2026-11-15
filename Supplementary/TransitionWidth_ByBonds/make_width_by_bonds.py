"""
SI figure: raw transition width (Boltzmann dx, degrees C) for 7-nt (4, 5 sites) and 5-nt (9, 10 sites) interfaces,
with the all-or-none width expected if the whole interface melted as one unit (dx = 1.5 R Tm^2 / |dH_NN|).
Source: Figure2_Updated/PanelD_CooperativeUnit/data.csv (Figure 2 robust set).  Writes figure.svg + data.csv here.
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(os.path.dirname(HERE))
D = pd.read_csv(os.path.join(REPO, 'Figure2_Updated/PanelD_CooperativeUnit/data.csv'), dtype={'code': str})
GRP = [(7, 4), (7, 5), (5, 9), (5, 10)]; BL = {5: '#E8A33D', 7: '#3B6AA0'}
D = D[[ (int(a), int(b)) in GRP for a, b in zip(D.nt, D.n_sites)]].copy()
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'svg.fonttype': 'none'})
fig, ax = plt.subplots(figsize=(5.4, 4.2)); rng = np.random.default_rng(0); med = {}
for i, (nt, N) in enumerate(GRP):
    s = D[(D.nt == nt) & (D.n_sites == N)]; v = s.width_dx_C.values
    ax.scatter(i + rng.uniform(-.17, .17, len(v)), v, s=13, color=BL[nt], alpha=.7, edgecolor='none', zorder=2)
    m = np.median(v); q = np.percentile(v, [25, 75]); med[(nt, N)] = m
    ax.plot([i - .3, i + .3], [m, m], color='#0b0b0b', lw=2.2, zorder=3); ax.plot([i, i], q, color='#0b0b0b', lw=1.1, zorder=3)
    a = s.allornone_dx_C.median(); ax.plot([i - .3, i + .3], [a, a], color='#c0392b', lw=1.4, ls=(0, (3, 2)), zorder=3)
    ax.text(i + .33, m, f'{m:.2f}', va='center', fontsize=7.5)
for (i, j) in ((0, 1), (2, 3)):
    d = med[GRP[j]] - med[GRP[i]]
    ax.annotate(f'Δ {d:+.2f} °C', xy=((i + j) / 2, max(med[GRP[i]], med[GRP[j]]) + 1.3), ha='center', fontsize=8)
ax.plot([], [], color='#c0392b', lw=1.4, ls=(0, (3, 2)), label='all-or-none width (median)')
ax.axvline(1.5, color='#cccccc', lw=1)
ax.set_xticks(range(4)); ax.set_xticklabels([f'{N} sites\n{nt}-nt\n(n={((D.nt == nt) & (D.n_sites == N)).sum()})' for nt, N in GRP])
ax.set(ylabel='transition width δ (°C)', ylim=(0, D.width_dx_C.max() * 1.15)); ax.legend(fontsize=7.5, frameon=False, loc='upper right')
ax.grid(True, axis='y', color='#e6e5e0', lw=.6, zorder=0)
ax.set_title('Transition width by bond length and number of sites', loc='left', fontsize=9.5, weight='bold')
fig.tight_layout(); fig.savefig(os.path.join(HERE, 'figure.svg'))
D[['code', 'nt', 'n_sites', 'width_dx_C', 'allornone_dx_C']].round(4).to_csv(os.path.join(HERE, 'data.csv'), index=False)
print('width medians:', {f'{k[0]}nt/{k[1]}': round(v, 2) for k, v in med.items()})
