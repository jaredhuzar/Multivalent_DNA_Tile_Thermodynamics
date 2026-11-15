"""
Why a fifth 7-nt site raises Tm but a tenth 5-nt site barely does (supports the Figure 2D -> 2F transition text).

For every add-one-site pair (lo -> hi) in the 4->5 (7-nt added) and 9->10 (5-nt added) regimes, docking readout:
  * C_eff is the value fitted to the lo construct's measured Tm (PanelH_BindingCooperativity/data.csv, docking rows);
  * each ADDED duplex i (the two duplexes of the new site) is evaluated at that C_eff:
        T*_i = dH_i / (dS_i + R ln C_eff)          melting temperature of the duplex alone with its partner at C_eff
        p_i  = K_i C_eff / (1 + K_i C_eff)  at the lo interface's Tm   (fraction of time formed while docked)
  * dTm at constant C_eff = Tm(hi, C_eff,lo) - Tm(lo)   vs the measured dTm (same data.csv).
Panels: A  T*_i - Tm(lo) of the added duplexes;  B  p_i of the added duplexes at Tm(lo);
        C  predicted (constant C_eff) vs measured dTm per pair.
Same NN thermo / salt / readout as make_uniform_ceff_panel.py (its header is executed).  Run after that script.
Writes PanelSI_AddedBondOccupancy/{figure.svg, duplexes.csv, pairs.csv}.
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(HERE)
_hdr = os.path.join(HERE, 'make_uniform_ceff_panel.py'); _src = open(_hdr).read(); _cut = '# ---- robust set'
assert _cut in _src, 'NN model header marker missing in make_uniform_ceff_panel.py'
ns = {'__file__': _hdr}; exec(_src[:_src.index(_cut)], ns); os.chdir(REPO)
R, duplexes = ns['R'], ns['duplexes']
OUT = os.path.join(HERE, 'PanelSI_AddedBondOccupancy'); os.makedirs(OUT, exist_ok=True)

P = pd.read_csv(os.path.join(HERE, 'PanelH_BindingCooperativity/data.csv'), dtype={'lo': str, 'hi': str})
P = P[(P.readout == 'docking') & P.regime.isin(['4->5', '9->10'])].reset_index(drop=True)
TM = pd.read_csv(os.path.join(HERE, 'PanelH_UniformCeff_Width/data.csv'), dtype={'code': str}).set_index('code').Tm_C
LAB = {'4->5': '4→5 sites\n(7-nt site added)', '9->10': '9→10 sites\n(5-nt site added)'}
rows = []
for k, r in P.iterrows():
    T = TM[r.lo] + 273.15; C = 10 ** r.log10Ceff_lo
    lo = [d[:2] for d in duplexes(r.lo)]; added = [d[:2] for d in duplexes(r.hi)]
    for d in lo: added.remove(d)
    for H, S in added:
        x = np.exp(-(H - T * S) / (R * T)) * C
        rows.append(dict(pair=k, lo=r.lo, hi=r.hi, regime=r.regime, Tm_lo_C=T - 273.15, log10Ceff_M=r.log10Ceff_lo,
                         dH_kcal=H, dS_kcal_K=S, Tstar_C=H / (S + R * np.log(C)) - 273.15,
                         Tstar_minus_Tm_C=H / (S + R * np.log(C)) - T, KCeff=x, occupancy=x / (1 + x)))
Dd = pd.DataFrame(rows); Dd.round(4).to_csv(os.path.join(OUT, 'duplexes.csv'), index=False)
P[['lo', 'hi', 'regime', 'log10Ceff_lo', 'dTm_constCeff_C', 'dTm_meas_C', 'shortfall_C', 'ddG_coop_kcal']].round(4) \
  .to_csv(os.path.join(OUT, 'pairs.csv'), index=False)

INK, INK2, GRID = '#0b0b0b', '#52514e', '#e6e5e0'; COL = {'4->5': '#3B6AA0', '9->10': '#E8A33D'}
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'svg.fonttype': 'none'})
fig, ax = plt.subplots(1, 3, figsize=(12.6, 4.1)); rng = np.random.default_rng(0); regs = ['4->5', '9->10']
def strip(a, col, ylabel, fmt, ref=None):
    for i, reg in enumerate(regs):
        v = Dd[Dd.regime == reg][col].values
        a.scatter(i + rng.uniform(-.18, .18, len(v)), v, s=11, color=COL[reg], alpha=.55, edgecolor='none', zorder=2)
        mn = v.mean(); a.plot([i - .3, i + .3], [mn, mn], color=INK, lw=2.2, zorder=3)
        a.text(i + .33, mn, fmt(mn), va='center', fontsize=8, color=INK)
    if ref is not None: a.axhline(ref, color=INK2, lw=1, ls='--', zorder=1)
    a.set_xticks(range(2)); a.set_xticklabels([f'{LAB[r]}\n(n = {(Dd.regime == r).sum()} duplexes)' for r in regs], fontsize=8)
    a.set_xlim(-.5, 1.75); a.set_ylabel(ylabel); a.grid(True, axis='y', color=GRID, lw=.6, zorder=0)
a = ax[0]; strip(a, 'Tstar_minus_Tm_C', 'added duplex: $T^*$ at $C_{eff}$ − interface $T_m$ (°C)', lambda m: f'mean {m:+.1f} °C', 0)
a.text(-.47, .4, '↑ melts above the interface', ha='left', va='bottom', fontsize=7, color=INK2)
a.text(-.47, -.4, '↓ melts below', ha='left', va='top', fontsize=7, color=INK2)
a.set_title('A   added duplexes alone at $C_{eff}$', loc='left', fontsize=9.5, weight='bold')
a = ax[1]; strip(a, 'occupancy', 'fraction of time formed at interface $T_m$\n$K_iC_{eff}/(1+K_iC_{eff})$', lambda m: f'mean {m:.2f}', .5)
a.set_ylim(0, 1.05); a.set_title('B   how often the added duplexes are formed', loc='left', fontsize=9.5, weight='bold')
a = ax[2]
for i, reg in enumerate(regs):
    s = P[P.regime == reg]
    for j, (col, mk, lab) in enumerate((('dTm_constCeff_C', 'o', 'predicted, $C_{eff}$ held fixed'), ('dTm_meas_C', 's', 'measured'))):
        v = s[col].values; x = i + (-.17 if j == 0 else .17)
        a.scatter(x + rng.uniform(-.07, .07, len(v)), v, s=10, color=COL[reg], alpha=.5 if j else .8,
                  marker=mk, edgecolor='none' if j == 0 else COL[reg], facecolor=COL[reg] if j == 0 else 'none', zorder=2)
        md = np.median(v); a.plot([x - .13, x + .13], [md, md], color=INK, lw=2.2, zorder=3)
        a.text(x, max(s.dTm_constCeff_C.max(), s.dTm_meas_C.max()) + .4, f'{md:+.1f}', ha='center', fontsize=8)
a.scatter([], [], marker='o', color=INK2, label='predicted, $C_{eff}$ held fixed'); a.scatter([], [], marker='s', facecolor='none', edgecolor=INK2, label='measured')
a.axhline(0, color=INK2, lw=.8); a.set_xticks(range(2)); a.set_xticklabels([f'{LAB[r]}\n(n = {(P.regime == r).sum()} pairs)' for r in regs], fontsize=8)
a.set(xlim=(-.5, 1.5), ylabel='Δ$T_m$ from the added site (°C)'); a.grid(True, axis='y', color=GRID, lw=.6, zorder=0)
a.legend(fontsize=7.5, frameon=False, loc='upper right')
a.set_title('C   predicted vs measured gain (medians)', loc='left', fontsize=9.5, weight='bold')
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'figure.svg'))
import sys
if '--png' in sys.argv: fig.savefig(sys.argv[sys.argv.index('--png') + 1], dpi=150)
g = Dd.groupby('regime'); q = P.groupby('regime')
print('added-bond occupancy: ' + '; '.join(
    f"{reg}: T*-Tm mean {g.get_group(reg).Tstar_minus_Tm_C.mean():+.1f} C (median {g.get_group(reg).Tstar_minus_Tm_C.median():+.1f}), "
    f"occupancy mean {g.get_group(reg).occupancy.mean():.2f} (median {g.get_group(reg).occupancy.median():.2f}), "
    f"dTm const-Ceff median {q.get_group(reg).dTm_constCeff_C.median():+.2f}, measured {q.get_group(reg).dTm_meas_C.median():+.2f}" for reg in regs))
