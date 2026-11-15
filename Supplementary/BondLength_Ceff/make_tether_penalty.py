"""
Per-bond tethering penalty relative to a free duplex at 1 M:  dG_teth = RT ln(1 M / C_eff)  at 37 C.

  A  every interface's own C_eff (multivalent model, one uniform C_eff per interface, intermediate readout;
     Supplementary/Ceff_Geometry_OLS/data.csv), grouped by the interface's bond length, vs the worm-like-chain (WLC)
     tether estimate for each length (Table S1 model B4: C_WLC = exp(dS_teth / R)).
  B  bond-type effective concentration fitted inside the Figure 4 model (BondLength_Ceff/fits.csv, no-kappa fit),
     evaluated at the median duplex mean distance of the 236 interfaces, with 95% CI from phi -- vs the WLC estimate.
Run after Ceff_Geometry_OLS/make_ceff_ols.py and BondLength_Ceff/make_bondlength_ceff.py.
Writes tether_penalty.svg + tether_penalty.csv here.
"""
import os, sys, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(REPO, 'Figure4_MultivalentModel')); import mv_model as mv
RT = 1.98720425864e-3 * 310.15                              # kcal/mol at 37 C
pen = lambda C: RT * np.log(1.0 / C)                        # C in M

D = pd.read_csv(os.path.join(REPO, 'Supplementary/Ceff_Geometry_OLS/data.csv'), dtype={'code': str})
D['penalty_kcal'] = pen(10 ** D.log10_Ceff_M)
F = pd.read_csv(os.path.join(HERE, 'fits.csv')); fn = F[F.model.str.contains('no kappa')].iloc[0]
WLC = {'5-nt': fn.WLC_Ceff_5nt_mM * 1e-3, '6-nt': fn.WLC_Ceff_6nt_mM * 1e-3, '7-nt': fn.WLC_Ceff_7nt_mM * 1e-3}
A = mv.arrays(mv.UM['Bl']); dmed = float(np.median(A['D'][A['mask']]))
base = 10 ** fn.log10_Cref_M * np.exp(-dmed / fn.lambda_nm)
FIT = {'7-nt': (base, base, base), '5-nt': (base * fn.phi5, base * fn.phi5_lo, base * fn.phi5_hi),
       '6-nt': (base * fn.phi6, base * fn.phi6_lo, base * fn.phi6_hi)}

GRP = ['5-nt', '6-nt', '7-nt', 'A/B']; LEN = {'5-nt': '5-nt', '6-nt': '6-nt', '7-nt': '7-nt', 'A/B': '7-nt'}
GC = {'7-nt': '#3B6AA0', '5-nt': '#E8A33D', '6-nt': '#4DA167', 'A/B': '#8a6fbf'}; INK, INK2, GRID, WC = '#0b0b0b', '#52514e', '#e6e5e0', '#e34948'
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'svg.fonttype': 'none'})
fig, ax = plt.subplots(1, 2, figsize=(10.4, 4.1), gridspec_kw={'width_ratios': [1.25, 1]}); rng = np.random.default_rng(0)
a = ax[0]; rows = []
for i, g in enumerate(GRP):
    v = D[D.group == g].penalty_kcal.values
    a.scatter(i + rng.uniform(-.17, .17, len(v)), v, s=12, color=GC[g], alpha=.6, edgecolor='none', zorder=2)
    md = np.median(v); a.plot([i - .28, i + .28], [md, md], color=INK, lw=2.2, zorder=3); a.text(i + .3, md, f'{md:.2f}', va='center', fontsize=7.5)
    w = pen(WLC[LEN[g]]); a.plot([i - .28, i + .28], [w, w], color=WC, lw=2.2, ls=(0, (3, 1.5)), zorder=3)
    rows.append(dict(panel='A', group=g, n=len(v), median_penalty_kcal=md, WLC_penalty_kcal=w))
a.plot([], [], color=INK, lw=2.2, label='median of fitted interfaces'); a.plot([], [], color=WC, lw=2.2, ls=(0, (3, 1.5)), label='worm-like-chain tether estimate')
a.set_xticks(range(len(GRP))); a.set_xticklabels([f'{g}\n(n={(D.group == g).sum()})' for g in GRP])
a.set(ylabel=r'per-bond penalty  $RT\,\ln(1\,\mathrm{M}/C_{eff})$  (kcal/mol, 37 °C)', ylim=(2, 6))
a.legend(fontsize=7.5, frameon=False, loc='upper left'); a.grid(True, axis='y', color=GRID, lw=.6, zorder=0)
a.set_title('A   each interface\'s own $C_{eff}$, by bond length', loc='left', fontsize=9.5, weight='bold')
b = ax[1]; L = ['5-nt', '6-nt', '7-nt']
for i, g in enumerate(L):
    c, clo, chi = FIT[g]; p, plo, phi = pen(c), pen(chi), pen(clo)
    b.errorbar(i - .12, p, yerr=[[p - plo], [phi - p]], fmt='o', color=GC[g], ms=8, capsize=3, lw=1.4, zorder=3)
    w = pen(WLC[g]); b.scatter(i + .12, w, marker='D', s=45, color=WC, zorder=3)
    rows.append(dict(panel='B', group=g, n=np.nan, fitted_penalty_kcal=p, fitted_lo=plo, fitted_hi=phi, WLC_penalty_kcal=w,
                     at_mean_distance_nm=dmed))
b.errorbar([], [], fmt='o', color=INK2, label=f'fitted in the multivalent model\n(at the median distance, {dmed:.0f} nm)')
b.scatter([], [], marker='D', color=WC, label='worm-like-chain estimate')
b.set_xticks(range(3)); b.set_xticklabels(L); b.set(ylabel=r'$RT\,\ln(1\,\mathrm{M}/C_{eff})$  (kcal/mol)', ylim=(2, 6), xlim=(-.5, 2.5))
b.legend(fontsize=7.5, frameon=False, loc='upper left'); b.grid(True, axis='y', color=GRID, lw=.6, zorder=0)
b.set_title('B   bond-length dependence: fit vs tether model', loc='left', fontsize=9.5, weight='bold')
fig.tight_layout(); fig.savefig(os.path.join(HERE, 'tether_penalty.svg'))
if '--png' in sys.argv: fig.savefig(sys.argv[sys.argv.index('--png') + 1], dpi=150)
pd.DataFrame(rows).round(4).to_csv(os.path.join(HERE, 'tether_penalty.csv'), index=False)
print('tether penalty:', pd.DataFrame(rows).round(2).to_string(index=False), f'\n  overall median {D.penalty_kcal.median():.2f} kcal/mol')
