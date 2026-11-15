"""
Figure 4 (multivalent model).  Run fit_model.py first (writes model_fit/).  Run from anywhere.  --png PATH: preview.

  A  model (schematic placeholder + equations + fitted parameters)
  B  out-of-fold predicted vs measured Tm, 236 interfaces (10-fold CV, folds over interfaces)
  C  what each ingredient buys: CV RMS of the nested variants (+ the SI phenomenological models for reference)
  D  THEORY  Tm vs number of sites (clustered, contiguous windows containing the reporter site 5), 7-nt and 5-nt
             bonds; dashed = same model without the valency term; markers = measured / fitted group medians
  E  THEORY  Tm vs mean gap for every four-site 7-nt arrangement on the 11-site lattice (site 5 = 5-nt reporter);
             measured medians overlaid; dashed = one uniform C_eff
  F  THEORY  melting curves theta(T) of SYNTHETIC four-site interfaces (identical average 7-nt bonds) at gap 6-20 nm
  G  THEORY  per-duplex occupancy K C_eff / (1 + K C_eff) at the interface Tm vs the duplex's mean distance
  H  THEORY  design map: Tm over number of sites x uniform spacing (7-nt bonds, 5-nt reporter at site 5)

Theory curves use the fitted main model (3 parameters, intermediate readout) and the real per-position sequences.
Writes figure.svg and theory/*.csv here.
"""
import os, sys, itertools, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mv_model as mv

HERE = mv.HERE; FIT = os.path.join(HERE, 'model_fit'); TH = os.path.join(HERE, 'theory'); os.makedirs(TH, exist_ok=True)
P = pd.read_csv(os.path.join(FIT, 'params.csv')).set_index('parameter')
CV = pd.read_csv(os.path.join(FIT, 'cv_comparison.csv'))
O = pd.read_csv(os.path.join(FIT, 'oof_predictions.csv'), dtype={'code': str})
p = P.point_estimate.values.astype(float)                       # [log10 C_ref, lambda, kappa]
p_k0 = np.array([p[0], p[1], 0.0])
_u = CV[(CV.readout == 'intermediate') & (CV.n_params == 1)].iloc[0]; p_uni = np.array([_u.log10_Cref_M])
C0 = 273.15
ci = 'boot' if 'boot_CI95_low' in P else 'GN'
lo_ci, hi_ci = P[f'{ci}_CI95_low'].values, P[f'{ci}_CI95_high'].values

def code_of(sites, ch):
    c = ['0'] * 11
    for s in sites: c[s] = ch
    c[5] = '5'; return ''.join(c)
def gap_nm(sites): s = sorted(sites); return (s[-1] - s[0]) / (len(s) - 1) * 6.0 if len(s) > 1 else np.nan

# ---------------- theory D: Tm vs number of sites ----------------
rowsD = []
for typ in ('7', '5'):
    for N in range(1, 12):
        wins = [list(range(s, s + N)) for s in range(0, 12 - N) if s <= 5 < s + N]
        A = mv.arrays_from_codes([code_of(w, typ) for w in wins])
        t, t0 = mv.tm(p, A) - C0, mv.tm(p_k0, A) - C0
        rowsD.append(dict(bond=f'{typ}-nt', n_sites=N, n_windows=len(wins), Tm_mean_C=t.mean(), Tm_min_C=t.min(),
                          Tm_max_C=t.max(), Tm_no_valency_mean_C=t0.mean()))
D = pd.DataFrame(rowsD); D.round(3).to_csv(os.path.join(TH, 'Tm_vs_NumSites.csv'), index=False)

# ---------------- theory E/F/G: all four-site 7-nt arrangements ----------------
ARR = [sorted([5, *c]) for c in itertools.combinations([s for s in range(11) if s != 5], 3)]
CODES4 = [code_of(a, '7') for a in ARR]; A4 = mv.arrays_from_codes(CODES4)
T4 = mv.tm(p, A4); T4u = mv.tm(p_uni, A4); W4 = mv.width(p, A4)
_, _, OCC = mv.state(T4, p, A4)
E = pd.DataFrame(dict(code=CODES4, mean_gap_nm=[gap_nm(a) for a in ARR], Tm_model_C=T4 - C0,
                      Tm_uniformCeff_C=T4u - C0, width_model_C=W4))
E.round(3).to_csv(os.path.join(TH, 'Tm_vs_MeanGap_4sites.csv'), index=False)
G = pd.DataFrame([dict(code=CODES4[i], duplex=k, site=int(round(A4['D'][i, k] * 0)), mean_distance_nm=A4['D'][i, k],
                       reporter=bool(A4['rep'][i, k]), occupancy_at_Tm=OCC[i, k])
                  for i in range(len(CODES4)) for k in range(8)]).drop(columns='site')
G.round(4).to_csv(os.path.join(TH, 'Occupancy_vs_Distance_4sites.csv'), index=False)
# Panel F is SYNTHETIC: four identical sites (every duplex = the mean salt-corrected NN dH/dS of the 22 7-nt duplexes),
# equally spaced at gap g (twin duplexes +-1.5 nm), so only geometry changes between curves; reporter = 2nd site.
_th = np.array([mv.UM['TH']['7'][h] for h in range(1, 23)])
H7 = _th[:, 0].mean(); S7 = (_th[:, 1] / 1000.0 + 0.368e-3 * 6 * mv.LN_NAEQ).mean()
def synth(g):
    z = np.array([k * g + t for k in range(4) for t in (-1.5, 1.5)])
    A = dict(H=np.full((1, 8), H7), S=np.full((1, 8), S7), mask=np.ones((1, 8), bool), N=np.array([4.0]),
             D=np.array([[np.mean(np.abs(np.delete(z, i) - z[i])) for i in range(8)]]),
             rep=np.array([[False, False, True, True, False, False, False, False]]))
    A['nrep'] = A['rep'].sum(1); A['hasrep'] = A['rep'].sum(1) > 0
    return A
GAPS = np.arange(6.0, 21.0, 2.0)
Tgrid = np.arange(20.0, 70.0, 0.1) + C0
SYN = {}
for g_ in GAPS:
    A1 = synth(g_); tmK = mv.tm(p, A1)[0]; wd = mv.width(p, A1)[0]
    SYN[g_] = (tmK, wd, mv.theta(Tgrid, p, mv._sub(A1, np.zeros(len(Tgrid), int))))
pd.DataFrame({'T_C': Tgrid - C0, **{f'theta_gap{int(g_)}nm': SYN[g_][2] for g_ in GAPS}}).round(5) \
  .to_csv(os.path.join(TH, 'MeltCurves_vs_Spacing.csv'), index=False)
pd.DataFrame(dict(gap_nm=GAPS, Tm_C=[SYN[g][0] - C0 for g in GAPS], width_C=[SYN[g][1] for g in GAPS],
                  bond_dH_kcal=H7, bond_dS_kcal_K=S7)).round(4).to_csv(os.path.join(TH, 'MeltCurves_vs_Spacing_summary.csv'), index=False)

# ---------------- theory H: design map ----------------
rowsH = []
for N in range(2, 12):
    for k in range(1, 6):
        sets = [[s + k * j for j in range(N)] for s in range(11)]
        sets = [q for q in sets if q[-1] <= 10 and 5 in q]
        if not sets: continue
        t = mv.tm(p, mv.arrays_from_codes([code_of(q, '7') for q in sets])) - C0
        rowsH.append(dict(n_sites=N, spacing_sites=k, spacing_nm=6.0 * k, n_placements=len(sets), Tm_mean_C=t.mean()))
H = pd.DataFrame(rowsH); H.round(3).to_csv(os.path.join(TH, 'DesignMap_Tm.csv'), index=False)

# ---------------- figure ----------------
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e6e5e0"
GC = {'7-nt': '#3B6AA0', '5-nt': '#E8A33D', '6-nt': '#4DA167', 'A/B': '#8a6fbf'}
MOD, UNI = '#1baf7a', '#8a8984'
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'svg.fonttype': 'none',
                     'axes.edgecolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2})
fig = plt.figure(figsize=(15.5, 14.2))
gs = GridSpec(3, 3, figure=fig, hspace=.48, wspace=.32, top=.95, bottom=.05, left=.055, right=.975)
def title(ax, s): ax.set_title(s, loc='left', fontsize=9.6, weight='bold')
def grid(ax): ax.grid(True, axis='y', color=GRID, lw=.6, zorder=0)

# A
a = fig.add_subplot(gs[0, 0:2]); a.axis('off'); title(a, 'A   multivalent binding model')
a.add_patch(plt.Rectangle((.0, .0), 1, 1, fill=False, ec='#bbb', ls='--', lw=1.1, transform=a.transAxes))
eq = [r'$K_a=\dfrac{1}{C_{\mathrm{ref}}}\left[\prod_{i=1}^{2N}\left(1+K_i\,C_{\mathrm{eff},i}\right)-1\right]$',
      r'$C_{\mathrm{eff},i}=C_{\mathrm{ref}}\;e^{-d_i/\lambda}\;e^{-\kappa\,(N-4)}$',
      r'$\theta=f\,\left[w+(1-w)\,p_{\mathrm{rep}}\right],\qquad p_i=\dfrac{K_iC_{\mathrm{eff},i}}{1+K_iC_{\mathrm{eff},i}}$']
for k, s in enumerate(eq): a.text(.04, .80 - .21 * k, s, fontsize=12.5, transform=a.transAxes, va='center')
pt = [f'$C_{{\\mathrm{{ref}}}}$ = {10 ** p[0] * 1e3:.0f} mM  ({10 ** lo_ci[0] * 1e3:.0f}–{10 ** hi_ci[0] * 1e3:.0f})',
      f'λ = {p[1]:.1f} nm  ({lo_ci[1]:.1f}–{hi_ci[1]:.1f})',
      f'κ = {p[2]:.3f} per site  ({lo_ci[2]:.3f}–{hi_ci[2]:.3f})',
      f'w = {mv.W_INT} (fixed, Fig. 2)']
a.text(.64, .80, 'fitted (95% CI, ' + ('bootstrap' if ci == 'boot' else 'Gauss–Newton') + ')', fontsize=8.5, color=INK2,
       transform=a.transAxes, weight='bold')
for k, s in enumerate(pt): a.text(.64, .69 - .1 * k, s, fontsize=9.5, transform=a.transAxes)
a.text(.64, .22, '$d_i$: mean distance of duplex $i$ to the other duplexes\n$N$: occupied sites;  $K_i$: salt-corrected NN\n'
       '[schematic of tile edge, bonds and $C_{\\mathrm{eff}}$ decay to be added]', fontsize=7.6, color=INK2, transform=a.transAxes, va='top')

# B
b = fig.add_subplot(gs[0, 2]); title(b, 'B   out-of-fold prediction, 236 interfaces')
for gname in ('7-nt', '5-nt', '6-nt', 'A/B'):
    s = O[O.group == gname]
    b.scatter(s.Tm_pred_oof_C, s.Tm_meas_C, s=16, color=GC[gname], alpha=.75, edgecolor='white', lw=.3, label=f'{gname} (n={len(s)})', zorder=3)
lim = [min(O.Tm_pred_oof_C.min(), O.Tm_meas_C.min()) - 1, max(O.Tm_pred_oof_C.max(), O.Tm_meas_C.max()) + 1]
b.plot(lim, lim, color=INK2, lw=1, zorder=1); b.set(xlim=lim, ylim=lim, xlabel='predicted $T_m$ (°C, held-out fold)', ylabel='measured $T_m$ (°C)')
main = CV[(CV.readout == 'intermediate') & (CV.n_params == 3)].iloc[0]
b.text(.04, .96, f'RMS {main.cvRMS:.2f} °C\nSpearman ρ = {main.cvrho:.3f}', transform=b.transAxes, va='top', fontsize=8.5)
b.legend(fontsize=7, frameon=False, loc='lower right'); grid(b)

# C
c = fig.add_subplot(gs[1, 0]); title(c, 'C   what each ingredient buys')
mi = CV[CV.readout == 'intermediate'].sort_values('n_params')          # main model readout only (docking + SI models in model_fit/)
mi[['model', 'n_params', 'cvRMS', 'cvrho', 'log10_Cref_M', 'lambda_nm', 'kappa_per_site']].round(4) \
  .to_csv(os.path.join(TH, 'NestedModelComparison.csv'), index=False)
x = np.arange(3)
c.bar(x, mi.cvRMS, width=.6, color=[UNI, '#7fd0b0', MOD], edgecolor='none', zorder=2)
for xi, (rm, rh) in enumerate(zip(mi.cvRMS, mi.cvrho)): c.text(xi, rm + .05, f'{rm:.2f} °C\nρ {rh:.2f}', ha='center', fontsize=7.4)
c.set_xticks(x); c.set_xticklabels(['one\nuniform\n$C_{eff}$', '+ distance\n(λ)', '+ valency\n(κ)'], fontsize=7.6)
c.set(ylabel='10-fold CV RMS error in $T_m$ (°C)', ylim=(0, mi.cvRMS.max() * 1.25), xlim=(-.6, 2.6))
grid(c)

# D
d = fig.add_subplot(gs[1, 1]); title(d, 'D   Tm vs number of sites (clustered)')
for typ in ('7-nt', '5-nt'):
    s = D[D.bond == typ]; col = GC[typ]
    d.fill_between(s.n_sites, s.Tm_min_C, s.Tm_max_C, color=col, alpha=.15, lw=0)
    d.plot(s.n_sites, s.Tm_mean_C, color=col, lw=2, label=f'{typ} model')
    d.plot(s.n_sites, s.Tm_no_valency_mean_C, color=col, lw=1.2, ls=(0, (3, 2)))
    m = O[O.group == typ].groupby('n_sites')
    for N_, q in m:
        if len(q) < 3: continue
        md_ = q.Tm_meas_C.median(); i25, i75 = np.percentile(q.Tm_meas_C, [25, 75])
        d.errorbar(N_, md_, yerr=[[md_ - i25], [i75 - md_]], fmt='o', color=col, mec=INK, mew=.8, ms=6, capsize=2.5, zorder=4)
        d.scatter(N_ + .18, q.Tm_pred_fit_C.median(), marker='o', s=34, facecolor='white', edgecolor=col, lw=1.4, zorder=4)
d.plot([], [], color=INK2, lw=1.2, ls=(0, (3, 2)), label='model without valency term (κ = 0)')
d.errorbar([], [], fmt='o', color=INK2, mec=INK, label='measured median (IQR), all spacings')
d.scatter([], [], facecolor='white', edgecolor=INK2, lw=1.4, label='model, same constructs')
d.set(xlabel='occupied binding sites  (site 5 = 5-nt reporter)', ylabel='$T_m$ (°C)', xticks=range(1, 12))
d.legend(fontsize=6.6, frameon=False, loc='lower right'); grid(d)

# E
e = fig.add_subplot(gs[1, 2]); title(e, 'E   Tm vs spacing (four sites, 7-nt)')
ug = np.sort(E.mean_gap_nm.round(2).unique())
med = lambda df, col: np.array([df[np.isclose(df.mean_gap_nm, u)][col].median() for u in ug])
e.scatter(E.mean_gap_nm, E.Tm_model_C, s=9, color=MOD, alpha=.35, edgecolor='none', zorder=2)
e.plot(ug, med(E, 'Tm_model_C'), color=MOD, lw=2.2, label='model (all 120 arrangements)', zorder=3)
e.plot(ug, med(E, 'Tm_uniformCeff_C'), color=UNI, lw=1.4, ls=(0, (3, 2)), label='one uniform $C_{eff}$', zorder=3)
m4 = O[(O.group == '7-nt') & (O.n_sites == 4)].copy(); m4['mean_gap_nm'] = [gap_nm([i for i, ch in enumerate(cc) if ch != '0']) for cc in m4.code]
ugm = np.sort(m4.mean_gap_nm.round(2).unique())
e.plot(ugm, [m4[np.isclose(m4.mean_gap_nm, u)].Tm_meas_C.median() for u in ugm], color=INK, lw=1.6, marker='o', ms=4,
       label=f'measured median (n = {len(m4)})', zorder=4)
e.set(xlabel='mean gap between bonded sites (nm)', ylabel='$T_m$ (°C)'); e.legend(fontsize=7, frameon=False, loc='upper right'); grid(e)

# F
f_ = fig.add_subplot(gs[2, 0]); title(f_, 'F   identical bonds: spreading lowers $T_m$')
cols = dict(zip(GAPS, plt.cm.viridis(np.linspace(0, .9, len(GAPS)))))   # dark = clustered, light = spread
for g_ in GAPS:
    tmK, wd, th = SYN[g_]
    f_.plot(Tgrid - C0, th, color=cols[g_], lw=2.0, label=f'{g_:.0f} nm: $T_m$ {tmK - C0:.1f} °C, width {wd:.2f} °C')
    f_.plot(tmK - C0, .5, 'o', color=cols[g_], mec=INK, ms=4.5)
_tms = [SYN[g][0] - C0 for g in GAPS]
f_.axhline(.5, color=GRID, lw=1); f_.set(xlabel='temperature (°C)', ylabel='bound fraction θ (FRET readout)',
                                          xlim=(min(_tms) - 9, max(_tms) + 7))
f_.text(.99, .97, 'synthetic: 4 identical 7-nt sites,\nequally spaced', transform=f_.transAxes, ha='right', va='top', fontsize=6.8, color=INK2)
f_.legend(fontsize=6.3, frameon=False, loc='lower left', title='mean gap', title_fontsize=6.5); grid(f_)

# G
g = fig.add_subplot(gs[2, 1]); title(g, 'G   distant bonds are formed less at the interface $T_m$')
nr = G[~G.reporter]; rp = G[G.reporter]
g.scatter(nr.mean_distance_nm, nr.occupancy_at_Tm, s=6, color=GC['7-nt'], alpha=.25, edgecolor='none', label='7-nt duplexes')
g.scatter(rp.mean_distance_nm, rp.occupancy_at_Tm, s=6, color=GC['5-nt'], alpha=.35, edgecolor='none', label='5-nt reporter duplexes')
bins = np.arange(3, 40, 3.0)
for s_, col in ((nr, GC['7-nt']), (rp, GC['5-nt'])):
    k_ = np.digitize(s_.mean_distance_nm, bins); xs, ys = [], []
    for b_ in np.unique(k_):
        q = s_[k_ == b_]
        if len(q) >= 8: xs.append(q.mean_distance_nm.median()); ys.append(q.occupancy_at_Tm.median())
    g.plot(xs, ys, color=col, lw=2.2, mec=INK)
g.set(xlabel='mean distance of the duplex to the other duplexes (nm)', ylabel='occupancy $K_iC_{eff,i}/(1+K_iC_{eff,i})$ at $T_m$', ylim=(0, 1.02))
g.legend(fontsize=7, frameon=False, loc='lower left', markerscale=2.5); grid(g)

# H
h = fig.add_subplot(gs[2, 2]); title(h, 'H   design map: predicted Tm (7-nt bonds)')
Ns = list(range(2, 12)); Ks = list(range(1, 6))
Z = np.full((len(Ks), len(Ns)), np.nan)
for _, r in H.iterrows(): Z[Ks.index(int(r.spacing_sites)), Ns.index(int(r.n_sites))] = r.Tm_mean_C
im = h.imshow(Z, origin='lower', aspect='auto', cmap='viridis', extent=(1.5, 11.5, .5, 5.5))
for (ki, ni), v in np.ndenumerate(Z):
    if np.isfinite(v): h.text(Ns[ni], Ks[ki], f'{v:.0f}', ha='center', va='center', fontsize=7,
                              color='white' if v < np.nanmean(Z) else INK)
h.set(xlabel='occupied binding sites', ylabel='spacing between adjacent sites (nm)', xticks=Ns, yticks=Ks,
      yticklabels=[f'{6 * k}' for k in Ks])
h.spines[['top', 'right']].set_visible(True)
cb = fig.colorbar(im, ax=h, fraction=.05, pad=.02); cb.set_label('$T_m$ (°C)')
h.text(.99, -.2, 'blank: does not fit the 11-site edge', transform=h.transAxes, ha='right', fontsize=6.8, color=INK2)

# ---- per-panel standalone SVGs (for figure assembly / pseudo-figures), cropped to each axes' tight bbox ----
from matplotlib.transforms import Bbox
PANELS = os.path.join(HERE, 'panels'); os.makedirs(PANELS, exist_ok=True)
fig.canvas.draw(); _rend = fig.canvas.get_renderer()
def _save_panel(axes, name, padin=0.08):
    bb = Bbox.union([ax.get_tightbbox(_rend) for ax in axes]).transformed(fig.dpi_scale_trans.inverted())
    bb = Bbox.from_extents(bb.x0 - padin, bb.y0 - padin, bb.x1 + padin, bb.y1 + padin)
    fig.savefig(os.path.join(PANELS, name + '.svg'), bbox_inches=bb)
for _axes, _nm in (([a], 'A'), ([b], 'B'), ([c], 'C'), ([d], 'D'), ([e], 'E'), ([f_], 'F'), ([g], 'G'), ([h, cb.ax], 'H')):
    _save_panel(_axes, _nm)

fig.savefig(os.path.join(HERE, 'figure.svg'))
if '--png' in sys.argv: fig.savefig(sys.argv[sys.argv.index('--png') + 1], dpi=110)
plt.close(fig)
print(f"Figure 4: params {np.round(p, 4)}; CV RMS {main.cvRMS:.3f} C, rho {main.cvrho:.3f}; "
      f"theory Tm(7-nt, N=4 clustered) {D[(D.bond == '7-nt') & (D.n_sites == 4)].Tm_mean_C.iloc[0]:.1f} C")
