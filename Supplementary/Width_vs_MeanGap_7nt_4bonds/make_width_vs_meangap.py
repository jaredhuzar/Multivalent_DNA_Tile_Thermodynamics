"""
SI Fig. S10 (rebuilt): melt SHAPE vs mean gap for four-site 7-nt interfaces.

Replaces HS_vs_MeanGap_7nt_4bonds (apparent two-state dH/dS vs mean gap).  Those fitted dH/dS are strongly
correlated (rho ~ 0.999) and individually poorly determined, so they showed no trend.  In a two-state description
dS_vH = dH_vH/Tm + R ln(4/C_T) is fixed by dH_vH and Tm, so the only shape information beyond Tm is the width.

Panels
  A  width-derived van't Hoff enthalpy |dH_vH| = 1.5 R Tm^2 / dx (bimolecular; same estimator as Figure 2 C)
  B  cooperative fraction f_coop = |dH_vH| / |dH_NN|   (1 = whole interface melts as one unit)
  C  transition width dx: measured vs model predictions
       - multivalent model, C_eff,i = C_ref exp(-d_i/lambda) (d_i = mean distance of site i to the other bonded
         sites), intermediate readout (w = 0.57), 2 parameters fitted to the measured Tm of these interfaces only
       - same model with one uniform C_eff (1 parameter)
       - Fig 4 model A2: two-state, dH_eff = sum_k g_k dH_k, g_k = exp(-d_k/lambda_A2), lambda_A2 read from
         Supplementary/TableS1_Model_Comparison_All/phenomenological_parameter_CIs.csv, pinned at the measured Tm
     Widths of model curves: logistic-equivalent 1/(4 |dtheta/dT|) at theta = 1/2 (as in Figure 2 Panel H).

Data: same selection as Fig 3E -- Filtered_DF_HS Reliable rows, median per construct, four occupied sites, codes
without A/B/6/8/9 or '55' (four-site 7-nt designs, which carry the 5-nt reporter site in some designs).
Writes figure.svg + data.csv here.   Optional: --png PATH.
"""
import os, sys
import numpy as np, pandas as pd
from scipy.optimize import least_squares
from scipy.stats import spearmanr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
_hdr = os.path.join(REPO, 'Figure2_Updated', 'make_uniform_ceff_panel.py')
_src = open(_hdr).read(); _cut = '# ---- robust set'
assert _cut in _src, 'NN model header marker missing in make_uniform_ceff_panel.py'
_ns = {'__file__': _hdr}; exec(_src[:_src.index(_cut)], _ns)
R, C_TILE, duplexes = _ns['R'], _ns['C_TILE'], _ns['duplexes']
os.chdir(REPO)
W_INT, PITCH, TWIN = 0.57, 6.0, 1.5

# ---- data ----
HS = pd.read_csv('Filtered_DF_HS.csv', dtype={'Code': str}); HS = HS[HS.Reliable == True].copy()
for c in ('Tm', 'dx', 'Predicted_deltaH', 'Enthalpy', 'Entropy'): HS[c] = pd.to_numeric(HS[c], errors='coerce')
HS['code'] = HS.Code.str.partition('.')[0]
G = HS.groupby('code').agg(Tm=('Tm', 'median'), dx=('dx', 'median'), dH_NN=('Predicted_deltaH', 'median'),
                           dH_fit=('Enthalpy', 'median'), dS_fit=('Entropy', 'median'), n_reps=('Tm', 'size')).reset_index()
G = G[~G.code.str.contains('A|B|6|8|9|55')]
G = G[G.code.apply(lambda c: sum(ch != '0' for ch in c)) == 4].dropna(subset=['Tm', 'dx']).reset_index(drop=True)
sites = [[i for i, ch in enumerate(c) if ch != '0'] for c in G.code]
G['mean_gap_nm'] = [(s[-1] - s[0]) / (len(s) - 1) * PITCH for s in sites]
TmK = G.Tm.values + 273.15
G['dH_vH_kcal'] = 1.5 * R * TmK ** 2 / G.dx                      # magnitude
G['f_coop'] = G.dH_vH_kcal / G.dH_NN.abs()
n = len(G)

# ---- per-duplex arrays (8 duplexes = 4 sites x 2) ----
M = 8
HH = np.zeros((n, M)); SS = np.zeros((n, M)); RP = np.zeros((n, M), bool); DS = np.zeros((n, M)); DZ = np.zeros((n, M))
for i, (c, s) in enumerate(zip(G.code, sites)):
    d = duplexes(c)
    HH[i] = [x[0] for x in d]; SS[i] = [x[1] for x in d]; RP[i] = [x[2] for x in d]
    sd = [np.mean([abs(p - q) for q in s if q != p]) * PITCH for p in s]              # site mean distance
    DS[i] = np.repeat(sd, 2)
    z = np.array([p * PITCH + sgn * TWIN for p in s for sgn in (-1, 1)])                # A2 duplex geometry
    DZ[i] = [np.mean(np.abs(np.delete(z, k) - z[k])) for k in range(M)]
NREP = np.maximum(RP.sum(1), 1)

def theta(T, lc, lam, w=W_INT):
    C = 10 ** lc * np.exp(-DS / lam)
    x = np.exp(np.clip(-(HH - T[:, None] * SS) / (R * T[:, None]), -700, 700)) * C
    a = np.expm1(np.sum(np.log1p(x), 1)) / 10 ** lc * C_TILE
    f = np.where(a > 1e-12, ((2 * a + 1) - np.sqrt(4 * a + 1)) / (2 * np.maximum(a, 1e-300)), a)
    return f * (w + (1 - w) * np.sum(np.where(RP, x / (1 + x), 0), 1) / NREP)

def tm_pred(lc, lam):
    lo = np.full(n, 230.0); hi = np.full(n, 400.0)
    for _ in range(46):
        mid = (lo + hi) / 2; up = theta(mid, lc, lam) > 0.5
        lo = np.where(up, mid, lo); hi = np.where(up, hi, mid)
    return (lo + hi) / 2

def width_at(fun, TK, h=0.05):
    return 1.0 / (4.0 * np.abs((fun(TK + h) - fun(TK - h)) / (2 * h)))

fits = {}
r2 = least_squares(lambda p: tm_pred(p[0], p[1]) - TmK, [-2.0, 15.0], bounds=([-8, 1], [1, 1e4]), diff_step=1e-4)
r1 = least_squares(lambda p: tm_pred(p[0], np.inf) - TmK, [-2.5], bounds=([-8], [1]), diff_step=1e-4)
for key, lc, lam in (('dist', r2.x[0], r2.x[1]), ('uniform', r1.x[0], np.inf)):
    Tp = tm_pred(lc, lam)
    G[f'Tm_pred_{key}_C'] = Tp - 273.15
    G[f'dx_pred_{key}_C'] = width_at(lambda T: theta(T, lc, lam), Tp)
    fits[key] = dict(lc=lc, lam=lam, rms=float(np.sqrt(np.mean((Tp - TmK) ** 2))))
CI = pd.read_csv('Supplementary/TableS1_Model_Comparison_All/phenomenological_parameter_CIs.csv').set_index('parameter')
LAM_A2 = float(CI.loc['lambda_nm', 'point_estimate'])
dH_eff = np.sum(np.exp(-DZ / LAM_A2) * HH, 1)                                         # kcal/mol (negative)
def twostate(T, dH, Tm):
    K = (2.0 / C_TILE) * np.exp(-dH / R * (1 / T - 1 / Tm)); a = K * C_TILE
    return ((2 * a + 1) - np.sqrt(4 * a + 1)) / (2 * a)
G['dx_pred_A2_C'] = width_at(lambda T: twostate(T, dH_eff, TmK), TmK)

out = G[['code', 'n_reps', 'mean_gap_nm', 'Tm', 'dx', 'dH_vH_kcal', 'dH_NN', 'f_coop', 'dH_fit', 'dS_fit',
         'Tm_pred_dist_C', 'dx_pred_dist_C', 'Tm_pred_uniform_C', 'dx_pred_uniform_C', 'dx_pred_A2_C']].copy()
out.columns = ['code', 'n_reps', 'mean_gap_nm', 'Tm_C', 'width_obs_C', 'dH_vH_width_kcal', 'dH_NN_kcal', 'f_coop',
               'dH_twostate_fit_kcal', 'dS_twostate_fit_cal', 'Tm_pred_multivalent_dist_C', 'width_pred_multivalent_dist_C',
               'Tm_pred_multivalent_uniform_C', 'width_pred_multivalent_uniform_C', 'width_pred_A2_C']
out.round(4).to_csv(os.path.join(HERE, 'data.csv'), index=False)

# ---- stats ----
gap = G.mean_gap_nm.values
def sp(v): r, p = spearmanr(gap, v); return r, p
S = {k: sp(G[k].values) for k in ('Tm', 'dx', 'dH_vH_kcal', 'f_coop', 'dH_fit', 'dS_fit',
                                    'dx_pred_dist_C', 'dx_pred_uniform_C', 'dx_pred_A2_C')}
q1, q2 = np.percentile(gap, [100 / 3, 200 / 3]); lo, hi = G[gap <= q1], G[gap >= q2]

# ---- figure ----
INK, INK2, GRID, PT = "#0b0b0b", "#52514e", "#e6e5e0", "#3B6AA0"
COL = {'meas': INK, 'dist': '#1baf7a', 'uniform': '#8a8984', 'A2': '#e34948'}
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'svg.fonttype': 'none',
                     'axes.edgecolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2})
fig, ax = plt.subplots(1, 3, figsize=(13.2, 4.2))
ug = np.unique(np.round(gap, 2)); rng = np.random.default_rng(0)
def binned(v): return np.array([np.median(v[np.isclose(np.round(gap, 2), u)]) for u in ug])
def scatter(a, v):
    a.scatter(gap + rng.uniform(-.18, .18, n), v, s=14, color=PT, alpha=.55, edgecolor='none', zorder=2)
    a.plot(ug, binned(v), color=INK, lw=1.8, marker='o', ms=3.5, zorder=3, label='measured (median per gap)')
    a.grid(True, axis='y', color=GRID, lw=.6, zorder=0); a.set_xlabel('mean gap between bonded sites (nm)')
a = ax[0]; scatter(a, G.dH_vH_kcal.values)
a.set(ylabel=r'$|\Delta H_{vH}| = 1.5\,RT_m^2/\delta$  (kcal/mol)')
a.set_title(f'A   melt sharpens as bonds spread apart\nSpearman ρ = {S["dH_vH_kcal"][0]:+.2f}, P = {S["dH_vH_kcal"][1]:.0e}, n = {n}',
            loc='left', fontsize=9.3, weight='bold')
a = ax[1]; scatter(a, G.f_coop.values)
a.axhline(1, color=COL['A2'], ls='--', lw=1.2); a.text(gap.max(), .97, 'all-or-none', color=COL['A2'], ha='right', va='top', fontsize=7)
a.set(ylabel=r'cooperative fraction  $|\Delta H_{vH}|/|\Delta H_{NN}|$', ylim=(0, 1.08))
a.text(.03, .05, f'clustered third (≤{q1:.0f} nm): {lo.f_coop.median():.2f}\ndispersed third (≥{q2:.0f} nm): {hi.f_coop.median():.2f}',
       transform=a.transAxes, fontsize=7.2, color=INK2, va='bottom')
a.set_title(f'B   a larger fraction melts as one unit\nSpearman ρ = {S["f_coop"][0]:+.2f}', loc='left', fontsize=9.3, weight='bold')
a = ax[2]; scatter(a, G.dx.values)
for key, lab in (('dist', f'multivalent, $C_{{eff}}$ decays with distance (λ = {fits["dist"]["lam"]:.0f} nm)'),
                 ('uniform', 'multivalent, one uniform $C_{eff}$'),
                 ('A2', f'Fig. 4 model A2 (λ = {LAM_A2:.0f} nm)')):
    v = G[f'dx_pred_{key}_C'].values
    a.plot(ug, binned(v), color=COL[key], lw=1.8, ls='-' if key != 'uniform' else (0, (3, 2)), zorder=4,
           label=f'{lab}: ρ = {S[f"dx_pred_{key}_C"][0]:+.2f}')
a.set(ylabel='transition width δ (°C)', ylim=(0, 5.6))
a.legend(fontsize=6.6, frameon=True, framealpha=.92, edgecolor='none', loc='upper right')
a.set_title(f'C   width vs mean gap: measured ρ = {S["dx"][0]:+.2f}', loc='left', fontsize=9.3, weight='bold')
fig.tight_layout(); fig.savefig(os.path.join(HERE, 'figure.svg'))
if '--png' in sys.argv: fig.savefig(sys.argv[sys.argv.index('--png') + 1], dpi=170)
plt.close(fig)

print(f"S10 width vs mean gap: n={n}; " + "; ".join(f"{k} rho={v[0]:+.2f} (P={v[1]:.1g})" for k, v in S.items()))
print(f"  f_coop clustered third {lo.f_coop.median():.2f} -> dispersed third {hi.f_coop.median():.2f}; "
      f"Tm {lo.Tm.median():.1f} -> {hi.Tm.median():.1f} C; dx {lo.dx.median():.2f} -> {hi.dx.median():.2f} C")
print(f"  multivalent dist fit: log10 C_ref={fits['dist']['lc']:.2f}, lambda={fits['dist']['lam']:.1f} nm, Tm RMS={fits['dist']['rms']:.2f} C; "
      f"uniform: log10 C={fits['uniform']['lc']:.2f}, RMS={fits['uniform']['rms']:.2f} C; A2 lambda={LAM_A2:.1f} nm")

# ---- standalone single panels (four-site 7-nt): width-derived van't Hoff enthalpy and transition width vs mean gap ----
PAN = os.path.join(HERE, 'panels'); os.makedirs(PAN, exist_ok=True)
for key, col, ylab, ttl in (
        ('dHvH_vs_gap', 'dH_vH_kcal', r"van't Hoff enthalpy $|\Delta H_{vH}| = 1.5\,RT_m^2/\delta$ (kcal/mol)",
         "van't Hoff enthalpy vs bond separation"),
        ('width_vs_gap', 'dx', 'transition width δ (°C)', 'Transition width vs bond separation')):
    v = G[col].values; r_, p_ = spearmanr(gap, v)
    fig, a = plt.subplots(figsize=(4.8, 3.9)); rng = np.random.default_rng(0)
    a.scatter(gap + rng.uniform(-.18, .18, n), v, s=15, color=PT, alpha=.55, edgecolor='none', zorder=2)
    a.plot(ug, binned(v), color=INK, lw=1.8, marker='o', ms=4, zorder=3, label='median per mean gap')
    a.set(xlabel='mean gap between bonded sites (nm)', ylabel=ylab, ylim=(None, v.max() + .22 * (v.max() - v.min())))
    a.text(.03, .97, f'Spearman ρ = {r_:+.2f}, P = {p_:.0e}\nn = {n} four-site 7-nt interfaces', transform=a.transAxes,
           va='top', fontsize=8)
    a.grid(True, axis='y', color=GRID, lw=.6, zorder=0); a.legend(fontsize=7.5, frameon=False, loc='lower right')
    a.set_title(ttl, loc='left', fontsize=9.5, weight='bold')
    fig.tight_layout(); fig.savefig(os.path.join(PAN, key + '.svg'))
    if '--png' in sys.argv: fig.savefig(os.path.join(PAN, key + '.png'), dpi=150)
    plt.close(fig)
    G[['code', 'mean_gap_nm', 'Tm', 'dx', 'dH_vH_kcal']].rename(columns={'Tm': 'Tm_C', 'dx': 'width_C'}).round(4) \
      .to_csv(os.path.join(PAN, key + '.csv'), index=False)
print("  standalone panels -> panels/dHvH_vs_gap.svg, panels/width_vs_gap.svg (+ csv)")

# ---- standalone: two-state model A2 (distance weight on both dH and dS) vs measured transition width ----
wA = G.dx_pred_A2_C.values; wM = G.dx.values; ratio = wA / wM
fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.9))
a = ax[0]; lim = [0, max(wM.max(), wA.max()) * 1.08]
a.scatter(wM, wA, s=15, color=COL['A2'], alpha=.6, edgecolor='none', zorder=2)
a.plot(lim, lim, color=INK2, lw=1, zorder=1); a.text(lim[1] * .97, lim[1] * .93, 'predicted = measured', ha='right', fontsize=7, color=INK2)
a.plot(lim, [0, lim[1] / 2], color=INK2, lw=1, ls=(0, (3, 2)), zorder=1); a.text(lim[1] * .97, lim[1] * .43, 'half as wide', ha='right', fontsize=7, color=INK2)
a.set(xlim=lim, ylim=lim, xlabel='measured transition width δ (°C)', ylabel='two-state model width (°C)')
a.text(.03, .97, f'median predicted/measured = {np.median(ratio):.2f}\n(IQR {np.percentile(ratio, 25):.2f}–{np.percentile(ratio, 75):.2f}), n = {n}',
       transform=a.transAxes, va='top', fontsize=8)
a.grid(True, color=GRID, lw=.6, zorder=0)
a.set_title('A   two-state melts are about half as wide', loc='left', fontsize=9.5, weight='bold')
a = ax[1]; rng = np.random.default_rng(0)
a.scatter(gap + rng.uniform(-.18, .18, n), wM, s=12, color=PT, alpha=.45, edgecolor='none', zorder=2)
a.scatter(gap + rng.uniform(-.18, .18, n), wA, s=12, color=COL['A2'], alpha=.45, edgecolor='none', zorder=2)
a.plot(ug, binned(wM), color=INK, lw=1.8, marker='o', ms=4, zorder=3, label=f'measured  (ρ = {S["dx"][0]:+.2f})')
a.plot(ug, binned(wA), color=COL['A2'], lw=1.8, marker='s', ms=4, zorder=3, label=f'two-state model  (ρ = {S["dx_pred_A2_C"][0]:+.2f})')
a.set(xlabel='mean gap between bonded sites (nm)', ylabel='transition width δ (°C)', ylim=(0, wM.max() * 1.1))
a.legend(fontsize=7.5, frameon=False, loc='upper right'); a.grid(True, axis='y', color=GRID, lw=.6, zorder=0)
a.set_title('B   model broadens with spacing; data narrow', loc='left', fontsize=9.5, weight='bold')
fig.tight_layout(); fig.savefig(os.path.join(PAN, 'twostate_width_vs_measured.svg'))
if '--png' in sys.argv: fig.savefig(os.path.join(PAN, 'twostate_width_vs_measured.png'), dpi=150)
plt.close(fig)
G[['code', 'mean_gap_nm', 'Tm', 'dx', 'dx_pred_A2_C']].assign(ratio_pred_over_meas=ratio) \
  .rename(columns={'Tm': 'Tm_C', 'dx': 'width_measured_C', 'dx_pred_A2_C': 'width_twostate_A2_C'}).round(4) \
  .to_csv(os.path.join(PAN, 'twostate_width_vs_measured.csv'), index=False)
print(f"  two-state A2 width: median ratio {np.median(ratio):.2f}; rho(gap, A2 width) {S['dx_pred_A2_C'][0]:+.2f} vs measured {S['dx'][0]:+.2f}")

# ---- merged S11: van't Hoff enthalpy (A) and transition width (B) vs mean gap in one figure ----
fig, axs = plt.subplots(1, 2, figsize=(9.6, 3.9)); rng = np.random.default_rng(0)
for a, (lab, col, ylab, ttl) in zip(axs, (
        ('A', 'dH_vH_kcal', r"$|\Delta H_{vH}| = 1.5\,RT_m^2/\delta$ (kcal/mol)", "van't Hoff enthalpy vs bond separation"),
        ('B', 'dx', 'transition width δ (°C)', 'Transition width vs bond separation'))):
    v = G[col].values; r_, p_ = spearmanr(gap, v)
    a.scatter(gap + rng.uniform(-.18, .18, n), v, s=14, color=PT, alpha=.55, edgecolor='none', zorder=2)
    a.plot(ug, binned(v), color=INK, lw=1.8, marker='o', ms=4, zorder=3, label='median per mean gap')
    a.set(xlabel='mean gap between bonded sites (nm)', ylabel=ylab, ylim=(None, v.max() + .22 * (v.max() - v.min())))
    a.text(.03, .97, f'Spearman ρ = {r_:+.2f}, P = {p_:.0e}', transform=a.transAxes, va='top', fontsize=8)
    a.grid(True, axis='y', color=GRID, lw=.6, zorder=0); a.legend(fontsize=7.5, frameon=False, loc='lower right')
    a.set_title(f'{lab}   {ttl}', loc='left', fontsize=9.5, weight='bold')
fig.suptitle(f'Four-site 7-nt interfaces (n = {n})', fontsize=9.5, y=1.0); fig.tight_layout()
fig.savefig(os.path.join(PAN, 'dHvH_and_width_vs_gap.svg')); plt.close(fig)
G[['code', 'mean_gap_nm', 'Tm', 'dx', 'dH_vH_kcal']].rename(columns={'Tm': 'Tm_C', 'dx': 'width_C'}).round(4) \
  .to_csv(os.path.join(PAN, 'dHvH_and_width_vs_gap.csv'), index=False)
