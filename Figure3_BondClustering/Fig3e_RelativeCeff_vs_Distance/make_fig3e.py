"""
Fig 3E (rebuilt): relative effective concentration of a bond vs the distance to its nearest bonded neighbour,
fitted to MEASURED Tm -- no van't Hoff dG.

Replaces Fig3d_RelativeOccupancy_vs_Distance (notebook), which decomposed the van't Hoff dG_exp over distance
bins.  That dG is the free energy of the small cooperative unit that controls the transition (Figure 2D), is
noise-biased, and compresses per-bond differences, so the decay it implied (0.31 at 30 nm) was too shallow.

MODEL (same multivalent model as Figure2_Updated Panel H, plus one distance factor per bin):
    C_i  = C_ref * phi[bin(i)],   bin = nearest occupied site 1,2,3,4,>=5 sites away (6 nm pitch), phi[6 nm] = 1
    K_a  = (1/C_ref) * [ prod_i (1 + K_i C_i) - 1 ]      K_i: per-duplex NN (SantaLucia 1998, salt-corrected)
    mass balance at 4.5 nM per tile;  readout theta = f * (w + (1-w) * p_reporter)
Five global parameters (log10 C_ref, log10 phi[12..30 nm]) are fitted by least squares to the median Tm of
every 4-site 7nt interface (same 109 constructs as Fig 3C).  phi_d = C_eff(d)/C_eff(6 nm) is the relative
effective association constant of a bond -- the quantity the old exp(-(beta_d - beta_d1)/RT) estimated; it
equals the relative occupancy when bonds are far from saturation.
This figure uses the DOCKING readout (w = 1) for BOTH curves it draws; the intermediate (w = 0.57) and
reporter (w = 0) fits are still computed and written to data.csv for robustness.  It plots two quantities
vs distance: (i) relative effective concentration phi_d = C_eff(d)/C_eff(6 nm), the geometric/tethering
factor (no saturation), and (ii) relative occupancy theta(d)/theta(6 nm) with theta_i = K_i C_eff,i /
(1 + K_i C_eff,i) evaluated at the MEDIAN Tm of the four-site 7nt interfaces shown -- i.e. C_eff passed
through Langmuir saturation, so it is compressed toward 1 and decays less steeply than phi_d.  The two
coincide only far from saturation.
95% CI: Gauss-Newton covariance (t, df = n - 5); constructs bootstrap (200x) for the intermediate readout;
occupancy CI maps the docking phi_d CI through the saturation at the median Tm.
Also fitted: an exponential phi = exp(-(d - 6 nm)/lambda) through the bin estimates (log-space least squares).

Writes figure.svg + data.csv here and Figure_Source_Data/Fig3E_RelativeCeff_vs_Distance.csv.
Run after the Filtered_DF tables exist; reuses the NN model header of Figure2_Updated/make_uniform_ceff_panel.py.
Optional: --png PATH  (preview).
"""
import os, sys
import numpy as np, pandas as pd
from scipy.optimize import least_squares
from scipy.stats import t as tdist
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
_hdr_file = os.path.join(REPO, 'Figure2_Updated', 'make_uniform_ceff_panel.py')
_src = open(_hdr_file).read(); _cut = '# ---- robust set'
assert _cut in _src, 'NN model header marker missing in make_uniform_ceff_panel.py'
_ns = {'__file__': _hdr_file}; exec(_src[:_src.index(_cut)], _ns)          # NN table, duplexes(), constants
R, C_TILE, duplexes = _ns['R'], _ns['C_TILE'], _ns['duplexes']
os.chdir(REPO)

# ---- data: 4-site 7nt interfaces (same selection as Fig 3C / the old 3E), median Tm over Reliable wells ----
HS = pd.read_csv('Filtered_DF_HS.csv', dtype={'Code': str}); HS = HS[HS.Reliable == True]
G = HS.groupby('Code').agg(Tm=('Tm', 'median'), n_wells=('Tm', 'size')).reset_index().rename(columns={'Code': 'code'})
G = G[~G.code.str.contains('A|B|6|8|9|55')]
G = G[G.code.apply(lambda c: sum(ch != '0' for ch in c)) == 4].reset_index(drop=True)
NM, NB = 6.0, 5                                         # site pitch (nm); bins 1..4 and >=5 sites

def site_bins(code):
    pos = [i for i, ch in enumerate(code) if ch != '0']
    return {p: min(min(abs(p - q) for q in pos if q != p), NB) for p in pos}

n = len(G); M = 8                                       # 4 sites x 2 duplexes
HH = np.zeros((n, M)); SS = np.zeros((n, M)); RP = np.zeros((n, M), bool); BN = np.zeros((n, M), int)
for i, c in enumerate(G.code):
    d = duplexes(c); b = site_bins(c)
    HH[i] = [x[0] for x in d]; SS[i] = [x[1] for x in d]; RP[i] = [x[2] for x in d]
    BN[i] = [b[p] for p, ch in enumerate(c) if ch != '0' for _ in (0, 1)]   # both duplexes of a site share its bin
NREP = RP.sum(1); Tobs = G.Tm.values

def theta(T, lcr, lphi, w):
    C = 10 ** (lcr + np.concatenate([[0.0], lphi])[BN - 1])
    x = np.exp(np.clip(-(HH - T[:, None] * SS) / (R * T[:, None]), -700, 700)) * C
    a = np.expm1(np.sum(np.log1p(x), axis=1)) / 10 ** lcr * C_TILE
    f = np.where(a > 1e-12, ((2 * a + 1) - np.sqrt(4 * a + 1)) / (2 * np.maximum(a, 1e-300)), a)
    return f * (w + (1 - w) * np.sum(np.where(RP, x / (1 + x), 0), 1) / NREP)

def tm_pred(p, w):                                      # vectorized bisection (theta falls with T)
    lo = np.full(n, 230.0); hi = np.full(n, 400.0)
    for _ in range(44):
        mid = (lo + hi) / 2; up = theta(mid, p[0], p[1:], w) > 0.5
        lo = np.where(up, mid, lo); hi = np.where(up, hi, mid)
    return (lo + hi) / 2 - 273.15

def fit(w, idx=None, nphi=NB - 1):
    idx = np.arange(n) if idx is None else idx
    full = lambda p: np.r_[p, np.zeros(NB - 1 - nphi)]
    return least_squares(lambda p: (tm_pred(full(p), w) - Tobs)[idx], x0=[-2.5] + [0.0] * nphi, diff_step=1e-3)

READ = [('intermediate', 0.57, '#1baf7a'), ('docking', 1.0, '#6250d6'), ('reporter', 0.0, '#e34948')]
d_nm = np.arange(1, NB + 1) * NM
rows, fits = [], {}
for name, w, _ in READ:
    r = fit(w); J = r.jac; dof = n - J.shape[1]
    cov = (r.fun @ r.fun) / dof * np.linalg.inv(J.T @ J); tc = tdist.ppf(.975, dof)
    boot = None
    if name == 'intermediate':
        rng = np.random.default_rng(1)
        boot = np.array([fit(w, rng.integers(0, n, n)).x[1:] for _ in range(200)])
    lphi = np.r_[0.0, r.x[1:]]; se = np.r_[0.0, np.sqrt(np.diag(cov))[1:]]
    xd = d_nm[1:] - NM; lam = (xd @ xd) / (-(xd @ (lphi[1:] * np.log(10))))     # exp through phi(6 nm) = 1
    fits[name] = dict(lphi=lphi, se=se, tc=tc, lam=lam, rms=np.sqrt(np.mean(r.fun ** 2)), lcr=r.x[0])
    for k in range(NB):
        bl, bh = ((10 ** np.percentile(boot[:, k - 1], [2.5, 97.5])) if (boot is not None and k) else
                  ((1.0, 1.0) if k == 0 else (np.nan, np.nan)))
        rows.append(dict(readout=name, w=w, nearest_bond_distance_nm=d_nm[k], rel_Ceff=round(10 ** lphi[k], 4),
                         ci_lo=round(10 ** (lphi[k] - tc * se[k]), 4), ci_hi=round(10 ** (lphi[k] + tc * se[k]), 4),
                         boot_ci_lo=round(bl, 4), boot_ci_hi=round(bh, 4),
                         n_sites_in_bin=int((BN == k + 1).sum() // 2), lambda_exp_nm=round(lam, 2),
                         log10_Cref_M=round(r.x[0], 3), Tm_rms_C=round(fits[name]['rms'], 3)))
r0 = fit(0.57, nphi=0); rms0 = np.sqrt(np.mean(r0.fun ** 2))

# ---- relative occupancy vs distance (DOCKING fit), evaluated at the median Tm of the interfaces shown ----
# occupancy theta_i = K_i C_eff,i / (1 + K_i C_eff,i) is C_eff passed through Langmuir saturation: it equals
# relative C_eff only far from saturation and is compressed toward 1 near it.  K_i = exp(-(dH - T dS)/RT) at
# T = median interface Tm; C_eff,i = 10**(lcr + lphi[bin]).  Both C_eff and occupancy use the docking (w=1) fit.
fd = fits['docking']; lcr_d, lphi_d, se_d, tc_d = fd['lcr'], fd['lphi'], fd['se'], fd['tc']
Tmed_C = float(np.median(Tobs)); Tmed_K = Tmed_C + 273.15
# representative per-bond K at the median Tm (median over the 7nt duplexes in the set).  K is a property of the
# duplex sequence, NOT of distance (Fig 3C moves a CONSTANT-sequence bond), so we hold it fixed and let only
# C_eff vary with distance -- otherwise the differing sequence composition of each distance bin leaks in.
Kbar = float(np.median(np.exp(np.clip(-(HH - Tmed_K * SS) / (R * Tmed_K), -700, 700))))
def _rel_occ(rel_by_bin):                                                      # rel_by_bin: C_eff(bin)/C_eff(6 nm)
    x = Kbar * 10 ** lcr_d * np.asarray(rel_by_bin); occ = x / (1 + x)         # Langmuir saturation of C_eff(bin)
    return occ / occ[0]
occ_pt = _rel_occ(10 ** lphi_d)                                                # point estimate
occ_lo = _rel_occ(10 ** (lphi_d - tc_d * se_d))                                # phi_d CI mapped through saturation
occ_hi = _rel_occ(10 ** (lphi_d + tc_d * se_d))

D = pd.DataFrame(rows)
D['eval_Tm_C'] = np.nan; D['rel_occupancy_medianTm'] = np.nan; D['rel_occ_ci_lo'] = np.nan; D['rel_occ_ci_hi'] = np.nan
_dm = D.readout == 'docking'; D.loc[_dm, 'eval_Tm_C'] = round(Tmed_C, 2)
for k in range(NB):
    _m = _dm & (D.nearest_bond_distance_nm == d_nm[k])
    D.loc[_m, 'rel_occupancy_medianTm'] = round(occ_pt[k], 4)
    D.loc[_m, 'rel_occ_ci_lo'] = round(occ_lo[k], 4); D.loc[_m, 'rel_occ_ci_hi'] = round(occ_hi[k], 4)
D.to_csv(os.path.join(HERE, 'data.csv'), index=False)
os.makedirs('Figure_Source_Data', exist_ok=True)
D.to_csv(os.path.join('Figure_Source_Data', 'Fig3E_RelativeCeff_vs_Distance.csv'), index=False)

# ---- figure ----
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e6e5e0"
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'svg.fonttype': 'none',
                     'axes.edgecolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2})
fig, ax = plt.subplots(figsize=(4.8, 4.0))
DOCK, OCC = '#6250d6', '#e6820a'; dd = np.linspace(6, 31, 200)
# (i) docking relative effective concentration (geometric; no saturation) + its exponential length scale
ax.plot(dd, np.exp(-(dd - NM) / fd['lam']), color=DOCK, lw=1, ls=(0, (3, 2)), zorder=1,
        label=f'exp(−(d − 6 nm)/λ), λ = {fd["lam"]:.0f} nm')
qd = D[D.readout == 'docking'].sort_values('nearest_bond_distance_nm'); yC = qd.rel_Ceff.values
ax.errorbar(qd.nearest_bond_distance_nm, yC, yerr=[yC - qd.ci_lo, qd.ci_hi - yC], fmt='o', ms=7,
            color=DOCK, mfc=DOCK, mew=1.5, capsize=3, lw=1.3, zorder=3, label='relative $C_{eff}$ (docking)')
# (ii) relative occupancy at the median interface Tm = C_eff after Langmuir saturation (flatter than C_eff)
ax.errorbar(d_nm + 0.7, occ_pt, yerr=[occ_pt - occ_lo, occ_hi - occ_pt], fmt='s', ms=6,
            color=OCC, mfc='white', mew=1.5, capsize=3, lw=1.3, zorder=3,
            label=f'relative occupancy at median $T_m$ = {Tmed_C:.0f} °C (docking)')
ax.axhline(1, color=GRID, lw=1, zorder=0)
ax.set(xlim=(3.5, 32.5), ylim=(0, 1.15), xticks=d_nm, xticklabels=['6', '12', '18', '24', '≥30'],
       xlabel='distance to nearest bonded site (nm)', ylabel='value relative to a clustered (6 nm) bond')
ax.legend(fontsize=7, frameon=False, loc='lower left')
ax.text(31.8, .99, f'fit to $T_m$ of {n} four-site\n7-nt interfaces; RMS {fd["rms"]:.1f} °C\n'
        f'(one uniform $C_{{eff}}$: {rms0:.1f} °C)', ha='right', va='top', fontsize=6.8, color=INK2)
ax.set_title('E   isolated bonds: lower $C_{eff}$, less time engaged at $T_m$', loc='left', fontsize=9.5, weight='bold')
ax.grid(True, axis='y', color=GRID, lw=.6, zorder=0)
fig.tight_layout(); fig.savefig(os.path.join(HERE, 'figure.svg'))
if '--png' in sys.argv: fig.savefig(sys.argv[sys.argv.index('--png') + 1], dpi=200)
plt.close(fig)
print(f"Fig 3E (Tm model): n={n}; " + "; ".join(
    f"{k}: phi=" + ",".join(f"{10 ** v:.2f}" for v in f['lphi']) + f" lambda={f['lam']:.1f} nm RMS={f['rms']:.2f} C"
    for k, f in fits.items()) + f"; uniform-C_eff RMS={rms0:.2f} C")
