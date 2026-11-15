import os as _os
def _find_repo():
    d = _os.path.dirname(_os.path.abspath(__file__))
    while d != '/' and not _os.path.exists(_os.path.join(d, 'Filtered_DF_Tm.csv')):
        d = _os.path.dirname(d)
    return d
REPO = _find_repo()
"""
Figure 2 (updated) -- Panel H: a single, uniform effective concentration cannot reproduce the melt.

MODEL (no all-or-none assumption).  Once the two tiles dock, each sticky-end duplex i binds independently
at a local concentration C_eff.  Summing every bound microstate (one bimolecular docking):

    K_d(T) = (1/C_eff) * [ prod_i (1 + K_i(T) * C_eff) - 1 ]
    K_i(T) = exp(-(dH_i - T*dS_i) / RT)          per-duplex SantaLucia-1998 NN, dS salt-corrected
                                                  (0.368*(L-1)*ln[Na+]eq, [Na+]eq = 3.79*sqrt(0.0113 M))
    mass balance (each tile at c = 4.5 nM):  K_d = f / ((1-f)^2 c)   ->  f = ((2a+1) - sqrt(4a+1)) / (2a), a = K_d*c

  The summed-shortfall C_eff formula is the all-or-none limit of this (keep only the fully-bound term).

READOUT (two limits of what the FRET signal reports):
    docking   theta = f                                  (signal tracks tile proximity)
    reporter  theta = f * mean_r[ K_r C_eff / (1 + K_r C_eff) ]   (the two duplexes at the central site)

TEST.  For each construct the ONE free parameter, C_eff (taken temperature-independent), is fitted so that
theta(Tm_measured) = 1/2.  The predicted transition width is then a genuine prediction, compared with the
measured Boltzmann width dx.  Reference: all-or-none bimolecular two-state with dH = sum NN, pinned at Tm.
Width of any model curve = logistic-equivalent width 1/(4 |dtheta/dT|) at theta = 1/2.

DATA.  The figure's robust set (same as make_updated_figure2.py): Filtered_DF_HS rows with Reliable==True &
Confidence=='high', median per construct, >= 2 replicates, 9-nt excluded (A/B designs excluded: no per-
position NN sequence model).  Panel (left) shows one representative 7nt 4-site construct: its raw replicate
curves normalized with the pipeline's droop-trimmed Boltzmann fit (dx >= 0.3 C), the measured fit, and the
three models; (right) predicted vs observed width for every construct.

PANEL I (readout weight).  Real FRET readout is likely between the two limits: when the central reporter
duplexes fray but the tiles stay docked, the dye/quencher are still held within a few nm and stay partly
quenched.  theta = f * (w + (1-w) * p_reporter), w in [0,1] (w=0 reporter, w=1 docking).  For each w, C_eff
is refit to the measured Tm; w_best is where predicted width = observed width.  A single w (~0.57) that is
flat across bond length and valency reproduces the widths, so the width test alone CANNOT show that bonds
are non-equivalent -- it shows the pure readouts bracket the data and all-or-none is excluded.

Writes Figure2_Updated/PanelH_UniformCeff_Width/{figure.svg, data.csv, example_curves.csv},
Figure2_Updated/PanelI_ReadoutWeight/{figure.svg, data.csv, width_ratio_vs_w.csv}, and
Figure_Source_Data/Fig2U_H_UniformCeff_Width.csv + Fig2U_I_ReadoutWeight.csv.
Run from anywhere (after the Filtered_DF tables exist).   Optional:  --png PATH   (PNG preview of H;
I is written next to it with suffix _I).
"""
import os, sys, math, numpy as np, pandas as pd
from scipy.optimize import brentq, curve_fit
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
os.chdir(REPO)

OUTDIR = os.path.join("Figure2_Updated", "PanelH_UniformCeff_Width"); os.makedirs(OUTDIR, exist_ok=True)
SRC = "Figure_Source_Data"; os.makedirs(SRC, exist_ok=True)
R = 1.987204e-3                      # kcal/(mol K)
C_TILE = 4.5e-9                      # each tile (M); total strands C_T = 9 nM
LN_NAEQ = math.log(3.79 * math.sqrt(0.0113))
# model colours chosen NOT to reuse Figure 2's bond-length hues (blue 7nt / orange 5nt / green 6nt);
# identity is also carried by line style (lines) and labels, so it is never colour-alone.
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e6e5e0"
MODEL = {'aon': ("#8a8984", (0, (4, 2)), 'all-or-none, ΔH = ΣNN'),
         'docking': ("#6250d6", '-', 'uniform $C_{eff}$, docking readout'),
         'reporter': ("#e34948", (0, (1.2, 1.2)), 'uniform $C_{eff}$, reporter readout')}
MARK = {7: 'o', 5: 's', 6: '^'}

# ---- SantaLucia 1998 unified NN (1 M NaCl), orientation-correct lookup (same table as Tile_Tile_Master) ----
NN_H = {"AA/TT": -7.9, "AT/TA": -7.2, "TA/AT": -7.2, "CA/GT": -8.5, "GT/CA": -8.4, "CT/GA": -7.8,
        "GA/CT": -8.2, "CG/GC": -10.6, "GC/CG": -9.8, "GG/CC": -8.0}
NN_S = {"AA/TT": -22.2, "AT/TA": -20.4, "TA/AT": -21.3, "CA/GT": -22.7, "GT/CA": -22.4, "CT/GA": -21.0,
        "GA/CT": -22.2, "CG/GC": -27.2, "GC/CG": -24.4, "GG/CC": -19.9}
INIT_H = {"AT": 2.3, "GC": 0.1}; INIT_S = {"AT": 4.1, "GC": -2.8}
_COMP = str.maketrans("ATCG", "TAGC")
def nn_thermo(seq):
    """(dH kcal/mol, dS kcal/mol/K) at 1 M NaCl for a non-self-complementary duplex."""
    seq = seq.upper(); H = S = 0.0
    for i in range(len(seq) - 1):
        top = seq[i:i + 2]; bot = top.translate(_COMP); k = f"{top}/{bot}"
        if k not in NN_H: k = f"{bot[::-1]}/{top[::-1]}"          # read from the other strand
        H += NN_H[k]; S += NN_S[k]
    for b in (seq[0], seq[-1]):
        t = "AT" if b in "AT" else "GC"; H += INIT_H[t]; S += INIT_S[t]
    return H, S / 1000.0

SEQ = {t: pd.read_csv(f"nt{t}.csv", index_col=0)["Sequence"].astype(str).tolist() for t in ("5", "6", "7")}
def duplexes(code):
    """per duplex: (dH, dS_salt, is_reporter) -- two duplexes per occupied site; reporter = central site (5)."""
    out = []
    for i, ch in enumerate(code):
        if ch == '0': continue
        t = '5' if ch == 'X' else ch
        for s in (SEQ[t][2 * i], SEQ[t][2 * i + 1]):
            H, S = nn_thermo(s); L = len(s)
            out.append((H, S + 0.368e-3 * (L - 1) * LN_NAEQ, i == 5))
    return out

def theta(TK, dup, logC, mode):
    C = 10.0 ** logC
    H = np.array([d[0] for d in dup]); S = np.array([d[1] for d in dup]); rep = np.array([d[2] for d in dup])
    x = np.exp(-(H[None, :] - TK[:, None] * S[None, :]) / (R * TK[:, None])) * C
    Kd = np.expm1(np.sum(np.log1p(x), axis=1)) / C
    a = Kd * C_TILE; f = np.where(a > 1e-12, ((2 * a + 1) - np.sqrt(4 * a + 1)) / (2 * a), a)
    return f if mode == 'docking' else f * (x[:, rep] / (1 + x[:, rep])).mean(axis=1)

def all_or_none(TK, dH, TmK):
    K = (2.0 / C_TILE) * np.exp(-dH / R * (1 / TK - 1 / TmK))     # K(Tm) gives f = 1/2
    a = K * C_TILE; return ((2 * a + 1) - np.sqrt(4 * a + 1)) / (2 * a)

def width(TK, th):
    i = int(np.argmin(np.abs(th - 0.5))); return float(1.0 / (4.0 * -np.gradient(th, TK)[i]))

# ---- robust set (identical definition to make_updated_figure2.load) ----
d = pd.read_csv('Filtered_DF_HS.csv')
d = d[(d['Reliable'] == True) & (d['Confidence'] == 'high')].copy()
for c in ['Tm', 'dx', 'Predicted_deltaH', 'Num_Bonds', 'NucleotideBonds']:
    d[c] = pd.to_numeric(d[c], errors='coerce')
d['code'] = d['Code'].astype(str).str.partition('.')[0]
g = d.groupby('code').agg(Tm=('Tm', 'median'), dx=('dx', 'median'), dHNN=('Predicted_deltaH', 'median'),
                          N=('Num_Bonds', 'first'), nt=('NucleotideBonds', 'first'), nrep=('Tm', 'size')).dropna()
g = g[(g.nrep >= 2) & (g.nt != 9)].copy(); g['nt'] = g['nt'].astype(int)

TK = np.arange(0.0, 90.0, 0.02) + 273.15
rows = []
for code, r in g.iterrows():
    dup = duplexes(code); TmK = r.Tm + 273.15
    rec = dict(code=code, nt=int(r.nt), n_sites=int(r.N), n_reps=int(r.nrep), Tm_C=r.Tm, width_obs_C=r.dx,
               width_allornone_C=width(TK, all_or_none(TK, -abs(r.dHNN), TmK)))
    for mode in ('docking', 'reporter'):
        lc = brentq(lambda v: theta(np.array([TmK]), dup, v, mode)[0] - 0.5, -14, 4)
        rec[f'log10_Ceff_{mode}_M'] = lc
        rec[f'log10_Ceff_over_Ctile_{mode}'] = lc - math.log10(C_TILE)
        rec[f'width_uniformCeff_{mode}_C'] = width(TK, theta(TK, dup, lc, mode))
    rows.append(rec)
D = pd.DataFrame(rows).sort_values(['nt', 'n_sites', 'code']).reset_index(drop=True)
D = D.round({'Tm_C': 3, 'width_obs_C': 3, 'width_allornone_C': 3, 'width_uniformCeff_docking_C': 3,
             'width_uniformCeff_reporter_C': 3, 'log10_Ceff_docking_M': 3, 'log10_Ceff_reporter_M': 3,
             'log10_Ceff_over_Ctile_docking': 3, 'log10_Ceff_over_Ctile_reporter': 3})
D.to_csv(os.path.join(OUTDIR, 'data.csv'), index=False)
D.to_csv(os.path.join(SRC, 'Fig2U_H_UniformCeff_Width.csv'), index=False)

# ---- representative construct: 7nt, 4 sites, >=3 replicates, width closest to the 7nt median ----
cand = D[(D.nt == 7) & (D.n_sites == 4) & (D.n_reps >= 3)].copy()
cand['_d'] = (cand.width_obs_C - D[D.nt == 7].width_obs_C.median()).abs()
ex = cand.sort_values(['_d', 'code']).iloc[0]; code = ex.code
def boltz(x, A1, A2, x0, dx): return A1 + (A2 - A1) / (1 + np.exp((x0 - x) / dx))
BND = ([-np.inf, -np.inf, -np.inf, 0.3], [np.inf] * 4)
pts = []
for exp in sorted(d[d.code == code].Experiment.unique()):
    raw = pd.read_csv(os.path.join('Melting_Results', exp), sep='\t' if exp.endswith('.txt') else ',', encoding='latin1')
    T = pd.to_numeric(raw.iloc[:, 1], errors='coerce').values
    for k in [k for k in raw.columns[2:] if str(k).partition('.')[0].replace('X', '5') == code]:
        y = pd.to_numeric(raw[k], errors='coerce').values; ok = np.isfinite(T) & np.isfinite(y)
        t, y = T[ok], y[ok]; o = np.argsort(t); t, y = t[o], y[o]
        p0 = [np.percentile(y, 5), np.percentile(y, 95), np.median(t), 3.0]
        p, _ = curve_fit(boltz, t, y, p0=p0, maxfev=20000, bounds=BND)
        sm = pd.Series(y).rolling(7, center=True, min_periods=1).median().values
        post = t > p[2] + 1.5 * p[3]
        onset = t[np.where(post)[0][np.argmax(sm[post])]] if post.sum() >= 3 else t.max()
        keep = t <= min(onset + 3.0, t.max())
        p, _ = curve_fit(boltz, t[keep], y[keep], p0=p0, maxfev=20000, bounds=BND)
        pts.append(pd.DataFrame({'well': k, 'plate': exp, 'T_C': t[keep], 'theta': (p[1] - y[keep]) / (p[1] - p[0])}))
P = pd.concat(pts, ignore_index=True)
TmK = ex.Tm_C + 273.15; dup = duplexes(code); Tp = np.arange(20.0, 65.0, 0.05) + 273.15
curves = {'measured_fit': 1 - 1 / (1 + np.exp((TmK - Tp) / ex.width_obs_C)),
          'aon': all_or_none(Tp, -abs(g.loc[code, 'dHNN']), TmK),
          'docking': theta(Tp, dup, ex.log10_Ceff_docking_M, 'docking'),
          'reporter': theta(Tp, dup, ex.log10_Ceff_reporter_M, 'reporter')}
pd.concat([P.assign(kind='raw_replicate'),
           pd.DataFrame({'T_C': Tp - 273.15, **{f'theta_{k}': v for k, v in curves.items()}}).assign(kind='model_curves')],
          ignore_index=True).to_csv(os.path.join(OUTDIR, 'example_curves.csv'), index=False)

# ---- figure ----
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'svg.fonttype': 'none',
                     'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2})
fig, ax = plt.subplots(1, 2, figsize=(9.8, 4.1), gridspec_kw={'width_ratios': [1.15, 1]})
a = ax[0]
a.scatter(P.T_C, P.theta, s=7, color=INK2, alpha=.35, edgecolor='none', zorder=1, label=f'measured ({P.well.nunique()} replicates)')
a.plot(Tp - 273.15, curves['measured_fit'], color=INK, lw=2, zorder=4, label=f'measured fit (width {ex.width_obs_C:.1f} °C)')
for key, w in (('aon', ex.width_allornone_C), ('docking', ex.width_uniformCeff_docking_C), ('reporter', ex.width_uniformCeff_reporter_C)):
    col, ls, lab = MODEL[key]
    a.plot(Tp - 273.15, curves[key], color=col, ls=ls, lw=2, zorder=3, label=f'{lab} (width {w:.1f} °C)')
a.axvline(ex.Tm_C, color=GRID, lw=1, zorder=0)
a.set(xlim=(25, 60), ylim=(-0.08, 1.08), xlabel='Temperature (°C)', ylabel='fraction bound  θ')
a.set_title(f'one construct ({code}; 7nt, 4 sites)', loc='left', fontsize=9.5, weight='bold')
a.legend(fontsize=7, frameon=False, loc='lower left')

b = ax[1]; lim = (0.3, 12)
for key, col_key in (('aon', 'width_allornone_C'), ('docking', 'width_uniformCeff_docking_C'), ('reporter', 'width_uniformCeff_reporter_C')):
    col, _, lab = MODEL[key]
    for nt in (7, 5, 6):
        s = D[D.nt == nt]
        if len(s): b.scatter(s.width_obs_C, s[col_key], s=14, marker=MARK[nt], color=col, alpha=.6, edgecolor='white',
                             linewidth=.4, zorder=2, label=lab if nt == 7 else None)
b.plot(lim, lim, color=INK2, lw=1, ls='--', zorder=1)
b.text(3.2, 2.6, 'model = data', color=INK2, fontsize=7.3, ha='left', va='top')
b.text(11, 0.34, 'below line: model too sharp\nabove line: model too broad', color=INK2, fontsize=7, ha='right', va='bottom')
b.set(xscale='log', yscale='log', xlim=lim, ylim=(0.3, 40), xlabel='observed width  dx (°C)', ylabel='model-predicted width (°C)')
b.set_xticks([0.5, 1, 2, 5, 10]); b.set_xticklabels(['0.5', '1', '2', '5', '10'])
b.set_yticks([0.5, 1, 2, 5, 10, 20]); b.set_yticklabels(['0.5', '1', '2', '5', '10', '20'])   # headroom above data for the legend
b.grid(True, which='major', color=GRID, lw=.6, zorder=0)
b.set_title(f'all constructs (n={len(D)}): pure readouts bracket the data', loc='left', fontsize=9.5, weight='bold')
b.legend(fontsize=7, frameon=False, loc='upper left', title='●  7nt   ■  5nt   ▲  6nt', title_fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(OUTDIR, 'figure.svg'))
if '--png' in sys.argv: fig.savefig(sys.argv[sys.argv.index('--png') + 1], dpi=200)
plt.close(fig)

ratio = {k: float((D[c] / D.width_obs_C).median()) for k, c in
         (('all-or-none', 'width_allornone_C'), ('docking', 'width_uniformCeff_docking_C'), ('reporter', 'width_uniformCeff_reporter_C'))}
print(f"Panel H: n={len(D)} constructs; example {code}; median predicted/observed width: "
      + ", ".join(f"{k} {v:.2f}" for k, v in ratio.items())
      + f"; median log10(Ceff/C_tile) docking {D.log10_Ceff_over_Ctile_docking.median():.2f}, "
        f"reporter {D.log10_Ceff_over_Ctile_reporter.median():.2f}")

# ============================ PANEL I: readout weight w ============================
OUT_I = os.path.join("Figure2_Updated", "PanelI_ReadoutWeight"); os.makedirs(OUT_I, exist_ok=True)
BL = {5: '#E8A33D', 6: '#4DA167', 7: '#3B6AA0'}          # Figure 2 bond-length colours (entity = bond length here)

def theta_w(TK, dup, logC, w):                            # theta = f * (w + (1-w) * p_reporter)
    C = 10.0 ** logC
    H = np.array([d[0] for d in dup]); S = np.array([d[1] for d in dup]); rep = np.array([d[2] for d in dup])
    x = np.exp(-(H[None, :] - TK[:, None] * S[None, :]) / (R * TK[:, None])) * C
    Kd = np.expm1(np.sum(np.log1p(x), axis=1)) / C; a = Kd * C_TILE
    f = np.where(a > 1e-12, ((2 * a + 1) - np.sqrt(4 * a + 1)) / (2 * a), a)
    return f * (w + (1 - w) * (x[:, rep] / (1 + x[:, rep])).mean(axis=1))

TKw = np.arange(0.0, 90.0, 0.05) + 273.15
W = np.round(np.arange(0.0, 1.0001, 0.05), 3)
def width_ratio(r, dup, w):
    TmK = r.Tm_C + 273.15
    lc = brentq(lambda v: theta_w(np.array([TmK]), dup, v, w)[0] - 0.5, -14, 4)
    return width(TKw, theta_w(TKw, dup, lc, w)) / r.width_obs_C

curves, wbest, at_bound = [], [], []
for _, r in D.iterrows():
    dup = duplexes(r.code); rr = np.array([width_ratio(r, dup, w) for w in W]); curves.append(rr)
    s = np.sign(np.log(rr)); idx = np.where(np.diff(s) != 0)[0]
    if len(idx):                                          # bracketed: refine the crossing
        i = idx[0]; wbest.append(brentq(lambda w: math.log(width_ratio(r, dup, w)), W[i], W[i + 1], xtol=1e-3)); at_bound.append(False)
    else:                                                 # never crosses: report the nearer bound, flagged
        wbest.append(1.0 if rr[-1] > 1 else 0.0); at_bound.append(True)
DI = D[['code', 'nt', 'n_sites', 'n_reps', 'Tm_C', 'width_obs_C']].copy()
DI['w_best'] = np.round(wbest, 3); DI['w_at_bound'] = at_bound
DI.to_csv(os.path.join(OUT_I, 'data.csv'), index=False); DI.to_csv(os.path.join(SRC, 'Fig2U_I_ReadoutWeight.csv'), index=False)
CR = pd.DataFrame(np.round(np.array(curves), 4), columns=[f'ratio_w={w:g}' for w in W]); CR.insert(0, 'code', D.code.values)
CR.to_csv(os.path.join(OUT_I, 'width_ratio_vs_w.csv'), index=False)

fig, ax = plt.subplots(1, 2, figsize=(9.8, 4.1), gridspec_kw={'width_ratios': [1.1, 1]})
a = ax[0]; Y = np.array(curves)
for nt, rr in zip(D.nt, Y): a.plot(W, rr, color=BL[nt], alpha=.12, lw=.8, zorder=1)
for nt in (7, 5, 6):
    msk = D.nt.values == nt
    if msk.any(): a.plot(W, np.median(Y[msk], axis=0), color=BL[nt], lw=2.4, zorder=3, label=f'{nt}nt median (n={int(msk.sum())})')
a.axhline(1, color=INK, lw=1, ls='--', zorder=2); a.text(0.01, 1.06, 'predicted = observed width', color=INK2, fontsize=7.3)
wc = float(np.median(DI.w_best)); a.axvline(wc, color=INK2, lw=1, ls=':', zorder=2)
a.text(wc + .01, 0.36, f'median best w = {wc:.2f}', color=INK2, fontsize=7.3)
a.set(yscale='log', ylim=(0.3, 6), xlim=(0, 1), xlabel='readout weight w', ylabel='predicted / observed width')
a.set_yticks([0.5, 1, 2, 5]); a.set_yticklabels(['0.5', '1', '2', '5'])
a.set_xticks([0, .25, .5, .75, 1]); a.set_xticklabels(['0\nreporter', '0.25', '0.5', '0.75', '1\ndocking'])
a.grid(True, which='major', axis='y', color=GRID, lw=.6, zorder=0)
a.set_title('width ratio vs readout weight (one line per construct)', loc='left', fontsize=9.5, weight='bold')
a.legend(fontsize=7.2, frameon=False, loc='upper right')

groups = sorted([(nt, N) for (nt, N), s in DI.groupby(['nt', 'n_sites']) if len(s) >= 4], key=lambda t: (t[0] != 7, t[0], t[1]))
def draw_bestfit_w(b):
    rng = np.random.default_rng(0)
    for i, (nt, N) in enumerate(groups):
        s = DI[(DI.nt == nt) & (DI.n_sites == N)]; v = s.w_best.values; j = rng.uniform(-.17, .17, len(v))
        b.scatter(i + j[~s.w_at_bound.values], v[~s.w_at_bound.values], s=13, color=BL[nt], alpha=.75, edgecolor='none', zorder=2)
        b.scatter(i + j[s.w_at_bound.values], v[s.w_at_bound.values], s=16, facecolor='none', edgecolor=BL[nt], lw=.9, zorder=2)
        md = np.median(v); b.plot([i - .3, i + .3], [md, md], color=INK, lw=2.2, zorder=3)
        b.plot([i, i], np.percentile(v, [25, 75]), color=INK, lw=1.1, zorder=3)
        b.text(i, 1.05, f'{md:.2f}', ha='center', fontsize=7.5, color=INK)
    b.axhline(wc, color=INK2, lw=1, ls=':', zorder=1)
    b.set_xticks(range(len(groups)))
    b.set_xticklabels([f'{N} sites\n{nt}nt (n={int(((DI.nt == nt) & (DI.n_sites == N)).sum())})' for nt, N in groups])
    b.set(ylim=(-0.05, 1.12), ylabel='best-fit readout weight w')
    b.set_yticks([0, .25, .5, .75, 1]); b.set_yticklabels(['0  reporter', '0.25', '0.5', '0.75', '1  docking'])
    b.grid(True, axis='y', color=GRID, lw=.6, zorder=0)
    b.set_title('best-fit w per construct', loc='left', fontsize=9.5, weight='bold')
    b.text(len(groups) - .5, -0.03, 'open = never reaches width match (shown at bound)', ha='right', va='bottom', fontsize=7, color=INK2)
draw_bestfit_w(ax[1])
fig.tight_layout(); fig.savefig(os.path.join(OUT_I, 'figure.svg'))
_f1, _b1 = plt.subplots(figsize=(5.2, 3.9)); draw_bestfit_w(_b1); _f1.tight_layout()   # SI: best-fit w only
_f1.savefig(os.path.join(OUT_I, 'figure_bestfit_w.svg')); plt.close(_f1)
if '--png' in sys.argv:
    p = sys.argv[sys.argv.index('--png') + 1]; fig.savefig(p.replace('.png', '_I.png'), dpi=200)
plt.close(fig)
from scipy.stats import spearmanr
print(f"Panel I: median best w {wc:.2f} (n={len(DI)}, {int(DI.w_at_bound.sum())} at a bound); by group "
      + ", ".join(f"{nt}nt/{N}: {DI[(DI.nt == nt) & (DI.n_sites == N)].w_best.median():.2f}" for nt, N in groups)
      + "; within-length Spearman(w, sites): "
      + ", ".join(f"{nt}nt {spearmanr(DI[DI.nt == nt].n_sites, DI[DI.nt == nt].w_best)[0]:+.2f}" for nt in (7, 5)))

# ============ MAIN-FIGURE PANEL H: binding cooperativity free energy (DDG_coop) ============
# For every add-one-bond pair (lo -> hi, hi = lo with one empty site filled), C_eff is fitted SEPARATELY to each
# construct's measured Tm (full multivalent model, no all-or-none):
#   DDG_coop = -RT ln(C_eff,hi / C_eff,lo)   at 37 C      (> 0: the added bond LOWERS the effective concentration felt
#                                                            by every bond -- it helps less than at constant C_eff)
# Same information as the constant-C_eff Tm shortfall (also kept: C_eff fitted to lo, hi's Tm predicted with it;
#   shortfall = measured dTm - predicted dTm, < 0 = helps less); the two are rank-equivalent (|Spearman| > 0.99).
# Reported for three readouts (docking w=1, intermediate w=0.57, reporter w=0) so the conclusion is readout-robust.
OUT_C = os.path.join("Figure2_Updated", "PanelH_BindingCooperativity"); os.makedirs(OUT_C, exist_ok=True)
from scipy.stats import wilcoxon
READ = [('docking', 1.0, '#6250d6'), ('intermediate', 0.57, '#1baf7a'), ('reporter', 0.0, '#e34948')]
RT37 = R * 310.15                                          # kcal/mol
Dc = D.set_index('code'); codes = set(Dc.index); nbd = lambda s: sum(ch != '0' for ch in s)
DUP = {c: duplexes(c) for c in codes}
LCF = {}
def logC_fit(c, w):                                        # log10 C_eff that reproduces construct c's measured Tm
    if (c, w) not in LCF:
        LCF[(c, w)] = brentq(lambda v: theta_w(np.array([Dc.loc[c, 'Tm_C'] + 273.15]), DUP[c], v, w)[0] - 0.5, -14, 4)
    return LCF[(c, w)]
crow = []
for hi in sorted(codes):
    for p_, ch in enumerate(hi):
        if ch == '0': continue
        lo = hi[:p_] + '0' + hi[p_ + 1:]
        if lo not in codes or nbd(lo) != nbd(hi) - 1: continue
        L, Hh = Dc.loc[lo], Dc.loc[hi]
        for name, w, _ in READ:
            lc, lh = logC_fit(lo, w), logC_fit(hi, w)
            try: Thi = brentq(lambda T: theta_w(np.array([T]), DUP[hi], lc, w)[0] - 0.5, 250, 400) - 273.15
            except ValueError: Thi = np.nan
            crow.append(dict(lo=lo, hi=hi, regime=f'{nbd(lo)}->{nbd(hi)}', nt_added={'5': 5, 'X': 5, '6': 6, '7': 7, '9': 9, 'A': 7, 'B': 7}.get(ch, np.nan),
                             readout=name, w=w, log10Ceff_lo=round(lc, 4), log10Ceff_hi=round(lh, 4),
                             Ceff_ratio_hi_lo=round(10 ** (lh - lc), 4), ddG_coop_kcal=round(-RT37 * np.log(10) * (lh - lc), 4),
                             dTm_meas_C=round(Hh.Tm_C - L.Tm_C, 4), dTm_constCeff_C=round(Thi - L.Tm_C, 4),
                             shortfall_C=round((Hh.Tm_C - L.Tm_C) - (Thi - L.Tm_C), 4)))
DC = pd.DataFrame(crow)
DC.to_csv(os.path.join(OUT_C, 'data.csv'), index=False); DC.to_csv(os.path.join(SRC, 'Fig2U_H_BindingCooperativity.csv'), index=False)
CREG = [r for r in ('4->5', '9->10', '10->11') if (DC.regime == r).any()]
summ = []
for reg in CREG:
    for name, w, _ in READ:
        q = DC[(DC.regime == reg) & (DC.readout == name)]; g_ = q.ddG_coop_kcal
        summ.append(dict(regime=reg, readout=name, n_pairs=len(q), median_ddG_coop_kcal=round(g_.median(), 3),
                         median_Ceff_ratio_hi_lo=round(10 ** (-g_.median() / (RT37 * np.log(10))), 3),
                         pct_cost=round(100 * (g_ > 0).mean(), 1), wilcoxon_p=(wilcoxon(g_).pvalue if len(q) >= 6 else np.nan),
                         median_shortfall_C=round(q.shortfall_C.median(), 3)))
pd.DataFrame(summ).to_csv(os.path.join(OUT_C, 'summary.csv'), index=False)

def draw_binding_cooperativity(ax, DC, compact=False):
    """Strip plot of DDG_coop by added site, one colour per readout (mirrored in make_updated_figure2.panel_H_coop)."""
    REGL = {'4->5': '4→5\n7nt added', '9->10': '9→10\n5nt added', '10->11': '10→11\n5nt added'}
    regs = [r for r in ('4->5', '9->10') if (DC.regime == r).any()]      # 10->11 dropped (shares one 11-site construct)
    rng_ = np.random.default_rng(0); fs = 6.2 if compact else 6.8
    plot_read = [(nm, w, col) for (nm, w, col) in READ if nm == 'docking']    # main figure: docking readout only
    off = list(np.linspace(-.27, .27, len(plot_read))) if len(plot_read) > 1 else [0.0]
    lim = np.abs(DC[DC.readout.isin([r[0] for r in plot_read])].ddG_coop_kcal).max() * 1.05
    ax.axhline(0, color=INK, lw=1, zorder=1)
    for i, reg in enumerate(regs):
        if reg == '10->11': ax.axvspan(i - .45, i + .45, color=GRID, alpha=.55, zorder=0)
        for k, (name, w, col) in enumerate(plot_read):
            s = DC[(DC.regime == reg) & (DC.readout == name)].ddG_coop_kcal.values
            if not len(s): continue
            x = i + off[k]
            ax.scatter(x + rng_.uniform(-.08, .08, len(s)), s, s=5 if compact else 9, color=col, alpha=.45, edgecolor='none', zorder=2,
                       label=f'{name}' + (f' (w={w:g})' if 0 < w < 1 else '') if i == 0 else None)
            md = np.median(s); ax.plot([x - .1, x + .1], [md, md], color=INK, lw=2, zorder=3)
            ax.plot([x, x], np.percentile(s, [25, 75]), color=INK, lw=1, zorder=3)
            ax.text(x, lim * 1.08, f'{md:+.2f}', ha='center', va='bottom', fontsize=fs, color=INK)
    n = DC[DC.readout == READ[0][0]].groupby('regime').size()
    ax.set_xticks(range(len(regs))); ax.set_xticklabels([f'{REGL[r]}\n(n={n.get(r, 0)})' for r in regs], fontsize=7 if compact else 8)
    ax.set(xlim=(-.55, len(regs) - .45), ylim=(-lim, lim * 1.25))
    ax.set_ylabel(r'$\Delta\Delta G_{coop} = -RT\,\ln(C_{eff,hi}/C_{eff,lo})$ (kcal/mol)')
    ax.text(1, lim * .62, '↑ added bond lowers $C_{eff}$\n(helps less)', ha='center', va='center', fontsize=fs, color=INK2)
    ax.text(len(regs) - 1, -lim * .5, '↓ added bond raises $C_{eff}$\n(helps more)', ha='center', va='center', fontsize=fs, color=INK2)
    ax.grid(True, axis='y', color=GRID, lw=.6, zorder=0)
    if len(plot_read) > 1:
        ax.legend(title='readout', fontsize=fs, title_fontsize=fs, frameon=False, loc='lower left', bbox_to_anchor=(.36, 0), markerscale=1.5)
    else:
        ax.text(.02, .02, f'{plot_read[0][0]} readout', transform=ax.transAxes, fontsize=fs, color=INK2)

fig, ax = plt.subplots(figsize=(6.6, 4.4)); draw_binding_cooperativity(ax, DC)
ax.set_title('cooperativity free energy of the added bond (>0: helps less)', loc='left', fontsize=9.5, weight='bold')
fig.tight_layout(); fig.savefig(os.path.join(OUT_C, 'figure.svg'))
if '--png' in sys.argv:
    p = sys.argv[sys.argv.index('--png') + 1]; fig.savefig(p.replace('.png', '_Hcoop.png'), dpi=200)
plt.close(fig)
print("Panel H (DDG_coop): " + "; ".join(f"{r['regime']} {r['readout']} {r['median_ddG_coop_kcal']:+.3f} kcal/mol (n={r['n_pairs']})"
                                                      for r in summ))
