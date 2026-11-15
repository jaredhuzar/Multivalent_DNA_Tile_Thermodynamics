"""
SI figure: melting curves of the 9-nt-bond interfaces that passed QC (wells retained in Filtered_DF_std.csv;
12 of the 48 designs on plate DNAExtensionThermo10_30, 3 wells each).

For each design (edge code), each replicate well's raw fluorescence is fitted with a Boltzmann sigmoid
    y = A1 + (A2 - A1) / (1 + exp((x0 - T)/dx))
using the same droop-trimmed fit as Figure2_Updated/make_uniform_ceff_panel.py (fit, find the post-transition droop
onset, refit up to onset + 3 C; dx >= 0.3 C).  Panels show the RAW fluorescence (points) and the fitted sigmoid (lines)
in raw units; a fit counts as valid only if its midpoint lies within the measured range and dx < 10 C (otherwise the
well shows no transition in range and no line is drawn).  fraction_melted in raw_points.csv uses the fitted baselines
for valid fits only.
Tm shown per panel: the pipeline's QC'd Tm (Filtered_DF_std.csv, median over the wells that passed QC) when the
design passed QC; designs that did not pass are marked "failed QC" with the median refit x0 shown for reference.
(9-nt designs are excluded from the main analyses because their NN thermodynamics are unreliable.)
Writes figure.svg, raw_points.csv (long format), fits.csv (per well), summary.csv (per design) here.
"""
import os, numpy as np, pandas as pd
from scipy.optimize import curve_fit
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(os.path.dirname(HERE)); os.chdir(REPO)
PLATE = 'DNAExtensionThermo10_30_Formatted.csv'
raw = pd.read_csv(os.path.join('Melting_Results', PLATE), encoding='latin1')
T = pd.to_numeric(raw.iloc[:, 1], errors='coerce').values
STD = pd.read_csv('Filtered_DF_std.csv', dtype={'Code': str}); STD = STD[STD.Experiment == PLATE]
QC = STD.set_index('Samples')
wells = [k for k in raw.columns[2:] if '9' in str(k).partition('.')[0] and k in QC.index]   # QC-passed wells only

def boltz(x, A1, A2, x0, dx): return A1 + (A2 - A1) / (1 + np.exp((x0 - x) / dx))
BND = ([-np.inf, -np.inf, -np.inf, 0.3], [np.inf] * 4)
def fit_well(t, y):
    p0 = [np.percentile(y, 5), np.percentile(y, 95), np.median(t), 3.0]
    p, _ = curve_fit(boltz, t, y, p0=p0, maxfev=20000, bounds=BND)
    sm = pd.Series(y).rolling(7, center=True, min_periods=1).median().values
    post = t > p[2] + 1.5 * p[3]
    onset = t[np.where(post)[0][np.argmax(sm[post])]] if post.sum() >= 3 else t.max()
    keep = t <= min(onset + 3.0, t.max())
    p, _ = curve_fit(boltz, t[keep], y[keep], p0=p0, maxfev=20000, bounds=BND)
    return p, keep

pts, fits = [], []
for w in wells:
    y = pd.to_numeric(raw[w], errors='coerce').values; ok = np.isfinite(T) & np.isfinite(y)
    t, yy = T[ok], y[ok]; o = np.argsort(t); t, yy = t[o], yy[o]
    code = w.partition('.')[0].replace('X', '5')
    try:
        p, keep = fit_well(t, yy); ok_fit = bool(p[1] != p[0]) and (t.min() - 1 <= p[2] <= t.max()) and p[3] < 10
    except Exception:
        p, keep, ok_fit = [np.nan] * 4, np.ones(len(t), bool), False
    norm = (yy - p[0]) / (p[1] - p[0]) if ok_fit else np.full(len(yy), np.nan)
    q = QC.loc[w] if w in QC.index else None
    fits.append(dict(code=code, well=w, refit_ok=ok_fit, A1=p[0], A2=p[1], refit_Tm_C=p[2] if ok_fit else np.nan, refit_dx_C=p[3],
                     fit_Tmax_C=float(t[keep].max()), passed_QC=q is not None,
                     pipeline_Tm_C=float(q.Tm) if q is not None else np.nan, pipeline_dx_C=float(q.dx) if q is not None else np.nan))
    for tt, yr, yn, kk in zip(t, yy, norm, keep):
        pts.append(dict(code=code, well=w, temperature_C=tt, raw_signal=yr, fraction_melted=yn, used_in_fit=bool(kk)))
P = pd.DataFrame(pts); Fw = pd.DataFrame(fits)
S = Fw.groupby('code').agg(n_wells=('well', 'size'), n_passed_QC=('passed_QC', 'sum'), Tm_C=('pipeline_Tm_C', 'median'),
                           Tm_sd_C=('pipeline_Tm_C', 'std'), refit_Tm_C=('refit_Tm_C', 'median')).reset_index()
S['passed_QC'] = S.n_passed_QC > 0
S = S.sort_values(['passed_QC', 'Tm_C', 'code'], ascending=[False, False, True]).reset_index(drop=True)
P.round(5).to_csv(os.path.join(HERE, 'raw_points.csv'), index=False)
Fw.round(4).to_csv(os.path.join(HERE, 'fits.csv'), index=False); S.round(3).to_csv(os.path.join(HERE, 'summary.csv'), index=False)

# ---- figure: one panel per design ----
plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False, 'svg.fonttype': 'none'})
ncol = 4; nrow = int(np.ceil(len(S) / ncol)); RC = ['#3B6AA0', '#E8A33D', '#1baf7a']; tt = np.linspace(np.nanmin(T), np.nanmax(T), 300)
fig, axs = plt.subplots(nrow, ncol, figsize=(3.0 * ncol, 2.35 * nrow), sharex=True)
for ax, (_, s) in zip(axs.flat, S.iterrows()):
    for k, (_, f) in enumerate(Fw[Fw.code == s.code].sort_values('well').iterrows()):
        d = P[P.well == f.well]; col = RC[k % 3]
        ax.scatter(d.temperature_C, d.raw_signal, s=4, color=col, alpha=.5 if f.passed_QC else .35, edgecolor='none')
        if f.refit_ok:
            ax.plot(tt, boltz(tt, f.A1, f.A2, f.refit_Tm_C, f.refit_dx_C), color=col, lw=1.1, ls='-' if f.passed_QC else (0, (3, 2)))
    if s.passed_QC:
        ax.axvline(s.Tm_C, color='#0b0b0b', lw=.8, ls=':')
        lab = f'$T_m$ {s.Tm_C:.1f} °C' + (f' ± {s.Tm_sd_C:.1f}' if np.isfinite(s.Tm_sd_C) else '') + f'  ({int(s.n_passed_QC)}/{int(s.n_wells)} wells)'
        ax.set_title(f'{s.code}\n{lab}', fontsize=7.5)
    else:
        sub = f'refit $x_0$ {s.refit_Tm_C:.1f} °C' if np.isfinite(s.refit_Tm_C) else 'no transition in range'
        ax.set_title(f'{s.code}\nfailed QC ({sub})', fontsize=7.5, color='#c0392b')
        ax.set_facecolor('#fbf6f5')
    ax.grid(True, color='#eeeeee', lw=.5)
for ax in list(axs.flat)[len(S):]: ax.axis('off')
for ax in axs[-1, :]: ax.set_xlabel('temperature (°C)')
for ax in axs[:, 0]: ax.set_ylabel('fluorescence (a.u.)')
fig.suptitle(f'9-nt-bond interfaces that passed QC ({len(S)} designs, {len(Fw)} wells): raw (points) and fitted (lines) melting curves',
             fontsize=10, y=1.0)
fig.tight_layout(); fig.savefig(os.path.join(HERE, 'figure.svg'))
import sys
if '--png' in sys.argv: fig.savefig(sys.argv[sys.argv.index('--png') + 1], dpi=110)
print(f'9-nt melting curves: {len(S)} designs, {len(Fw)} wells, {int(S.passed_QC.sum())} designs passed QC; '
      f'refit vs pipeline Tm (QC wells) median |diff| {np.nanmedian(np.abs(Fw.refit_Tm_C - Fw.pipeline_Tm_C)):.2f} C')
