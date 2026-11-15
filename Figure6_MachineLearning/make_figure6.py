"""
Figure 6 (updated 2026-11-12 for the multivalent model of Figure 4).  Run after Figure4_MultivalentModel/fit_model.py.

Manuscript panels produced here (folder names keep their historical letters):
  B  XGBoost out-of-fold predicted vs measured Tm                       -> Fig6b_XGBoost_StabilityScatter/
  C  held-out Spearman rho (bars: out-of-fold rho; circles: fold mean +- s.d.) and RMS of six models, the
     physical comparator now being the multivalent model of Figure 4      -> Fig6d_ModelComparison_SpearmanRho/
  E  edge-site effect is a spacing effect: four-site 7-nt interfaces with vs without an occupied edge site
     (measured Tm; Tm after regressing out mean gap; measured minus multivalent-model Tm) -> Fig6e_EdgeSites_Spacing/
  combined (A and D as placeholders: A = definitions artwork, D = SHAP from Supplementary/S6_SHAP_Values)
                                                                         -> figure6_combined.svg
Data: the canonical 236-interface set and folds of Table S1 / Figure 4 (KFold(10, shuffle=True, random_state=0)).
ML features (as in TmPredictor_NewTm.ipynb): 11 per-site values = mean salt-corrected NN dG at 37 C of the site's
two duplexes, 0 for an empty site.  Multivalent-model out-of-fold predictions are read from
Figure4_MultivalentModel/model_fit/oof_predictions.csv (same folds).   Optional: --png PATH (combined preview).
"""
import os, sys, numpy as np, pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
from sklearn.model_selection import KFold
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures
from xgboost import XGBRegressor
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(REPO, 'Figure4_MultivalentModel'))
import mv_model as mv
os.chdir(REPO)

CODES = list(mv.CODES); Y = np.asarray(mv.Y_C, float); n = len(Y); TH = mv.UM['TH']
feat = mv.site_dG_features                         # shared with Supplementary/S6_SHAP_Values (same 236, same features)
X = np.vstack([feat(c) for c in CODES])
FOLDS = list(KFold(10, shuffle=True, random_state=0).split(np.arange(n)))

class SumModel:
    def fit(self, X, y): self.m = LinearRegression().fit(X.sum(1, keepdims=True), y); return self
    def predict(self, X): return self.m.predict(X.sum(1, keepdims=True))
MODELS = {'XGBoost': lambda: XGBRegressor(n_estimators=500, learning_rate=0.05, max_depth=4, subsample=0.8,
                                          colsample_bytree=0.8, tree_method='hist', random_state=0),
          'Random forest': lambda: RandomForestRegressor(n_estimators=300, random_state=0),
          'Quadratic regression': lambda: make_pipeline(PolynomialFeatures(2, include_bias=False, interaction_only=True), LinearRegression()),
          'Linear regression': LinearRegression,
          'Summed NN ΔG': SumModel}
OOF = {}
for name, mk in MODELS.items():
    o = np.full(n, np.nan)
    for tr, te in FOLDS:
        m = mk().fit(X[tr], Y[tr]); o[te] = m.predict(X[te])
    OOF[name] = o
MVO = pd.read_csv(os.path.join(REPO, 'Figure4_MultivalentModel/model_fit/oof_predictions.csv'), dtype={'code': str}).set_index('code')
assert list(MVO.index) == CODES, 'Figure 4 out-of-fold predictions are not in the canonical construct order'
OOF['Multivalent model (Fig. 4)'] = MVO.Tm_pred_oof_C.values
GROUP = MVO.group.values

rows = []
for name, o in OOF.items():
    fr = [spearmanr(Y[te], o[te])[0] for _, te in FOLDS]
    rows.append(dict(Model=name, OOF_Rho=spearmanr(Y, o)[0], Mean_Fold_Rho=np.mean(fr), Std_Fold_Rho=np.std(fr),
                     OOF_RMS_C=np.sqrt(np.mean((o - Y) ** 2)), n=n))
PERF_ALL = pd.DataFrame(rows).sort_values('OOF_Rho', ascending=False).reset_index(drop=True)
PERF = PERF_ALL[PERF_ALL.Model != 'Summed NN ΔG'].reset_index(drop=True)   # panel C drops the additive-NN series
PERF_ALL.round(4).to_csv(os.path.join(HERE, 'Fig6d_ModelComparison_SpearmanRho', 'all_models_incl_summedNN.csv'), index=False)
print(PERF_ALL.round(3).to_string(index=False))

# ---- edge-site analysis: four-site 7-nt interfaces ----
sites = [[i for i, ch in enumerate(c) if ch != '0'] for c in CODES]
sel = np.array([(GROUP[i] == '7-nt') and len(s) == 4 for i, s in enumerate(sites)])
E = pd.DataFrame(dict(code=np.array(CODES)[sel], Tm_meas_C=Y[sel], Tm_mv_oof_C=OOF['Multivalent model (Fig. 4)'][sel],
                      edge=[(0 in s) or (10 in s) for s, k in zip(sites, sel) if k],
                      mean_gap_nm=[(s[-1] - s[0]) / 3 * 6.0 for s, k in zip(sites, sel) if k]))
b = np.polyfit(E.mean_gap_nm, E.Tm_meas_C, 1); E['Tm_resid_gap_C'] = E.Tm_meas_C - np.polyval(b, E.mean_gap_nm)
E['Tm_resid_model_C'] = E.Tm_meas_C - E.Tm_mv_oof_C
def mw(col): return mannwhitneyu(E[E.edge][col], E[~E.edge][col], alternative='two-sided').pvalue
ES = dict(n_edge=int(E.edge.sum()), n_noedge=int((~E.edge).sum()),
          Tm_edge=E[E.edge].Tm_meas_C.median(), Tm_noedge=E[~E.edge].Tm_meas_C.median(), P_Tm=mw('Tm_meas_C'),
          gap_edge=E[E.edge].mean_gap_nm.mean(), gap_noedge=E[~E.edge].mean_gap_nm.mean(), P_resid_gap=mw('Tm_resid_gap_C'),
          model_edge=E[E.edge].Tm_mv_oof_C.median(), model_noedge=E[~E.edge].Tm_mv_oof_C.median(), P_model=mw('Tm_mv_oof_C'),
          P_resid_model=mw('Tm_resid_model_C'))
print('edge-site:', {k: (round(v, 4) if isinstance(v, float) else v) for k, v in ES.items()})

# ---- write panel data ----
def wr(folder, df, note):
    d = os.path.join(HERE, folder); os.makedirs(d, exist_ok=True)
    df.to_csv(os.path.join(d, 'data.csv'), index=False); open(os.path.join(d, 'README.txt'), 'w').write(note)
    return d
dB = wr('Fig6b_XGBoost_StabilityScatter', pd.DataFrame(dict(code=CODES, group=GROUP, exp_Tm=Y, pred_Tm_oof=OOF['XGBoost'])).round(4),
        'Fig 6B. XGBoost out-of-fold (10-fold CV, canonical 236 interfaces, Table S1 folds) predicted vs measured Tm.\n'
        'Written by ../make_figure6.py.\n')
dC = wr('Fig6d_ModelComparison_SpearmanRho', PERF.round(4),
        'Fig 6C (folder keeps its historical "6d" name). Held-out Spearman rho and RMS of five models on the same 236\n'
        'interfaces and folds (the additive summed-NN baseline is dropped from the panel; all six incl. it are in\n'
        'all_models_incl_summedNN.csv); the physical comparator is the multivalent model of Figure 4 (out-of-fold predictions\n'
        'from Figure4_MultivalentModel/model_fit/). The former comparator (phenomenological A2, rho 0.933) is in\n'
        'Supplementary/SI_PhenomenologicalModel_A2/. Written by ../make_figure6.py.\n')
dE = wr('Fig6e_EdgeSites_Spacing', E.round(4),
        'Fig 6E. Four-site 7-nt interfaces with vs without an occupied edge site (site 0 or 10): measured Tm, Tm after\n'
        'regressing out mean gap, and measured minus multivalent-model (out-of-fold) Tm. Two-sided Mann-Whitney U.\n'
        + ''.join(f'  {k} = {v:.4g}\n' if isinstance(v, float) else f'  {k} = {v}\n' for k, v in ES.items())
        + 'Written by ../make_figure6.py.\n')

# ---- figure ----
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e6e5e0"
GC = {'7-nt': '#3B6AA0', '5-nt': '#E8A33D', '6-nt': '#4DA167', 'A/B': '#8a6fbf'}
MV_COL, ML_COL = '#1baf7a', '#005CAB'; EDGE_C, NOEDGE_C = '#1f3f66', '#9dbbe0'
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'svg.fonttype': 'none',
                     'axes.edgecolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2})
def grid(ax, axis='y'): ax.grid(True, axis=axis, color=GRID, lw=.6, zorder=0)

def panel_B(ax):
    o = OOF['XGBoost']; r = PERF.set_index('Model').loc['XGBoost']
    for gname in ('7-nt', '5-nt', '6-nt', 'A/B'):
        m = GROUP == gname
        ax.scatter(o[m], Y[m], s=16, color=GC[gname], alpha=.75, edgecolor='white', lw=.3, label=f'{gname} (n={m.sum()})', zorder=3)
    lim = [min(o.min(), Y.min()) - 1, max(o.max(), Y.max()) + 1]
    ax.plot(lim, lim, color=INK2, lw=1, zorder=1); ax.set(xlim=lim, ylim=lim, xlabel='XGBoost predicted $T_m$ (°C, held-out fold)', ylabel='measured $T_m$ (°C)')
    ax.text(.04, .96, f'RMS {r.OOF_RMS_C:.2f} °C\nSpearman ρ = {r.OOF_Rho:.3f}', transform=ax.transAxes, va='top', fontsize=8.5)
    ax.legend(fontsize=7, frameon=False, loc='lower right'); grid(ax)
    ax.set_title('B   XGBoost, out-of-fold', loc='left', fontsize=9.6, weight='bold')

def panel_C(ax):
    P = PERF.iloc[::-1].reset_index(drop=True); yy = np.arange(len(P))
    cols = [MV_COL if 'Multivalent' in m else ML_COL for m in P.Model]
    ax.barh(yy, P.OOF_Rho, color=cols, height=.62, zorder=2)
    ax.errorbar(P.Mean_Fold_Rho, yy, xerr=P.Std_Fold_Rho, fmt='o', color=INK, ms=4.5, capsize=3, lw=1, zorder=3)
    for yi, (rh, rm) in enumerate(zip(P.OOF_Rho, P.OOF_RMS_C)):
        ax.text(1.04, yi, f'ρ {rh:.2f}\nRMS {rm:.2f} °C', va='center', fontsize=7.2, color=INK)
    ax.set_yticks(yy); ax.set_yticklabels(P.Model, fontsize=8)
    ax.set(xlim=(0, 1.3), xticks=[0, .2, .4, .6, .8, 1.0], xlabel='held-out Spearman ρ\n(bar: out-of-fold; circle: fold mean ± s.d.)')
    ax.axvline(0, color=INK2, lw=.8); grid(ax, 'x')
    ax.set_title('C   held-out performance (same folds)', loc='left', fontsize=9.6, weight='bold')

def panel_E(sub):
    axs = [plt.subplot(s) for s in sub]
    rng = np.random.default_rng(0)
    spec = [('Tm_meas_C', 'measured $T_m$ (°C)', ES['P_Tm'], 'measured'),
            ('Tm_resid_gap_C', '$T_m$ residual after mean gap (°C)', ES['P_resid_gap'], 'mean gap regressed out'),
            ('Tm_resid_model_C', 'measured − model $T_m$ (°C)', ES['P_resid_model'], 'minus multivalent model')]
    for a, (col, lab, pv, ttl) in zip(axs, spec):
        for k, (flag, c, nm) in enumerate(((True, EDGE_C, 'edge'), (False, NOEDGE_C, 'no edge'))):
            v = E[E.edge == flag][col].values
            a.scatter(k + rng.uniform(-.16, .16, len(v)), v, s=10, color=c, alpha=.7, edgecolor='none', zorder=2)
            md = np.median(v); q = np.percentile(v, [25, 75])
            a.plot([k - .28, k + .28], [md, md], color=INK, lw=2.2, zorder=3); a.plot([k, k], q, color=INK, lw=1.1, zorder=3)
        a.set_xticks([0, 1]); a.set_xticklabels([f'edge\n(n={ES["n_edge"]})', f'no edge\n(n={ES["n_noedge"]})'], fontsize=7.4)
        a.set_ylabel(lab, fontsize=8); a.set_title(f'{ttl}\nP = {pv:.1g}', fontsize=8, color=INK2); grid(a)
        a.set_xlim(-.6, 1.6)
    axs[0].text(-.05, 1.25, 'E   the edge-site deficit tracks spacing', transform=axs[0].transAxes,
                fontsize=9.6, weight='bold')
    axs[0].annotate(f'mean gap {ES["gap_edge"]:.1f} vs {ES["gap_noedge"]:.1f} nm', xy=(.5, .02), xycoords='axes fraction',
                    ha='center', fontsize=6.8, color=INK2)
    return axs

def placeholder(ax, ttl, txt):
    ax.axis('off'); ax.add_patch(plt.Rectangle((0, 0), 1, 1, fill=False, ec='#bbb', ls='--', lw=1.1, transform=ax.transAxes))
    ax.set_title(ttl, loc='left', fontsize=9.6, weight='bold'); ax.text(.5, .5, txt, ha='center', va='center', fontsize=8, color=INK2, transform=ax.transAxes)

# per-panel svgs
f, a = plt.subplots(figsize=(4.8, 4.6)); panel_B(a); f.tight_layout(); f.savefig(os.path.join(dB, 'figure.svg')); plt.close(f)
f, a = plt.subplots(figsize=(6.4, 3.9)); panel_C(a); f.tight_layout(); f.savefig(os.path.join(dC, 'figure.svg')); plt.close(f)
f = plt.figure(figsize=(8.4, 3.6)); gE = GridSpec(1, 3, figure=f, wspace=.55, top=.72, bottom=.16, left=.08, right=.98)
panel_E([gE[0, 0], gE[0, 1], gE[0, 2]]); f.savefig(os.path.join(dE, 'figure.svg')); plt.close(f)

# combined
fig = plt.figure(figsize=(15, 10.2))
gs = GridSpec(2, 3, figure=fig, hspace=.5, wspace=.42, top=.94, bottom=.08, left=.07, right=.98, width_ratios=[1, 1.25, 1.1])
placeholder(fig.add_subplot(gs[0, 0]), 'A   definitions', 'interface code → 11 per-site NN ΔG values\n(0 = empty site)\n[existing artwork]')
panel_B(fig.add_subplot(gs[0, 1])); panel_C(fig.add_subplot(gs[0, 2]))
placeholder(fig.add_subplot(gs[1, 0]), 'D   SHAP values (XGBoost)', 'unchanged:\nSupplementary/S6_SHAP_Values/figure.svg')
sub = GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[1, 1:], wspace=.6)
panel_E([sub[0, 0], sub[0, 1], sub[0, 2]])
fig.savefig(os.path.join(HERE, 'figure6_combined.svg'))
if '--png' in sys.argv: fig.savefig(sys.argv[sys.argv.index('--png') + 1], dpi=120)
plt.close(fig)
print('wrote Figure6_MachineLearning panels B, C (Fig6d folder), E (Fig6e_EdgeSites_Spacing) + figure6_combined.svg')
