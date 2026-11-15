"""
Assemble SI_updates.docx: every SI figure and table that changed on 2026-11-12, with captions and numbers read from
the current outputs, plus a list of retired and unchanged SI items.  Native Word tables; figures as PNG (SVG sources
converted with LibreOffice).  Run after the SI scripts, Figure 4 and Table S1/S4 builders (see run_all.sh).
Needs pandoc and LibreOffice (soffice) on PATH.
"""
import os, shutil, subprocess, tempfile, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(os.path.dirname(HERE))
IMG = os.path.join(HERE, 'img'); os.makedirs(IMG, exist_ok=True)
J = lambda *p: os.path.join(REPO, *p)

def png(src, name):
    dst = os.path.join(IMG, name + '.png')
    if src.endswith('.png'): shutil.copyfile(src, dst); return dst
    with tempfile.TemporaryDirectory() as td:
        subprocess.run(['soffice', '--headless', '--convert-to', 'png', '--outdir', td, src], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        shutil.move(os.path.join(td, os.path.splitext(os.path.basename(src))[0] + '.png'), dst)
    return dst

# ---- numbers ----
ols = pd.read_csv(J('Supplementary/Ceff_Geometry_OLS/partial_correlations.csv')).set_index('predictor')
olsc = pd.read_csv(J('Supplementary/Ceff_Geometry_OLS/ols_coefficients.csv'))
ced = pd.read_csv(J('Supplementary/Ceff_Geometry_OLS/data.csv'))
bl = pd.read_csv(J('Supplementary/BondLength_Ceff/fits.csv')); blk, bln = bl.iloc[0], bl.iloc[1]
pv = pd.read_csv(J('Supplementary/PredictedVsMeasured_NN_Multivalent/data.csv'))
from scipy.stats import spearmanr
rms = lambda a, b: float(np.sqrt(np.mean((a - b) ** 2)))
w = pd.read_csv(J('Supplementary/Width_vs_MeanGap_7nt_4bonds/data.csv'))
q1, q2 = np.percentile(w.mean_gap_nm, [100 / 3, 200 / 3]); lo_, hi_ = w[w.mean_gap_nm <= q1], w[w.mean_gap_nm >= q2]
cv = pd.read_csv(J('Figure4_MultivalentModel/model_fit/cv_comparison.csv'))
a2 = cv[cv.model.str.contains('A2')].iloc[0]
coop = pd.read_csv(J('Figure2_Updated/PanelH_BindingCooperativity/summary.csv'))
ck = lambda reg, rd, col: coop[(coop.regime == reg) & (coop.readout == rd)][col].iloc[0]
T1 = pd.read_csv(J('Supplementary/TableS1_Model_Comparison_All/TableS1_all_models.csv'))
T4 = pd.read_csv(J('Supplementary/TableS4_Parameter_CIs/TableS4_parameters.csv'), keep_default_na=False)
p = lambda x: f'{x:.1g}'

F = []   # (title, image, caption)
F.append(('Effective concentration versus interface geometry (replaces the avidity-α regression; SI Fig. S2 of the 11/1 draft, cited as Fig. S8 in the main text)',
          png(J('Supplementary/Ceff_Geometry_OLS/figure.png'), 'ceff_ols'),
          f"Each interface's effective concentration was fitted on its own to its measured T_m with the multivalent binding model "
          f"(one uniform C_eff per interface, intermediate FRET readout) and regressed on the number of occupied sites, bond span, "
          f"mean bond length and GC fraction. Panels are added-variable plots: log10 C_eff and each factor after regressing out the "
          f"other three (OLS R² = {olsc.R2.iloc[0]:.2f}, n = {len(ced)}; median C_eff {10 ** ced.log10_Ceff_M.median() * 1e3:.1f} mM). "
          f"Bond span dominates (partial ρ = {ols.loc['span_nm', 'partial_rho']:+.2f}, P = {p(ols.loc['span_nm', 'partial_p'])}); bond length "
          f"(ρ = {ols.loc['bond_len', 'partial_rho']:+.2f}, P = {p(ols.loc['bond_len', 'partial_p'])}) and GC fraction "
          f"(ρ = {ols.loc['GC', 'partial_rho']:+.2f}, P = {p(ols.loc['GC', 'partial_p'])}) have no detectable effect; at fixed span, more sites "
          f"(a denser interface) slightly raise C_eff (ρ = {ols.loc['n_sites', 'partial_rho']:+.2f}, P = {p(ols.loc['n_sites', 'partial_p'])}). "
          f"Number of sites and bond length are collinear in this library (VIF ≈ 7–9). The avidity coefficient α is no longer used: it is "
          f"referenced to the 1 M standard state, so α < 1 and its decline with bond number are expected even when C_eff is constant."))
F.append(('Transition width versus bond spacing (replaces van’t Hoff ΔH and ΔS versus mean gap; SI Fig. S4 of the 11/1 draft, cited as Fig. S10 in the main text)',
          png(J('Supplementary/Width_vs_MeanGap_7nt_4bonds/figure.svg'), 'width_gap'),
          f"Four-site 7-nt interfaces (n = {len(w)}). (A) Width-derived van’t Hoff enthalpy |ΔH_vH| = 1.5 RT_m²/δ rises with mean gap "
          f"(ρ = {spearmanr(w.mean_gap_nm, w.dH_vH_width_kcal)[0]:+.2f}, P = {p(spearmanr(w.mean_gap_nm, w.dH_vH_width_kcal)[1])}). "
          f"(B) The cooperative fraction |ΔH_vH|/|ΔH_NN| rises from {lo_.f_coop.median():.2f} (most clustered third, ≤{q1:.0f} nm) to "
          f"{hi_.f_coop.median():.2f} (most dispersed third, ≥{q2:.0f} nm) while T_m falls from {lo_.Tm_C.median():.1f} to {hi_.Tm_C.median():.1f} °C. "
          f"(C) Measured width (ρ = {spearmanr(w.mean_gap_nm, w.width_obs_C)[0]:+.2f}) and predictions fitted to T_m alone: multivalent model with "
          f"distance-dependent C_eff (ρ = {spearmanr(w.mean_gap_nm, w.width_pred_multivalent_dist_C)[0]:+.2f}); one uniform C_eff "
          f"(ρ = {spearmanr(w.mean_gap_nm, w.width_pred_multivalent_uniform_C)[0]:+.2f}); two-state model A2 "
          f"(ρ = {spearmanr(w.mean_gap_nm, w.width_pred_A2_C)[0]:+.2f}, broadening). In a two-state description ΔS_vH = ΔH_vH/T_m + R ln(4/C_T) "
          f"is fixed by ΔH_vH and T_m, so only the width carries shape information; fitted two-state ΔH° and ΔS° are correlated "
          f"(ρ ≈ 0.999) and show no trend (ρ = {spearmanr(w.mean_gap_nm, w.dH_twostate_fit_kcal)[0]:+.2f}, "
          f"{spearmanr(w.mean_gap_nm, w.dS_twostate_fit_cal)[0]:+.2f})."))
F.append(('Bond length and effective concentration (replaces entropic penalty by bond length; SI Fig. S11 of the 11/1 draft)',
          png(J('Supplementary/BondLength_Ceff/figure.png'), 'bondlength'),
          f"(A) Effective concentration of 5- and 6-nt bonds relative to 7-nt bonds, fitted as an extra factor φ in the multivalent "
          f"model (C_eff,i = C_ref φ_type e^(−d_i/λ) e^(−κ(N−4))), with the valency term (circles; φ5 = {blk.phi5:.2f}, 95% CI "
          f"{blk.phi5_lo:.2f}–{blk.phi5_hi:.2f}; φ6 = {blk.phi6:.2f}) and without it (squares; φ5 = {bln.phi5:.2f}, {bln.phi5_lo:.2f}–{bln.phi5_hi:.2f}; "
          f"φ6 = {bln.phi6:.2f}). A worm-like-chain estimate of the tether effective concentration (red; Table S1 model B4) predicts the "
          f"opposite and much larger dependence (φ5 ≈ {blk.WLC_phi5:.0f}, φ6 ≈ {blk.WLC_phi6:.1f}). (B) Cross-validated RMS error is unchanged "
          f"({blk.base_cvRMS:.2f}, {blk.cvRMS:.2f} and {bln.cvRMS:.2f} °C). φ5 and κ are strongly correlated (r = {blk.corr_phi5_kappa:.2f}) "
          f"because 5-nt bonds dominate the many-site interfaces, so a valency cost and slightly weaker 5-nt bonds cannot be separated "
          f"in this library; the paired comparison of Fig. 2F tests the valency cost within one bond length. Replaces the earlier figure, "
          f"which plotted the noise-biased two-state van’t Hoff entropy."))
F.append(('Sequence-predicted versus measured melting temperature (replaces SI Fig. S12 of the 11/1 draft)',
          png(J('Supplementary/PredictedVsMeasured_NN_Multivalent/figure.png'), 'pred_meas'),
          f"(A) Additive nearest-neighbor two-state T_m, summing all duplexes with each referenced to the 1 M standard state "
          f"(RMS {rms(pv.Tm_NN_additive_C, pv.Tm_meas_C):.0f} °C, ρ = {spearmanr(pv.Tm_NN_additive_C, pv.Tm_meas_C)[0]:.2f}); it predicts "
          f"the many-site 5-nt interfaces to be the most stable (median {pv[pv.group == '5-nt'].Tm_NN_additive_C.median():.1f} versus "
          f"{pv[pv.group == '7-nt'].Tm_NN_additive_C.median():.1f} °C for 7-nt) although they melt lowest. (B) The multivalent model of Fig. 4, "
          f"same sequences, out-of-fold (RMS {rms(pv.Tm_multivalent_oof_C, pv.Tm_meas_C):.2f} °C, ρ = {spearmanr(pv.Tm_multivalent_oof_C, pv.Tm_meas_C)[0]:.3f}), "
          f"n = {len(pv)}."))
F.append(('Uniform effective concentration and the transition width (moved from main-text Fig. 2)',
          png(J('Figure2_Updated/PanelH_UniformCeff_Width/figure.svg'), 'uniform_width'),
          "One uniform C_eff per interface, fitted to its measured T_m with the multivalent model, predicts the transition width. With "
          "the FRET signal taken to report tile docking, predicted widths are 0.61 of measured (too sharp); with it taken to report "
          "the central reporter duplexes, 2.7 (too broad); an all-or-none transition with ΔH = ΣΔH_NN gives 0.37. The two pure readouts "
          "bracket the data and an all-or-none melt is excluded (n = 205; left, a representative 7-nt four-site interface)."))
F.append(('FRET readout weight (moved from main-text Fig. 2)',
          png(J('Figure2_Updated/PanelI_ReadoutWeight/figure.svg'), 'readout_w'),
          "Readout θ = f[w + (1 − w)p_rep]; for each w the uniform C_eff is refit to T_m and w_best is where the predicted width equals "
          "the measured width. Median w ≈ 0.57, flat across bond length and number of sites; this value is fixed in Figs. 2–4."))
F.append(('Cooperative unit in duplexes (moved from main-text Fig. 2; the cooperative fraction is Fig. 2C)',
          png(J('Figure2_Updated/PanelD_CooperativeUnit/figure.svg'), 'coop_unit'),
          "Number of duplexes that melt as one unit, n_coop = ΔH_vH/(ΔH_NN/2N), from the width-derived van’t Hoff enthalpy. About 3 "
          "duplexes at four to five sites rising to 6–7 at nine to ten sites, always fewer than the 2N duplexes of the interface "
          "(dashed); the melt is about 2.4 times broader than an all-or-none transition."))
F.append(('Cooperativity free energy for three FRET readouts (Fig. 2F shows the docking readout)',
          png(J('Figure2_Updated/PanelH_BindingCooperativity/figure.svg'), 'ddg_readouts'),
          f"ΔΔG_coop = −RT ln(C_eff,hi/C_eff,lo) for every pair of interfaces differing by one site, with C_eff fitted separately to each "
          f"construct, for docking, intermediate (w = 0.57) and reporter readouts. 4→5 (7-nt, n = {int(ck('4->5', 'docking', 'n_pairs'))}): "
          f"{ck('4->5', 'docking', 'median_ddG_coop_kcal'):+.2f}, {ck('4->5', 'intermediate', 'median_ddG_coop_kcal'):+.2f} and "
          f"{ck('4->5', 'reporter', 'median_ddG_coop_kcal'):+.2f} kcal mol⁻¹; 9→10 (5-nt, n = {int(ck('9->10', 'docking', 'n_pairs'))}): "
          f"{ck('9->10', 'docking', 'median_ddG_coop_kcal'):+.2f}, {ck('9->10', 'intermediate', 'median_ddG_coop_kcal'):+.2f} and "
          f"{ck('9->10', 'reporter', 'median_ddG_coop_kcal'):+.2f} kcal mol⁻¹. 10→11 (shaded) comprises 8 pairs sharing one 11-site "
          f"construct. Only 9→10 is positive under all three readouts."))
F.append(('Two-state model A2: predicted versus measured T_m (former main-text Fig. 4B)',
          png(J('Supplementary/SI_PhenomenologicalModel_A2/Phenomenological_Model/pred_vs_exp.svg'), 'a2_pred'),
          f"Two-state model T_m = ΣgΔH_NN / (ΣgΔS_NN + NΔS_Bond + ΔS_0 + R ln(C/2)), g = e^(−d/λ), out-of-fold (10-fold CV over the same "
          f"236 interfaces): RMS {a2.cvRMS:.2f} °C, ρ = {a2.cvrho:.3f}. It predicts T_m as well as the multivalent model but, being "
          f"two-state, implies transitions about half as wide as measured and broadening with spacing (Fig. S10). Parameters in Table S4."))
F.append(('Two-state model A2: theoretical curves (former main-text Fig. 4C–F)',
          png(J('Supplementary/SI_PhenomenologicalModel_A2/theory_fig4_cdef.svg'), 'a2_theory'),
          "Enthalpy, entropy and bound fraction at 40 °C versus number of bonds (clustered and spread) and versus spacing, from model A2 "
          "with the fitted parameters of Table S4. Superseded in the main text by the multivalent-model curves of Fig. 4D–H."))

def mdtable(df, cols, heads):
    out = ['| ' + ' | '.join(heads) + ' |', '|' + '|'.join(['---'] * len(heads)) + '|']
    for _, r in df.iterrows(): out.append('| ' + ' | '.join(str(r[c]) for c in cols) + ' |')
    return '\n'.join(out)
t1 = T1.copy(); t1['cvRMS'] = t1.cvRMS.map('{:.2f}'.format); t1['cvR2'] = t1.cvR2.map('{:.3f}'.format); t1['cvrho'] = t1.cvrho.map('{:.3f}'.format)
t1['form'] = t1.form.str.replace('|', '/', regex=False)
t4 = T4.copy(); t4['estimate'] = [f'{v:.1f}' if pn == 'C_ref' else f'{v:.3g}' for v, pn in zip(t4.estimate, t4.parameter)]

md = ['# SI updates (2026-11-12)', '',
      'Figures and tables of the Supplementary Information that change with the multivalent-model revision. Numbers are read '
      'from the current repository outputs. Figure numbers refer to the 11/1/2026 SI draft; the main text cites a different '
      'numbering, so SI citations in the main text need renumbering once the order below is fixed.', '']
for k, (t, im, cap) in enumerate(F, 1):
    md += [f'## Figure U{k}. {t}', '', f'![]({im}){{width=6.5in}}', '', f'**Figure U{k}.** {cap}', '']
md += ['## Table S1. Tested models of interface T_m', '',
       mdtable(t1, ['id', 'model', 'n_params', 'cvRMS', 'cvrho', 'params'],
               ['ID', 'Model', 'k', 'CV RMS (°C)', 'CV ρ', 'Fitted parameters (all data)']), '',
       'All models were fitted to the median T_m of the same 236 interfaces and evaluated by 10-fold cross-validation over interfaces '
       'with identical folds; CV statistics are out-of-fold (CV R² and model forms in the formatted table). M: multivalent binding models (M3 = Fig. 4; "d" = docking readout). '
       'A: two-state, distance-weighted. B: two-state with mechanical-coupling or worm-like-chain tether terms. A formatted version is '
       'in Supplementary/TableS1_Model_Comparison_All/TableS1_all_models.docx.', '',
       '## Table S4. Fitted parameters and uncertainties', '',
       mdtable(t4, ['model', 'parameter', 'units', 'estimate', 'CI', 'meaning'],
               ['Model', 'Parameter', 'Units', 'Estimate', 'Bootstrap 95% CI', 'Meaning']), '',
       'Standard errors are in the formatted table. Multivalent model: construct bootstrap (200 resamples). Two-state A2: jackknife SE and bootstrap (1000). '
       'C_ref SEs are in log10 units. Formatted version: Supplementary/TableS4_Parameter_CIs/TableS4_parameters.docx.', '',
       '## Retired SI items', '',
       '- Avidity α regression (11/1 SI Fig. S2): α is referenced to 1 M; replaced by Figure U1.',
       '- Regression of the experimental free energy ΔG_exp on geometry (11/1 SI Fig. S10): ΔG_exp is the noise-biased van’t Hoff free energy of the cooperative unit; removed (Figure U1 covers the geometry dependence).',
       '- Van’t Hoff ΔH and ΔS versus mean gap (11/1 SI Fig. S4): replaced by Figure U2.',
       '- Entropic penalty by bond length (11/1 SI Fig. S11): replaced by Figure U3.',
       '- NN-predicted versus measured T_m alone (11/1 SI Fig. S12): replaced by Figure U4.',
       '- ΔT_m per added site versus the additive NN prediction (former Fig. 2E draft panel): the additive prediction is referenced to 1 M; superseded by Fig. 2F and Figure U8.',
       '- Old Table S1 (two-state models only, 11/1 numbers): replaced by the Table S1 above.', '',
       '## Unchanged SI items', '',
       '- AFM of tile dimers (S1); geometry-controlled regression of T_m (S3); coarse-grained simulations (S5); SHAP values (S6); oxDNA fluctuations (S7); T_m by bond position (S8); transition steepness by bond length and number (S9); FRET pair design (S13); example interface (S14); DNA sequences.', '']
mdf = os.path.join(HERE, 'SI_updates.md'); open(mdf, 'w').write('\n'.join(md))
out = os.path.join(HERE, 'SI_updates.docx')
subprocess.run(['pandoc', mdf, '-o', out, '--resource-path', HERE], check=True)
import zipfile                                   # declare png by extension (pandoc only adds per-part overrides)
tmp = out + '.tmp'
with zipfile.ZipFile(out) as zi, zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as zo:
    for it in zi.infolist():
        if it.filename.endswith('/'): continue
        data = zi.read(it.filename)
        if it.filename == '[Content_Types].xml' and b'Extension="png"' not in data:
            data = data.replace(b'<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">',
                                b'<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="png" ContentType="image/png"/>', 1)
        zo.writestr(it, data)
os.replace(tmp, out)
print('wrote', out, f'({len(F)} figures, Tables S1 and S4)')
