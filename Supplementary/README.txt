SI figure regeneration (S2-S5, S7)
==================================
Regenerated from the CURRENT analysis pipeline, source table Filtered_DF_std.csv
(post droop-corrected fit + std reliability filter; still retains 9nt).
Per-construct medians are the plotted observations; barplot error bars = 95% bootstrap CI.

DECISIONS (per user): current pipeline; 9nt shown in S5 only (excluded from S3/S4); S6/S8 separate.

S2_BondLengthCount_Steepness  Steepness vs Num_Bonds, panels 5nt (9,10) and 7nt (4,5).  Reliable set.
S3_BondPosition_Entropy       Entropy (kcal/mol/K) of all interfaces containing a bond at each
                              position 0-10; red dashed = overall mean (= mean at the always-present
                              center position 5). Physical-dH/dS set, 9nt excluded.
S4_BondPosition_Melting       Same, Tm (deg C). Reliable set, 9nt excluded.
S5_BondLength_Entropy         Entropy (kcal/mol/K) vs bond length 5/7/9. Physical set, 9nt INCLUDED.
S7_PredictedVsMeasured_Tm     NN-predicted Tm vs measured Tm (per construct), y=x line. 9nt excluded.

NOTES
- Entropy is shown in kcal/mol/K (native column is cal/mol/K, divided by 1000) to match the SI axis.
- Values differ somewhat from the printed SI because the fitting pipeline changed since (droop
  correction, reliability/confidence gates, per-bond salt fix). Trends/shapes are preserved.
- S6 (spacer vs bond length) was built from a bespoke in-memory frame (dfA2) not saved in the repo;
  regenerate once the 5 matched 7nt/5nt edge pairs are specified.  S8 is a separate model (TmPredictor).

Code: ../SI_Figures_regenerate.py

ADDED (OLS regressions + mean-gap HS)
-------------------------------------
OLS_Tm_geometry/      geometry-controlled OLS, response = Tm (n=250, 9nt excl).
OLS_dG_geometry/      geometry-controlled OLS, response = experimental dG (n=245).
OLS_alpha_geometry/   geometry-controlled OLS, response = avidity alpha (n=213).
   Each: figure.svg (added-variable panels, alpha-residual vs factor-residual controlling the other
   three), partial_correlations.csv (marginal vs partial Spearman), ols_coefficients.csv
   (standardized beta, SE, t, p, VIF, with R^2 in the figure title), data.csv (residual coords).
   Predictors (all three): number of bonds, bond span (nm), bond length (nt), GC fraction.
HS_vs_MeanGap_7nt_4bonds/  van't Hoff Enthalpy & Entropy vs MEAN GAP (= span/(k-1), nm) for 7nt /
   4-bond interfaces, high-confidence set (n=105). Both ~flat (dH rho=-0.10, dS rho=-0.12).
Code: ../SI_Figures/add_ols_and_meangap.py
   SUPERSEDED (2026-11-12) for SI Fig. S10 by Width_vs_MeanGap_7nt_4bonds/ (width-derived dH_vH, cooperative
   fraction, and model-predicted widths vs mean gap); see that folder's README.


==== SI UPDATE 2026-11-12 (multivalent-model revision) ====
Assembled document: SI_Updates/SI_updates.docx (make_si_updates.py) -- every changed SI figure with caption and numbers,
Tables S1 and S4, and lists of retired / unchanged items.  Figure numbers below refer to the 11/1/2026 SI draft.
NEW / REPLACEMENT
  Ceff_Geometry_OLS/                  per-interface log10 C_eff vs N, span, bond length, GC (replaces the alpha OLS,
                                      S2_OLS_alpha_geometry/ and Alpha_Geometry_Controlled/; 11/1 Fig S2)
  Width_vs_MeanGap_7nt_4bonds/        width-derived dH_vH, cooperative fraction and model widths vs mean gap
                                      (replaces S4_HS_vs_MeanGap / HS_vs_MeanGap; 11/1 Fig S4, main-text S10)
  BondLength_Ceff/                    bond-type effective concentration vs WLC tether estimate (replaces
                                      S11_BondLength_Entropy; 11/1 Fig S11)
  PredictedVsMeasured_NN_Multivalent/ additive NN (1 M) vs multivalent prediction (replaces S12_PredictedVsMeasured_Tm)
  TableS1_Model_Comparison_All/TableS1_all_models.{csv,docx}  families M (multivalent), A, B; same 236 + folds
                                      (make_tableS1_combined.py; the older phenomenological_models_table.docx and
                                      .pdf cover families A/B only)
  TableS4_Parameter_CIs/TableS4_parameters.{csv,docx}  multivalent (GN + bootstrap) and A2 (jackknife + bootstrap)
  SI_PhenomenologicalModel_A2/        former main-text Figure 4 (two-state model A2 and its theory curves)
  ../Figure2_Updated/Panel{H_UniformCeff_Width, I_ReadoutWeight, D_CooperativeUnit, H_BindingCooperativity}  SI panels
RETIRED (kept on disk, not for the SI): S2_OLS_alpha_geometry/, Alpha_Geometry_Controlled/ (alpha is 1 M-referenced);
  S10_OLS_dG_geometry/ (dG_exp is the van't Hoff dG of the cooperative unit); S4_HS_vs_MeanGap_7nt_4bonds/;
  S11_BondLength_Entropy/; S12_PredictedVsMeasured_Tm/; Figure2_Updated/PanelE_dTm_per_added_site/ (1 M additive ref).
UNCHANGED: AFM (S1), S3_OLS_Tm_geometry (S3), MD (S5), S6_SHAP_Values (S6), oxDNA (S7), S8_BondPosition_Melting (S8),
  S9_BondLengthCount_Steepness (S9), FRET design (S13), scadnano interface (S14).

==== SI FIGURE NUMBERING (2026-11-13) ====
The numbered SI figures (S4, S6-S18) are assembled into ../Pseudo_Figures/SI/ by `python3 make_pseudo_figures.py SI`
(re-runs every generator, then rebuilds the folder); ../Pseudo_Figures/SI/README.md maps each number to its SVG/CSV
files and source.  Generators added for the numbering: NineNt_MeltingCurves/ (S4), TransitionWidth_ByBonds/ (S6),
Width_vs_MeanGap_7nt_4bonds/panels/ (S11 merged + S11a/b, S14), SI_PhenomenologicalModel_A2/make_twostate_pred_rank.py
(S13: prediction + rank), BondLength_Ceff/make_kappa_confound.py (S15), oxDNA/make_rmsf_panel.py (S18);
S12 = Figure2_Updated/PanelI_ReadoutWeight/figure_bestfit_w.svg; S16 = Figure 5 MD panels e/f (Tm only) + melting example;
S17 = S6_SHAP_Values/figure.svg + mean_shap_by_site.svg.
