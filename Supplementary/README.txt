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
