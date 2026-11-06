Fig 6e. SHAP values of the XGBoost melting-temperature model.

DATA / MODEL
  Target   : Tm (°C), per-construct median.
  Samples  : the 235 dimer interfaces of the Tm table (Filtered_DF_Tm, reliability-gated, 9nt excluded).
  Features : 11 inputs, one per edge position (Helix 0..10). Each = the salt-corrected nearest-neighbor
             predicted dG (G_corrected) of the bond at that position (0 if no bond). From make_Energy_Array.
  Model    : XGBRegressor(n_estimators=500, lr=0.05, max_depth=4, subsample=0.8, colsample_bytree=0.8,
             tree_method='hist', random_state=0), trained on ALL n=235 samples (not a single fold).
  SHAP     : shap.TreeExplainer on the full-data model; beeswarm fixed in Helix 0 (top) -> Helix 10 (bottom)
             order (NOT importance-sorted). Color = helix feature value (low=blue .. high=red).

FILES
  figure.svg                  beeswarm (Helix-ordered).
  data.csv                    per-sample SHAP_helix0..10 and featval_helix0..10 (235 rows).
  recompute_shap_fulltrain.py self-contained: rebuilds X/y, trains on all 235, recomputes SHAP, redraws.

NOTE: Helix 5 (the always-present central reporter) has near-zero mean|SHAP| (~0.05) because it is
constant across interfaces, so it carries no predictive variance. In-sample fit is high (rho~0.98) by
construction (full-data fit, no holdout); SHAP here explains the trained model, not out-of-sample skill
(use Fig 6b/6c for predictive performance).
