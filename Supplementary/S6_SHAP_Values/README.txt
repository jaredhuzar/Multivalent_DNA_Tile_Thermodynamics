Fig 6e. SHAP values of the XGBoost melting-temperature model.

DATA / MODEL
  Target   : Tm (°C), per-construct median.
  Samples  : the canonical 236 interfaces (Table S1 / Figures 4 and 6; updated 2026-11-13).
  Features : 11 inputs, one per site (Helix 0..10) = mean salt-corrected NN dG at 37 C of the site's two duplexes
             (0 if empty); identical to Figure 6 (mv_model.site_dG_features).
  Model    : XGBRegressor(n_estimators=500, lr=0.05, max_depth=4, subsample=0.8, colsample_bytree=0.8,
             tree_method='hist', random_state=0), trained on ALL n=236 samples (not a single fold).
  SHAP     : shap.TreeExplainer on the full-data model; beeswarm fixed in Helix 0 (top) -> Helix 10 (bottom)
             order (NOT importance-sorted). Color = helix feature value (low=blue .. high=red).

FILES
  figure.svg                  beeswarm (Helix-ordered).
  data.csv                    per-sample SHAP_helix0..10 and featval_helix0..10 (236 rows).
  recompute_shap_fulltrain.py self-contained: rebuilds X/y, trains on all 236, recomputes SHAP, redraws.

NOTE: Helix 5 (the always-present central reporter) has near-zero mean|SHAP| (~0.05) because it is
constant across interfaces, so it carries no predictive variance. In-sample fit is high (rho~0.98) by
construction (full-data fit, no holdout); SHAP here explains the trained model, not out-of-sample skill
(use Fig 6b/6c for predictive performance).
