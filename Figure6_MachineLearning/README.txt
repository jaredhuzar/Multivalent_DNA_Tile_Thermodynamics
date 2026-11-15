Figure 6 -- machine learning vs the physical model (updated 2026-11-12: comparator = multivalent model of Figure 4).
Generator: make_figure6.py (run after Figure4_MultivalentModel/fit_model.py).  Same 236 interfaces and the same
10-fold folds (KFold shuffle, random_state 0) for every model.  Manuscript letters -> folders:
  6A  definitions (artwork, not generated)
  6B  XGBoost out-of-fold predicted vs measured Tm          Fig6b_XGBoost_StabilityScatter/
  6C  held-out Spearman rho + RMS of six models             Fig6d_ModelComparison_SpearmanRho/  (historical "6d")
  6D  SHAP values (unchanged)                               ../Supplementary/S6_SHAP_Values/
  6E  edge-site effect vs spacing (+ multivalent model)     Fig6e_EdgeSites_Spacing/
  combined preview (A, D placeholders)                      figure6_combined.svg
  Fig6c_XGBoost_RankRank/ (rank-rank, not in the manuscript) is still exported by Tile_Tile_Master.ipynb.

RESULTS (out-of-fold, n = 236)
  model                        rho    fold mean +- sd   RMS (C)
  Multivalent model (Fig. 4)   0.934  0.924 +- 0.052    1.34
  XGBoost                      0.902  0.895 +- 0.075    1.72
  Random forest                0.843  0.845 +- 0.061    2.10
  Quadratic regression         0.809  0.810 +- 0.104    2.98
  Linear regression            0.659  0.681 +- 0.097    3.25
  Summed NN dG (linear in sum) 0.091  0.112 +- 0.245    4.10
  -> the 3-parameter physical model now ranks AND predicts better than every learner (the former A2 comparator had
     rho 0.933, RMS 1.30 C; see Supplementary/SI_PhenomenologicalModel_A2/).
  Edge sites (four-site 7-nt, n = 109; 58 with an occupied edge site 0 or 10): measured median 39.8 vs 42.2 C
  (Mann-Whitney P = 1e-4); mean gap 15.3 vs 11.3 nm; after regressing Tm on mean gap P = 0.81.  The multivalent
  model (out-of-fold) predicts 40.5 vs 41.7 C (about half of the 2.4 C deficit, P = 0.003); measured - model
  residuals still differ slightly (medians -0.16 vs +0.31 C, P = 0.047), so edge sites are a little weaker than the
  model's distance term alone predicts (consistent with, not proof of, extra edge mobility, cf. oxDNA SI).
