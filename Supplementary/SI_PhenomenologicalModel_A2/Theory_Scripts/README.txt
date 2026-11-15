Theory_Scripts — generators for the theory-curve CSVs/SVGs (run from Tile_TileResults).

  entropy_penalty_vs_spacing.py   -> EntropyPenalty_vs_Spacing_4bonds.csv + .svg
      dG and its decomposition vs bond spacing (4-bond, 7nt), from the occupancy model.
      The Nscale parameters (currently lam=90.1 nm, dS_int=-8.38, base=7.31) are READ AT RUNTIME
      from Supplementary/TableS1_Model_Comparison_All/phenomenological_parameter_CIs.csv, so they
      track the fit automatically. Geometry: 6 nm site pitch, 3 nm twin spacing (unified).
      Includes dG_no_occ_weighting (g=1).
  theory_fracbound_noocc.py       -> TheoreticalBondSpacing_FracBound_noOccWeight.svg
      Fraction-bound vs T at several spacings (4-bond) PLUS the g=1 (no occupancy) baseline line.

Other CSV generators already live with their outputs (current repo paths):
  Figure2_NegativeCooperativity/Fig2c_Alpha_PerBond_Efficiency/compute_alpha.py -> alpha_per_interface.csv
  Extras/Effective_Molarity/effective_molarity.py                      -> effective_molarity.csv
  Figure2_NegativeCooperativity/Fig2e_Avidity/{plot_avidity.py,em_avidity_report.py} -> avidity_*.csv
  Supplementary/SI_Figures_regenerate.py                               -> SI data.csv
  Supplementary/add_ols_and_meangap.py                                 -> OLS + mean-gap data.csv
  Supplementary/TableS1_Model_Comparison_All/unified_models.py         -> model_comparison_all_235.csv
