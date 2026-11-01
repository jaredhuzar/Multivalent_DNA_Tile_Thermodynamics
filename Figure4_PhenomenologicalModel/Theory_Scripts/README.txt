Theory_Scripts — generators for the theory-curve CSVs/SVGs (run from Tile_TileResults).

  entropy_penalty_vs_spacing.py   -> EntropyPenalty_vs_Spacing_4bonds.csv + .svg
      dG and its decomposition vs bond spacing (4-bond, 7nt), from the occupancy model
      (Nscale: lam=139.82, dS_int=-8.131, base=5.237). Includes dG_no_occ_weighting (g=1).
  theory_fracbound_noocc.py       -> TheoreticalBondSpacing_FracBound_noOccWeight.svg
      Fraction-bound vs T at several spacings (4-bond) PLUS the g=1 (no occupancy) baseline line.

Other CSV generators already live with their outputs:
  Paper_Figures/Alpha_PerBond_Efficiency/compute_alpha.py        -> alpha_per_interface.csv
  Paper_Figures/Effective_Molarity/effective_molarity.py         -> effective_molarity.csv
  Paper_Figures/Avidity/{plot_avidity.py,em_avidity_report.py}   -> avidity_*.csv
  Paper_Figures/SI_Figures/SI_Figures_regenerate.py              -> S2..S7 data.csv
  Paper_Figures/SI_Figures/add_ols_and_meangap.py                -> OLS + mean-gap data.csv
  Paper_Figures/SI_Figures/Model_Comparison_All/unified_models.py-> model_comparison_all_235.csv
