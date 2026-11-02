Fig 5 E-H binder-position summaries (HOOMD), copied from plots_robust/mechanism/.
  heatmap_grid.csv          4x4 matrices: per chain_bend (linker) x rod_bend (scaffold) x temperature:
                            tipgap_mean, tipgap_sigma, tipgap_n_unbound, anchorspacing_mean/sigma  (-> G,H)
  registration_summary.csv  per rod_bend x temperature: intra_mean/sigma, dS_reg, transverse_sigma  (-> F)
  registration_raw.csv      per-frame swarm: rod_bend,temperature,rep,metric,value (43k rows)        (-> E,F)
  panelC_histograms_*.csv   tip-gap distributions (raw + binned)                                     (-> E)
  *.svg                     rendered panels.
