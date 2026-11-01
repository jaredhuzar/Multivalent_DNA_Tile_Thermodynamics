Alpha_PerBond_Efficiency
========================
Per-bond avidity coefficient alpha for each multivalent interface.

Model:
  dG_avg^poly = alpha * dG^mono          (per bond)
  N * dG_avg^poly = dG_N^poly = alpha * N * dG^mono

  EACH BINDING SITE = TWO duplex bonds (add_Pairs in Tile_Tile_Master sums two helices per code
  position).  So N = 2 * (number of non-zero code positions), and dG^mono is taken over the
  INDIVIDUAL duplex bonds (two per site, which differ in sequence).

  dG_avg^poly  = experimental interface dG (37 C) / N_bonds   (realized per-bond dG; N = 2*sites)
  dG^mono      = dG of a single monovalent bond; PRIMARY = median of the interface's individual
                 bond dGs (mean, std, min, max also reported — bonds are heterogeneous).
  alpha        = dG_avg^poly / dG^mono   (alpha<1 => each bond realizes less than its monovalent
                 dG, i.e. an avidity / negative-cooperativity penalty)

Per-bond dG: each bond's energy comes ONLY from its own sequence and its own salt correction --
  dG_bond = calc_unified_nn(sequence)["Delta_G_37 (kcal/mol)"] - 0.114*L*ln(sqrt(ionic))
where ionic = 3.3*0.0114 (same saltCorrection as Tile_Tile_Master) and L is that bond's length.
No coaxial-stacking term, no half-splitting: each bond is an isolated, salt-corrected duplex at 37 C
(kcal/mol), matching the experimental buffer/temperature.  compute_alpha.py prints a provenance
check listing each bond's sequence, NN dG37, salt term, and final dG.

Robust experimental-dG set: Filtered_DF_HS rows with Reliable & high-confidence, >=2 replicates,
aggregated per construct by median, 9nt bonds excluded  ->  213 interfaces.

RESULTS (each site = 2 bonds; each bond = NN(sequence) - full salt correction; 9nt excluded)
  alpha (median monovalent dG): median = 0.388, mean = 0.388, std = 0.113, IQR = [0.328, 0.435]
  alpha (mean monovalent dG):   median = 0.401, mean = 0.384
  alpha_global (through-origin slope of dG_avg^poly vs dG^mono_median): 0.389
  => each bond realizes ~0.39 of its monovalent (isolated, salt-corrected) free energy.
  NOTE: alpha is ~invariant to the 2-bonds-per-site count itself (doubling N halves dG_avg^poly and
  splitting each site halves dG^mono); the shift from ~0.33 to ~0.39 is from using the FULL per-bond
  salt correction and dropping the site-level 9nt stacking, i.e. a cleaner monovalent reference.

FILES
  compute_alpha.py          self-contained code (rebuilds NN per-bond dG from Tile_Tile_Master,
                            computes alpha, writes the csv + svg). Run from Tile_TileResults.
  alpha_per_interface.csv   one row per interface: code, n_sites, N_bonds, n_reps, Exp_dG,
                            dG_avg_poly, dG_mono_{median,mean,std,min,max}, alpha_median, alpha_mean
  alpha_figure.svg          left: dG_avg^poly vs dG^mono (slope = global alpha, with per-interface
                            bond-spread error bars; alpha=1 reference); right: distribution of alpha.
