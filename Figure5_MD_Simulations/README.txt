Figure 5 — Coarse-grained MD (HOOMD-blue). Self-contained, canonical path:

  Figure_AH_SelfContained/
    Analysis_Tm_trends.ipynb   generator for every Fig 5 / SI S11 panel
    std_results_v2/*.npz        processed per-run data (nbound, potential, temperatures)
    plots_robust_v2/{b..f,f_matched}/   Tm panels (CANONICAL output)
    plots_robust/mechanism/     E-H: tip-gap histograms, registration, chain x rod heatmaps

plots_robust_v2/ is mirrored up to:  Figure5_MD_Simulations/plots_robust_v2/  and
../Figure_Source_Data/plots_robust_v2/   (and the exact plotted points: ../Figure_Source_Data/Fig5*.csv)

Conventions:
  - Simulated Tm = temperature where mean formed bonds crosses n_max/2  (1.0 for 2-bond, 1.5 for 3-bond).
  - Panel b (Fig 5B, 3-bond spacing series): CENTER series only; edge variants (_Edge) are excluded.
  - Panel f = raw valency (spacing-confounded); f_matched = valency at matched adjacent spacing s=4 (preferred).
  - Error bars = mean +/- s.d. over 3 seeds (11/22/33).

legacy figure_panels/ holds an older manuscript-panel copy (kept in sync); the v2 path above is authoritative.
