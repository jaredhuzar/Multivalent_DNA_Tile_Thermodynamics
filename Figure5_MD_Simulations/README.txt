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
    B(T) = EQUAL-WEIGHT mean of the 3 seeds' per-replicate means (each seed averaged first, then averaged
    across seeds), NOT a pool of thinned samples -- pooling weights each seed by its surviving sample count
    (replicates that stay bound keep more samples), which over-weights them and biases Tm upward.
  - Panel b (Fig 5B, 3-bond spacing series): CENTER series only; edge variants (_Edge) are excluded.
  - Panel f = raw valency (spacing-confounded); f_matched = valency at matched adjacent spacing s=4 (preferred).
    At matched spacing there is NO valency effect (Flex 2->3 bond 1.46->1.44, Rigid 1.51->1.51); the apparent
    gain in the raw panel was a spacing confound plus a wrong B=1.5 crossing applied to the 2-bond systems.
  - Error bars = mean +/- s.d. over 3 seeds (11/22/33).

legacy figure_panels/ holds an older manuscript-panel copy (kept in sync); the v2 path above is authoritative.
