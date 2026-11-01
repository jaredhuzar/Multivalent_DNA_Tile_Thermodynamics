Alpha_Geometry_Controlled
=========================
Geometry-controlled regression of the avidity coefficient alpha (from Alpha_PerBond_Efficiency),
isolating each structural factor's INDEPENDENT effect -- same method as the Tm geometry analysis.

Set: robust-dG alpha interfaces, 9nt excluded, n = 213.
Response: alpha_median (per-interface avidity = dG_avg^poly / dG^mono).
Predictors: number of bonds, bond span (nm), mean bond length (nt), GC fraction.

OLS (standardized predictors), R^2 = 0.817:
  number of bonds : beta=-0.114  t=-14.5  p~0      (VIF 5.5)
  bond span (nm)  : beta=-0.047  t=-10.6  p~0      (VIF 1.7)
  bond length (nt): beta=-0.065  t= -8.3  p=1.3e-14 (VIF 5.5)
  GC fraction     : beta=-0.014  t= -3.3  p=1.3e-3  (VIF 1.5)

Partial Spearman (controlling the other three):
  number of bonds : marginal -0.79 -> partial -0.67 ***
  bond span       : marginal -0.74 -> partial -0.55 ***
  bond length     : marginal +0.11 -> partial -0.40 ***   (marginal was +; partial is - once
                     GC and span are controlled: longer/stronger bonds underdeliver more)
  GC fraction     : marginal -0.13 -> partial -0.12  ns

INTERPRETATION
  - Number of bonds is the dominant independent driver of reduced avidity (crowding / saturation).
    CAVEAT: alpha = dG_exp/(N * dG^mono) has 1/N built in, and dG_exp saturates with N, so part of
    this is mechanical -- it is the quantification of diminishing returns, not a new effect.
  - Bond span independently lowers avidity (partial rho=-0.55): at fixed bond count, spreading bonds
    out reduces per-bond efficiency (colocalization / effective-concentration effect).
  - Bond length independently lowers avidity (partial rho=-0.40): bonds with a larger monovalent
    expectation realize a smaller fraction of it.
  - GC has no independent effect.

FILES
  alpha_geometry_ols.py                       self-contained code (OLS + partial corr + plots).
  alpha_geometry_data.csv                     per-interface: code, alpha, predictors, and the
                                              residual coords (xres_*, yres_alpha_vs_*) for each panel.
  alpha_geometry_partial_correlations.csv     marginal & partial rho per factor.
  alpha_geometry_controlled.svg               added-variable plots: alpha-residual vs each factor,
                                              controlling the other three.
