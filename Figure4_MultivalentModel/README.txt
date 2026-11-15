Figure 4 -- multivalent binding model (rebuilt 2026-11-12).  The former Figure 4 (phenomenological two-state
model A2, Tm = sum g dH / [sum g dS + N dS_bond + dS0 + R ln(C/2)]) is now in
../Supplementary/SI_PhenomenologicalModel_A2/ and stays in Table S1 as a comparator.

MODEL (mv_model.py)
  K_a = (1/C_ref) [ prod_i (1 + K_i C_eff,i) - 1 ]        one docking step, then each duplex binds independently
  C_eff,i = C_ref exp(-d_i/lambda) exp(-kappa (N - 4))     d_i: mean distance of duplex i to the other duplexes
                                                           (6 nm site pitch, twin duplexes +-1.5 nm); N: occupied sites
  K_i: SantaLucia 1998 NN, dS salt-corrected (0.368 (L-1) ln[Na+]eq, [Na+]eq = 3.79 sqrt(0.0113 M))
  mass balance at 4.5 nM per tile; readout theta = f [w + (1-w) p_reporter], w = 0.57 (Figure 2 Panel I), reporter =
  site 5 duplexes; Tm at theta = 1/2.
  Same 236 interfaces, sequences, geometry and 10-fold folds (KFold shuffle, random_state 0) as Table S1.

RESULTS (model_fit/, fit_model.py)
  main model (3 parameters, intermediate readout): CV RMS 1.34 C, Spearman rho 0.934 (out-of-fold)
     log10 C_ref = -1.573 (C_ref 27 mM; bootstrap 95% CI 22-31 mM)
     lambda      = 12.8 nm (11.8-14.3)
     kappa       = 0.079 per site (0.065-0.092): each added site lowers every bond's C_eff by ~8%
  nested (intermediate / docking readout): one uniform C_eff 2.40 / 2.22 C; + distance 1.57 / 1.51 C;
     + valency 1.34 / 1.41 C.   SI phenomenological A2 1.30 C (rho 0.933); B1 (no distance) 2.26 C.
  C_eff at typical spacings is a few mM (C_ref e^{-d/lambda}), consistent with Figure 2E; lambda agrees with the
  Fig 3E decay (9-10 nm nearest-bond, ~14-16 nm mean-distance fits); kappa agrees in sign and size with Fig 2F.
  CAVEATS: w was set from melt widths (Figure 2 Panel I), so absolute width agreement is partly by construction (the
  width TREND with spacing in SI Fig S10 is a prediction).  N and bond length are confounded in the library (4-5 sites
  are 7-nt, 9-10 sites are 5-nt), so kappa can absorb 5- vs 7-nt NN errors; the paired DDG_coop of Fig 2F is the
  cleaner test of the valency cost.  Every construct has a 5-nt reporter site at position 5.

PANELS (make_figure4.py -> figure.svg; theory/*.csv)
  A  equations + fitted parameters (schematic artwork to be added)
  B  out-of-fold predicted vs measured Tm (236)
  C  what each ingredient buys (CV RMS of nested variants; SI A2/B1 for reference)
  D  THEORY Tm vs occupied sites, clustered (all contiguous windows containing site 5; band = min-max over windows),
     7-nt and 5-nt; dashed = kappa 0; filled markers = measured group medians (IQR, all spacings), open = model for
     the same constructs.  Clustered theory sits above the measured medians because the library mixes spacings.
  E  THEORY Tm vs mean gap for all 120 four-site 7-nt arrangements; black = measured medians (n = 109); dashed = one
     uniform C_eff (flat: no spacing dependence without lambda).
  F  THEORY theta(T) for typical arrangements (Tm and width nearest the medians) at mean gap 6 / 12 / 20 nm:
     Tm 43.9 / 41.8 / 38.6 C, width 1.78 / 1.72 / 1.57 C (median width over all arrangements 1.84 -> 1.66 C).
  G  THEORY per-duplex occupancy K C_eff/(1 + K C_eff) at the interface Tm vs mean distance (120 arrangements):
     7-nt ~0.96 -> ~0.88, 5-nt reporter ~0.5 -> ~0.4 with distance; weak 7-nt duplexes form the lower cloud.
  H  THEORY design map: Tm over occupied sites x uniform spacing (7-nt bonds, reporter at site 5), averaged over
     placements on the 11-site edge; blank = does not fit.

Run order: fit_model.py (~2 min + ~0.5 s per bootstrap replicate; --boot 0 to skip), then make_figure4.py (--png PATH).
