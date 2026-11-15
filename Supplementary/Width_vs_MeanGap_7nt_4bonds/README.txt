SI Fig. S10 (rebuilt 2026-11-12): melt shape vs mean gap, four-site 7-nt interfaces (n = 109).
Replaces HS_vs_MeanGap_7nt_4bonds (apparent two-state dH/dS vs mean gap; both flat, rho -0.06/-0.07, because the
fitted dH and dS are correlated at rho ~ 0.999 and individually poorly determined).  In a two-state description
dS_vH = dH_vH/Tm + R ln(4/C_T), so only the width carries shape information beyond Tm.

Code: make_width_vs_meangap.py  (run from anywhere; --png PATH for a preview).  Outputs figure.svg, data.csv.
Selection: same as Fig 3E (Filtered_DF_HS Reliable, median per construct, 4 sites, no A/B/6/8/9/'55').

Results
  A  |dH_vH| = 1.5 R Tm^2 / dx rises with mean gap: Spearman rho = +0.46, P = 6e-7.
  B  cooperative fraction |dH_vH|/|dH_NN|: 0.38 (clustered third, <= 12 nm) -> 0.49 (dispersed third, >= 16 nm);
     Tm 42.2 -> 38.6 C; width 2.34 -> 1.81 C.
  C  width vs mean gap, measured rho = -0.48 (melt sharpens).  Predictions:
       multivalent model, C_eff = C_ref exp(-d/lambda) (mean distance), intermediate readout, fitted to Tm only:
           rho = -0.48, lambda = 16 nm, Tm RMS 1.22 C
       multivalent model, one uniform C_eff: rho = -0.16 (P = 0.09), Tm RMS 2.13 C  -> the Tm^2 dependence of the
           width alone does not produce the trend; the distance-dependent C_eff does
       Fig 4 model A2 (two-state, dH_eff = sum g_k dH_k, lambda = 90 nm): rho = +0.97 -- predicts BROADER melts with
           separation (down-weighting the enthalpy of distant bonds), the opposite of the data.
Interpretation: an elastic loss of enthalpy (or A2's enthalpy down-weighting) predicts broadening; a lower effective
concentration for dispersed bonds requires more bonds to form together, so the melt sharpens, as observed.
