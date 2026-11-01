Model comparison on the IDENTICAL 235-interface set (Filtered_DF_Tm medians).
All variants of both families refit and 10-fold cross-validated on the same data.

SHARED CORE (two-state, closed form):
  Tm(C) = 1000*dH_eff / dS_eff - 273.15
  dH_eff = sum_k g_k * dH_NN,k            (kcal/mol)
  g_k    = exp(-dbar_k / lambda)          dbar_k = mean distance of bond k to the others (A4/B use variants)
  every dS_eff includes + R*ln(C_T/2) ~ -39.6 cal/mol/K   (bimolecular term, C_T=4.5 nM)

FAMILY A (occupancy-weighted phenomenological):
  A1 mean      dS = sum g_k (dS_NN,k + dS_Bond) + base
  A2 Nscale    dS = sum g_k dS_NN,k + N*dS_Bond + base          <-- best (RMS 1.98, rho 0.884)
  A3 hOnly     dS = sum dS_NN,k + N*dS_Bond + base ; g on dH only
  A4 nn        A1 with g_k = exp(-dmin_k/lambda)
  A5 pertype   dS = sum g_k (dS_NN,k + dS_Bond[short/long]) + base
  A6 splitlam  dH uses g(lamH), dS uses g(lamS)
  A7 occ+rot   A2 - k_rot*R*ln(1 + sum_k (z_k - zbar)^2)         (elastic Green's-fn term added)

FAMILY B (elastic / mechanical-coupling):
  B1 g=1       dS = sum dS_NN,k + N*dS_Bond + base   (no arrangement weighting)
  B2 g=1+rot   B1 - k_rot*R*ln(1 + inertia)          (rotational Green's fn; best mechanical)
  B3 gw+teth   dS = sum g_k (dS_NN,k + dS_teth[type]) + base   (WLC tethering, derived)
  B4 teth only dS = sum (dS_NN,k + dS_teth[type]) + base

  Rotational (elastic Green's fn):  dS_rot = -k_rot*R*ln(1 + sum_k (z_k-zbar)^2)
  WLC tethering (derived, no fit):  dS_teth[type] = R*ln(V2/(V1^2 * C0)),  V=(4/3)pi<r^2>^1.5,
      <r^2> from a worm-like chain (l_p=1 nm, rise 0.63 nm/nt) at each type's ss states
      SS_STATES = {5:(7,2), 6:(7,1), 7:(8,1), 9:(9,0), A/B:(8,1)}.

RESULTS (10-fold CV, n=235): see model_comparison_all_235.csv.
  Occupancy best A2: RMS 1.98 C, R2 0.805, rho 0.884.
  Best mechanical B2 (rotational, g=1): RMS 2.20, rho 0.860.
  Derived elastic tethering B3/B4: RMS 4.1-6.0, rho 0.54-0.63.
TAKEAWAYS:
  - Occupancy weighting exp(-dbar/lambda) is the single most important ingredient
    (g vs g=1: RMS 2.56 -> 2.08; rho 0.748 -> 0.876).
  - The elastic Green's-function rotational term helps ONLY without g (B1->B2) and is
    degenerate once g is present (A2->A7: no gain, k_rot -> ~0). Arrangement is better
    captured by occupancy weighting than by the mechanical kernel.
  - Physically-derived WLC tethering entropy badly overestimates -> poor fit.
  - NOTE: the previously reported rho=0.795 for the mechanical model came from a different
    (now-missing interfaces.csv) dataset; on the common 235 set the gap is ~0.88 vs ~0.86.

Code: unified_models.py
