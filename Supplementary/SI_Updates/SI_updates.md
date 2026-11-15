# SI updates (2026-11-12)

Figures and tables of the Supplementary Information that change with the multivalent-model revision. Numbers are read from the current repository outputs. Figure numbers refer to the 11/1/2026 SI draft; the main text cites a different numbering, so SI citations in the main text need renumbering once the order below is fixed.

## Figure U1. Effective concentration versus interface geometry (replaces the avidity-α regression; SI Fig. S2 of the 11/1 draft, cited as Fig. S8 in the main text)

![](/Users/jaredhuzar/hoomd-blue/build/Claude/Tile_TileResults/Figure_Data_and_Code/Supplementary/SI_Updates/img/ceff_ols.png){width=6.5in}

**Figure U1.** Each interface's effective concentration was fitted on its own to its measured T_m with the multivalent binding model (one uniform C_eff per interface, intermediate FRET readout) and regressed on the number of occupied sites, bond span, mean bond length and GC fraction. Panels are added-variable plots: log10 C_eff and each factor after regressing out the other three (OLS R² = 0.61, n = 236; median C_eff 4.0 mM). Bond span dominates (partial ρ = -0.70, P = 0); bond length (ρ = +0.05, P = 0.5) and GC fraction (ρ = -0.02, P = 0.7) have no detectable effect; at fixed span, more sites (a denser interface) slightly raise C_eff (ρ = +0.13, P = 0.04). Number of sites and bond length are collinear in this library (VIF ≈ 7–9). The avidity coefficient α is no longer used: it is referenced to the 1 M standard state, so α < 1 and its decline with bond number are expected even when C_eff is constant.

## Figure U2. Transition width versus bond spacing (replaces van’t Hoff ΔH and ΔS versus mean gap; SI Fig. S4 of the 11/1 draft, cited as Fig. S10 in the main text)

![](/Users/jaredhuzar/hoomd-blue/build/Claude/Tile_TileResults/Figure_Data_and_Code/Supplementary/SI_Updates/img/width_gap.png){width=6.5in}

**Figure U2.** Four-site 7-nt interfaces (n = 109). (A) Width-derived van’t Hoff enthalpy |ΔH_vH| = 1.5 RT_m²/δ rises with mean gap (ρ = +0.46, P = 6e-07). (B) The cooperative fraction |ΔH_vH|/|ΔH_NN| rises from 0.38 (most clustered third, ≤12 nm) to 0.49 (most dispersed third, ≥16 nm) while T_m falls from 42.2 to 38.6 °C. (C) Measured width (ρ = -0.48) and predictions fitted to T_m alone: multivalent model with distance-dependent C_eff (ρ = -0.48); one uniform C_eff (ρ = -0.16); two-state model A2 (ρ = +0.97, broadening). In a two-state description ΔS_vH = ΔH_vH/T_m + R ln(4/C_T) is fixed by ΔH_vH and T_m, so only the width carries shape information; fitted two-state ΔH° and ΔS° are correlated (ρ ≈ 0.999) and show no trend (ρ = -0.06, -0.07).

## Figure U3. Bond length and effective concentration (replaces entropic penalty by bond length; SI Fig. S11 of the 11/1 draft)

![](/Users/jaredhuzar/hoomd-blue/build/Claude/Tile_TileResults/Figure_Data_and_Code/Supplementary/SI_Updates/img/bondlength.png){width=6.5in}

**Figure U3.** (A) Effective concentration of 5- and 6-nt bonds relative to 7-nt bonds, fitted as an extra factor φ in the multivalent model (C_eff,i = C_ref φ_type e^(−d_i/λ) e^(−κ(N−4))), with the valency term (circles; φ5 = 0.75, 95% CI 0.57–0.99; φ6 = 0.90) and without it (squares; φ5 = 0.59, 0.52–0.65; φ6 = 0.74). A worm-like-chain estimate of the tether effective concentration (red; Table S1 model B4) predicts the opposite and much larger dependence (φ5 ≈ 10, φ6 ≈ 3.8). (B) Cross-validated RMS error is unchanged (1.34, 1.34 and 1.34 °C). φ5 and κ are strongly correlated (r = 0.92) because 5-nt bonds dominate the many-site interfaces, so a valency cost and slightly weaker 5-nt bonds cannot be separated in this library; the paired comparison of Fig. 2F tests the valency cost within one bond length. Replaces the earlier figure, which plotted the noise-biased two-state van’t Hoff entropy.

## Figure U4. Sequence-predicted versus measured melting temperature (replaces SI Fig. S12 of the 11/1 draft)

![](/Users/jaredhuzar/hoomd-blue/build/Claude/Tile_TileResults/Figure_Data_and_Code/Supplementary/SI_Updates/img/pred_meas.png){width=6.5in}

**Figure U4.** (A) Additive nearest-neighbor two-state T_m, summing all duplexes with each referenced to the 1 M standard state (RMS 24 °C, ρ = 0.36); it predicts the many-site 5-nt interfaces to be the most stable (median 65.4 versus 64.4 °C for 7-nt) although they melt lowest. (B) The multivalent model of Fig. 4, same sequences, out-of-fold (RMS 1.34 °C, ρ = 0.934), n = 236.

## Figure U5. Uniform effective concentration and the transition width (moved from main-text Fig. 2)

![](/Users/jaredhuzar/hoomd-blue/build/Claude/Tile_TileResults/Figure_Data_and_Code/Supplementary/SI_Updates/img/uniform_width.png){width=6.5in}

**Figure U5.** One uniform C_eff per interface, fitted to its measured T_m with the multivalent model, predicts the transition width. With the FRET signal taken to report tile docking, predicted widths are 0.61 of measured (too sharp); with it taken to report the central reporter duplexes, 2.7 (too broad); an all-or-none transition with ΔH = ΣΔH_NN gives 0.37. The two pure readouts bracket the data and an all-or-none melt is excluded (n = 205; left, a representative 7-nt four-site interface).

## Figure U6. FRET readout weight (moved from main-text Fig. 2)

![](/Users/jaredhuzar/hoomd-blue/build/Claude/Tile_TileResults/Figure_Data_and_Code/Supplementary/SI_Updates/img/readout_w.png){width=6.5in}

**Figure U6.** Readout θ = f[w + (1 − w)p_rep]; for each w the uniform C_eff is refit to T_m and w_best is where the predicted width equals the measured width. Median w ≈ 0.57, flat across bond length and number of sites; this value is fixed in Figs. 2–4.

## Figure U7. Cooperative unit in duplexes (moved from main-text Fig. 2; the cooperative fraction is Fig. 2C)

![](/Users/jaredhuzar/hoomd-blue/build/Claude/Tile_TileResults/Figure_Data_and_Code/Supplementary/SI_Updates/img/coop_unit.png){width=6.5in}

**Figure U7.** Number of duplexes that melt as one unit, n_coop = ΔH_vH/(ΔH_NN/2N), from the width-derived van’t Hoff enthalpy. About 3 duplexes at four to five sites rising to 6–7 at nine to ten sites, always fewer than the 2N duplexes of the interface (dashed); the melt is about 2.4 times broader than an all-or-none transition.

## Figure U8. Cooperativity free energy for three FRET readouts (Fig. 2F shows the docking readout)

![](/Users/jaredhuzar/hoomd-blue/build/Claude/Tile_TileResults/Figure_Data_and_Code/Supplementary/SI_Updates/img/ddg_readouts.png){width=6.5in}

**Figure U8.** ΔΔG_coop = −RT ln(C_eff,hi/C_eff,lo) for every pair of interfaces differing by one site, with C_eff fitted separately to each construct, for docking, intermediate (w = 0.57) and reporter readouts. 4→5 (7-nt, n = 216): +0.06, +0.04 and -0.22 kcal mol⁻¹; 9→10 (5-nt, n = 40): +0.11, +0.10 and +0.06 kcal mol⁻¹. 10→11 (shaded) comprises 8 pairs sharing one 11-site construct. Only 9→10 is positive under all three readouts.

## Figure U9. Two-state model A2: predicted versus measured T_m (former main-text Fig. 4B)

![](/Users/jaredhuzar/hoomd-blue/build/Claude/Tile_TileResults/Figure_Data_and_Code/Supplementary/SI_Updates/img/a2_pred.png){width=6.5in}

**Figure U9.** Two-state model T_m = ΣgΔH_NN / (ΣgΔS_NN + NΔS_Bond + ΔS_0 + R ln(C/2)), g = e^(−d/λ), out-of-fold (10-fold CV over the same 236 interfaces): RMS 1.30 °C, ρ = 0.933. It predicts T_m as well as the multivalent model but, being two-state, implies transitions about half as wide as measured and broadening with spacing (Fig. S10). Parameters in Table S4.

## Figure U10. Two-state model A2: theoretical curves (former main-text Fig. 4C–F)

![](/Users/jaredhuzar/hoomd-blue/build/Claude/Tile_TileResults/Figure_Data_and_Code/Supplementary/SI_Updates/img/a2_theory.png){width=6.5in}

**Figure U10.** Enthalpy, entropy and bound fraction at 40 °C versus number of bonds (clustered and spread) and versus spacing, from model A2 with the fitted parameters of Table S4. Superseded in the main text by the multivalent-model curves of Fig. 4D–H.

## Table S1. Tested models of interface T_m

| ID | Model | k | CV RMS (°C) | CV ρ | Fitted parameters (all data) |
|---|---|---|---|---|---|
| M1 | one uniform C_eff | 1 | 2.40 | 0.824 | C_ref = 3.9 mM |
| M2 | + distance | 2 | 1.57 | 0.933 | C_ref = 27.6 mM, λ = 11.6 nm |
| M3 | + valency (Fig. 4 model) | 3 | 1.34 | 0.934 | C_ref = 26.7 mM, λ = 12.8 nm, κ = 0.079 |
| M4 | + bond-type C_eff (no valency) | 4 | 1.34 | 0.935 | C_ref = 30.6 mM, λ = 12.6 nm, φ5 = 0.59, φ6 = 0.74 |
| M5 | + bond-type C_eff + valency | 5 | 1.34 | 0.935 | C_ref = 28.8 mM, λ = 12.7 nm, κ = 0.041, φ5 = 0.75, φ6 = 0.90 |
| M1d | one uniform C_eff [docking readout] | 1 | 2.22 | 0.831 | C_ref = 2.8 mM |
| M2d | + distance [docking readout] | 2 | 1.51 | 0.932 | C_ref = 19.3 mM, λ = 12.6 nm |
| M3d | + valency (Fig. 4 model) [docking readout] | 3 | 1.41 | 0.931 | C_ref = 19.8 mM, λ = 13.1 nm, κ = 0.052 |
| A1 | mean (g-weighted ΔS_Bond) | 3 | 1.33 | 0.930 | λ = 28.14, ΔS_Bond = -10.83, base = 21.91 |
| A2 | N-scaled ΔS_Bond | 3 | 1.30 | 0.933 | λ = 90.09, ΔS_Bond = -8.38, base = 7.31 |
| A3 | g on ΔH only | 3 | 1.31 | 0.934 | λ→∞ (unidentified), ΔS_Bond = -8.25, base = 8.07 |
| A4 | nearest-neighbor g | 3 | 1.98 | 0.840 | λ = 32.74, ΔS_Bond = -11.69, base = 18.4 |
| A5 | per-type ΔS_Bond | 4 | 1.33 | 0.932 | λ = 26.54, ΔS_short = -10.91, ΔS_long = -11.16, base = 23.67 |
| A6 | split λ_H/λ_S | 4 | 1.29 | 0.932 | λ_H = 55.48, λ_S = 57.84, ΔS_Bond = -9.37, base = 15.87 |
| A7 | occ. + rotational (Green's fn) | 4 | 1.27 | 0.937 | λ = 138.55, ΔS_Bond = -8.68, base = 29.09, k_rot = 1.85 |
| B1 | g=1 (no arrangement weighting) | 2 | 2.26 | 0.783 | ΔS_Bond = -11.27, base = 7.54 |
| B2 | g=1 + rotational | 3 | 1.37 | 0.933 | ΔS_Bond = -9.45, base = 69.88, k_rot = 5.32 |
| B3 | g-weighted + WLC tethering | 2 | 3.18 | 0.602 | λ→∞ (unidentified), base = 49.22 |
| B4 | WLC tethering only (derived) | 2 | 3.17 | 0.613 | λ→∞ (unidentified), base = 51.19 |

All models were fitted to the median T_m of the same 236 interfaces and evaluated by 10-fold cross-validation over interfaces with identical folds; CV statistics are out-of-fold (CV R² and model forms in the formatted table). M: multivalent binding models (M3 = Fig. 4; "d" = docking readout). A: two-state, distance-weighted. B: two-state with mechanical-coupling or worm-like-chain tether terms. A formatted version is in Supplementary/TableS1_Model_Comparison_All/TableS1_all_models.docx.

## Table S4. Fitted parameters and uncertainties

| Model | Parameter | Units | Estimate | Bootstrap 95% CI | Meaning |
|---|---|---|---|---|---|
| Multivalent (Fig. 4) | C_ref | mM | 26.7 | 22.0 – 31.3 | effective concentration of a bond at zero distance in a four-site interface |
| Multivalent (Fig. 4) | λ | nm | 12.8 | 11.8 – 14.3 | decay length of the effective concentration with mean distance |
| Multivalent (Fig. 4) | κ | per site | 0.0791 | 0.0652 – 0.0924 | fractional loss of effective concentration per occupied site |
| Multivalent (Fig. 4) | w |  | 0.57 | — | FRET readout weight, set from transition widths (Fig. 2, SI) |
| Two-state A2 (SI) | λ | nm | 90.1 | 82.6 – 99.1 | decay length of the two-state weight g = e^(−d/λ) |
| Two-state A2 (SI) | ΔS_Bond | cal mol⁻¹ K⁻¹ | -8.38 | -8.63 – -8.14 | per-bond entropy (≈ R ln C_eff) |
| Two-state A2 (SI) | ΔS_0 | cal mol⁻¹ K⁻¹ | 7.31 | 5.69 – 8.87 | dimerization entropy offset |

Standard errors are in the formatted table. Multivalent model: construct bootstrap (200 resamples). Two-state A2: jackknife SE and bootstrap (1000). C_ref SEs are in log10 units. Formatted version: Supplementary/TableS4_Parameter_CIs/TableS4_parameters.docx.

## Retired SI items

- Avidity α regression (11/1 SI Fig. S2): α is referenced to 1 M; replaced by Figure U1.
- Regression of the experimental free energy ΔG_exp on geometry (11/1 SI Fig. S10): ΔG_exp is the noise-biased van’t Hoff free energy of the cooperative unit; removed (Figure U1 covers the geometry dependence).
- Van’t Hoff ΔH and ΔS versus mean gap (11/1 SI Fig. S4): replaced by Figure U2.
- Entropic penalty by bond length (11/1 SI Fig. S11): replaced by Figure U3.
- NN-predicted versus measured T_m alone (11/1 SI Fig. S12): replaced by Figure U4.
- ΔT_m per added site versus the additive NN prediction (former Fig. 2E draft panel): the additive prediction is referenced to 1 M; superseded by Fig. 2F and Figure U8.
- Old Table S1 (two-state models only, 11/1 numbers): replaced by the Table S1 above.

## Unchanged SI items

- AFM of tile dimers (S1); geometry-controlled regression of T_m (S3); coarse-grained simulations (S5); SHAP values (S6); oxDNA fluctuations (S7); T_m by bond position (S8); transition steepness by bond length and number (S9); FRET pair design (S13); example interface (S14); DNA sequences.
