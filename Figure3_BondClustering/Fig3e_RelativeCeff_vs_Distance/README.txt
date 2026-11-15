Fig 3E (rebuilt from Tm; replaces Fig3d_RelativeOccupancy_vs_Distance).   Generate: python3 make_fig3e.py

Relative effective concentration of a bond, phi(d) = C_eff(d)/C_eff(6 nm), vs the distance d to its nearest bonded
site, for the 109 four-site 7-nt interfaces of Fig 3C.  NO van't Hoff dG is used.

Model (Figure2_Updated Panel H multivalent model + one factor per distance bin):
  C_i = C_ref * phi[bin(i)];  K_a = (1/C_ref)[prod_i(1 + K_i C_i) - 1];  K_i = per-duplex salt-corrected NN;
  mass balance at 4.5 nM per tile; FRET readout theta = f*(w + (1-w)*p_reporter).
  5 global parameters (log10 C_ref, phi at 12/18/24/>=30 nm) least-squares fitted to median Tm.
  phi_d is the relative effective association constant of a bond (what the old exp(-(beta_d-beta_d1)/RT) estimated
  from dG_exp); it equals the relative occupancy only when bonds are far from saturation.

Result (intermediate readout w = 0.57; docking in brackets):
  phi = 1, 0.61 [0.58], 0.33 [0.36], 0.16 [0.20], 0.06 [0.08] at 6, 12, 18, 24, >=30 nm
  ~ exp(-(d - 6 nm)/lambda), lambda ~ 9-10 nm.   Tm RMS 1.37 C vs 2.13 C for one uniform C_eff.
  95% CI: Gauss-Newton (t, df = 104); bootstrap over constructs (200x) for w = 0.57 agrees (data.csv boot_ci_*).
  Reporter readout (w = 0; in data.csv only, its widths are ~2.7x too broad) gives a similar decay with a
  flatter 12-nm point (0.97) and 0.02 at >=30 nm.

Compared with the old dG-based panel (0.91, 0.80, 0.47, 0.31): the decay is steeper.  The van't Hoff dG is the
free energy of the cooperative unit (Figure 2D) and compresses per-bond differences, so it understated the loss.

Assumptions: T-independent C_eff; isolated-duplex NN energies; bins by nearest occupied site (both duplexes of a
site share the bin); a single bond's own bulk association is scaled by phi too (same form as the joint fit).
Outputs: figure.svg, data.csv (all three readouts), ../../Figure_Source_Data/Fig3E_RelativeCeff_vs_Distance.csv.

Descriptor caveat: nearest-bond and mean distance to the other bonds are highly correlated here (r = 0.93).
With the same model, an exponential in MEAN distance fits these Tm better (RMS 1.21 C, lambda ~14 nm) than one in
nearest-bond distance (1.37 C, ~10 nm); fitting both drives the nearest-bond lambda to ~76 nm.  The Fig 4 model
agrees (10-fold CV RMS 1.30 C mean vs 1.64 C nearest-bond).  So this panel describes the decay along the
nearest-bond axis; it does not show that the nearest neighbour is the controlling distance.
