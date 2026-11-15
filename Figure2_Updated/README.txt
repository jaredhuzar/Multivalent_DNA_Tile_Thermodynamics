Figure 2 (UPDATED / Tm-anchored reframe)   -- v2, rewritten after review
========================================================================

Replaces the van't Hoff-dG panels of the original Figure2_NegativeCooperativity (old A alpha/EM defs,
B, D).  The van't Hoff dG (Exp_DeltaG) is the free energy of the small cooperative unit that controls
the observed transition, NOT the whole-interface assembly free energy, so the original per-bond
"efficiency" alpha, the measured-vs-additive dG scatter, and the molar EM were confounded.

Generate:   python3 Figure2_Updated/make_uniform_ceff_panel.py     (run FIRST: main panel H + SI width/readout panels)
            python3 Figure2_Updated/make_updated_figure2.py        (panels B-F + combined MAIN svg; run from anywhere)

MAIN figure (Figure2_Updated.svg):  B, C, D (top) / F, H, G-placeholder (bottom).
SI:  E (PanelE_dTm_per_added_site), PanelH_UniformCeff_Width (width test), PanelI_ReadoutWeight.

ONE consistent robust set for EVERY panel
-----------------------------------------
  Filtered_DF_HS rows, Reliable==True & Confidence=='high', median per construct,
  >= 2 replicates, 9-nt bonds excluded.   -> n = 205 constructs (165 7nt, 35 5nt, 5 6nt).
  All panels (including F) draw from this same set; n is printed at run time.

Panels (each sub-folder: figure.svg + data.csv; exact points also in ../Figure_Source_Data/Fig2U_*.csv)
  B  PanelB_StrongestBond       interface Tm vs the NN Tm of its STRONGEST single bond (per bond length).
  C  PanelC_Tm_vs_Sites         Tm vs binding sites, by bond length.
  E  PanelE_dTm_per_added_site  dTm from adding one bond, MEASURED vs the ADDITIVE NN prediction.
  D  PanelD_CooperativeUnit     measured transition width (dx) vs valency, vs the all-or-none width.
  F  PanelF_Avidity             avidity log10(EM / C_partner), C_partner = 4.5 nM (one tile).
  (G external 3D mechanism render, unchanged, not produced here.)
  H  PanelH_BindingCooperativity  MAIN.  Cooperativity free energy of each added bond.  For every add-one-bond pair,
                                C_eff is fitted SEPARATELY to each construct's measured Tm (full multivalent model, no
                                all-or-none) and DDG_coop = -RT ln(C_eff,hi / C_eff,lo) at 37 C (>0: the added bond lowers
                                the effective concentration felt by every bond -> helps less), for docking / intermediate
                                (w=0.57) / reporter readouts.  9->10 (5nt) is a cost under every readout (+0.11 / +0.11 /
                                +0.05 kcal/mol, i.e. C_eff x0.84-0.92; 70-80% of pairs; p<=0.02); 4->5 (7nt) is small and
                                readout-dependent (+0.06 / +0.04 / -0.22); 10->11 (8 pairs, one shared 11-site construct)
                                is shaded.  The penalty is ~0.1 kcal/mol, small next to one bond's dG (several kcal/mol):
                                mild negative cooperativity on a highly avid interface.  data.csv also keeps the
                                equivalent constant-C_eff Tm shortfall (rank-identical, |Spearman| > 0.99); summary.csv.
  SI PanelH_UniformCeff_Width   uniform-C_eff model fitted to each Tm: the two pure readouts (docking /
                                reporter) bracket the measured melt width; all-or-none is too sharp
                                (make_uniform_ceff_panel.py; also writes example_curves.csv).
  SI PanelI_ReadoutWeight       readout weight w that matches each construct's width (one uniform C_eff);
                                w ~0.57, flat across bond length and valency (same script; also writes
                                width_ratio_vs_w.csv).

HONEST SCOPE  -- what each panel does and does NOT show
-------------------------------------------------------
Panel B.  The cleanest, most assumption-light avidity statement: every interface melts far above the
  NN Tm of its single strongest bond (7nt: ~4.6 C best bond vs ~42 C interface, gain +38.5, min +31;
  5nt: ~ -19 C vs ~34 C, gain +53, min +49), so no single bond could hold the two tiles together at
  room temperature.  Uses ONLY standard nearest-neighbor thermodynamics for a 5-7 bp duplex (that
  model's best-validated regime), at the same 9 nM / 4.5 nM-per-tile and salt -- no two-state fit, no
  cooperative-unit, no effective-molarity assumption.  "Strongest" = the bond with the highest own Tm.

Panel C.  Direct measured Tm.  Uses the figure's robust set (>=2 reps, HS QC), so it is NOT identical
  to the original Fig 2b, which groups the broader Filtered_DF_Tm table (236 constructs): e.g. 7nt
  4-site here is n=103, median 40.6 C vs 109 / 40.97 C there.  This is a deliberate consistency choice.

Panel E.  "Each added bond helps less" is MOSTLY EXPECTED even with zero cooperativity (Tm saturates
  as bonds are added, and 5nt bonds are weaker).  So the panel plots measured dTm against the ADDITIVE
  nearest-neighbor prediction (grey).  Result: 4->5 realizes ~95% of additivity (NO extra cooperativity
  there); only 9->10 falls short (~12%) -- the genuine high-density negative-cooperativity signal.
  10->11 (+0.9 C, n=8) exceeds 9->10, so the effect is non-monotonic, not a simple "helps less".

Panel D.  The van't Hoff 'Enthalpy' column (ln K vs 1/T fit) is NOISE-BIASED here: errors in theta near
  the 0.15/0.85 window edges are amplified on the log scale and flatten the slope (on synthetic melts at
  the measured plate noise it returns 111 for a true 147 kcal/mol, and 50 at twice the noise), it barely
  tracks the fitted width (per-curve Spearman(|Enthalpy|,dx) ~ -0.08), and only 0.22 of its variance is
  real construct differences.  So Panel D does NOT use that column.  It plots the directly-measured
  transition width dx, and derives enthalpy from the WIDTH: dH_vH = 1.5*R*Tm^2/dx (bimolecular), a
  validated estimator (recovers 145 for a true 147 kcal/mol; 13% replicate scatter; 0.55 variance from
  signal).  The all-or-none reference is the matching bimolecular width dx = 1.5*R*Tm^2/ΔH_NN (~0.9 C).
  Result: observed ~2.1 C is ~2.4x broader than all-or-none -> the interface does NOT melt as one unit;
  the cooperative unit is ~3 helices at 4-5 sites rising to ~7 at 9-10 sites, a FRACTION of the interface
  (< 2N).  Molecularity caveat: a unimolecular sub-transition (unit melting while tiles stay docked)
  would use R*Tm^2/dx, a factor 1.5 smaller -- the >=2.4x-broad conclusion holds, the helix count shifts.
  Reporter caveat (not excluded by the data): the FRET dye sits on the central bond and one site is
  exactly 2 duplexes, so a small apparent unit could partly reflect what the reporter senses locally.
  data.csv reports dx, the width-derived dH_vH, the all-or-none dx, the ratio, and coop_unit_helices.

Panel F.  avidity = log10(EM / C_partner), C_partner = 4.5 nM (one tile) -- the bulk concentration the
  bond's partner would otherwise see; this matches the original Fig 2e convention (dividing by 9 nM
  total gives ~0.3 lower).  EM is derived from dG_total via the measured-minus-additive shortfall, so it
  inherits the dG_total assumption below and rests on attributing that shortfall to effective molarity.

Panel H.  Tests the effective-molarity picture WITHOUT the all-or-none assumption.  Model: once the tiles
  dock, each duplex binds independently at one uniform local concentration C_eff; summing every bound
  microstate, K_d = (1/C_eff)*[prod_i(1 + K_i*C_eff) - 1], with K_i from SantaLucia-1998 NN (salt-corrected)
  and mass balance at 4.5 nM per tile.  The summed-shortfall C_eff formula is the all-or-none limit of this.
  C_eff (T-independent) is fitted per construct to the measured Tm; the transition width is then a free
  prediction.  Two readout limits: docking (theta = fraction docked) and reporter (theta = docked AND the
  two central-site duplexes closed).  Result (n=205): predicted/observed width median 0.61 (docking, too
  sharp) and 2.7 (reporter, too broad); all-or-none 0.37.  So the two PURE readouts bracket the measured
  width and an all-or-none melt is excluded.  This does NOT by itself show that bonds are non-equivalent:
  see Panel I -- an intermediate readout with ONE uniform C_eff reproduces the widths.
  The fitted C_eff itself (log10(C_eff/4.5 nM) ~ 5.8 docking, ~6.4 reporter) supports enhancement >> 1.
  Assumptions: T-independent C_eff; reporter = the two duplexes at the central site; isolated-duplex NN
  energies (no coaxial stacking with the tile helices).  9-nt excluded; A/B designs excluded (no per-site
  NN sequence model).  Left panel: one representative 7nt 4-site construct (raw replicates normalized with
  the pipeline's droop-trimmed Boltzmann fit).

Panel I.  Readout weight.  The real FRET readout is likely between H's two limits: when the central
  reporter duplexes fray but the tiles stay docked, dye and quencher are still held within a few nm and
  stay partly quenched.  theta = f*(w + (1-w)*p_reporter), w=0 reporter, w=1 docking; for each w, C_eff is
  refit to Tm, and w_best is where predicted width = observed width.  Result: median w ~0.57, flat across
  bond length (7nt/5nt) and valency (within-length Spearman ~0); a few constructs never cross (shown open,
  at a bound).  Interpretation: the widths are CONSISTENT with one uniform C_eff if the reporter keeps
  ~57% of its quenching when frayed-but-docked -- a testable prediction (dye/quencher geometry, or a
  docking-only FRET control with the pair away from the bonds).  The scatter of w tracks each construct's
  measured width (Spearman ~ -0.76) and is of the size expected from ~13% replicate scatter in dx, so it
  cannot separate measurement noise from real bond non-equivalence.

Definitions
-----------
  dG_total(37C) = ΔH_NN*(1 - T/Tm) + R*T*ln(C_T/4)
      State-function whole-interface dG.  Because dG is path-independent, this does NOT require the melt
      to be two-state in SHAPE.  It DOES assume ΔH_total = ΔH_NN (Predicted_deltaH; all bonds formed, no
      strain reduction) and dCp ~ 0.  Using the (unreliable) measured van't Hoff ΔH instead would shift
      the median dG_total from about -17 to -13 kcal/mol, so the number is sensitive to this choice.
      R*T*ln(C_T/4) = -12.27 kcal/mol is a CONSTANT standard-state offset (identical for every construct,
      fixed by defining Tm at theta = 1/2); it is NOT a per-construct "docking/translational cost" and
      makes up ~70% of dG_total's magnitude, so construct-to-construct variation comes from the ΔH term.

  width dx is the Boltzmann logistic width from the droop-corrected fit (robust; used directly in D).
  van't Hoff ΔH is fit with the correct BIMOLECULAR, non-self-complementary K = 2*theta/((1-theta)^2*C_T),
  C_T = 9 nM total = two tiles at 4.5 nM each -- but see Panel D: that column is unreliable here.

Data source:  Filtered_DF_HS.csv  and  Extras/Effective_Molarity/effective_molarity.csv (run first,
  for the per-bond NN reference in F).  9-nt bonds excluded throughout (unreliable predicted thermo).

==== MAIN FIGURE LAYOUT v3 (2026-11-12) -- Figure2_Updated.svg ====
  A  interface + FRET readout schematic         (placeholder; existing artwork)
  B  interface Tm vs NN Tm of its strongest single bond   PanelB_StrongestBond   (+38 C 7nt, +53 C 5nt)
  C  cooperative-unit FRACTION f = dH_vH(width)/dH_NN      PanelC_CoopFraction    (medians 0.43/0.30/0.37/0.30)
  D  Tm vs binding sites, by bond length                   PanelC_Tm_vs_Sites (folder keeps old letter)
  E  avidity log10(C_eff/C_partner), C_eff = multivalent-model fit per interface (docking readout)
                                                           PanelF_Avidity         (median 10^5.79, C_eff 2.8 mM; n=205)
  F (standalone, docking only, 4->5 and 9->10)             PanelF_DDGcoop_Docking
  F  DDG_coop per added bond, DOCKING readout only        PanelH_BindingCooperativity
     (medians +0.06 / +0.11 kcal/mol at 4->5 / 9->10; 10->11 not shown (8 pairs, one construct); all readouts
      and regimes stay in
      data.csv and in the standalone figure.svg written by make_uniform_ceff_panel.py)
  G  mechanism schematic                                   (placeholder; existing artwork)
Panel letters are set by the LET dict in make_updated_figure2.py; folder names keep their historical letters.
Moved to SI: PanelD_CooperativeUnit (helix count), PanelE_dTm_per_added_site, PanelH_UniformCeff_Width,
PanelI_ReadoutWeight.

SI PanelSI_AddedBondOccupancy (make_added_bond_occupancy.py; after make_uniform_ceff_panel.py), docking readout:
  added duplexes in 4->5 (7-nt; 432 duplexes / 216 pairs) vs 9->10 (5-nt; 80 / 40), each at the C_eff fitted to the
  lower construct: T* - Tm mean +8.7 C (median +11.8) vs -2.0 C (median +1.9); fraction formed at Tm mean 0.76
  (median 0.93) vs 0.47 (median 0.58); dTm at constant C_eff median +3.2 vs +0.9 C; measured +2.9 vs +0.1 C.
  Both distributions are wide and bimodal (sequence-dependent): the contrast is in the typical added bond, not a
  property of bond length alone.
