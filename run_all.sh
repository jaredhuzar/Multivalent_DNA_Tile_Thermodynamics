#!/usr/bin/env bash
# Regenerate every table and figure panel from scratch, in dependency order.
# Self-contained: notebooks write panels directly into the organized FigureN/ folders (no copying),
# and each standalone script auto-detects the repo root and writes back into its own folder.
cd "$(dirname "$0")"
LOG=run_all.log; : > "$LOG"
pass=0; fail=0; failed=()
run(){
  printf '%-62s ' "$1"
  if python3 "$1" >>"$LOG" 2>&1; then echo "OK"; pass=$((pass+1));
  else echo "FAIL"; fail=$((fail+1)); failed+=("$1"); fi
}
run_nb(){
  printf '%-62s ' "$1"
  if jupyter execute "$1" >>"$LOG" 2>&1; then echo "OK"; pass=$((pass+1));
  else echo "FAIL"; fail=$((fail+1)); failed+=("$1"); fi
}
echo "=== Tables (notebook) -> Filtered_DF_*.csv ==="
run_nb Tile_Tile_MasterRegenerate.ipynb
echo "=== Authoritative construct table (Tm/dG from Filtered_DF; source for the geometry OLS) ==="
run Experimental_Trends/build_construct_medians.py
echo "=== Notebook figure panels -> Fig2a/2b/2d, 3b/3c/3e, 6b/6c/6d ==="
run_nb Tile_Tile_Master.ipynb
echo "=== Figure 1C: fitted melting curves (svg + raw data) ==="
run Figure1_Platform_AFM/make_melting_curves.py
echo "=== Figure 2 ==="
run Figure2_NegativeCooperativity/Fig2a_NNdG_ExpdG_saltfixed/regen_nndg_vs_exp.py
run Figure2_NegativeCooperativity/Fig2d_ddG/regen_fig2d_ddg.py
run Figure2_NegativeCooperativity/Fig2c_Alpha_PerBond_Efficiency/compute_alpha.py
run Extras/Effective_Molarity/effective_molarity.py
run Figure2_NegativeCooperativity/Fig2e_Avidity/plot_avidity.py
run Figure2_NegativeCooperativity/Fig2e_Avidity/em_avidity_report.py
echo "=== Figure 2 (updated / Tm-anchored reframe: C,E,D,F + combined) ==="
run Figure2_Updated/make_uniform_ceff_panel.py      # Panels H (cooperativity, main) + width test + I (SI)
run Figure2_Updated/make_updated_figure2.py         # Panels B-F + combined main figure (B,C,D,F,H,G)
run Figure2_Updated/make_added_bond_occupancy.py      # SI: why a 5th 7-nt site helps but a 10th 5-nt site barely does
echo "=== Figure 3 ==="
run Extras/HS_vs_Distance_7nt_4bonds/hs_vs_distance.py
run Extras/HS_vs_Distance_7nt_4bonds/linear_vs_quadratic.py
run Figure3_BondClustering/Fig3e_RelativeCeff_vs_Distance/make_fig3e.py   # Fig 3E from Tm (replaces dG-based Fig3d occupancy)
echo "=== SI: former Figure 4 (phenomenological model A2; notebook also writes Fig 6 ML comparison data) ==="
run_nb TmPredictor_NewTm.ipynb
run Supplementary/SI_PhenomenologicalModel_A2/Theory_Scripts/entropy_penalty_vs_spacing.py
run Supplementary/SI_PhenomenologicalModel_A2/Theory_Scripts/entropy_penalty_vs_numbonds.py
run Supplementary/SI_PhenomenologicalModel_A2/Theory_Scripts/theory_fracbound_noocc.py
echo "=== Figure 6 ==="
run Supplementary/S6_SHAP_Values/recompute_shap_fulltrain.py
echo "=== Supplementary ==="
run Extras/CooperativeUnit_vs_BondCount/coop_unit.py
run Supplementary/Alpha_Geometry_Controlled/alpha_geometry_ols.py
run Supplementary/add_ols_and_meangap.py
run Supplementary/Width_vs_MeanGap_7nt_4bonds/make_width_vs_meangap.py   # SI Fig S10 (replaces HS_vs_MeanGap)
run Figure3_BondClustering/Fig3b_Tm_vs_BondSpan_OLSResiduals/make_fig3b.py
run Supplementary/SI_Figures_regenerate.py
run Supplementary/TableS1_Model_Comparison_All/unified_models.py
echo "=== Figure 4 (multivalent model: fit + 10-fold CV + bootstrap CIs, then figure with theory curves; needs Table S1 above) ==="
run Figure4_MultivalentModel/fit_model.py
run Figure4_MultivalentModel/make_figure4.py
echo "=== Figure 6 B, C, E (ML vs multivalent model; after Figure 4) ==="
run Figure6_MachineLearning/make_figure6.py
echo "=== SI figures/tables built on the multivalent model (after Figure 4) ==="
run Supplementary/Ceff_Geometry_OLS/make_ceff_ols.py                       # replaces alpha OLS (S2_OLS_alpha / Alpha_Geometry_Controlled)
run Supplementary/BondLength_Ceff/make_bondlength_ceff.py                  # replaces S11_BondLength_Entropy
run Supplementary/BondLength_Ceff/make_tether_penalty.py                   # per-bond RT ln(1 M/C_eff) by bond length vs WLC
run Supplementary/BondLength_Ceff/make_kappa_confound.py                   # SI Fig S11: kappa vs bond-length confound
run Supplementary/NineNt_MeltingCurves/make_9nt_melting_curves.py         # SI: raw + fitted melts of all 48 9-nt designs
run Supplementary/SI_PhenomenologicalModel_A2/make_twostate_pred_rank.py  # SI: two-state model pred + rank
run Supplementary/TransitionWidth_ByBonds/make_width_by_bonds.py          # SI: width by bond length/count
run Supplementary/oxDNA/make_rmsf_panel.py                                # SI: tile RMSF
run Supplementary/PredictedVsMeasured_NN_Multivalent/make_pred_vs_meas.py  # replaces S12_PredictedVsMeasured_Tm
run Supplementary/TableS1_Model_Comparison_All/make_tableS1_combined.py    # Table S1 (families M, A, B)
run Supplementary/TableS4_Parameter_CIs/make_tableS4.py                    # Table S4
run Supplementary/SI_Tables/make_si_tables.py                              # Tables S1 + S4 CSVs -> Pseudo_Figures/SI (checked vs clean SI docx)
run Supplementary/SI_Tables/make_dataset_s1.py                             # Dataset S1 / Table S3 interface summary (236) -> Pseudo_Figures/SI (checked vs construct_medians)
run Supplementary/SI_Updates/make_si_updates.py                           # SI_updates.docx (changed SI figures + Tables S1/S4)
run Supplementary/TableS1_Model_Comparison_All/pdf_report.py
run Supplementary/TableS1_Model_Comparison_All/make_table_docx.py
echo "=== Figure 5 b-f Tm panels (from processed npz; syncs mirror copies) ==="
run Figure5_MD_Simulations/Figure_AH_SelfContained/regen_fig5_tm_panels.py
echo "=== Figure_Source_Data: exact plotted points per panel ==="
run Figure_Source_Data/make_figure_points.py
echo "=== Pseudo-figures: per-panel SVG + CSV and composite layouts (panels already regenerated above) ==="
printf '%-62s ' make_pseudo_figures.py; if python3 make_pseudo_figures.py --no-regen >>"$LOG" 2>&1; then echo OK; pass=$((pass+1)); else echo FAIL; fail=$((fail+1)); failed+=(make_pseudo_figures.py); fi
echo "-------------------------------------"
echo "PASS=$pass FAIL=$fail"
[ $fail -gt 0 ] && { printf 'FAILED: %s\n' "${failed[@]}"; exit 1; }
echo "all panels regenerated (see $LOG)"
