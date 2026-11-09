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
echo "=== Figure 3 ==="
run Extras/HS_vs_Distance_7nt_4bonds/hs_vs_distance.py
run Extras/HS_vs_Distance_7nt_4bonds/linear_vs_quadratic.py
echo "=== Figure 4 (phenomenological model: Fig 4b pred-vs-exp, model performances, param CIs) ==="
run_nb TmPredictor_NewTm.ipynb
run Figure4_PhenomenologicalModel/Theory_Scripts/entropy_penalty_vs_spacing.py
run Figure4_PhenomenologicalModel/Theory_Scripts/entropy_penalty_vs_numbonds.py
run Figure4_PhenomenologicalModel/Theory_Scripts/theory_fracbound_noocc.py
echo "=== Figure 6 ==="
run Supplementary/S6_SHAP_Values/recompute_shap_fulltrain.py
echo "=== Supplementary ==="
run Extras/CooperativeUnit_vs_BondCount/coop_unit.py
run Supplementary/Alpha_Geometry_Controlled/alpha_geometry_ols.py
run Supplementary/add_ols_and_meangap.py
run Figure3_BondClustering/Fig3b_Tm_vs_BondSpan_OLSResiduals/make_fig3b.py
run Supplementary/SI_Figures_regenerate.py
run Supplementary/TableS1_Model_Comparison_All/models_and_pdf.py
run Supplementary/TableS1_Model_Comparison_All/pdf_report.py
run Supplementary/TableS1_Model_Comparison_All/make_table_docx.py
echo "=== Figure_Source_Data: exact plotted points per panel ==="
run Figure_Source_Data/make_figure_points.py
echo "-------------------------------------"
echo "PASS=$pass FAIL=$fail"
[ $fail -gt 0 ] && { printf 'FAILED: %s\n' "${failed[@]}"; exit 1; }
echo "all panels regenerated (see $LOG)"
