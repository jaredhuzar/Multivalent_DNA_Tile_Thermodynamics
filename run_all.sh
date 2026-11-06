#!/usr/bin/env bash
# Regenerate every script-produced figure panel, in dependency order.
# Each script auto-detects the repo root and writes back into its own folder.
cd "$(dirname "$0")"
LOG=run_all.log; : > "$LOG"
pass=0; fail=0; failed=()
run(){
  printf '%-62s ' "$1"
  if python3 "$1" >>"$LOG" 2>&1; then echo "OK"; pass=$((pass+1));
  else echo "FAIL"; fail=$((fail+1)); failed+=("$1"); fi
}
echo "Renegerate Tables"
jupyter execute Tile_Tile_MasterRegenerate.ipynb
echo "=== Figure 2 ==="
run Figure2_NegativeCooperativity/NNdG_ExpdG_saltfixed/regen_nndg_vs_exp.py
run Figure2_NegativeCooperativity/Fig2d_ddG/regen_fig2d_ddg.py
run Figure2_NegativeCooperativity/Alpha_PerBond_Efficiency/compute_alpha.py
run Figure2_NegativeCooperativity/Effective_Molarity/effective_molarity.py
run Figure2_NegativeCooperativity/Avidity/plot_avidity.py
run Figure2_NegativeCooperativity/Avidity/em_avidity_report.py
echo "=== Figure 3 ==="
run Figure3_BondClustering/HS_vs_Distance_7nt_4bonds/hs_vs_distance.py
run Figure3_BondClustering/HS_vs_Distance_7nt_4bonds/linear_vs_quadratic.py
echo "=== Figure 4 ==="
run Figure4_PhenomenologicalModel/Theory_Scripts/entropy_penalty_vs_spacing.py
run Figure4_PhenomenologicalModel/Theory_Scripts/entropy_penalty_vs_numbonds.py
run Figure4_PhenomenologicalModel/Theory_Scripts/theory_fracbound_noocc.py
echo "=== Figure 6 ==="
run Figure6_MachineLearning/Fig6e_SHAP_Values/recompute_shap_fulltrain.py
echo "=== Supplementary ==="
run Supplementary/CooperativeUnit_vs_BondCount/coop_unit.py
run Supplementary/Alpha_Geometry_Controlled/alpha_geometry_ols.py
run Supplementary/add_ols_and_meangap.py
run Supplementary/SI_Figures_regenerate.py
run Supplementary/Model_Comparison_All/models_and_pdf.py
run Supplementary/Model_Comparison_All/pdf_report.py
run Supplementary/Model_Comparison_All/make_table_docx.py
echo "-------------------------------------"
echo "PASS=$pass FAIL=$fail"
[ $fail -gt 0 ] && { printf 'FAILED: %s\n' "${failed[@]}"; exit 1; }
echo "all panels regenerated (see $LOG)"
