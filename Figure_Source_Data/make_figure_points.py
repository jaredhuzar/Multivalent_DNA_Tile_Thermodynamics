import os as _os
def _find_repo():
    d=_os.path.dirname(_os.path.abspath(__file__))
    while d!='/' and not _os.path.exists(_os.path.join(d,'Filtered_DF_Tm.csv')):
        d=_os.path.dirname(d)
    return d
__REPO__=_find_repo()
"""
Figure_Source_Data: one CSV per manuscript panel holding ONLY the exact points that are plotted.

Several panels plot a CONDITIONAL SUBSET of the full dataframe (e.g. Fig 2b/2c use only 7-nt bonds at
4/5 sites and 5-nt bonds at 9/10 sites). This script re-derives each panel's plotted points by reading
the current committed panel data and applying that panel's exact filter, so the folder stays in sync
with any code/data change (it is the last step of run_all.sh). Run from anywhere.
"""
import os, numpy as np, pandas as pd
PROJ=__REPO__; os.chdir(PROJ)
OUT="Figure_Source_Data"; os.makedirs(OUT,exist_ok=True)
F2="Figure2_NegativeCooperativity"; F3="Figure3_BondClustering"
F4="Figure4_PhenomenologicalModel"; F6="Figure6_MachineLearning"
readme=[]
def save(name, df, cond):
    df=df.reset_index(drop=True); df.to_csv(f"{OUT}/{name}.csv", index=False)
    readme.append((name, len(df), cond)); print(f"  {name:32} n={len(df):4d}  | {cond}")
def nt_of(code):                       # dominant bond length from the 11-char code
    for ch in ('A','B','6','7','9','5'):
        if ch in str(code): return 7 if ch in 'AB' else int(ch)
    return None

print("Figure_Source_Data (exact plotted points):")
# ---- Figure 2 ----
d=pd.read_csv(f"{F2}/Fig2a_NNdG_ExpdG_saltfixed/NNdG_ExpdG_saltfixed.csv")
save("Fig2a_NNdG_vs_ExpdG", d, "all robust interfaces: NN-predicted vs measured dG (no subset)")

d=pd.read_csv(f"{F2}/Fig2b_Tm_vs_BindingSites_x_Length/data.csv"); d["NucleotideBonds"]=d["NucleotideBonds"].astype(str)
m=((d.NucleotideBonds=="7")&(d.Num_Bonds.isin([4,5])))|((d.NucleotideBonds=="5")&(d.Num_Bonds.isin([9,10])))
save("Fig2b_Tm_by_Sites", d[m], "SUBSET: 7-nt bonds at 4/5 sites + 5-nt bonds at 9/10 sites")

d=pd.read_csv(f"{F2}/Fig2c_Alpha_PerBond_Efficiency/alpha_per_interface.csv"); d["nt"]=d["code"].map(nt_of)
m=((d.nt==7)&(d.n_sites.isin([4,5])))|((d.nt==5)&(d.n_sites.isin([9,10])))
save("Fig2c_Alpha_by_Sites", d[m][["code","nt","n_sites","N_bonds","n_reps","alpha_median","alpha_mean"]],
     "SUBSET: cooperativity alpha for the Fig 2b interfaces (7-nt 4/5 sites + 5-nt 9/10 sites)")

d=pd.read_csv(f"{F2}/Fig2d_ddG/test2_add_one_bond_pairs.csv")
save("Fig2d_ddG_addbond", d[d.regime.isin(["4->5","9->10"])],
     "SUBSET: one-added-site ddG, regimes 4->5 (7-nt) and 9->10 (5-nt) only")

d=pd.read_csv(f"{F2}/Fig2e_Avidity/avidity_per_interface.csv")
d["bond_length"]=d["nt"].astype(str).map({"5":"5nt","7":"7nt","A":"7nt","B":"7nt"})   # 7nt incl A/B anchors; 6nt dropped
d=d[d["bond_length"].notna()]
save("Fig2e_Avidity", d[["code","nt","bond_length","N_bonds","n_reps","avidity_enh_perbond","EM_per_bond_M"]],
     "SUBSET: avidity enhancement grouped 7-nt (7/A/B) vs 5-nt bonds (6-nt excluded)")

# ---- Figure 3 (panels plot their full data.csv) ----
save("Fig3a_Tm_vs_MeanGap", pd.read_csv(f"{F3}/Fig3a_Tm_vs_MeanGap_7nt_4bonds/data.csv"),
     "SUBSET: 7-nt bonds, 4 binding sites (Tm vs mean bond gap)")
save("Fig3b_Tm_vs_SpanResidual", pd.read_csv(f"{F3}/Fig3b_Tm_vs_BondSpan_OLSResiduals/data.csv"),
     "all interfaces: Tm residual vs bond-span residual (other covariates controlled)")
save("Fig3c_Tm_vs_BondPosition", pd.read_csv(f"{F3}/Fig3c_Tm_vs_BondPosition_ConstantSeq/data.csv"),
     "SUBSET: constant-sequence positional series (6 constructs)")
save("Fig3d_Occupancy_vs_Distance", pd.read_csv(f"{F3}/Fig3d_RelativeOccupancy_vs_Distance/data.csv"),
     "SUBSET: 7-nt interfaces, inferred occupancy vs nearest-bond distance")

# ---- Figure 4b + Figure 6 (model panels, plot full data.csv) ----
save("Fig4b_Pred_vs_Exp", pd.read_csv(f"{F4}/Phenomenological_Model/pred_vs_exp.csv"),
     "all interfaces: phenomenological out-of-fold predicted vs experimental Tm")
save("Fig6b_XGB_Pred_vs_Exp", pd.read_csv(f"{F6}/Fig6b_XGBoost_StabilityScatter/data.csv"),
     "all interfaces: XGBoost out-of-fold predicted vs experimental Tm")
save("Fig6c_XGB_RankRank", pd.read_csv(f"{F6}/Fig6c_XGBoost_RankRank/data.csv"),
     "all interfaces: XGBoost predicted vs experimental Tm rank")
save("Fig6d_ModelComparison", pd.read_csv(f"{F4}/Phenomenological_Model/combined_model_performance.csv"),
     "per-model Spearman rho, 6 models incl. Phenomenological (the panel's full bar set)")

# ---- Figure 4c-f (model-generated theory curves, not experimental points) ----
save("Fig4c_MeltingCurves_vs_NumBonds", pd.read_csv(f"{F4}/theory_vs_Nbonds.csv"),
     "THEORY: predicted melt (Tm, theta@40C) per bond count 2-12")
save("Fig4d_Decomp_vs_NumBonds", pd.read_csv(f"{F4}/EntropyPenalty_vs_NumBonds_spacing4.csv"),
     "THEORY: free-energy decomposition vs number of bonds")
save("Fig4e_MeltingCurves_vs_Spacing", pd.read_csv(f"{F4}/theory_vs_spacing.csv"),
     "THEORY: predicted melt (Tm, theta@40C) vs bond spacing")
save("Fig4f_Decomp_vs_Spacing", pd.read_csv(f"{F4}/EntropyPenalty_vs_Spacing_4bonds.csv"),
     "THEORY: free-energy decomposition vs bond spacing")

# ---- Figure 5 (coarse-grained MD) — canonical self-contained source ----
_AH="Figure5_MD_Simulations/Figure_AH_SelfContained"
MD=f"{_AH}/plots_robust_v2"                 # thermodynamic panels b/c/d/e/f (v2 rerun)
MECH=f"{_AH}/plots_robust/mechanism"        # D swarm, F registration, G/H heatmaps
_tm=["panel","series","x","Tm","bonds"]
save("Fig5b_Tm_vs_BondSpacing", pd.read_csv(f"{MD}/b/panel_b.csv")[_tm],
     "MD: simulated Tm (kT) vs bond spacing; flexible vs rigid scaffold")
d=pd.read_csv(f"{MD}/d/panel_d.csv"); save("Fig5c_Tm_vs_LinkerRigidity", d[d.series.isin(["sp2","sp6"])][_tm],
     "MD: simulated Tm (kT) vs linker rigidity; bond spacings 2 & 6 (rigid scaffold rod80)")
c=pd.read_csv(f"{MD}/c/panel_c.csv"); save("Fig5e_Tm_vs_ScaffoldRigidity", c[c.series.isin(["sp2","sp6"])][_tm],
     "MD: simulated Tm (kT) vs scaffold rigidity (rod 4/20/80/240); bond spacings 2 & 6")
save("Fig5d_BinderDistance_byLinker", pd.read_csv(f"{MECH}/panelC_histograms_raw.csv"),
     "MD: binder (tip-tip) distance per linker bending stiffness 4/10/20/40 (kT=2 swarm)")
rs=pd.read_csv(f"{MECH}/registration_summary.csv"); rs=rs[rs.rod_bend.isin([4,20,80]) & (rs.temperature==1.5)]
save("Fig5f_ExtensionSpread_vs_RodBending", rs[[c for c in ["rod_bend","temperature","intra_mean","intra_sigma","intra_sigma_reperr"] if c in rs.columns]],
     "MD: extension-position spread (intra_sigma) vs rod bending stiffness 4/20/80 (sp4, kT=1.5)")
hg=pd.read_csv(f"{MECH}/heatmap_grid.csv"); hg2=hg[hg.temperature==2.0]
save("Fig5g_MeanTipGap_heatmap", hg2[["chain_bend","rod_bend","tipgap_mean"]],
     "MD: mean tip gap over linker (chain_bend) x scaffold (rod_bend), kT=2")
save("Fig5h_TipGapSpread_heatmap", hg2[["chain_bend","rod_bend","tipgap_sigma"]],
     "MD: tip-gap spread (sigma) over linker x scaffold, kT=2")

with open(f"{OUT}/README.txt","w") as fh:
    fh.write("Exact plotted points per manuscript panel (regenerated by make_figure_points.py, last step of run_all).\n"
             "SUBSET = panel plots a conditional subset of the full dataframe; the filter is applied here.\n\n")
    for n,k,c in readme: fh.write(f"{n+'.csv':36} n={k:<5} {c}\n")
print(f"wrote {len(readme)} panel point-files + README to {OUT}/")
