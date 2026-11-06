import os as _os
def _find_repo():
    d=_os.path.dirname(_os.path.abspath(__file__))
    while d!='/' and not _os.path.exists(_os.path.join(d,'Filtered_DF_Tm.csv')):
        d=_os.path.dirname(d)
    return d
__REPO__=_find_repo()
"""
Report the two standard multivalency quantities per interface:
  (A) effective molarity  EM  (vs the 1 M bimolecular standard state; EM<1M => negative cooperativity)
  (B) avidity enhancement EM / C_exp  (vs the 4.5 nM working strand concentration; >>1 => strong avidity)

EM is obtained from  dG_Exp - dG_Seq(1 init) = (N-1) RT ln(1 M / EM)  (one bimolecular initiation for
the whole tile-tile docking; bonds 2..N intramolecular).  Robust experimental-dG set, 9nt excluded.
Writes avidity_raw_datapoints.csv and EM_and_avidity_enhancement.svg.  Run from Tile_TileResults.
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
PROJ=__REPO__; os.chdir(PROJ)
OUT=os.path.join(PROJ,"Figure2_NegativeCooperativity","Fig2e_Avidity")
RT=1.987204e-3*310.15; C_EXP=4.5e-9

em=pd.read_csv("Extras/Effective_Molarity/effective_molarity.csv")
al=pd.read_csv("Figure2_NegativeCooperativity/Fig2c_Alpha_PerBond_Efficiency/alpha_per_interface.csv")
hs=pd.read_csv("Filtered_DF_HS.csv"); hs["code"]=hs["Code"].astype(str).str.partition(".")[0]
nt=hs.groupby("code")["NucleotideBonds"].first().astype(str)

d=em.merge(al[["code","dG_mono_median_kcal_mol","alpha_median"]],on="code",how="left")
d["nt"]=d["code"].map(nt)
d["EM_mM"]=d["EM_per_bond_M"]*1e3
d["avidity_enhancement_fold"]=d["EM_per_bond_M"]/C_EXP
d["avidity_dG_perbond_kcal"]=RT*np.log(d["avidity_enhancement_fold"])

raw=d[["code","nt","n_sites","N_bonds","n_reps","Exp_dG","dG_mono_median_kcal_mol",
       "dG_Seq_1init","shortfall_kcal","EM_per_bond_M","EM_mM","alpha_median",
       "avidity_enhancement_fold","avidity_dG_perbond_kcal"]].rename(
       columns={"Exp_dG":"dG_Exp_kcal","dG_mono_median_kcal_mol":"dG_mono_median_kcal",
                "alpha_median":"cooperativity_alpha"})
raw.to_csv(os.path.join(OUT,"avidity_raw_datapoints.csv"),index=False)

emM=d["EM_per_bond_M"]; enh=d["avidity_enhancement_fold"]
print(f"n={len(d)}")
print(f"effective molarity : median={emM.median()*1e3:.2f} mM  IQR=[{emM.quantile(.25)*1e3:.2f},{emM.quantile(.75)*1e3:.2f}] mM  (100% < 1 M)")
print(f"avidity enhancement: median={enh.median():.2e}x  (= {RT*np.log(enh.median()):.1f} kcal/mol per extra bond vs 4.5 nM)")

fig,ax=plt.subplots(1,2,figsize=(13.2,5.1))
ax[0].hist(np.log10(emM),bins=24,color="#4C78A8",alpha=.85,edgecolor="white")
ax[0].axvline(0,color="k",ls="--",lw=1.3,label="1 M standard state")
ax[0].axvline(np.log10(emM.median()),color="C3",lw=2,label=f"median = {emM.median()*1e3:.2f} mM")
ax[0].set(xlabel="log$_{10}$(effective molarity / M)",ylabel="interfaces",
          title="(A) Effective molarity  (vs 1 M)\nEM < 1 M → negative cooperativity")
ax[0].legend(fontsize=9)
ax[1].hist(np.log10(enh),bins=24,color="#54A24B",alpha=.85,edgecolor="white")
ax[1].axvline(0,color="k",ls="--",lw=1.3,label="C$_{exp}$ = 4.5 nM (no gain)")
ax[1].axvline(np.log10(enh.median()),color="C3",lw=2,label=f"median = {enh.median():.0e}×")
ax[1].set(xlabel="log$_{10}$(avidity enhancement, EM / C$_{exp}$)",ylabel="interfaces",
          title="(B) Avidity enhancement  (vs 4.5 nM bulk)\nEM ≫ C$_{exp}$ → strong avidity")
ax[1].legend(fontsize=9)
plt.tight_layout(); plt.savefig(os.path.join(OUT,"EM_and_avidity_enhancement.svg"),bbox_inches="tight")
print("saved avidity_raw_datapoints.csv + EM_and_avidity_enhancement.svg in",OUT)
