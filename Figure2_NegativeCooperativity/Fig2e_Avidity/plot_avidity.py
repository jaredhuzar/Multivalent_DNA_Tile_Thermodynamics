import os as _os
def _find_repo():
    d=_os.path.dirname(_os.path.abspath(__file__))
    while d!='/' and not _os.path.exists(_os.path.join(d,'Filtered_DF_Tm.csv')):
        d=_os.path.dirname(d)
    return d
__REPO__=_find_repo()
"""
Avidity (companion to the cooperativity alpha).

Two faces of the SAME interfaces (robust dG set, 9nt excluded):
  (A) Avidity GROWS with valency: total measured interface free energy |dG_Exp| vs number of bonds,
      against the flat single-bond monovalent |dG_mono|.  More bonds => stronger overall interface,
      even though each bond's per-bond efficiency (alpha) falls.
  (B) Avidity enhancement over the working concentration: per-bond effective molarity relative to the
      4.5 nM bulk strand concentration, EM / C_exp.  This is the Whitesides avidity axis (reference =
      bulk), distinct from cooperativity (reference = 1 M standard state).
        EM < 1 M  -> negative cooperativity (alpha<1)
        EM >> C   -> strong avidity          (this panel)

Inputs: Effective_Molarity/effective_molarity.csv (EM per interface),
        Alpha_PerBond_Efficiency/alpha_per_interface.csv (dG_mono, alpha).
Run from Tile_TileResults.  Writes avidity_per_interface.csv + avidity_figure.svg here.
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
PROJ=__REPO__; os.chdir(PROJ)
OUT=os.path.join(PROJ,"Figure2_NegativeCooperativity","Fig2e_Avidity")
RT=1.987204e-3*310.15; C_EXP=4.5e-9

em=pd.read_csv("Extras/Effective_Molarity/effective_molarity.csv")
al=pd.read_csv("Figure2_NegativeCooperativity/Fig2c_Alpha_PerBond_Efficiency/alpha_per_interface.csv")
d=em.merge(al[["code","dG_mono_median_kcal_mol","alpha_median"]],on="code",how="left")
# nt bond type per construct
hs=pd.read_csv("Filtered_DF_HS.csv"); hs["code"]=hs["Code"].astype(str).str.partition(".")[0]
nt=hs.groupby("code")["NucleotideBonds"].first().astype(str)
d["nt"]=d["code"].map(nt)

d["avidity_enh_perbond"]=d["EM_per_bond_M"]/C_EXP            # fold, local vs bulk
d["avidity_dG_perbond_kcal"]=RT*np.log(d["EM_per_bond_M"]/C_EXP)  # +kcal/mol favourable gain per extra bond vs bulk
d["dG_Exp_abs"]=d["Exp_dG"].abs()
d.to_csv(os.path.join(OUT,"avidity_per_interface.csv"),index=False)

print(f"interfaces: {len(d)}")
print(f"avidity enhancement EM/C_exp:  median={d.avidity_enh_perbond.median():.2e}  "
      f"(= {RT*np.log(d.avidity_enh_perbond.median()):.1f} kcal/mol per extra bond vs 4.5 nM bulk)")
print(f"total interface |dG_Exp|: {d.dG_Exp_abs.median():.1f} kcal/mol (median) vs single-bond |dG_mono| "
      f"{d.dG_mono_median_kcal_mol.abs().median():.1f}")

cmap={"7":"#4C78A8","5":"#E45756","6":"#72B7B2","A":"#54A24B","B":"#B279A2"}
fig,ax=plt.subplots(1,2,figsize=(13.5,5.3))

# Panel A: avidity (total interface strength) grows with valency
for t,sub in d.groupby("nt"):
    ax[0].scatter(sub.N_bonds,sub.dG_Exp_abs,s=26,alpha=.45,color=cmap.get(t,"0.5"),label=f"{t}nt")
gb=d.groupby("N_bonds")["dG_Exp_abs"].median()
ax[0].plot(gb.index,gb.values,"k-o",ms=6,lw=1.6,label="median")
mono=d.dG_mono_median_kcal_mol.abs().median()
ax[0].axhline(mono,color="0.4",ls="--",lw=1.4,label=f"single monovalent bond ({mono:.1f})")
ax[0].set(xlabel="number of bonds  (= 2 x sites)",ylabel="interface |$\\Delta$G$_{Exp}$| (kcal/mol)",
          title="Avidity: total interface strength grows with valency")
ax[0].legend(fontsize=8.3,loc="upper left")

# Panel B: avidity enhancement over the 4.5 nM bulk
le=np.log10(d.avidity_enh_perbond)
ax[1].hist(le,bins=22,color="#54A24B",alpha=.85,edgecolor="white")
ax[1].axvline(le.median(),color="C3",lw=2,label=f"median = {d.avidity_enh_perbond.median():.0e}×")
ax[1].axvline(0,color="k",ls="--",lw=1.2,label="EM = C$_{exp}$ (no enhancement)")
ax[1].set(xlabel="log$_{10}$( EM / C$_{exp}$ ),  C$_{exp}$=4.5 nM",ylabel="interfaces",
          title="Avidity enhancement: local conc. vs 4.5 nM bulk")
ax[1].legend(fontsize=8.6)
plt.tight_layout(); plt.savefig(os.path.join(OUT,"avidity_figure.svg"),bbox_inches="tight")
print("saved avidity_per_interface.csv + avidity_figure.svg in",OUT)
