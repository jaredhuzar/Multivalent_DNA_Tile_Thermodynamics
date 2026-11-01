import os as _os
def _find_repo():
    d=_os.path.dirname(_os.path.abspath(__file__))
    while d!='/' and not _os.path.exists(_os.path.join(d,'Filtered_DF_Tm.csv')):
        d=_os.path.dirname(d)
    return d
__REPO__=_find_repo()
import os, numpy as np, pandas as pd
from scipy.stats import spearmanr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
PROJ=__REPO__; os.chdir(PROJ)
OUT=os.path.join(PROJ,"Supplementary","CooperativeUnit_vs_BondCount")

hs=pd.read_csv("Filtered_DF_HS.csv"); hs["code"]=hs["Code"].astype(str).str.partition(".")[0]
hs=hs[~hs["code"].str.contains("9")]                      # 9nt excluded
cm=pd.read_csv("Experimental_Trends/construct_medians.csv"); cm["code"]=cm["code"].astype(str).str.partition(".")[0]
g=hs.groupby("code").agg(dH=("Enthalpy","median"),dS=("Entropy","median"),sites=("Num_Bonds","first")).reset_index()
g=g.merge(cm[["code","bond_len"]],on="code",how="left").dropna(subset=["dH","dS","sites"])
g["N_bonds"]=2*g["sites"].astype(int)                     # each binding site = 2 bonds
g.to_csv(os.path.join(OUT,"coop_unit_vs_bondcount.csv"),index=False)

rH=spearmanr(g.N_bonds,g.dH); rS=spearmanr(g.N_bonds,g.dS)
fig,ax=plt.subplots(1,2,figsize=(13,5.2))
for a,(col,lab,rr,ylab) in zip(ax,[("dH","ΔH",rH,"van't Hoff ΔH (kcal/mol)"),("dS","ΔS",rS,"van't Hoff ΔS (cal/mol/K)")]):
    sc=a.scatter(g.N_bonds+np.random.normal(0,0.15,len(g)),g[col],c=g.bond_len,cmap="viridis",s=26,alpha=0.7,edgecolor="none")
    med=g.groupby("N_bonds")[col].median()
    a.plot(med.index,med.values,"o-",color="#E4572E",lw=2.5,ms=8,label="median per bond count")
    a.axhline(g[col].median(),color="0.4",ls="--",lw=1.2,label=f"overall median = {g[col].median():.0f}")
    a.set(xlabel="number of bonds (2 per site)",ylabel=ylab,
          title=f"{lab} vs bond count — flat (Spearman ρ={rr[0]:+.2f}, p={rr[1]:.1e})")
    a.legend(fontsize=8.5,loc="upper right")
cb=fig.colorbar(sc,ax=ax,fraction=0.025,pad=0.02); cb.set_label("bond length (nt)")
ax[0].text(0.03,0.06,"≈ one helix (cooperative unit),\nnot the sum of all bonds",transform=ax[0].transAxes,fontsize=9,color="0.25")
fig.suptitle("Cooperative melting unit is ~1 helix: ΔH and ΔS barely change from 4 to 20 bonds (9nt excluded, n=%d)"%len(g),y=1.02,fontsize=12)
plt.savefig(os.path.join(OUT,"coop_unit_vs_bondcount.svg"),bbox_inches="tight")
print(f"n={len(g)} | ΔH: Spearman ρ={rH[0]:+.3f} p={rH[1]:.1e}, median={g.dH.median():.1f} kcal/mol")
print(f"        ΔS: Spearman ρ={rS[0]:+.3f} p={rS[1]:.1e}, median={g.dS.median():.1f} cal/mol/K")
print("median ΔH, ΔS by bond count:")
print(g.groupby("N_bonds").agg(n=("dH","size"),dH=("dH","median"),dS=("dS","median")).round(1).to_string())
print("saved coop_unit_vs_bondcount.{csv,svg} in",OUT)
