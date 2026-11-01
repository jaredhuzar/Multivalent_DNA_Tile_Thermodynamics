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
OUT=os.path.join(PROJ,"Figure3_BondClustering","HS_vs_Distance_7nt_4bonds")

# high-confidence dH/dS dataset, 7nt bonds, 4 binding sites
d=pd.read_csv("Filtered_DF_HS.csv"); d["code"]=d["Code"].astype(str).str.partition(".")[0]
sub=d[(d["NucleotideBonds"]=="7") & (d["Num_Bonds"]==4) & (d["Reliable"]==True) & (d["Confidence"]=="high")]
g=(sub.groupby("code")
      .agg(dH=("Enthalpy","median"), dS=("Entropy","median"),
           Distances=("Distances","first"), n_reps=("Enthalpy","size")).reset_index())
g["bond_distance_nm"]=g["Distances"]*8.8          # mean pairwise inter-bond distance
g.to_csv(os.path.join(OUT,"hs_vs_distance_7nt_4bonds.csv"), index=False)

x=g["bond_distance_nm"].values
rH=spearmanr(x,g.dH); rS=spearmanr(x,g.dS)
fig,ax=plt.subplots(1,2,figsize=(13,5.2))
for a,(col,rr,ylab,c) in zip(ax,[("dH",rH,"van't Hoff ΔH (kcal/mol)","#4C78A8"),
                                  ("dS",rS,"van't Hoff ΔS (cal/mol/K)","#E4572E")]):
    a.scatter(x,g[col],alpha=0.5,s=28,color=c,edgecolor="none")
    b=np.polyfit(x,g[col],1); xs=np.linspace(x.min(),x.max(),30); a.plot(xs,np.polyval(b,xs),"k--",lw=1.6)
    a.set(xlabel="mean inter-bond distance (nm)",ylabel=ylab,
          title=f"{'ΔH' if col=='dH' else 'ΔS'} vs bond distance — flat (Spearman ρ={rr[0]:+.2f}, p={rr[1]:.2f})")
fig.suptitle("7nt bonds, 4 binding sites, high-confidence ΔH/ΔS set (n=%d constructs)"%len(g),y=1.02,fontsize=12)
plt.tight_layout(); plt.savefig(os.path.join(OUT,"hs_vs_distance_7nt_4bonds.svg"),bbox_inches="tight")
print(f"n={len(g)} constructs | distance {x.min():.1f}-{x.max():.1f} nm")
print(f"ΔH: Spearman ρ={rH[0]:+.3f} p={rH[1]:.2f}, median={g.dH.median():.1f} kcal/mol")
print(f"ΔS: Spearman ρ={rS[0]:+.3f} p={rS[1]:.2f}, median={g.dS.median():.1f} cal/mol/K")
print("saved hs_vs_distance_7nt_4bonds.{csv,svg}")
