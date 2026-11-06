import os as _os
def _find_repo():
    d=_os.path.dirname(_os.path.abspath(__file__))
    while d!='/' and not _os.path.exists(_os.path.join(d,'Filtered_DF_Tm.csv')):
        d=_os.path.dirname(d)
    return d
__REPO__=_find_repo()
"""
Fig 2a panel: NN-predicted dG37 (salt-corrected) vs experimental dG37, one point per interface.

Construct set / experimental dG: Filtered_DF_HS rows with Reliable & Confidence=='high',
aggregated per construct code by MEDIAN Exp_DeltaG (n_reps = #replicates).  These are experimental
and are NOT affected by the salt correction.

NN predicted dG = the salt-corrected additive sequence prediction Predicted_deltaG (sum of per-bond
SantaLucia NN dG37 + per-duplex salt correction, computed in Tile_Tile_Master).  After the salt fix
the per-bond and per-site salt conventions coincide (each L-bp duplex contributes (L-1) phosphate
steps; a binding site has two duplexes), so there is a single NN_Predicted column.

Writes NNdG_ExpdG_saltfixed.csv + .svg here.  Run from anywhere.
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
PROJ=__REPO__; os.chdir(PROJ)
OUT="Figure2_NegativeCooperativity/Fig2a_NNdG_ExpdG_saltfixed"

hs=pd.read_csv("Filtered_DF_HS.csv"); hs["code"]=hs["Code"].astype(str).str.partition(".")[0]
sub=hs[(hs["Reliable"]==True)&(hs["Confidence"]=="high")]
g=(sub.groupby("code")
      .agg(Exp_dG_37C=("Exp_DeltaG","median"),
           n_reps=("Exp_DeltaG","size"),
           NN_Predicted_dG_37C=("Predicted_deltaG","median"))
      .reset_index())
g=g[g["n_reps"]>=2].reset_index(drop=True)          # >=2 replicates (matches the original panel set)
g.to_csv(os.path.join(OUT,"NNdG_ExpdG_saltfixed.csv"),index=False)

x=g["NN_Predicted_dG_37C"].values; y=g["Exp_dG_37C"].values
r=pearsonr(x,y)[0]; rho=spearmanr(x,y)[0]
sl,ic=np.polyfit(x,y,1)
print(f"interfaces: {len(g)}  Pearson r={r:.3f}  Spearman rho={rho:.3f}  slope={sl:.3f} intercept={ic:.2f}")

fig,ax=plt.subplots(figsize=(6.2,6.0))
ax.scatter(x,y,s=26,alpha=0.55,color="#4C78A8",edgecolor="white",linewidth=0.4)
lo=min(x.min(),y.min())-2; hi=max(x.max(),y.max())+2
ax.plot([lo,hi],[lo,hi],"k--",lw=1.0,label="y = x")
xs=np.linspace(x.min(),x.max(),50)
ax.plot(xs,sl*xs+ic,"r-",lw=1.8,label=f"OLS: slope={sl:.2f}  (r={r:.2f})")
ax.set(xlabel="NN-predicted $\\Delta G_{37}$ (salt-corrected, kcal/mol)",
       ylabel="Experimental $\\Delta G_{37}$ (kcal/mol)",
       title="Predicted vs measured interface free energy")
ax.legend(fontsize=9,loc="upper left"); plt.tight_layout()
plt.savefig(os.path.join(OUT,"NNdG_ExpdG_saltfixed.svg"),bbox_inches="tight")
print("saved NNdG_ExpdG_saltfixed.csv + .svg in",OUT)
