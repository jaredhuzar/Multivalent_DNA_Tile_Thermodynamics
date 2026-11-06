import os as _os
def _find_repo():
    d=_os.path.dirname(_os.path.abspath(__file__))
    while d!='/' and not _os.path.exists(_os.path.join(d,'Filtered_DF_Tm.csv')):
        d=_os.path.dirname(d)
    return d
__REPO__=_find_repo()
import os, math, numpy as np, pandas as pd, nbformat, warnings
warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import shap
from xgboost import XGBRegressor
PROJ=__REPO__; os.chdir(PROJ)
OUT="Supplementary/S6_SHAP_Values"

# ---- rebuild NN tables + make_Energy_Array from the notebook ----
nb=nbformat.read("Tile_Tile_Master.ipynb",as_version=4)
G={"pd":pd,"np":np,"math":math}
for i in (3,4,6,7,8): exec(nb.cells[i].source,G)
for i in (20,21): exec(nb.cells[i].source,G)
for t in ("nt5","nt6","nt7","nt9"): G[t]=pd.read_csv(f"{t}.csv",index_col=0)
G["nt9_stack"]=pd.read_csv("nt9_stack.csv",index_col=0)
exec(nb.cells[30].source,G)        # builds nt5_2 etc.
exec(nb.cells[31].source,G)        # make_Energy_Array
make_Energy_Array=G["make_Energy_Array"]

# ---- 235-interface Tm table: X = per-position NN dG, y = Tm ----
tm=pd.read_csv("Filtered_DF_Tm.csv"); tm["code"]=tm["Samples"].astype(str).str.partition(".")[0]
agg=tm.groupby("code").agg(Tm=("Tm","median")).reset_index().rename(columns={"code":"Samples"})
agg=agg.reset_index(drop=True)
G_ns={**G,"df":agg}
X=make_Energy_Array(agg); y=agg["Tm"].values.astype(float)
print(f"X shape {X.shape}  (samples x helix positions);  target = Tm;  n={len(y)}")

# ---- train XGBoost on ALL n=235 (no fold) ----
model=XGBRegressor(n_estimators=500,learning_rate=0.05,max_depth=4,subsample=0.8,
                   colsample_bytree=0.8,tree_method="hist",random_state=0)
model.fit(X,y)
expl=shap.TreeExplainer(model); sv=expl(X)
SH=np.asarray(sv.values); FV=np.asarray(sv.data)
nfeat=SH.shape[1]; names=[f"Helix {i}" for i in range(nfeat)]

# ---- save data ----
df=pd.DataFrame(SH,columns=[f"SHAP_helix{i}" for i in range(nfeat)])
for j in range(nfeat): df[f"featval_helix{j}"]=FV[:,j]
df.to_csv(os.path.join(OUT,"data.csv"),index=False)

# ---- beeswarm, fixed order Helix 0 (top) -> Helix 10 (bottom) ----
cmap=getattr(getattr(shap,"plots",None),"colors",None)
cmap=getattr(cmap,"red_blue",None) or plt.cm.RdBu_r
fig,ax=plt.subplots(figsize=(8.2,6.4))
rng=np.random.default_rng(0)
for i in range(nfeat):
    x=SH[:,i].astype(float); v=FV[:,i].astype(float)
    lo,hi=np.nanpercentile(v,5),np.nanpercentile(v,95)
    col=np.clip((v-lo)/(hi-lo),0,1) if hi>lo else np.full_like(v,0.5)
    # density-aware vertical jitter
    N=len(x); nbins=60
    if x.max()>x.min(): edges=np.linspace(x.min(),x.max(),nbins+1)
    else: edges=np.array([x.min()-.5,x.max()+.5])
    wb=np.clip(np.digitize(x,edges)-1,0,len(edges)-2)
    maxc=max((wb==b).sum() for b in range(len(edges)-1))
    yoff=np.zeros(N)
    for b in range(len(edges)-1):
        m=np.where(wb==b)[0]; k=len(m)
        if k==0: continue
        s=np.linspace(-1,1,k) if k>1 else np.array([0.0])
        yoff[m]=s*0.42*min(1.0,k/max(1,0.5*maxc))
    sc=ax.scatter(x, i+yoff, c=col, cmap=cmap, s=14, alpha=0.85, linewidths=0, vmin=0, vmax=1)
ax.axvline(0,color="0.6",lw=0.8)
ax.set_yticks(range(nfeat)); ax.set_yticklabels(names)
ax.invert_yaxis()                                  # Helix 0 at TOP, Helix 10 at BOTTOM
ax.set_xlabel("SHAP value (impact on predicted Tm, °C)")
cb=fig.colorbar(plt.cm.ScalarMappable(cmap=cmap),ax=ax,ticks=[0,1],pad=0.02)
cb.ax.set_yticklabels(["low","high"]); cb.set_label("helix NN ΔG feature value")
plt.tight_layout(); plt.savefig(os.path.join(OUT,"figure.svg"),bbox_inches="tight"); plt.close()

# quick self-check: OOF-free fit quality
from scipy.stats import spearmanr
print("in-sample rho(pred,Tm)=%.3f"%spearmanr(model.predict(X),y)[0])
print("mean|SHAP| by helix:", {i:round(np.abs(SH[:,i]).mean(),3) for i in range(nfeat)})
print("saved figure.svg + data.csv in",OUT)
