import os as _os
def _find_repo():
    d=_os.path.dirname(_os.path.abspath(__file__))
    while d!='/' and not _os.path.exists(_os.path.join(d,'Filtered_DF_Tm.csv')):
        d=_os.path.dirname(d)
    return d
__REPO__=_find_repo()
import os, sys, numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import shap
from xgboost import XGBRegressor
PROJ=__REPO__; os.chdir(PROJ)
OUT="Supplementary/S6_SHAP_Values"
# ---- canonical 236-interface set and the SAME features as Figure 6 (updated 2026-11-13) ----
# X = per-site mean salt-corrected NN dG (37 C) of the site's two duplexes (0 = empty), y = median Tm;
# identical constructs/features to Figure6_MachineLearning/make_figure6.py (mv_model.site_dG_features).
sys.path.insert(0, os.path.join(PROJ, "Figure4_MultivalentModel")); import mv_model as mv
os.chdir(PROJ)
codes=list(mv.CODES); y=np.asarray(mv.Y_C,float)
X=np.vstack([mv.site_dG_features(c) for c in codes])
print(f"X shape {X.shape}  (samples x site positions);  target = Tm;  n={len(y)}")

# ---- train XGBoost on ALL interfaces (no fold) ----
model=XGBRegressor(n_estimators=500,learning_rate=0.05,max_depth=4,subsample=0.8,
                   colsample_bytree=0.8,tree_method="hist",random_state=0)
model.fit(X,y)
expl=shap.TreeExplainer(model); sv=expl(X)
SH=np.asarray(sv.values); FV=np.asarray(sv.data)
nfeat=SH.shape[1]; names=[f"Helix {i}" for i in range(nfeat)]

# ---- save data ----
df=pd.DataFrame(SH,columns=[f"SHAP_helix{i}" for i in range(nfeat)]); df.insert(0,"code",codes); df.insert(1,"Tm_C",y)
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
cb.ax.set_yticklabels(["low","high"]); cb.set_label("site NN ΔG feature value")
plt.tight_layout(); plt.savefig(os.path.join(OUT,"figure.svg"),bbox_inches="tight"); plt.close()

# quick self-check: OOF-free fit quality
from scipy.stats import spearmanr
print("in-sample rho(pred,Tm)=%.3f"%spearmanr(model.predict(X),y)[0])
print("mean|SHAP| by helix:", {i:round(np.abs(SH[:,i]).mean(),3) for i in range(nfeat)})
print("saved figure.svg + data.csv in",OUT)

# ---- mean SHAP at each site: occupied vs empty (SI companion panel) ----
rows=[]
for i in range(nfeat):
    occ=FV[:,i]!=0
    rows.append(dict(site=i, n_occupied=int(occ.sum()), mean_SHAP_occupied=float(SH[occ,i].mean()) if occ.any() else np.nan,
                     mean_SHAP_empty=float(SH[~occ,i].mean()) if (~occ).any() else np.nan, mean_abs_SHAP=float(np.abs(SH[:,i]).mean())))
MS=pd.DataFrame(rows); MS.round(4).to_csv(os.path.join(OUT,"mean_shap_by_site.csv"),index=False)
fig,ax=plt.subplots(figsize=(7.2,3.8)); xs=np.arange(nfeat); w=.38
ax.bar(xs-w/2,MS.mean_SHAP_occupied,width=w,color="#3B6AA0",label="site occupied")
ax.bar(xs+w/2,MS.mean_SHAP_empty,width=w,color="#c9c9c9",label="site empty")
ax.axhline(0,color="0.3",lw=.8)
for e in (0,10): ax.axvspan(e-.5,e+.5,color="#f6e3d6",zorder=0)
ax.set_xticks(xs); ax.set_xticklabels([f"{i}" for i in xs]); ax.set(xlabel="site (helix) position; 5 = reporter, 0/10 = edges",ylabel="mean SHAP value (°C)")
ax.legend(fontsize=8,frameon=False); ax.set_title("Mean SHAP contribution of each site (XGBoost, all 236 interfaces)",fontsize=9.5,loc="left")
plt.tight_layout(); plt.savefig(os.path.join(OUT,"mean_shap_by_site.svg"),bbox_inches="tight"); plt.close()
print("saved mean_shap_by_site.svg + .csv")
