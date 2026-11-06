import os as _os
def _find_repo():
    d=_os.path.dirname(_os.path.abspath(__file__))
    while d!='/' and not _os.path.exists(_os.path.join(d,'Filtered_DF_Tm.csv')):
        d=_os.path.dirname(d)
    return d
__REPO__=_find_repo()
"""
Geometry-controlled regression of the avidity coefficient alpha.

alpha (per-interface avidity, = dG_avg^poly / dG^mono, from Alpha_PerBond_Efficiency) is regressed
on the interface structural factors to isolate each one's INDEPENDENT effect, exactly as done for Tm:
  - OLS multiple regression (standardized predictors): coefficient, t, p, VIF, R^2
  - partial Spearman (rank, controlling the other three)
  - partial-regression (added-variable) plots: alpha-residual vs factor-residual, controlling the rest

Predictors: number of bonds (N), bond span (nm, clustering), mean bond length (nt), GC fraction.
Set: robust-dG alpha interfaces, 9nt excluded (n=213).

Run from Tile_TileResults.  Writes alpha_geometry_data.csv and alpha_geometry_controlled.svg here.
"""
import os, numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata, t as tdist
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
PROJ=__REPO__; os.chdir(PROJ)
OUT=os.path.join(PROJ,"Supplementary","Alpha_Geometry_Controlled")

# ---- assemble alpha + geometry ----
a=pd.read_csv("Figure2_NegativeCooperativity/Fig2c_Alpha_PerBond_Efficiency/alpha_per_interface.csv")
cm=pd.read_csv("Experimental_Trends/construct_medians.csv"); cm["code"]=cm["code"].astype(str).str.partition(".")[0]
d=a.merge(cm[["code","span_nm","bond_len","GC","MeanGap_nm"]],on="code",how="left")
d=d.rename(columns={"N_bonds":"Num_Bonds"})
preds=["Num_Bonds","span_nm","bond_len","GC"]; resp="alpha_median"
labels={"Num_Bonds":"number of bonds","span_nm":"bond span (nm)","bond_len":"bond length (nt)","GC":"GC fraction"}
d=d.dropna(subset=[resp]+preds).reset_index(drop=True)
print(f"n interfaces = {len(d)}")

# ---- predictor collinearity ----
print("\n=== predictor collinearity (Spearman) ===")
print(d[preds].corr(method="spearman").round(2).to_string())

# ---- OLS (standardized) with t-tests and VIF ----
def z(v): return (v-v.mean())/v.std()
X=np.column_stack([np.ones(len(d))]+[z(d[c]).values for c in preds]); y=d[resp].values
b,*_=np.linalg.lstsq(X,y,rcond=None); dof=len(d)-X.shape[1]
s2=np.sum((y-X@b)**2)/dof; cov=s2*np.linalg.inv(X.T@X); se=np.sqrt(np.diag(cov)); tval=b/se
pval=2*(1-tdist.cdf(np.abs(tval),dof)); R2=1-np.sum((y-X@b)**2)/np.sum((y-y.mean())**2)
print(f"\n=== OLS: alpha ~ (standardized) {' + '.join(preds)}   R2={R2:.3f}, n={len(d)} ===")
for nm,cf,s,tt,pp in zip(["intercept"]+preds,b,se,tval,pval):
    print(f"  {nm:12s}: beta={cf:+.4f}  SE={s:.4f}  t={tt:+.2f}  p={pp:.2e}")
print("  VIF:",end=" ")
for c in preds:
    o=[k for k in preds if k!=c]; Xo=np.column_stack([np.ones(len(d))]+[d[k].values for k in o])
    bo,*_=np.linalg.lstsq(Xo,d[c].values,rcond=None); r2=1-np.sum((d[c].values-Xo@bo)**2)/np.sum((d[c].values-d[c].mean())**2)
    print(f"{c}={1/(1-r2):.1f}",end="  ")
print()

# ---- partial Spearman (on controlled residuals) + plots ----
def pcorr_resid(y_, x_, ctrl):   # returns residual arrays + partial Spearman
    m=d[[y_,x_]+ctrl].dropna()
    def ores(t):
        Xc=np.column_stack([np.ones(len(m))]+[m[k].values for k in ctrl]); bb=np.linalg.lstsq(Xc,m[t].values,rcond=None)[0]
        return m[t].values-Xc@bb
    rx,ry=ores(x_),ores(y_); rr,pp=spearmanr(rx,ry); return rx,ry,rr,pp

print("\n=== marginal vs partial Spearman (alpha ~ factor) ===")
rows=[]; fig,ax=plt.subplots(1,4,figsize=(19,4.6))
for j,c in enumerate(preds):
    ctrl=[k for k in preds if k!=c]
    mr,_=spearmanr(d[resp],d[c]); rx,ry,pr,pp=pcorr_resid(resp,c,ctrl)
    flag="***" if pp<0.001 else ("**" if pp<0.01 else ("*" if pp<0.05 else "ns"))
    print(f"  {labels[c]:16s}: marginal rho={mr:+.2f} -> partial rho={pr:+.2f} (p={pp:.1e}) {flag}")
    rows.append(dict(response="alpha",predictor=c,marginal_rho=round(mr,3),partial_rho=round(pr,3),partial_p=f"{pp:.1e}",n=len(d)))
    A=ax[j]; A.scatter(rx,ry,alpha=0.45,s=16,color="#4C78A8")
    bb=np.polyfit(rx,ry,1); xs=np.linspace(rx.min(),rx.max(),40); A.plot(xs,np.polyval(bb,xs),"r-",lw=1.8)
    A.set(xlabel=f"{labels[c]} (residual)",ylabel="alpha (residual)",title=f"alpha vs {labels[c]}\npartial rho={pr:+.2f}, p={pp:.1e}")
    # store residual coords
    d[f"xres_{c}"]=np.nan; d.loc[d.index[:len(rx)],f"xres_{c}"]=rx; d[f"yres_alpha_vs_{c}"]=np.nan; d.loc[d.index[:len(ry)],f"yres_alpha_vs_{c}"]=ry
fig.suptitle("Avidity alpha vs geometry, each factor controlling the other three (9nt excluded, n=%d)"%len(d),y=1.02,fontsize=12)
plt.tight_layout(); plt.savefig(os.path.join(OUT,"alpha_geometry_controlled.svg"),bbox_inches="tight")

pd.DataFrame(rows).to_csv(os.path.join(OUT,"alpha_geometry_partial_correlations.csv"),index=False)
keep=["code","alpha_median","alpha_mean","Num_Bonds","span_nm","bond_len","GC","MeanGap_nm",
      "Exp_dG_37C_kcal_mol","dG_avg_poly_kcal_mol","dG_mono_median_kcal_mol"]+ \
     [f"xres_{c}" for c in preds]+[f"yres_alpha_vs_{c}" for c in preds]
d[keep].to_csv(os.path.join(OUT,"alpha_geometry_data.csv"),index=False)
print("\nsaved alpha_geometry_data.csv, alpha_geometry_partial_correlations.csv, alpha_geometry_controlled.svg in",OUT)
