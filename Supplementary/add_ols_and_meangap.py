import os as _os
def _find_repo():
    d=_os.path.dirname(_os.path.abspath(__file__))
    while d!='/' and not _os.path.exists(_os.path.join(d,'Filtered_DF_Tm.csv')):
        d=_os.path.dirname(d)
    return d
__REPO__=_find_repo()
"""
Add to Supplementary/:
  (A) The three geometry-controlled OLS regressions (response = Tm, Exp_dG, avidity alpha), each over
      the four factors {Num_Bonds, span_nm, bond_len, GC}.  9nt excluded.  Added-variable (partial-
      regression) panels + partial Spearman + standardized-OLS coefficient table (beta,t,p,VIF,R2).
  (B) van't Hoff Enthalpy and Entropy vs MEAN GAP for 7nt / 4-bond interfaces (high-confidence set).
"""
import os, ast, numpy as np, pandas as pd
from scipy.stats import spearmanr, t as tdist
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
PROJ=__REPO__; os.chdir(PROJ)
BASE="Supplementary"
def out(sub): d=os.path.join(BASE,sub); os.makedirs(d,exist_ok=True); return d
PREDS=["Num_Bonds","span_nm","bond_len","GC"]
LABELS={"Num_Bonds":"number of bonds","span_nm":"bond span (nm)","bond_len":"bond length (nt)","GC":"GC fraction"}

def run_ols(df, resp, folder, title, disp=None):
    disp=disp or resp
    d=df.dropna(subset=[resp]+PREDS).reset_index(drop=True); n=len(d)
    def z(v): return (v-v.mean())/v.std()
    X=np.column_stack([np.ones(n)]+[z(d[c]).values for c in PREDS]); y=d[resp].values
    b,*_=np.linalg.lstsq(X,y,rcond=None); dof=n-X.shape[1]
    s2=np.sum((y-X@b)**2)/dof; cov=s2*np.linalg.inv(X.T@X); se=np.sqrt(np.diag(cov)); tval=b/se
    pval=2*(1-tdist.cdf(np.abs(tval),dof)); R2=1-np.sum((y-X@b)**2)/np.sum((y-y.mean())**2)
    vif={}
    for c in PREDS:
        o=[k for k in PREDS if k!=c]; Xo=np.column_stack([np.ones(n)]+[d[k].values for k in o])
        bo,*_=np.linalg.lstsq(Xo,d[c].values,rcond=None); r2=1-np.sum((d[c].values-Xo@bo)**2)/np.sum((d[c].values-d[c].mean())**2)
        vif[c]=1/(1-r2)
    coef=pd.DataFrame([dict(predictor=nm,beta_std=round(cf,4),SE=round(s,4),t=round(tt,2),
                            p=f"{pp:.2e}",VIF=round(vif.get(nm,np.nan),2) if nm in vif else "")
                       for nm,cf,s,tt,pp in zip(["intercept"]+PREDS,b,se,tval,pval)])
    coef.attrs={}; 
    O=out(folder); coef.to_csv(os.path.join(O,"ols_coefficients.csv"),index=False)
    # partial residuals + added-variable plots
    def ores(m,t,ctrl):
        Xc=np.column_stack([np.ones(len(m))]+[m[k].values for k in ctrl]); bb=np.linalg.lstsq(Xc,m[t].values,rcond=None)[0]
        return m[t].values-Xc@bb
    rows=[]; fig,axs=plt.subplots(2,2,figsize=(11,9.4)); ax=axs.flatten(); dres=d.copy()
    for j,c in enumerate(PREDS):
        ctrl=[k for k in PREDS if k!=c]; m=d[[resp,c]+ctrl]
        rx=ores(m,c,ctrl); ry=ores(m,resp,ctrl)
        mr,_=spearmanr(d[resp],d[c]); pr,pp=spearmanr(rx,ry)
        rows.append(dict(response=resp,predictor=c,marginal_rho=round(mr,3),partial_rho=round(pr,3),
                         partial_p=f"{pp:.1e}",n=n))
        A=ax[j]; A.scatter(rx,ry,s=16,alpha=.45,color="#4C78A8")
        bb=np.polyfit(rx,ry,1); xs=np.linspace(rx.min(),rx.max(),40); A.plot(xs,np.polyval(bb,xs),"r-",lw=1.8)
        A.set(xlabel=f"{LABELS[c]} (residual)",ylabel=f"{disp} (residual)",
              title=f"{disp} vs {LABELS[c]}\npartial ρ={pr:+.2f}, p={pp:.1e}")
        dres[f"xres_{c}"]=rx; dres[f"yres_{resp}_vs_{c}"]=ry
    fig.suptitle(f"{title}  —  each factor controlling the other three (9nt excluded, n={n}, R²={R2:.2f})",y=1.0,fontsize=12)
    plt.tight_layout(); plt.savefig(os.path.join(O,"figure.svg"),bbox_inches="tight"); plt.close()
    pd.DataFrame(rows).to_csv(os.path.join(O,"partial_correlations.csv"),index=False)
    keep=["code"] if "code" in dres.columns else []
    keep=[k for k in (keep+[resp]+PREDS+[f"xres_{c}" for c in PREDS]+[f"yres_{resp}_vs_{c}" for c in PREDS]) if k in dres.columns]
    dres[keep].to_csv(os.path.join(O,"data.csv"),index=False)
    print(f"  {folder}: n={n}, R2={R2:.3f}")
    return R2,n

# ---- sources ----
cm=pd.read_csv("Experimental_Trends/construct_medians.csv"); cm["code"]=cm["code"].astype(str)
cm=cm[~cm["code"].str.contains("9")].reset_index(drop=True)   # exclude 9nt
cm_dg=cm.rename(columns={"dG37":"Exp_dG"})
ag=pd.read_csv("Supplementary/Alpha_Geometry_Controlled/alpha_geometry_data.csv")
if "code" not in ag.columns: ag["code"]=np.arange(len(ag))

print("OLS regressions:")
run_ols(cm, "Tm", "S3_OLS_Tm_geometry", "Melting Temperature vs geometry", disp="Melting Temperature (°C)")
run_ols(cm_dg, "Exp_dG", "S10_OLS_dG_geometry", "Experimental Free Energy vs geometry", disp="Experimental Free Energy (kcal/mol)")
run_ols(ag.rename(columns={"alpha_median":"alpha"}), "alpha", "S2_OLS_alpha_geometry", "Avidity α vs geometry", disp="avidity α")

# ===================== (B) dH / dS vs MEAN GAP, 7nt / 4-bond, high-confidence =====================
hs=pd.read_csv("Filtered_DF_HS.csv"); hs["code"]=hs["Code"].astype(str).str.partition(".")[0]
sub=hs[(hs.NucleotideBonds.astype(str)=="7")&(hs.Num_Bonds==4)&(hs.Reliable==True)&(hs.Confidence=="high")].copy()
def meangap_nm(posstr):
    p=sorted(ast.literal_eval(posstr) if isinstance(posstr,str) else list(posstr))
    if len(p)<2: return np.nan
    span=p[-1]-p[0]; return (span/(len(p)-1))*8.8     # mean gap = span/(k-1) sites, ×8.8 nm/site
g=(sub.groupby("code").agg(dH=("Enthalpy","median"),dS=("Entropy","median"),
                           Positions=("Positions","first"),n_reps=("Enthalpy","size")).reset_index())
g["MeanGap_nm"]=g["Positions"].apply(meangap_nm)
g=g.dropna(subset=["MeanGap_nm"]).reset_index(drop=True)
O=out("S4_HS_vs_MeanGap_7nt_4bonds")
g[["code","dH","dS","MeanGap_nm","n_reps"]].to_csv(os.path.join(O,"data.csv"),index=False)
fig,ax=plt.subplots(1,2,figsize=(12,4.8)); x=g.MeanGap_nm.values
for a,(col,lab,u) in zip(ax,[("dH","Enthalpy","kcal/mol"),("dS","Entropy","cal/mol/K")]):
    y=g[col].values; rho,p=spearmanr(x,y)
    a.scatter(x,y,s=24,alpha=.5,color="#4C78A8")
    bb=np.polyfit(x,y,1); xs=np.linspace(x.min(),x.max(),40); a.plot(xs,np.polyval(bb,xs),"r-",lw=1.6)
    a.set(xlabel="mean gap (nm)",ylabel=f"{lab} ({u})",
          title=f"{lab} vs mean gap (7nt, 4 bonds, high-conf, n={len(g)})\nSpearman ρ={rho:+.2f}, p={p:.2f}")
plt.tight_layout(); plt.savefig(os.path.join(O,"figure.svg"),bbox_inches="tight"); plt.close()
print(f"S4_HS_vs_MeanGap_7nt_4bonds: n={len(g)}, gap {x.min():.1f}-{x.max():.1f} nm")
print("  dH rho=%.3f  dS rho=%.3f"%(spearmanr(x,g.dH)[0],spearmanr(x,g.dS)[0]))
