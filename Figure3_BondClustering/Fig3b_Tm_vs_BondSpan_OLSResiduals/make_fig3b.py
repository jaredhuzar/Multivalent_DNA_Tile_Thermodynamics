import os as _os
def _find_repo():
    d=_os.path.dirname(_os.path.abspath(__file__))
    while d!='/' and not _os.path.exists(_os.path.join(d,'Filtered_DF_Tm.csv')):
        d=_os.path.dirname(d)
    return d
__REPO__=_find_repo()
"""
Fig 3b: added-variable (partial-regression) plot of melting temperature vs bond span (g*).

Residual of Tm and residual of bond span after ordinary-least-squares regression of each on the other
three covariates (bond count, bond length, GC fraction). A nonzero slope here is the INDEPENDENT effect
of bond clustering (span) on Tm, controlling for the rest. Dashed red line is a linear fit.

This is the bond-span sub-panel of the four-covariate Tm OLS (SI Fig S3). The residuals are read from
`Supplementary/S3_OLS_Tm_geometry/data.csv` (columns xres_span_nm, yres_Tm_vs_span_nm), which
`add_ols_and_meangap.py` writes from the authoritative Filtered_DF_Tm set (n=235). Run after it.
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import spearmanr
PROJ=__REPO__; os.chdir(PROJ)
OUT="Figure3_BondClustering/Fig3b_Tm_vs_BondSpan_OLSResiduals"; os.makedirs(OUT,exist_ok=True)

src=pd.read_csv("Supplementary/S3_OLS_Tm_geometry/data.csv")
x=src["xres_span_nm"].values; y=src["yres_Tm_vs_span_nm"].values
rho,p=spearmanr(x,y); sl,ic=np.polyfit(x,y,1)
out=pd.DataFrame({"bond_span_residual_nm":x,"Tm_residual_C":y})
if "code" in src.columns: out.insert(0,"code",src["code"].values)
out.to_csv(os.path.join(OUT,"data.csv"),index=False)

fig,ax=plt.subplots(figsize=(5.4,5.0))
ax.scatter(x,y,s=18,alpha=0.45,color="#808080",edgecolor="none")
ax.axhline(0,color="0.75",lw=0.8,ls=":")
xs=np.linspace(x.min(),x.max(),50); ax.plot(xs,sl*xs+ic,"r--",lw=1.8)
ax.set(xlabel="Bond Span Residual (nm)",ylabel="Experimental T$_m$ Residual (°C)",
       title=f"Tm vs bond span (g*), other covariates controlled\npartial ρ={rho:+.2f}, p={p:.1e}, n={len(x)}")
plt.tight_layout(); plt.savefig(os.path.join(OUT,"figure.svg"),bbox_inches="tight"); plt.close()
print(f"Fig3b: n={len(x)}  partial rho={rho:+.3f}  slope={sl:+.3f}  -> {OUT}")
