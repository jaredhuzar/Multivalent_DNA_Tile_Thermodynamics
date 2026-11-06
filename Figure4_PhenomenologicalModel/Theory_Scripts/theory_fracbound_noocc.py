import os as _os
def _find_repo():
    d=_os.path.dirname(_os.path.abspath(__file__))
    while d!='/' and not _os.path.exists(_os.path.join(d,'Filtered_DF_Tm.csv')):
        d=_os.path.dirname(d)
    return d
__REPO__=_find_repo()
import os, numpy as np, pandas as pd, nbformat, warnings
warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, seaborn as sns
os.chdir(__REPO__)
nb=nbformat.read("TmPredictor_NewTm.ipynb",as_version=4)
G={"np":np,"pd":pd}
exec(nb.cells[1].source,G); exec(nb.cells[2].source,G)   # constants + TH (Excel)
R_CAL=G["R_CAL"]; C_TILE=G["C_TILE"]; EDGE=G["EDGE"]; TH=G["TH"]
_ci=pd.read_csv("Supplementary/TableS1_Model_Comparison_All/phenomenological_parameter_CIs.csv").set_index("parameter")["point_estimate"]
FIT=dict(lam=float(_ci["lambda_nm"]), dS_int=float(_ci["dS_Bond_cal_mol_K"]), base=float(_ci["dS_0_cal_mol_K"]))   # current best Nscale fit

def _rep_thermo(L):
    t={5:'5',7:'7',9:'9'}.get(int(L),'7'); v=list(TH[t].values())
    return np.mean([x[0] for x in v]), np.mean([x[1] for x in v])
def _meand(z):
    z=np.asarray(z,float)
    if len(z)<=1: return np.zeros(len(z))
    return np.array([np.mean(np.abs(np.delete(z,i)-z[i])) for i in range(len(z))])
def make_interface(n_bonds, duplex_len=7, spacing_nm=None, center_nm=None, edge=EDGE):
    if center_nm is None: center_nm=edge/2
    if spacing_nm is None:
        span=0.8*edge; xs=np.linspace(center_nm-span/2,center_nm+span/2,n_bonds) if n_bonds>1 else [center_nm]
    else:
        total=(n_bonds-1)*spacing_nm; xs=[center_nm-total/2+k*spacing_nm for k in range(n_bonds)]
    dH,dS=_rep_thermo(duplex_len); bonds=[]
    for x in xs: bonds+=[(x,dH,dS),(x,dH,dS)]
    return bonds
def fracbound(bonds, occ=True, C=None):
    C=C_TILE if C is None else C
    g=np.exp(-_meand([b[0] for b in bonds])/FIT['lam']) if occ else np.ones(len(bonds))
    dH=sum(gi*b[1] for gi,b in zip(g,bonds)); dS=sum(gi*b[2] for gi,b in zip(g,bonds))
    dS_eff=dS+len(bonds)*FIT['dS_int']+FIT['base']+R_CAL*np.log(C/2.0)
    T_C=np.linspace(25,55,500); T_K=T_C+273.15
    dG=dH*1000.0-T_K*dS_eff; K=np.exp(-dG/(R_CAL*T_K)); return T_C, K/(1+K)

sns.set_theme(style="white", context="talk")
sps=[1,3,6,12,18,36]; cmap=plt.cm.viridis; norm=plt.Normalize(min(sps),max(sps))
fig,ax=plt.subplots(figsize=(8,6))
for sp in sps:
    T,f=fracbound(make_interface(4,duplex_len=7,spacing_nm=sp),occ=True)
    sns.lineplot(x=T,y=f,ax=ax,color=cmap(norm(sp)),linewidth=2.5,label=f"{sp}nm")
# --- the one extra line: dG with NO occupancy weighting (g=1); spacing-independent ---
T,f=fracbound(make_interface(4,duplex_len=7,spacing_nm=6),occ=False)
ax.plot(T,f,color="black",lw=3,ls="--",label="no occupancy\nweighting (g=1)")
ax.set_xlabel("Temperature (°C)"); ax.set_ylabel("Fraction Bound"); ax.set_ylim(-0.05,1.05)
ax.legend(title="Bond Spacing", bbox_to_anchor=(1.02,1), loc="upper left")
plt.tight_layout()
plt.savefig("Figure4_PhenomenologicalModel/TheoreticalBondSpacing_FracBound_noOccWeight.svg", dpi=600)
print("saved TheoreticalBondSpacing_FracBound_noOccWeight.svg (original untouched)")
print("g=1 baseline Tm (midpoint) ~ %.1f C"%(T[np.argmin(np.abs(f-0.5))]))
