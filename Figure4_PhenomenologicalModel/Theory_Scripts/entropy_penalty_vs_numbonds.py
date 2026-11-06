import os as _os
def _find_repo():
    d=_os.path.dirname(_os.path.abspath(__file__))
    while d!='/' and not _os.path.exists(_os.path.join(d,'Filtered_DF_Tm.csv')):
        d=_os.path.dirname(d)
    return d
__REPO__=_find_repo()
"""
Generate EntropyPenalty_vs_NumBonds_spacing4.csv (+ .svg): dG and its decomposition vs number of
binding sites (fixed 4 nm spacing, 7nt), from the occupancy model (best Nscale params).
entropic_penalty is PER BOND (N_slots = 2 * sites), consistent with dG_net's denominator.
Num_Bonds column = number of binding sites; each site = 2 doubled slots (duplex bonds).
Run from Tile_TileResults.
"""
import numpy as np, pandas as pd, nbformat, warnings
warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import os; os.chdir(__REPO__)

nb=nbformat.read("TmPredictor_NewTm.ipynb",as_version=4); G={"np":np,"pd":pd}
exec(nb.cells[1].source,G); exec(nb.cells[2].source,G)
R_CAL,C_TILE,EDGE,RLN,TH = G["R_CAL"],G["C_TILE"],G["EDGE"],G["RLN"],G["TH"]
_ci=pd.read_csv("Supplementary/TableS1_Model_Comparison_All/phenomenological_parameter_CIs.csv").set_index("parameter")["point_estimate"]
LAM,DS_INT,BASE,T = float(_ci["lambda_nm"]),float(_ci["dS_Bond_cal_mol_K"]),float(_ci["dS_0_cal_mol_K"]),310.15   # current best Nscale fit
SPACING=4.0
v=list(TH['7'].values()); DH7=np.mean([x[0] for x in v]); DS7=np.mean([x[1] for x in v])

def meand(z):
    z=np.asarray(z,float)
    return np.array([np.mean(np.abs(np.delete(z,i)-z[i])) for i in range(len(z))]) if len(z)>1 else np.zeros(len(z))
def slots(n_sites,sp=SPACING,c=EDGE/2):
    total=(n_sites-1)*sp; xs=[c-total/2+k*sp for k in range(n_sites)]; z=[]
    for x in xs: z+=[x,x]          # 2 doubled slots (duplex bonds) per site
    return np.array(z)

rows=[]
for n in range(1,23):
    z=slots(n); g=np.exp(-meand(z)/LAM)
    dH=float(np.sum(g*DH7)); dS_eff=float(np.sum(g*DS7)+len(z)*DS_INT+BASE+RLN)
    TdS=T*dS_eff/1000.0; dG_net=dH-TdS
    ent_pen=-len(z)*DS_INT*T/1000.0            # PER-BOND: N_slots = 2*n (matches dG_net)
    rows.append(dict(Num_Bonds=n, dH_kcal_mol=dH, TdS_37C_kcal_mol=TdS,
                     dG_no_entropic_penalty_kcal_mol=dG_net-ent_pen,
                     entropic_penalty_kcal_mol=ent_pen, dG_net_kcal_mol=dG_net))
d=pd.DataFrame(rows); d.to_csv("Figure4_PhenomenologicalModel/EntropyPenalty_vs_NumBonds_spacing4.csv",index=False)

fig,ax=plt.subplots(figsize=(8,6))
ax.plot(d.Num_Bonds,d.dG_net_kcal_mol,'o-',color="#4C78A8",lw=2.5,label="ΔG net (with entropic penalty)")
ax.plot(d.Num_Bonds,d.dG_no_entropic_penalty_kcal_mol,'s--',color="#E45756",lw=2.0,label="ΔG, no entropic penalty")
ax.set(xlabel="Number of binding sites",ylabel="ΔG (kcal/mol)",
       title="Entropic penalty vs valency (per-bond penalty; 4 nm spacing)")
ax.legend(frameon=False); plt.tight_layout(); plt.savefig("Figure4_PhenomenologicalModel/EntropyPenalty_vs_NumBonds_spacing4.svg",dpi=600)
print(f"wrote EntropyPenalty_vs_NumBonds_spacing4.csv (+svg), n={len(d)}; per-bond entropic_penalty (N=4 sites)={-8*DS_INT*T/1000:.3f}")
