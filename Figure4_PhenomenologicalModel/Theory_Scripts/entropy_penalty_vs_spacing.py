import os as _os
def _find_repo():
    d=_os.path.dirname(_os.path.abspath(__file__))
    while d!='/' and not _os.path.exists(_os.path.join(d,'Filtered_DF_Tm.csv')):
        d=_os.path.dirname(d)
    return d
__REPO__=_find_repo()
"""
Generate EntropyPenalty_vs_Spacing_4bonds.csv (+ .svg): dG and its decomposition vs bond spacing
for a 4-bond, 7nt interface, from the phenomenological occupancy model (best Nscale params).
Includes dG WITHOUT occupancy weighting (g=1), the spacing-independent additive baseline.

Model (per predict_Tm):  g_k = exp(-dbar_k/lam);  dH = sum g_k*dH7;
  dS_eff = sum g_k*dS7 + N_slots*dS_int + base + R*ln(C_T/2);  dG = dH - T*dS_eff/1000  (37 C).
Reported decomposition columns use the per-BOND penalty (N=4): entropic_penalty = -4*dS_int*T/1000.
Run from Tile_TileResults.
"""
import numpy as np, pandas as pd, nbformat, warnings
warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import os; os.chdir(__REPO__)

# physical constants + per-type NN thermo (from TmPredictor_NewTm cells 1-2)
nb=nbformat.read("TmPredictor_NewTm.ipynb",as_version=4); G={"np":np,"pd":pd}
exec(nb.cells[1].source,G); exec(nb.cells[2].source,G)
R_CAL,C_TILE,EDGE,RLN,TH = G["R_CAL"],G["C_TILE"],G["EDGE"],G["RLN"],G["TH"]
_ci=pd.read_csv("Supplementary/TableS1_Model_Comparison_All/phenomenological_parameter_CIs.csv").set_index("parameter")["point_estimate"]
LAM,DS_INT,BASE,T = float(_ci["lambda_nm"]),float(_ci["dS_Bond_cal_mol_K"]),float(_ci["dS_0_cal_mol_K"]),310.15   # current best Nscale fit
N_BONDS=4
v=list(TH['7'].values()); DH7=np.mean([x[0] for x in v]); DS7=np.mean([x[1] for x in v])

def meand(z):
    z=np.asarray(z,float)
    return np.array([np.mean(np.abs(np.delete(z,i)-z[i])) for i in range(len(z))]) if len(z)>1 else np.zeros(len(z))
def slots(sp,n=N_BONDS,c=EDGE/2):
    total=(n-1)*sp; xs=[c-total/2+k*sp for k in range(n)]; z=[]
    for x in xs: z+=[x,x]          # 2 doubled slots per bond
    return np.array(z)

ent_pen = -len(slots(1))*DS_INT*T/1000.0    # PER-BOND penalty: N_slots = 2*sites (matches dG_net)
dG_no_occ = (len(slots(1))*DH7) - T*(len(slots(1))*DS7 + len(slots(1))*DS_INT + BASE + RLN)/1000.0  # g=1

rows=[]
for sp in np.arange(1,51,1.0):
    z=slots(sp); g=np.exp(-meand(z)/LAM)
    dH=float(np.sum(g*DH7)); dS_eff=float(np.sum(g*DS7)+len(z)*DS_INT+BASE+RLN)
    TdS=T*dS_eff/1000.0; dG_net=dH-TdS
    rows.append(dict(Spacing_nm=sp, dH_kcal_mol=dH, TdS_37C_kcal_mol=TdS,
                     dG_no_entropic_penalty_kcal_mol=dG_net-ent_pen,
                     entropic_penalty_kcal_mol=ent_pen, dG_net_kcal_mol=dG_net,
                     dG_no_occ_weighting_kcal_mol=dG_no_occ))
d=pd.DataFrame(rows); d.to_csv("Figure4_PhenomenologicalModel/EntropyPenalty_vs_Spacing_4bonds.csv",index=False)

fig,ax=plt.subplots(figsize=(8,6))
ax.plot(d.Spacing_nm,d.dG_net_kcal_mol,'-',color="#4C78A8",lw=2.5,label="ΔG net (occupancy-weighted)")
ax.axhline(dG_no_occ,color="black",lw=2.5,ls="--",label="ΔG, no occupancy weighting (g=1)")
ax.set(xlabel="Bond Spacing (nm)",ylabel="ΔG (kcal/mol)",title="4-bond interface: ΔG vs spacing")
ax.legend(frameon=False); plt.tight_layout(); plt.savefig("Figure4_PhenomenologicalModel/EntropyPenalty_vs_Spacing_4bonds.svg",dpi=600)
print(f"wrote EntropyPenalty_vs_Spacing_4bonds.csv (+svg); dG_no_occ={dG_no_occ:.3f} kcal/mol, n={len(d)}")
