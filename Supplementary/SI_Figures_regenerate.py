import os as _os
def _find_repo():
    d=_os.path.dirname(_os.path.abspath(__file__))
    while d!='/' and not _os.path.exists(_os.path.join(d,'Filtered_DF_Tm.csv')):
        d=_os.path.dirname(d)
    return d
__REPO__=_find_repo()
"""
Regenerate SI figures S2-S5, S7 from the CURRENT pipeline data (Filtered_DF_std.csv, which still
retains 9nt). One subfolder per figure under Supplementary/, each with figure.svg + data.csv.
Decisions: current pipeline; 9nt shown in S5 only (not S3/S4); S6/S8 handled separately.
"""
import os, ast, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import seaborn as sns
PROJ=__REPO__; os.chdir(PROJ)
BASE="Supplementary"
def out(sub):
    d=os.path.join(BASE,sub); os.makedirs(d,exist_ok=True); return d

d=pd.read_csv("Filtered_DF_std.csv")
d["Code"]=d["Code"].astype(str); d["nt"]=d["NucleotideBonds"].astype(str)
for c in ["Tm","Steepness","Entropy","Enthalpy","Predicted_tm","Num_Bonds"]:
    d[c]=pd.to_numeric(d[c],errors="coerce")
d["Entropy_kcal"]=d["Entropy"]/1000.0   # cal/mol/K -> kcal/mol/K (SI axis scale)

def per_construct(df):
    g=df.groupby("Code")
    m=g.agg(Tm=("Tm","median"),Steepness=("Steepness","median"),Entropy=("Entropy","median"),
            Entropy_kcal=("Entropy_kcal","median"),Enthalpy=("Enthalpy","median"),
            Predicted_tm=("Predicted_tm","median"),nt=("nt","first"),
            Num_Bonds=("Num_Bonds","first"),Positions=("Positions","first"),
            Experiment=("Experiment","first")).reset_index()
    return m

rel = d[d["Reliable"]==True]
phys = d[(d["Entropy"]<0)&(d["Enthalpy"]<0)]            # physical dH/dS = HS dataset
mTm   = per_construct(rel)
mPhys = per_construct(phys)

# ===================== S2: Steepness vs Num_Bonds, 5nt(9,10) / 7nt(4,5) =====================
sub=mTm[mTm["nt"].isin(["5","7"])].copy()
fig,ax=plt.subplots(1,2,figsize=(9,4.6),sharey=True)
for a,(nt,nbs,pal,lab) in zip(ax,[("5",[9,10],["#2ca25f","#d8721c"],"5nts"),
                                  ("7",[4,5],["#1f77b4","#d4a017"],"7nts")]):
    s=sub[(sub.nt==nt)&(sub.Num_Bonds.isin(nbs))]
    sns.barplot(data=s,x="Num_Bonds",y="Steepness",ax=a,palette=pal,errorbar=("ci",95))
    a.set_title(lab,fontsize=14); a.set_xlabel("Num_Bonds"); a.set_ylabel("Steepness" if nt=="5" else "")
plt.tight_layout(); plt.savefig(os.path.join(out("S2_BondLengthCount_Steepness"),"figure.svg"),bbox_inches="tight")
plt.close()
sub[["Code","nt","Num_Bonds","Steepness"]].rename(columns={"nt":"NucleotideBonds"}).to_csv(
    os.path.join(out("S2_BondLengthCount_Steepness"),"data.csv"),index=False)

# ===================== S3 / S4: Bond Position (Entropy / Tm), 9nt EXCLUDED =====================
def pos_table(m, col):
    m=m[~m["nt"].isin(["9"])].reset_index(drop=True)
    rows=[]
    for i,p in enumerate(m["Positions"].values):
        try: plist=ast.literal_eval(p) if isinstance(p,str) else list(p)
        except Exception: continue
        v=m.loc[i,col]
        if pd.isna(v): continue
        for pos in plist: rows.append({"Bond Position":int(pos), col:v})
    return pd.DataFrame(rows)

# S3 Entropy (physical set)
s3=pos_table(mPhys,"Entropy_kcal")
fig,axp=plt.subplots(figsize=(6.2,4.4))
sns.barplot(data=s3,x="Bond Position",y="Entropy_kcal",ax=axp,color="#2b7bba",errorbar=("ci",95))
axp.axhline(s3[s3["Bond Position"]==5]["Entropy_kcal"].mean(),color="red",ls="--",lw=2)
axp.set_ylabel("Entropy"); axp.set_xlabel("Bond Position")
plt.tight_layout(); plt.savefig(os.path.join(out("S3_BondPosition_Entropy"),"figure.svg"),bbox_inches="tight"); plt.close()
s3.to_csv(os.path.join(out("S3_BondPosition_Entropy"),"data.csv"),index=False)

# S4 Tm (reliable set)
s4=pos_table(mTm,"Tm")
fig,axp=plt.subplots(figsize=(6.2,4.4))
sns.barplot(data=s4,x="Bond Position",y="Tm",ax=axp,color="#2b7bba",errorbar=("ci",95))
axp.axhline(s4[s4["Bond Position"]==5]["Tm"].mean(),color="red",ls="--",lw=2)
axp.set_ylim(35,45); axp.set_ylabel("Melting Temperature"); axp.set_xlabel("Bond Position")
plt.tight_layout(); plt.savefig(os.path.join(out("S4_BondPosition_Melting"),"figure.svg"),bbox_inches="tight"); plt.close()
s4.rename(columns={"Tm":"Melting Temperature"}).to_csv(os.path.join(out("S4_BondPosition_Melting"),"data.csv"),index=False)

# ===================== S5: Entropy vs bond length (5,7,9) =====================
s5=mPhys[mPhys["nt"].isin(["5","7","9"])].copy()
s5["len"]=s5["nt"].astype(int)
fig,axp=plt.subplots(figsize=(6.2,4.4))
sns.barplot(data=s5,x="len",y="Entropy_kcal",ax=axp,color="#2b7bba",errorbar=("ci",95))
axp.set_xlabel("Length of Bonds (nt)"); axp.set_ylabel("Entropy")
plt.tight_layout(); plt.savefig(os.path.join(out("S5_BondLength_Entropy"),"figure.svg"),bbox_inches="tight"); plt.close()
s5[["Code","len","Entropy_kcal","Entropy"]].rename(columns={"len":"Length_nt"}).to_csv(
    os.path.join(out("S5_BondLength_Entropy"),"data.csv"),index=False)

# ===================== S7: NN-predicted Tm vs measured Tm =====================
s7=mTm[~mTm["nt"].isin(["9"])].dropna(subset=["Predicted_tm","Tm"]).copy()
from scipy.stats import spearmanr,pearsonr
rho=spearmanr(s7.Predicted_tm,s7.Tm)[0]; r=pearsonr(s7.Predicted_tm,s7.Tm)[0]
fig,axp=plt.subplots(figsize=(5.6,5.4))
axp.scatter(s7.Predicted_tm,s7.Tm,s=22,alpha=.5,color="#4C78A8")
lim=[min(s7.Predicted_tm.min(),s7.Tm.min())-1,max(s7.Predicted_tm.max(),s7.Tm.max())+1]
axp.plot(lim,lim,"k--",lw=1.2,label="y = x")
axp.set(xlabel="Predicted Tm (°C, nearest-neighbor)",ylabel="Measured Tm (°C)",
        title=f"Spearman ρ={rho:.2f}, Pearson r={r:.2f}, n={len(s7)}")
axp.legend(); plt.tight_layout()
plt.savefig(os.path.join(out("S7_PredictedVsMeasured_Tm"),"figure.svg"),bbox_inches="tight"); plt.close()
s7[["Code","nt","Num_Bonds","Predicted_tm","Tm"]].rename(columns={"nt":"NucleotideBonds"}).to_csv(
    os.path.join(out("S7_PredictedVsMeasured_Tm"),"data.csv"),index=False)

print("S2 constructs:",len(sub),"| S3 rows:",len(s3),"| S4 rows:",len(s4))
print("S5 by length:"); print(s5.groupby("len")["Entropy_kcal"].agg(["median","count"]))
print("S7 n:",len(s7),"rho=%.2f r=%.2f"%(rho,r))
print("done -> Supplementary/{S2..S5,S7}")
