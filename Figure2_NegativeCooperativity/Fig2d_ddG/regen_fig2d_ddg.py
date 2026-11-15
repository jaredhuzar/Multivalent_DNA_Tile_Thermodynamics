import os as _os
def _find_repo():
    d=_os.path.dirname(_os.path.abspath(__file__))
    while d!='/' and not _os.path.exists(_os.path.join(d,'Filtered_DF_Tm.csv')):
        d=_os.path.dirname(d)
    return d
__REPO__=_find_repo()
"""
Fig 2d/2E backing data AND figure: 'add one bond' pairs.

For every measured construct 'hi' that has exactly one more bond than a measured construct 'lo'
(hi == lo with one empty position filled), compare the measured dDeltaG to the sequence-predicted
dDeltaG of the added bond:

  ddg_meas  = Exp_dG(hi) - Exp_dG(lo)                      [experimental; salt/C_T-independent]
  ddg_pred  = Predicted_deltaG(hi) - Predicted_deltaG(lo)  [salt-corrected NN dG of the added bond]
  efficiency = ddg_meas / ddg_pred

FIXES (round-2 review item 5):
  - Pairs are ENUMERATED from the data (every valid lo->hi), not read from a frozen 125-pair list
    that covered only two plates.  The 4->5 (added 7nt) set now uses every measured pair.
  - The measured dG uses the SAME robust set as Fig 2a/2c/2e: Filtered_DF_HS rows with
    Reliable==True & Confidence=='high', grouped by X->5 code, median Exp_DeltaG.
  - This script now GENERATES test2_measured_only.svg (both subplots) so the panel is reproducible
    (it was previously a static orphan with stale counts).

Run from anywhere.  Writes test2_add_one_bond_pairs.csv, test2_measured_only_summary.csv,
test2_measured_only.svg.
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
PROJ=__REPO__; os.chdir(PROJ)
OUT="Figure2_NegativeCooperativity/Fig2d_ddG"
LEN={'5':5,'X':5,'6':6,'7':7,'9':9,'A':7,'B':7}   # added-bond nucleotide length by position label
def _nb(c): return sum(ch!='0' for ch in c)

# ---- robust per-construct measured & predicted dG (identical filter to Fig 2a/2c/2e) ----
hs=pd.read_csv("Filtered_DF_HS.csv"); hs["code"]=hs["Code"].astype(str).str.partition(".")[0]
sub=hs[(hs["Reliable"]==True)&(hs["Confidence"]=="high")]
g=(sub.groupby("code").agg(meas=("Exp_DeltaG","median"), pred=("Predicted_deltaG","median"),
                           n_reps=("Exp_DeltaG","size")).reset_index())
meas=g.set_index("code")["meas"]; pred=g.set_index("code")["pred"]; codes=set(g["code"])

# ---- enumerate EVERY add-one-bond pair (lo has one fewer bond than hi, identical elsewhere) ----
rows=[]
for hi in codes:
    for p,ch in enumerate(hi):
        if ch=='0': continue
        lo=hi[:p]+'0'+hi[p+1:]
        if lo in codes and _nb(lo)==_nb(hi)-1:
            rows.append(dict(lo=lo, hi=hi, regime=f"{_nb(lo)}->{_nb(hi)}", added_nt=LEN.get(ch,np.nan),
                             ddg_meas=float(meas[hi]-meas[lo]), ddg_pred=float(pred[hi]-pred[lo])))
d=pd.DataFrame(rows).sort_values(["regime","added_nt","lo","hi"]).reset_index(drop=True)
d["efficiency"]=d["ddg_meas"]/d["ddg_pred"]
d.to_csv(os.path.join(OUT,"test2_add_one_bond_pairs.csv"),index=False)

# ---- summary per (regime, added_nt) ----
from scipy.stats import wilcoxon
srows=[]
for (reg,ant),s in d.groupby(["regime","added_nt"]):
    x=s["ddg_meas"].dropna().values
    try: wp=wilcoxon(x, alternative="less").pvalue if len(x)>=6 else np.nan
    except Exception: wp=np.nan
    srows.append(dict(regime=reg, added_nt=ant, n_pairs=len(s), median_ddG=round(float(np.median(x)),3),
                      mean_ddG=round(float(np.mean(x)),3), pct_stabilizing=round(100*float(np.mean(x<0)),1),
                      wilcoxon_p_stabilizing=wp))
pd.DataFrame(srows).to_csv(os.path.join(OUT,"test2_measured_only_summary.csv"),index=False)

# ---- figure: left = measured ddG for BOTH regimes (4->5 stabilizes vs 9->10 no gain); right = 4->5 counts ----
g45 = d[(d.regime=="4->5") & (d.added_nt==7)]["ddg_meas"].values
g910= d[(d.regime=="9->10") & (d.added_nt==5)]["ddg_meas"].values
x=g45
fig,ax=plt.subplots(1,2,figsize=(11,4.4),gridspec_kw={"width_ratios":[1.4,1]})
rng=np.random.default_rng(0); ax[0].axhline(0,color="#888",lw=1,ls="--")
groups=[(0,g45),(1,g910)]
for xi,vals in groups:
    ax[0].scatter(xi+rng.uniform(-0.2,0.2,len(vals)), vals, s=24,
                  c=np.where(vals<0,"#4C78A8","#E45756"), alpha=0.75, edgecolor="none")
    ax[0].plot([xi-0.3,xi+0.3],[np.median(vals)]*2, color="k", lw=2.5)
_all=np.concatenate([g45,g910])
for xi,vals in groups:
    ax[0].text(xi, _all.max()+0.3, f"n={len(vals)}\nmed={np.median(vals):+.2f}", ha="center", va="bottom", fontsize=8)
ax[0].set_ylim(_all.min()-0.3, _all.max()+1.1)
ax[0].set_xticks([0,1]); ax[0].set_xticklabels(["4→5\n(add 7nt)","9→10\n(add 5nt)"])
ax[0].set_ylabel("measured ΔΔG of the added bond (kcal/mol)")
ax[0].set_title("adding a bond: 4→5 stabilizes (ΔΔG<0); 9→10 does not")
# right: 4->5 binned counts (reproducible, labelled)
bins=[-np.inf,-1.0,0.0,1.0,np.inf]; labels=["< −1","[−1,0)","[0,1)","≥ 1"]
_cats=pd.cut(x,bins=bins); cnt=pd.Series(_cats).value_counts().reindex(_cats.categories).values
ax[1].bar(range(4),cnt,color=["#2E5E8C","#4C78A8","#E99","#E45756"])
ax[1].set_xticks(range(4)); ax[1].set_xticklabels(labels)
ax[1].set_xlabel("4→5 measured ΔΔG bin (kcal/mol)"); ax[1].set_ylabel("number of pairs")
for i,v in enumerate(cnt): ax[1].text(i,v+0.5,str(int(v)),ha="center",fontsize=9)
ax[1].set_title("4→5 counts by ΔΔG bin")
fig.tight_layout(); fig.savefig(os.path.join(OUT,"test2_measured_only.svg")); plt.close(fig)

print(f"Fig2d: {len(d)} pairs; 4->5 added7 n={len(g45)} med={np.median(g45):+.3f}; "
      f"9->10 added5 n={len(g910)} med={np.median(g910):+.3f}; counts(4->5 by bin)={list(map(int,cnt))}")
for _,r in pd.DataFrame(srows).iterrows():
    print(f"  {r.regime} added{int(r.added_nt)}: n={r.n_pairs} median={r.median_ddG:+.3f} %stab={r.pct_stabilizing}")
