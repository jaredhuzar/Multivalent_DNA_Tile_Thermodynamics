import os as _os
def _find_repo():
    d=_os.path.dirname(_os.path.abspath(__file__))
    while d!='/' and not _os.path.exists(_os.path.join(d,'Filtered_DF_Tm.csv')):
        d=_os.path.dirname(d)
    return d
__REPO__=_find_repo()
"""
Fig 2d backing data: 'add one bond' pairs.  For each measured pair (lo, hi) that differs by a single
added bond, compare the measured dDeltaG to the sequence-predicted dDeltaG of the added bond.

  ddg_meas  = Exp_dG(hi) - Exp_dG(lo)                         [experimental; UNCHANGED by salt]
  ddg_pred  = Predicted_deltaG(hi) - Predicted_deltaG(lo)     [= salt-corrected NN dG of the added bond]
  efficiency = ddg_meas / ddg_pred

The pairing (lo, hi, regime, added_nt) and ddg_meas are experimental metadata and are read from the
committed pairs file unchanged; only the salt-dependent columns ddg_pred and efficiency are recomputed
from the corrected additive prediction (Predicted_deltaG in Filtered_DF).  After the salt fix the
additive Predicted_deltaG uses a single consistent per-duplex convention, so ddg_pred is simply the
difference of the two interfaces' predictions.  The displayed panel (test2_measured_only.svg) uses
ddg_meas only and is unaffected.

Run from anywhere.  Rewrites test2_add_one_bond_pairs.csv in place.
"""
import os, numpy as np, pandas as pd
PROJ=__REPO__; os.chdir(PROJ)
OUT="Figure2_NegativeCooperativity/Fig2d_ddG"
PAIRS=os.path.join(OUT,"test2_add_one_bond_pairs.csv")

# corrected additive NN prediction per construct code (sequence-derived; identical across replicates)
hs=pd.read_csv("Filtered_DF_HS.csv"); hs["code"]=hs["Code"].astype(str).str.partition(".")[0]
pred=hs.groupby("code")["Predicted_deltaG"].median()

d=pd.read_csv(PAIRS)
missing=sorted(set(d.lo)|set(d.hi)) ; missing=[c for c in missing if c not in pred.index]
if missing:
    raise SystemExit(f"codes missing from Filtered_DF_HS: {missing[:6]} ...")

d["ddg_pred"]=d["hi"].map(pred).values - d["lo"].map(pred).values
d["efficiency"]=d["ddg_meas"]/d["ddg_pred"]
d.to_csv(PAIRS,index=False)

for reg,sub in d.groupby("regime"):
    print(f"regime {reg}: n={len(sub)}  median ddg_pred={sub.ddg_pred.median():+.3f}  "
          f"median efficiency={sub.efficiency.median():+.3f}")
print("rewrote",PAIRS,"(ddg_pred, efficiency recomputed; ddg_meas unchanged)")
