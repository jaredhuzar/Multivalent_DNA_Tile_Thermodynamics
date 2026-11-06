import os as _os
def _find_repo():
    d=_os.path.dirname(_os.path.abspath(__file__))
    while d!='/' and not _os.path.exists(_os.path.join(d,'Filtered_DF_Tm.csv')):
        d=_os.path.dirname(d)
    return d
__REPO__=_find_repo()
"""
Rebuild Experimental_Trends/construct_medians.csv from the AUTHORITATIVE droop-corrected fits
(the same fit that makes Filtered_DF), so the geometry-controlled OLS panels and coop_unit use the
SAME Tm / dG / construct set as the rest of the paper -- NOT the separate Boltzmann-sigmoid pipeline
(Experimental_Trends_Analysis.ipynb) that produced the old construct_medians.

  Tm            : per-construct median from Filtered_DF_Tm   (fit_Not_Normalized: droop-corrected
                  sigmoid midpoint + SE/Reliable/plateau QC)
  dG37/dH/dS    : per-construct median from Filtered_DF_HS   (droop-corrected van't Hoff)
  Num_Bonds, MeanGap_nm, span_nm, bond_len, GC : functions of the construct code (helix pitch 8.8 nm;
                  GC/length from the repo sequence files) -- identical values to the old file.

Self-contained: reads only repo files (Filtered_DF_Tm.csv, Filtered_DF_HS.csv, DNA Sequences.xlsx,
DNA Sequences - 6nt bonds.csv).  Run before the Supplementary OLS / coop_unit scripts.
"""
import os, numpy as np, pandas as pd
PROJ=__REPO__; os.chdir(PROJ)

# ---- sequence table -> GC fraction and bond length (same definition as Experimental_Trends) ----
xl=pd.ExcelFile("DNA Sequences.xlsx")
def _hdr(raw):
    for r in range(len(raw)):
        row=[str(x) for x in raw.iloc[r].tolist()]
        if 'Tile 1' in row: return r,{v:i for i,v in enumerate(row)}
def _parse(raw):
    hdr,cols=_hdr(raw); d={}
    for r in range(hdr+1,len(raw)):
        v=raw.iloc[r]
        try: h=int(float(v[cols['#']]))
        except: continue
        t=str(v[cols['Tile 1']]).strip().upper()
        if t and t!='NAN': d[h]=t
    return d
SEQ={k:_parse(pd.read_excel(xl,sheet_name=s,header=None)) for k,s in
     [('5','5nt'),('7','7nt'),('9','9nt'),('A','A'),('B','B')]}
SEQ['6']=_parse(pd.read_csv("DNA Sequences - 6nt bonds.csv",header=None))
LEN={'5':5,'X':5,'6':6,'7':7,'9':9,'A':7,'B':7}
def gc_and_len(code):
    gcs,lens=[],[]
    for i,ch in enumerate(code):
        if ch=='0': continue
        c='5' if ch=='X' else ch; lens.append(LEN.get(ch,np.nan))
        s=SEQ.get(c,{}).get(2*i+1)                 # odd physical helix
        if s: gcs.append((s.count('G')+s.count('C'))/len(s))
    return (np.mean(gcs) if gcs else np.nan),(np.mean(lens) if lens else np.nan)
def geom(code):
    pos=[i for i,ch in enumerate(code) if ch!='0']; nb=len(pos)
    mg=np.mean(np.diff(sorted(pos)))*8.8 if nb>1 else np.nan
    sp=(max(pos)-min(pos))*8.8 if nb>1 else 0.0
    gc,bl=gc_and_len(code)
    return nb,mg,sp,gc,bl

# ---- authoritative per-construct medians ----
tm=pd.read_csv("Filtered_DF_Tm.csv"); tm["code"]=tm["Code"].astype(str).str.partition(".")[0]
tmG=tm.groupby("code").agg(Tm=("Tm","median")).reset_index()
hs=pd.read_csv("Filtered_DF_HS.csv"); hs["code"]=hs["Code"].astype(str).str.partition(".")[0]
hsG=hs.groupby("code").agg(dG37=("Exp_DeltaG","median"),Entropy=("Entropy","median"),
                           Enthalpy=("Enthalpy","median")).reset_index()
# Base the table on the authoritative Tm construct set (Filtered_DF_Tm). dG/dH/dS are looked up from
# Filtered_DF_HS ONLY for those constructs (left join), so no OLS can exceed the Tm set: a construct
# with a valid H/S fit but outside the Tm set is NOT added.  dG37 is NaN where a Tm construct lacks a
# valid van't Hoff dG (the dG-OLS drops those via dropna).
m=pd.merge(tmG,hsG,on="code",how="left")
_g=m["code"].apply(geom)
m["Num_Bonds"]=[g[0] for g in _g]; m["MeanGap_nm"]=[g[1] for g in _g]
m["span_nm"]  =[g[2] for g in _g]; m["GC"]       =[g[3] for g in _g]; m["bond_len"]=[g[4] for g in _g]
m=m[["code","Tm","dG37","Entropy","Enthalpy","Num_Bonds","MeanGap_nm","span_nm","bond_len","GC"]]
os.makedirs("Experimental_Trends",exist_ok=True)
m.to_csv("Experimental_Trends/construct_medians.csv",index=False)
print(f"rebuilt Experimental_Trends/construct_medians.csv on Filtered_DF: {len(m)} constructs "
      f"(Tm from {len(tmG)} Filtered_DF_Tm, dG/dH/dS from {len(hsG)} Filtered_DF_HS)")
