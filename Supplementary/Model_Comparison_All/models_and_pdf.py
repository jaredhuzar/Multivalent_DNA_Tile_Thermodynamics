import os as _os
def _find_repo():
    d=_os.path.dirname(_os.path.abspath(__file__))
    while d!='/' and not _os.path.exists(_os.path.join(d,'Filtered_DF_Tm.csv')):
        d=_os.path.dirname(d)
    return d
__REPO__=_find_repo()
import os, numpy as np, pandas as pd, nbformat, warnings
from scipy.optimize import least_squares
from scipy.stats import spearmanr
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
warnings.filterwarnings("ignore")
os.chdir(__REPO__)
OUT="Supplementary/Model_Comparison_All"

# ---- constants / NN thermo / TH ----
R_CAL=1.98720425864; C_TILE=4.5e-9; EDGE=88.0; RLN=R_CAL*np.log(C_TILE/2.0)
ELL_P_SS,RISE_SS=1.0,0.63; NA=6.02214076e23; C0_nm3=1e-24*NA; T=310.15
NN_dH={'AA':-7.6,'AT':-7.2,'AC':-8.4,'AG':-7.8,'TA':-7.2,'TT':-7.6,'TC':-8.2,'TG':-8.5,'CA':-8.5,'CT':-7.8,'CC':-8.0,'CG':-10.6,'GA':-8.2,'GT':-8.4,'GC':-9.8,'GG':-8.0}
NN_dS={'AA':-22.2,'AT':-20.4,'AC':-22.4,'AG':-21.0,'TA':-21.3,'TT':-22.2,'TC':-22.2,'TG':-22.7,'CA':-22.7,'CT':-21.0,'CC':-19.9,'CG':-27.2,'GA':-22.2,'GT':-22.4,'GC':-24.4,'GG':-19.9}
dH_iGC,dH_iAT,dS_iGC,dS_iAT=0.1,2.3,-2.8,4.1
def duplex(seq):
    seq=str(seq).upper(); dH=sum(NN_dH[seq[i:i+2]] for i in range(len(seq)-1)); dS=sum(NN_dS[seq[i:i+2]] for i in range(len(seq)-1))
    for e in (seq[0],seq[-1]): dH+=(dH_iGC if e in 'GC' else dH_iAT); dS+=(dS_iGC if e in 'GC' else dS_iAT)
    return dH,dS
xl=pd.ExcelFile('DNA Sequences.xlsx')
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
SEQ={'5':_parse(pd.read_excel(xl,'5nt',header=None)),'7':_parse(pd.read_excel(xl,'7nt',header=None)),
     '9':_parse(pd.read_excel(xl,'9nt',header=None)),'A':_parse(pd.read_excel(xl,'A',header=None)),
     'B':_parse(pd.read_excel(xl,'B',header=None)),'6':_parse(pd.read_csv('DNA Sequences - 6nt bonds.csv',header=None))}
TH={t:{h:duplex(s) for h,s in d.items()} for t,d in SEQ.items()}

# ---- WLC tethering entropy; 6nt is ASYMMETRIC (8nt one tile, 7nt other) ----
def _msq(nt):
    if nt<=0: return 1e-9
    Lc=nt*RISE_SS; x=Lc/ELL_P_SS; return 2*ELL_P_SS*Lc*(1-(ELL_P_SS/Lc)*(1-np.exp(-x)))
def _V(nt):
    if nt<=0: nt=0.5
    return (4/3)*np.pi*_msq(nt)**1.5
SS={'5':(7,2),'X':(7,2),'7':(8,1),'9':(9,0),'A':(8,1),'B':(8,1)}   # symmetric (n_unbound, tail)
def dS_teth(lab):
    if lab=='6':                       # asymmetric: strands 8nt & 7nt -> tails 2nt & 1nt
        return R_CAL*np.log(_V(1.5)/(_V(8)*_V(7)*C0_nm3))   # tail=mean(2,1); two distinct free ends
    nu,nb=SS.get(lab,(7,2)); return R_CAL*np.log(_V(nb)/(_V(nu)**2*C0_nm3))

# ---- 235-interface dataset ----
tm=pd.read_csv('Filtered_DF_Tm.csv'); tm['code']=tm['Samples'].astype(str).str.partition('.')[0]
agg=tm.groupby('code').agg(Tm=('Tm','median')).reset_index()
def bonds(code):
    labs=[]; [labs.extend(['5' if ch=='X' else ch]*2) for ch in code]
    z=(np.arange(len(labs))+0.5)/len(labs)*EDGE; out=[]
    for s,lab in enumerate(labs):
        if lab=='0': continue
        th=TH.get(lab,{}).get(s+1)
        if th is None: return None
        out.append((z[s],th[0],th[1],0 if lab in '7AB' else 1,lab))
    return out
agg['B']=agg['code'].apply(bonds); agg=agg[agg['B'].notna()&agg['B'].apply(lambda b:len(b)>0)].reset_index(drop=True)
Y=agg['Tm'].values.astype(float); n=len(agg)
# precompute per construct
D=[]
for b in agg['B']:
    z=np.array([x[0] for x in b]); N=len(z)
    md=np.array([np.mean(np.abs(np.delete(z,i)-z[i])) for i in range(N)]) if N>1 else np.zeros(1)
    mn=np.array([np.min(np.abs(np.delete(z,i)-z[i])) for i in range(N)]) if N>1 else np.zeros(1)
    D.append(dict(md=md,mn=mn,dH=np.array([x[1] for x in b]),dS=np.array([x[2] for x in b]),
                  flag=np.array([x[3] for x in b]),teth=np.array([dS_teth(x[4]) for x in b]),N=N,
                  I=float(np.sum((z-z.mean())**2))))
def tmK(dH,dS): return dH*1000.0/dS-273.15

def P_mean(d,p):    lam,s,b=p; g=np.exp(-d['md']/lam); return tmK((g*d['dH']).sum(),(g*(d['dS']+s)).sum()+b+RLN)
def P_Nscale(d,p):  lam,s,b=p; g=np.exp(-d['md']/lam); return tmK((g*d['dH']).sum(),(g*d['dS']).sum()+d['N']*s+b+RLN)
def P_hOnly(d,p):   lam,s,b=p; g=np.exp(-d['md']/lam); return tmK((g*d['dH']).sum(),d['dS'].sum()+d['N']*s+b+RLN)
def P_nn(d,p):      lam,s,b=p; g=np.exp(-d['mn']/lam); return tmK((g*d['dH']).sum(),(g*(d['dS']+s)).sum()+b+RLN)
def P_pertype(d,p): lam,ss,sl,b=p; g=np.exp(-d['md']/lam); pen=np.where(d['flag'],ss,sl); return tmK((g*d['dH']).sum(),(g*(d['dS']+pen)).sum()+b+RLN)
def P_split(d,p):   lh,ls,s,b=p; gH=np.exp(-d['md']/lh); gS=np.exp(-d['md']/ls); return tmK((gH*d['dH']).sum(),(gS*(d['dS']+s)).sum()+b+RLN)
def P_rot(d,p):     lam,s,b,k=p; g=np.exp(-d['md']/lam); return tmK((g*d['dH']).sum(),(g*d['dS']).sum()+d['N']*s+b+RLN-k*R_CAL*np.log(1+d['I']))
def P_g1(d,p):      s,b=p; return tmK(d['dH'].sum(),d['dS'].sum()+d['N']*s+b+RLN)
def P_g1rot(d,p):   s,b,k=p; return tmK(d['dH'].sum(),d['dS'].sum()+d['N']*s+b+RLN-k*R_CAL*np.log(1+d['I']))
def P_gwteth(d,p):  lam,b=p; g=np.exp(-d['md']/lam); return tmK((g*d['dH']).sum(),(g*(d['dS']+d['teth'])).sum()+b+RLN)
def P_teth(d,p):    lam,b=p; g=np.exp(-d['md']/lam); return tmK((g*d['dH']).sum(),(d['dS']+d['teth']).sum()+b+RLN)

MODELS=[("A1 mean (g-weighted ΔS_Bond)",P_mean,[120,-8,5],([1,-60,-200],[1000,5,400])),
 ("A2 N-scaled ΔS_Bond",P_Nscale,[120,-8,5],([1,-60,-200],[1000,5,400])),
 ("A3 g on ΔH only",P_hOnly,[120,-8,5],([1,-60,-200],[1000,5,400])),
 ("A4 nearest-neighbor g",P_nn,[120,-8,5],([1,-60,-200],[1000,5,400])),
 ("A5 per-type ΔS_Bond",P_pertype,[120,-8,-8,5],([1,-60,-60,-200],[1000,5,5,400])),
 ("A6 split λ_H/λ_S",P_split,[120,120,-8,5],([1,1,-60,-200],[1000,1000,5,400])),
 ("A7 occ.+rotational",P_rot,[120,-8,5,.1],([1,-60,-200,-5],[1000,5,400,20])),
 ("B1 g=1 (no weighting)",P_g1,[-8,5],([-60,-200],[5,400])),
 ("B2 g=1 + rotational",P_g1rot,[-8,5,.1],([-60,-200,-5],[5,400,20])),
 ("B3 g-weighted + WLC tethering",P_gwteth,[120,5],([1,-200],[1000,400])),
 ("B4 WLC tethering only",P_teth,[120,5],([1,-200],[1000,400]))]
def fit(P,p0,bn,idx=None):
    I=range(n) if idx is None else idx
    return least_squares(lambda p:[P(D[i],p)-Y[i] for i in I],p0,bounds=bn).x
def ev(P,p0,bn):
    p=fit(P,p0,bn); fp=np.array([P(D[i],p) for i in range(n)])
    oof=np.full(n,np.nan)
    for tr,te in KFold(10,shuffle=True,random_state=0).split(np.arange(n)):
        q=fit(P,p0,bn,idx=tr)
        for j in te: oof[j]=P(D[j],q)
    return dict(cvRMS=np.sqrt(np.mean((oof-Y)**2)),cvR2=r2_score(Y,oof),cvrho=spearmanr(Y,oof)[0],
                fitR2=r2_score(Y,fp),params=np.round(p,2))
rows=[]
for nm,P,p0,bn in MODELS:
    r=ev(P,p0,bn); rows.append(dict(model=nm,cvRMS=round(r['cvRMS'],3),cvR2=round(r['cvR2'],3),
                                    cvrho=round(r['cvrho'],3),fitR2=round(r['fitR2'],3),params=str(r['params'])))
    print(f"{nm:34s} cvRMS={r['cvRMS']:.2f} cvrho={r['cvrho']:.3f}")
df=pd.DataFrame(rows); df.to_csv(os.path.join(OUT,"model_comparison_all_235.csv"),index=False)
df.to_pickle('/tmp/mc_rows.pkl')
print("saved CSV")
