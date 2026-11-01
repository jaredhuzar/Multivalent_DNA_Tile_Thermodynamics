import os as _os
def _find_repo():
    d=_os.path.dirname(_os.path.abspath(__file__))
    while d!='/' and not _os.path.exists(_os.path.join(d,'Filtered_DF_Tm.csv')):
        d=_os.path.dirname(d)
    return d
__REPO__=_find_repo()
import numpy as np, pandas as pd, warnings
from scipy.optimize import least_squares
from scipy.stats import spearmanr
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
warnings.filterwarnings("ignore")
import os; os.chdir(__REPO__)

# ---------- constants ----------
R_CAL=1.98720425864; C_TILE=4.5e-9; EDGE=88.0; NSLOT=22; RLN=R_CAL*np.log(C_TILE/2.0)
ELL_P_SS, RISE_SS = 1.0, 0.63; NA=6.02214076e23; C0_nm3=1e-24*NA

NN_dH={'AA':-7.6,'AT':-7.2,'AC':-8.4,'AG':-7.8,'TA':-7.2,'TT':-7.6,'TC':-8.2,'TG':-8.5,'CA':-8.5,'CT':-7.8,'CC':-8.0,'CG':-10.6,'GA':-8.2,'GT':-8.4,'GC':-9.8,'GG':-8.0}
NN_dS={'AA':-22.2,'AT':-20.4,'AC':-22.4,'AG':-21.0,'TA':-21.3,'TT':-22.2,'TC':-22.2,'TG':-22.7,'CA':-22.7,'CT':-21.0,'CC':-19.9,'CG':-27.2,'GA':-22.2,'GT':-22.4,'GC':-24.4,'GG':-19.9}
dH_iGC,dH_iAT,dS_iGC,dS_iAT=0.1,2.3,-2.8,4.1
def duplex_thermo(seq):
    seq=str(seq).upper(); dH=sum(NN_dH[seq[i:i+2]] for i in range(len(seq)-1)); dS=sum(NN_dS[seq[i:i+2]] for i in range(len(seq)-1))
    for e in (seq[0],seq[-1]):
        if e in 'GC': dH+=dH_iGC; dS+=dS_iGC
        else: dH+=dH_iAT; dS+=dS_iAT
    return dH,dS

# ---------- sequences / per-type thermo ----------
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
TH={t:{h:duplex_thermo(s) for h,s in d.items()} for t,d in SEQ.items()}

# ---------- canonical 235-interface dataset ----------
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
agg['B']=agg['code'].apply(bonds)
agg=agg[agg['B'].notna() & agg['B'].apply(lambda b:len(b)>0)].reset_index(drop=True)
Bl=list(agg['B']); Y=agg['Tm'].values.astype(float); n=len(Bl)
print(f"canonical set: n={n} interfaces\n")

# ---------- weights / physical terms ----------
def gk(b,lam,mode='mean'):
    z=np.array([x[0] for x in b])
    if len(z)==1: return np.array([1.0])
    if mode=='mean': d=np.array([np.mean(np.abs(np.delete(z,i)-z[i])) for i in range(len(z))])
    else:            d=np.array([np.min(np.abs(np.delete(z,i)-z[i])) for i in range(len(z))])
    return np.exp(-d/lam)
def inertia(b): z=np.array([x[0] for x in b]); return np.sum((z-z.mean())**2)
def _msq(nt):
    if nt<=0: return 1e-9
    Lc=nt*RISE_SS; x=Lc/ELL_P_SS; return 2*ELL_P_SS*Lc*(1-(ELL_P_SS/Lc)*(1-np.exp(-x)))
def _Vsph(nt):
    if nt<=0: nt=0.5
    return (4/3)*np.pi*_msq(nt)**1.5
SS_STATES={'5':(7,2),'X':(7,2),'6':(7,1),'7':(8,1),'9':(9,0),'A':(8,1),'B':(8,1)}
def dS_teth(lab):
    nu,nb=SS_STATES.get(lab,(7,2)); return R_CAL*np.log(_Vsph(nb)/(_Vsph(nu)**2*C0_nm3))

TK=273.15
def _tm(dH,dS): return (dH*1000)/dS - TK          # dH kcal, dS cal/K; returns C

# ---------- model predictors (b, p) -> Tm(C) ----------
def P_mean(b,p):      lam,dSb,base=p; g=gk(b,lam);   return _tm(sum(gi*x[1] for gi,x in zip(g,b)), sum(gi*(x[2]+dSb) for gi,x in zip(g,b))+base+RLN)
def P_Nscale(b,p):    lam,dSb,base=p; g=gk(b,lam);   return _tm(sum(gi*x[1] for gi,x in zip(g,b)), sum(gi*x[2] for gi,x in zip(g,b))+len(b)*dSb+base+RLN)
def P_hOnly(b,p):     lam,dSb,base=p; g=gk(b,lam);   return _tm(sum(gi*x[1] for gi,x in zip(g,b)), sum(x[2] for x in b)+len(b)*dSb+base+RLN)
def P_nn(b,p):        lam,dSb,base=p; g=gk(b,lam,'nn'); return _tm(sum(gi*x[1] for gi,x in zip(g,b)), sum(gi*(x[2]+dSb) for gi,x in zip(g,b))+base+RLN)
def P_pertype(b,p):   lam,ds,dl,base=p; g=gk(b,lam); return _tm(sum(gi*x[1] for gi,x in zip(g,b)), sum(gi*(x[2]+(ds if x[3] else dl)) for gi,x in zip(g,b))+base+RLN)
def P_splitlam(b,p):  lamH,lamS,dSb,base=p; gH=gk(b,lamH); gS=gk(b,lamS); return _tm(sum(gi*x[1] for gi,x in zip(gH,b)), sum(gi*(x[2]+dSb) for gi,x in zip(gS,b))+base+RLN)
def P_Nscale_rot(b,p):lam,dSb,base,kr=p; g=gk(b,lam); return _tm(sum(gi*x[1] for gi,x in zip(g,b)), sum(gi*x[2] for gi,x in zip(g,b))+len(b)*dSb+base+RLN-kr*R_CAL*np.log(1+inertia(b)))
def P_g1(b,p):        dSb,base=p;                     return _tm(sum(x[1] for x in b), sum(x[2] for x in b)+len(b)*dSb+base+RLN)
def P_g1_rot(b,p):    dSb,base,kr=p;                  return _tm(sum(x[1] for x in b), sum(x[2] for x in b)+len(b)*dSb+base+RLN-kr*R_CAL*np.log(1+inertia(b)))
def P_gw_tether(b,p): lam,base=p; g=gk(b,lam);        return _tm(sum(gi*x[1] for gi,x in zip(g,b)), sum(gi*(x[2]+dS_teth(x[4])) for gi,x in zip(g,b))+base+RLN)
def P_tether_only(b,p):lam,base=p; g=gk(b,lam);       return _tm(sum(gi*x[1] for gi,x in zip(g,b)), sum(x[2]+dS_teth(x[4]) for x in b)+base+RLN)

MODELS=[
 ("A1 mean (g-weighted ΔS_Bond)",        P_mean,     [120,-8,5],  ([1,-60,-200],[1000,5,400])),
 ("A2 N-scaled ΔS_Bond  [best occ.]",    P_Nscale,   [120,-8,5],  ([1,-60,-200],[1000,5,400])),
 ("A3 g on ΔH only",                     P_hOnly,    [120,-8,5],  ([1,-60,-200],[1000,5,400])),
 ("A4 nearest-neighbor g",               P_nn,       [120,-8,5],  ([1,-60,-200],[1000,5,400])),
 ("A5 per-type ΔS_Bond",                 P_pertype,  [120,-8,-8,5],([1,-60,-60,-200],[1000,5,5,400])),
 ("A6 split λ_H/λ_S",                    P_splitlam, [120,120,-8,5],([1,1,-60,-200],[1000,1000,5,400])),
 ("A7 occ. + rotational (Green's fn)",   P_Nscale_rot,[120,-8,5,.1],([1,-60,-200,-5],[1000,5,400,20])),
 ("B1 g=1 (no arrangement weighting)",   P_g1,       [-8,5],      ([-60,-200],[5,400])),
 ("B2 g=1 + rotational",                 P_g1_rot,   [-8,5,.1],   ([-60,-200,-5],[5,400,20])),
 ("B3 g-weighted + WLC tethering",       P_gw_tether,[120,5],     ([1,-200],[1000,400])),
 ("B4 WLC tethering only (derived)",     P_tether_only,[120,5],   ([1,-200],[1000,400])),
]

def fit(pred,p0,bnds,idx=None):
    I=range(n) if idx is None else idx
    return least_squares(lambda p:[pred(Bl[i],p)-Y[i] for i in I],p0,bounds=bnds).x
def evaluate(pred,p0,bnds):
    p=fit(pred,p0,bnds); fp=np.array([pred(b,p) for b in Bl])
    frms=np.sqrt(np.mean((fp-Y)**2)); fr2=r2_score(Y,fp); frho=spearmanr(Y,fp)[0]
    oof=np.full(n,np.nan)
    for tr,te in KFold(10,shuffle=True,random_state=0).split(np.arange(n)):
        q=fit(pred,p0,bnds,idx=tr)
        for j in te: oof[j]=pred(Bl[j],q)
    cvr=np.sqrt(np.mean((oof-Y)**2)); cvr2=r2_score(Y,oof); cvrho=spearmanr(Y,oof)[0]
    return dict(fitRMS=frms,fitR2=fr2,fitrho=frho,cvRMS=cvr,cvR2=cvr2,cvrho=cvrho,params=np.round(p,2))

rows=[]
print(f"{'model':38s} {'cvRMS':>6} {'cvR2':>6} {'cvrho':>6} {'fitR2':>6}")
print("-"*72)
for name,pred,p0,bnds in MODELS:
    r=evaluate(pred,p0,bnds)
    rows.append(dict(model=name,**{k:(round(v,3) if np.isscalar(v) else str(v)) for k,v in r.items()}))
    print(f"{name:38s} {r['cvRMS']:6.2f} {r['cvR2']:6.3f} {r['cvrho']:6.3f} {r['fitR2']:6.3f}")
out=pd.DataFrame(rows)
os.makedirs("Supplementary/Model_Comparison_All",exist_ok=True)
out.to_csv("Supplementary/Model_Comparison_All/model_comparison_all_235.csv",index=False)
print("\nsaved Supplementary/Model_Comparison_All/model_comparison_all_235.csv")
