import numpy as np, pandas as pd
from scipy import stats

d=pd.read_csv("Figure3_BondClustering/HS_vs_Distance_7nt_4bonds/hs_vs_distance_7nt_4bonds.csv")
x=d["bond_distance_nm"].values.astype(float)
n=len(d)
print(f"n={n}, distance range {x.min():.1f}-{x.max():.1f} nm\n")

def fit(y,deg):
    X=np.vander(x,deg+1)                      # columns high->low power
    beta,_,_,_=np.linalg.lstsq(X,y,rcond=None)
    yhat=X@beta; rss=np.sum((y-yhat)**2); k=deg+1
    r2=1-rss/np.sum((y-y.mean())**2)
    adj=1-(1-r2)*(n-1)/(n-k)
    aic=n*np.log(rss/n)+2*k; bic=n*np.log(rss/n)+k*np.log(n)
    return beta,rss,r2,adj,aic,bic,k

def loocv(y,deg):
    err=[]
    for i in range(n):
        m=np.ones(n,bool); m[i]=False
        X=np.vander(x[m],deg+1); b,_,_,_=np.linalg.lstsq(X,y[m],rcond=None)
        err.append((y[i]-np.vander(x[i:i+1],deg+1)@b)[0])
    return np.sqrt(np.mean(np.array(err)**2))

for name,y in [("Delta H (kcal/mol)",d["dH"].values.astype(float)),
               ("Delta S (cal/mol/K)",d["dS"].values.astype(float))]:
    print("="*68); print(name); print("="*68)
    b1,rss1,r2_1,adj1,aic1,bic1,k1=fit(y,1)
    b2,rss2,r2_2,adj2,aic2,bic2,k2=fit(y,2)
    # nested F-test: does quadratic term help?
    F=((rss1-rss2)/(k2-k1))/(rss2/(n-k2)); pF=1-stats.f.cdf(F,k2-k1,n-k2)
    cv1,cv2=loocv(y,1),loocv(y,2)
    print(f"  LINEAR    : R2={r2_1:+.4f} adjR2={adj1:+.4f} AIC={aic1:.1f} BIC={bic1:.1f} LOO-RMSE={cv1:.2f}")
    print(f"  QUADRATIC : R2={r2_2:+.4f} adjR2={adj2:+.4f} AIC={aic2:.1f} BIC={bic2:.1f} LOO-RMSE={cv2:.2f}")
    print(f"  dAIC(quad-lin)={aic2-aic1:+.2f}  dBIC={bic2-bic1:+.2f}   (negative favors quadratic)")
    print(f"  F-test quad term: F={F:.3f}, p={pF:.3f}")
    # vertex of parabola y=b2[0]x^2+b2[1]x+b2[2]
    a,bb,c=b2
    if abs(a)>0: vx=-bb/(2*a); print(f"  parabola vertex at x={vx:.1f} nm ({'INSIDE' if x.min()<=vx<=x.max() else 'OUTSIDE'} data range {x.min():.0f}-{x.max():.0f}); curvature a={a:+.4f} ({'down (∩)' if a<0 else 'up (∪)'})")
    print()
