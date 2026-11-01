import os as _os
def _find_repo():
    d=_os.path.dirname(_os.path.abspath(__file__))
    while d!='/' and not _os.path.exists(_os.path.join(d,'Filtered_DF_Tm.csv')):
        d=_os.path.dirname(d)
    return d
__REPO__=_find_repo()
"""Render phenomenological-model formulas + performance table to a PDF."""
import os, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
os.chdir(__REPO__)
OUT="Supplementary/Model_Comparison_All"
df=pd.read_csv(os.path.join(OUT,"model_comparison_all_235.csv"))

CORE=(r"$T_m(^\circ\mathrm{C})=\frac{1000\,\sum_k g_k\,\Delta H_{NN,k}}{\Delta S_{\mathrm{eff}}}-273.15,"
      r"\qquad g_k=e^{-\bar d_k/\lambda}$")
NOTE=(r"$\bar d_k$ = mean distance of bond $k$ to the others;  every $\Delta S_{\mathrm{eff}}$ also includes "
      r"$+R\ln(C_T/2)\approx-39.6$ cal/mol$\cdot$K.  $N=$ # bonds (2 per site).")
FORMS=[
 (r"A1  mean",                    r"$\Delta S_{\mathrm{eff}}=\sum_k g_k(\Delta S_{NN,k}+\Delta S_{Bond})+\mathrm{base}$"),
 (r"A2  N-scaled $\Delta S_{Bond}$ (best)", r"$\Delta S_{\mathrm{eff}}=\sum_k g_k\Delta S_{NN,k}+N\,\Delta S_{Bond}+\mathrm{base}$"),
 (r"A3  g on $\Delta H$ only",     r"$\Delta S_{\mathrm{eff}}=\sum_k \Delta S_{NN,k}+N\,\Delta S_{Bond}+\mathrm{base}$  (g applied only in $\Delta H$)"),
 (r"A4  nearest-neighbor g",       r"same as A1 with $g_k=e^{-d_{\min,k}/\lambda}$"),
 (r"A5  per-type $\Delta S_{Bond}$", r"$\Delta S_{\mathrm{eff}}=\sum_k g_k(\Delta S_{NN,k}+\Delta S_{Bond}^{\,\mathrm{type}})+\mathrm{base}$"),
 (r"A6  split $\lambda_H/\lambda_S$", r"$\Delta H$ uses $g(\lambda_H)$, $\Delta S$ uses $g(\lambda_S)$"),
 (r"A7  occ. + rotational",        r"A2 $-\,k_{rot}R\ln(1+\sum_k(z_k-\bar z)^2)$"),
 (r"B1  g=1 (no weighting)",       r"$\Delta S_{\mathrm{eff}}=\sum_k \Delta S_{NN,k}+N\,\Delta S_{Bond}+\mathrm{base}$   ($g_k\equiv1$)"),
 (r"B2  g=1 + rotational",         r"B1 $-\,k_{rot}R\ln(1+\sum_k(z_k-\bar z)^2)$"),
 (r"B3  g-weighted + WLC tethering", r"$\Delta S_{\mathrm{eff}}=\sum_k g_k(\Delta S_{NN,k}+\Delta S_{teth}^{\,\mathrm{type}})+\mathrm{base}$"),
 (r"B4  WLC tethering only",       r"$\Delta S_{\mathrm{eff}}=\sum_k(\Delta S_{NN,k}+\Delta S_{teth}^{\,\mathrm{type}})+\mathrm{base}$"),
]
ELAST=[r"Rotational (elastic Green's fn):  $\Delta S_{rot}=-k_{rot}\,R\ln(1+\sum_k(z_k-\bar z)^2)$",
 r"WLC tethering (derived):  $\Delta S_{teth}^{\,\mathrm{type}}=R\ln\frac{V_2}{V_1^{2}\,C_0}$,  $V=\frac{4}{3}\pi\langle r^2\rangle^{3/2}$  (WLC, $\ell_p$=1 nm, 0.63 nm/nt)",
 r"6nt is asymmetric (8 nt / 7 nt strands):  $\Delta S_{teth}^{6}=R\ln\frac{V(1.5)}{V(8)\,V(7)\,C_0}$"]

with PdfPages(os.path.join(OUT,"phenomenological_models.pdf")) as pdf:
    # page 1: formulas
    fig=plt.figure(figsize=(8.5,11)); fig.text(0.5,0.955,"Phenomenological $T_m$ models — formulas",ha="center",fontsize=15,weight="bold")
    fig.text(0.5,0.915,"Shared two-state core",ha="center",fontsize=11,style="italic")
    fig.text(0.5,0.885,CORE,ha="center",fontsize=13)
    fig.text(0.5,0.855,NOTE,ha="center",fontsize=8.5)
    fig.text(0.07,0.815,"Family A — occupancy-weighted  |  Family B — elastic / mechanical-coupling",fontsize=10,weight="bold")
    y=0.785
    for name,form in FORMS:
        fig.text(0.07,y,name,fontsize=10,weight="bold"); fig.text(0.40,y,form,fontsize=10.5); y-=0.047
    fig.text(0.07,y-0.01,"Mechanical terms",fontsize=10,weight="bold"); y-=0.055
    for e in ELAST: fig.text(0.07,y,e,fontsize=9.5); y-=0.042
    pdf.savefig(fig); plt.close(fig)
    # page 2: table
    fig,ax=plt.subplots(figsize=(11,8.5)); ax.axis("off")
    ax.set_title("Performance (10-fold CV, n = 235 interfaces)",fontsize=14,weight="bold",pad=18)
    cols=["model","cvRMS","cvR2","cvrho","fitR2","params"]
    tbl=ax.table(cellText=df[cols].values,colLabels=["Model","cv RMS (°C)","cv R²","cv ρ","fit R²","fitted params"],
                 loc="center",cellLoc="center",colWidths=[0.27,0.09,0.08,0.08,0.08,0.40])
    tbl.auto_set_font_size(False); tbl.set_fontsize(8.5); tbl.scale(1,1.6)
    for (r,c),cell in tbl.get_celld().items():
        if r==0: cell.set_facecolor("#4C78A8"); cell.set_text_props(color="white",weight="bold")
        elif df.iloc[r-1]["model"].startswith("A2"): cell.set_facecolor("#DCE8F5")
    pdf.savefig(fig,bbox_inches="tight"); plt.close(fig)
print("wrote",os.path.join(OUT,"phenomenological_models.pdf"))
