import os as _os
def _find_repo():
    d=_os.path.dirname(_os.path.abspath(__file__))
    while d!='/' and not _os.path.exists(_os.path.join(d,'Filtered_DF_Tm.csv')):
        d=_os.path.dirname(d)
    return d
__REPO__=_find_repo()
"""Tested-models performance table (.docx).  Every number is read from
model_comparison_all_235.csv (written by unified_models.py) -- nothing is hard-coded, so the
docx can never desync from the CSV / PDF.  The fitted-parameter column is labelled per model
from a fixed parameter-name schema and the numeric values parsed out of the CSV `params` field."""
import os, re, numpy as np, pandas as pd
os.chdir(__REPO__)
from docx import Document
from docx.shared import Pt, RGBColor, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.section import WD_ORIENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

OUT="Supplementary/TableS1_Model_Comparison_All"
df=pd.read_csv(os.path.join(OUT,"model_comparison_all_235.csv"))

# per-model fitted-parameter labels, in the same order unified_models.py fits them
PNAMES={
 "A1":["λ","ΔS_Bond","base"], "A2":["λ","ΔS_Bond","base"], "A3":["λ","ΔS_Bond","base"],
 "A4":["λ","ΔS_Bond","base"], "A5":["λ","ΔS_short","ΔS_long","base"],
 "A6":["λ_H","λ_S","ΔS_Bond","base"], "A7":["λ","ΔS_Bond","base","k_rot"],
 "B1":["ΔS_Bond","base"], "B2":["ΔS_Bond","base","k_rot"],
 "B3":["λ","base"], "B4":["λ","base"],
}
def fmt_params(model, raw):
    vals=[float(x) for x in re.findall(r"-?\d+\.?\d*(?:[eE][-+]?\d+)?", str(raw))]
    names=PNAMES.get(model.split()[0], [f"p{i}" for i in range(len(vals))])
    def _fmt(n,v):
        if n.startswith("λ") and v>=500: return f"{n}→∞ (g→1, unidentified)"   # >> 60 nm interface
        return f"{n}={v:g}"
    return ", ".join(_fmt(n,v) for n,v in zip(names, vals))
_best = str(df.loc[df["cvRMS"].astype(float).idxmin(),"model"])   # data-driven "best" = lowest cv RMS

headers=["Model","cv RMS (°C)","cv R²","cv ρ","fit R²","Fitted parameters"]
widths=[Inches(2.5),Inches(1.0),Inches(0.8),Inches(0.8),Inches(0.8),Inches(3.6)]
rows=[[str(r["model"]), f"{r['cvRMS']:.2f}", f"{r['cvR2']:.3f}", f"{r['cvrho']:.3f}",
       f"{r['fitR2']:.3f}", fmt_params(r["model"], r["params"])] for _,r in df.iterrows()]

def shade(cell,hexcolor):
    tcPr=cell._tc.get_or_add_tcPr(); sh=OxmlElement('w:shd')
    sh.set(qn('w:val'),'clear'); sh.set(qn('w:color'),'auto'); sh.set(qn('w:fill'),hexcolor); tcPr.append(sh)
def setcell(cell,text,w,bold=False,white=False,center=False):
    cell.width=w; p=cell.paragraphs[0]; p.alignment=WD_ALIGN_PARAGRAPH.CENTER if center else WD_ALIGN_PARAGRAPH.LEFT
    r=p.add_run(text); r.font.size=Pt(8.5); r.font.bold=bold
    if white: r.font.color.rgb=RGBColor(0xFF,0xFF,0xFF)

doc=Document()
sec=doc.sections[0]; sec.orientation=WD_ORIENT.LANDSCAPE
sec.page_width,sec.page_height=Inches(11),Inches(8.5)
sec.left_margin=Inches(0.6); sec.right_margin=Inches(0.6); sec.top_margin=Inches(0.7)

doc.add_heading("Phenomenological Tm models — performance",level=1)
sub=doc.add_paragraph(); r=sub.add_run("10-fold cross-validated on the canonical interface set. "
    "Family A = occupancy-weighted; Family B = elastic / mechanical-coupling. ΔS and k_rot in cal/mol·K, λ in nm. "
    "SantaLucia 1998 NN table; 6 nm site pitch.")
r.italic=True; r.font.size=Pt(9); r.font.color.rgb=RGBColor(0x55,0x55,0x55)

tbl=doc.add_table(rows=1,cols=len(headers)); tbl.style="Table Grid"; tbl.autofit=False
for j,hname in enumerate(headers):
    c=tbl.rows[0].cells[j]; shade(c,"4C78A8"); setcell(c,hname,widths[j],bold=True,white=True,center=(j not in (0,5)))
for row in rows:
    cells=tbl.add_row().cells; best=(row[0]==_best)
    for j,val in enumerate(row):
        if best: shade(cells[j],"DCE8F5")
        setcell(cells[j],val,widths[j],bold=(best and j==0),center=(j not in (0,5)))

doc.add_paragraph()
_Arows=df[df["model"].str.startswith(("A1","A2","A3","A5","A6","A7"))]["cvRMS"].astype(float)
foot=doc.add_paragraph(); r=foot.add_run(
    f"Lowest cv RMS: {_best} (highlighted). The occupancy-weighted models (A1-A3, A5-A7) are statistically "
    f"comparable (cv RMS {_Arows.min():.2f}-{_Arows.max():.2f}); A2 (N-scaled per-bond penalty) is the parsimonious "
    "occupancy model used for the Fig 4 theory curves. The elastic Green's-function (rotational) term adds only a "
    "small gain (A2 to A7). For the WLC-tethering models (B3/B4), lambda is unidentified (g->1) and sits at/near the "
    "bound. B3/B4 use the corrected asymmetric 6nt bond (8 nt / 7 nt strands). All values read directly from "
    "model_comparison_all_235.csv.")
r.font.size=Pt(8); r.font.color.rgb=RGBColor(0x55,0x55,0x55)

doc.save(os.path.join(OUT,"phenomenological_models_table.docx"))
print("wrote phenomenological_models_table.docx from model_comparison_all_235.csv (",len(rows),"rows )")
