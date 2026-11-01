from docx import Document
from docx.shared import Pt, RGBColor, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.section import WD_ORIENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

rows=[
 ["A1 mean (g-weighted ΔS_Bond)","2.08","0.785","0.876","0.793","λ=53, ΔS_Bond=−10.07, base=16.05"],
 ["A2 N-scaled ΔS_Bond (best)","1.98","0.805","0.884","0.811","λ=139.8, ΔS_Bond=−8.13, base=5.24"],
 ["A3 g on ΔH only","2.05","0.792","0.886","0.795","λ→∞, ΔS_Bond=−7.27, base=4.23"],
 ["A4 nearest-neighbor g","2.56","0.676","0.748","0.679","λ→∞, ΔS_Bond=−10.6, base=5.93"],
 ["A5 per-type ΔS_Bond","2.04","0.793","0.873","0.801","λ=57.5; short=−9.94, long=−9.19; base=11.82"],
 ["A6 split λ_H/λ_S","2.04","0.793","0.878","0.804","λ_H=484, λ_S=754, ΔS_Bond=−8.1, base=5.52"],
 ["A7 occ. + rotational","2.00","0.802","0.881","0.813","λ=113, ΔS_Bond=−7.98, base=−5.34, k_rot=−0.86"],
 ["B1 g=1 (no weighting)","2.56","0.676","0.748","0.679","ΔS_Bond=−10.6, base=5.8"],
 ["B2 g=1 + rotational","2.20","0.761","0.860","0.775","ΔS_Bond=−9.16, base=48.1, k_rot=3.49"],
 ["B3 g-weighted + WLC tethering","4.05","0.190","0.547","0.201","λ=73.8, base=48.3"],
 ["B4 WLC tethering only","5.89","−0.715","0.636","−0.699","λ→∞, base=77.1"],
]
headers=["Model","cv RMS (°C)","cv R²","cv ρ","fit R²","Fitted parameters"]
widths=[Inches(2.5),Inches(1.0),Inches(0.8),Inches(0.8),Inches(0.8),Inches(3.6)]

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
for m in (sec.left_margin,): pass
sec.left_margin=Inches(0.6); sec.right_margin=Inches(0.6); sec.top_margin=Inches(0.7)

h=doc.add_heading("Phenomenological Tm models — performance",level=1)
sub=doc.add_paragraph(); r=sub.add_run("10-fold cross-validated, n = 235 interfaces. Family A = occupancy-weighted; Family B = elastic / mechanical-coupling. ΔS_Bond and k_rot in cal/mol·K, λ in nm.")
r.italic=True; r.font.size=Pt(9); r.font.color.rgb=RGBColor(0x55,0x55,0x55)

tbl=doc.add_table(rows=1,cols=len(headers)); tbl.style="Table Grid"; tbl.autofit=False
for j,hname in enumerate(headers):
    c=tbl.rows[0].cells[j]; shade(c,"4C78A8"); setcell(c,hname,widths[j],bold=True,white=True,center=(j not in (0,5)))
for row in rows:
    cells=tbl.add_row().cells; best=row[0].startswith("A2")
    for j,val in enumerate(row):
        if best: shade(cells[j],"DCE8F5")
        setcell(cells[j],val,widths[j],bold=(best and j==0),center=(j not in (0,5)))

doc.add_paragraph()
foot=doc.add_paragraph(); r=foot.add_run("Best model: A2 (occupancy-weighted, N-scaled per-bond penalty). The elastic Green's-function (rotational) term is degenerate once g is present (A2→A7, no gain); the physically-derived WLC tethering (B3/B4) performs worst. B3/B4 use the corrected asymmetric 6nt bond (8 nt / 7 nt strands).")
r.font.size=Pt(8); r.font.color.rgb=RGBColor(0x55,0x55,0x55)

doc.save("Supplementary/Model_Comparison_All/phenomenological_models_table.docx")
print("wrote phenomenological_models_table.docx")
