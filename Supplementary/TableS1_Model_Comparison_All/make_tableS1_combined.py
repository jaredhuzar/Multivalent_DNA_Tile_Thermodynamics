"""
Table S1 (updated 2026-11-12): every tested Tm model on the canonical 236 interfaces, 10-fold CV with identical folds.
  Family M  multivalent binding models (Figure 4): Figure4_MultivalentModel/model_fit/cv_comparison.csv and
            Supplementary/BondLength_Ceff/fits.csv
  Family A  two-state, distance-weighted ("occupancy") models     } model_comparison_all_235.csv (unified_models.py)
  Family B  two-state, elastic / mechanical-coupling / tethering   }
cv R^2 = 1 - MSE / var(Tm) for every row.  Writes TableS1_all_models.csv and TableS1_all_models.docx here.
Run after unified_models.py, Figure4_MultivalentModel/fit_model.py and BondLength_Ceff/make_bondlength_ceff.py.
"""
import os, re, numpy as np, pandas as pd
from docx import Document
from docx.shared import Pt, RGBColor, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.section import WD_ORIENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(os.path.dirname(HERE))
Y = pd.read_csv(os.path.join(REPO, 'Figure4_MultivalentModel/model_fit/oof_predictions.csv')).Tm_meas_C.values
VAR = np.var(Y)
r2 = lambda rms: 1 - rms ** 2 / VAR

rows = []
CV = pd.read_csv(os.path.join(REPO, 'Figure4_MultivalentModel/model_fit/cv_comparison.csv'))
FORM = {1: 'C_eff,i = C_ref', 2: 'C_eff,i = C_ref e^(−d_i/λ)', 3: 'C_eff,i = C_ref e^(−d_i/λ) e^(−κ(N−4))'}
NAME = {1: 'one uniform C_eff', 2: '+ distance', 3: '+ valency (Fig. 4 model)'}
for _, r in CV[CV.readout.isin(['intermediate', 'docking'])].iterrows():
    pr = [f"C_ref = {10 ** r.log10_Cref_M * 1e3:.1f} mM"]
    if r.n_params > 1: pr.append(f"λ = {r.lambda_nm:.1f} nm")
    if r.n_params > 2: pr.append(f"κ = {r.kappa_per_site:.3f}")
    tag = f"M{int(r.n_params)}" + ('' if r.readout == 'intermediate' else 'd')
    rows.append(dict(id=tag, family='M multivalent', model=f"{NAME[int(r.n_params)]}" + (' [docking readout]' if r.readout == 'docking' else ''),
                     form=FORM[int(r.n_params)], n_params=int(r.n_params), cvRMS=r.cvRMS, cvR2=r2(r.cvRMS), cvrho=r.cvrho,
                     params=', '.join(pr)))
BL = pd.read_csv(os.path.join(REPO, 'Supplementary/BondLength_Ceff/fits.csv'))
for _, r in BL.iterrows():
    kap = 'with' in r.model
    pr = [f"C_ref = {10 ** r.log10_Cref_M * 1e3:.1f} mM", f"λ = {r.lambda_nm:.1f} nm"] + ([f"κ = {r.kappa_per_site:.3f}"] if kap else []) + \
         [f"φ5 = {r.phi5:.2f}", f"φ6 = {r.phi6:.2f}"]
    rows.append(dict(id='M5' if kap else 'M4', family='M multivalent', model='+ bond-type C_eff' + (' + valency' if kap else ' (no valency)'),
                     form='C_eff,i = C_ref φ_type e^(−d_i/λ)' + (' e^(−κ(N−4))' if kap else ''), n_params=int(r.n_params),
                     cvRMS=r.cvRMS, cvR2=r2(r.cvRMS), cvrho=r.cvrho, params=', '.join(pr)))
T = pd.read_csv(os.path.join(HERE, 'model_comparison_all_235.csv'))
PN = {"A1": ["λ", "ΔS_Bond", "base"], "A2": ["λ", "ΔS_Bond", "base"], "A3": ["λ", "ΔS_Bond", "base"], "A4": ["λ", "ΔS_Bond", "base"],
      "A5": ["λ", "ΔS_short", "ΔS_long", "base"], "A6": ["λ_H", "λ_S", "ΔS_Bond", "base"], "A7": ["λ", "ΔS_Bond", "base", "k_rot"],
      "B1": ["ΔS_Bond", "base"], "B2": ["ΔS_Bond", "base", "k_rot"], "B3": ["λ", "base"], "B4": ["λ", "base"]}
FORMS = {"A1": "ΔS = Σg(ΔS_NN+ΔS_Bond)+base", "A2": "ΔS = Σg ΔS_NN + NΔS_Bond + base", "A3": "g on ΔH only",
         "A4": "A1 with nearest-bond g", "A5": "A1 with ΔS_Bond per bond length", "A6": "separate λ for ΔH and ΔS",
         "A7": "A2 − k_rot R ln(1+Σ(z−z̄)²)", "B1": "g = 1: ΣΔH / (ΣΔS_NN + NΔS_Bond + base)", "B2": "B1 − k_rot R ln(1+Σ(z−z̄)²)",
         "B3": "Σg(ΔS_NN + ΔS_teth,WLC) + base", "B4": "Σ(ΔS_NN + ΔS_teth,WLC) + base"}
for _, r in T.iterrows():
    k = r.model.split()[0]; vals = [float(x) for x in re.findall(r"-?\d+\.?\d*(?:[eE][-+]?\d+)?", str(r.params))]
    pr = ', '.join((f"{nm}→∞ (unidentified)" if nm.startswith('λ') and v >= 500 else f"{nm} = {v:g}") for nm, v in zip(PN[k], vals))
    rows.append(dict(id=k, family='A two-state, distance-weighted' if k[0] == 'A' else 'B two-state, mechanical / tethering',
                     model=' '.join(r.model.split()[1:]), form=FORMS[k], n_params=len(PN[k]), cvRMS=r.cvRMS, cvR2=r2(r.cvRMS),
                     cvrho=r.cvrho, params=pr))
ORDER = ['M1', 'M2', 'M3', 'M4', 'M5', 'M1d', 'M2d', 'M3d']
S = pd.DataFrame(rows); S['_o'] = [ORDER.index(i) if i in ORDER else 100 + k for k, i in enumerate(S.id)]
S = S.sort_values('_o').drop(columns='_o').reset_index(drop=True); S.round(4).to_csv(os.path.join(HERE, 'TableS1_all_models.csv'), index=False)

def shade(cell, hexcolor):
    tcPr = cell._tc.get_or_add_tcPr(); sh = OxmlElement('w:shd')
    sh.set(qn('w:val'), 'clear'); sh.set(qn('w:color'), 'auto'); sh.set(qn('w:fill'), hexcolor); tcPr.append(sh)
def setcell(cell, text, w, bold=False, white=False, center=False):
    cell.width = w; p = cell.paragraphs[0]; p.alignment = WD_ALIGN_PARAGRAPH.CENTER if center else WD_ALIGN_PARAGRAPH.LEFT
    r = p.add_run(text); r.font.size = Pt(8); r.font.bold = bold
    if white: r.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)
doc = Document(); sec = doc.sections[0]; sec.orientation = WD_ORIENT.LANDSCAPE
sec.page_width, sec.page_height = Inches(11), Inches(8.5)
for m in ('left_margin', 'right_margin'): setattr(sec, m, Inches(0.5))
sec.top_margin = Inches(0.6)
doc.add_heading('Table S1. Tested models of interface Tm', level=1)
H = ['ID', 'Model', 'Form', 'k', 'CV RMS (°C)', 'CV R²', 'CV ρ', 'Fitted parameters (all data)']
W = [Inches(.45), Inches(2.1), Inches(2.6), Inches(.3), Inches(.75), Inches(.6), Inches(.6), Inches(2.6)]
tb = doc.add_table(rows=1, cols=len(H)); tb.style = 'Table Grid'; tb.autofit = False
for j, h in enumerate(H): c = tb.rows[0].cells[j]; shade(c, '3B6AA0'); setcell(c, h, W[j], bold=True, white=True, center=j in (3, 4, 5, 6))
fam = None
for _, r in S.iterrows():
    if r.family != fam:
        fam = r.family; c = tb.add_row().cells; m = c[0].merge(c[-1]); shade(m, 'E6E5E0')
        setcell(m, {'M': 'M. Multivalent binding models (this work, Fig. 4)', 'A': 'A. Two-state models, distance-weighted bond contributions',
                    'B': 'B. Two-state models, mechanical coupling and tether entropy'}[fam[0]], Inches(10), bold=True)
    c = tb.add_row().cells; main = r.id == 'M3'
    vals = [r.id, r.model, r.form, str(r.n_params), f'{r.cvRMS:.2f}', f'{r.cvR2:.3f}', f'{r.cvrho:.3f}', r.params]
    for j, v in enumerate(vals):
        if main: shade(c[j], 'D7F0E5')
        setcell(c[j], v, W[j], bold=main and j < 2, center=j in (3, 4, 5, 6))
doc.add_paragraph()
p = doc.add_paragraph(); r = p.add_run(
    'All models fitted by least squares to the median Tm of the same 236 interfaces and evaluated by 10-fold cross-validation over '
    'interfaces with identical folds (KFold, shuffle, random_state 0); CV statistics are out-of-fold. Family M: binding polynomial '
    'K_a = C_ref⁻¹[Π(1 + K_i C_eff,i) − 1] with mass balance at 4.5 nM per tile; FRET readout θ = f[w + (1 − w)p_rep], w = 0.57 '
    '(intermediate) or 1 (docking); d_i = mean distance of duplex i to the other duplexes; N = occupied sites; φ_type = effective '
    'concentration of 5- or 6-nt bonds relative to 7-nt bonds. φ5 and κ are strongly correlated (r = 0.92 when fitted together) '
    'because 5-nt bonds dominate the many-site interfaces. Families A and B: two-state Tm = ΣgΔH_NN / (ΔS_eff + R ln(C/2)), '
    'g = e^(−d/λ); ΔS in cal mol⁻¹ K⁻¹, λ in nm. Shaded: model shown in Fig. 4. Two-state models imply melts about half as wide '
    'as measured (SI Appendix, Fig. S10).')
r.font.size = Pt(8); r.font.color.rgb = RGBColor(0x52, 0x51, 0x4E)
_z = doc.settings.element.find(qn('w:zoom'))          # python-docx template omits the required w:percent
if _z is not None and _z.get(qn('w:percent')) is None: _z.set(qn('w:percent'), '100')
doc.save(os.path.join(HERE, 'TableS1_all_models.docx'))
print(S[['id', 'model', 'n_params', 'cvRMS', 'cvrho']].round(3).to_string(index=False))
