"""
Table S4 (updated 2026-11-12): fitted parameters and uncertainties.
  Multivalent model (Figure 4; Figure4_MultivalentModel/model_fit/params.csv): point estimate, Gauss-Newton SE,
      construct bootstrap SE and 95% percentile CI (200 resamples), n = 236.
  Two-state model A2 (SI; Supplementary/TableS1_Model_Comparison_All/phenomenological_parameter_CIs.csv):
      point estimate, jackknife SE, bootstrap SE and 95% CI (1000 resamples).
Writes TableS4_parameters.csv and TableS4_parameters.docx here.
"""
import os, numpy as np, pandas as pd
from docx import Document
from docx.shared import Pt, RGBColor, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(os.path.dirname(HERE))
M = pd.read_csv(os.path.join(REPO, 'Figure4_MultivalentModel/model_fit/params.csv')).set_index('parameter')
A = pd.read_csv(os.path.join(REPO, 'Supplementary/TableS1_Model_Comparison_All/phenomenological_parameter_CIs.csv')).set_index('parameter')
rows = []
lc = M.loc['log10_Cref_M']
rows.append(dict(model='Multivalent (Fig. 4)', parameter='C_ref', units='mM', estimate=10 ** lc.point_estimate * 1e3,
                 SE_a=f"log10 SE {lc.GN_SE:.3f} (GN)", SE_b=f"log10 SE {lc.boot_SE:.3f}",
                 CI=f"{10 ** lc.boot_CI95_low * 1e3:.1f} – {10 ** lc.boot_CI95_high * 1e3:.1f}",
                 meaning='effective concentration of a bond at zero distance in a four-site interface'))
for key, nm, un, mean in (('lambda_nm', 'λ', 'nm', 'decay length of the effective concentration with mean distance'),
                          ('kappa_per_site', 'κ', 'per site', 'fractional loss of effective concentration per occupied site')):
    r = M.loc[key]
    rows.append(dict(model='Multivalent (Fig. 4)', parameter=nm, units=un, estimate=r.point_estimate, SE_a=f"{r.GN_SE:.3g} (GN)",
                     SE_b=f"{r.boot_SE:.3g}", CI=f"{r.boot_CI95_low:.3g} – {r.boot_CI95_high:.3g}", meaning=mean))
rows.append(dict(model='Multivalent (Fig. 4)', parameter='w', units='', estimate=0.57, SE_a='fixed', SE_b='fixed', CI='—',
                 meaning='FRET readout weight, set from transition widths (Fig. 2, SI)'))
for key, nm, un, mean in (('lambda_nm', 'λ', 'nm', 'decay length of the two-state weight g = e^(−d/λ)'),
                          ('dS_Bond_cal_mol_K', 'ΔS_Bond', 'cal mol⁻¹ K⁻¹', 'per-bond entropy (≈ R ln C_eff)'),
                          ('dS_0_cal_mol_K', 'ΔS_0', 'cal mol⁻¹ K⁻¹', 'dimerization entropy offset')):
    r = A.loc[key]
    rows.append(dict(model='Two-state A2 (SI)', parameter=nm, units=un, estimate=r.point_estimate, SE_a=f"{r.jackknife_SE:.3g} (jackknife)",
                     SE_b=f"{r.bootstrap_SE:.3g}", CI=f"{r.boot_CI95_low:.3g} – {r.boot_CI95_high:.3g}", meaning=mean))
T = pd.DataFrame(rows); T.to_csv(os.path.join(HERE, 'TableS4_parameters.csv'), index=False)

def shade(cell, hexcolor):
    tcPr = cell._tc.get_or_add_tcPr(); sh = OxmlElement('w:shd')
    sh.set(qn('w:val'), 'clear'); sh.set(qn('w:color'), 'auto'); sh.set(qn('w:fill'), hexcolor); tcPr.append(sh)
def setcell(cell, text, w, bold=False, white=False, center=False):
    cell.width = w; p = cell.paragraphs[0]; p.alignment = WD_ALIGN_PARAGRAPH.CENTER if center else WD_ALIGN_PARAGRAPH.LEFT
    r = p.add_run(text); r.font.size = Pt(8.5); r.font.bold = bold
    if white: r.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)
doc = Document(); sec = doc.sections[0]
for m in ('left_margin', 'right_margin'): setattr(sec, m, Inches(0.6))
doc.add_heading('Table S4. Fitted parameters and uncertainties', level=1)
H = ['Model', 'Parameter', 'Units', 'Estimate', 'SE (analytic)', 'SE (bootstrap)', 'Bootstrap 95% CI', 'Meaning']
W = [Inches(.95), Inches(.6), Inches(.7), Inches(.6), Inches(.85), Inches(.7), Inches(.85), Inches(1.75)]
tb = doc.add_table(rows=1, cols=len(H)); tb.style = 'Table Grid'; tb.autofit = False
for j, h in enumerate(H): c = tb.rows[0].cells[j]; shade(c, '3B6AA0'); setcell(c, h, W[j], bold=True, white=True, center=j in (3, 4, 5, 6))
for _, r in T.iterrows():
    c = tb.add_row().cells
    est = f"{r.estimate:.3g}" if r.parameter != 'C_ref' else f"{r.estimate:.1f}"
    for j, v in enumerate([r.model, r.parameter, r.units, est, r.SE_a, r.SE_b, r.CI, r.meaning]):
        if r.model.startswith('Multivalent'): shade(c[j], 'D7F0E5')
        setcell(c[j], str(v), W[j], center=j in (3, 4, 5, 6))
doc.add_paragraph()
p = doc.add_paragraph(); q = p.add_run(
    f"Multivalent model: n = {int(M.n_constructs.iloc[0])} interfaces; Gauss–Newton (GN) SE from the least-squares Jacobian; bootstrap "
    f"over constructs ({int(M.n_bootstrap.iloc[0])} resamples, percentile CI). C_ref was fitted as log10 C_ref (M); its SEs are in log10 units. "
    "Two-state model A2 (former Fig. 4, now SI): jackknife and bootstrap (1000 resamples) over the same interfaces.")
q.font.size = Pt(8); q.font.color.rgb = RGBColor(0x52, 0x51, 0x4E)
_z = doc.settings.element.find(qn('w:zoom'))          # python-docx template omits the required w:percent
if _z is not None and _z.get(qn('w:percent')) is None: _z.set(qn('w:percent'), '100')
doc.save(os.path.join(HERE, 'TableS4_parameters.docx'))
print(T[['model', 'parameter', 'estimate', 'CI']].to_string(index=False))
