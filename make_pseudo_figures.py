import os as _os
def _find_repo():
    d = _os.path.dirname(_os.path.abspath(__file__))
    while d != '/' and not _os.path.exists(_os.path.join(d, 'Filtered_DF_Tm.csv')):
        d = _os.path.dirname(d)
    return d
REPO = _find_repo()
"""
Pseudo-figures: one composite SVG per manuscript figure, tiling that figure's panels into a labelled grid in the
CURRENT paper layout (panel map below), plus every individual panel's SVG and CSV.  A pseudo-figure is a layout
mock-up (every panel/slot at a glance), NOT the final publication figure.  Real data panels are NESTED as vector SVG
(ids namespaced per panel so clip-paths / glyph defs do not collide).  Slots without a standalone SVG are drawn as
labelled placeholders:
  - dashed grey  = schematic / external art supplied separately (e.g. Greg's cartoons)
  - dashed amber = data panel that still needs script work
  - dashed red   = a file that was expected but is missing

EVERY RUN REGENERATES EVERYTHING:
  1. re-runs the scripts that produce the panels (fast generators; see GENERATORS) so the panels reflect the current
     code and data -- skip with --no-regen (e.g. right after run_all.sh, which already ran them);
     --refit additionally re-fits the Figure 4 model (fit_model.py, ~3 min incl. bootstrap);
     --notebooks additionally re-executes the notebooks behind Fig 3B/3E (Tile_Tile_Master.ipynb, slow).
  2. deletes and rebuilds Pseudo_Figures/FigureN/ (no stale files survive): for each panel <FigN><letter>.svg and
     <FigN><letter>.csv (the data behind it; extra tables get a suffix), and MANIFEST.txt (slot -> sources).
  3. writes the composite Pseudo_Figures/FigureN_pseudo.svg from those copies.

Usage:  python3 make_pseudo_figures.py [--no-regen] [--refit] [--notebooks] [Figure2 Figure4 ...]
All panel sources use the updated (2026-11-12/13) analyses: Tm-anchored Figure 2 with multivalent-model avidity and
DDG_coop, Tm-based Fig 3E, the multivalent Figure 4, and Figure 6 with the multivalent comparator.
"""
import os, re, sys, math, shutil, subprocess, textwrap, time
os.chdir(REPO)
OUT = "Pseudo_Figures"; os.makedirs(OUT, exist_ok=True)

SCHEM = ("schematic", None)                       # external art (Greg / AFM); supplied separately
def TODO(note): return ("todo", note)             # data panel that still needs script work

# ================= generators (run in this order before assembling) =================
F2 = "Figure2_Updated"; F3 = "Figure3_BondClustering"; F4 = "Figure4_MultivalentModel"
F5 = "Figure5_MD_Simulations/Figure_AH_SelfContained"; F6 = "Figure6_MachineLearning"
GENERATORS = {   # figure -> scripts (python) that rebuild its panels; order matters (dependencies first)
 "Figure1": ["Figure1_Platform_AFM/make_melting_curves.py"],
 "Figure2": [f"{F2}/make_uniform_ceff_panel.py", f"{F2}/make_updated_figure2.py"],
 "Figure3": [f"{F3}/Fig3b_Tm_vs_BondSpan_OLSResiduals/make_fig3b.py", f"{F3}/Fig3e_RelativeCeff_vs_Distance/make_fig3e.py"],
 "Figure4": [f"{F4}/make_figure4.py"],
 "Figure5": [f"{F5}/regen_fig5_tm_panels.py"],
 "Figure6": [f"{F6}/make_figure6.py", "Supplementary/S6_SHAP_Values/recompute_shap_fulltrain.py"],
 "SI": [f"{F2}/make_uniform_ceff_panel.py", f"{F2}/make_updated_figure2.py", f"{F2}/make_added_bond_occupancy.py",
        "Supplementary/Width_vs_MeanGap_7nt_4bonds/make_width_vs_meangap.py", "Supplementary/Ceff_Geometry_OLS/make_ceff_ols.py",
        "Supplementary/BondLength_Ceff/make_bondlength_ceff.py", "Supplementary/BondLength_Ceff/make_tether_penalty.py", "Supplementary/BondLength_Ceff/make_kappa_confound.py", "Supplementary/NineNt_MeltingCurves/make_9nt_melting_curves.py",
        "Supplementary/PredictedVsMeasured_NN_Multivalent/make_pred_vs_meas.py", "Supplementary/TransitionWidth_ByBonds/make_width_by_bonds.py",
        "Supplementary/add_ols_and_meangap.py", f"{F5}/regen_fig5_tm_panels.py", "Supplementary/S6_SHAP_Values/recompute_shap_fulltrain.py",
        "Supplementary/oxDNA/make_rmsf_panel.py", "Supplementary/SI_PhenomenologicalModel_A2/make_twostate_pred_rank.py"],
}
REFIT = [f"{F4}/fit_model.py"]                    # --refit: run before Figure 4 / Figure 6 generators
NOTEBOOKS = ["Tile_Tile_Master.ipynb"]            # --notebooks: Fig 3B (mean gap) and Fig 3E (positional series)

# ================= current paper panel map: (label, svg, csv | [csv, ...]) or (label, placeholder) =================
FIGS = {
 "Figure1": ("Figure 1 — Platform & melting readout", 2, [
    ("A–B  platform schematic", SCHEM),
    ("C  fraction melted — melting curves", "Figure1_Platform_AFM/Melting Curves.svg", "Figure1_Platform_AFM/Melting_Curves_data.csv")]),

 "Figure2": ("Figure 2 — Multivalency, cooperative unit, avidity", 3, [
    ("A  interface + FRET schematic", SCHEM),
    ("B  melts far above its best single bond", f"{F2}/PanelB_StrongestBond/figure.svg", f"{F2}/PanelB_StrongestBond/data.csv"),
    ("C  size of the cooperative unit", f"{F2}/PanelD_CooperativeUnit/figure.svg", f"{F2}/PanelD_CooperativeUnit/data.csv"),
    ("D  melting temperature vs sites", f"{F2}/PanelC_Tm_vs_Sites/figure.svg", f"{F2}/PanelC_Tm_vs_Sites/data.csv"),
    ("E  avidity (multivalent-model effective concentration)", f"{F2}/PanelF_Avidity/figure.svg", f"{F2}/PanelF_Avidity/data.csv"),
    ("F  added-bond ΔΔG_coop (docking readout)", f"{F2}/PanelF_DDGcoop_Docking/figure.svg", f"{F2}/PanelF_DDGcoop_Docking/data.csv"),
    ("G  mechanism schematic", SCHEM)]),

 "Figure3": ("Figure 3 — Bond clustering / spatial arrangement", 3, [
    ("A  interface schematic (Greg)", SCHEM),
    ("B  Tm vs mean gap (7nt, 4 bonds)", f"{F3}/Fig3a_Tm_vs_MeanGap_7nt_4bonds/figure.svg", f"{F3}/Fig3a_Tm_vs_MeanGap_7nt_4bonds/data.csv"),
    ("C  residual Tm vs bond span", f"{F3}/Fig3b_Tm_vs_BondSpan_OLSResiduals/figure.svg", f"{F3}/Fig3b_Tm_vs_BondSpan_OLSResiduals/data.csv"),
    ("D  fixed bonds, varying position schematic (Greg)", SCHEM),
    ("E  moving same bond away from fixed pair (melting result)", f"{F3}/Fig3c_Tm_vs_BondPosition_ConstantSeq/figure.svg", f"{F3}/Fig3c_Tm_vs_BondPosition_ConstantSeq/data.csv"),
    ("F  C_eff decays with spacing", f"{F3}/Fig3e_RelativeCeff_vs_Distance/figure.svg", f"{F3}/Fig3e_RelativeCeff_vs_Distance/data.csv")]),

 "Figure4": ("Figure 4 — Multivalent model", 3, [
    ("A  model schematic (Greg)", SCHEM),
    ("B  out-of-fold predictions (new model)", f"{F4}/panels/B.svg", f"{F4}/model_fit/oof_predictions.csv"),
    ("C  nested model comparison", f"{F4}/panels/C.svg", f"{F4}/theory/NestedModelComparison.csv"),
    ("D  theoretical Tm vs binding sites", f"{F4}/panels/D.svg", [f"{F4}/theory/Tm_vs_NumSites.csv", f"{F4}/model_fit/oof_predictions.csv"]),
    ("E  theory: identical bonds, spacing lowers Tm", f"{F4}/panels/F.svg", [f"{F4}/theory/MeltCurves_vs_Spacing_summary.csv", f"{F4}/theory/MeltCurves_vs_Spacing.csv"]),
    ("F  design map: predicted-Tm heatmap", f"{F4}/panels/H.svg", f"{F4}/theory/DesignMap_Tm.csv")]),

 "Figure5": ("Figure 5 — Coarse-grained MD", 3, [
    ("A  dimer schematic (Greg)", SCHEM),
    ("B  decay with spacing (rigid vs flexible centre)", f"{F5}/plots_robust_v2/b/panel_b.svg", f"{F5}/plots_robust_v2/b/panel_b.csv"),
    ("C  Tm vs linker bending stiffness (Sp 2 & 6)", f"{F5}/plots_robust_v2/d/panel_d.svg", f"{F5}/plots_robust_v2/d/panel_d.csv"),
    ("D  Tm vs rod bending stiffness (Sp 2 & 6)", f"{F5}/plots_robust_v2/c/panel_c.svg", f"{F5}/plots_robust_v2/c/panel_c.csv"),
    ("E  linker stiffness × unbound tip gap", f"{F5}/plots_robust/mechanism_panels/e/panel_e.svg", [f"{F5}/plots_robust/mechanism_panels/e/panel_e.csv", f"{F5}/plots_robust/mechanism_panels/e/panel_e_summary.csv"]),
    ("F  anchor spacing × scaffold bending (3 T)", f"{F5}/plots_robust/mechanism_panels/f/panel_f.svg", f"{F5}/plots_robust/mechanism_panels/f/panel_f.csv"),
    ("G  mean tip-gap heatmap (linker × scaffold)", f"{F5}/plots_robust/mechanism_panels/g/panel_g.svg", f"{F5}/plots_robust/mechanism_panels/g/panel_g.csv"),
    ("H  tip-gap SD heatmap (linker × scaffold)", f"{F5}/plots_robust/mechanism_panels/h/panel_h.svg", f"{F5}/plots_robust/mechanism_panels/h/panel_h.csv")]),

 "Figure6": ("Figure 6 — Machine learning", 3, [
    ("A  model schematic (Greg)", SCHEM),
    ("B  XGBoost vs experimental Tm", f"{F6}/Fig6b_XGBoost_StabilityScatter/figure.svg", f"{F6}/Fig6b_XGBoost_StabilityScatter/data.csv"),
    ("C  Spearman ρ OOF (drop Additive-NN series)", f"{F6}/Fig6d_ModelComparison_SpearmanRho/figure.svg", f"{F6}/Fig6d_ModelComparison_SpearmanRho/data.csv"),
    ("D  SHAP values", "Supplementary/S6_SHAP_Values/figure.svg", "Supplementary/S6_SHAP_Values/data.csv"),
    ("E  Tm vs edge-site occupied (raw + mean-gap OLS; 7nt 4-bond)", f"{F6}/Fig6e_EdgeSites_Spacing/figure.svg", f"{F6}/Fig6e_EdgeSites_Spacing/data.csv")]),
}
SW = "Supplementary/Width_vs_MeanGap_7nt_4bonds"; SA2 = "Supplementary/SI_PhenomenologicalModel_A2"
F5M = f"{F5}/plots_robust/mechanism_panels"; F5V = f"{F5}/plots_robust_v2"; SH = "Supplementary/S6_SHAP_Values"
# SI figures: (SI number or None, file name, description, svg, csv | [csv, ...]).  Numbered -> Pseudo_Figures/SI/FigS<nn><x>_<name>.*,
# None -> Pseudo_Figures/SI/Extra/.  README.md in Pseudo_Figures/SI maps each number to its files.
SI_FIGS = [
 ("4", "NineNt_MeltingCurves", "9-nt bonds: raw and fitted melting curves of the 9-nt interfaces that passed QC (edge code and Tm per panel)",
  "Supplementary/NineNt_MeltingCurves/figure.svg",
  ["Supplementary/NineNt_MeltingCurves/summary.csv", "Supplementary/NineNt_MeltingCurves/fits.csv", "Supplementary/NineNt_MeltingCurves/raw_points.csv"]),
 ("6", "TransitionWidth_ByBonds", "Transition width (raw, °C) of 5- and 7-nt interfaces at different bond counts",
  "Supplementary/TransitionWidth_ByBonds/figure.svg", "Supplementary/TransitionWidth_ByBonds/data.csv"),
 ("7", "UniformCeff_Width", "Example melting curve with uniform C_eff: all-or-none, reporter readout and docking readout (PanelH_UniformCeff_Width)",
  f"{F2}/PanelH_UniformCeff_Width/figure.svg", [f"{F2}/PanelH_UniformCeff_Width/data.csv", f"{F2}/PanelH_UniformCeff_Width/example_curves.csv"]),
 ("8", "CooperativeFraction", "Cooperative fraction ΔH_vH/ΔH_NN (PanelC_CoopFraction)", f"{F2}/PanelC_CoopFraction/figure.svg", f"{F2}/PanelC_CoopFraction/data.csv"),
 ("9", "AddedBondOccupancy", "Added duplexes at C_eff: T* relative to the interface Tm, fraction formed, and predicted vs measured ΔTm (4→5 7-nt vs 9→10 5-nt)",
  f"{F2}/PanelSI_AddedBondOccupancy/figure.svg", [f"{F2}/PanelSI_AddedBondOccupancy/duplexes.csv", f"{F2}/PanelSI_AddedBondOccupancy/pairs.csv"]),
 ("10", "OLS_Tm_Geometry", "OLS of melting temperature on number of bonds, bond span, bond length and GC fraction",
  "Supplementary/S3_OLS_Tm_geometry/figure.svg",
  ["Supplementary/S3_OLS_Tm_geometry/data.csv", "Supplementary/S3_OLS_Tm_geometry/ols_coefficients.csv", "Supplementary/S3_OLS_Tm_geometry/partial_correlations.csv"]),
 ("11", "vantHoff_and_Width_vs_Gap", "van't Hoff enthalpy (A) and transition width (B) vs separation, four-site 7-nt interfaces (merged)",
  f"{SW}/panels/dHvH_and_width_vs_gap.svg", f"{SW}/panels/dHvH_and_width_vs_gap.csv"),
 ("11a", "vantHoffEnthalpy_vs_Gap", "panel A alone: van't Hoff enthalpy change with separation", f"{SW}/panels/dHvH_vs_gap.svg", f"{SW}/panels/dHvH_vs_gap.csv"),
 ("11b", "Width_vs_Gap", "panel B alone: transition width with separation", f"{SW}/panels/width_vs_gap.svg", f"{SW}/panels/width_vs_gap.csv"),
 ("12", "ReadoutWeight_BestFit", "Best-fit FRET readout weight w per construct (PanelI_ReadoutWeight, best-fit panel only)",
  f"{F2}/PanelI_ReadoutWeight/figure_bestfit_w.svg", f"{F2}/PanelI_ReadoutWeight/data.csv"),
 ("13", "TwoState_PredAndRank", "Two-state model with distance weight on enthalpy and entropy, out-of-fold: (A) predicted vs measured Tm (RMSD 1.30 °C, ρ 0.933), (B) rank vs rank",
  f"{SA2}/Phenomenological_Model/pred_and_rank.svg",
  [f"{SA2}/Phenomenological_Model/pred_and_rank.csv", f"{SA2}/Phenomenological_Model/pred_vs_exp.csv", f"{SA2}/Phenomenological_Model/rank_rank.csv"]),
 ("14", "TwoState_Width", "Two-state model implies transitions about half as wide as measured and broadening with spacing",
  f"{SW}/panels/twostate_width_vs_measured.svg", f"{SW}/panels/twostate_width_vs_measured.csv"),
 ("15", "KappaConfound", "κ cannot be separated from a 5-nt bond-length effect: (A) equal CV error, (B) κ–φ5 trade-off at constant error",
  "Supplementary/BondLength_Ceff/kappa_confound.svg",
  ["Supplementary/BondLength_Ceff/kappa_confound_cv.csv", "Supplementary/BondLength_Ceff/kappa_confound_profile.csv"]),
 ("16a", "MD_Tm_vs_BondStrength", "Simulated Tm versus bond strength (Morse depth)", f"{F5V}/e/panel_e.svg", [f"{F5V}/e/panel_e.csv", f"{F5V}/e/panel_e_replicates.csv"]),
 ("16b", "MD_Tm_vs_BondsPerInterface", "Simulated Tm vs bonds per interface, flexible and rigid scaffold", f"{F5V}/f/panel_f.svg", [f"{F5V}/f/panel_f.csv", f"{F5V}/f/panel_f_replicates.csv"]),
 ("16c", "MD_MeltingExample", "Simulated melting example: mean bonds formed vs temperature (3 bonds, spacing 3)",
  f"{F5M}/melting_example/panel_melting_example.svg", f"{F5M}/melting_example/panel_melting_example.csv"),
 ("17a", "SHAP_Beeswarm", "SHAP values for XGBoost (all 236 interfaces)", f"{SH}/figure.svg", f"{SH}/data.csv"),
 ("17b", "SHAP_MeanBySite", "Mean SHAP value at each site, occupied vs empty", f"{SH}/mean_shap_by_site.svg", f"{SH}/mean_shap_by_site.csv"),
 ("18", "RMSF_Tile", "RMSF of the tile (oxDNA monomer)", "Supplementary/oxDNA/rmsf_tile.svg", "Supplementary/oxDNA/rmsf_tile.csv"),
 (None, "WidthShape_vs_Gap", "melt shape vs mean gap with model widths (former S10 composite)", f"{SW}/figure.svg", f"{SW}/data.csv"),
 (None, "TetherPenalty", "per-bond RT ln(1 M/C_eff) by bond length vs WLC", "Supplementary/BondLength_Ceff/tether_penalty.svg", "Supplementary/BondLength_Ceff/tether_penalty.csv"),
 (None, "BondLength_Ceff", "bond-type effective concentration vs WLC", "Supplementary/BondLength_Ceff/figure.svg", "Supplementary/BondLength_Ceff/fits.csv"),
 (None, "Ceff_Geometry", "per-interface C_eff vs geometry (OLS)", "Supplementary/Ceff_Geometry_OLS/figure.svg",
  ["Supplementary/Ceff_Geometry_OLS/data.csv", "Supplementary/Ceff_Geometry_OLS/ols_coefficients.csv", "Supplementary/Ceff_Geometry_OLS/partial_correlations.csv"]),
 (None, "PredVsMeasured", "additive NN (1 M) vs multivalent predicted Tm", "Supplementary/PredictedVsMeasured_NN_Multivalent/figure.svg", "Supplementary/PredictedVsMeasured_NN_Multivalent/data.csv"),
 (None, "TwoState_Theory", "two-state model A2 theory curves (former Fig 4C-F)", f"{SA2}/theory_fig4_cdef.svg", [f"{SA2}/theory_vs_Nbonds.csv", f"{SA2}/theory_vs_spacing.csv"]),
]
FIGS_SI = {"SI": ("Supplementary figures", 3, None)}
FIGS.update(FIGS_SI)

# ================= SVG composition =================
def load_panel(path, idx):
    txt = open(path, encoding='utf-8', errors='replace').read()
    m = re.search(r'<svg\b[^>]*>', txt, re.S)
    root = m.group(0)
    vb = re.search(r'viewBox="([^"]+)"', root)
    if vb:
        viewbox = vb.group(1)
    else:
        w = re.search(r'width="([\d.]+)', root); h = re.search(r'height="([\d.]+)', root)
        viewbox = f"0 0 {float(w.group(1))} {float(h.group(1))}"
    inner = txt[m.end(): txt.rstrip().rfind('</svg>')]
    pre = f"q{idx}_"
    inner = re.sub(r'id="([^"]+)"', lambda mm: f'id="{pre}{mm.group(1)}"', inner)
    inner = re.sub(r'url\(#([^)]+)\)', lambda mm: f'url(#{pre}{mm.group(1)})', inner)
    inner = re.sub(r'\b(xlink:href|href)="#([^"]+)"', lambda mm: f'{mm.group(1)}="#{pre}{mm.group(2)}"', inner)
    return viewbox, inner

def esc(s): return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')

PH_COLOR = {'schematic': '#9a9a9a', 'todo': '#c98a00', 'missing': '#c0392b'}
PH_KIND  = {'schematic': 'schematic — supplied separately', 'todo': 'needs script work', 'missing': 'file missing'}

def placeholder(x, y, w, h, kind, note):
    col = PH_COLOR[kind]
    P = [f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="#fbfbfa" stroke="{col}" stroke-width="1.1" '
         f'stroke-dasharray="6,4"/>',
         f'<text x="{x+w/2}" y="{y+h/2-18}" font-size="11" font-style="italic" text-anchor="middle" '
         f'font-family="sans-serif" fill="{col}">[{PH_KIND[kind]}]</text>']
    for i, ln in enumerate(textwrap.wrap(note, 34)[:4]):
        P.append(f'<text x="{x+w/2}" y="{y+h/2+i*14}" font-size="10.5" text-anchor="middle" '
                 f'font-family="sans-serif" fill="#555">{esc(ln)}</text>')
    return '\n'.join(P)

def compose(fig_name, title, slots, cols=3, cellw=360, cellh=270, out=None):
    """slots: list of (label, src) where src is a placeholder tuple or an SVG path (already copied)."""
    if not slots:
        print(f"  {fig_name}: no panels, skipped"); return None
    n = len(slots); cols = min(cols, n); rows = math.ceil(n / cols)
    pad, labh, titleh = 14, 22, 44
    W = cols * cellw + (cols + 1) * pad
    H = titleh + rows * (cellh + labh) + (rows + 1) * pad
    P = [f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
         f'width="{W}pt" height="{H}pt" viewBox="0 0 {W} {H}">',
         f'<rect width="{W}" height="{H}" fill="white"/>',
         f'<text x="{W/2}" y="28" font-size="19" font-weight="bold" text-anchor="middle" '
         f'font-family="sans-serif">{esc(title)}</text>']
    for i, (label, src) in enumerate(slots):
        r, c = divmod(i, cols)
        x = pad + c * (cellw + pad)
        y = titleh + pad + r * (cellh + labh + pad)
        P.append(f'<text x="{x+3}" y="{y+15}" font-size="11.5" font-weight="bold" font-family="sans-serif" '
                 f'fill="#222">{esc(label)}</text>')
        if isinstance(src, tuple):
            P.append(placeholder(x, y + labh, cellw, cellh, src[0], src[1] or ''))
        elif os.path.exists(src):
            viewbox, inner = load_panel(src, i)
            P.append(f'<rect x="{x}" y="{y+labh}" width="{cellw}" height="{cellh}" fill="none" stroke="#ddd" '
                     f'stroke-width="0.7"/>')
            P.append(f'<clipPath id="cell{i}"><rect x="{x}" y="{y+labh}" width="{cellw}" height="{cellh}"/></clipPath>'
                     f'<g clip-path="url(#cell{i})">'      # clip panels cropped from a larger figure (all renderers)
                     f'<svg x="{x}" y="{y+labh}" width="{cellw}" height="{cellh}" viewBox="{viewbox}" overflow="hidden" '
                     f'preserveAspectRatio="xMidYMid meet">{inner}</svg></g>')
        else:
            P.append(placeholder(x, y + labh, cellw, cellh, 'missing', src))
    P.append('</svg>')
    out = out or os.path.join(OUT, f"{fig_name}_pseudo.svg")
    open(out, 'w', encoding='utf-8').write('\n'.join(P))
    return out

# ================= regenerate, collect, compose =================
def run(cmd, label):
    t0 = time.time()
    r = subprocess.run(cmd, cwd=REPO, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    ok = r.returncode == 0
    print(f"    {'OK  ' if ok else 'FAIL'} {label}  ({time.time() - t0:.0f} s)")
    if not ok: print(textwrap.indent('\n'.join(r.stdout.strip().splitlines()[-12:]), '        '))
    return ok

def regenerate(figs, refit, notebooks):
    print("Regenerating panel sources:")
    done, failed = set(), []
    if notebooks:
        for nb in NOTEBOOKS:
            if not run(['jupyter', 'execute', nb], nb): failed.append(nb)
    if refit and ({'Figure4', 'Figure6'} & set(figs)):
        for s in REFIT:
            if not run([sys.executable, s], s): failed.append(s)
    for f in figs:
        for s in GENERATORS.get(f, []):
            if s in done: continue
            done.add(s)
            if not run([sys.executable, s], s): failed.append(s)
    return failed

# Which columns of each panel's CSV(s) are drawn in its SVG (written into each Figure<N>/MANIFEST.txt).
# File names are the copies in Pseudo_Figures/Figure<N>/. Columns not listed are carried for reference only.
# Keep in sync with the plotting code named in each entry when a panel changes.
PLOTTED = {
 "Figure1": {
  "C": ["Fig1C.csv: one line per `code` (one interface, ~200 curves): x = `temperature_C`, y = `normalized_signal`",
        "`Tm_C`, `dx` = the per-interface median sigmoid parameters that generate each curve (not drawn separately)",
        "(Figure1_Platform_AFM/make_melting_curves.py)"],
 },
 "Figure2": {
  "B": ["Fig2B.csv, split by `nt` (7-nt left, 5-nt right):",
        "grey swarm = `strongest_bond_tm_C`; coloured swarm = `interface_tm_C`; black bars = median and IQR of each",
        "arrow label (+38 °C, +53 °C) = median of `gain_C` per `nt`",
        "(Figure2_Updated/make_updated_figure2.py, panel_B)"],
  "C": ["Fig2C.csv: x = `n_sites`, y = `coop_unit_helices` (swarm, coloured by `nt`); black bars = median and IQR per `n_sites`",
        "red dashed line = 2N (all helices as one unit), computed from `n_sites`, not a column",
        "text box = range of the per-`n_sites` medians of `coop_unit_helices`",
        "not drawn: `width_dx_C`, `dHvh_width_kcal`, `allornone_dx_C`, `width_ratio_obs_over_aon` (inputs to coop_unit_helices)",
        "(make_updated_figure2.py, panel_D)"],
  "D": ["Fig2D.csv: y = `Tm_C` (swarm), one column per (`nt`, `n_sites`) group: 7-nt 4, 5 sites; 5-nt 9, 10 sites",
        "black bars = median and IQR of `Tm_C` per group; colour = `nt`",
        "(make_updated_figure2.py, panel_C)"],
  "E": ["Fig2E.csv: histogram of `avidity_log10_Ceff_over_Cpartner` (bin width 0.1), one colour per `nt`",
        "dashed line + label = median of `avidity_log10_Ceff_over_Cpartner`; label C_eff = 10^median(`log10_Ceff_M`) in mM",
        "not drawn: `Tm_C` (the input the per-interface C_eff is fitted to), `n_sites`",
        "(make_updated_figure2.py, panel_F; docking readout)"],
  "F": ["Fig2F.csv: x = `regime` (4->5, 9->10), y = `ddG_coop_kcal` (jittered points, one per add-one-site pair)",
        "black bars + label = median and IQR of `ddG_coop_kcal` per `regime`; all rows have `readout` = docking",
        "not drawn: `lo`, `hi`, `log10Ceff_lo`, `log10Ceff_hi`, `Ceff_ratio_hi_lo` (ddG_coop = -RT ln Ceff_ratio_hi_lo)",
        "(make_updated_figure2.py, panel_H_coop; data from make_uniform_ceff_panel.py)"],
 },
 "Figure3": {
  "B": ["Fig3B.csv: x = `MeanGap_nm`, y = `Tm` (one point per four-site 7-nt interface)",
        "red line = least-squares fit of `Tm` on `MeanGap_nm`; shaded band = seaborn regplot 95% bootstrap CI (unseeded, redraws slightly each run)",
        "not drawn: `Tm_SE`",
        "(Tile_Tile_Master.ipynb, Fig3a export cell; sns.regplot)"],
  "C": ["Fig3C.csv: x = `bond_span_residual_nm`, y = `Tm_residual_C` (OLS added-variable residuals)",
        "red dashed line = least-squares fit of y on x; title partial ρ = Spearman(x, y)",
        "(Figure3_BondClustering/Fig3b_Tm_vs_BondSpan_OLSResiduals/make_fig3b.py)"],
  "E": ["Fig3E.csv: x = `Bond_Position`, y = `Tm`, error bars = `Tm_SE`",
        "(Tile_Tile_Master.ipynb, Fig3c export cell)"],
  "F": ["Fig3F.csv, rows with `readout` = docking only (intermediate and reporter rows are not drawn):",
        "x = `nearest_bond_distance_nm` (30 = the >=30 nm bin); filled circles = `rel_Ceff`, error bars `ci_lo`–`ci_hi`",
        "open squares (x + 0.7 nm) = `rel_occupancy_medianTm`, error bars `rel_occ_ci_lo`–`rel_occ_ci_hi`",
        "dashed curve = exp(-(d - 6 nm)/`lambda_exp_nm`); text = `Tm_rms_C` and the uniform-C_eff RMS",
        "not drawn: `w`, `boot_ci_lo`, `boot_ci_hi`, `n_sites_in_bin`, `log10_Cref_M`, `eval_Tm_C`",
        "(Figure3_BondClustering/Fig3e_RelativeCeff_vs_Distance/make_fig3e.py)"],
 },
 "Figure4": {
  "B": ["Fig4B.csv: x = `Tm_pred_oof_C`, y = `Tm_meas_C`, colour = `group`; grey line = identity",
        "RMS and ρ text come from model_fit/cv_comparison.csv (intermediate readout, 3 parameters), not from this file",
        "not drawn: `n_sites`, `Tm_pred_fit_C`, `width_pred_C`, `oof_*` (other models' out-of-fold predictions)",
        "(Figure4_MultivalentModel/make_figure4.py, panel B)"],
  "C": ["Fig4C.csv: one bar per `model` (ordered by `n_params`), bar height = `cvRMS`; bar labels = `cvRMS` and `cvrho`",
        "not drawn: `log10_Cref_M`, `lambda_nm`, `kappa_per_site`",
        "(make_figure4.py, panel C)"],
  "D": ["Fig4D.csv, one colour per `bond` (7-nt, 5-nt): x = `n_sites`; solid line = `Tm_mean_C`;",
        "shaded band = `Tm_min_C`–`Tm_max_C`; dashed line = `Tm_no_valency_mean_C` (κ = 0). Not drawn: `n_windows`",
        "Fig4D_2_oof_predictions.csv, `group` 7-nt and 5-nt, `n_sites` with >= 3 interfaces: filled circle = median `Tm_meas_C`,",
        "error bars = its IQR; open circle (x + 0.18) = median `Tm_pred_fit_C` of the same interfaces",
        "(make_figure4.py, panel D)"],
  "E": ["Fig4E_2_MeltCurves_vs_Spacing.csv: x = `T_C`, y = `theta_gap<g>nm`, one curve per gap (6–20 nm)",
        "Fig4E.csv: marker at (`Tm_C`, 0.5) on each curve; legend = `gap_nm`, `Tm_C`, `width_C`",
        "not drawn: `bond_dH_kcal`, `bond_dS_kcal_K` (the identical-bond parameters used)",
        "(make_figure4.py, panel F; inner title reads F)"],
  "F": ["Fig4F.csv: heatmap, x = `n_sites`, y = `spacing_sites` (tick labels = `spacing_nm`), colour and cell number = `Tm_mean_C`",
        "blank cells = no placement fits the 11-site edge. Not drawn: `n_placements`",
        "(make_figure4.py, panel H; inner title reads H)"],
 },
 "Figure5": {
  "B": ["Fig5B.csv: x = `x` (bond spacing), y = `Tm` (kT), error bars = `Tm_sd` (s.d. over seeds); one line per `series`",
        "(Analysis_Tm_trends.ipynb cell 8, _build_tm_panel 'b')"],
  "C": ["Fig5C.csv: x = `x` (linker bending stiffness, chain_bend), y = `Tm` (kT), error bars = `Tm_sd`; one line per `series` (sp2, sp6)",
        "(Analysis_Tm_trends.ipynb cell 8, _build_tm_panel 'd')"],
  "D": ["Fig5D.csv: x = `x` (scaffold bending stiffness, rod_bend; log axis), y = `Tm` (kT), error bars = `Tm_sd`; one line per `series` (sp2, sp6)",
        "(Analysis_Tm_trends.ipynb cell 8, _build_tm_panel 'c')"],
  "E": ["Fig5E.csv: x = `unbound_tip_gap_beads` (horizontal swarm), y row = `linker_bending_stiffness` (4, 10, 20, 40); all rows `temperature` = 2",
        "Fig5E_2_panel_e_summary.csv: vertical bar = `median`, horizontal bar = `q1`–`q3`",
        "not drawn: Fig5E.csv `pair_index`; summary `n`, `mean`, `sd`",
        "(Figure_AH_SelfContained/make_fig5_mechanism_panels.py, E)"],
  "F": ["Fig5F.csv: x = `rod_bend`, y = `intra_sigma`, error bars = `intra_sigma_reperr`; one line per `temperature` (1.25, 1.5, 1.75)",
        "not drawn: `n_reps`",
        "(make_fig5_mechanism_panels.py, F)"],
  "G": ["Fig5G.csv: heatmap, x = `rod_bend`, y = `chain_bend`, colour and cell number = `tipgap_mean` (`temperature` = 2)",
        "not drawn: `tipgap_n_unbound`",
        "(make_fig5_mechanism_panels.py, G)"],
  "H": ["Fig5H.csv: heatmap, x = `rod_bend`, y = `chain_bend`, colour and cell number = `tipgap_sigma` (`temperature` = 2)",
        "(make_fig5_mechanism_panels.py, H)"],
 },
 "Figure6": {
  "B": ["Fig6B.csv: x = `pred_Tm_oof`, y = `exp_Tm`, colour = `group`; grey line = identity",
        "RMS and ρ text = the XGBoost row of Fig6C.csv (`OOF_RMS_C`, `OOF_Rho`)",
        "(Figure6_MachineLearning/make_figure6.py, panel_B)"],
  "C": ["Fig6C.csv: one bar per `Model`, bar length = `OOF_Rho`; circle + error bar = `Mean_Fold_Rho` ± `Std_Fold_Rho`",
        "right-hand labels = `OOF_Rho` and `OOF_RMS_C`; green bar = the multivalent model. Not drawn: `n`",
        "(make_figure6.py, panel_C)"],
  "D": ["Fig6D.csv: one row per site i = 0–10: x = `SHAP_helix<i>` (density-jittered swarm),",
        "colour = `featval_helix<i>` scaled between its 5th and 95th percentiles (blue low, red high). Not drawn: `code`, `Tm_C`",
        "(Supplementary/S6_SHAP_Values/recompute_shap_fulltrain.py)"],
  "E": ["Fig6E.csv, three sub-panels, each split by `edge` (True = edge, False = no edge), jittered points with median and IQR bars:",
        "left y = `Tm_meas_C`; middle y = `Tm_resid_gap_C`; right y = `Tm_resid_model_C`",
        "P values in the sub-titles are computed in make_figure6.py (not stored here). Not drawn: `code`, `Tm_mv_oof_C`, `mean_gap_nm`",
        "(make_figure6.py, panel_E)"],
 },
}

def collect(fig_name, short, slots):
    """rebuild Pseudo_Figures/<fig_name>/ from scratch; return composite slots + manifest lines + problems."""
    d = os.path.join(OUT, fig_name)
    # clear only the generated panel files (<short><letter>*.svg/.csv) and MANIFEST.txt; keep anything added by hand
    os.makedirs(d, exist_ok=True)
    for f in os.listdir(d):
        p = os.path.join(d, f)
        if os.path.isfile(p) and (f == 'MANIFEST.txt' or (f.startswith(short) and f.endswith(('.svg', '.csv')))): os.remove(p)
    comp, man, probs = [], [f"{fig_name}: panel -> source (regenerated by make_pseudo_figures.py)", ""], []
    for slot in slots:
        label, src = slot[0], slot[1]
        letter = label.split()[0].replace('–', '-')
        if isinstance(src, tuple):
            comp.append((label, src)); man.append(f"{letter:5s} {label}\n      [{src[0]}] {src[1] or ''}"); continue
        stem = f"{short}{letter}"
        if not os.path.exists(src):
            comp.append((label, src)); probs.append(f"{stem}: missing SVG {src}"); man.append(f"{letter:5s} {label}\n      MISSING {src}"); continue
        svg_dst = os.path.join(d, stem + '.svg'); shutil.copyfile(src, svg_dst); comp.append((label, svg_dst))
        csvs = slot[2] if len(slot) > 2 else []
        csvs = [csvs] if isinstance(csvs, str) else list(csvs)
        lines = [f"{letter:5s} {label}", f"      svg: {src}"]
        if not csvs: probs.append(f"{stem}: no CSV listed")
        for k, c in enumerate(csvs):
            if not os.path.exists(c): probs.append(f"{stem}: missing CSV {c}"); lines.append(f"      csv: MISSING {c}"); continue
            dst = os.path.join(d, stem + ('' if k == 0 else f'_{k + 1}_' + os.path.splitext(os.path.basename(c))[0]) + '.csv')
            shutil.copyfile(c, dst); lines.append(f"      csv: {c}  ->  {os.path.basename(dst)}")
        pl = PLOTTED.get(fig_name, {}).get(letter, [])
        if pl: lines += ["      plotted:"] + [f"        {p}" for p in pl]
        else: probs.append(f"{stem}: no plotted-column description in PLOTTED")
        man += lines
    open(os.path.join(d, 'MANIFEST.txt'), 'w', encoding='utf-8').write('\n'.join(man) + '\n')
    return comp, probs

def _copy(src, dst_dir, stem, csvs, probs):
    shutil.copyfile(src, os.path.join(dst_dir, stem + '.svg')); names = []
    csvs = [csvs] if isinstance(csvs, str) else list(csvs)
    for k, c in enumerate(csvs):
        if not os.path.exists(c): probs.append(f"{stem}: missing CSV {c}"); continue
        nm = stem + ('' if k == 0 else '_' + os.path.splitext(os.path.basename(c))[0]) + '.csv'
        shutil.copyfile(c, os.path.join(dst_dir, nm)); names.append(nm)
    return names

def collect_si():
    """rebuild Pseudo_Figures/SI/ (numbered SI figures) and SI/Extra/ (everything else) + README.md mapping."""
    d = os.path.join(OUT, 'SI'); ex = os.path.join(d, 'Extra')
    # clear only what this function generates (FigS* files and Extra/); anything else placed in SI/ by hand
    # (e.g. the oxDNA/ folder) is kept
    if os.path.isdir(d):
        for f in os.listdir(d):
            if f.startswith('FigS') and os.path.isfile(os.path.join(d, f)): os.remove(os.path.join(d, f))
    os.makedirs(ex, exist_ok=True)
    for f in os.listdir(ex):                       # only generated Extra_* files; hand-added files are kept
        if f.startswith('Extra_') and os.path.isfile(os.path.join(ex, f)): os.remove(os.path.join(ex, f))
    probs, comp, comp_x = [], [], []
    md = ['# Supplementary figures — file map', '',
          'Regenerated by `make_pseudo_figures.py SI` (run `python3 make_pseudo_figures.py SI` after changing any analysis).',
          'Each SI figure number maps to its SVG and CSV file(s) in this folder; letters (a, b, c) are panels of one figure.',
          'Figures not in the numbered list are in `Extra/`.', '',
          '| SI Fig. | Description | SVG | CSV | Source |', '|---|---|---|---|---|']
    def key(num): return (int(''.join(ch for ch in num if ch.isdigit())), num)
    for num, name, desc, svg, csvs in sorted([x for x in SI_FIGS if x[0]], key=lambda x: key(x[0])):
        stem = f"FigS{int(''.join(ch for ch in num if ch.isdigit())):02d}{''.join(ch for ch in num if ch.isalpha())}_{name}"
        if not os.path.exists(svg):
            probs.append(f"{stem}: missing SVG {svg}"); md.append(f"| S{num} | {desc} | MISSING | | `{svg}` |"); continue
        cs = _copy(svg, d, stem, csvs, probs); comp.append((f"S{num}  {name}", os.path.join(d, stem + '.svg')))
        md.append(f"| S{num} | {desc} | `{stem}.svg` | " + '<br>'.join(f'`{c}`' for c in cs) + f" | `{svg}` |")
    md += ['', '## Extra/', '', '| File | Description | Source |', '|---|---|---|']
    for _, name, desc, svg, csvs in [x for x in SI_FIGS if not x[0]]:
        stem = f"Extra_{name}"
        if not os.path.exists(svg): probs.append(f"{stem}: missing SVG {svg}"); continue
        _copy(svg, ex, stem, csvs, probs); comp_x.append((name, os.path.join(ex, stem + '.svg')))
        md.append(f"| `{stem}.svg` | {desc} | `{svg}` |")
    # the generated map goes to README.md only if README.md is absent or is itself a generated map; a hand-written
    # README.md is never overwritten (the map then goes to SI_FILE_MAP.md)
    rd = os.path.join(d, 'README.md')
    ours = not os.path.exists(rd) or open(rd, encoding='utf-8', errors='ignore').read().startswith(md[0])
    if not ours: probs.append("SI/README.md is hand-written: generated map written to SI/SI_FILE_MAP.md instead")
    open(rd if ours else os.path.join(d, 'SI_FILE_MAP.md'), 'w', encoding='utf-8').write('\n'.join(md) + '\n')
    compose('SI', 'Supplementary figures (numbered)', comp, cols=3)
    compose('SI_Extra', 'Supplementary — extra figures', comp_x, cols=3, out=os.path.join(ex, 'Extra_pseudo.svg'))
    return comp, comp_x, probs

if __name__ == "__main__":
    args = sys.argv[1:]
    only = [a for a in args if not a.startswith('--')] or list(FIGS)
    bad = [a for a in only if a not in FIGS]
    if bad: sys.exit(f"unknown figure(s): {bad}; choose from {list(FIGS)}")
    failed = [] if '--no-regen' in args else regenerate(only, '--refit' in args, '--notebooks' in args)
    print("Pseudo-figures (current paper layout):")
    problems = []
    for name in only:
        title, cols, slots = FIGS[name]
        if name == 'SI':
            comp, comp_x, probs = collect_si(); problems += probs
            print(f"  SI: {len(comp)} numbered panels -> {OUT}/SI/ (+ {len(comp_x)} in SI/Extra/), map in {OUT}/SI/README.md")
            continue
        comp, probs = collect(name, 'Fig' + name.replace('Figure', ''), slots)
        out = compose(name, title, comp, cols=cols)
        nreal = sum(1 for _, s in comp if not isinstance(s, tuple) and os.path.exists(s))
        print(f"  {name}: {len(comp)} slots ({nreal} real panels) -> {out}, panels in {OUT}/{name}/")
        problems += probs
    if failed: print("\nGENERATOR FAILURES (panels kept from the previous run):\n  " + "\n  ".join(failed))
    if problems: print("\nPROBLEMS:\n  " + "\n  ".join(problems))
    if failed or problems: sys.exit(1)
