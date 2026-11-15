import os as _os
def _find_repo():
    d = _os.path.dirname(_os.path.abspath(__file__))
    while d != '/' and not _os.path.exists(_os.path.join(d, 'Filtered_DF_Tm.csv')):
        d = _os.path.dirname(d)
    return d
REPO = _find_repo()
"""
Figure 2 (UPDATED / Tm-anchored reframe).  v2 -- rewritten after review.

Replaces the van't Hoff-dG panels of the original Figure2_NegativeCooperativity (old A alpha/EM defs,
B, D).  Every panel is a direct measured Tm, a Tm difference vs an additive null, a measured
transition width, or the state-function dG_total.

ONE consistent robust set is used for every panel:
    Filtered_DF_HS rows, Reliable==True & Confidence=='high', median per construct,
    >= 2 replicates, 9-nt bonds excluded.
n is printed per panel; all panels draw from the same set (Panel F additionally needs a per-bond
NN reference and reports any construct it cannot place).

Panels
  B  interface Tm vs the NN Tm of its strongest single bond (pure nearest-neighbor; no fit)
  C  Tm vs binding sites, by bond length
  E  dTm per added site vs the ADDITIVE nearest-neighbor prediction (measured / predicted)
  D  measured transition width vs valency, vs the all-or-none (ΔH=ΣNN) prediction
  F  avidity enhancement log10(EM / C_partner), C_partner = 4.5 nM (one tile), from dG_total
  (G mechanism cartoon is an external 3D render, unchanged, not produced here)

HONEST SCOPE / what each panel does and does not show -- see README.txt.  In brief:
  * dG_total = ΔH_NN*(1 - T/Tm) + R*T*ln(C_T/4) is a STATE-FUNCTION extrapolation: it does not need
    the melt to be two-state in SHAPE, but it DOES assume ΔH_total = ΔH_NN (all bonds formed, dCp~0).
    R*T*ln(C_T/4) = -12.27 kcal/mol is a constant standard-state offset (identical for every construct,
    fixed by defining Tm at theta=1/2); it is NOT a per-construct "docking cost".
  * Panel D does NOT use the van't Hoff 'Enthalpy' column (ln K vs 1/T), which is noise-biased here
    (returns 111 for a true 147 kcal/mol on synthetic melts; per-curve Spearman(|Enthalpy|,dx) ~ -0.08;
    0.22 variance from signal).  It plots the directly-fit width dx and derives enthalpy from the width:
    dH_vH = 1.5*R*Tm^2/dx (bimolecular -- two tiles), a validated estimator (recovers 145 for true 147;
    0.55 variance from signal).  The all-or-none reference is the matching bimolecular 1.5*R*Tm^2/ΔH_NN.
    The melt is ~2.4x broader than all-or-none -> cooperative unit ~3 helices (4-5 sites) rising to ~7
    (9-10 sites), a FRACTION of the interface (< 2N); molecularity shifts the count by 1.5x, not the sign.
  * Panel E: 4->5 realizes ~95% of the additive prediction (no cooperativity there); the signal is the
    9->10 shortfall.  Diminishing dTm alone is expected from Tm saturation + weaker 5nt bonds.

Writes figure.svg + data.csv per panel folder, exact points to Figure_Source_Data/, combined svg.
Run from anywhere.
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
os.chdir(REPO)
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'svg.fonttype': 'none'})

R, CT, T37 = 1.987, 9e-9, 310.15
C_PARTNER = 4.5e-9                                  # one tile; avidity reference (matches original Fig2e)
DOCK = R * T37 * np.log(CT / 4) / 1000.0            # -12.27 kcal/mol: CONSTANT standard-state offset
BL = {5: '#E8A33D', 6: '#4DA167', 7: '#3B6AA0'}
DEG = u'℃'; ARR = u'→'; APPROX = u'≈'; ARR_UP = u'↑'; ARR_DN = u'↓'
# Panel letters of the MAIN figure (v3 layout, 2026-11-12):
#   A schematic | B strongest bond vs interface | C cooperative-unit fraction | D Tm vs sites
#   E avidity   | F DDG_coop                    | G mechanism schematic
# Folder names keep their historical letters (PanelC_Tm_vs_Sites = D, PanelF_Avidity = E, etc.).
LET = dict(strongest='B', frac='C', sites='D', avid='E', coop='F', width='SI', dTm='SI')   # see make_pseudo_figures.py for the assembled layout
F2 = "Figure2_Updated"; SRC = "Figure_Source_Data"; os.makedirs(SRC, exist_ok=True)

def load():
    d = pd.read_csv('Filtered_DF_HS.csv')
    d = d[(d['Reliable'] == True) & (d['Confidence'] == 'high')].copy()
    for c in ['Tm', 'dx', 'Predicted_deltaH', 'Predicted_deltaS', 'Num_Bonds', 'NucleotideBonds']:
        d[c] = pd.to_numeric(d[c], errors='coerce')
    d['code'] = d['Code'].astype(str).str.partition('.')[0]
    g = d.groupby('code').agg(Tm=('Tm', 'median'), dx=('dx', 'median'), dH=('Predicted_deltaH', 'median'),
                              dS=('Predicted_deltaS', 'median'), N=('Num_Bonds', 'first'),
                              nt=('NucleotideBonds', 'first'), nrep=('Tm', 'size')).dropna()
    g = g[(g['nrep'] >= 2) & (g['nt'] != 9)].copy()          # >=2 replicates; 9nt excluded
    g['nt'] = g['nt'].astype(int)
    g['dG_total'] = g['dH'] * (1 - T37 / (g['Tm'] + 273.15)) + DOCK
    return g

def swarm(ax, x, v, colors, rng, w=.3):
    v = np.asarray(v, float); j = rng.uniform(-.16, .16, len(v))
    ax.scatter(x + j, v, s=13, c=colors, alpha=.78, edgecolor='none', zorder=2)
    md = np.median(v); ax.plot([x - w, x + w], [md, md], 'k-', lw=2.3, zorder=3)
    ax.plot([x, x], [np.percentile(v, 25), np.percentile(v, 75)], 'k-', lw=1.1, zorder=3)
    return md

def save(df, panel_dir, src_name):
    df.to_csv(os.path.join(F2, panel_dir, 'data.csv'), index=False)
    df.to_csv(os.path.join(SRC, src_name), index=False)

g = load()
print("robust set: n=%d constructs (>=2 reps, 9nt excluded); nt counts %s"
      % (len(g), dict(g.nt.value_counts().sort_index())))
GRP = [(7, 4), (7, 5), (5, 9), (5, 10)]

# ---- strongest-single-bond NN Tm per construct (for Panel B) ----
# For each construct, the NN melting temperature of its STRONGEST single bond alone, at the same 9 nM
# (two tiles, 4.5 nM each) and salt.  Uses ONLY standard nearest-neighbor thermodynamics for a 5-7 bp
# duplex -- no two-state fit, no cooperative-unit, no effective-molarity assumption.  dH/dS from
# calc_unified_nn (Tile_Tile_Master); dS salt-corrected (von Ahsen, [Na+]eq = 3.79*sqrt(0.0113)).
def _strongest_single_bond_tm():
    import math, nbformat
    nb = nbformat.read('Tile_Tile_Master.ipynb', as_version=4)
    GG = {'pd': pd, 'np': np, 'math': math}
    def runc(i): exec(nb.cells[i].source, GG)
    for i in (3, 4, 6, 7, 8, 20, 21): runc(i)
    for t in ('nt7', 'nt6', 'nt5', 'nt9'): GG[t] = pd.read_csv(f'{t}.csv', index_col=0)
    GG['nt9_stack'] = pd.read_csv('nt9_stack.csv', index_col=0); runc(30)
    LEN = {'5': 5, 'X': 5, '6': 6, '7': 7, '9': 9, 'A': 7, 'B': 7}
    nt1 = {'5': GG['nt5_1'], '6': GG['nt6_1'], '7': GG['nt7_1'], '9': GG['nt9_1']}
    const1 = GG['constant_1']; calc = GG['calc_unified_nn']; lnI = math.log(3.79 * math.sqrt(0.0113))
    def seqs(code):
        out = []
        for i, ch in enumerate(code):
            if ch == '0': continue
            L = LEN.get(ch)
            if L is None: return None
            if ch in '5X679':
                t = '5' if ch == 'X' else ch
                out += [(str(nt1[t]['Sequence'][2 * i]), L), (str(nt1[t]['Sequence'][2 * i + 1]), L)]
            elif ch == 'A': out += [(str(const1['Sequence'][0]), L), (str(const1['Sequence'][1]), L)]
            elif ch == 'B': out += [(str(const1['Sequence'][2]), L), (str(const1['Sequence'][3]), L)]
        return out
    def bond_tm(seq, L):
        d = calc(seq); dS = d['Delta_S (cal/K/mol)'] + 0.3676 * (L - 1) * lnI      # salt-corrected dS
        return d['Delta_H (kcal/mol)'] * 1000 / (dS + R * np.log(CT / 4)) - 273.15
    def best(code):
        s = seqs(code)
        return max(bond_tm(q, L) for q, L in s) if s else np.nan
    return best
_best_bond = _strongest_single_bond_tm()
g['best_bond_tm'] = [_best_bond(c) for c in g.index]

# ---------------- Panel B : interface Tm vs strongest single bond ----------------
def panel_B(ax, rng):
    gb = g[g.nt.isin([5, 7])].dropna(subset=['best_bond_tm']).copy()
    xpos = {(7, 'b'): 0.0, (7, 'i'): 0.85, (5, 'b'): 2.15, (5, 'i'): 3.0}
    ax.axhline(22, color='#999', ls=':', lw=1); ax.text(1.5, 23, 'room temp', fontsize=6.5, color='#999', ha='center', va='bottom')
    rows = []
    for nt in (7, 5):
        s = gb[gb.nt == nt]
        mb = swarm(ax, xpos[(nt, 'b')], s.best_bond_tm.values, ['#b8b8b8'] * len(s), rng, .3)
        mi = swarm(ax, xpos[(nt, 'i')], s.Tm.values, [BL[nt]] * len(s), rng, .3)
        xa, xb = xpos[(nt, 'b')], xpos[(nt, 'i')]
        ax.annotate('', xy=(xb, mi - 2), xytext=(xa, mb + 2), arrowprops=dict(arrowstyle='-|>', color='#333', lw=1.2))
        ax.text((xa + xb) / 2, (mb + mi) / 2, f'+{(s.Tm - s.best_bond_tm).median():.0f}' + DEG,
                ha='center', va='center', fontsize=8, weight='bold', bbox=dict(boxstyle='round,pad=.15', fc='white', ec='none', alpha=.85))
    ax.set_xticks([xpos[(7, 'b')], xpos[(7, 'i')], xpos[(5, 'b')], xpos[(5, 'i')]])
    ax.set_xticklabels(['strongest\nbond', '7-nt\ninterface', 'strongest\nbond', '5-nt\ninterface'], fontsize=7.5)
    ax.set_ylabel(r'$T_m$ (' + DEG + ')'); ax.set_ylim(-28, 50); ax.set_xlim(-.5, 3.5)
    ax.set_title(LET['strongest'] + '   melts far above its best single bond', fontsize=9.3, loc='left', weight='bold')
    out = gb.reset_index()[['code', 'nt', 'best_bond_tm', 'Tm']].copy()
    out['nt'] = out['nt'].astype(int); out['gain_C'] = out['Tm'] - out['best_bond_tm']
    out.columns = ['code', 'nt', 'strongest_bond_tm_C', 'interface_tm_C', 'gain_C']
    return out.round({'strongest_bond_tm_C': 2, 'interface_tm_C': 2, 'gain_C': 2})

# ---------------- Panel C : Tm vs binding sites ----------------
def panel_C(ax, rng):
    rows = []
    for i, (nt, N) in enumerate(GRP):
        s = g[(g.nt == nt) & (g.N == N)]
        swarm(ax, i, s.Tm.values, [BL[nt]] * len(s), rng)
        for code, tm in s.Tm.items():
            rows.append(dict(code=code, nt=int(nt), n_sites=int(N), Tm_C=round(float(tm), 2)))
    ax.set_xticks(range(4)); ax.set_xticklabels([f'{N}\n(n={len(g[(g.nt==nt)&(g.N==N)])})' for nt, N in GRP])
    ax.axvline(1.5, color='#ccc', lw=1); ax.text(.5, 50.5, '7-nt', color=BL[7], ha='center'); ax.text(2.5, 50.5, '5-nt', color=BL[5], ha='center')
    ax.set_xlabel('Binding sites (interfaces)'); ax.set_ylabel(r'$T_m$ (' + DEG + ')'); ax.set_ylim(28, 52)
    ax.set_title(LET['sites'] + '   more bonds ' + ARR + r' higher $T_m$', fontsize=9.5, loc='left', weight='bold')
    return pd.DataFrame(rows)

# ---------------- Panel E : dTm per added site vs additive NN ----------------
def panel_E(ax, rng):
    g2 = g.copy()
    g2['Tm_add'] = g2['dH'] * 1000 / (g2['dS'] + R * np.log(CT / 4)) - 273.15   # additive-NN two-state Tm
    tm, tp = g2['Tm'], g2['Tm_add']; codes = set(g2.index); nb = lambda s: sum(ch != '0' for ch in s)
    rows = []
    for hi in sorted(codes):                                      # sorted -> deterministic
        for p, ch in enumerate(hi):
            if ch == '0': continue
            lo = hi[:p] + '0' + hi[p + 1:]
            if lo in codes and nb(lo) == nb(hi) - 1:
                rows.append(dict(lo=lo, hi=hi, regime=f'{nb(lo)}->{nb(hi)}', nt=int(g2.loc[hi, 'nt']),
                                 dTm_meas_C=round(float(tm[hi] - tm[lo]), 3),
                                 dTm_additive_C=round(float(tp[hi] - tp[lo]), 3)))
    pp = pd.DataFrame(rows).sort_values(['regime', 'lo', 'hi']).reset_index(drop=True)
    regimes = [('4->5', 7), ('9->10', 5)]
    regimes = [(r, nt) for r, nt in regimes if (pp.regime == r).sum() >= 3]
    ax.axhline(0, color='#888', ls='--', lw=1)
    top = 0
    for i, (reg, nt) in enumerate(regimes):
        s = pp[pp.regime == reg]
        md = swarm(ax, i - .17, s.dTm_meas_C.values, [BL[nt]] * len(s), rng, .14)
        mp = np.median(s.dTm_additive_C.values)
        ax.plot([i + .03, i + .31], [mp, mp], color='#555', lw=2.3, ls='-')         # additive prediction
        ax.scatter([i + .17] * len(s), s.dTm_additive_C.values, s=10, c='#bbb', alpha=.6, edgecolor='none', zorder=1)
        top = max(top, s.dTm_meas_C.max(), s.dTm_additive_C.max())
        pct = 100 * md / mp if mp else np.nan
        ax.text(i, -4.3, f'meas {md:+.1f}\nadd {mp:+.1f}\n({pct:.0f}%)', ha='center', fontsize=7.3)
    ax.set_xticks(range(len(regimes)))
    ax.set_xticklabels([f'{r.replace("->", ARR)}\n{nt}-nt' for r, nt in regimes])
    ax.set_xlabel('Binding sites (pairs)'); ax.set_ylabel(r'$\Delta T_m$ per added site (' + DEG + ')')
    ax.set_ylim(-5.5, top + 1.5)
    ax.plot([], [], color=BL[7], marker='o', ls='', label='measured'); ax.plot([], [], color='#555', lw=2.3, label='additive NN')
    ax.legend(fontsize=7, loc='upper right')
    ax.set_title(LET['dTm'] + '   additive except at high density (9' + ARR + '10)', fontsize=9.5, loc='left', weight='bold')
    return pp[pp.regime.isin([r for r, _ in regimes])].reset_index(drop=True)   # only plotted regimes

# ---------------- Panel D : transition width vs valency (vs all-or-none) ----------------
# BIMOLECULAR (two tiles): width-derived van't Hoff enthalpy dH = 1.5*R*Tm^2/dx, and an all-or-none
# melt of ΔH = ΣNN would have width dx = 1.5*R*Tm^2/ΔH_NN.  This width-derived ΔH is a validated,
# less-noisy estimator (recovers true ΔH on synthetic melts; 13% vs 20% replicate scatter; variance-
# from-signal 0.55 vs 0.22) -- used here INSTEAD of the noise-biased 'Enthalpy' (ln K vs 1/T) column.
def panel_D(ax, rng):
    gg = g.copy(); TmK = gg['Tm'] + 273.15
    gg['dHvh_width'] = 1.5 * R * TmK ** 2 / (gg['dx'] * 1000.0)              # width-derived vanHoff ΔH (bimol) kcal/mol
    gg['w_allornone_C'] = 1.5 * R * TmK ** 2 / (gg['dH'].abs() * 1000.0)    # dx if ΔH = ΣNN (bimolecular)
    gg['width_ratio'] = gg['dx'] / gg['w_allornone_C']
    gg['coop_unit_helices'] = gg['dHvh_width'] / (gg['dH'].abs() / (2 * gg['N']))
    Ns = sorted(gg.N.unique())
    for N in Ns:
        s = gg[gg.N == N]
        swarm(ax, N, s.coop_unit_helices.values, [BL[int(nt)] for nt in s.nt], rng, .42)
    xp = np.array(Ns)
    ax.plot(xp, 2 * xp, 'r--', lw=1.7, marker='s', ms=4)
    ax.text(xp[-2], 2 * xp[-2] + .4, 'all helices as\none unit (2N)', color='r', fontsize=6.8, ha='center', va='bottom')
    lo = gg[gg.N.isin([4, 5])].coop_unit_helices.median(); hi = gg[gg.N.isin([9, 10])].coop_unit_helices.median()
    wobs = gg.dx.median(); wpred = gg.w_allornone_C.median()
    ax.text(.03, .97, 'cooperative unit ' + APPROX + f' {lo:.0f}' + u'→' + f'{hi:.0f} helices\n(< 2N, rises with valency)\n'
            + f'[melt {wobs/wpred:.1f}' + u'×' + ' broader than all-or-none]',
            transform=ax.transAxes, ha='left', va='top', fontsize=7.3, color='#333')
    ax.set_xlabel('Binding sites'); ax.set_ylabel('cooperative unit (# helices)')
    ax.set_ylim(0, 2 * max(Ns) + 2)
    ax.set_title(LET['width'] + '   melts as a partial cooperative unit', fontsize=9.5, loc='left', weight='bold')
    out = gg.reset_index()[['code', 'nt', 'N', 'dx', 'dHvh_width', 'w_allornone_C', 'width_ratio', 'coop_unit_helices']].copy()
    out['nt'] = out['nt'].astype(int)
    out.columns = ['code', 'nt', 'n_sites', 'width_dx_C', 'dHvh_width_kcal', 'allornone_dx_C', 'width_ratio_obs_over_aon', 'coop_unit_helices']
    return out.round({'width_dx_C': 3, 'dHvh_width_kcal': 1, 'allornone_dx_C': 3, 'width_ratio_obs_over_aon': 2, 'coop_unit_helices': 2})

# ---------------- Panel C (main) : cooperative-unit FRACTION of the interface ----------------
# f_coop = dH_vH(width) / dH_NN(all bonds) = fraction of the interface's base-pairing enthalpy that melts as one
# unit (bimolecular width estimator, as in the SI width panel).  1 = all-or-none; observed ~0.3-0.45.
def panel_C_frac(ax, rng):
    gg = g.copy(); TmK = gg['Tm'] + 273.15
    gg['dHvh_width'] = 1.5 * R * TmK ** 2 / (gg['dx'] * 1000.0)
    gg['f_coop'] = gg['dHvh_width'] / gg['dH'].abs()
    rows = []
    ax.axhline(1, color='#c0392b', ls='--', lw=1.4, zorder=1)
    ax.text(3.45, 1.02, 'all-or-none (whole interface)', color='#c0392b', fontsize=6.8, ha='right', va='bottom')
    for i, (nt, N) in enumerate(GRP):
        s = gg[(gg.nt == nt) & (gg.N == N)]
        md = swarm(ax, i, s.f_coop.values, [BL[nt]] * len(s), rng)
        ax.text(i + .33, md, f'{md:.2f}', fontsize=7, va='center', ha='left')
        for code, r in s.iterrows():
            rows.append(dict(code=code, nt=int(nt), n_sites=int(N), width_dx_C=round(float(r.dx), 3),
                             dHvh_width_kcal=round(float(r.dHvh_width), 1), dH_NN_kcal=round(float(abs(r.dH)), 1),
                             f_coop=round(float(r.f_coop), 3)))
    ax.set_xticks(range(4)); ax.set_xticklabels([f'{N}\n(n={len(gg[(gg.nt==nt)&(gg.N==N)])})' for nt, N in GRP])
    ax.axvline(1.5, color='#ccc', lw=1)
    ax.text(.5, 1.13, '7-nt', color=BL[7], ha='center'); ax.text(2.5, 1.13, '5-nt', color=BL[5], ha='center')
    ax.set_xlabel('Binding sites'); ax.set_ylim(0, 1.2); ax.set_xlim(-.5, 3.5)
    ax.set_ylabel(r'cooperative fraction  $\Delta H_{vH}/\Delta H_{NN}$')
    ax.set_title(LET['frac'] + '   only part of the interface melts as one unit', fontsize=9.5, loc='left', weight='bold')
    return pd.DataFrame(rows)

# ---------------- Panel F (main E) : avidity = multivalent-model effective concentration ----------------
# Updated 2026-11-13: each interface's C_eff is fitted on its own to its measured Tm with the multivalent model of
# Figure 4 (one uniform C_eff per construct, DOCKING readout theta = f, as in panel F; Figure4_MultivalentModel/mv_model.py).
# Avidity enhancement = C_eff / C_partner (4.5 nM).  Replaces the dG_total / one-initiation EM estimate, which put the
# initiation penalty in a different place and gave ~10x lower values.
def panel_F(ax):
    import sys as _s; _s.path.insert(0, os.path.join(REPO, 'Figure4_MultivalentModel')); import mv_model as _mv
    codes = [c.replace('X', '5') for c in g.index]
    lc = _mv.fit_ceff_per_construct(_mv.arrays_from_codes(codes), g.Tm.values + 273.15, w=1.0)   # docking readout (Figure 2 uses docking throughout)
    gg = g.copy(); gg['log10_Ceff_M'] = lc; gg['avid'] = lc - np.log10(C_PARTNER)
    bins = np.arange(np.floor(gg.avid.min() * 5) / 5, np.ceil(gg.avid.max() * 5) / 5 + .1, .1)
    for nt in (7, 5, 6):
        v = gg[gg.nt == nt].avid.values
        if len(v): ax.hist(v, bins=bins, color=BL[nt], alpha=.8, label=f'{nt}-nt (n={len(v)})')
    med = gg.avid.median(); cmed = 10 ** gg.log10_Ceff_M.median() * 1e3
    ax.axvline(med, color='k', ls='--', lw=1)
    ax.text(med, ax.get_ylim()[1] * .9, f' median {med:.2f}\n ($C_{{eff}}$ {cmed:.1f} mM)', fontsize=7.5, va='top')
    ax.set_xlabel(r'avidity enhancement  $\log_{10}(C_{eff}/C_{\mathrm{partner}})$   ($C_{\mathrm{partner}}=4.5$ nM)')
    ax.set_ylabel('interfaces'); ax.legend(fontsize=7.5, loc='upper left')
    ax.set_title(LET['avid'] + r'   avidity, median $10^{%.1f}$' % med, fontsize=9.5, loc='left', weight='bold')
    print("Panel F (avidity, multivalent C_eff, docking): n=%d, median log10(C_eff/C) %.2f, C_eff %.2f mM" % (len(gg), med, cmed))
    out = gg.reset_index()[['code', 'nt', 'N', 'Tm', 'log10_Ceff_M', 'avid']].copy()
    out['nt'] = out['nt'].astype(int); out['N'] = out['N'].astype(int)
    out.columns = ['code', 'nt', 'n_sites', 'Tm_C', 'log10_Ceff_M', 'avidity_log10_Ceff_over_Cpartner']
    return out.round({'log10_Ceff_M': 4, 'avidity_log10_Ceff_over_Cpartner': 4})

# ---------------- Panel H : binding cooperativity free energy (DDG_coop) ----------------
# Drawn from PanelH_BindingCooperativity/data.csv (written by make_uniform_ceff_panel.py, run BEFORE this script).
# Per add-one-bond pair: C_eff fitted SEPARATELY to each construct's measured Tm (full multivalent model, no
# all-or-none); DDG_coop = -RT ln(C_eff,hi / C_eff,lo) at 37 C (> 0: the added bond lowers the effective
# concentration felt by every bond -> helps less), for three FRET readouts.  Rank-equivalent to the constant-C_eff
# Tm shortfall (also in data.csv).
def panel_H_coop(ax):
    fn = os.path.join(F2, 'PanelH_BindingCooperativity', 'data.csv')
    if not os.path.exists(fn):
        ax.axis('off'); ax.text(.5, .5, 'Panel H data missing:\nrun make_uniform_ceff_panel.py first', ha='center', transform=ax.transAxes)
        return pd.DataFrame()
    DC = pd.read_csv(fn, dtype={'lo': str, 'hi': str})
    if 'ddG_coop_kcal' not in DC:
        raise SystemExit('PanelH_BindingCooperativity/data.csv is stale (no ddG_coop_kcal): rerun make_uniform_ceff_panel.py')
    READ = [('docking', '#6250d6', 'docking')]          # main figure: docking FRET readout only (others in data.csv)
    DC = DC[DC.readout == 'docking']
    REGL = {'4->5': '4' + ARR + '5\n7nt added', '9->10': '9' + ARR + '10\n5nt added', '10->11': '10' + ARR + '11\n5nt added'}
    regs = [r for r in ('4->5', '9->10') if (DC.regime == r).any()]       # 10->11 (8 pairs, one construct) left out
    DC = DC[DC.regime.isin(regs)]
    rng_ = np.random.default_rng(0); off = [0.0]; INK2 = '#52514e'
    lim = DC.ddG_coop_kcal.abs().max() * 1.05
    ax.axhline(0, color='#0b0b0b', lw=1, zorder=1)
    for i, reg in enumerate(regs):
        if reg == '10->11': ax.axvspan(i - .45, i + .45, color='#e6e5e0', alpha=.55, zorder=0)
        for k, (name, col, lab) in enumerate(READ):
            v = DC[(DC.regime == reg) & (DC.readout == name)].ddG_coop_kcal.values
            if not len(v): continue
            x = i + off[k]
            ax.scatter(x + rng_.uniform(-.2, .2, len(v)), v, s=7, color=col, alpha=.5, edgecolor='none', zorder=2)
            md = np.median(v); ax.plot([x - .28, x + .28], [md, md], color='#0b0b0b', lw=2.2, zorder=3)
            ax.plot([x, x], np.percentile(v, [25, 75]), color='#0b0b0b', lw=1, zorder=3)
            ax.text(x, lim * 1.04, f'median {md:+.2f}', ha='center', va='bottom', fontsize=6.2, color='#0b0b0b')
    n = DC[DC.readout == 'docking'].groupby('regime').size()
    ax.set_xticks(range(len(regs))); ax.set_xticklabels([f'{REGL[r]}\n(n={n.get(r, 0)})' for r in regs], fontsize=7)
    ax.set(xlim=(-.55, len(regs) - .45), ylim=(-lim, lim * 1.2),
           ylabel=r'$\Delta\Delta G_{coop} = -RT\,\ln(C_{eff,hi}/C_{eff,lo})$ (kcal/mol)')
    ax.text(1, lim * .62, ARR_UP + ' added bond lowers $C_{eff}$\n(helps less)', ha='center', va='center', fontsize=6.2, color=INK2)
    ax.text(len(regs) - 1, -lim * .5, ARR_DN + ' added bond raises $C_{eff}$\n(helps more)', ha='center', va='center', fontsize=6.2, color=INK2)
    ax.grid(True, axis='y', color='#e6e5e0', lw=.6, zorder=0)
    ax.text(.02, .02, 'docking FRET readout', transform=ax.transAxes, fontsize=6.5, color=INK2)
    ax.set_title(LET['coop'] + r'   added bonds lower $C_{eff}$ at high valency', fontsize=9.5, loc='left', weight='bold')
    return DC[['lo', 'hi', 'regime', 'readout', 'log10Ceff_lo', 'log10Ceff_hi', 'Ceff_ratio_hi_lo', 'ddG_coop_kcal']].reset_index(drop=True)

def render(fn, panel_dir, src_name, figsize=(4.6, 4.3)):
    fig, ax = plt.subplots(figsize=figsize); rng = np.random.default_rng(0)
    df = fn(ax, rng) if fn.__code__.co_argcount == 2 else fn(ax)
    fig.tight_layout(); fig.savefig(os.path.join(F2, panel_dir, 'figure.svg')); plt.close(fig)
    save(df, panel_dir, src_name)

render(panel_B, 'PanelB_StrongestBond', 'Fig2U_B_StrongestBond.csv')
render(panel_C, 'PanelC_Tm_vs_Sites', 'Fig2U_C_Tm_by_Sites.csv')
render(panel_E, 'PanelE_dTm_per_added_site', 'Fig2U_E_dTm_per_added_site.csv')
render(panel_D, 'PanelD_CooperativeUnit', 'Fig2U_D_TransitionWidth.csv')
render(panel_F, 'PanelF_Avidity', 'Fig2U_F_Avidity.csv')

os.makedirs(os.path.join(F2, 'PanelC_CoopFraction'), exist_ok=True)
render(panel_C_frac, 'PanelC_CoopFraction', 'Fig2U_C_CoopFraction.csv')
os.makedirs(os.path.join(F2, 'PanelF_DDGcoop_Docking'), exist_ok=True)      # main-figure F as a standalone panel
render(panel_H_coop, 'PanelF_DDGcoop_Docking', 'Fig2U_F_DDGcoop_Docking.csv')

# combined MAIN figure (v3):  A schematic (wide) | B
#                             C (coop fraction) | D (Tm vs sites) | E
#                             F | G schematic (wide)
# SI (standalone svgs): PanelE_dTm_per_added_site, PanelD_CooperativeUnit (helix count), PanelH_UniformCeff_Width,
# PanelI_ReadoutWeight.
def placeholder(ax, title, lines):
    ax.axis('off')
    ax.add_patch(plt.Rectangle((.02, .06), .96, .84, fill=False, ec='#bbb', ls='--', lw=1.2, transform=ax.transAxes))
    ax.set_title(title, fontsize=9.5, loc='left', weight='bold')
    ax.text(.5, .55, lines, ha='center', va='center', fontsize=7.5, transform=ax.transAxes, color='#666')
    ax.text(.5, .16, '[keep existing artwork]', ha='center', fontsize=7, style='italic', transform=ax.transAxes, color='#999')

fig = plt.figure(figsize=(14, 12))
gs = GridSpec(3, 3, figure=fig, hspace=.5, wspace=.33, top=.94, bottom=.06, left=.06, right=.98)
rng = np.random.default_rng(0)
placeholder(fig.add_subplot(gs[0, 0:2]), 'A   tile-tile interface and FRET melt readout',
            'interface schematic (binding sites, 5/7-nt bonds, FRET pair)\n'
            '+ example melt curve: docking ' + ARR + ' reporter')
panel_B(fig.add_subplot(gs[0, 2]), rng)
panel_C_frac(fig.add_subplot(gs[1, 0]), rng)
panel_C(fig.add_subplot(gs[1, 1]), rng)
panel_F(fig.add_subplot(gs[1, 2]))
panel_H_coop(fig.add_subplot(gs[2, 0]))
placeholder(fig.add_subplot(gs[2, 1:3]), 'G   mechanism',
            'docking pays the translational cost once ' + ARR + ' later bonds see their partner at ' + APPROX + ' mM\n'
            'only a few neighbouring helices melt together; at high valency each added bond lowers $C_{eff}$')
fig.suptitle('Figure 2: multivalency stabilizes the interface far beyond any single bond; '
             'the melt is partially cooperative and added bonds cost effective concentration',
             fontsize=10.8, weight='bold', y=.995)
fig.savefig(os.path.join(F2, 'Figure2_Updated.svg'))
import sys as _sys
if '--png' in _sys.argv: fig.savefig(_sys.argv[_sys.argv.index('--png') + 1], dpi=130)
plt.close(fig)
print("wrote Figure2_Updated/ panels + combined MAIN svg (A-G v3 layout) + Figure_Source_Data CSVs")
