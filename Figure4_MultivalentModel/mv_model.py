"""
Multivalent binding model for tile-tile interfaces (Figure 4).  Imported by fit_model.py and make_figure4.py.

    K_i(T)    = exp(-(dH_i - T dS_i) / RT)        per-duplex SantaLucia-1998 NN, dS salt-corrected
                                                   (+0.368 (L-1) ln[Na+]eq, [Na+]eq = 3.79 sqrt(0.0113 M))
    C_eff,i   = C_ref * exp(-d_i / lambda) * exp(-kappa (N - 4))
                d_i = mean distance of duplex i to the other duplexes (6 nm site pitch, twin duplexes +-1.5 nm;
                      identical geometry to the Fig 4 / Table S1 phenomenological models), N = occupied sites
    K_a       = (1 / C_ref) [ prod_i (1 + K_i C_eff,i) - 1 ]           (one bimolecular docking step)
    mass balance (each tile at C = 4.5 nM):  f = ((2a + 1) - sqrt(4a + 1)) / (2a),  a = K_a C
    readout   theta = f [ w + (1 - w) <p_rep> ],  p_rep = K_r C_eff,r / (1 + K_r C_eff,r) of the reporter-site
                (site 5) duplexes;  w = 0.57 (intermediate; Figure 2 Panel I), w = 1 docking
    Tm: theta(Tm) = 1/2.

Parameter vector p = [log10 C_ref (M), lambda (nm), kappa (per site)]; shorter vectors drop the later terms
(length 1: one uniform C_eff; length 2: + distance).  Data: the canonical 236-interface set of
Supplementary/TableS1_Model_Comparison_All/unified_models.py (same constructs, sequences, geometry and folds).
"""
import os, numpy as np
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
_UM = os.path.join(REPO, 'Supplementary', 'TableS1_Model_Comparison_All', 'unified_models.py')
_src = open(_UM).read(); _cut = 'rows=[]'
assert _cut in _src, 'unified_models.py layout changed: marker "rows=[]" missing'
UM = {'__file__': _UM}
_cwd = os.getcwd(); exec(_src[:_src.index(_cut)], UM); os.chdir(_cwd)
bonds, AGG = UM['bonds'], UM['agg']
CODES = AGG.code.values; Y_C = UM['Y']; N_ALL = UM['n']

R = 1.987204e-3                        # kcal/(mol K)
C_TILE = 4.5e-9                        # each tile (M)
W_INT = 0.57
LN_NAEQ = np.log(3.79 * np.sqrt(0.0113))
LEN = {'5': 5, '6': 6, '7': 7, '9': 9, 'A': 7, 'B': 7}
BOUNDS = ([-8.0, 1.0, -2.0], [1.0, 1e4, 2.0])
P0 = [-2.0, 13.0, 0.05]

def arrays(blist):
    """padded per-duplex arrays for a list of bonds() outputs."""
    n = len(blist); M = max(len(b) for b in blist)
    A = {k: np.zeros((n, M)) for k in ('H', 'S', 'D')}
    A['mask'] = np.zeros((n, M), bool); A['rep'] = np.zeros((n, M), bool); A['N'] = np.zeros(n)
    for i, b in enumerate(blist):
        z = np.array([x[0] for x in b])
        for k, x in enumerate(b):
            A['H'][i, k] = x[1]
            A['S'][i, k] = x[2] / 1000.0 + 0.368e-3 * (LEN.get(x[4], 7) - 1) * LN_NAEQ
            A['mask'][i, k] = True
            A['rep'][i, k] = int(round(x[0] / 6.0)) == 5
            A['D'][i, k] = np.mean(np.abs(np.delete(z, k) - z[k])) if len(z) > 1 else 0.0
        A['N'][i] = len(b) / 2
    A['nrep'] = np.maximum(A['rep'].sum(1), 1); A['hasrep'] = A['rep'].sum(1) > 0
    return A

def arrays_from_codes(codes):
    bl = [bonds(c) for c in codes]
    assert all(b is not None for b in bl), 'a code has no sequence for one of its positions'
    return arrays(bl)

def _sub(A, idx):
    return {k: (v[idx] if isinstance(v, np.ndarray) else v) for k, v in A.items()}

def ceff(p, A):
    lc = p[0]; lam = p[1] if len(p) > 1 else np.inf; kap = p[2] if len(p) > 2 else 0.0
    return 10 ** lc * np.exp(-A['D'] / lam) * np.exp(-kap * (A['N'] - 4))[:, None]

def state(T, p, A, w=W_INT):
    """theta, docked fraction f, per-duplex conditional occupancy x/(1+x) at temperatures T (K, one per row)."""
    T = np.asarray(T, float)
    x = np.exp(np.clip(-(A['H'] - T[:, None] * A['S']) / (R * T[:, None]), -700, 700)) * ceff(p, A)
    x = np.where(A['mask'], x, 0.0)
    a = np.expm1(np.sum(np.log1p(x), 1)) / 10 ** p[0] * C_TILE
    f = np.where(a > 1e-12, ((2 * a + 1) - np.sqrt(4 * a + 1)) / (2 * np.maximum(a, 1e-300)), a)
    occ = x / (1 + x)
    prep = np.where(A['hasrep'], np.sum(np.where(A['rep'], occ, 0), 1) / A['nrep'], 1.0)
    return f * (w + (1 - w) * prep), f, occ

def theta(T, p, A, w=W_INT): return state(T, p, A, w)[0]

def tm(p, A, w=W_INT):
    """Tm in K for every row of A (vectorized bisection; theta falls with T)."""
    n = len(A['N']); lo = np.full(n, 230.0); hi = np.full(n, 420.0)
    for _ in range(48):
        mid = (lo + hi) / 2; up = theta(mid, p, A, w) > 0.5
        lo = np.where(up, mid, lo); hi = np.where(up, hi, mid)
    return (lo + hi) / 2

def width(p, A, w=W_INT, h=0.05):
    """logistic-equivalent width 1/(4|dtheta/dT|) at Tm (K = C)."""
    T = tm(p, A, w)
    return 1.0 / (4.0 * np.abs((theta(T + h, p, A, w) - theta(T - h, p, A, w)) / (2 * h)))

def fit(A, Y_K, npar=3, w=W_INT, idx=None, p0=None):
    idx = np.arange(len(Y_K)) if idx is None else np.asarray(idx)
    As, Ys = _sub(A, idx), Y_K[idx]
    lb, ub = BOUNDS[0][:npar], BOUNDS[1][:npar]
    return least_squares(lambda p: tm(p, As, w) - Ys, (P0 if p0 is None else p0)[:npar], bounds=(lb, ub), diff_step=1e-4)

def group(code):
    s = set(code) - {'0'}
    if s & {'A', 'B'}: return 'A/B'
    if '6' in s: return '6-nt'
    if '7' in s: return '7-nt'
    if s <= {'5'}: return '5-nt'
    return 'other'

def n_sites(code): return sum(ch != '0' for ch in code)

def fit_ceff_per_construct(A, Y_K, w=W_INT):
    """log10 C_eff (M) for each row separately: one uniform C_eff per construct (no distance/valency terms),
    solved so that theta(measured Tm) = 1/2 (vectorized bisection; theta increases with C_eff)."""
    T = np.asarray(Y_K, float); n = len(T); lo = np.full(n, -14.0); hi = np.full(n, 3.0)
    for _ in range(60):
        lc = (lo + hi) / 2; C = 10 ** lc
        x = np.exp(np.clip(-(A['H'] - T[:, None] * A['S']) / (R * T[:, None]), -700, 700)) * C[:, None]
        x = np.where(A['mask'], x, 0.0)
        a = np.expm1(np.sum(np.log1p(x), 1)) / C * C_TILE
        f = np.where(a > 1e-12, ((2 * a + 1) - np.sqrt(4 * a + 1)) / (2 * np.maximum(a, 1e-300)), a)
        occ = x / (1 + x)
        prep = np.where(A['hasrep'], np.sum(np.where(A['rep'], occ, 0), 1) / A['nrep'], 1.0)
        up = f * (w + (1 - w) * prep) > 0.5
        hi = np.where(up, lc, hi); lo = np.where(up, lo, lc)
    return (lo + hi) / 2

def nn_twostate_tm(A, C_tile=C_TILE):
    """additive NN two-state Tm (all duplexes, 1 M reference): sum dH / (sum dS + R ln(C/2))  [K]"""
    return np.sum(np.where(A['mask'], A['H'], 0), 1) / (np.sum(np.where(A['mask'], A['S'], 0), 1) + R * np.log(C_tile / 2))

def site_dG_features(code):
    """ML features (Figure 6 and SHAP): 11 per-site values = mean salt-corrected NN dG at 37 C (kcal/mol) of the
    site's two duplexes, 0 for an empty site.  Same TH table / slot convention as bonds()."""
    TH = UM['TH']; x = np.zeros(11)
    for i, ch in enumerate(code):
        if ch == '0': continue
        ch = '5' if ch == 'X' else ch; dg = []
        for h in (2 * i + 1, 2 * i + 2):
            dH, dS = TH[ch][h]; dS = dS + 0.368 * (LEN[ch] - 1) * LN_NAEQ
            dg.append(dH - 310.15 * dS / 1000.0)
        x[i] = np.mean(dg)
    return x
