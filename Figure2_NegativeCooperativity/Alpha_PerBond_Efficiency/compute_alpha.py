import os as _os
def _find_repo():
    d=_os.path.dirname(_os.path.abspath(__file__))
    while d!='/' and not _os.path.exists(_os.path.join(d,'Filtered_DF_Tm.csv')):
        d=_os.path.dirname(d)
    return d
__REPO__=_find_repo()
"""
Per-bond avidity coefficient alpha for each multivalent interface.

Model:   dG_avg^poly = alpha * dG^mono
         N * dG_avg^poly = dG_N^poly = alpha * N * dG^mono

  dG_avg^poly = (experimental interface dG at 37 C) / (number of bonds N)   -- realized per-bond dG
  dG^mono     = the dG of a single (monovalent) bond.  Primary estimate = MEDIAN of the interface's
                individual bond dG; mean and spread (std/min/max) are also reported because the
                bonds are heterogeneous (different lengths / sequences).
  alpha       = dG_avg^poly / dG^mono     (per interface; alpha<1 => each bond realizes less than its
                monovalent dG, i.e. negative cooperativity / avidity penalty)

Monovalent bond dG = salt-corrected SantaLucia dG_37 (G_corrected), i.e. the same conditions /
temperature as the experimental dG, so units (kcal/mol) and T (37 C) match.

Robust experimental-dG set: Filtered_DF_HS rows with Reliable & high-confidence, >=2 replicates,
aggregated per construct by median.

Run from the Tile_TileResults directory.  Writes alpha_per_interface.csv and alpha_figure.svg here.
"""
import os, math, warnings, numpy as np, pandas as pd, nbformat
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
PROJ = __REPO__
OUT = os.path.join(PROJ, "Figure2_NegativeCooperativity", "Alpha_PerBond_Efficiency")
os.chdir(PROJ)

# ---- rebuild the NN per-bond thermo tables (salt-corrected G_corrected) from Tile_Tile_Master ----
nb = nbformat.read("Tile_Tile_Master.ipynb", as_version=4)
G = {"pd": pd, "np": np, "math": math}
def run(i): exec(nb.cells[i].source, G)
for i in (3, 4, 6, 7, 8):          # NN dicts, calc_unified_nn, saltCorrection, add_Pairs, add_nn_thermo_columns
    run(i)
for i in (20, 21):                 # constant / constant_2 (A,B helices)
    run(i)
G["nt7"] = pd.read_csv("nt7.csv", index_col=0); G["nt6"] = pd.read_csv("nt6.csv", index_col=0)
G["nt5"] = pd.read_csv("nt5.csv", index_col=0); G["nt9"] = pd.read_csv("nt9.csv", index_col=0)
G["nt9_stack"] = pd.read_csv("nt9_stack.csv", index_col=0)
run(30)                            # builds nt5_1/2, nt6_1/2, nt7_1/2, nt9_1/2

# ---- EACH BINDING SITE = TWO duplex bonds (add_Pairs sums two helices per code position).
# Each bond's monovalent energy comes ONLY from its OWN sequence (SantaLucia NN dG37, computed by
# calc_unified_nn) and its OWN full salt correction  dG_salt = dG_1M - 0.114*(L-1)*ln(3.3*sqrt([Mg])).
# No coaxial-stacking term and no half-splitting: each bond is an isolated, salt-corrected duplex.
LEN = {"5": 5, "X": 5, "6": 6, "7": 7, "9": 9, "A": 7, "B": 7}
nt1 = {"5": G["nt5_1"], "6": G["nt6_1"], "7": G["nt7_1"], "9": G["nt9_1"]}  # per-helix rows: Sequence + dG37
const1 = G["constant_1"]
calc = G["calc_unified_nn"]
ionic = 3.3 * math.sqrt(0.0113)   # [Na+]_eq = 3.3*sqrt([Mg2+]); Mg = 11.3 mM
def _corr(L): return 0.114 * (L-1) * math.log(ionic)     # salt term for one duplex of length L

def _bond_seqs(code):
    """the two helix SEQUENCES (+length) for each bonded code position; None if unknown char."""
    sites = []
    for i, ch in enumerate(code):
        if ch == "0": continue
        L = LEN.get(ch)
        if L is None: return None
        if ch in "5X679":
            t = "5" if ch == "X" else ch
            s1 = str(nt1[t]["Sequence"][2*i]); s2 = str(nt1[t]["Sequence"][2*i+1])
        elif ch == "A":
            s1 = str(const1["Sequence"][0]); s2 = str(const1["Sequence"][1])
        elif ch == "B":
            s1 = str(const1["Sequence"][2]); s2 = str(const1["Sequence"][3])
        sites.append((s1, s2, L))
    return sites

def per_bond_dG(code):
    """salt-corrected monovalent dG (kcal/mol) for EACH individual bond (2 per site),
    computed directly from each bond's sequence via calc_unified_nn + full salt correction."""
    sites = _bond_seqs(code)
    if sites is None: return None
    out = []
    for s1, s2, L in sites:
        out += [calc(s1)["Delta_G_37 (kcal/mol)"] - _corr(L),
                calc(s2)["Delta_G_37 (kcal/mol)"] - _corr(L)]
    return out

# ---- provenance check: each bond's energy = NN(its sequence) - full salt correction ----
_ex = "00000X09990"
print("PROVENANCE CHECK for", _ex, "(each bond = calc_unified_nn(sequence) - 0.114*(L-1)*ln(3.3*sqrt([Mg]))):")
for s1, s2, L in _bond_seqs(_ex):
    for s in (s1, s2):
        nn = calc(s)["Delta_G_37 (kcal/mol)"]
        print(f"   seq={s:8s} L={L}  NN dG37={nn:+.3f}  - salt({_corr(L):+.3f}) = {nn-_corr(L):+.3f} kcal/mol")

# ---- robust experimental dG set ----
d = pd.read_csv("Filtered_DF_HS.csv"); d["code"] = d["Code"].astype(str).str.partition(".")[0]
sub = d[(d["Reliable"] == True) & (d["Confidence"] == "high")]
g = (sub.groupby("code").agg(Exp_dG=("Exp_DeltaG", "median"), N=("Num_Bonds", "first"),
                             n_reps=("Exp_DeltaG", "size")).reset_index())
g = g[g["n_reps"] >= 2].dropna(subset=["Exp_dG", "N"]).reset_index(drop=True)
g = g[~g["code"].str.contains("9")].reset_index(drop=True)   # ignore 9nt bonds (stacking handled separately)

rows = []
for _, r in g.iterrows():
    mono = per_bond_dG(r["code"])
    if not mono: continue
    mono = np.array(mono)
    n_sites = int(r["N"])            # non-zero code positions
    N = len(mono)                    # TOTAL bonds = 2 per site
    dG_avg_poly = r["Exp_dG"] / N
    med, mean, sd = np.median(mono), mono.mean(), mono.std(ddof=0)
    rows.append(dict(code=r["code"], n_sites=n_sites, N_bonds=N, n_reps=int(r["n_reps"]),
                     Exp_dG_37C_kcal_mol=round(r["Exp_dG"], 4),
                     dG_avg_poly_kcal_mol=round(dG_avg_poly, 4),
                     dG_mono_median_kcal_mol=round(med, 4),
                     dG_mono_mean_kcal_mol=round(mean, 4),
                     dG_mono_std_kcal_mol=round(sd, 4),
                     dG_mono_min_kcal_mol=round(mono.min(), 4),
                     dG_mono_max_kcal_mol=round(mono.max(), 4),
                     alpha_median=round(dG_avg_poly / med, 4),
                     alpha_mean=round(dG_avg_poly / mean, 4)))
a = pd.DataFrame(rows)
a.to_csv(os.path.join(OUT, "alpha_per_interface.csv"), index=False)

# global alpha = slope of dG_avg_poly vs dG_mono through the origin
x = a["dG_mono_median_kcal_mol"].values; y = a["dG_avg_poly_kcal_mol"].values
alpha_global = float(np.sum(x * y) / np.sum(x * x))
am = a["alpha_median"]
print(f"interfaces: {len(a)}")
print(f"alpha (per-interface, using MEDIAN monovalent dG): median={am.median():.3f}  mean={am.mean():.3f}  "
      f"std={am.std():.3f}  IQR=[{am.quantile(.25):.3f},{am.quantile(.75):.3f}]")
print(f"alpha (using MEAN monovalent dG): median={a['alpha_mean'].median():.3f}  mean={a['alpha_mean'].mean():.3f}")
print(f"alpha_global (through-origin slope, median mono): {alpha_global:.3f}")

# ---- figure ----
fig, ax = plt.subplots(1, 2, figsize=(13, 5.3))
ax[0].errorbar(x, y, xerr=a["dG_mono_std_kcal_mol"].values, fmt="o", ms=5, alpha=0.5,
               color="#4C78A8", ecolor="#B0C4DE", elinewidth=0.8, capsize=0, label="interfaces (±bond spread)")
xs = np.linspace(min(x.min(), 0), 0, 50)
ax[0].plot(xs, alpha_global * xs, "r-", lw=2, label=f"dG_avg^poly = alpha · dG^mono,  alpha={alpha_global:.2f}")
ax[0].plot(xs, xs, "k--", lw=1.1, label="alpha = 1 (full additivity)")
ax[0].set(xlabel="monovalent bond dG^mono (median, kcal/mol)",
          ylabel="realized per-bond dG_avg^poly = dG_exp / N (kcal/mol)",
          title="Per-bond avidity: realized vs monovalent dG (37 °C)")
ax[0].legend(fontsize=8.5, loc="lower right")
ax[1].hist(am, bins=24, color="#4C78A8", alpha=0.8, edgecolor="white")
ax[1].axvline(am.median(), color="C3", lw=2, label=f"median alpha = {am.median():.2f}")
ax[1].axvline(1.0, color="k", ls="--", lw=1.1, label="alpha = 1")
ax[1].set(xlabel="alpha  (per interface, = dG_avg^poly / dG^mono_median)", ylabel="interfaces",
          title="Distribution of per-interface alpha")
ax[1].legend(fontsize=8.5)
plt.tight_layout(); plt.savefig(os.path.join(OUT, "alpha_figure.svg"), bbox_inches="tight")
print("saved alpha_per_interface.csv and alpha_figure.svg in", OUT)
print(a.head(6).to_string(index=False))
