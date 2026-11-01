import os as _os
def _find_repo():
    d=_os.path.dirname(_os.path.abspath(__file__))
    while d!='/' and not _os.path.exists(_os.path.join(d,'Filtered_DF_Tm.csv')):
        d=_os.path.dirname(d)
    return d
__REPO__=_find_repo()
"""
Effective molarity (EM) of the avidity shortfall.

Attribute the gap between the measured interface free energy and the ONE-INITIATION sequence
reference to bonds 2..N each hybridizing intramolecularly at an effective local concentration EM
instead of the 1 M bimolecular standard:

    dG_Exp - dG_Seq(1 init)  =  (N - 1) * RT * ln(1 M / EM)
    =>  EM = exp( -(dG_Exp - dG_Seq(1init)) / ((N-1) * RT) )   [M]

  dG_Seq(1 init) = sum of per-bond (stacking + salt) + ONE initiation (not N).
                   The two tiles dock as ONE body: only the nucleating bond pays the bimolecular
                   initiation (translational/rotational entropy); bonds 2..N form intramolecularly.
  N              = total duplex bonds = 2 per binding site (add_Pairs).
  EM > 1 M  => favourable preorganization (positive avidity);
  EM < 1 M  => bonds effectively dilute (negative cooperativity).

Each bond's stacking+salt+initiation come ONLY from its own sequence (calc_unified_nn + full salt
correction), identical to compute_alpha.py.  Robust experimental-dG set, 9nt excluded.
Run from Tile_TileResults.  Writes effective_molarity.csv + .svg here.
"""
import os, math, warnings, numpy as np, pandas as pd, nbformat
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
PROJ = __REPO__
OUT = os.path.join(PROJ, "Figure2_NegativeCooperativity", "Effective_Molarity"); os.chdir(PROJ)

# ---- rebuild NN per-bond thermo (same as compute_alpha.py) ----
nb = nbformat.read("Tile_Tile_Master.ipynb", as_version=4)
G = {"pd": pd, "np": np, "math": math}
for i in (3, 4, 6, 7, 8):  run = exec(nb.cells[i].source, G)
for i in (20, 21):         exec(nb.cells[i].source, G)
G["nt7"] = pd.read_csv("nt7.csv", index_col=0); G["nt6"] = pd.read_csv("nt6.csv", index_col=0)
G["nt5"] = pd.read_csv("nt5.csv", index_col=0); G["nt9"] = pd.read_csv("nt9.csv", index_col=0)
G["nt9_stack"] = pd.read_csv("nt9_stack.csv", index_col=0)
exec(nb.cells[30].source, G)

LEN = {"5": 5, "X": 5, "6": 6, "7": 7, "9": 9, "A": 7, "B": 7}
nt1 = {"5": G["nt5_1"], "6": G["nt6_1"], "7": G["nt7_1"], "9": G["nt9_1"]}
const1 = G["constant_1"]; calc = G["calc_unified_nn"]; init_G = G["init_G"]
ionic = 3.3 * 0.0114
def _corr(L): return 0.114 * L * math.log(math.sqrt(ionic))
def _init(seq):                                   # per-bond bimolecular initiation dG (both ends)
    seq = seq.upper(); return sum(init_G["AT" if b in "AT" else "GC"] for b in (seq[0], seq[-1]))

def _bond_seqs(code):
    sites = []
    for i, ch in enumerate(code):
        if ch == "0": continue
        L = LEN.get(ch)
        if L is None: return None
        if ch in "5X679":
            t = "5" if ch == "X" else ch
            s1 = str(nt1[t]["Sequence"][2*i]); s2 = str(nt1[t]["Sequence"][2*i+1])
        elif ch == "A": s1, s2 = str(const1["Sequence"][0]), str(const1["Sequence"][1])
        elif ch == "B": s1, s2 = str(const1["Sequence"][2]), str(const1["Sequence"][3])
        sites.append((s1, s2, L))
    return sites

def bond_table(code):
    """per bond: salt-corrected dG (incl init) and the init part alone."""
    sites = _bond_seqs(code)
    if sites is None: return None
    dG, ini = [], []
    for s1, s2, L in sites:
        for s in (s1, s2):
            dG.append(calc(s)["Delta_G_37 (kcal/mol)"] - _corr(L))
            ini.append(_init(s))
    return np.array(dG), np.array(ini)

RT = 1.987204e-3 * 310.15     # kcal/mol at 37 C

# ---- robust experimental dG set (same filter as alpha) ----
d = pd.read_csv("Filtered_DF_HS.csv"); d["code"] = d["Code"].astype(str).str.partition(".")[0]
sub = d[(d["Reliable"] == True) & (d["Confidence"] == "high")]
g = (sub.groupby("code").agg(Exp_dG=("Exp_DeltaG", "median"), n_sites=("Num_Bonds", "first"),
                             n_reps=("Exp_DeltaG", "size")).reset_index())
g = g[g["n_reps"] >= 2].dropna(subset=["Exp_dG", "n_sites"])
g = g[~g["code"].str.contains("9")].reset_index(drop=True)

rows = []
for _, r in g.iterrows():
    bt = bond_table(r["code"])
    if bt is None: continue
    dG, ini = bt; N = len(dG)
    if N < 2: continue
    dG_seq_Ninit = dG.sum()
    dG_seq_1init = dG_seq_Ninit - ini.sum() * (N - 1) / N     # keep ONE (mean) initiation
    shortfall = r["Exp_dG"] - dG_seq_1init                     # >0 => measured weaker
    EM_bonds = math.exp(-shortfall / ((N - 1) * RT))
    EM_sites = math.exp(-shortfall / ((r["n_sites"] - 1) * RT)) if r["n_sites"] > 1 else np.nan
    rows.append(dict(code=r["code"], n_sites=int(r["n_sites"]), N_bonds=N, n_reps=int(r["n_reps"]),
                     Exp_dG=round(r["Exp_dG"], 3), dG_Seq_Ninit=round(dG_seq_Ninit, 3),
                     dG_Seq_1init=round(dG_seq_1init, 3), shortfall_kcal=round(shortfall, 3),
                     EM_per_bond_M=EM_bonds, EM_per_site_M=EM_sites))
a = pd.DataFrame(rows); a.to_csv(os.path.join(OUT, "effective_molarity.csv"), index=False)

em = a["EM_per_bond_M"]
print(f"interfaces: {len(a)}   (robust dG, 9nt excluded)")
print(f"shortfall dG_Exp - dG_Seq(1init): median={a.shortfall_kcal.median():+.2f} kcal/mol "
      f"(>0 => measured weaker than one-init reference)")
print(f"EM per bond (M):  median={em.median():.3e}  IQR=[{em.quantile(.25):.3e}, {em.quantile(.75):.3e}]")
print(f"                  = median {em.median()*1e3:.2f} mM ;  log10(EM) median={np.log10(em).median():+.2f}")
print(f"EM per site (M):  median={a['EM_per_site_M'].median():.3e}")
frac_below = (em < 1).mean()*100
print(f"interfaces with EM < 1 M (dilute / negative coop): {frac_below:.0f}%")

# ---- figure ----
fig, ax = plt.subplots(1, 2, figsize=(13, 5.1))
ax[0].hist(np.log10(em), bins=24, color="#4C78A8", alpha=0.85, edgecolor="white")
ax[0].axvline(0, color="k", ls="--", lw=1.2, label="EM = 1 M (no avidity)")
ax[0].axvline(np.log10(em.median()), color="C3", lw=2, label=f"median = {em.median()*1e3:.1f} mM")
ax[0].set(xlabel="log$_{10}$(effective molarity / M)", ylabel="interfaces",
          title="Effective molarity of bonds 2..N")
ax[0].legend(fontsize=9)
ax[1].scatter(a["N_bonds"], np.log10(em), s=22, alpha=0.5, color="#4C78A8")
ax[1].axhline(0, color="k", ls="--", lw=1.1)
ax[1].set(xlabel="N bonds (= 2 x sites)", ylabel="log$_{10}$(EM / M)",
          title="EM vs valency")
plt.tight_layout(); plt.savefig(os.path.join(OUT, "effective_molarity.svg"), bbox_inches="tight")
print("saved effective_molarity.csv + .svg in", OUT)
