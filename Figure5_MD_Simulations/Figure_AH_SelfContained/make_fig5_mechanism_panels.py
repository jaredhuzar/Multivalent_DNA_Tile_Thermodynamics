"""
Standalone Figure 5 E-H mechanism panels, built from the mechanism CSVs already in plots_robust/mechanism/
(written by Analysis_Tm_trends.ipynb cells 4, 6, 10 from the GSD trajectories in this folder; no GSD needed here).

  E  unbound tip-tip gap per linker bending stiffness (chain_bend 4/10/20/40; rod 80; kT = 2): horizontal swarm of every
       sample with median (vertical bar) and IQR (horizontal bar)  <- panelC_histograms_raw.csv
  F  anchor-spacing spread sigma vs scaffold bending stiffness (rod_bend 4/20/80), three temperatures 1.25/1.5/1.75
       <- registration_summary.csv (error bars = s.d. over replicates)
  G  mean unbound tip gap, linker x scaffold grid (kT = 2)      <- heatmap_grid.csv
  H  tip-gap spread sigma, linker x scaffold grid (kT = 2)      <- heatmap_grid.csv
Writes plots_robust/mechanism_panels/<e|f|g|h>/panel_<x>.svg + panel_<x>.csv (exactly the plotted values).
Run from anywhere; called by regen_fig5_tm_panels.py.
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "plots_robust", "mechanism"); OUT = os.path.join(HERE, "plots_robust", "mechanism_panels")
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'svg.fonttype': 'none'})
INK2, GRID = '#52514e', '#e6e5e0'
def folder(x): d = os.path.join(OUT, x); os.makedirs(d, exist_ok=True); return d

# ---- E: unbound tip gap vs linker stiffness (horizontal swarm, median + IQR) ----
raw = pd.read_csv(os.path.join(SRC, "panelC_histograms_raw.csv"))
chains = [4, 10, 20, 40]; cols = plt.cm.Oranges(np.linspace(.35, .85, len(chains)))
def swarm_offsets(x, width=.38, nbins=40):
    """density-aware vertical offsets so points in dense x-ranges spread wider (beeswarm-like)."""
    edges = np.linspace(x.min(), x.max() + 1e-9, nbins + 1); b = np.clip(np.digitize(x, edges) - 1, 0, nbins - 1)
    cnt = np.bincount(b, minlength=nbins); mx = cnt.max(); off = np.zeros(len(x)); rng = np.random.default_rng(0)
    for k in range(nbins):
        m = np.where(b == k)[0]
        if len(m): off[m] = rng.permutation(np.linspace(-1, 1, len(m))) * width * min(1, len(m) / (.6 * mx))
    return off
fig, ax = plt.subplots(figsize=(2.9, 2.5)); rowsE = []
for i, (c, col) in enumerate(zip(chains, cols)):
    g = raw[raw.chain_bend == c].tip_gap.values
    ax.scatter(g, i + swarm_offsets(g), s=5, color=col, alpha=.75, edgecolor='none', zorder=2)
    md = np.median(g); q1, q3 = np.percentile(g, [25, 75])
    ax.plot([q1, q3], [i, i], color='k', lw=1.4, zorder=3); ax.plot([md, md], [i - .3, i + .3], color='k', lw=1.8, zorder=3)
    rowsE.append(dict(linker_bending_stiffness=c, n=len(g), median=round(md, 4), q1=round(q1, 4), q3=round(q3, 4), mean=round(g.mean(), 4), sd=round(g.std(), 4)))
ax.set_yticks(range(len(chains))); ax.set_yticklabels(chains)
ax.set(xlabel="Unbound tip gap (beads)", ylabel="Linker bending\nstiffness", xlim=(0, 8), xticks=[0, 2, 4, 6, 8], ylim=(-.6, len(chains) - .4))
fig.tight_layout(); d = folder("e"); fig.savefig(os.path.join(d, "panel_e.svg")); plt.close(fig)
raw[raw.chain_bend.isin(chains)][["chain_bend", "temperature", "pair_index", "tip_gap"]].rename(
    columns={"chain_bend": "linker_bending_stiffness", "tip_gap": "unbound_tip_gap_beads"}).to_csv(os.path.join(d, "panel_e.csv"), index=False)
pd.DataFrame(rowsE).to_csv(os.path.join(d, "panel_e_summary.csv"), index=False)

# ---- F: anchor-spacing spread vs scaffold bending, three temperatures ----
rs = pd.read_csv(os.path.join(SRC, "registration_summary.csv")); rs = rs[rs.rod_bend.isin([4, 20, 80])].copy()
Ts = sorted(rs.temperature.unique()); tc = plt.cm.plasma(np.linspace(.15, .8, len(Ts)))
fig, ax = plt.subplots(figsize=(5.0, 3.8))
for T, col in zip(Ts, tc):
    s = rs[rs.temperature == T].sort_values("rod_bend")
    ax.errorbar(s.rod_bend, s.intra_sigma, yerr=s.intra_sigma_reperr, marker='o', color=col, capsize=3, lw=1.6, label=f"kT = {T:g}")
ax.set_xscale('log'); ax.set_xticks([4, 20, 80]); ax.set_xticklabels(["4", "20", "80"])
ax.set(xlabel="scaffold (rod) bending stiffness", ylabel="anchor-spacing spread σ (± s.d. over replicates)")
ax.legend(fontsize=8, frameon=False); ax.grid(True, axis='y', color=GRID, lw=.6)
ax.set_title("Anchor-spacing fluctuation vs scaffold bending stiffness (spacing 4)", fontsize=9)
fig.tight_layout(); d = folder("f"); fig.savefig(os.path.join(d, "panel_f.svg")); plt.close(fig)
rs[["rod_bend", "temperature", "n_reps", "intra_sigma", "intra_sigma_reperr"]].sort_values(["temperature", "rod_bend"]) \
  .to_csv(os.path.join(d, "panel_f.csv"), index=False)

# ---- G / H: linker x scaffold heatmaps at kT = 2 ----
hg = pd.read_csv(os.path.join(SRC, "heatmap_grid.csv")); hg = hg[hg.temperature == 2.0]
CH, RO = sorted(hg.chain_bend.unique()), sorted(hg.rod_bend.unique())
for key, col, lab, cmap, ttl in (("g", "tipgap_mean", "mean unbound tip gap", "viridis_r", "Mean unbound tip gap (kT = 2)"),
                                 ("h", "tipgap_sigma", "tip-gap spread σ", "viridis_r", "Tip-gap spread σ (kT = 2)")):
    Z = hg.pivot(index="chain_bend", columns="rod_bend", values=col).reindex(index=CH, columns=RO).values
    fig, ax = plt.subplots(figsize=(4.6, 3.8))
    im = ax.imshow(Z, origin='lower', cmap=cmap, aspect='auto')
    for (i, j), v in np.ndenumerate(Z):
        ax.text(j, i, f"{v:.2f}", ha='center', va='center', fontsize=8, color='white' if v > np.nanmean(Z) else 'black')
    ax.set_xticks(range(len(RO))); ax.set_xticklabels(RO); ax.set_yticks(range(len(CH))); ax.set_yticklabels(CH)
    ax.set(xlabel="scaffold (rod) bending stiffness", ylabel="linker bending stiffness"); ax.spines[['top', 'right']].set_visible(True)
    cb = fig.colorbar(im, ax=ax, fraction=.05, pad=.03); cb.set_label(lab)
    ax.set_title(ttl, fontsize=9)
    fig.tight_layout(); d = folder(key); fig.savefig(os.path.join(d, f"panel_{key}.svg")); plt.close(fig)
    hg[["chain_bend", "rod_bend", "temperature", col] + (["tipgap_n_unbound"] if key == "g" else [])] \
      .sort_values(["chain_bend", "rod_bend"]).to_csv(os.path.join(d, f"panel_{key}.csv"), index=False)
print("Fig 5 E-H mechanism panels -> plots_robust/mechanism_panels/{e,f,g,h}/")
