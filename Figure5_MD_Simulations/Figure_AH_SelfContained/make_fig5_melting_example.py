"""
SI: example simulated melting curves -- mean number of bonds formed vs temperature for the 3-bond interface at
spacing 3 with a flexible and a rigid scaffold (FlexScaffold3bonds / RigidScaffold3bonds, sp3), from std_results_v2/*.npz.
Per seed: mean <bonds> at each temperature (crashed runs skipped); curve = mean over seeds, error = s.d. over seeds.
Tm = crossing of n_max/2 = 1.5 bonds (same convention as the Figure 5 Tm panels).
Writes plots_robust/mechanism_panels/melting_example/panel_melting_example.{svg,csv}.
"""
import os, glob, re, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "plots_robust", "mechanism_panels", "melting_example")
os.makedirs(OUT, exist_ok=True)
SER = [("Flexible scaffold", "FlexScaffold3bonds", '#1f77b4'), ("Rigid scaffold", "RigidScaffold3bonds", '#ff7f0e')]; SP = 3; HALF = 1.5
rows = []
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'svg.fonttype': 'none'})
fig, ax = plt.subplots(figsize=(5.0, 3.8))
for lab, stem, col in SER:
    per = {}
    for f in sorted(glob.glob(os.path.join(HERE, "std_results_v2", f"std_{stem}_sp{SP}_rep*_seed*.npz"))):
        seed = re.search(r"seed(\d+)", f).group(1); d = np.load(f, allow_pickle=True)
        for T, nb, cr in zip(d["temperatures"], d["nbound"], d["crashed"]):
            s = np.asarray(nb, float)
            if cr or s.size < 10: continue
            per.setdefault(float(T), {})[seed] = float(np.mean(s))
    Ts = np.array(sorted(per)); m = np.array([np.mean(list(per[t].values())) for t in Ts]); e = np.array([np.std(list(per[t].values())) for t in Ts])
    tm = np.nan
    for i in range(len(Ts) - 1):
        if (m[i] - HALF) * (m[i + 1] - HALF) <= 0 and m[i] != m[i + 1]:
            tm = Ts[i] + (HALF - m[i]) * (Ts[i + 1] - Ts[i]) / (m[i + 1] - m[i]); break
    ax.errorbar(Ts, m, yerr=e, marker='o', color=col, capsize=3, lw=1.6, label=f"{lab}  ($T_m$ {tm:.2f} kT)")
    ax.plot(tm, HALF, 'o', color=col, mec='k', ms=6, zorder=4)
    for t, mm, ee in zip(Ts, m, e): rows.append(dict(series=lab, stem=stem, spacing=SP, temperature_kT=t, mean_bonds=round(mm, 4), sd_over_seeds=round(ee, 4), n_seeds=len(per[t]), Tm_kT=round(tm, 4)))
ax.axhline(HALF, color='#999', lw=1, ls='--'); ax.text(ax.get_xlim()[1], HALF + .05, 'n_max / 2', ha='right', fontsize=7, color='#666')
ax.set(xlabel='temperature (kT)', ylabel='mean bonds formed ⟨n⟩', ylim=(0, 3.1))
ax.legend(fontsize=7.5, frameon=False, loc='lower left'); ax.grid(alpha=.3)
ax.set_title('Simulated melting: 3 bonds, spacing 3', fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "panel_melting_example.svg")); plt.close(fig)
pd.DataFrame(rows).to_csv(os.path.join(OUT, "panel_melting_example.csv"), index=False)
print("melting example ->", OUT, {r['series']: r['Tm_kT'] for r in rows})
