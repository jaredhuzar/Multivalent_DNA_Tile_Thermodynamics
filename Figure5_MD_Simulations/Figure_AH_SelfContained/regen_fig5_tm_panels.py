"""
Regenerate the Fig 5 b-f Tm panels from the processed std_results_v2/*.npz, and sync the two mirror
copies of plots_robust_v2/.  This runs ONLY the Tm/thermo panel cells (8, 12, 14) of
Analysis_Tm_trends.ipynb -- NOT the mechanism/heatmap cells that need the (gitignored) GSD
trajectories -- so Fig 5 b-f are reproducible in run_all.sh from the committed npz alone.

Tm convention (round-2 review fixes): B(T) = equal-weight mean of the 3 seeds' per-seed means
(full series, no thinning); error = s.d. over the per-seed Tm crossings (non-crossing seeds censored
at the ladder bound); per-bond dH/dG/dS are per-seed (spread over seeds), not pooled.
"""
import os, shutil, numpy as np, nbformat, matplotlib
matplotlib.use("Agg")
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
def _repo(d):
    while d != "/" and not os.path.exists(os.path.join(d, "Filtered_DF_Tm.csv")):
        d = os.path.dirname(d)
    return d
REPO = _repo(HERE)

nb = nbformat.read("Analysis_Tm_trends.ipynb", as_version=4)
G = {"np": np}
for ci in (8, 12, 14):                       # panel builder, per-replicate CSVs, f_matched
    exec(nb.cells[ci].source, G)

# standalone Fig 5 E-H mechanism panels from the mechanism CSVs (plots_robust/mechanism_panels/)
import subprocess, sys as _sys
subprocess.run([_sys.executable, os.path.join(HERE, "make_fig5_mechanism_panels.py")], check=True)
subprocess.run([_sys.executable, os.path.join(HERE, "make_fig5_melting_example.py")], check=True)

# sync canonical plots_robust_v2/ -> the two mirror copies (byte-identical)
CANON = os.path.join(HERE, "plots_robust_v2")
for mir in (os.path.join(REPO, "Figure5_MD_Simulations", "plots_robust_v2"),
            os.path.join(REPO, "Figure_Source_Data", "plots_robust_v2")):
    for root, _dirs, files in os.walk(CANON):
        rel = os.path.relpath(root, CANON); dst = os.path.join(mir, rel)
        os.makedirs(dst, exist_ok=True)
        for f in files:
            shutil.copy2(os.path.join(root, f), os.path.join(dst, f))
print("Fig 5 b-f panels regenerated (cells 8/12/14) and mirrors synced")
