"""
SI: root-mean-square fluctuation (RMSF) of the simulated DNA tile monomer (oxDNA, 300M steps; edge code 55555555555).
SVG: the oxView mean structure coloured by RMSF (Monomer_oxDNA/Monomer.png) with its colour bar, embedded as images
(the 3D render is raster), plus a per-nucleotide RMSF trace.  CSV: per-nucleotide RMSF (nm) from Monomer_deviations.json.
Writes rmsf_tile.svg + rmsf_tile.csv here.
"""
import os, json, base64, io, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); M = os.path.join(HERE, "Monomer_oxDNA")
r = np.asarray(json.load(open(os.path.join(M, "Monomer_deviations.json")))["RMSF (nm)"], float)
pd.DataFrame(dict(nucleotide_index=np.arange(len(r)), RMSF_nm=np.round(r, 4))).to_csv(os.path.join(HERE, "rmsf_tile.csv"), index=False)
fig = plt.figure(figsize=(10, 4.2))
a1 = fig.add_axes([0.0, 0.05, 0.42, 0.9]); a1.imshow(plt.imread(os.path.join(M, "Monomer.png"))); a1.axis('off')
a1.set_title('mean structure coloured by RMSF (oxView)', fontsize=9)
a2 = fig.add_axes([0.43, 0.15, 0.03, 0.7]); a2.imshow(plt.imread(os.path.join(M, "Monomer_colorbar.png"))); a2.axis('off')
a3 = fig.add_axes([0.55, 0.15, 0.42, 0.7]); a3.plot(np.arange(len(r)), r, lw=.4, color='#3B6AA0')
a3.set(xlabel='nucleotide index', ylabel='RMSF (nm)'); a3.set_title(f'per-nucleotide RMSF (median {np.median(r):.2f} nm)', fontsize=9)
a3.spines[['top', 'right']].set_visible(False)
fig.savefig(os.path.join(HERE, "rmsf_tile.svg"), dpi=150); plt.close(fig)
print(f"RMSF tile: {len(r)} nucleotides, median {np.median(r):.2f} nm")
