import os as _os
def _find_repo():
    d=_os.path.dirname(_os.path.abspath(__file__))
    while d!='/' and not _os.path.exists(_os.path.join(d,'Filtered_DF_Tm.csv')):
        d=_os.path.dirname(d)
    return d
__REPO__=_find_repo()
"""
Fig 1C "Fitted Melting Curves": the normalized two-state (Boltzmann) fit of every construct's
droop-corrected melting transition, overlaid.

Each curve is  fU(T) = (y - A1)/(A2 - A1)  with  y = A1 + (A2 - A1)/(1 + exp((Tm - T)/dx)),
i.e. the normalized fitted sigmoid (A1/A2 are the fitted baselines, stored normalized to 0/1;
Tm = x0 midpoint, dx = transition width) from the authoritative Filtered_DF_Tm fits
(fit_Not_Normalized). One curve per construct = per-`Code` median of (A1, A2, Tm, dx).

Same selection as the original notebook cell (Tile_Tile_MasterRegenerate.ipynb): keep constructs
with dx >= 1, 30 <= Tm <= 50, A1 <= 1; evaluate on 500 points over 26-55 C.

Outputs (self-contained; run from anywhere):
  Figure1_Platform_AFM/Melting Curves.svg                 (the figure)
  Figure1_Platform_AFM/Melting_Curves_data.csv            (raw plotted data, long format)
  Figure_Source_Data/Fig1C_MeltingCurves.csv              (same, for the source-data folder)
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
PROJ=__REPO__; os.chdir(PROJ)
OUT="Figure1_Platform_AFM"; os.makedirs(OUT, exist_ok=True)

def boltzmann(x, A1, A2, x0, dx):
    return A1 + (A2 - A1) / (1 + np.exp((x0 - x) / dx))

# per-construct median fit parameters (one curve per interface)
tm = pd.read_csv("Filtered_DF_Tm.csv"); tm["code"] = tm["Code"].astype(str).str.partition(".")[0]
med = tm.groupby("code").agg(A1=("A1","median"), A2=("A2","median"),
                             Tm=("Tm","median"), dx=("dx","median")).reset_index()

x_fit = np.linspace(26, 55, 500)
plt.figure(figsize=(6, 6))
rows = []
for _, r in med.iterrows():
    A1, A2, Tm, dx = r["A1"], r["A2"], r["Tm"], r["dx"]
    if dx < 1 or Tm < 30 or Tm > 50 or A1 > 1:        # original selection
        continue
    y = boltzmann(x_fit, A1, A2, Tm, dx)
    fU = (y - A1) / (A2 - A1)                          # normalized signal
    plt.plot(x_fit, fU, alpha=0.2)
    for T, f in zip(x_fit, fU):
        rows.append((r["code"], round(float(Tm), 4), round(float(dx), 4),
                     round(float(T), 4), round(float(f), 6)))
plt.xlim(25, 55); plt.ylim(-0.05, 1.05)
plt.xlabel("Temperature (°C)"); plt.ylabel("Normalized Signal")
plt.title("Fitted Melting Curves"); plt.tight_layout()
plt.savefig(os.path.join(OUT, "Melting Curves.svg"))

data = pd.DataFrame(rows, columns=["code", "Tm_C", "dx", "temperature_C", "normalized_signal"])
data.to_csv(os.path.join(OUT, "Melting_Curves_data.csv"), index=False)
os.makedirs("Figure_Source_Data", exist_ok=True)
data.to_csv("Figure_Source_Data/Fig1C_MeltingCurves.csv", index=False)
n_curves = data["code"].nunique()
print(f"Fig1C: {n_curves} fitted melting curves x {len(x_fit)} points -> "
      f"{OUT}/Melting Curves.svg + Melting_Curves_data.csv + Figure_Source_Data/Fig1C_MeltingCurves.csv")
