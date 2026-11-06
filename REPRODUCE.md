# How to reproduce the figures

This repository is self-contained: raw plate-reader melting curves, processed tables, the analysis
notebooks, and one script per figure are all included. Below is everything needed to go from a fresh
clone to the published figures.

## 1. Setup
```bash
git clone git@github.com:jaredhuzar/Multivalent_DNA_TIle_Thermodynamics.git
cd Multivalent_DNA_TIle_Thermodynamics
pip install -r requirements.txt
```
Requires Python 3.9+ with: numpy, pandas, scipy, scikit-learn, matplotlib, seaborn, xgboost, shap,
nbformat, openpyxl (and python-docx for the model-comparison table).

## 2. Regenerate everything (one command)
```bash
bash run_all.sh
```
This runs every script-produced panel in dependency order and writes each `figure.svg` + `data.csv`
back into its own folder. It prints a PASS/FAIL line per script and exits non-zero on any failure
(expected: **17/17 PASS**; log in `run_all.log`). The two heavy steps are the SHAP retrain and the
11-model cross-validation (each ~1–3 min).

## 3. How the data flows
```
Melting_Results/  (raw FRET curves) + Tile_TileResultsKey*.csv
        │   Tile_Tile_Master.ipynb  (droop-corrected two-state fits, QC gates)
        ▼
Filtered_DF_Tm.csv / Filtered_DF_HS.csv / Filtered_DF_std.csv   (processed, committed)
        │   per-figure scripts (auto-detect repo root, read these tables)
        ▼
FigureN_*/<panel>/figure.svg + data.csv
```
The processed tables are committed, so the figure scripts run without re-fitting. To rebuild the tables
from the raw curves, run `Tile_Tile_Master.ipynb` top to bottom.

## 4. Reproduce a single panel
Each script auto-detects the repo root and writes into its own folder, so run it from anywhere:

| panel | command / source |
|---|---|
| Fig 2a | `Tile_Tile_Master.ipynb` → `Figure2_NegativeCooperativity/Fig2a_NNdG_ExpdG_saltfixed/` |
| Fig 2b | filter `Filtered_DF_Tm.csv` (`NucleotideBonds`, `Num_Bonds`) |
| Fig 2c | `python3 Figure2_NegativeCooperativity/Fig2c_Alpha_PerBond_Efficiency/compute_alpha.py` |
| Fig 2d | `Tile_Tile_Master.ipynb` → `Figure2_NegativeCooperativity/Fig2d_ddG/` |
| Fig 2e | `python3 Extras/Effective_Molarity/effective_molarity.py` then `.../Avidity/plot_avidity.py` |
| Fig 3b/3c/3e | `Tile_Tile_Master.ipynb` → `Figure3_BondClustering/<panel>/` |
| Fig 3 ΔH/ΔS | `python3 Extras/HS_vs_Distance_7nt_4bonds/hs_vs_distance.py` |
| Fig 4b | `TmPredictor_NewTm.ipynb` → `Figure4_PhenomenologicalModel/Phenomenological_Model/` |
| Fig 4c–4f | `python3 Figure4_PhenomenologicalModel/Theory_Scripts/*.py` |
| Fig 5 | `Figure5_MD_Simulations/simulation_notebooks/*.ipynb`; panels in `figure_panels/` (E–H binder-position data in `figure_panels/mechanism_EFGH/`) |
| Fig 6b/6c/6d | `Tile_Tile_Master.ipynb` (XGBoost cells) → `Figure6_MachineLearning/` |
| Fig 6e | `python3 Supplementary/S6_SHAP_Values/recompute_shap_fulltrain.py` |
| SI S2/S3/S10 | `python3 Supplementary/Alpha_Geometry_Controlled/alpha_geometry_ols.py`; `python3 Supplementary/add_ols_and_meangap.py` |
| SI S4 | `python3 Extras/HS_vs_Distance_7nt_4bonds/hs_vs_distance.py`; `.../add_ols_and_meangap.py` |
| SI S7 | oxDNA tile RMSF — analysis not included (external oxDNA trajectory) |
| SI S8/S9/S11/S12 | `python3 Supplementary/SI_Figures_regenerate.py` |
| Tested-models table | `python3 Supplementary/TableS1_Model_Comparison_All/models_and_pdf.py` then `pdf_report.py` |

## 5. Verify you got it right
Key values the regenerated data should reproduce:

| quantity | expected |
|---|---|
| Fig 2b Tm (7nt 4 / 5 sites; 5nt 9 / 10 sites) | 40.8 / 44.3 ; 33.9 / 34.1 °C |
| Fig 2c cooperativity α (same groups) | 0.42 / 0.36 ; 0.22 / 0.19 |
| Fig 2e avidity enhancement, median log₁₀(EM/C) (7nt ; 5nt) | 4.94 ; 5.40 |
| Fig 4b occupancy model (held-out) | RMSD 1.98 °C, ρ 0.88 |
| Fig 6 XGBoost (held-out) | R² 0.81, ρ 0.87 |
| Supp. best model (A2, N-scaled ΔS_Bond) | cv RMS 1.98 °C, ρ 0.884 |

## 6. Notes
- **Raw MD trajectories** (`*.gsd`, ~0.7 GB) are not included; Fig 5 reproduces from the processed
  `.npz` in `Figure5_MD_Simulations/processed_data_npz/` and the simulation notebooks.
- `Tile_Tile_Master.ipynb` is trimmed to the cells that build the tables and export the figures
  (exploratory/diagnostic cells removed). Its in-notebook save paths target the original layout; the
  standalone scripts above are the maintained, verified reproduction path.
- See `README.md` for the full folder map.
