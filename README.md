# Multivalent DNA Thermodynamics — figure data & code

Self-contained reproduction of every figure in *Characterization and Prediction of Multivalent DNA
Thermodynamics* (main text + SI, 11/1/2026). **Data and code are organized one folder per figure.**

## Layout
- **Top level = shared inputs** (used by many figures; exec'd notebook cells read these by name, so they
  stay at the root): `Filtered_DF_Tm.csv` / `_high` / `Filtered_DF_HS.csv` / `Filtered_DF_std.csv`,
  `Experimental_Trends/construct_medians.csv`, `nt{5,6,7,9}.csv`, `nt9_stack.csv`,
  `DNA Sequences.xlsx`, `DNA Sequences - 6nt bonds.csv`, `Melting_Results/` (raw plate-reader curves) +
  `Tile_TileResultsKey*.csv`, and the notebooks `Tile_Tile_Master.ipynb`, `TmPredictor_NewTm.ipynb`
  (outputs stripped).
- **One folder per figure:** `Figure1_Platform_AFM/`, `Figure2_NegativeCooperativity/`,
  `Figure3_BondClustering/`, `Figure4_PhenomenologicalModel/`, `Figure5_MD_Simulations/`,
  `Figure6_MachineLearning/`, and `Supplementary/`. Each panel subfolder holds its `data.csv`,
  `figure.svg`, a `README`, and the script that makes it.

## Run it
Scripts auto-detect the repo root (they locate `Filtered_DF_Tm.csv`) and write back into their own
folder, so run any of them from anywhere:
```
python3 Figure2_NegativeCooperativity/Alpha_PerBond_Efficiency/compute_alpha.py
python3 Figure4_PhenomenologicalModel/Theory_Scripts/entropy_penalty_vs_spacing.py
python3 Supplementary/add_ols_and_meangap.py
```
Deps: python3 + numpy, pandas, scipy, scikit-learn, matplotlib, seaborn, xgboost, shap, nbformat, openpyxl.
Raw MD trajectories (`*.gsd`, ~0.7 GB) are **not** included; Fig 5 reproduces from the processed `.npz`.

## Figure 2 — Negative cooperativity  → `Figure2_NegativeCooperativity/`
`NNdG_ExpdG_saltfixed/` (2a) · `Alpha_PerBond_Efficiency/` (2c, compute_alpha.py) · `Fig2d_ddG/` (2d) ·
`Effective_Molarity/` + `Avidity/` (2e) · Tm panel (2b): filter `../Filtered_DF_Tm.csv`.

**Salt correction (SantaLucia).** The NN-predicted ΔG uses a monovalent-equivalent salt model:
`[Na⁺]_eq = 3.3·√[Mg²⁺]` with `[Mg²⁺] = 11.3 mM`, and `ΔG_salt = ΔG_1M − 0.114·n_steps·ln[Na⁺]_eq`,
where `n_steps = n_duplex·(L−1)` phosphate steps (each sticky-end duplex of length `L` bp contributes
`L−1`; a binding site has `n_duplex = 2` duplexes). This per-duplex convention is applied consistently
in the notebook (`saltCorrection`), the α/EM scripts, and the per-interface table. The salt term shifts
only the sequence-predicted ΔG (α, EM, Fig 2a/2d backing data, Fig 6 features); it does not touch any
measured Tm/ΔH/ΔS, the phenomenological model, or model/figure rankings.

## Figure 3 — Bond clustering  → `Figure3_BondClustering/`
`Fig3b_Tm_vs_MeanGap_7nt_4bonds/` · `Fig3c_RelativeOccupancy_vs_Distance/` · `Fig3e_Tm_vs_BondPosition_ConstantSeq/` · `HS_vs_Distance_7nt_4bonds/`.

## Figure 4 — Phenomenological model  → `Figure4_PhenomenologicalModel/`
`Phenomenological_Model/` (pred-vs-exp, rank-rank) · `Theory_Scripts/` + `theory_*`, `EntropyPenalty_*`,
`Theoretical*` curves · model fit: `../TmPredictor_NewTm.ipynb`. Model comparison: `../Supplementary/Model_Comparison_All/`.

## Figure 5 — Coarse-grained MD  → `Figure5_MD_Simulations/`
(Panels E–H binder-position data — tip-gap/anchor 4×4 matrices + per-frame swarm — are in `figure_panels/mechanism_EFGH/`.)
`simulation_notebooks/` (24 HOOMD notebooks) · `figure_panels/` (b–h) · `processed_data_npz/` (291 runs).

## Figure 6 — Machine learning  → `Figure6_MachineLearning/`
`Fig6b_XGBoost_RankRank/` · `Fig6c_XGBoost_StabilityScatter/` · `Fig6d_ModelComparison_SpearmanRho/` ·
`Fig6e_SHAP_Values/` (recompute_shap_fulltrain.py). Features: `../Tile_Tile_Master.ipynb`.

## Supplementary (11/1 numbering)  → `Supplementary/`
S2 OLS avidity `Alpha_Geometry_Controlled/`,`OLS_alpha_geometry/` · S3 OLS Tm `OLS_Tm_geometry/` ·
S4 ΔH/ΔS vs distance `HS_vs_MeanGap_7nt_4bonds/` (+ Fig3 `HS_vs_Distance`) · S5 MD (Figure5) ·
S6 SHAP (Figure6) · S7 oxDNA tile RMSF — analysis not included (needs the oxDNA trajectory, kept outside the repo) · S8 `S4_BondPosition_Melting/` ·
S9 `S2_BondLengthCount_Steepness/` · S10 OLS dG `OLS_dG_geometry/` · S11 `S5_BondLength_Entropy/` ·
S12 `S7_PredictedVsMeasured_Tm/` · Tested-models table `Model_Comparison_All/`.
Regenerators: `SI_Figures_regenerate.py`, `add_ols_and_meangap.py`. Also `CooperativeUnit_vs_BondCount/`,
`Extra_ReproducibleTm_UnphysicalHS/`, `Cooperativity_Tests_source/`.

## Reproducing each panel
`bash run_all.sh` regenerates every script-produced panel (verified: 17/17 pass). Individually:

| panel | how to generate |
|---|---|
| 2a | `python3 Figure2_NegativeCooperativity/NNdG_ExpdG_saltfixed/regen_nndg_vs_exp.py` |
| 2b | filter `Filtered_DF_Tm.csv` by `NucleotideBonds` & `Num_Bonds` (7nt 4/5, 5nt 9/10) |
| 2c | `python3 Figure2_NegativeCooperativity/Alpha_PerBond_Efficiency/compute_alpha.py` |
| 2d | `python3 Figure2_NegativeCooperativity/Fig2d_ddG/regen_fig2d_ddg.py` (backing data); panel `test2_measured_only.svg` is measured-only |
| 2e | `python3 Figure2_NegativeCooperativity/Effective_Molarity/effective_molarity.py` then `…/Avidity/plot_avidity.py` |
| 3b/3c/3e | `Tile_Tile_Master.ipynb` → `Figure3_BondClustering/<panel>/` |
| 3 ΔH/ΔS vs distance | `python3 Figure3_BondClustering/HS_vs_Distance_7nt_4bonds/hs_vs_distance.py` |
| 4b | `TmPredictor_NewTm.ipynb` → `Figure4…/Phenomenological_Model/` |
| 4c–4f | `python3 Figure4_PhenomenologicalModel/Theory_Scripts/*.py` |
| 5 | `Figure5_MD_Simulations/simulation_notebooks/*.ipynb` (plots in `figure_panels/`) |
| 6b/6c/6d | `Tile_Tile_Master.ipynb` (XGBoost cells) → `Figure6…/` |
| 6e | `python3 Figure6_MachineLearning/Fig6e_SHAP_Values/recompute_shap_fulltrain.py` |
| S2 | `python3 Supplementary/Alpha_Geometry_Controlled/alpha_geometry_ols.py`, `…/add_ols_and_meangap.py` |
| S3, S10 | `python3 Supplementary/add_ols_and_meangap.py` (OLS_Tm_geometry, OLS_dG_geometry) |
| S4 | `python3 Figure3_BondClustering/HS_vs_Distance_7nt_4bonds/hs_vs_distance.py`, `…/add_ols_and_meangap.py` |
| S6 | `python3 Figure6_MachineLearning/Fig6e_SHAP_Values/recompute_shap_fulltrain.py` |
| S7 | oxDNA tile RMSF — not included (external oxDNA trajectory) |
| S8, S9, S11, S12 | `python3 Supplementary/SI_Figures_regenerate.py` |
| tested-models table | `python3 Supplementary/Model_Comparison_All/models_and_pdf.py` then `pdf_report.py` |

Install deps with `pip install -r requirements.txt`.
