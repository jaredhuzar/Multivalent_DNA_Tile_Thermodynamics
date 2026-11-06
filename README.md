# Multivalent DNA Thermodynamics — figure data & code

Self-contained reproduction of every figure in *Characterization and Prediction of Multivalent DNA
Thermodynamics* (main text + SI, 11/1/2026). **Data and code are organized one folder per figure.**

## Layout
- **Top level = shared inputs** (used by many figures; exec'd notebook cells read these by name, so they
  stay at the root): `Filtered_DF_Tm.csv` / `_high` / `Filtered_DF_HS.csv` / `Filtered_DF_std.csv`,
  `Experimental_Trends/construct_medians.csv` (**rebuilt from the authoritative tables** by
  `Experimental_Trends/build_construct_medians.py`: Tm from `Filtered_DF_Tm`, dG/dH/dS from
  `Filtered_DF_HS`, geometry from the code — this is the source for the geometry OLS, so they use the
  same droop-corrected Tm and construct set as the rest of the paper), `nt{5,6,7,9}.csv`, `nt9_stack.csv`,
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
python3 Figure2_NegativeCooperativity/Fig2c_Alpha_PerBond_Efficiency/compute_alpha.py
python3 Figure4_PhenomenologicalModel/Theory_Scripts/entropy_penalty_vs_spacing.py
python3 Supplementary/add_ols_and_meangap.py
```
Deps: python3 + numpy, pandas, scipy, scikit-learn, matplotlib, seaborn, xgboost, shap, nbformat, openpyxl.
Raw MD trajectories (`*.gsd`, ~0.7 GB) are **not** included; Fig 5 reproduces from the processed `.npz`.

## Figure 2 — Negative cooperativity  → `Figure2_NegativeCooperativity/`
`Fig2a_NNdG_ExpdG_saltfixed/` (2a) · `Fig2b_Tm_vs_BindingSites_x_Length/` (2b) ·
`Fig2c_Alpha_PerBond_Efficiency/` (2c, compute_alpha.py) · `Fig2d_ddG/` (2d) ·
`Fig2e_Avidity/` (2e, uses `../Extras/Effective_Molarity/`).

**Salt correction (SantaLucia).** The NN-predicted ΔG uses a monovalent-equivalent salt model:
`[Na⁺]_eq = 3.3·√[Mg²⁺]` with `[Mg²⁺] = 11.3 mM`, and `ΔG_salt = ΔG_1M − 0.114·n_steps·ln[Na⁺]_eq`,
where `n_steps = n_duplex·(L−1)` phosphate steps (each sticky-end duplex of length `L` bp contributes
`L−1`; a binding site has `n_duplex = 2` duplexes). This per-duplex convention is applied consistently
in the notebook (`saltCorrection`), the α/EM scripts, and the per-interface table. The salt term shifts
only the sequence-predicted ΔG (α, EM, Fig 2a/2d backing data, Fig 6 features); it does not touch any
measured Tm/ΔH/ΔS, the phenomenological model, or model/figure rankings.

## Figure 3 — Bond clustering  → `Figure3_BondClustering/`
`Fig3a_Tm_vs_MeanGap_7nt_4bonds/` (3a) · `Fig3b_Tm_vs_BondSpan_OLSResiduals/` (3b, Tm vs bond span g*
with bond count/length/GC controlled — the span sub-panel of the SI S3 OLS, `make_fig3b.py`) ·
`Fig3c_Tm_vs_BondPosition_ConstantSeq/` (3c) · `Fig3d_RelativeOccupancy_vs_Distance/` (3d).

## Figure 4 — Phenomenological model  → `Figure4_PhenomenologicalModel/`
`Phenomenological_Model/` (pred-vs-exp, rank-rank) · `Theory_Scripts/` + `theory_*`, `EntropyPenalty_*`,
`Theoretical*` curves · model fit: `../TmPredictor_NewTm.ipynb`. Model comparison: `../Supplementary/TableS1_Model_Comparison_All/`.

## Figure 5 — Coarse-grained MD  → `Figure5_MD_Simulations/`
(Panels E–H binder-position data — tip-gap/anchor 4×4 matrices + per-frame swarm — are in `figure_panels/mechanism_EFGH/`.)
`simulation_notebooks/` (24 HOOMD notebooks) · `figure_panels/` (b–h) · `processed_data_npz/` (291 runs).

## Figure 6 — Machine learning  → `Figure6_MachineLearning/`
`Fig6b_XGBoost_StabilityScatter/` (6b, Exp vs pred Tm) · `Fig6c_XGBoost_RankRank/` (6c, rank–rank) ·
`Fig6d_ModelComparison_SpearmanRho/` (6d). Features: `../Tile_Tile_Master.ipynb`.
(SHAP is **SI Fig S6**, in `../Supplementary/S6_SHAP_Values/`.)

## Supplementary (11/1 numbering)  → `Supplementary/`
S2 OLS avidity/α `S2_OLS_alpha_geometry/` (+ `Alpha_Geometry_Controlled/` intermediate) · S3 OLS Tm `S3_OLS_Tm_geometry/` ·
S4 ΔH/ΔS vs mean gap `S4_HS_vs_MeanGap_7nt_4bonds/` · S5 MD (Figure5) ·
S6 SHAP `S6_SHAP_Values/` · S7 oxDNA tile RMSF — analysis not included (needs the oxDNA trajectory, kept outside the repo) · S8 `S8_BondPosition_Melting/` ·
S9 `S9_BondLengthCount_Steepness/` · S10 OLS dG `S10_OLS_dG_geometry/` · S11 `S11_BondLength_Entropy/` ·
S12 `S12_PredictedVsMeasured_Tm/` · Tested-models **Table S1** `TableS1_Model_Comparison_All/`.
Regenerators: `SI_Figures_regenerate.py`, `add_ols_and_meangap.py`.

## Extras (not in the paper or SI)  → `../Extras/`
Panels/analyses kept for reference but not in the manuscript: `PredicteddG_vs_Tm/`,
`NormalizedBindingStrength/`, `BondPosition_Entropy/`, `HS_vs_Distance_7nt_4bonds/`,
`CooperativeUnit_vs_BondCount/`, `Effective_Molarity/` (feeds Fig 2e avidity),
`ReproducibleTm_UnphysicalHS/`, `Cooperativity_Tests_source/`.

## Reproducing each panel
`bash run_all.sh` regenerates **everything from scratch, self-contained** (verified: 21/21 pass): it runs
`Tile_Tile_MasterRegenerate.ipynb` (→ `Filtered_DF_*` tables), then `Tile_Tile_Master.ipynb` (→ the
notebook panels Fig 2a/2b/2d, 3a/3c/3d, 6b/6c/6d, written directly into their `FigureN/` folders — no
copying), then the standalone scripts below. Individually:

| panel | how to generate |
|---|---|
| 2a | `python3 Figure2_NegativeCooperativity/Fig2a_NNdG_ExpdG_saltfixed/regen_nndg_vs_exp.py` |
| 2b | filter `Filtered_DF_Tm.csv` by `NucleotideBonds` & `Num_Bonds` (7nt 4/5, 5nt 9/10) |
| 2c | `python3 Figure2_NegativeCooperativity/Fig2c_Alpha_PerBond_Efficiency/compute_alpha.py` |
| 2d | `python3 Figure2_NegativeCooperativity/Fig2d_ddG/regen_fig2d_ddg.py` (backing data); panel `test2_measured_only.svg` is measured-only |
| 2e | `python3 Extras/Effective_Molarity/effective_molarity.py` then `Figure2…/Fig2e_Avidity/plot_avidity.py` |
| 3a/3c/3d | `Tile_Tile_Master.ipynb` → `Figure3_BondClustering/<panel>/` |
| 3b | `python3 Figure3_BondClustering/Fig3b_Tm_vs_BondSpan_OLSResiduals/make_fig3b.py` (after `add_ols_and_meangap.py`) |
| 3 ΔH/ΔS vs distance | `python3 Extras/HS_vs_Distance_7nt_4bonds/hs_vs_distance.py` |
| 4b | `TmPredictor_NewTm.ipynb` → `Figure4…/Phenomenological_Model/` |
| 4c–4f | `python3 Figure4_PhenomenologicalModel/Theory_Scripts/*.py` |
| 5 | `Figure5_MD_Simulations/simulation_notebooks/*.ipynb` (plots in `figure_panels/`) |
| 6b/6c/6d | `Tile_Tile_Master.ipynb` (XGBoost cells) → `Figure6…/` |
| 6e | `python3 Supplementary/S6_SHAP_Values/recompute_shap_fulltrain.py` |
| S2 | `python3 Supplementary/Alpha_Geometry_Controlled/alpha_geometry_ols.py`, `…/add_ols_and_meangap.py` |
| construct table | `python3 Experimental_Trends/build_construct_medians.py` (Tm/dG from `Filtered_DF`; run before the OLS below) |
| S3, S10 | `python3 Supplementary/add_ols_and_meangap.py` (→ S3_OLS_Tm_geometry, S10_OLS_dG_geometry) |
| S4 | `python3 Extras/HS_vs_Distance_7nt_4bonds/hs_vs_distance.py`, `…/add_ols_and_meangap.py` |
| S6 | `python3 Supplementary/S6_SHAP_Values/recompute_shap_fulltrain.py` |
| S7 | oxDNA tile RMSF — not included (external oxDNA trajectory) |
| S8, S9, S11, S12 | `python3 Supplementary/SI_Figures_regenerate.py` |
| tested-models table | `python3 Supplementary/TableS1_Model_Comparison_All/models_and_pdf.py` then `pdf_report.py` |

Install deps with `pip install -r requirements.txt`.
