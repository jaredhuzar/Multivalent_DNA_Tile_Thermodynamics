# Figure label mapping — repo ↔ final manuscript / SI

Cross-referenced from **PNAS_main_manuscript.docx** and **PNAS_SI_Appendix.docx** (SI figs renumbered to S1–S15).
The repo folder/label names predate the final numbering, so many differ. Legend: ✓ already matches · ⚠ **mismatch, needs relabel** · ? ambiguous/verify.

## Main figures

| repo path / source-data label | content | final manuscript panel | |
|---|---|---|---|
| `Figure1_Platform_AFM/` AFM images | tile dimers (AFM) | **Fig 1B** | ✓ |
| `Figure1_Platform_AFM/Melting Curves.svg` → `Fig1C_MeltingCurves.csv` | fitted melting curves | **Fig 1C** | ✓ |
| `Fig2a_NNdG_ExpdG_saltfixed` | NN-predicted vs measured ΔG | **Fig 2A, 2B** | ✓ |
| `Fig2b_Tm_vs_BindingSites_x_Length` | Tm by sites × bond length | **Fig 2C** | ⚠ 2b→2C |
| `Fig2c_Alpha_PerBond_Efficiency` | cooperativity α | **Fig 2D** | ⚠ 2c→2D |
| `Fig2d_ddG` | ΔΔG of adding one site | **Fig 2E** | ⚠ 2d→2E |
| `Fig2e_Avidity` | avidity enhancement | **Fig 2F** | ⚠ 2e→2F |
| `Fig3a_Tm_vs_MeanGap_7nt_4bonds` | Tm vs mean gap (7nt, 4 sites) | **Fig 3A** (+3C) | ✓/? |
| `Fig3c_Tm_vs_BondPosition_ConstantSeq` | Tm vs bond position | **Fig 3D** (+3B schematic) | ⚠ 3c→3D |
| `Fig3d_RelativeOccupancy_vs_Distance` | inferred occupancy vs distance | **Fig 3E, 3G** | ⚠ 3d→3E/3G |
| `Fig3b_Tm_vs_BondSpan_OLSResiduals` | Tm vs bond span (covariates out) | **Fig 3F** | ⚠ 3b→3F |
| `Figure4…/Phenomenological_Model/pred_vs_exp` | predicted vs measured Tm (OOF) | **Fig 4B** | ✓ |
| `EntropyPenalty_vs_NumBonds_spacing4.csv` | net ΔG°37 vs bond number | **Fig 4C** | ⚠ (source-data labeled 4d) |
| `theory_vs_spacing.csv` | Tm vs spacing | **Fig 4D** | ⚠ (source-data labeled 4e) |
| `theory_vs_Nbonds.csv` | Tm vs bond number (diminishing) | **Fig 4E** (+4G plateau) | ⚠ (source-data labeled 4c) |
| `EntropyPenalty_vs_Spacing_4bonds.csv` | ΔG decomposition vs spacing | **Fig 4F** (map) | ? verify (4F is a 2-D map) |
| `plots_robust_v2/b` → `Fig5b_Tm_vs_BondSpacing` | Tm vs spacing, rigid/flex | **Fig 5B** | ✓ |
| `plots_robust_v2/d` → `Fig5c_Tm_vs_LinkerRigidity` | Tm vs linker stiffness | **Fig 5C** | ✓ |
| `plots_robust_v2/c` → `Fig5e_Tm_vs_ScaffoldRigidity` | Tm vs scaffold rigidity | **Fig 5D** | ⚠ source-data 5e→5D |
| `mechanism panelC` → `Fig5d_BinderDistance_byLinker` | tip–tip distance by linker | **Fig 5E** | ⚠ source-data 5d→5E |
| `registration` → `Fig5f_ExtensionSpread_vs_RodBending` | tether spread vs scaffold | **Fig 5F** | ✓ |
| `heatmap_grid` → `Fig5g/Fig5h` | tip-gap mean / spread heatmaps | **Fig 5G, 5H** | ✓ |
| `Fig6b_XGBoost_StabilityScatter` | XGBoost predicted vs measured | **Fig 6B** | ✓ |
| `Fig6d_ModelComparison_SpearmanRho` | per-model Spearman ρ | **Fig 6C** | ⚠ 6d→6C |
| `S6_SHAP_Values` (main copy) | SHAP per site | **Fig 6D** (+ SI S12) | ⚠ |
| edge Mann–Whitney panel | Tm by edge occupancy | **Fig 6E** | ✓ |
| `Fig6c_XGBoost_RankRank` | rank–rank | *no main panel* — SI or dropped | ? |

## SI figures (final SI number ← repo folder)

| final SI | content | repo folder | |
|---|---|---|---|
| **S1** | entropic penalty (ΔS°) by bond length | `S11_BondLength_Entropy` | ⚠ S11→S1 |
| S2 | scadnano design (all sequences) | *(notebook/design)* | — |
| S3 | interface regions, every base | *(design)* | — |
| S4 | fluorophore–quencher reporter | *(design)* | — |
| **S5** | width dx of the melting transition | `S9_BondLengthCount_Steepness` | ⚠ S9→S5 ? |
| S6 | AFM of tile dimers | `Figure1_Platform_AFM` (images) | — |
| **S7** | NN-predicted vs measured Tm | `S12_PredictedVsMeasured_Tm` | ⚠ S12→S7 |
| **S8** | OLS of α (added-variable) | `S2_OLS_alpha_geometry` (+`Alpha_Geometry_Controlled`) | ⚠ S2→S8 |
| **S9** | OLS of Tm | `S3_OLS_Tm_geometry` | ⚠ S3→S9 |
| **S10** | ΔH°/ΔS° vs mean gap | `S4_HS_vs_MeanGap_7nt_4bonds` | ⚠ S4→S10 |
| S11 | simulated stability vs strength/count/T | `Figure5_MD_Simulations` | — |
| **S12** | SHAP values (XGBoost) | `S6_SHAP_Values` | ⚠ S6→S12 |
| **S13** | Tm by occupied site | `S8_BondPosition_Melting` | ⚠ S8→S13 |
| S14 | RMSF of simulated tile monomer | *(oxDNA, external)* | — |
| **S15** | OLS of ΔG°37 | `S10_OLS_dG_geometry` | ⚠ S10→S15 |
| — | *(no final-SI match)* | `S3_BondPosition_Entropy` | ? likely deprecated |

Tables: repo `TableS1_Model_Comparison_All` → **Table S1** ✓; fitted params → **Table S4**; CG parameters → **Table S5**.
