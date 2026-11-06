# Trajectory generators

These scripts/notebooks produce the raw simulation data read by `../Analysis_Tm_trends.ipynb`
to make figure panels A–H. All use HOOMD-blue (Brownian dynamics, γ=5, dt·0.1, Morse bonds
D0=8), box 30×28×28, 2-bond interfaces unless noted. `dimer_defs.py` holds the shared
definitions (layout, forces, relax/equilibrate, bridge_state) extracted from the base notebook.

## What generates what

| Generator | Output (in the figure folder) | Feeds panel(s) |
|---|---|---|
| `run_batch.py` + the interface notebooks (`Base.ipynb`, `RigidScaffold.ipynb`, `*Chain.ipynb`, `*3bonds.ipynb`, …) | `std_results_v2/` (npz) + `gsd_results/` (GSD trajectories) | **B** (3-bond spacing), **C** (linker rigidity), **E** (scaffold rigidity), **F** (sp4 registration) |
| `linker_confirm.py` | `confirm_results/linker_chain{4,10,20,40}.npz` | **D** (binder distance vs linker stiffness) |
| `grid_confirm.py` | `grid_results/gsd_grid_c{chain}_r{rod}_T{T}.gsd` | **G, H** (chain×rod heatmaps) |

Panel **A** is a hand-drawn schematic (no generator).

## Parameters (what each sweep varies; everything else fixed)
- **run_batch / interface notebooks** — each notebook fixes its own `chain_bend` (linker) and
  `rod_bend` (scaffold); the injected `standardized_melt` cell runs the uniform ladder
  (kT 0.5–2.0 step 0.25), 1.5M burn-in, 300 samples, 3 seeds, and writes a GSD per temperature.
  - B: `FlexScaffold3bonds` (rod 4) vs `RigidScaffold3bonds` (rod 80), 3-bond, chain 4, D0 8, spacings swept
  - C: `RigidScaffold_*Chain` (rod 80, chain 3→40), D0 8, spacings 2 & 6
  - E/F: `Base`/`ModerateScaffold`/`RigidScaffold`/`VeryRigidScaffold` (rod 4/20/80/240), chain 4, D0 8
- **linker_confirm.py** — rod 80 fixed, chain ∈ {4,10,20,40}, sp4, kT {1.4,1.6,1.8,2.0}.
- **grid_confirm.py** — chain ∈ {4,10,20,40} × rod ∈ {4,20,80,240}, sp4, kT {1.5,2.0}.

## Running (regenerates the data dirs the analysis notebook reads)
```
python3 linker_confirm.py        # -> ../confirm_results/   (panel D)
python3 grid_confirm.py          # -> ../grid_results/      (panels G,H); GRID_NW sets workers
# main batch (panels B,C,E,F): run the dispatcher from the FIGURE ROOT so paths resolve:
cd ..  &&  python3 trajectory_generators/run_batch.py     # -> std_results_v2/, gsd_results/
```
Note: `grid_confirm.py` and `linker_confirm.py` were reconstructed from the project history
(their originals were deleted); their geometry was verified to match the shipped data exactly.
`light_rerun.py` (the scaffold-rigidity force-decomposition run that made `lightrun_rod*.npz`)
is **not** included — it feeds the non-figure `mechanism.svg`, not panels A–H.
