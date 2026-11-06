"""Generates the chain x rod GRID trajectories behind panels G and H.
2-bond, spacing=4. For every (chain_bend, rod_bend) in {4,10,20,40} x {4,20,80,240}, at a few
temperatures, run a melt-from-bound trajectory and write a GSD at the sampling cadence to
../grid_results/gsd_grid_c{chain}_r{rod}_T{T}.gsd. The heatmap cell then computes tip-gap
mean (G) and spread (H) post-hoc. Reconstructed from the project history.
Run:  python3 grid_confirm.py          (uses a small process Pool; env GRID_NW sets workers)
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1"); os.environ.setdefault("MPLBACKEND", "Agg")
import time, numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.abspath(os.path.join(HERE, "..", "grid_results")); os.makedirs(OUT, exist_ok=True)

GEOM = [9, 13]                       # sp4, 2-bond
NUMBONDS, BPC, BDIAM, GRAFT_SP, TIP_GAP = 4, 3, 1, 1, 6
CHAIN_BOND, ROD_BOND = 50, 60
D0, GAMMA = 8, 5.0
TEMPS = [float(x) for x in os.environ.get("GRID_TEMPS", "1.5,2.0").split(",")]
BURN = int(os.environ.get("GRID_BURN", "150000"))
NSAMP = int(os.environ.get("GRID_NSAMP", "300"))
SINT = int(os.environ.get("GRID_SINT", "3000"))
CHAINS = [int(x) for x in os.environ.get("GRID_CHAINS", "4,10,20,40").split(",")]
RODS = [int(x) for x in os.environ.get("GRID_RODS", "4,20,80,240").split(",")]
BOX = (30.0, 28.0, 28.0)
NW = int(os.environ.get("GRID_NW", "4"))

def run_cell(args):
    chain, rod = args
    import hoomd
    exec(open(os.path.join(HERE, "dimer_defs.py")).read(), globals())
    cp = [chain, CHAIN_BOND]; rp = [rod, ROD_BOND]
    todo = [T for T in TEMPS if not os.path.exists(os.path.join(OUT, f"gsd_grid_c{chain}_r{rod}_T{T:g}.gsd"))]
    if not todo:
        return f"skip c{chain} r{rod} (done)"
    snap, bondsArr, anchor_tags, chain1_tags, signs = layoutSimulation(
        NUMBONDS, BPC, BDIAM, GRAFT_SP, TIP_GAP, GEOM, GEOM)
    snap.configuration.box = hoomd.Box(Lx=BOX[0], Ly=BOX[1], Lz=BOX[2])
    bridge_pairs = getBondPairs(bondsArr[0], bondsArr[1])
    sim = hoomd.Simulation(device=cpu, seed=1); sim.create_state_from_snapshot(snap)
    U, F, rmn, rmx = morse_Curve_Table(d0_ratio=D0)
    f = forces_MD([U, F, rmn, rmx], cp, rp, anchor_tags, chain1_tags, GRAFT_SP, signs)
    sim = relax_MD(sim, list(f))
    U1, F1, r1, r2 = morse_Curve_Table_assosciation(d0_ratio=100)
    f = forces_MD([U1, F1, r1, r2], cp, rp, anchor_tags, chain1_tags, GRAFT_SP, signs)
    sim = Equilibrate_MD(sim, list(f), bridge_pairs)
    bound = sim.state.get_snapshot()
    for T in todo:
        gpath = os.path.join(OUT, f"gsd_grid_c{chain}_r{rod}_T{T:g}.gsd")
        sim = hoomd.Simulation(device=cpu, seed=7); sim.create_state_from_snapshot(bound)
        U, F, rmn, rmx = morse_Curve_Table(d0_ratio=D0)
        lj, angle_force, fene, rod_bond, bridge_bond, angle, dih, cust = forces_MD(
            [U, F, rmn, rmx], cp, rp, anchor_tags, chain1_tags, GRAFT_SP, signs)
        integ = hoomd.md.Integrator(dt=dt * 0.1)
        integ.forces = [lj, angle_force, fene, rod_bond, bridge_bond, angle, cust]
        integ.methods = [hoomd.md.methods.Brownian(filter=hoomd.filter.All(), kT=T, default_gamma=GAMMA)]
        sim.operations.integrator = integ
        sim.run(BURN)
        w = hoomd.write.GSD(filename=gpath, mode="wb",
                            trigger=hoomd.trigger.Periodic(period=SINT, phase=sim.timestep % SINT),
                            filter=hoomd.filter.All(), dynamic=["property"])
        sim.operations.writers.append(w)
        sim.run(SINT * NSAMP)
        try: w.flush()
        except Exception: pass
        sim.operations.writers.remove(w)
    return f"done c{chain} r{rod} ({time.strftime('%H:%M:%S')})"

if __name__ == "__main__":
    tasks = [(c, r) for c in CHAINS for r in RODS]
    print(f"=== chain x rod grid: {len(tasks)} designs x {len(TEMPS)} temps, pool={NW} -> {OUT} ===", flush=True)
    t0 = time.time()
    with Pool(processes=NW) as pool:
        for msg in pool.imap_unordered(run_cell, tasks):
            print(f"  {msg}   [{int(time.time()-t0)}s]", flush=True)
    print(f"=== GRID DONE {time.strftime('%Y-%m-%d %H:%M:%S')} ({int(time.time()-t0)}s) ===", flush=True)
