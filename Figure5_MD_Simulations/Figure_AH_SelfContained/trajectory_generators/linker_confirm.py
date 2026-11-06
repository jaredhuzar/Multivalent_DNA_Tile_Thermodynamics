"""Generates the linker-stiffness confirmation data behind panel D (binder distance).
Rigid scaffold (rod_bend=80), 2-bond, spacing=4, chain_bend in {4,10,20,40}. Logs, per
production sample: n_bound, every force's potential energy, the bridge tip-tip distances,
the bridge-tip positions, and the two rod COMs. Writes ../confirm_results/linker_chain{chain}.npz;
panel D is the distribution of the unbound tip-tip (binder) distance vs chain_bend.
Reconstructed from the project history.   Run:  python3 linker_confirm.py
"""
import os, time, collections, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.abspath(os.path.join(HERE, "..", "confirm_results")); os.makedirs(OUT, exist_ok=True)
exec(open(os.path.join(HERE, "dimer_defs.py")).read(), globals())
import hoomd, matplotlib; matplotlib.use("Agg")

GEOM = [11 - 4 // 2, 11 + 4 // 2]      # sp4, 2-bond
NUMBONDS, BPC, BDIAM, GRAFT_SP, TIP_GAP = 4, 3, 1, 1, 6
ROD_BEND, ROD_BOND, CHAIN_BOND = 80, 60, 50
D0 = 8; GAMMA = 5.0; R_BIND = 1.4
TEMPS = [1.4, 1.6, 1.8, 2.0]; BURN = 150000; NSAMP = 300; SINT = 3000
CHAINS = [4, 10, 20, 40]; BOX = (30.0, 28.0, 28.0)

def rod_groups(snap):
    bt = list(snap.bonds.types); rid = bt.index('rod'); grp = snap.bonds.group; tid = snap.bonds.typeid
    adj = collections.defaultdict(set); nodes = set()
    for k in range(len(grp)):
        if tid[k] == rid:
            i, j = int(grp[k][0]), int(grp[k][1]); adj[i].add(j); adj[j].add(i); nodes |= {i, j}
    seen = set(); comps = []
    for n in nodes:
        if n in seen: continue
        st = [n]; c = []
        while st:
            x = st.pop()
            if x in seen: continue
            seen.add(x); c.append(x); st.extend(adj[x] - seen)
        comps.append(sorted(c))
    comps.sort(key=len, reverse=True); return comps[0], comps[1]

def run_design(chain):
    cp = [chain, CHAIN_BOND]; rp = [ROD_BEND, ROD_BOND]
    snap, bondsArr, anchor_tags, chain1_tags, signs = layoutSimulation(
        NUMBONDS, BPC, BDIAM, GRAFT_SP, TIP_GAP, GEOM, GEOM)
    snap.configuration.box = hoomd.Box(Lx=BOX[0], Ly=BOX[1], Lz=BOX[2])
    bridge_pairs = getBondPairs(bondsArr[0], bondsArr[1])
    tips = sorted(set([i for p in bridge_pairs for i in p]))
    sim = hoomd.Simulation(device=cpu, seed=1); sim.create_state_from_snapshot(snap)
    g0, g1 = rod_groups(sim.state.get_snapshot())
    U, F, rmn, rmx = morse_Curve_Table(d0_ratio=D0)
    f = forces_MD([U, F, rmn, rmx], cp, rp, anchor_tags, chain1_tags, GRAFT_SP, signs); sim = relax_MD(sim, list(f))
    U1, F1, r1, r2 = morse_Curve_Table_assosciation(d0_ratio=100)
    f = forces_MD([U1, F1, r1, r2], cp, rp, anchor_tags, chain1_tags, GRAFT_SP, signs)
    sim = Equilibrate_MD(sim, list(f), bridge_pairs)
    bound = sim.state.get_snapshot()
    out = {"chain": chain, "temps": TEMPS, "bridge_pairs": np.array(bridge_pairs), "tips": np.array(tips)}
    for T in TEMPS:
        sim = hoomd.Simulation(device=cpu, seed=7); sim.create_state_from_snapshot(bound)
        U, F, rmn, rmx = morse_Curve_Table(d0_ratio=D0)
        lj, angle_force, fene, rod_bond, bridge_bond, angle, dih, cust = forces_MD(
            [U, F, rmn, rmx], cp, rp, anchor_tags, chain1_tags, GRAFT_SP, signs)
        integ = hoomd.md.Integrator(dt=dt * 0.1)
        integ.forces = [lj, angle_force, fene, rod_bond, bridge_bond, angle, cust]
        integ.methods = [hoomd.md.methods.Brownian(filter=hoomd.filter.All(), kT=T, default_gamma=GAMMA)]
        sim.operations.integrator = integ; sim.run(BURN)
        named = dict(bridge=bridge_bond, bend=angle_force, lj=lj, fene=fene, rodbond=rod_bond, graft=angle, zalign=cust)
        E = {k: np.empty(NSAMP) for k in named}; nb = np.empty(NSAMP, int)
        dist = np.empty((NSAMP, len(bridge_pairs))); tippos = np.empty((NSAMP, len(tips), 3)); rcom = np.empty((NSAMP, 2, 3))
        inv = {t: k for k, t in enumerate(tips)}
        for s in range(NSAMP):
            sim.run(SINT)
            n, _, _ = bridge_state(sim, bridge_pairs=bridge_pairs, r_bind=R_BIND); nb[s] = n
            for k, fo in named.items():
                try: E[k][s] = float(fo.energy)
                except Exception: E[k][s] = np.nan
            ss = sim.state.get_snapshot(); pos = np.array(ss.particles.position)
            dist[s] = [np.linalg.norm(pos[i] - pos[j]) for (i, j) in bridge_pairs]
            for t in tips: tippos[s, inv[t]] = pos[t]
            rcom[s, 0] = pos[g0].mean(0); rcom[s, 1] = pos[g1].mean(0)
        out[f"T{T}_n"] = nb; out[f"T{T}_dist"] = dist; out[f"T{T}_tippos"] = tippos; out[f"T{T}_rcom"] = rcom
        for k in named: out[f"T{T}_E_{k}"] = E[k]
        print(f"  chain={chain} T={T}: <n>={nb.mean():.2f} ({time.time():.0f})", flush=True)
    np.savez(os.path.join(OUT, f"linker_chain{chain}.npz"), **out)
    print(f"saved linker_chain{chain}.npz", flush=True)

if __name__ == "__main__":
    print(f"=== linker-stiffness confirmation (rod=80, 2-bond, sp=4) -> {OUT} ===", flush=True)
    for ch in CHAINS:
        print(f"--- chain_bend={ch} ---", flush=True); run_design(ch)
    print("DONE", flush=True)
