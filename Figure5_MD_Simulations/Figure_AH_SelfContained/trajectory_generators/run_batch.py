#!/usr/bin/env python3
"""Parallel batch executor for the patched Base notebooks.
Each notebook runs in its OWN working directory (batch_run/<stem>/) so same-named
scratch files (MorseTestEquil.gsd, Summary*.csv, ...) from parallel jobs cannot
collide. All per-sample results go to one shared std_results/ via STD_OUTDIR.
Resumable: skips a notebook whose executed copy already exists, and (via the
injected function) skips any replicate whose .npz already exists."""
import subprocess, os, glob, time, shutil
from concurrent.futures import ThreadPoolExecutor

PROJ = "/Users/jaredhuzar/hoomd-blue/build/Claude"
os.chdir(PROJ)
# self-protect against idle/system sleep on macOS (releases when this process exits)
import platform
if platform.system() == "Darwin":
    try:
        subprocess.Popen(["caffeinate", "-is", "-w", str(os.getpid())])
    except Exception:
        pass
RESULTS = os.environ.get("STD_OUTDIR", os.path.join(PROJ, "std_results_v2"))
GSDDIR = os.environ.get("STD_GSDDIR", os.path.join(PROJ, "gsd_results"))
BATCHDIR = os.environ.get("BATCH_RUNDIR", "batch_run_v2")
LOGDIR = os.environ.get("BATCH_LOGDIR", "batch_logs_v2")
for d in (RESULTS, GSDDIR, BATCHDIR, LOGDIR):
    os.makedirs(d, exist_ok=True)
files = [os.path.basename(f) for f in sorted(glob.glob("*.ipynb"))
         if "_Energy-Copy1" not in f and "_Edge_LongRun" not in f and ".exec." not in f and os.path.basename(f)!="Analysis_Tm_trends.ipynb"]
BASEENV = dict(os.environ, MPLBACKEND="Agg", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
               OPENBLAS_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1",
               STD_OUTDIR=RESULTS, STD_GSDDIR=GSDDIR)
NPAR = int(os.environ.get("NPAR", "12"))

def run(f):
    stem = f[:-6]
    rundir = os.path.join(PROJ, BATCHDIR, stem)
    os.makedirs(rundir, exist_ok=True)
    exec_copy = os.path.join(rundir, f"{f}.exec.ipynb")
    if os.path.exists(exec_copy):                      # resume: whole notebook done
        print(f"SKIP  {f} (already executed)", flush=True)
        return 0
    shutil.copy(os.path.join(PROJ, f), rundir)         # isolate source in its own dir
    t0 = time.time()
    print(f"START {f} {time.strftime('%H:%M:%S')}", flush=True)
    with open(os.path.join(PROJ, LOGDIR, f"{f}.log"), "w") as log:
        rc = subprocess.call(
            ["jupyter", "nbconvert", "--to", "notebook", "--execute", "--allow-errors",
             "--ExecutePreprocessor.timeout=0", "--output", f"{f}.exec.ipynb", f],
            cwd=rundir, stdout=log, stderr=subprocess.STDOUT, env=BASEENV)
    npz = len(glob.glob(os.path.join(RESULTS, f"std_{stem}_*.npz")))
    print(f"DONE  {f} rc={rc} {int(time.time()-t0)}s npz={npz}", flush=True)
    return rc

print(f"launching {len(files)} notebooks, {NPAR} parallel, isolated dirs, "
      f"results -> {RESULTS}", flush=True)
with ThreadPoolExecutor(max_workers=NPAR) as ex:
    rcs = list(ex.map(run, files))
print(f"=== BATCH COMPLETE {time.strftime('%Y-%m-%d %H:%M:%S')} ===", flush=True)
print(f"nonzero rc: {sum(1 for r in rcs if r)}/{len(files)} | "
      f"total npz: {len(glob.glob(os.path.join(RESULTS, '*.npz')))}", flush=True)
