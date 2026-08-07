
import os, glob
import numpy as np
import pyemu


def g07_evaluate(X):
    # X is (n_run, 10).  Returns a dict of observation name -> (n_run,) array.
    #
    # obj/g1..g8 are the standard g07 problem.  obj2 is an ADDED second objective - the squared
    # distance from the origin - which trades off against obj and gives pestpp-mou a genuine
    # two-objective front to sort.  With one objective mou degenerates to a constrained EA and
    # never exercises dominance/crowding/archive, which is most of what could differ by
    # platform.  The other four tools carry obj2 as a zero-weight observation, unaffected.
    x = [X[:, i] for i in range(10)]
    obj = (x[0]**2 + x[1]**2 + x[0]*x[1] - 14*x[0] - 16*x[1] + (x[2] - 10)**2
           + 4*(x[3] - 5)**2 + (x[4] - 3)**2 + 2*(x[5] - 1)**2 + 5*x[6]**2
           + 7*(x[7] - 11)**2 + 2*(x[8] - 10)**2 + (x[9] - 7)**2 + 45)
    obj2 = sum(xi**2 for xi in x)
    out = {"obj": obj, "obj2": obj2}
    out["g1"] = -105 + 4*x[0] + 5*x[1] - 3*x[6] + 9*x[7]
    out["g2"] = 10*x[0] - 8*x[1] - 17*x[6] + 2*x[7]
    out["g3"] = -8*x[0] + 2*x[1] + 5*x[8] - 2*x[9] - 12
    out["g4"] = 3*(x[0]-2)**2 + 4*(x[1]-3)**2 + 2*x[2]**2 - 7*x[3] - 120
    out["g5"] = 5*x[0]**2 + 8*x[1] + (x[2]-6)**2 - 2*x[3] - 40
    out["g6"] = x[0]**2 + 2*(x[1]-2)**2 - 2*x[0]*x[1] + 14*x[4] - 6*x[5]
    out["g7"] = 0.5*(x[0] - 8)**2 + 2*(x[1]-4)**2 + 3*x[4]**2 - x[5] - 30
    out["g8"] = -3*x[0] + 6*x[1] + 12*(x[8] - 8)**2 - 7*x[9]
    return out


def forward_g07_ensemble(rns_file=None):
    if rns_file is None:
        best, nbest = None, -1
        for c in [f for f in glob.glob("*.rns") if "fail" not in f.lower()]:
            h, _, _ = pyemu.helpers.RunStor.file_info(c)
            if h["n_runs"] > nbest:
                best, nbest = c, h["n_runs"]
        rns_file = best
    header, par_names, obs_names = pyemu.helpers.RunStor.file_info(rns_file)
    npar, nobs, nrun = len(par_names), len(obs_names), int(header["n_runs"])
    run_start, run_size = int(header["run_start"]), int(header["run_size"])
    par_off = 1 + 1001 + 8            # r_status(int8) + info_txt(1001) + info_val(float64)
    obs_off = par_off + npar * 8

    par = np.empty((nrun, npar))
    with open(rns_file, "rb") as f:
        for i in range(nrun):
            f.seek(run_start + i * run_size + par_off)
            par[i, :] = np.fromfile(f, dtype=np.float64, count=npar)

    pmap = {n.lower(): j for j, n in enumerate(par_names)}
    X = np.column_stack([par[:, pmap["x{0}".format(k)]] for k in range(1, 11)])
    sim = g07_evaluate(X)

    omap = {n.lower(): j for j, n in enumerate(obs_names)}
    obs = np.zeros((nrun, nobs))
    for name, vals in sim.items():
        j = omap.get(name)
        if j is not None:
            obs[:, j] = vals

    # obs values, run_status=1 (completed), and the trailing buffer byte
    with open(rns_file, "r+b") as f:
        for i in range(nrun):
            base = run_start + i * run_size
            f.seek(base); f.write(np.int8(1).tobytes())
            f.seek(base + obs_off); f.write(obs[i].astype(np.float64).tobytes())
            f.seek(base + obs_off + nobs * 8); f.write(np.int8(0).tobytes())
    print("forward_g07_ensemble: solved {0} runs from {1}".format(nrun, rns_file))


def forward_g07_single():
    # the CLASSIC path: read par.dat written from the template file, write obs.dat for the
    # instruction file.  Some tools run the model DIRECTLY rather than through the run manager -
    # pestpp-opt does it for "running the model once with optimal decision variables" - and that
    # run reads and writes files, not the run-storage.  Without this the direct runs fail with
    # "output file './obs.dat' not found" and the tool carries on with a missing result.
    if not os.path.exists("par.dat"):
        return
    vals = {}
    for line in open("par.dat"):
        t = line.strip().split()
        if len(t) >= 2:
            try:
                vals[t[0].strip().lower()] = float(t[1])
            except ValueError:
                pass
    if len(vals) < 10:
        return
    X = np.array([[vals["x{0}".format(k)] for k in range(1, 11)]])
    sim = g07_evaluate(X)
    with open("obs.dat", "w") as f:
        for name in ["obj", "obj2"] + ["g{0}".format(i) for i in range(1, 9)]:
            f.write("{0}  {1:20.12E}\n".format(name, float(sim[name][0])))


if __name__ == "__main__":
    # BOTH, unconditionally: the batch when there is one, and the single-run files when the
    # tool wrote par.dat.  One script that serves the external run manager and a direct run
    # means no tool needs a different model command.
    forward_g07_single()
    forward_g07_ensemble()
