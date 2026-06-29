import os, sys, glob, re
import numpy as np
import pyemu

def evolve_lorenz(X0, n_steps, dt_step=0.01, F=8.0):
    """vectorized fixed-step RK4 evolution of a lorenz96 ensemble.  X0 is (M,N) (or (N,)); returns a
    (n_steps, M, N) trajectory where out[k] is X0 advanced (k+1) steps of dt_step.  fixed-step (not
    adaptive) so each realization integrates identically regardless of the rest of the batch - the
    batched ensemble solve then reproduces a single-realization solve exactly (clean twin)."""
    import numpy as np
    def fv(X):
        return (np.roll(X, -1, axis=-1) - np.roll(X, 2, axis=-1)) * np.roll(X, 1, axis=-1) - X + F
    x = np.atleast_2d(np.asarray(X0, dtype=float))
    out = np.empty((n_steps,) + x.shape)
    for k in range(n_steps):
        k1 = fv(x); k2 = fv(x + 0.5*dt_step*k1); k3 = fv(x + 0.5*dt_step*k2); k4 = fv(x + dt_step*k3)
        x = x + (dt_step/6.0)*(k1 + 2*k2 + 2*k3 + k4)
        out[k] = x
    return out


def forward_run_ensemble(rns_file=None, F=8.0):
    """external ('/e') run-manager forward run: read the FULL ensemble of parameter sets from the
    pest++ binary run-storage file, solve lorenz96 for all realizations in one vectorized pass, and
    pack the simulated states back into the run-storage file as observations.  uses pyemu's RunStor
    to parse the file layout and writes results positionally (the time-0 obs and the state-IC pars
    share names, so name-based access can't be used)."""
    import os, re, glob
    import numpy as np
    import pyemu

    # locate the run-storage file (case.rns), skipping the failed-runs file
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
    par_off = 1 + 1001 + 8          # r_status(int8) + info_txt(1001) + info_val(float64)
    obs_off = par_off + npar * 8    # obs block follows the par block (matches pyemu RunStor.update)

    # read every parameter set (positional read avoids the par/obs name collision at time 0)
    par = np.empty((nrun, npar))
    with open(rns_file, "rb") as f:
        for i in range(nrun):
            f.seek(run_start + i * run_size + par_off)
            par[i] = np.fromfile(f, dtype=np.float64, count=npar)
    pmap = {n: j for j, n in enumerate(par_names)}
    omap = {n: j for j, n in enumerate(obs_names)}

    # initial states are the par names x<i>_00.000, ordered by state index
    state_cols = sorted([n for n in par_names if re.match(r"x\d+_0*0\.000$", n)],
                        key=lambda n: int(n.split("_")[0][1:]))
    snames = [c.split("_")[0] for c in state_cols]
    M, N = nrun, len(state_cols)
    X0 = par[:, [pmap[c] for c in state_cols]]
    delt = float(par[0, pmap["delt"]]); t_start = float(par[0, pmap["t_start"]]); t_end = float(par[0, pmap["t_end"]])
    times = np.arange(0.0, t_end - t_start, delt)

    # vectorized, batched, fixed-step evolution of the whole (M,N) ensemble.  NOTE: the integration
    # step is 0.01 (matching evolve_lorenz/the original rand_evolve default) - 'delt' is ONLY the obs
    # time-label spacing, not the integration step.  out[k] = state after (k+1) steps.
    traj = evolve_lorenz(X0, len(times), dt_step=0.01, F=F)
    obs = np.zeros((nrun, nobs))
    for k, otime in enumerate(times):
        for i, sn in enumerate(snames):
            j = omap.get("{0}_{1:06.3f}".format(sn, otime))
            if j is not None:
                obs[:, j] = traj[k, :, i]

    # pack results back into the run-storage file: obs values, run_status=1 (completed), buffer=0
    with open(rns_file, "r+b") as f:
        for i in range(nrun):
            base = run_start + i * run_size
            f.seek(base); f.write(np.int8(1).tobytes())
            f.seek(base + obs_off); f.write(obs[i].astype(np.float64).tobytes())
            f.seek(base + obs_off + nobs * 8); f.write(np.int8(0).tobytes())
    print("forward_run_ensemble: solved {0} runs from {1}".format(nrun, rns_file))

if __name__ == '__main__':
    forward_run_ensemble()
