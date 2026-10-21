"""anafield: an analytical, nonlinear, scalable benchmark for the ensemble methods.

the parameter field is a log10 property on an nrow by ncol grid drawn from an
exponential variogram.  there is no model to run - the observations are closed
form functionals of the field, so the problem goes to hundreds of thousands of
parameters and the forward run is a numpy call.  three kinds of observation,
each with its own flavour of nonlinearity:

    point   y = sinh(a * x) / a at a cell.  a is the nonlinearity dial: a -> 0
            is linear, a = 2 or 3 is strongly convex in |x|
    path    y = log10 of the harmonic mean of 10**x along a straight line of
            cells.  a soft minimum - the low cells on the line control it
    block   y = log10 of the arithmetic mean of 10**x over a window.  a soft
            maximum - the high cells in the window control it

all observations are on the log10 scale so one noise level fits them all.
the truth is a draw from the same prior, so it is reachable and the fit can
go to the noise floor.  the closed form jacobian is here too for checking
against finite differences or for pestpp-glm.

setup() writes a template dir with the control file, the forward script, the
prior ensemble, the prior covariance (dense, so only up to cov_max_pars) and
a rook lattice graph for enif.  anafield_dashboard.py lets you play with the
same inputs before you commit to a size.
"""
import os
import json
import shutil
import numpy as np
import pandas as pd
import pyemu


def make_config(nrow=20, ncol=20, n_points=20, n_paths=10, n_windows=5,
                window=4, a=2.0, seed=12345):
    """pick where the observations sit.  everything is drawn from one rng so
    the same seed gives the same layout at any field size."""
    rng = np.random.default_rng(seed)
    pts = np.column_stack([rng.integers(0, nrow, n_points),
                           rng.integers(0, ncol, n_points)])
    paths = []
    for _ in range(n_paths):
        i0, j0 = rng.integers(0, nrow), rng.integers(0, ncol)
        i1, j1 = rng.integers(0, nrow), rng.integers(0, ncol)
        paths.append(path_cells(i0, j0, i1, j1))
    window = int(max(1, min(window, nrow, ncol)))
    wins = np.column_stack([rng.integers(0, nrow - window + 1, n_windows),
                            rng.integers(0, ncol - window + 1, n_windows)])
    return {"nrow": int(nrow), "ncol": int(ncol), "a": float(a), "window": window,
            "points": pts.tolist(),
            "paths": [p.tolist() for p in paths],
            "windows": wins.tolist()}


def path_cells(i0, j0, i1, j1):
    """the cells on the straight line between two cells, each counted once"""
    n = int(max(abs(i1 - i0), abs(j1 - j0))) + 1
    ii = np.rint(np.linspace(i0, i1, n)).astype(int)
    jj = np.rint(np.linspace(j0, j1, n)).astype(int)
    cells = np.column_stack([ii, jj])
    _, keep = np.unique(cells, axis=0, return_index=True)
    return cells[np.sort(keep)]


def obs_names(cfg):
    names = ["point_{0:04d}".format(k) for k in range(len(cfg["points"]))]
    names += ["path_{0:04d}".format(k) for k in range(len(cfg["paths"]))]
    names += ["block_{0:04d}".format(k) for k in range(len(cfg["windows"]))]
    return names


def forward(x, cfg):
    """x is the log10 field (nrow, ncol).  returns the observation vector in
    the order of obs_names()"""
    x = np.asarray(x, dtype=float)
    a = cfg["a"]
    out = []
    for i, j in cfg["points"]:
        xv = x[i, j]
        out.append(xv if a == 0.0 else np.sinh(a * xv) / a)
    for cells in cfg["paths"]:
        cells = np.asarray(cells)
        xv = x[cells[:, 0], cells[:, 1]]
        # log10 of the harmonic mean of 10**x, written to stay finite
        out.append(np.log10(len(xv)) - np.log10(np.sum(10.0 ** (-xv))))
    w = cfg["window"]
    for i, j in cfg["windows"]:
        xv = x[i:i + w, j:j + w]
        out.append(np.log10(np.mean(10.0 ** xv)))
    return np.array(out)


def jacobian(x, cfg):
    """d obs / d x, dense (nobs, nrow*ncol), x flattened row major.  the path
    and block rows are the softmin / softmax weights of their cells"""
    x = np.asarray(x, dtype=float)
    nrow, ncol = x.shape
    a = cfg["a"]
    rows = []
    for i, j in cfg["points"]:
        r = np.zeros(nrow * ncol)
        r[i * ncol + j] = 1.0 if a == 0.0 else np.cosh(a * x[i, j])
        rows.append(r)
    for cells in cfg["paths"]:
        cells = np.asarray(cells)
        r = np.zeros(nrow * ncol)
        wgt = 10.0 ** (-x[cells[:, 0], cells[:, 1]])
        r[cells[:, 0] * ncol + cells[:, 1]] = wgt / wgt.sum()
        rows.append(r)
    w = cfg["window"]
    for i, j in cfg["windows"]:
        r = np.zeros((nrow, ncol))
        wgt = 10.0 ** x[i:i + w, j:j + w]
        r[i:i + w, j:j + w] = wgt / wgt.sum()
        rows.append(r.ravel())
    return np.array(rows)


def par_names(nrow, ncol):
    return ["k_{0:04d}_{1:04d}".format(i, j) for i in range(nrow) for j in range(ncol)]


def geostruct(range_cells=8.0, sill=0.25):
    """exponential variogram in cell units.  pyemu's a is the e-folding length,
    the practical range is nearer 3a.  sill 0.25 is a log10 std of 0.5"""
    return pyemu.geostats.GeoStruct(variograms=[pyemu.geostats.ExpVario(contribution=sill, a=range_cells)])


def draw_fields(nrow, ncol, num_reals, range_cells=8.0, sill=0.25, seed=12345):
    """log10 fields of shape (num_reals, nrow, ncol), mean zero, by spectral
    simulation so the size does not matter"""
    gs = geostruct(range_cells, sill)
    ss = pyemu.geostats.SpecSim2d(delx=np.ones(ncol), dely=np.ones(nrow), geostruct=gs)
    return ss.draw_arrays(num_reals=num_reals, mean_value=0.0, rng=np.random.RandomState(seed))


def model_precision(nrow, ncol, range_cells=8.0, sill=0.25, radius=3):
    """a sparse prior precision from the variogram itself, no draws involved.

    each cell is regressed on the cells within radius that come before it in
    raster order (the half plane above and to the left), with the regression
    coefficients taken from the exact exponential covariance.  that gives one
    row of a unit lower triangular factor L and a conditional variance d, and
    the precision is L^T diag(1/d) L: symmetric positive definite by
    construction, with a bandwidth set by radius.  returns the parameter names
    and the precision as a scipy csr matrix in that order.

    this is the input the enif graph path could not get from the ensemble at
    npar >> N: a precision that is right because the prior is known, not
    because it was estimated."""
    import scipy.sparse as sp
    names = par_names(nrow, ncol)
    offs = [(di, dj) for di in range(-radius, 1) for dj in range(-radius, radius + 1)
            if ((di < 0) or (di == 0 and dj < 0)) and np.hypot(di, dj) <= radius]
    offs = np.array(offs)

    def cov(d):
        return sill * np.exp(-d / range_cells)

    rows, cols, vals = [], [], []
    for i in range(nrow):
        for j in range(ncol):
            node = i * ncol + j
            pi, pj = i + offs[:, 0], j + offs[:, 1]
            ok = (pi >= 0) & (pi < nrow) & (pj >= 0) & (pj < ncol)
            if ok.sum() == 0:
                d = sill
                b = np.zeros(0)
                preds = np.zeros(0, dtype=int)
            else:
                pi, pj = pi[ok], pj[ok]
                preds = pi * ncol + pj
                dd = np.hypot(pi[:, None] - pi[None, :], pj[:, None] - pj[None, :])
                cnn = cov(dd)
                cni = cov(np.hypot(pi - i, pj - j))
                b = np.linalg.solve(cnn, cni)
                d = max(sill - cni @ b, 1.0e-8 * sill)
            sc = 1.0 / np.sqrt(d)
            rows.append(node)
            cols.append(node)
            vals.append(sc)
            for k, pr in enumerate(preds):
                rows.append(node)
                cols.append(pr)
                vals.append(-b[k] * sc)
    n = nrow * ncol
    L = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
    prec = (L.T @ L).tocsr()
    return names, prec


def distance_localizer(cfg, nrow, ncol, radius_cells):
    """a gaspari-cohn taper on the distance from each parameter cell to the
    nearest cell of each observation's footprint: the cell for a point, the
    line for a path, the window for a block.  returns the nobs x npar array
    in obs_names() by par_names() order"""
    w = cfg["window"]
    foot = []
    for i, j in cfg["points"]:
        foot.append(np.array([[i, j]]))
    for cells in cfg["paths"]:
        foot.append(np.asarray(cells))
    for i, j in cfg["windows"]:
        ii, jj = np.meshgrid(np.arange(i, i + w), np.arange(j, j + w), indexing="ij")
        foot.append(np.column_stack([ii.ravel(), jj.ravel()]))
    gi, gj = np.meshgrid(np.arange(nrow), np.arange(ncol), indexing="ij")
    gi, gj = gi.ravel(), gj.ravel()
    loc = np.zeros((len(foot), nrow * ncol))
    c = float(radius_cells)
    for k, f in enumerate(foot):
        d = np.full(nrow * ncol, np.inf)
        for (fi, fj) in f:
            d = np.minimum(d, np.hypot(gi - fi, gj - fj))
        z = d / c
        t = np.zeros_like(z)
        m = z <= 1.0
        t[m] = -0.25 * z[m] ** 5 + 0.5 * z[m] ** 4 + 0.625 * z[m] ** 3 - (5.0 / 3.0) * z[m] ** 2 + 1.0
        m = (z > 1.0) & (z <= 2.0)
        t[m] = (1.0 / 12.0) * z[m] ** 5 - 0.5 * z[m] ** 4 + 0.625 * z[m] ** 3 + (5.0 / 3.0) * z[m] ** 2 - 5.0 * z[m] + 4.0 - 2.0 / (3.0 * z[m])
        loc[k] = np.clip(t, 0.0, 1.0)
    return loc


def block_localizer(cfg, nrow, ncol, tile, radius_cells):
    """a 0/1 localizer on parameter tiles: parameters are grouped into tile x
    tile blocks and an observation may inform a block when any cell of its
    footprint is within radius_cells of the block.  pestpp then solves one
    local problem per block instead of one per parameter, which is what makes
    a localized run affordable at thousands of parameters.  returns the
    nobs x nblocks array, the block names in column order, and the block name
    of every parameter in par_names() order"""
    nbi = int(np.ceil(nrow / tile))
    nbj = int(np.ceil(ncol / tile))
    bnames = ["blk_{0:03d}_{1:03d}".format(bi, bj) for bi in range(nbi) for bj in range(nbj)]
    par_block = ["blk_{0:03d}_{1:03d}".format(i // tile, j // tile) for i in range(nrow) for j in range(ncol)]
    w = cfg["window"]
    foot = []
    for i, j in cfg["points"]:
        foot.append(np.array([[i, j]]))
    for cells in cfg["paths"]:
        foot.append(np.asarray(cells))
    for i, j in cfg["windows"]:
        ii, jj = np.meshgrid(np.arange(i, i + w), np.arange(j, j + w), indexing="ij")
        foot.append(np.column_stack([ii.ravel(), jj.ravel()]))
    loc = np.zeros((len(foot), len(bnames)))
    for b, name in enumerate(bnames):
        bi, bj = int(name.split("_")[1]), int(name.split("_")[2])
        i0, i1 = bi * tile, min((bi + 1) * tile, nrow) - 1
        j0, j1 = bj * tile, min((bj + 1) * tile, ncol) - 1
        for k, f in enumerate(foot):
            # distance from each footprint cell to the block, zero inside it
            di = np.maximum(0, np.maximum(i0 - f[:, 0], f[:, 0] - i1))
            dj = np.maximum(0, np.maximum(j0 - f[:, 1], f[:, 1] - j1))
            if np.hypot(di, dj).min() <= radius_cells:
                loc[k, b] = 1.0
    return loc, bnames, par_block


def run_forward(ws="."):
    """the model command.  reads the config and the k array the template
    wrote, writes obs.dat"""
    with open(os.path.join(ws, "anafield_config.json")) as f:
        cfg = json.load(f)
    k = np.loadtxt(os.path.join(ws, "k.dat")).reshape(cfg["nrow"], cfg["ncol"])
    y = forward(np.log10(k), cfg)
    with open(os.path.join(ws, "obs.dat"), "w") as f:
        for n, v in zip(obs_names(cfg), y):
            f.write("{0} {1:20.12E}\n".format(n, v))


def run_forward_rns(ws=".", case="anafield"):
    """the model command for pestpp's external run manager (the /e switch).
    pestpp writes every run of the batch to <case>.rns and calls the model
    once; this fills in the observations for all of them in one numpy pass
    and marks them complete.  parameters in the file are model values (k),
    not log values"""
    from pyemu.utils.helpers import RunStor
    with open(os.path.join(ws, "anafield_config.json")) as f:
        cfg = json.load(f)
    fname = os.path.join(ws, case + ".rns")
    _, pnames, onames = RunStor.file_info(fname)
    rs = RunStor(fname)
    df = rs.get_data()
    # the file keeps its own name order (hashed), so pick columns by name
    lower = {n.lower(): n for n in pnames}
    pcols = [lower[n] for n in par_names(cfg["nrow"], cfg["ncol"])]
    lower = {n.lower(): n for n in onames}
    ocols = [lower[n] for n in obs_names(cfg)]
    kvals = df.loc[:, pcols].values.astype(float)
    y = np.array([forward(np.log10(k).reshape(cfg["nrow"], cfg["ncol"]), cfg) for k in kvals])
    df.loc[:, ocols] = y
    df.loc[:, "run_status"] = 1
    rs.update(df)


def setup(new_d="anafield_template", nrow=20, ncol=20, n_points=20, n_paths=10,
          n_windows=5, window=4, a=2.0, range_cells=8.0, sill=0.25,
          noise_std=0.05, num_reals=100, seed=12345, cov_max_pars=5000,
          graph_kind="rook"):
    """write the whole problem.  one extra realization is drawn and becomes
    the truth, so the truth is in the prior but not in the ensemble.
    returns the pst.

    the prior covariance is dense, nrow*ncol squared, so it is only written
    when the parameter count is at or below cov_max_pars.  the prior ensemble
    and the graph are written at any size."""
    if os.path.exists(new_d):
        shutil.rmtree(new_d)
    os.makedirs(new_d)
    cfg = make_config(nrow, ncol, n_points, n_paths, n_windows, window, a, seed)
    with open(os.path.join(new_d, "anafield_config.json"), "w") as f:
        json.dump(cfg, f)

    names = par_names(nrow, ncol)
    with open(os.path.join(new_d, "k.dat.tpl"), "w") as f:
        f.write("ptf ~\n")
        for i in range(nrow):
            f.write(" ".join("~ {0} ~".format(names[i * ncol + j]) for j in range(ncol)) + "\n")
    onames = obs_names(cfg)
    with open(os.path.join(new_d, "obs.dat.ins"), "w") as f:
        f.write("pif ~\n")
        for n in onames:
            f.write("l1 w !{0}!\n".format(n))

    fields = draw_fields(nrow, ncol, num_reals + 1, range_cells, sill, seed)
    truth = fields[0]
    np.savetxt(os.path.join(new_d, "truth_logk.dat"), truth, fmt="%15.6E")
    np.savetxt(os.path.join(new_d, "k.dat"), 10.0 ** truth, fmt="%15.6E")

    # the forward script is this module, run as the model command.  full path
    # so it works when ci runs the tests from another folder
    shutil.copy2(os.path.abspath(__file__), os.path.join(new_d, "anafield_setup.py"))
    with open(os.path.join(new_d, "forward_run.py"), "w") as f:
        f.write("import anafield_setup\nif __name__ == '__main__':\n    anafield_setup.run_forward()\n")
    # the batch version for the external run manager: model command
    # "python forward_run_rns.py" and run pestpp with /e
    with open(os.path.join(new_d, "forward_run_rns.py"), "w") as f:
        f.write("import anafield_setup\nif __name__ == '__main__':\n    anafield_setup.run_forward_rns()\n")
    run_forward(new_d)

    pst = pyemu.Pst.from_io_files([os.path.join(new_d, "k.dat.tpl")], [os.path.join(new_d, "k.dat")],
                                  [os.path.join(new_d, "obs.dat.ins")], [os.path.join(new_d, "obs.dat")],
                                  pst_path=".")
    par = pst.parameter_data
    par.loc[:, "partrans"] = "log"
    par.loc[:, "parval1"] = 1.0
    par.loc[:, "parlbnd"] = 1.0e-3
    par.loc[:, "parubnd"] = 1.0e+3
    par.loc[:, "pargp"] = "k"
    par.loc[:, "i"] = [int(n.split("_")[1]) for n in par.parnme]
    par.loc[:, "j"] = [int(n.split("_")[2]) for n in par.parnme]

    rng = np.random.default_rng(seed + 1)
    obs = pst.observation_data
    y = forward(truth, cfg)
    obs.loc[onames, "obsval"] = y + noise_std * rng.standard_normal(len(y))
    obs.loc[onames, "weight"] = 1.0 / noise_std
    obs.loc[onames, "obgnme"] = [n.split("_")[0] for n in onames]
    pd.DataFrame({"obsnme": onames, "truth": y}).to_csv(os.path.join(new_d, "truth_obs.csv"), index=False)

    pst.model_command = "python forward_run.py"
    pst.control_data.noptmax = 0
    pst.pestpp_options["ies_par_en"] = "prior_pe.jcb"
    pst.pestpp_options["ies_num_reals"] = num_reals

    reals = fields[1:].reshape(num_reals, nrow * ncol)
    pe = pyemu.ParameterEnsemble(pst=pst, df=pd.DataFrame(10.0 ** reals, columns=names,
                                 index=["real_{0}".format(r) for r in range(num_reals)]))
    pe.to_binary(os.path.join(new_d, "prior_pe.jcb"))

    if len(names) <= cov_max_pars:
        pd_ = pst.parameter_data
        cov = geostruct(range_cells, sill).covariance_matrix(pd_["j"].values.astype(float),
                                                             pd_["i"].values.astype(float),
                                                             names=names)
        cov.to_coo(os.path.join(new_d, "prior_cov.jcb"))
        pst.pestpp_options["parcov_filename"] = "prior_cov.jcb"
    else:
        print("  {0} parameters > cov_max_pars {1}: prior covariance not written".format(
            len(names), cov_max_pars))

    import synth2d_graph
    graph = synth2d_graph.lattice_graph(pst, kind=graph_kind)
    graph.to_coo(os.path.join(new_d, "graph_{0}.jcb".format(graph_kind)))

    pst.write(os.path.join(new_d, "anafield.pst"), version=2)
    print("anafield: {0}x{1} grid, {2} parameters, {3} observations, a={4}, range {5} cells".format(
        nrow, ncol, pst.npar_adj, pst.nnz_obs, a, range_cells))
    return pst


if __name__ == "__main__":
    setup()
