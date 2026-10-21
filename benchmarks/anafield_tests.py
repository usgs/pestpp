"""tests on the anafield benchmark - the analytical, nonlinear, scalable
problem in anafield_setup.py.  no model binaries needed, the forward run is
numpy, so these run anywhere pestpp-ies does."""
import os
import json
import shutil
import platform
import numpy as np
import pandas as pd
import pyemu

import anafield_setup

bin_path = os.path.join("test_bin")
if "linux" in platform.platform().lower():
    bin_path = os.path.join(bin_path, "linux")
elif "darwin" in platform.platform().lower() or "macos" in platform.platform().lower():
    bin_path = os.path.join(bin_path, "mac")
else:
    bin_path = os.path.join(bin_path, "win")
bin_path = os.path.abspath("test_bin")
os.environ["PATH"] += os.pathsep + bin_path

# same convention as the other benchmark scripts: the exe path is relative to a
# run dir two levels below this one, and the repo dir is called pestpp
if os.path.exists(os.path.join("pestpp", "bin")):
    bin_path = os.path.join("..", "..", "pestpp", "bin")
else:
    bin_path = os.path.join("..", "..", "..", "..", "pestpp", "bin")
if "windows" in platform.platform().lower():
    exe_path = os.path.join(bin_path, "win", "pestpp-ies.exe")
elif "darwin" in platform.platform().lower() or "macos" in platform.platform().lower():
    exe_path = os.path.join(bin_path, "mac", "pestpp-ies")
else:
    exe_path = os.path.join(bin_path, "linux", "pestpp-ies")


def anafield_forward_test():
    """the closed form jacobian agrees with finite differences for every
    observation type, the a=0 point obs is linear, and the path and block
    obs are bracketed by the min and max of their cells"""
    cfg = anafield_setup.make_config(nrow=12, ncol=15, n_points=6, n_paths=4, n_windows=3,
                                     window=3, a=2.0, seed=7)
    x = anafield_setup.draw_fields(12, 15, 1, range_cells=4.0, seed=7)[0]
    y = anafield_setup.forward(x, cfg)
    jac = anafield_setup.jacobian(x, cfg)
    assert jac.shape == (len(y), 12 * 15)
    fd = np.zeros_like(jac)
    h = 1.0e-6
    for c in range(12 * 15):
        xp = x.ravel().copy()
        xp[c] += h
        fd[:, c] = (anafield_setup.forward(xp.reshape(12, 15), cfg) - y) / h
    err = np.abs(fd - jac).max()
    assert err < 1.0e-4, "jacobian vs finite difference: {0}".format(err)

    # path rows sum to one (softmin weights), block rows too
    npt = len(cfg["points"])
    assert np.allclose(jac[npt:].sum(axis=1), 1.0)

    for k, cells in enumerate(cfg["paths"]):
        cells = np.asarray(cells)
        xv = x[cells[:, 0], cells[:, 1]]
        yv = y[npt + k]
        assert xv.min() - 1.0e-9 <= yv <= xv.mean() + 1.0e-9
    w = cfg["window"]
    for k, (i, j) in enumerate(cfg["windows"]):
        xv = x[i:i + w, j:j + w]
        yv = y[npt + len(cfg["paths"]) + k]
        assert xv.mean() - 1.0e-9 <= yv <= xv.max() + 1.0e-9

    cfg0 = dict(cfg)
    cfg0["a"] = 0.0
    y0 = anafield_setup.forward(x, cfg0)
    pts = np.asarray(cfg["points"])
    assert np.allclose(y0[:npt], x[pts[:, 0], pts[:, 1]])
    # sinh grows faster than linear away from zero
    assert np.all(np.abs(y[:npt]) >= np.abs(y0[:npt]) - 1.0e-12)


def anafield_ies_test(nrow=20, ncol=20, num_reals=50, noptmax=3):
    """build a small instance, check the written files reproduce the truth,
    run pestpp-ies in serial and check phi drops well below the prior"""
    model_d = "anafield"
    template_d = os.path.join(model_d, "template")
    os.makedirs(model_d, exist_ok=True)
    anafield_setup.setup(new_d=template_d, nrow=nrow, ncol=ncol, n_points=30, n_paths=15,
                         n_windows=8, window=4, a=2.0, range_cells=6.0, num_reals=num_reals)
    pst = pyemu.Pst(os.path.join(template_d, "anafield.pst"))
    assert pst.npar_adj == nrow * ncol
    assert pst.nnz_obs == 30 + 15 + 8
    for f in ("prior_pe.jcb", "prior_cov.jcb", "graph_rook.jcb", "truth_logk.dat"):
        assert os.path.exists(os.path.join(template_d, f)), f

    # the forward script run on the truth gives the noise free truth obs
    pst.control_data.noptmax = 0
    truth = np.loadtxt(os.path.join(template_d, "truth_logk.dat"))
    pst.parameter_data.loc[:, "parval1"] = 10.0 ** truth.ravel()
    m_d = os.path.join(model_d, "master_truth")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    shutil.copytree(template_d, m_d)
    pst.write(os.path.join(m_d, "anafield.pst"), version=2)
    pyemu.os_utils.run("{0} anafield.pst".format(exe_path), cwd=m_d)
    res = pyemu.Pst(os.path.join(m_d, "anafield.pst")).res
    tru = pd.read_csv(os.path.join(template_d, "truth_obs.csv"), index_col=0)
    diff = np.abs(res.loc[tru.index, "modelled"].values - tru["truth"].values).max()
    assert diff < 1.0e-4, "truth run vs truth obs: {0}".format(diff)

    # the prior ensemble has the variance the variogram says it should
    pe = pyemu.ParameterEnsemble.from_binary(pst, os.path.join(template_d, "prior_pe.jcb"))
    sd = np.log10(pe.values).std()
    assert 0.35 < sd < 0.65, "prior log10 std {0}".format(sd)

    # the ies run goes through the external run manager: pestpp writes each
    # batch to anafield.rns, forward_run_rns.py fills it in one numpy pass
    pst = pyemu.Pst(os.path.join(template_d, "anafield.pst"))
    pst.control_data.noptmax = noptmax
    pst.pestpp_options["ies_num_reals"] = num_reals
    pst.pestpp_options["ies_include_base"] = "false"
    pst.pestpp_options["random_seed"] = 112358
    pst.model_command = "python forward_run_rns.py"
    m_d = os.path.join(model_d, "master_ies")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    shutil.copytree(template_d, m_d)
    pst.write(os.path.join(m_d, "anafield.pst"), version=2)
    pyemu.os_utils.run("{0} anafield.pst /e".format(exe_path), cwd=m_d)

    # the batch path gives the same numbers as calling forward() on the prior
    pe0 = pyemu.ParameterEnsemble.from_csv(pst, os.path.join(m_d, "anafield.0.par.csv"))
    oe0 = pyemu.ObservationEnsemble.from_csv(pst, os.path.join(m_d, "anafield.0.obs.csv"))
    names = anafield_setup.par_names(nrow, ncol)
    with open(os.path.join(template_d, "anafield_config.json")) as f:
        cfg = json.load(f)
    for real in oe0.index[:5]:
        x = np.log10(pe0.loc[real, names].values.astype(float)).reshape(nrow, ncol)
        y = anafield_setup.forward(x, cfg)
        d = np.abs(oe0.loc[real, anafield_setup.obs_names(cfg)].values - y).max()
        assert d < 1.0e-6, "external run manager vs forward() on {0}: {1}".format(real, d)

    df = pd.read_csv(os.path.join(m_d, "anafield.phi.actual.csv"))
    phi = df["mean"].values
    print("ies mean phi by iteration:", np.round(phi, 1))
    assert len(phi) == noptmax + 1
    assert phi[-1] < 0.25 * phi[0], "phi did not drop enough: {0:.1f} -> {1:.1f}".format(phi[0], phi[-1])

    # the posterior mean field is closer to the truth than the prior mean was
    pe1 = pyemu.ParameterEnsemble.from_csv(pst, os.path.join(m_d, "anafield.{0}.par.csv".format(noptmax)))
    rmse0 = np.sqrt(np.mean((np.log10(pe0.loc[:, names].values).mean(axis=0) - truth.ravel()) ** 2))
    rmse1 = np.sqrt(np.mean((np.log10(pe1.loc[:, names].values).mean(axis=0) - truth.ravel()) ** 2))
    print("rmse of the mean log10 field vs truth: prior {0:.3f}, posterior {1:.3f}".format(rmse0, rmse1))
    assert rmse1 < rmse0
    return m_d


def _run_measured(cmd, cwd, log="run_stdout.txt"):
    """run a command, return wall seconds and peak rss in gb.  the command
    runs under a one-line python wrapper so the child rusage is this run's
    alone; stdout goes to log in cwd"""
    import subprocess
    import sys
    import time
    wrapper = ("import resource,subprocess,sys;"
               "r=subprocess.run(sys.argv[1],shell=True);"
               "print('MAXRSS',resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss);"
               "sys.exit(r.returncode)")
    t0 = time.time()
    with open(os.path.join(cwd, log), "w") as f:
        proc = subprocess.run([sys.executable, "-c", wrapper, cmd], cwd=cwd, stdout=f,
                              stderr=subprocess.STDOUT)
    wall = time.time() - t0
    if proc.returncode != 0:
        raise Exception("'{0}' failed with code {1} in {2}, see {3}".format(cmd, proc.returncode, cwd, log))
    peak = 0.0
    for line in open(os.path.join(cwd, log)):
        if line.startswith("MAXRSS"):
            peak = float(line.split()[1])
            # linux reports kb, mac bytes
            peak = peak / 1.0e6 if platform.system() == "Linux" else peak / 1.0e9
    return wall, peak


def anafield_scaling_compare(nrow=224, ncol=224, obs_sizes=(2000, 5000, 10000, 50000),
                             enif_max_obs=10000, num_reals=50, noptmax=3, range_cells=8.0, exe=None,
                             methods=("ies", "enif"), lasso_frac=0.01, num_threads=8):
    """ies against the enif woodbury form as the observation count climbs, at
    about 50k parameters.  not a ci test (no _test suffix): the full size
    takes a while and the point is the numbers, not a pass or fail.

    the parameter field is the same at every size (same seed) and the
    observations are split 30/50/20 between points, paths and blocks.  no
    prior covariance file - it is dense and 50k squared is 20 gb - so enif
    runs with ies_enif_current_prec, the hessian from the current ensemble,
    which is the variant that matched ies on the earlier problems.

    the woodbury path forms three dense matrices sized by the observation
    count: H^T and C H^T (parameters x observations) and the innovation
    covariance G (observations x observations) which it then factors.  at
    50k observations that is 20 gb each and an o(n^3) factor, where ies stays
    ensemble-sized throughout.  enif_max_obs stops the enif run above that
    count; raise it to try anyway.

    methods picks what runs: "ies", "enif" (woodbury, implicit H) and
    "enif_lasso" (woodbury with the explicit lasso H, ies_enif_h_lasso =
    lasso_frac, the lasso threaded over num_threads).  rows for methods not
    run this time are kept from an earlier csv, so the variants can be run
    one at a time and still land in one table and one pdf.

    writes anafield_scaling.csv and anafield_scaling.pdf in the work dir and
    returns the summary dataframe."""
    import time
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    if exe is None:
        exe = exe_path
    model_d = "anafield_scaling"
    os.makedirs(model_d, exist_ok=True)
    all_methods = ("ies", "enif", "enif_lasso")
    csv_name = os.path.join(model_d, "anafield_scaling.csv")
    rows = []
    if os.path.exists(csv_name):
        old = pd.read_csv(csv_name)
        rows = old.loc[~old.method.isin(methods)].to_dict("records")
    phis = {}
    for nobs in obs_sizes:
        template_d = os.path.join(model_d, "template_{0}".format(nobs))
        n_points = int(0.3 * nobs)
        n_paths = int(0.5 * nobs)
        n_windows = nobs - n_points - n_paths
        if not os.path.exists(os.path.join(template_d, "anafield.pst")):
            t0 = time.time()
            anafield_setup.setup(new_d=template_d, nrow=nrow, ncol=ncol, n_points=n_points, n_paths=n_paths,
                                 n_windows=n_windows, window=4, a=2.0, range_cells=range_cells,
                                 num_reals=num_reals, cov_max_pars=0)
            print("setup {0} obs took {1:.0f} s".format(nobs, time.time() - t0))
        for method in all_methods:
            if method not in methods:
                continue
            if method.startswith("enif") and nobs > enif_max_obs:
                print("skipping {0} at {1} obs (enif_max_obs {2})".format(method, nobs, enif_max_obs))
                continue
            pst = pyemu.Pst(os.path.join(template_d, "anafield.pst"))
            pst.control_data.noptmax = noptmax
            pst.pestpp_options["ies_num_reals"] = num_reals
            pst.pestpp_options["ies_include_base"] = "false"
            pst.pestpp_options["random_seed"] = 112358
            pst.pestpp_options["save_binary"] = "true"
            pst.model_command = "python forward_run_rns.py"
            if method.startswith("enif"):
                pst.pestpp_options["ies_use_enif"] = "true"
                pst.pestpp_options["ies_enif_current_prec"] = "true"
            if method == "enif_lasso":
                pst.pestpp_options["ies_enif_h_lasso"] = lasso_frac
                pst.pestpp_options["ies_num_threads"] = num_threads
            m_d = os.path.join(model_d, "master_{0}_{1}".format(method, nobs))
            if os.path.exists(m_d):
                shutil.rmtree(m_d)
            shutil.copytree(template_d, m_d)
            pst.write(os.path.join(m_d, "anafield.pst"), version=2)
            wall, peak = _run_measured("{0} anafield.pst /e".format(exe), m_d)
            df = pd.read_csv(os.path.join(m_d, "anafield.phi.actual.csv"))
            phi = df["mean"].values
            phis[(method, nobs)] = pd.read_csv(os.path.join(m_d, "anafield.phi.actual.csv"), index_col=0)
            # how much of the wall was the forward runs: the batch calls are echoed to stdout
            fwd = 0.0
            for line in open(os.path.join(m_d, "run_stdout.txt")):
                if "forward run command(s) finished, took" in line:
                    fwd += float(line.split("took")[1].split()[0])
            rows.append({"method": method, "nobs": nobs, "npar": nrow * ncol, "wall_s": wall,
                         "forward_s": fwd, "solve_s": wall - fwd, "peak_rss_gb": peak,
                         "phi0": phi[0], "phi_final": phi[-1], "phi_final_per_obs": phi[-1] / nobs})
            print("{0} {1} obs: phi {2:.0f} -> {3:.0f}, wall {4:.0f} s (forward {5:.0f} s), peak rss {6:.2f} gb".format(
                method, nobs, phi[0], phi[-1], wall, fwd, peak))
            pd.DataFrame(rows).to_csv(csv_name, index=False)

    summ = pd.DataFrame(rows).sort_values(["method", "nobs"])
    # phi traces of methods run earlier, from their master dirs
    for method in all_methods:
        for nobs in obs_sizes:
            f = os.path.join(model_d, "master_{0}_{1}".format(method, nobs), "anafield.phi.actual.csv")
            if (method, nobs) not in phis and os.path.exists(f):
                phis[(method, nobs)] = pd.read_csv(f, index_col=0)
    colors = {"ies": "tab:blue", "enif": "tab:red", "enif_lasso": "tab:green"}
    with PdfPages(os.path.join(model_d, "anafield_scaling.pdf")) as pdf:
        labels = {"ies": "ies", "enif": "enif woodbury, implicit H", "enif_lasso": "enif woodbury, lasso H"}
        fig, axes = plt.subplots(1, 4, figsize=(18, 4.5))
        for method in all_methods:
            d = summ.loc[summ.method == method]
            if d.shape[0] == 0:
                continue
            axes[0].plot(d.nobs, d.phi_final_per_obs, "o-", color=colors[method], label=labels[method])
            axes[1].plot(d.nobs, d.solve_s, "o-", color=colors[method], label=labels[method])
            axes[2].plot(d.nobs, d.wall_s, "o-", color=colors[method], label=labels[method])
            axes[3].plot(d.nobs, d.peak_rss_gb, "o-", color=colors[method], label=labels[method])
        for ax, lab in zip(axes, ("final mean phi per observation", "solve time, forward runs excluded (s)",
                                  "total wall time (s)", "peak memory (gb)")):
            ax.set_xscale("log")
            ax.set_xlabel("observations")
            ax.set_title(lab, fontsize=10)
            ax.grid(alpha=0.3)
            ax.legend(fontsize=8)
        axes[1].set_yscale("log")
        axes[2].set_yscale("log")
        fig.suptitle("anafield {0}x{1} = {2} parameters, {3} reals, {4} iterations".format(
            nrow, ncol, nrow * ncol, num_reals, noptmax))
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
        # phi per realization by iteration, one page per size
        for nobs in obs_sizes:
            fig, ax = plt.subplots(1, 1, figsize=(8, 5))
            for method in all_methods:
                if (method, nobs) not in phis:
                    continue
                df = phis[(method, nobs)]
                reals = df.iloc[:, 5:]
                for c in reals.columns:
                    ax.plot(df.index, reals[c].values / nobs, color=colors[method], lw=0.4, alpha=0.4)
                ax.plot(df.index, df["mean"].values / nobs, color=colors[method], lw=2.5, label=labels[method] + " mean")
            ax.set_yscale("log")
            ax.set_xlabel("iteration")
            ax.set_ylabel("phi per observation")
            ax.set_title("{0} observations".format(nobs))
            ax.grid(alpha=0.3)
            ax.legend()
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
    return summ


def anafield_model_precision_test():
    """the variogram-based precision is spd, its inverse is close to the exact
    covariance, and the localizer is one on the footprint and zero far away"""
    nrow, ncol = 14, 14
    names, Q = anafield_setup.model_precision(nrow, ncol, range_cells=6.0, sill=0.25, radius=4)
    Qd = Q.toarray()
    assert np.allclose(Qd, Qd.T)
    assert np.linalg.eigvalsh(Qd).min() > 0.0
    C = np.linalg.inv(Qd)
    ii, jj = np.meshgrid(np.arange(nrow), np.arange(ncol), indexing="ij")
    ii, jj = ii.ravel(), jj.ravel()
    Ct = 0.25 * np.exp(-np.hypot(ii[:, None] - ii[None, :], jj[:, None] - jj[None, :]) / 6.0)
    rms = np.sqrt(((C - Ct) ** 2).mean()) / 0.25
    assert rms < 0.08, "model precision inverse vs exact covariance, rms {0}".format(rms)
    assert np.abs(np.diag(C) - 0.25).max() / 0.25 < 0.2

    cfg = anafield_setup.make_config(nrow, ncol, 3, 2, 2, window=3, a=1.0, seed=3)
    loc = anafield_setup.distance_localizer(cfg, nrow, ncol, radius_cells=3.0)
    assert loc.shape == (7, nrow * ncol)
    i, j = cfg["points"][0]
    assert loc[0, i * ncol + j] == 1.0
    assert (loc >= 0.0).all() and (loc <= 1.0).all()
    far = np.hypot(ii - i, jj - j) > 6.0
    assert np.all(loc[0, far] == 0.0)


def anafield_reach_compare(nrow=64, ncol=64, nobs=1000, num_reals=500, noptmax=8, range_cells=8.0,
                           prec_radius=4, loc_radius=None, loc_tile=8, lasso_frac=0.01, num_threads=8, exe=None,
                           variants=None, model_d="anafield_reach"):
    """the experiment the enif graph path was built for: can a prior precision
    that is NOT estimated from the ensemble move the solution outside the
    ensemble span, and does that beat ies with localization?  not a ci test.

    the variants, all with the same prior ensemble, noise and seed:
      ies                 plain
      ies_loc             block localizer: loc_tile x loc_tile parameter blocks, an
                          observation informs a block when its footprint is within
                          loc_radius cells.  pestpp solves per block (64 of them at
                          64x64, tile 8); a per-parameter taper would be 4096 sequential
                          solves per upgrade, which took 50 minutes each here
      ies_aal             automatic adaptive localization - same per-parameter cost, so
                          not in the default list
      enif_wood_cur       woodbury, implicit H, current covariance (the control = ies)
      enif_graph_est      rook graph, lasso H, prior precision ESTIMATED from the prior ensemble
      enif_graph_model    lasso H, prior precision from the variogram (ies_enif_prec_file)
      enif_graph_model_hybrid   as above with lambda on the current-ensemble precision
      enif_graph_est_hybrid     estimated prior precision with the hybrid damping

    scoring: phi, rmse of the mean log field against the truth, spread, the
    fraction of the mean update that lies outside the span of the prior
    anomalies, and wall time.  writes anafield_reach.csv and anafield_reach.pdf
    in the work dir."""
    import time
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    if exe is None:
        exe = exe_path
    if loc_radius is None:
        loc_radius = 2.0 * range_cells
    template_d = os.path.join(model_d, "template")
    os.makedirs(model_d, exist_ok=True)
    n_points = int(0.3 * nobs)
    n_paths = int(0.5 * nobs)
    n_windows = nobs - n_points - n_paths
    if not os.path.exists(os.path.join(template_d, "anafield.pst")):
        anafield_setup.setup(new_d=template_d, nrow=nrow, ncol=ncol, n_points=n_points, n_paths=n_paths,
                             n_windows=n_windows, window=4, a=2.0, range_cells=range_cells,
                             num_reals=num_reals, cov_max_pars=0)
        pst = pyemu.Pst(os.path.join(template_d, "anafield.pst"))
        names, Q = anafield_setup.model_precision(nrow, ncol, range_cells, 0.25, radius=prec_radius)
        pyemu.Matrix(x=Q.toarray(), row_names=names, col_names=names).to_coo(
            os.path.join(template_d, "prec_model.jcb"))
        with open(os.path.join(template_d, "anafield_config.json")) as f:
            cfg = json.load(f)
        loc, bnames, par_block = anafield_setup.block_localizer(cfg, nrow, ncol, loc_tile, loc_radius)
        pyemu.Matrix(x=loc, row_names=anafield_setup.obs_names(cfg), col_names=bnames).to_coo(
            os.path.join(template_d, "loc.jcb"))
        pd.Series(par_block, index=names).to_csv(os.path.join(template_d, "par_blocks.csv"), header=False)
        print("model precision {0} non-zeros, block localizer {1} blocks, mean {2:.3f}".format(
            Q.nnz, len(bnames), loc.mean()))

    enif = {"ies_use_enif": "true"}
    graph = dict(enif, ies_enif_h_lasso=lasso_frac, ies_num_threads=num_threads)
    all_variants = {
        "ies": {},
        "ies_loc": {"ies_localizer": "loc.jcb", "ies_num_threads": num_threads},
        "ies_aal": {"ies_autoadaloc": "true", "ies_num_threads": num_threads},
        "enif_wood_cur": dict(enif, ies_enif_current_prec="true"),
        "enif_graph_est": dict(graph, ies_enif_graph="graph_rook.jcb"),
        "enif_graph_model": dict(graph, ies_enif_prec_file="prec_model.jcb"),
        "enif_graph_model_hybrid": dict(graph, ies_enif_prec_file="prec_model.jcb", ies_enif_hybrid_damp="true"),
        "enif_graph_est_hybrid": dict(graph, ies_enif_graph="graph_rook.jcb", ies_enif_hybrid_damp="true"),
    }
    if variants is None:
        variants = [v for v in all_variants if v != "ies_aal"]

    csv_name = os.path.join(model_d, "anafield_reach.csv")
    rows = []
    if os.path.exists(csv_name):
        old = pd.read_csv(csv_name)
        rows = old.loc[~old.variant.isin(variants)].to_dict("records")

    truth = np.loadtxt(os.path.join(template_d, "truth_logk.dat")).ravel()
    pst = pyemu.Pst(os.path.join(template_d, "anafield.pst"))
    names = anafield_setup.par_names(nrow, ncol)
    pe0 = pyemu.ParameterEnsemble.from_binary(pst, os.path.join(template_d, "prior_pe.jcb"))
    X0 = np.log10(pe0.loc[:, names].values.astype(float))
    m0 = X0.mean(axis=0)
    # orthonormal basis of the prior anomaly span
    Ub, _, _ = np.linalg.svd((X0 - m0).T, full_matrices=False)

    for tag in variants:
        opts = all_variants[tag]
        pst = pyemu.Pst(os.path.join(template_d, "anafield.pst"))
        pst.control_data.noptmax = noptmax
        pst.control_data.phiredstp = 0.0
        pst.control_data.nphistp = 100
        pst.control_data.nphinored = 100
        pst.pestpp_options["ies_num_reals"] = num_reals
        pst.pestpp_options["ies_include_base"] = "false"
        pst.pestpp_options["random_seed"] = 112358
        pst.pestpp_options["save_binary"] = "true"
        pst.pestpp_options.update(opts)
        if "ies_localizer" in opts:
            # the localizer columns are parameter blocks, so the blocks are the groups
            blocks = pd.read_csv(os.path.join(template_d, "par_blocks.csv"), header=None, index_col=0)[1]
            pst.parameter_data.loc[blocks.index, "pargp"] = blocks.values
            pst.rectify_pgroups()
        pst.model_command = "python forward_run_rns.py"
        m_d = os.path.join(model_d, "master_" + tag)
        if os.path.exists(m_d):
            shutil.rmtree(m_d)
        shutil.copytree(template_d, m_d)
        pst.write(os.path.join(m_d, "anafield.pst"), version=2)
        try:
            wall, peak = _run_measured("{0} anafield.pst /e".format(exe), m_d)
        except Exception as ex:
            print("{0} FAILED: {1}".format(tag, ex))
            rows.append({"variant": tag, "failed": True})
            pd.DataFrame(rows).to_csv(csv_name, index=False)
            continue
        df = pd.read_csv(os.path.join(m_d, "anafield.phi.actual.csv"))
        last = int(df["iteration"].values[-1])
        pe = pyemu.ParameterEnsemble.from_binary(pst, os.path.join(m_d, "anafield.{0}.par.jcb".format(last)))
        X = np.log10(pe.loc[:, names].values.astype(float))
        m = X.mean(axis=0)
        upd = m - m0
        inside = Ub @ (Ub.T @ upd)
        frac_out = np.linalg.norm(upd - inside) / max(np.linalg.norm(upd), 1.0e-30)
        rmse_out = np.sqrt(np.mean((m0 + inside - truth) ** 2))
        lam = np.nan
        for line in open(os.path.join(m_d, "anafield.rec")):
            if "current lambda" in line:
                try:
                    lam = float(line.split(":")[-1])
                except ValueError:
                    pass
        rows.append({"variant": tag, "failed": False, "iterations": last, "wall_s": wall, "peak_rss_gb": peak,
                     "phi0": df["mean"].values[0], "phi_final": df["mean"].values[-1],
                     "phi_final_per_obs": df["mean"].values[-1] / nobs,
                     "rmse_prior_mean": np.sqrt(np.mean((m0 - truth) ** 2)),
                     "rmse_post_mean": np.sqrt(np.mean((m - truth) ** 2)),
                     "rmse_post_mean_inspan_only": rmse_out,
                     "spread_prior": X0.std(axis=0).mean(), "spread_post": X.std(axis=0).mean(),
                     "frac_update_outside_span": frac_out, "last_lambda": lam})
        print("{0}: phi {1:.0f} -> {2:.0f} in {3} iters, rmse {4:.3f} -> {5:.3f}, spread {6:.3f} -> {7:.3f}, "
              "outside span {8:.3f}, wall {9:.0f} s".format(
                  tag, df["mean"].values[0], df["mean"].values[-1], last,
                  rows[-1]["rmse_prior_mean"], rows[-1]["rmse_post_mean"],
                  rows[-1]["spread_prior"], rows[-1]["spread_post"], frac_out, wall))
        pd.DataFrame(rows).to_csv(csv_name, index=False)

    summ = pd.DataFrame(rows)
    order = [v for v in all_variants if v in set(summ.variant)]
    cmap = plt.get_cmap("tab10")
    colors = {v: cmap(k) for k, v in enumerate(all_variants)}
    with PdfPages(os.path.join(model_d, "anafield_reach.pdf")) as pdf:
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        for tag in order:
            f = os.path.join(model_d, "master_" + tag, "anafield.phi.actual.csv")
            if not os.path.exists(f):
                continue
            df = pd.read_csv(f, index_col=0)
            reals = df.iloc[:, 5:]
            for c in reals.columns[::10]:
                ax.plot(df.index, reals[c].values / nobs, color=colors[tag], lw=0.3, alpha=0.3)
            ax.plot(df.index, df["mean"].values / nobs, color=colors[tag], lw=2.5, label=tag)
        ax.set_yscale("log")
        ax.set_xlabel("iteration")
        ax.set_ylabel("phi per observation (every 10th realization thin, mean thick)")
        ax.set_title("anafield {0}x{1}, {2} obs, {3} reals".format(nrow, ncol, nobs, num_reals))
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        ok = summ.loc[~summ.failed.astype(bool)].set_index("variant").loc[[v for v in order if v in summ.loc[~summ.failed.astype(bool), "variant"].values]]
        fig, axes = plt.subplots(2, 3, figsize=(16, 8))
        panels = [("phi_final_per_obs", "final phi per obs"), ("rmse_post_mean", "rmse of mean log field vs truth"),
                  ("spread_post", "mean posterior spread (log10)"), ("frac_update_outside_span", "fraction of mean update outside prior span"),
                  ("wall_s", "wall time (s)"), ("last_lambda", "last lambda")]
        for ax, (col, lab) in zip(axes.ravel(), panels):
            ax.bar(range(len(ok)), ok[col].values, color=[colors[v] for v in ok.index])
            ax.set_xticks(range(len(ok)))
            ax.set_xticklabels(ok.index, rotation=60, ha="right", fontsize=7)
            ax.set_title(lab, fontsize=10)
            ax.grid(alpha=0.3, axis="y")
            if col in ("wall_s", "last_lambda", "phi_final_per_obs"):
                ax.set_yscale("log")
        if "rmse_prior_mean" in ok.columns:
            axes[0, 1].axhline(ok["rmse_prior_mean"].values[0], color="k", ls="--", lw=1, label="prior mean")
            axes[0, 1].legend(fontsize=8)
            axes[0, 2].axhline(ok["spread_prior"].values[0], color="k", ls="--", lw=1)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # fields: truth, prior mean, posterior mean per variant
        n = len(order) + 2
        ncolf = 4
        nrowf = int(np.ceil(n / ncolf))
        fig, axes = plt.subplots(nrowf, ncolf, figsize=(4 * ncolf, 3.6 * nrowf))
        axes = np.atleast_1d(axes).ravel()
        lim = np.abs(truth).max()
        axes[0].imshow(truth.reshape(nrow, ncol), vmin=-lim, vmax=lim, cmap="viridis")
        axes[0].set_title("truth", fontsize=9)
        axes[1].imshow(m0.reshape(nrow, ncol), vmin=-lim, vmax=lim, cmap="viridis")
        axes[1].set_title("prior mean", fontsize=9)
        for k, tag in enumerate(order):
            f = os.path.join(model_d, "master_" + tag)
            files = sorted([x for x in os.listdir(f) if x.endswith(".par.jcb") and x.startswith("anafield.")
                            and x.split(".")[1].isdigit()], key=lambda x: int(x.split(".")[1]))
            if not files:
                continue
            pe = pyemu.ParameterEnsemble.from_binary(pst, os.path.join(f, files[-1]))
            m = np.log10(pe.loc[:, names].values.astype(float)).mean(axis=0)
            axes[k + 2].imshow(m.reshape(nrow, ncol), vmin=-lim, vmax=lim, cmap="viridis")
            axes[k + 2].set_title("{0} mean, rmse {1:.3f}".format(tag, np.sqrt(np.mean((m - truth) ** 2))), fontsize=9)
        for ax in axes:
            ax.set_xticks([])
            ax.set_yticks([])
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
    return summ


if __name__ == "__main__":
    anafield_forward_test()
    anafield_model_precision_test()
    anafield_ies_test()
