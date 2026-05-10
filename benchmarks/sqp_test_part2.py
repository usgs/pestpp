import os
import sys
import shutil
import platform
import re
import numpy as np
import pandas as pd
import platform
import pyemu

bin_path = os.path.join("test_bin")
if "linux" in platform.platform().lower():
    bin_path = os.path.join(bin_path,"linux")
elif "darwin" in platform.platform().lower():
    bin_path = os.path.join(bin_path,"mac")
else:
    bin_path = os.path.join(bin_path,"win")

bin_path = os.path.abspath("test_bin")
os.environ["PATH"] += os.pathsep + bin_path

bin_path = os.path.join("..","..","..","bin", "win")
exe = ""
if "windows" in platform.platform().lower():
    exe = ".exe"
exe_path = os.path.join(bin_path, "pestpp-sqp" + exe)

noptmax = 3
num_reals = 15
port = 4025

def _read_rec(m_d, pst_name):
    base = os.path.splitext(pst_name)[0]
    rec_path = os.path.join(m_d, base + ".rec")
    assert os.path.exists(rec_path), \
        "rec file missing: {0}".format(rec_path)
    with open(rec_path, "r") as fp:
        return rec_path, fp.read()

def _setup_rosenbrock_workdir(case_name, initial_decvars=(0.1, -1.6),
                              constraint_rhs=-1.5, noptmax=3, num_reals=15, subset_size=6,
                              random_seed=1234):

    w_d = os.path.join("rosenbrock", "part2_" + case_name)
    if os.path.exists(w_d):
        shutil.rmtree(w_d)
    shutil.copytree(os.path.join("rosenbrock", "template"), w_d)

    in_file = "par.dat"
    tpl_file = in_file + ".tpl"
    out_files = ["obs.dat", "constraints.dat"]
    ins_files = [f + ".ins" for f in out_files]

    cwd = os.getcwd()
    os.chdir(w_d)
    try:
        pst = pyemu.helpers.pst_from_io_files(tpl_file, in_file,
                                              ins_files, out_files)

        par = pst.parameter_data
        par.loc[:, "partrans"] = "none"
        par.loc[par.parnme[0], "parval1"] = initial_decvars[0]
        par.loc[par.parnme[1], "parval1"] = initial_decvars[1]
        par.loc[:, "parlbnd"] = -2.2
        par.loc[:, "parubnd"] = 2.2
        par.loc[:, "parchglim"] = "relative"

        obs = pst.observation_data
        obs.loc["obs", "obgnme"] = "obj_fn"
        obs.loc["obs", "obsval"] = 0.0
        obs.loc["constraint", "obgnme"] = "l_constraint"
        obs.loc["constraint", "obsval"] = constraint_rhs
        obs.loc[:, "weight"] = 1.0

        pst.pestpp_options["opt_obj_func"] = "obs"
        pst.pestpp_options["sqp_num_reals"] = num_reals
        pst.pestpp_options["sqp_update_hessian"] = "true"
        pst.pestpp_options["sqp_subset_size"] = subset_size
        pst.pestpp_options["par_sigma_range"] = 10
        pst.pestpp_options["random_seed"] = random_seed

        pst.model_command = ["python " + "rosenbrock_2par_one_linear_constrained.py"]
        pst.control_data.noptmax = noptmax
        pst_name = "part2_" + case_name + ".pst"
        pst.write(pst_name)
    finally:
        os.chdir(cwd)
    return os.path.abspath(w_d), pst_name

def _reset_master_dir(model_d, name):
    m_d = os.path.join(model_d, name)
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    return m_d


def _run_sqp_parallel(t_d, model_d, pst_name, m_d, port, num_workers=8):
    pyemu.os_utils.start_workers(t_d, exe_path, pst_name,num_workers=num_workers, master_dir=m_d,
                                 worker_root=model_d, port=port)

def _parse_filter_rec(rec_path):
    header_re = re.compile(r"SQP filter members \((?P<count>\d+)\) for iteration (?P<iter>-?\d+):")
    pass_re = re.compile(r"number of candidate realizations passing filter:\s*(?P<count>\d+)")
    none_re = re.compile(r"no realizations passed the filter")

    members = {}
    passed_per_event = []
    none_events = 0
    with open(rec_path, "r") as fp:
        for line in fp:
            m = header_re.search(line)
            if m:
                it = int(m.group("iter"))
                members.setdefault(it, []).append(int(m.group("count")))
                continue
            m = pass_re.search(line)
            if m:
                passed_per_event.append(int(m.group("count")))
                continue
            if none_re.search(line):
                none_events += 1
    return members, passed_per_event, none_events


def _candidate_csv(m_d, base, iteration):
    return os.path.join(m_d, "{0}.{1}.dv_candidates.csv".format(base, iteration))


def _cma_cov_path(m_d, base, iteration):
    return os.path.join(m_d, "{0}.{1}.CMA.cov".format(base, iteration))


def sqp_filter_test():
    # check filter populates and admits at least one candidate per iteration
    w_d, pst_name = _setup_rosenbrock_workdir("filter", noptmax=3, num_reals=15)
    base = os.path.splitext(pst_name)[0]

    m_d = _reset_master_dir("rosenbrock", "master_filter")
    _run_sqp_parallel(w_d, "rosenbrock", pst_name, m_d, port + 1)

    rec_path, _ = _read_rec(m_d, pst_name)
    members, passed_events, none_events = _parse_filter_rec(rec_path)

    assert len(members) >= 1, "no SQP filter membership reported in rec file"
    for it, counts in members.items():
        assert all(c >= 1 for c in counts), "filter is empty at iteration {0}: {1}".format(it, counts)
    assert any(c >= 1 for c in passed_events), "no candidate ever passed the filter (events={0})".format(passed_events)
    assert none_events == 0 or any(c >= 1 for c in passed_events), "filter rejected every candidate at every iteration"

    #check iter-specific dv/oe candidate output files exist
    pst = pyemu.Pst(os.path.join(w_d, pst_name))
    for i in range(1, pst.control_data.noptmax + 1):
        for ext in (".dv_candidates.csv", ".oe_candidates.csv"):
            f = os.path.join(m_d, "{0}.{1}{2}".format(base, i, ext))
            assert os.path.exists(f), "missing line-search candidate output: {0}".format(f)


def sqp_filter_tol_test():
    pass_counts = {}
    for label, tol in [("strict", 1.0e-6), ("relaxed", 5.0e-2)]:
        w_d, pst_name = _setup_rosenbrock_workdir(
            "filter_tol_" + label, noptmax=3, num_reals=15)
        base = os.path.splitext(pst_name)[0]

        pst = pyemu.Pst(os.path.join(w_d, pst_name))
        pst.pestpp_options["sqp_filter_tol"] = tol
        pst.write(os.path.join(w_d, pst_name))

        m_d = _reset_master_dir("rosenbrock", "master_filter_tol_" + label)
        _run_sqp_parallel(w_d, "rosenbrock", pst_name, m_d, port + 2)

        rec_path, _ = _read_rec(m_d, pst_name)
        _, passed_events, _ = _parse_filter_rec(rec_path)
        pass_counts[label] = sum(passed_events)

    assert pass_counts["strict"] >= pass_counts["relaxed"], "strict filter tol passed fewer candidates than relaxed: {0}".format(pass_counts)

def sqp_line_search_alpha_test():
    # the dv_candidates row count equals (subset_size + 1) * len(alpha_mults)
    alpha_mults = [0.05, 0.25, 0.5, 1.0]
    subset_size = 6
    w_d, pst_name = _setup_rosenbrock_workdir("ls_alpha", noptmax=2, num_reals=15)
    base = os.path.splitext(pst_name)[0]

    pst = pyemu.Pst(os.path.join(w_d, pst_name))
    pst.pestpp_options["sqp_alpha_mults"] = ",".join("{0}".format(a) for a in alpha_mults)
    pst.pestpp_options["sqp_subset_size"] = subset_size
    pst.pestpp_options["sqp_num_refined_search_pts"] = 0
    pst.write(os.path.join(w_d, pst_name))

    m_d = _reset_master_dir("rosenbrock", "master_ls_alpha")
    _run_sqp_parallel(w_d, "rosenbrock", pst_name, m_d, port + 4)

    cand = _candidate_csv(m_d, base, 1)
    assert os.path.exists(cand), "iter-1 dv_candidates.csv missing: {0}".format(cand)

    df = pd.read_csv(cand, index_col=0)
    expected = (subset_size + 1) * len(alpha_mults)
    assert df.shape[0] == expected, "iter-1 dv_candidates row count: expected {0}, got {1} (shape={2})".format(expected, df.shape[0], df.shape)

    # check companion line-search outputs are also present
    for ext in (".oe_candidates.csv", ".pcs.csv"):
        f = os.path.join(m_d, "{0}.1{1}".format(base, ext))
        assert os.path.exists(f), "missing line-search output: {0}".format(f)


def sqp_line_search_subset_test():
    # check increasing sqp_subset_size increases the dv_candidates row count
    alpha_mults = "0.1, 0.5, 1.0"
    n_alpha = len(alpha_mults.split(","))
    rows = {}
    for subset_size in (3, 8):
        w_d, pst_name = _setup_rosenbrock_workdir(
            "ls_subset_{0}".format(subset_size), noptmax=1, num_reals=15)
        base = os.path.splitext(pst_name)[0]

        pst = pyemu.Pst(os.path.join(w_d, pst_name))
        pst.pestpp_options["sqp_alpha_mults"] = alpha_mults
        pst.pestpp_options["sqp_subset_size"] = subset_size
        pst.pestpp_options["sqp_num_refined_search_pts"] = 0
        pst.write(os.path.join(w_d, pst_name))

        m_d = _reset_master_dir("rosenbrock", "master_ls_subset_{0}".format(subset_size))
        _run_sqp_parallel(w_d, "rosenbrock", pst_name, m_d, port + 5)

        df = pd.read_csv(_candidate_csv(m_d, base, 1), index_col=0)
        rows[subset_size] = df.shape[0]

        # exact count check too: (subset_size + 1) * n_alpha
        expected = (subset_size + 1) * n_alpha
        assert df.shape[0] == expected, "subset={0}: expected {1} candidates, got {2}".format(subset_size, expected, df.shape[0])

    assert rows[3] < rows[8], "increasing subset size did not increase candidate rows: {0}".format(rows)


# def sqp_line_search_refined_test():
#     # sqp_num_refined_search_pts > 0 should add candidate rows
#     counts = {}
#     for refined in (0, 2):
#         w_d, pst_name = _setup_rosenbrock_workdir(
#             "ls_refine_{0}".format(refined), noptmax=1, num_reals=15)
#         base = os.path.splitext(pst_name)[0]

#         pst = pyemu.Pst(os.path.join(w_d, pst_name))
#         pst.pestpp_options["sqp_alpha_mults"] = "0.1, 0.5, 1.0"
#         pst.pestpp_options["sqp_subset_size"] = 6
#         pst.pestpp_options["sqp_num_refined_search_pts"] = refined
#         pst.write(os.path.join(w_d, pst_name))

#         m_d = _reset_master_dir("rosenbrock", "master_ls_refine_{0}".format(refined))
#         _run_sqp_parallel(w_d, "rosenbrock", pst_name, m_d, port + 6)

#         df = pd.read_csv(_candidate_csv(m_d, base, 1), index_col=0)
#         counts[refined] = df.shape[0]

#     assert counts[2] > counts[0], \
#         "refined search did not produce extra candidates: {0}".format(counts)


def sqp_cma_save_cov_test():
    case_noptmax = 3
    w_d, pst_name = _setup_rosenbrock_workdir("cma_save_cov", noptmax=case_noptmax, num_reals=12)
    base = os.path.splitext(pst_name)[0]

    pst = pyemu.Pst(os.path.join(w_d, pst_name))
    pst.pestpp_options["sqp_save_cov_every"] = 1
    pst.write(os.path.join(w_d, pst_name))

    m_d = _reset_master_dir("rosenbrock", "master_cma_save_cov")
    _run_sqp_parallel(w_d, "rosenbrock", pst_name, m_d, port + 7)

    expected = [_cma_cov_path(m_d, base, i) for i in range(1, case_noptmax + 1)]
    missing = [p for p in expected if not os.path.exists(p)]
    assert not missing, "CMA cov files missing: {0}".format(missing)

    cov = pyemu.Cov.from_ascii(expected[0])
    n_dv = pst.npar_adj
    assert cov.shape == (n_dv, n_dv), "CMA cov shape mismatch: expected ({0}, {0}), got {1}".format(n_dv, cov.shape)


def sqp_cma_evolves_test():
    #check that the CMA covariance changes between iter 1 and iter noptmax
    #and remains symmetric / positive semi-definite throughout
    case_noptmax = 3
    w_d, pst_name = _setup_rosenbrock_workdir("cma_evolves", noptmax=case_noptmax, num_reals=14)
    base = os.path.splitext(pst_name)[0]

    pst = pyemu.Pst(os.path.join(w_d, pst_name))
    pst.pestpp_options["sqp_save_cov_every"] = 1
    pst.write(os.path.join(w_d, pst_name))

    m_d = _reset_master_dir("rosenbrock", "master_cma_evolves")
    _run_sqp_parallel(w_d, "rosenbrock", pst_name, m_d, port + 9)

    cov_1 = pyemu.Cov.from_ascii(_cma_cov_path(m_d, base, 1)).x
    cov_n = pyemu.Cov.from_ascii(_cma_cov_path(m_d, base, case_noptmax)).x
    assert not np.allclose(cov_1, cov_n), "CMA covariance is unchanged from iter 1 to iter {0}".format(case_noptmax)

    # check covariance is symmetric and PSD at each saved iter
    for i in (1, case_noptmax):
        cov = pyemu.Cov.from_ascii(_cma_cov_path(m_d, base, i)).x
        assert np.allclose(cov, cov.T, atol=1e-8), "CMA covariance at iter {0} is not symmetric".format(i)
        eigvals = np.linalg.eigvalsh((cov + cov.T) * 0.5)
        assert (eigvals > -1e-8).all(), "CMA covariance at iter {0} has negative eigenvalues: {1}".format(i, eigvals)

def sqp_cma_custom_rates_test():
    # user-supplied CMA learning rates should not crash the optimisation
    case_noptmax = 3
    w_d, pst_name = _setup_rosenbrock_workdir("cma_custom_rates", noptmax=case_noptmax, num_reals=14)
    base = os.path.splitext(pst_name)[0]

    pst = pyemu.Pst(os.path.join(w_d, pst_name))
    pst.pestpp_options["sqp_cma_c1"] = 0.05
    pst.pestpp_options["sqp_cma_cmu"] = 0.10
    pst.pestpp_options["sqp_cma_cc"] = 0.20
    pst.pestpp_options["sqp_cma_parent_num"] = 4
    pst.pestpp_options["sqp_save_cov_every"] = 1
    pst.write(os.path.join(w_d, pst_name))

    m_d = _reset_master_dir("rosenbrock", "master_cma_custom_rates")
    _run_sqp_parallel(w_d, "rosenbrock", pst_name, m_d, port + 10)

    for i in range(0, case_noptmax + 1):
        for ext in (".par.csv", ".obs.csv"):
            f = os.path.join(m_d, "{0}.{1}{2}".format(base, i, ext))
            assert os.path.exists(f), "missing per-iter output: {0}".format(f)
    assert os.path.exists(_cma_cov_path(m_d, base, case_noptmax)), "missing final CMA cov file"

def sqp_cma_reinflation_test():
    max_cond_num = 100.0 #use tight reinflation threhold to easily trigger reinflation and test its effect on CMA covariance condition number
    case_noptmax = 3
    w_d, pst_name = _setup_rosenbrock_workdir("cma_reinflate", noptmax=case_noptmax, num_reals=14)
    base = os.path.splitext(pst_name)[0]

    pst = pyemu.Pst(os.path.join(w_d, pst_name))
    pst.pestpp_options["sqp_cma_reinflation_factor"] = 1.5
    pst.pestpp_options["sqp_max_reinflation_cond_num"] = max_cond_num
    pst.pestpp_options["sqp_save_cov_every"] = 1
    pst.pestpp_options["sqp_debug_cma"] = "True"
    pst.write(os.path.join(w_d, pst_name))

    m_d = _reset_master_dir("rosenbrock", "master_cma_reinflate")
    _run_sqp_parallel(w_d, "rosenbrock", pst_name, m_d, port + 11)

    # check iter-specific output files exist
    cov_paths = []
    for i in range(1, case_noptmax + 1):
        cov_path = _cma_cov_path(m_d, base, i)
        assert os.path.exists(cov_path), "missing CMA cov for iter {0}".format(i)
        cov_paths.append(cov_path)
        for ext in (".par.csv", ".obs.csv"):
            f = os.path.join(m_d, "{0}.{1}{2}".format(base, i, ext))
            assert os.path.exists(f), "missing per-iter output: {0}".format(f)

    # parse CMA condition-number metrics written by report_cma_metrics:
    _, rec_text = _read_rec(m_d, pst_name)
    cond_re = re.compile(r"Condition No\.\s+([\d.eE+\-]+)\s+([\d.eE+\-]+)")
    prior_conds, post_conds = [], []
    for line in rec_text.splitlines():
        m = cond_re.search(line)
        if m:
            prior_conds.append(float(m.group(1)))
            post_conds.append(float(m.group(2)))

    assert len(prior_conds) > 0, \
        "no CMA condition-number metrics found in rec file (sqp_debug_cma must be true)"

    # at least one pre-update condition number must have exceeded the threshold,
    # confirming the covariance was actually ill-conditioned enough to need reinflation
    assert any(c > max_cond_num for c in prior_conds), (
        "no pre-update condition number exceeded threshold {0}; "
        "reinflation was never needed: prior_conds={1}".format(max_cond_num, prior_conds))

    # every post-update condition number must be within the threshold,
    # confirming reinflation successfully capped the condition number
    assert all(c <= max_cond_num * 1.05 for c in post_conds), ("post-update condition numbers exceed threshold {0}: {1}".format(max_cond_num, post_conds))

    for i, cov_path in enumerate(cov_paths, start=1):
        mat = pyemu.Cov.from_ascii(cov_path).x
        sym = (mat + mat.T) * 0.5
        eigvals = np.linalg.eigvalsh(sym)
        cond_num = eigvals.max() / (eigvals.min() + 1e-10)
        assert cond_num <= max_cond_num * 1.05, (
            "condition number at iter {0} ({1:.2f}) exceeds threshold ({2}); "
            "reinflation did not cap it".format(i, cond_num, max_cond_num))

if __name__ == "__main__":
    # sqp_filter_test()
    # sqp_filter_tol_test()
    # sqp_line_search_alpha_test()
    # sqp_line_search_subset_test()
    # sqp_line_search_refined_test()
    # sqp_cma_save_cov_test()
    # sqp_cma_evolves_test()
    # sqp_cma_custom_rates_test()
    sqp_cma_reinflation_test()

