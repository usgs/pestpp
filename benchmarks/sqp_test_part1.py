import os
import sys
import shutil
import platform
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

bin_path = os.path.join("..","..","..","bin")
exe = ""
if "windows" in platform.platform().lower():
    exe = ".exe"
exe_path = os.path.join(bin_path, "pestpp-sqp" + exe)


noptmax = 4
num_reals = 20
port = 4021


def _read_rec(m_d, pst_name):
    base = os.path.splitext(pst_name)[0]
    rec_path = os.path.join(m_d, base + ".rec")
    assert os.path.exists(rec_path), \
        "rec file missing: {0}".format(rec_path)
    with open(rec_path, "r") as fp:
        return rec_path, fp.read()

def _setup_rosenbrock_workdir(case_name, initial_decvars=(0.1, -1.6),
                              constraint_rhs=-1.5, noptmax=3, num_reals=15, subset_size= 10,
                              random_seed=8):

    w_d = os.path.join("rosenbrock", "part1_" + case_name)
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
        pst.pestpp_options["sqp_subset_size"] = 10
        pst.pestpp_options["par_sigma_range"] = 10
        pst.pestpp_options["random_seed"] = random_seed

        pst.model_command = ["python " + "rosenbrock_2par_one_linear_constrained.py"]
        pst.control_data.noptmax = noptmax
        pst_name = "part1_" + case_name + ".pst"
        pst.write(pst_name)
    finally:
        os.chdir(cwd)
    return os.path.abspath(w_d), pst_name


def _run_sqp_serial(w_d, pst_name):
    pyemu.os_utils.run("{0} {1}".format(exe_path, pst_name), cwd=w_d)

def _prep_freyberg_sqp_pst():
    model_d = "freyberg_sqp"
    t_d = os.path.join(model_d, "template")
    pst = pyemu.Pst(os.path.join(t_d, "freyberg_run_mou.pst"))

    for key in ["mou_dv_population_file", "mou_objectives", "mou_population_size",
                "mou_save_population_every", "forecasts", "opt_objective_function"]:
        pst.pestpp_options.pop(key, None)

    obj_obs = "oname:cum_otype:lst_usecol:wel_totim:4383.5"
    pst.observation_data.loc[obj_obs, "obgnme"] = "obj_fn"
    pst.pestpp_options["opt_obj_func"] = obj_obs

    return model_d, t_d, pst


def _reset_master_dir(model_d, name):
    m_d = os.path.join(model_d, name)
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    return m_d


def _run_sqp_parallel(t_d, model_d, pst_name, m_d, port, num_workers=8):
    pyemu.os_utils.start_workers(t_d, exe_path, pst_name, num_workers=num_workers, master_dir=m_d,
                                 worker_root=model_d, port=port)

def basic_sqp_init_test():
    # verify the initial BASE par/rei
    w_d, pst_name = _setup_rosenbrock_workdir("init_n0", noptmax=0, num_reals=10)
    base = os.path.splitext(pst_name)[0]
    _run_sqp_serial(w_d, pst_name)

    assert os.path.exists(os.path.join(w_d, base + ".base.par")), "noptmax=0 run did not produce {0}.base.par".format(base)
    assert os.path.exists(os.path.join(w_d, base + ".base.rei")), "noptmax=0 run did not produce {0}.base.rei".format(base)

    w_d, pst_name = _setup_rosenbrock_workdir("init_ensemble", noptmax=-1, num_reals=10)
    base = os.path.splitext(pst_name)[0]
    _run_sqp_serial(w_d, pst_name)

    # check rec and log are produced
    for ext in (".rec", ".log"):
        assert os.path.exists(os.path.join(w_d, base + ext)), "expected {0}{1} after noptmax=-1 run".format(base, ext)

    # check initial ensemble par/obs 
    par_csv = os.path.join(w_d, base + ".0.par.csv")
    obs_csv = os.path.join(w_d, base + ".0.obs.csv")
    assert os.path.exists(par_csv), "noptmax=-1 run did not produce initial parameter ensemble {0}.0.par.csv".format(base)
    assert os.path.exists(obs_csv), "noptmax=-1 run did not produce initial observation ensemble {0}.0.obs.csv".format(base)

    pst = pyemu.Pst(os.path.join(w_d, pst_name))
    n_reals = int(pst.pestpp_options["sqp_num_reals"])

    par_df = pd.read_csv(par_csv, index_col=0)
    assert par_df.shape == (n_reals + 1, pst.npar), "unexpected {0}.0.par.csv shape: expected ({1}, {2}), got {3}".format(base, n_reals + 1, pst.npar, par_df.shape)

    #check if BASE exists in initial pars
    base_real_names = [r for r in par_df.index.astype(str) if r.upper() == "BASE"]
    assert base_real_names, "BASE realization missing from initial par.csv"
    
    #check adjustable pars
    base_row = par_df.loc[base_real_names[0]]
    adj = pst.parameter_data[pst.parameter_data.partrans.isin(["none", "log"])]
    common = [p for p in adj.parnme if p in base_row.index]
    assert len(common) > 0, "no adjustable parameters found in par.csv columns"
    
    expected = adj.set_index("parnme").loc[common, "parval1"].astype(float)
    actual = base_row[common].astype(float)
    assert np.allclose(expected.values, actual.values, rtol=1e-6, atol=1e-8), "BASE realization parameters do not match parval1"

    # check obs shape
    obs_df = pd.read_csv(obs_csv, index_col=0)
    assert obs_df.shape == (n_reals + 1, pst.nobs + pst.nprior), "unexpected {0}.0.obs.csv shape: expected ({1}, {2}), got {3}".format(
            base, n_reals + 1, pst.nobs + pst.nprior, obs_df.shape)


def basic_sqp_iter_test():
    w_d, pst_name = _setup_rosenbrock_workdir("iter", noptmax=2, num_reals=10)
    base = os.path.splitext(pst_name)[0]

    alpha_mults = [0.1, 0.5, 1.0]
    pst = pyemu.Pst(os.path.join(w_d, pst_name))
    pst.pestpp_options["sqp_subset_size"] = 5
    pst.pestpp_options["sqp_alpha_mults"] = ", ".join(map(str, alpha_mults))
    pst.pestpp_options["sqp_num_refined_search_pts"] = 0
    pst.write(os.path.join(w_d, pst_name))

    m_d = _reset_master_dir("rosenbrock", "master_iter")
    _run_sqp_parallel(w_d, "rosenbrock", pst_name, m_d, port + 1)

    #check iter-specific output files
    for i in range(0, pst.control_data.noptmax + 1):
        for ext in (".par.csv", ".obs.csv", ".base.par", ".base.rei"):
            f = os.path.join(m_d, "{0}.{1}{2}".format(base, i, ext))
            assert os.path.exists(f), "missing per-iter output files: {0}".format(f)

    #check line search outputs
    for i in range(1, pst.control_data.noptmax + 1):
        for ext in (".dv_candidates.csv", ".oe_candidates.csv", ".pcs.csv"):
            f = os.path.join(m_d, "{0}.{1}{2}".format(base, i, ext))
            assert os.path.exists(f), "missing line-search candidate output files: {0}".format(f)

    #check candidate counts
    cand = pd.read_csv(os.path.join(m_d, "{0}.1.dv_candidates.csv".format(base)), index_col=0)
    expected = (pst.pestpp_options["sqp_subset_size"] + 1) * len(alpha_mults)
    assert cand.shape[0] == expected, "iter-1 dv_candidates row count: expected {0}, got {1} (shape={2})".format(expected, cand.shape[0], cand.shape)


def _read_base_obj(m_d, base, iter_i, obj_obs_name="obs"):
    rei_path = os.path.join(m_d, "{0}.{1}.base.rei".format(base, iter_i))
    rei = pyemu.pst_utils.read_resfile(rei_path)
    return float(rei.loc[obj_obs_name, "modelled"])


def basic_sqp_direction_test():
    base_pars = {}
    for label, direction in (("min", "min"), ("max", "max")):
        w_d, pst_name = _setup_rosenbrock_workdir("dir_" + label, noptmax=2, num_reals=10)
        base = os.path.splitext(pst_name)[0]

        pst = pyemu.Pst(os.path.join(w_d, pst_name))
        pst.pestpp_options["opt_direction"] = direction
        pst.pestpp_options["sqp_subset_size"] = 4
        pst.pestpp_options["sqp_alpha_mults"] = "0.1, 0.5, 1.0"
        pst.write(os.path.join(w_d, pst_name))

        m_d = _reset_master_dir("rosenbrock", "master_dir_" + label)
        _run_sqp_parallel(w_d, "rosenbrock", pst_name, m_d, port + 2)

        _, rec_text = _read_rec(m_d, pst_name)
        expected_token = "opt_direction: {0}".format(1 if direction == "min" else -1)
        assert expected_token in rec_text, "rec did not record '{0}' for direction='{1}'".format(expected_token, direction)

        noptmax_val = pst.control_data.noptmax
        par_path = os.path.join(m_d, base + ".{0}.base.par".format(noptmax_val))
        assert os.path.exists(par_path), "missing iter-{0} base.par for direction='{1}': {2}".format(noptmax_val, direction, par_path)
        base_pars[label] = pyemu.pst_utils.read_parfile(par_path)

        # check if obj progresses in the right direction
        obj0 = _read_base_obj(m_d, base, 0)
        obj_final = _read_base_obj(m_d, base, noptmax_val)
        if direction == "min":
            assert obj_final <= obj0 + 1e-6, "min run did not decrease obj: iter0={0:.6g}, iter{1}={2:.6g}".format(obj0, noptmax_val, obj_final)
        else:
            assert obj_final >= obj0 - 1e-6, "max run did not increase obj: iter0={0:.6g}, iter{1}={2:.6g}".format(obj0, noptmax_val, obj_final)

    #check if pars changed
    common = [p for p in base_pars["min"].index if p in base_pars["max"].index]
    assert len(common) > 0, "no common parameters between min and max base.par files"
    diffs = (base_pars["min"].loc[common, "parval1"] - base_pars["max"].loc[common, "parval1"]).abs()
    assert diffs.max() > 0.0, "min and max runs produced identical BASE parameters: max diff=0"

def basic_sqp_rosenbrock_chance_test():
    case_name = "rosenbrock_chance"
    w_d = os.path.join("rosenbrock", "part1_" + case_name)
    if os.path.exists(w_d):
        shutil.rmtree(w_d)
    shutil.copytree(os.path.join("rosenbrock", "template"), w_d)

    tpl_file = "par_chance.dat.tpl"
    in_file = "par.dat"
    out_files = ["obs.dat", "constraints.dat"]
    ins_files = [f + ".ins" for f in out_files]

    n_reals = 8
    stack_size = 10

    cwd = os.getcwd()
    os.chdir(w_d)
    try:
        pst = pyemu.helpers.pst_from_io_files(tpl_file, in_file, ins_files, out_files)

        par = pst.parameter_data
        for pname, val in [("par1", 0.1), ("par2", -1.6)]:
            par.loc[pname, "partrans"] = "none"
            par.loc[pname, "parval1"] = val
            par.loc[pname, "parlbnd"] = -2.2
            par.loc[pname, "parubnd"] = 2.2
            par.loc[pname, "parchglim"] = "relative"
            par.loc[pname, "pargp"] = "decvar"

        par.loc["theta", "partrans"] = "none"
        par.loc["theta", "parval1"] = 0.25
        par.loc["theta", "parlbnd"] = -1.0
        par.loc["theta", "parubnd"] = 1.5
        par.loc["theta", "parchglim"] = "relative"
        par.loc["theta", "pargp"] = "uncertain"

        obs = pst.observation_data
        obs.loc["obs", "obgnme"] = "obj_fn"
        obs.loc["obs", "obsval"] = 0.0
        obs.loc["constraint", "obgnme"] = "l_constraint"
        obs.loc["constraint", "obsval"] = -1.5
        obs.loc[:, "weight"] = 1.0

        pst.pestpp_options["opt_obj_func"] = "obs"
        pst.pestpp_options["opt_dec_var_groups"] = "decvar"
        pst.pestpp_options["sqp_risk"] = 0.95
        pst.pestpp_options["sqp_num_reals"] = n_reals
        pst.pestpp_options["sqp_subset_size"] = 5
        pst.pestpp_options["opt_stack_size"] = stack_size
        pst.pestpp_options["opt_chance_points"] = "SINGLE"
        pst.pestpp_options["sqp_alpha_mults"] = "0.1, 0.5, 1.0"
        pst.pestpp_options["par_sigma_range"] = 10
        pst.pestpp_options["sqp_update_hessian"] = "true"
        pst.model_command = ["python rosenbrock_2par_one_nonlinear_chance_constrained.py"]
        pst.control_data.noptmax = 1
        pst_name = "part1_" + case_name + ".pst"
        pst.write(pst_name)
    finally:
        os.chdir(cwd)

    m_d = _reset_master_dir("rosenbrock", "master_" + case_name)
    _run_sqp_parallel(w_d, "rosenbrock", pst_name, m_d, port + 6)

    base = os.path.splitext(pst_name)[0]
    assert os.path.exists(os.path.join(m_d, base + ".1.base.par")), "iter-1 base.par missing after rosenbrock chance run"

    # stack size checks
    par_stack_csv = os.path.join(m_d, base + ".1.par_stack.csv")
    obs_stack_csv = os.path.join(m_d, base + ".1.obs_stack.csv")
    assert os.path.exists(par_stack_csv), "iter-1 par_stack.csv missing"
    assert os.path.exists(obs_stack_csv), "iter-1 obs_stack.csv missing"

    par_stack_df = pd.read_csv(par_stack_csv, index_col=0)
    obs_stack_df = pd.read_csv(obs_stack_csv, index_col=0)

    assert par_stack_df.shape[0] == stack_size, "par_stack row count mismatch: expected {0} (opt_stack_size), got {1}".format(stack_size, par_stack_df.shape[0])
    assert obs_stack_df.shape[0] == stack_size, "obs_stack row count mismatch: expected {0} (opt_stack_size), got {1}".format(stack_size, obs_stack_df.shape[0])

    adj_pars = pst.parameter_data.loc[pst.parameter_data.partrans.isin(["none", "log"]), "parnme"].tolist()
    par_stack_cols = [c for c in adj_pars if c in par_stack_df.columns]
    assert len(par_stack_cols) == len(adj_pars), "par_stack column count mismatch: expected {0} adjustable pars, found {1} in csv".format(
            len(adj_pars), len(par_stack_cols))

    expected_obs_cols = pst.nobs + pst.nprior
    assert obs_stack_df.shape[1] == expected_obs_cols, "obs_stack column count mismatch: expected {0} (nobs+nprior), got {1}".format(
            expected_obs_cols, obs_stack_df.shape[1])

    # confirm non-decvar uncertainty was sampled
    assert "theta" in par_stack_df.columns, "non-DV parameter was not included in stack"
    assert par_stack_df["theta"].std() > 0, "non-DV uncertainty was not sampled"

    # checks for ALL chance points
    pst.pestpp_options["opt_chance_points"] = "ALL"
    pst_name_all = "part1_" + case_name + "_all.pst"
    cwd = os.getcwd()
    os.chdir(w_d)
    try:
        pst.write(pst_name_all)
    finally:
        os.chdir(cwd)

    m_d_all = _reset_master_dir("rosenbrock", "master_" + case_name + "_all")
    _run_sqp_parallel(w_d, "rosenbrock", pst_name_all, m_d_all, port + 7)

    base_all = os.path.splitext(pst_name_all)[0]
    assert os.path.exists(os.path.join(m_d_all, base_all + ".1.base.par")), "iter-1 base.par missing after ALL chance-points run"

    # check for stack size: total rows = n_reals * stack_size
    par_stack_all_csv = os.path.join(m_d_all, base_all + ".1.nested.par_stack.csv")
    obs_stack_all_csv = os.path.join(m_d_all, base_all + ".1.nested.obs_stack.csv")
    assert os.path.exists(par_stack_all_csv), "iter-1 par_stack.csv missing for ALL chance-points run"
    assert os.path.exists(obs_stack_all_csv), "iter-1 obs_stack.csv missing for ALL chance-points run"

    par_stack_all_df = pd.read_csv(par_stack_all_csv, index_col=0)
    obs_stack_all_df = pd.read_csv(obs_stack_all_csv, index_col=0)

    expected_all_rows = (n_reals + 1) * stack_size
    assert par_stack_all_df.shape[0] == expected_all_rows, "ALL par_stack row count mismatch: expected {0} (n_reals={1} * stack_size={2}), got {3}".format(
            expected_all_rows, n_reals, stack_size, par_stack_all_df.shape[0])
    assert obs_stack_all_df.shape[0] == expected_all_rows, "ALL obs_stack row count mismatch: expected {0} (n_reals={1} * stack_size={2}), got {3}".format(
            expected_all_rows, n_reals, stack_size, obs_stack_all_df.shape[0])

    # ALL chance points should produce more stack rows than SINGLE
    assert par_stack_all_df.shape[0] > par_stack_df.shape[0], "ALL par_stack ({0} rows) should be larger than SINGLE par_stack ({1} rows)".format(par_stack_all_df.shape[0], par_stack_df.shape[0])

    # stack dim checks
    par_stack_all_cols = [c for c in adj_pars if c in par_stack_all_df.columns]
    assert len(par_stack_all_cols) == len(adj_pars), "ALL par_stack column count mismatch: expected {0} adjustable pars, found {1}".format(
            len(adj_pars), len(par_stack_all_cols))
    assert obs_stack_all_df.shape[1] == expected_obs_cols, "ALL obs_stack column count mismatch: expected {0} (nobs+nprior), got {1}".format(
            expected_obs_cols, obs_stack_all_df.shape[1])

    # theta should still vary across all stack members
    assert "theta" in par_stack_all_df.columns, "theta not found in ALL par_stack columns"
    assert par_stack_all_df["theta"].std() > 0, "theta has zero spread in ALL par_stack"

# def basic_sqp_chance_test_highdim():
#     model_d, t_d, pst = _prep_freyberg_sqp_pst()

#     n_reals = 10
#     stack_size = 15
#     pst.pestpp_options["sqp_risk"] = 0.95
#     pst.pestpp_options["sqp_num_reals"] = n_reals
#     pst.pestpp_options["sqp_subset_size"] = 5
#     pst.pestpp_options["opt_stack_size"] = stack_size
#     pst.pestpp_options["opt_chance_points"] = "SINGLE"
#     pst.pestpp_options["sqp_alpha_mults"] = "0.1, 0.5, 1.0"
#     pst.pestpp_options["sqp_enforce_bounds"] = "true"   
#     pst.pestpp_options["random_seed"] = 8

#     pst.control_data.noptmax = 1
#     pst_name = "freyberg_run_sqp_chance.pst"
#     pst.write(os.path.join(t_d, pst_name), version=2)

#     m_d = _reset_master_dir(model_d, "master_sqp_chance")
#     _run_sqp_parallel(t_d, model_d, pst_name, m_d, port + 5)

#     _, rec_text = _read_rec(m_d, pst_name)
#     base = os.path.splitext(pst_name)[0]

#     assert os.path.exists(os.path.join(m_d, base + ".1.base.par")), "iter-1 base.par missing after chance run"
#     assert os.path.exists(os.path.join(m_d, base + ".1.par.csv")), "iter-1 par.csv missing after chance run"

#     # check stack generation
#     par_stack_csv = os.path.join(m_d, base + ".1.par_stack.csv")
#     obs_stack_csv = os.path.join(m_d, base + ".1.obs_stack.csv")
#     assert os.path.exists(par_stack_csv), "iter-1 par_stack.csv missing"
#     assert os.path.exists(obs_stack_csv), "iter-1 obs_stack.csv missing"

#     par_stack_df = pd.read_csv(par_stack_csv, index_col=0)
#     obs_stack_df = pd.read_csv(obs_stack_csv, index_col=0)

#     # check stack size
#     assert par_stack_df.shape[0] == stack_size, "par_stack row count mismatch: expected {0} (opt_stack_size), got {1}".format(stack_size, par_stack_df.shape[0])
#     assert obs_stack_df.shape[0] == stack_size, "obs_stack row count mismatch: expected {0} (opt_stack_size), got {1}".format( stack_size, obs_stack_df.shape[0])

#     adj_pars = pst.parameter_data.loc[pst.parameter_data.partrans.isin(["none", "log"]), "parnme"].tolist()
#     par_stack_cols = [c for c in adj_pars if c in par_stack_df.columns]
#     assert len(par_stack_cols) == len(adj_pars), "par_stack column count mismatch: expected {0} adjustable pars, found {1} in csv".format(len(adj_pars), len(par_stack_cols))

#     expected_obs_cols = pst.nobs + pst.nprior
#     assert obs_stack_df.shape[1] == expected_obs_cols, "obs_stack column count mismatch: expected {0} (nobs+nprior), got {1}".format(expected_obs_cols, obs_stack_df.shape[1])

#     constraint_groups = [g for g in pst.observation_data.obgnme.unique()
#                          if g.startswith(("l_", "g_", "e_","less_than", "greater_than", "equal_to"))]
#     assert len(constraint_groups) > 0, "no constraint observation groups found in pst"
#     constraint_obs = pst.observation_data.loc[pst.observation_data.obgnme.isin(constraint_groups), "obsnme"].tolist()
#     constraint_cols = [c for c in constraint_obs if c in obs_stack_df.columns]
#     assert len(constraint_cols) > 0, "no constraint obs found in obs_stack.csv columns"
#     constraint_stds = obs_stack_df[constraint_cols].std()
#     assert (constraint_stds > 0).any(), "all constraint obs have zero spread across stack members"

def basic_sqp_bounds_test():
    w_d, pst_name = _setup_rosenbrock_workdir("bounds", noptmax=1, num_reals=12, initial_decvars=(0.1, 0.1))
    base = os.path.splitext(pst_name)[0]

    pst = pyemu.Pst(os.path.join(w_d, pst_name))
    pst.parameter_data.loc[:, "parlbnd"] = 0.0
    pst.parameter_data.loc[:, "parubnd"] = 0.15
    pst.pestpp_options["sqp_enforce_bounds"] = "true"
    pst.write(os.path.join(w_d, pst_name))

    m_d = _reset_master_dir("rosenbrock", "master_bounds")
    _run_sqp_parallel(w_d, "rosenbrock", pst_name, m_d, port + 8)

    cand_csv = os.path.join(m_d, base + ".1.dv_candidates.csv")
    assert os.path.exists(cand_csv), "dv_candidates.csv missing after iter 1: {0}".format(cand_csv)
    cand_df = pd.read_csv(cand_csv, index_col=0)

    par = pst.parameter_data
    dv_pars = par.loc[par.partrans != "fixed", "parnme"].tolist()
    dv_pars = [p for p in dv_pars if p in cand_df.columns]
    assert len(dv_pars) > 0, "no DV columns found in dv_candidates.csv"

    lb = par.set_index("parnme").loc[dv_pars, "parlbnd"].astype(float)
    ub = par.set_index("parnme").loc[dv_pars, "parubnd"].astype(float)
    sub = cand_df[dv_pars].astype(float)

    tol = 1.0E-6
    too_low = (sub.subtract(lb, axis=1) < -tol).values
    too_high = (sub.subtract(ub, axis=1) > tol).values
    assert not too_low.any(), "dv_candidates below parlbnd with sqp_enforce_bounds=true"
    assert not too_high.any(), "dv_candidates above parubnd with sqp_enforce_bounds=true"

def _parse_hessians_from_rec(rec_path):
    import re
    hessians = []
    with open(rec_path) as f:
        lines = f.readlines()
    i = 0
    float_re = re.compile(r'^[\s\d.eE+\-]+$')
    while i < len(lines):
        if re.match(r'\s*hessian \(iter \d+\):', lines[i]):
            i += 1
            rows = []
            while i < len(lines) and float_re.match(lines[i]) and lines[i].strip():
                rows.append([float(v) for v in lines[i].split()])
                i += 1
            if rows:
                hessians.append(np.array(rows))
        else:
            i += 1
    return hessians


def basic_sqp_hessian_test():
    for label, update in (("bfgs", "true"), ("identity", "false")):
        w_d, pst_name = _setup_rosenbrock_workdir("hess_" + label, noptmax=2, num_reals=10)
        base = os.path.splitext(pst_name)[0]

        pst = pyemu.Pst(os.path.join(w_d, pst_name))
        pst.pestpp_options["sqp_update_hessian"] = update
        pst.pestpp_options["sqp_subset_size"] = 4
        pst.pestpp_options["sqp_alpha_mults"] = "0.1, 0.5, 1.0"
        if update == "true":
            pst.pestpp_options["sqp_debug_hessian"] = "true"
        pst.write(os.path.join(w_d, pst_name))

        m_d = _reset_master_dir("rosenbrock", "master_hess_" + label)
        _run_sqp_parallel(w_d, "rosenbrock", pst_name, m_d, port + 9)

        for i in range(0, pst.control_data.noptmax + 1):
            f = os.path.join(m_d, "{0}.{1}.par.csv".format(base, i))
            assert os.path.exists(f), "missing iter-specific par.csv for hessian='{0}': {1}".format(update, f)

        last_par = os.path.join(m_d, "{0}.{1}.base.par".format(base, pst.control_data.noptmax))
        assert os.path.exists(last_par), "missing final base.par for hessian='{0}': {1}".format(update, last_par)

        if update == "true":
            rec_path = os.path.join(m_d, base + ".rec")
            assert os.path.exists(rec_path), "missing rec file for BFGS run: " + rec_path
            hessians = _parse_hessians_from_rec(rec_path)
            assert len(hessians) > 0, "no hessian in rec file -- "
            n = hessians[0].shape[0]
            identity = np.eye(n)
            # assert np.allclose(hessians[0], identity), "initial Hessian is not the identity matrix:\n{}".format(hessians[0])
            final_H = hessians[-1]
            off_diag_diff = np.abs(final_H - identity)
            np.fill_diagonal(off_diag_diff, 0.0)
            diag_diff = np.abs(np.diag(final_H) - 1.0)
            assert off_diag_diff.max() > 1e-10 or diag_diff.max() > 1e-10, "the BFGS update path may not be exercised"


if __name__ == "__main__":
    # basic_sqp_init_test()
    # basic_sqp_iter_test()
    # basic_sqp_direction_test()
    basic_sqp_rosenbrock_chance_test()
    # basic_sqp_chance_test_highdim()
    # basic_sqp_bounds_test()
    # basic_sqp_hessian_test()

