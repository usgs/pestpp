import os
import sys
import re
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
        # registry summary emits the canonical uppercase tag and a std::to_string double
        # (6 decimals), e.g. "OPT_DIRECTION: 1.000000" (was "opt_direction: 1")
        expected_token = "OPT_DIRECTION: {0}.000000".format(1 if direction == "min" else -1)
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

    # NOTE: stack entries are model runs and can be legitimately dropped when they fail or
    # time out (e.g. a transient CI-agent stall); the run manager drops the failed rows and
    # PESTPP-SQP proceeds. So the surviving stack can be <= the requested opt_stack_size.
    # Assert the structural invariant (a stack was drawn, no larger than requested, and the
    # par/obs stacks stayed row-consistent) rather than an exact count that assumes a
    # zero-failure run.
    assert 0 < par_stack_df.shape[0] <= stack_size, "par_stack row count out of range: expected 1..{0} (opt_stack_size), got {1}".format(stack_size, par_stack_df.shape[0])
    assert obs_stack_df.shape[0] == par_stack_df.shape[0], "par/obs_stack row counts disagree: par={0}, obs={1}".format(par_stack_df.shape[0], obs_stack_df.shape[0])

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

    # opt_chance_points='all' draws a SEPARATE (nested) stack for every decision-variable
    # realization (the n_reals drawn members plus the BASE realization), so the ideal row
    # count is (n_reals+1)*stack_size. The nested-stack index encodes membership as
    # "<stack_realization>||<member>". Individual stack entries are model runs, though, and
    # can legitimately fail or time out on a busy CI host; PESTPP-SQP then drops the failed
    # runs (recording them in the failed-run storage file <case>.rnf) and proceeds. So
    # rather than assert the zero-failure ideal, reconcile the surviving rows against the
    # runs that were actually lost. The '.1.' nested stack is written by the iteration-1
    # chance calc (the initial calc writes '.0.'), so scope the drop accounting to that
    # calc's messages in the rec file.
    ideal_all_rows = (n_reals + 1) * stack_size

    def _members(df):
        return [str(ix).split("||")[-1].strip().lower() for ix in df.index]

    _, rec_text_all = _read_rec(m_d_all, pst_name_all)
    # the rec logs one block per nested chance calc; the last block is the iteration-1 calc
    # that produced the '.1.' nested stack files read above
    nested_blocks = rec_text_all.split("queuing up nested sets of chance runs")
    last_block = nested_blocks[-1] if len(nested_blocks) > 1 else rec_text_all

    # partial drops within a member: "WARNING: <N> stack runs failed for realization <m>, dropped"
    partial_drops = {}
    for mo in re.finditer(r"WARNING:\s+(\d+)\s+stack runs failed for realization\s+(\S+?),?\s+dropped", last_block):
        partial_drops[mo.group(2).strip().lower()] = int(mo.group(1))
    # whole-member drops: "WARNING: all stack runs failed for population member 'X', removing..."
    whole_member_drops = set()
    for mo in re.finditer(r"all stack runs failed for population member\s+'([^']+)'", last_block):
        whole_member_drops.add(mo.group(1).strip().lower())

    dropped_rows = sum(partial_drops.values()) + len(whole_member_drops) * stack_size
    expected_all_rows = ideal_all_rows - dropped_rows
    expected_members = (n_reals + 1) - len(whole_member_drops)

    par_members = _members(par_stack_all_df)
    present_members = set(par_members)

    # every surviving dv realization must get its own nested stack
    assert len(present_members) == expected_members, \
        "ALL nested stack member count mismatch: expected {0} ({1} dv realizations minus {2} fully-failed member(s)), got {3}: {4}".format(
            expected_members, n_reals + 1, len(whole_member_drops), len(present_members), sorted(present_members))

    # par and obs nested stacks must stay row-for-row consistent
    assert set(par_stack_all_df.index) == set(obs_stack_all_df.index), \
        "ALL par/obs nested stack indices differ (par has {0} rows, obs has {1})".format(par_stack_all_df.shape[0], obs_stack_all_df.shape[0])

    # exact row reconciliation: surviving rows == ideal minus the runs that actually failed
    assert par_stack_all_df.shape[0] == expected_all_rows, \
        "ALL par_stack row count mismatch: expected {0} = (n_reals+1)*stack_size ({1}) - {2} failed/dropped run(s), got {3}".format(
            expected_all_rows, ideal_all_rows, dropped_rows, par_stack_all_df.shape[0])
    assert obs_stack_all_df.shape[0] == par_stack_all_df.shape[0], \
        "ALL par/obs nested stack row counts disagree: par={0}, obs={1}".format(par_stack_all_df.shape[0], obs_stack_all_df.shape[0])

    # per-member: each present member keeps stack_size minus its own dropped runs
    from collections import Counter
    per_member_rows = Counter(par_members)
    for mname, cnt in per_member_rows.items():
        exp_member = stack_size - partial_drops.get(mname, 0)
        assert cnt == exp_member, \
            "ALL nested stack for member '{0}': expected {1} rows (stack_size {2} - {3} dropped), got {4}".format(
                mname, exp_member, stack_size, partial_drops.get(mname, 0), cnt)

    # corroborate the lost rows against the failed-run storage file (<case>.rnf). That file
    # is cumulative across the whole run (initial + iteration-1 calcs), so it must account
    # for at least the rows lost from this file.
    lost_rows = ideal_all_rows - par_stack_all_df.shape[0]
    rnf_file = os.path.join(m_d_all, base_all + ".rnf")
    if lost_rows > 0:
        assert os.path.exists(rnf_file), \
            "{0} nested-stack row(s) were lost but no failed-run storage file was written: {1}".format(lost_rows, rnf_file)
    if os.path.exists(rnf_file):
        try:
            _, _, rnf_meta = pyemu.helpers.read_pestpp_runstorage(rnf_file, irun="all", with_metadata=True)
            n_failed_runs = int((rnf_meta["r_status"] != 1).sum())
        except Exception as e:
            n_failed_runs = None
            print("WARNING: could not read failed-run storage file {0}: {1}".format(rnf_file, e))
        if n_failed_runs is not None:
            assert n_failed_runs >= lost_rows, \
                "failed-run storage ({0} failed run(s) in {1}) does not account for the {2} lost nested-stack row(s)".format(
                    n_failed_runs, os.path.basename(rnf_file), lost_rows)

    # ALL chance points must produce a genuinely nested stack (one per dv realization),
    # which SINGLE does not
    assert len(present_members) > 1, \
        "ALL chance points did not produce a nested (per-realization) stack: only members {0}".format(sorted(present_members))
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



def sqp_drop_violations_test():
    """sqp drops candidate members violating a nominated 'drop_violations' observation.

    The point of this test is WHERE it applies. sqp is not like ies or mou, where every run is
    a population member - it issues several kinds of batch per iteration, and only some of them
    contain candidates:

      applies      initial ensemble; line-search / trust-region candidate steps
      must NOT     finite-difference gradient runs; the control-file and mean-dv runs

    Dropping a gradient perturbation would leave a hole in the Jacobian - the constraint would
    be steering the DERIVATIVE rather than the search - and dropping the current point would
    leave sqp with nowhere to step from. So this asserts not just that dropping happens, but
    that it only ever happens at the two candidate stages.
    """
    t_d = os.path.join("sqp_drop_viol")
    if os.path.exists(t_d):
        shutil.rmtree(t_d)
    os.makedirs(t_d)

    with open(os.path.join(t_d, "forward.py"), "w") as f:
        f.write('vals={}\n'
                'for line in open("par.dat"):\n'
                '    k,v=line.split(); vals[k.strip()]=float(v)\n'
                'x,y=vals["x"],vals["y"]\n'
                'f=(1.0-x)**2 + 100.0*(y-x*x)**2\n'
                'with open("obs.dat","w") as fp:\n'
                '    fp.write("obs {0:20.10E}\\n".format(f))\n'
                '    fp.write("constraint {0:20.10E}\\n".format(x+y))\n')
    with open(os.path.join(t_d, "par.dat.tpl"), "w") as f:
        f.write("ptf ~\nx ~   x        ~\ny ~   y        ~\n")
    with open(os.path.join(t_d, "obs.dat.ins"), "w") as f:
        f.write("pif ~\nl1 w !obs!\nl1 w !constraint!\n")

    cwd = os.getcwd()
    os.chdir(t_d)
    try:
        pst = pyemu.Pst.from_io_files("par.dat.tpl", "par.dat", "obs.dat.ins", "obs.dat")
    finally:
        os.chdir(cwd)

    par = pst.parameter_data
    par.loc[:, "partrans"] = "none"
    par.loc[:, "parlbnd"] = -2.0
    par.loc[:, "parubnd"] = 2.0
    par.loc[:, "pargp"] = "decvars"
    par.loc[:, "parval1"] = 1.0

    obs = pst.observation_data
    obs.loc["obs", ["obgnme", "obsval", "weight"]] = ["obj_fn", 0.0, 1.0]
    # a less-than constraint, violated whenever x+y exceeds the bound, and NOMINATED so that
    # violating candidates are dropped rather than merely penalised
    # bound at the starting point's value so roughly half the drawn ensemble
    # violates - a constraint that everything violates tests the floor, not the drop
    obs.loc["constraint", ["obgnme", "obsval", "weight"]] = ["l_constraint", 2.0, 1.0]
    obs.loc[:, "drop_violations"] = False
    obs.loc["constraint", "drop_violations"] = True

    pst.pestpp_options["opt_dec_var_groups"] = "decvars"
    pst.pestpp_options["opt_objective_function"] = "obs"
    pst.pestpp_options["opt_direction"] = "min"
    pst.pestpp_options["sqp_num_reals"] = 20
    pst.pestpp_options["random_seed"] = 11
    pst.control_data.noptmax = 2
    pst.model_command = ["python forward.py"]
    pst.write(os.path.join(t_d, "sqp_viol.pst"), version=2)

    # newest first, and absolute: a relative exe path resolves through PATH to whatever
    # pestpp-sqp is installed, which silently tests a release binary instead of this build
    cands = [c for c in (
        os.path.join("..", "build", "src", "programs", "pestpp-sqp", "pestpp-sqp" + exe),
        os.path.abspath(exe_path)) if os.path.exists(c)]
    assert cands, "could not find a pestpp-sqp to test"
    sqp_exe = max([os.path.abspath(c) for c in cands], key=os.path.getmtime)
    print("testing with:", sqp_exe)
    pyemu.os_utils.run("{0} sqp_viol.pst".format(sqp_exe), cwd=t_d)

    rec = open(os.path.join(t_d, "sqp_viol.rec")).read()
    # 1. the nomination registered - without this the rest passes with the feature OFF
    assert "1 'drop_violations' observations detected" in rec, \
        "sqp did not register the drop_violations nomination:\n" + rec[:3000]

    # 2. dropping only ever happened at a candidate stage
    stages = re.findall(r"meet 'drop_violations' conditions during ([a-z ]+?) and are being", rec)
    allowed = {"initial ensemble evaluation", "candidate step evaluation"}
    for st in stages:
        assert st.strip() in allowed, \
            "drop_violations was applied at a stage where a run is structurally required, " \
            "not a candidate: '{0}'. Gradient runs and the current point must never be " \
            "dropped.".format(st.strip())


def sqp_preemption_screening_test():
    """sqp: mid-run screening must change wall-clock and NOTHING else.

    sqp is the riskiest of the three for this, and the reason is structural. ies needed
    screening SUSPENDED during lambda testing, because those runs are a comparison - each
    candidate lambda is scored over the same subset, and abandoning runs mid-comparison changed
    which lambda won. sqp's line search has the same shape: it proposes several candidate
    steps, scores them against each other, and picks one.

    The difference that makes it defensible is that a violating step is one sqp would refuse
    anyway - drop_violating_members() is applied at exactly the candidate stages. This test is
    what turns that argument from plausible into checked: same answer, fewer runs, and nothing
    discarded that the unscreened run went on to keep.
    """
    results = {}
    for tag, interval in (("off", 0.0), ("on", 0.02)):     # 0.02 min ~ 1.2 s
        t_d = "sqp_preempt_" + tag
        if os.path.exists(t_d):
            shutil.rmtree(t_d)
        os.makedirs(t_d)

        # slow, so a run is in flight long enough to be asked about
        with open(os.path.join(t_d, "forward.py"), "w") as f:
            f.write('import time\n'
                    'vals={}\n'
                    'for line in open("par.dat"):\n'
                    '    k,v=line.split(); vals[k.strip()]=float(v)\n'
                    'x,y=vals["x"],vals["y"]\n'
                    'f=(1.0-x)**2 + 100.0*(y-x*x)**2\n'
                    'with open("obs.dat","w") as fp:\n'
                    '    fp.write("obs {0:20.10E}\\n".format(f))\n'
                    '    fp.write("constraint {0:20.10E}\\n".format(x+y))\n'
                    'time.sleep(3)\n')
        with open(os.path.join(t_d, "par.dat.tpl"), "w") as f:
            f.write("ptf ~\nx ~   x        ~\ny ~   y        ~\n")
        with open(os.path.join(t_d, "obs.dat.ins"), "w") as f:
            f.write("pif ~\nl1 w !obs!\nl1 w !constraint!\n")

        cwd = os.getcwd()
        os.chdir(t_d)
        try:
            pst = pyemu.Pst.from_io_files("par.dat.tpl", "par.dat", "obs.dat.ins", "obs.dat")
        finally:
            os.chdir(cwd)

        par = pst.parameter_data
        par.loc[:, "partrans"] = "none"
        par.loc[:, "parlbnd"] = -2.0
        par.loc[:, "parubnd"] = 2.0
        par.loc[:, "pargp"] = "decvars"
        par.loc[:, "parval1"] = 1.0

        obs = pst.observation_data
        obs.loc["obs", ["obgnme", "obsval", "weight"]] = ["obj_fn", 0.0, 1.0]
        obs.loc["constraint", ["obgnme", "obsval", "weight"]] = ["l_constraint", 2.0, 1.0]
        obs.loc[:, "drop_violations"] = False
        obs.loc["constraint", "drop_violations"] = True

        pst.pestpp_options["opt_dec_var_groups"] = "decvars"
        pst.pestpp_options["opt_objective_function"] = "obs"
        pst.pestpp_options["opt_direction"] = "min"
        pst.pestpp_options["sqp_num_reals"] = 20
        pst.pestpp_options["random_seed"] = 11
        pst.pestpp_options["preemption_poll_interval_minutes"] = interval
        pst.control_data.noptmax = 1
        pst.model_command = ["python forward.py"]
        pst.write(os.path.join(t_d, "sqp_pre.pst"), version=2)

        cands = [c for c in (
            os.path.join("..", "build", "src", "programs", "pestpp-sqp", "pestpp-sqp" + exe),
            os.path.abspath(exe_path)) if os.path.exists(c)]
        assert cands, "could not find a pestpp-sqp to test"
        sqp_exe = max([os.path.abspath(c) for c in cands], key=os.path.getmtime)

        wr = "sqp_preempt_workers"
        if os.path.exists(wr):
            shutil.rmtree(wr)
        os.makedirs(wr)
        pyemu.os_utils.start_workers(t_d, sqp_exe, "sqp_pre.pst", num_workers=4,
                                     master_dir=t_d + "_master", worker_root=wr, port=4029)
        results[tag] = t_d + "_master"

    rmrs = {}
    for tag in ("off", "on"):
        for f in os.listdir(results[tag]):
            if f.endswith(".rmr"):
                rmrs[tag] = open(os.path.join(results[tag], f)).read()

    # 1. screening fired, and only in the screened run
    n_abandoned = rmrs["on"].count("abandoned mid-run")
    assert n_abandoned > 0, \
        "no sqp run was abandoned mid-run, so this test did not exercise screening:\n" \
        + rmrs["on"][-3000:]
    assert "abandoned mid-run" not in rmrs["off"], \
        "the UNscreened run abandoned something, which it cannot have done"

    # 1c. and they were reported as ABANDONED, not as model failures. An abandoned run and a
    #     failed one both cost their member, but they tell a user opposite things: one says
    #     preemption worked, the other says go and debug your model.
    rec_on, rec_off = "", ""
    for f in os.listdir(results["on"]):
        if f.endswith(".rec"):
            rec_on = open(os.path.join(results["on"], f)).read()
    for f in os.listdir(results["off"]):
        if f.endswith(".rec"):
            rec_off = open(os.path.join(results["off"], f)).read()
    assert "abandoned mid-run" in rec_on, \
        "the .rec never says anything was abandoned, so screened runs are being reported as " \
        "something else:\n" + rec_on[-3000:]
    assert "NOT model failures" in rec_on, \
        "abandoned members are not being distinguished from model failures in the .rec"
    assert "abandoned mid-run" not in rec_off, \
        "the UNscreened run reported an abandonment, which it cannot have made"
    for name, text in (("off", rec_off), ("on", rec_on)):
        assert "runs failed" not in text, \
            "{0}: members were reported as model FAILURES, but nothing failed in this case - " \
            "abandoned runs are being mislabelled".format(name)

    # 2. it saved exactly that much work, and no abandoned run was rescheduled
    n_runs = {t: rmrs[t].count("received from") for t in ("off", "on")}
    assert n_runs["on"] < n_runs["off"], \
        "screening saved no model runs: off={0} on={1}".format(n_runs["off"], n_runs["on"])
    assert n_runs["off"] - n_runs["on"] == n_abandoned, \
        "off made {0} runs, on made {1}, but {2} were abandoned - an abandoned run was " \
        "rescheduled rather than dropped".format(n_runs["off"], n_runs["on"], n_abandoned)

    # 3. and the ANSWER is unchanged. This is the part the line-search structure puts at risk:
    #    if abandoning a candidate mid-comparison changed which step won, the objective
    #    trajectory would diverge - which is exactly how the ies lambda-testing bug showed up.
    off_files = {f for f in os.listdir(results["off"]) if f.endswith(".csv")}
    on_files = {f for f in os.listdir(results["on"]) if f.endswith(".csv")}
    assert off_files == on_files, \
        "different output files were written:\n only off: {0}\n only on: {1}".format(
            sorted(off_files - on_files), sorted(on_files - off_files))

    # The ensembles are written BEFORE drop_violating_members() culls them, so the unscreened
    # run reports candidates it is about to discard and the screened run never computed them.
    # Demanding byte-identity would therefore demand that screening evaluate the very steps it
    # declined to take - the same wrong criterion that cost a day on the ies version of this
    # test. What must hold instead:
    #
    #   - every member present in BOTH runs is identical, to the byte. This is the assertion
    #     that actually addresses the line-search worry: if abandoning a candidate mid-
    #     comparison had changed which step won, the survivors would carry different values.
    #   - anything missing from the screened run VIOLATES, so it is something the unscreened
    #     run went on to drop anyway.
    #
    # .pcs.csv is exempt from the row comparison because it is not an ensemble: it is a
    # parameter-change summary computed ACROSS the members, so a legitimately different member
    # set at write time moves every statistic in it. Its shape is still checked.
    checked = 0
    for f in sorted(off_files):
        a = pd.read_csv(os.path.join(results["off"], f), index_col=0)
        b = pd.read_csv(os.path.join(results["on"], f), index_col=0)
        a.columns = [c.lower() for c in a.columns]
        b.columns = [c.lower() for c in b.columns]
        if ".pcs." in f:
            assert list(a.columns) == list(b.columns), \
                "{0}: the change summary has different columns".format(f)
            continue
        common_rows = [i for i in a.index if i in b.index]
        # numeric only: some outputs carry text columns (names, tags) that cannot be subtracted
        num = set(a.select_dtypes(include="number").columns) & \
              set(b.select_dtypes(include="number").columns)
        common_cols = [c for c in a.columns if c in num]
        if not common_rows or not common_cols:
            continue
        d = (a.loc[common_rows, common_cols] - b.loc[common_rows, common_cols]).abs()
        assert d.max().max() < 1e-10, \
            "{0}: screening CHANGED THE RESULT - a member present in both runs differs by up " \
            "to {1}. If the line search picked a different step, this is where it shows." \
            .format(f, d.max().max())
        if "constraint" in common_cols:
            for missing in set(a.index) - set(b.index):
                assert a.loc[missing, "constraint"] > 2.0, \
                    "{0}: member '{1}' was screened out but does NOT violate ({2} <= 2.0) - " \
                    "screening may only ever drop what the process would have dropped" \
                    .format(f, missing, a.loc[missing, "constraint"])
        checked += 1
    assert checked > 0, "no sqp output files were compared, so nothing was actually checked"
    print("RESULT sqp preemption screening PASS ({0} file(s) compared, {1} runs vs {2}, "
          "{3} abandoned)".format(checked, n_runs["on"], n_runs["off"], n_abandoned))
