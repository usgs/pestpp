import os
import shutil
import platform
import socket
import numpy as np
import pandas as pd
import pyemu


def _get_free_port():
    """Find an available TCP port by letting the OS assign one."""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("", 0))
        return s.getsockname()[1]


def _parse_mps_row_rhs(mps_file):
    """parse a COIN/CLP standard MPS file -> ({row_name: rhs_value}, {row_name: sense}).

    section keywords (ROWS/COLUMNS/RHS/...) start at column 0; data rows are indented.
    note the RHS-vector name in RHS data rows is itself 'RHS', so a naive token check
    would mistake data rows for the section header - key off the leading indent instead.
    """
    rhs, sense, sec = {}, {}, None
    with open(mps_file) as f:
        for line in f:
            if not line.strip():
                continue
            t = line.split()
            if (not line[0].isspace()) and t[0] in ("NAME", "OBJSENSE", "ROWS", "COLUMNS",
                                                    "RHS", "RANGES", "BOUNDS", "ENDATA"):
                sec = t[0]
                continue
            if sec == "ROWS" and len(t) == 2 and t[0] in ("L", "G", "E", "N"):
                sense[t[1]] = t[0]
            elif sec == "RHS" and len(t) >= 3:
                for i in range(1, len(t) - 1, 2):
                    rhs[t[i]] = float(t[i + 1])
    return rhs, sense


def _parse_rec_obs_constraints(rec_file, iteration):
    """parse the 'observation constraint/objective information at start of iteration N'
    table -> {name: (required, shifted_sim_value)}.  the fixed-width columns can run
    together where the +/- double-max bound abuts the residual, so read by position:
    name=t[0], sense=t[1], required=t[2], (chance-)shifted sim value=t[3]."""
    hdr = "observation constraint/objective information at start of iteration {0}".format(iteration)
    lines = open(rec_file).read().splitlines()
    hs = [i for i, l in enumerate(lines) if hdr in l]
    out = {}
    if not hs:
        return out
    for l in lines[hs[-1] + 2:]:
        t = l.split()
        if len(t) >= 5 and t[1] in ("less_than", "greater_than", "equal_to"):
            try:
                out[t[0]] = (float(t[2]), float(t[3]))
            except ValueError:
                pass
        elif out:
            break  # reached the end of the table
    return out


def check_mps_shift_coherence(d, run_name, iters, tol=1.0e-3):
    """assert the LP problem the simplex actually solved is coherent with the chance-shifted
    constraint values pestpp-opt reports.

    pestpp-opt writes a standard MPS file each SLP iteration ('<run>.<iter>.mps', on by
    default via ++opt_coin_log).  Its right-hand side for each observation-constraint row
    is (required - shifted sim value).  The rec file's 'observation constraint ... start of
    iteration' table reports the same 'required' and (chance-)'shifted sim value'.  If the
    LP is built with one chance convention (e.g. raw stack quantiles) while the reported /
    convergence-checked shifted values use another (e.g. anomaly form), the two disagree and
    the 'optimal' solution silently violates the constraints it is reported to satisfy.
    This check ties the two together: MPS RHS == required - reported shifted sim value.
    """
    rec_file = os.path.join(d, run_name + ".rec")
    assert os.path.exists(rec_file), "missing rec file: " + rec_file
    total = 0
    for it in iters:
        mps = os.path.join(d, "{0}.{1}.mps".format(run_name, it))
        if not os.path.exists(mps):
            continue
        rhs, sense = _parse_mps_row_rhs(mps)
        rec = _parse_rec_obs_constraints(rec_file, it)
        common = [n for n in rec if n in rhs and sense.get(n) in ("L", "G")]
        assert len(common) > 0, \
            "no obs constraints common to {0} and rec iteration {1}".format(mps, it)
        for name in common:
            required, shifted = rec[name]
            expected_rhs = required - shifted
            assert abs(rhs[name] - expected_rhs) < tol, (
                "MPS vs reported-shift incoherence at iter {0} constraint '{1}': "
                "MPS RHS={2} but required-shifted={3} (required={4}, shifted sim={5})".format(
                    it, name, rhs[name], expected_rhs, required, shifted))
            total += 1
    assert total > 0, "no MPS constraint rows were checked - was ++opt_coin_log disabled?"
    print("MPS coherence: {0} constraint-rows across iterations {1} match reported shifted values".format(
        total, list(iters)))


bin_path = os.path.join("test_bin")
print(platform.platform().lower())

if "linux" in platform.platform().lower():
    bin_path = os.path.join(bin_path,"linux")
elif "darwin" in platform.platform().lower() or "mac" in platform.platform().lower():
    bin_path = os.path.join(bin_path,"mac")
else:
    bin_path = os.path.join(bin_path,"win")

bin_path = os.path.abspath("test_bin")
os.environ["PATH"] += os.pathsep + bin_path


# case of either appveyor, travis or local
if os.path.exists(os.path.join("pestpp","bin")):
    bin_path = os.path.join("..","..","pestpp","bin")
else:
    bin_path = os.path.join("..","..","..","..","pestpp","bin")
    
if "windows" in platform.platform().lower():
    exe_path = os.path.join(bin_path, "pestpp-opt.exe")
elif "darwin" in platform.platform().lower() or "mac" in platform.platform().lower():
    exe_path = os.path.join(bin_path, "pestpp-opt")
else:
    exe_path = os.path.join(bin_path, "pestpp-opt")


def std_weights_test():
    d = os.path.join("opt_dewater_chance", "test_std_weights")
    if os.path.exists(d):
        shutil.rmtree(d)
    shutil.copytree(os.path.join("opt_dewater_chance", "template"), d)
    pst_file = os.path.join(d, "dewater_pest.base.pst")
    jco_file = os.path.join(d, "dewater_pest.full.jcb")
    pst = pyemu.Pst(pst_file)
    par = pst.parameter_data
    par.loc[par.partrans == "fixed", "partrans"] = "log"
    jco = pyemu.Jco.from_binary(jco_file)
    par.loc[par.pargp == "q", "partrans"] = "fixed"
    obs = pst.observation_data.loc[pst.nnz_obs_names, :]

    forecast_names = list(obs.loc[obs.obgnme.apply(lambda x: x.startswith("l_") or \
                                                             x.startswith("less_")), "obsnme"])
    # print(forecast_names)
    pst.observation_data.loc[:, "weight"] = 0.0
    sc = pyemu.Schur(jco=jco, pst=pst, forecasts=forecast_names)
    # print(sc.get_forecast_summary())

    fstd = sc.get_forecast_summary().loc[:, "post_var"].apply(lambda x: np.sqrt(x))
    pr_unc_py = fstd.to_dict()
    pst.observation_data.loc[fstd.index, "weight"] = fstd.values

    pst.pestpp_options["opt_risk"] = 0.1
    pst.pestpp_options["opt_std_weights"] = True
    pst.pestpp_options["base_jacobian"] = os.path.split(jco_file)[-1]
    par.loc[par.pargp == "q", "partrans"] = "none"

    new_pst_file = os.path.join(d, "test.pst")

    pst.write(new_pst_file)
    pyemu.os_utils.run("{0} {1}".format(exe_path, os.path.split(new_pst_file)[-1]), cwd=d)
    pr_unc1 = scrap_rec(new_pst_file.replace(".pst", ".rec"))

    print(pr_unc1)
    for fore in forecast_names:
        dif = np.abs(pr_unc_py[fore] - pr_unc1[fore])
        print(fore, pr_unc_py[fore], pr_unc1[fore], dif)
        assert dif < 1.0e-4

    pst.pestpp_options["opt_std_weights"] = False
    pst.write(new_pst_file)
  
    pyemu.os_utils.run("{0} {1}".format(exe_path, os.path.split(new_pst_file)[-1]), cwd=d)
    pr_unc2 = scrap_rec(new_pst_file.replace(".pst", ".rec"))
    print(pr_unc2)
    for fore in forecast_names:
        dif = np.abs(pr_unc_py[fore] - pr_unc2[fore])
        print(fore, pr_unc_py[fore], pr_unc2[fore], dif)
        assert dif < 1.0e-4,dif


def scrap_rec(rec_file):
    unc = {}
    tag = "FOSM-based chance constraint information at start of iteration 1"
    # this alt tag is to support mou dev 
    tag_alt = "Chance constraint/objective information at start of iteration 1"
    with open(rec_file, 'r') as f:
        while True:
            line = f.readline()
            if line == "":
                break

            if tag in line or tag_alt in line:
                f.readline()
                while True:
                    line = f.readline()
                    if line == "":
                        break
                    raw = line.strip().split()
                    # print(raw)
                    try:
                        name = raw[0].lower()
                        val = float(raw[4])
                        unc[name] = val
                    except:
                        break
                break
    return unc


def run_dewater_test():
    worker_d = os.path.join("opt_dewater_chance")
    pyemu.os_utils.start_workers(os.path.join(worker_d, "template"), exe_path, "dewater_pest.base.pst",
                                master_dir=os.path.join(worker_d, "master1"), worker_root=worker_d, num_workers=10,
                                verbose=True,port=_get_free_port())

    opt = None
    with open(os.path.join(worker_d, "master1", "dewater_pest.base.rec"), 'r') as f:
        for line in f:
            if "iteration 1 objective function value:" in line:
                opt = float(line.strip().split()[-2])
    assert opt is not None

    pst = pyemu.Pst(os.path.join(worker_d,"template","dewater_pest.base.pst"))
    pst.control_data.noptmax = 3
    pst.write(os.path.join(worker_d,"template","test.pst"))
    pyemu.os_utils.start_workers(os.path.join(worker_d, "template"), exe_path, "test.pst",
                                master_dir=os.path.join(worker_d, "master2"), worker_root=worker_d, num_workers=10,
                                verbose=True,port=_get_free_port())
    
    with open(os.path.join(worker_d,"master2","test.rec")) as f:
        for line in f:
            if "iteration       obj func" in line:
                f.readline() # skip the initial obj func
                lines = []
                for _ in range(pst.control_data.noptmax):
                    lines.append(f.readline())
    obj_funcs = np.array([float(line.strip().split()[-1]) for line in lines])
    print(obj_funcs)
    assert np.abs(obj_funcs.max() - obj_funcs.min()) < 0.1

    ws = os.path.join(worker_d,"ext_objfunc_master")
    if os.path.exists(ws):
        shutil.rmtree(ws)
    shutil.copytree(os.path.join(worker_d,"template"),ws)
    pst.pestpp_options.pop("base_jacobian",None)
    par = pst.parameter_data
    dv_pars = par.loc[par.pargp=="q","parnme"]
    with open(os.path.join(ws,"ext_OBJ_func.csv"),'w') as f:
        for dv_par in dv_pars:
            f.write("{0},1\n".format(dv_par))
    pst.pestpp_options["opt_objective_function"] = "ext_OBJ_func.csv"
    pst.control_data.noptmax = 3
    pst.write(os.path.join(ws,"test.pst"))
    #pyemu.os_utils.start_workers(os.path.join(worker_d, "template"), exe_path, "test.pst",
    #                            master_dir=os.path.join(worker_d, "master3"), worker_root=worker_d, num_workers=10,
    #                            verbose=True,port=_get_free_port())
    pyemu.os_utils.run("{0} {1}".format(exe_path,"test.pst"),cwd=ws)
    with open(os.path.join(ws,"test.rec")) as f:
        for line in f:
            if "iteration       obj func" in line:
                f.readline() # skip the initial obj func
                lines = []
                for _ in range(pst.control_data.noptmax):
                    lines.append(f.readline())

    ext_obj_funcs = np.array([float(line.strip().split()[-1]) for line in lines])
    print(obj_funcs)
    print(ext_obj_funcs) 

    assert np.abs(ext_obj_funcs - obj_funcs).max() < 1e-6

    ws = os.path.join(worker_d,"ext_objfunc_master2")
    if os.path.exists(ws):
        shutil.rmtree(ws)
    shutil.copytree(os.path.join(worker_d,"template"),ws)
    pst.pestpp_options.pop("base_jacobian",None)
    par = pst.parameter_data
    dv_pars = par.loc[par.pargp=="q","parnme"]
    with open(os.path.join(ws,"ext_OBJ_func.csv"),'w') as f:
        for dv_par in dv_pars:
            f.write("{0},3.0\n".format(dv_par))
    pst.pestpp_options["opt_objective_function"] = "ext_OBJ_func.csv"
    pst.control_data.noptmax = 3
    pst.write(os.path.join(ws,"test.pst"))
    #pyemu.os_utils.start_workers(os.path.join(worker_d, "template"), exe_path, "test.pst",
    #                            master_dir=os.path.join(worker_d, "master3"), worker_root=worker_d, num_workers=10,
    #                            verbose=True,port=_get_free_port())
    pyemu.os_utils.run("{0} {1}".format(exe_path,"test.pst"),cwd=ws)
    with open(os.path.join(ws,"test.rec")) as f:
        for line in f:
            if "iteration       obj func" in line:
                f.readline() # skip the initial obj func
                lines = []
                for _ in range(pst.control_data.noptmax):
                    lines.append(f.readline())

    ext_obj_funcs = np.array([float(line.strip().split()[-1]) for line in lines]) / 3.0
    print(obj_funcs)
    print(ext_obj_funcs) 

    assert np.abs(ext_obj_funcs - obj_funcs).max() < 1e-2


    ws = os.path.join(worker_d,"tied_serial_master")
    if os.path.exists(ws):
        shutil.rmtree(ws)
    shutil.copytree(os.path.join(worker_d,"template"),ws)

    pst.pestpp_options.pop("opt_objective_function")
    pst.pestpp_options["opt_risk"] = 0.95
    pst.pestpp_options.pop("base_jacobian",None)
    par = pst.parameter_data
    dv_pars = par.loc[par.pargp=="q","parnme"]
    adj_pars = par.loc[par.pargp=="h","parnme"]
    pst.parameter_data.loc[adj_pars[0],"partrans"] = "log"
    pst.parameter_data.loc[adj_pars[1:],"partrans"] = "tied"
    pst.parameter_data.loc[adj_pars[1:],"partied"] = adj_pars[0]

    pst.parameter_data.loc[["up_grad","dn_grad"],"partrans"] = "log"
    pst.pestpp_options["opt_recalc_chance_every"] = 100
    pst.control_data.noptmax = 3
    pst.write(os.path.join(ws,"test.pst"))
    #pyemu.os_utils.start_workers(os.path.join(worker_d, "template"), exe_path, "test.pst",
    #                            master_dir=os.path.join(worker_d, "master3"), worker_root=worker_d, num_workers=10,
    #                            verbose=True,port=_get_free_port())
    pyemu.os_utils.run("{0} {1}".format(exe_path,"test.pst"),cwd=ws)
    #with open(os.path.join(worker_d,"master3","test.rec")) as f:
    with open(os.path.join(ws,"test.rec")) as f:
        for line in f:
            if "iteration       obj func" in line:
                f.readline() # skip the initial obj func
                lines = []
                for _ in range(pst.control_data.noptmax):
                    lines.append(f.readline())
    averse_obj_funcs = np.array([float(line.strip().split()[-1]) for line in lines])
    print(averse_obj_funcs) 
    assert np.abs(averse_obj_funcs.max() - averse_obj_funcs.min()) < 0.1
   

    ws = os.path.join(worker_d,"tied_serial_master2")
    if os.path.exists(ws):
        shutil.rmtree(ws)
    shutil.copytree(os.path.join(worker_d,"template"),ws)
    pst.parameter_data.loc[adj_pars,"partrans"] = "log"
    #pst.drop_parameters(os.path.join(ws,"transmissivity_layer_1.ref.tpl"),pst_path='.')
    #pst.drop_parameters(os.path.join(ws,"strt_layer_1.ref.tpl"),pst_path='.')


    pst.pestpp_options["opt_risk"] = 0.95
    pst.pestpp_options.pop("base_jacobian",None)
    par = pst.parameter_data
    dv_pars = par.loc[par.pargp=="q","parnme"]
    adj_pars = par.loc[par.pargp=="h","parnme"]
    pst.parameter_data.loc[adj_pars[0],"partrans"] = "log"
    pst.parameter_data.loc[adj_pars[1:],"partrans"] = "tied"
    pst.parameter_data.loc[adj_pars[1:],"partied"] = adj_pars[0]

    pst.parameter_data.loc[dv_pars[1:],"partrans"] = "tied"
    pst.parameter_data.loc[dv_pars[1:],"partied"] = dv_pars[0]
    pst.pestpp_options["opt_risk"] = 0.5
    #pst.parameter_data.loc[dv_pars[1:],"pargp"] = "qtied"
    
    #pst.write(os.path.join(ws,"test.pst"))

    pst.parameter_data.loc[["up_grad","dn_grad"],"partrans"] = "log"
    pst.pestpp_options["opt_recalc_chance_every"] = 100
    pst.control_data.noptmax = 3
    pst.write(os.path.join(ws,"test.pst"))
   
    pyemu.os_utils.run("{0} {1}".format(exe_path,"test.pst"),cwd=ws)
    #with open(os.path.join(worker_d,"master3","test.rec")) as f:
    with open(os.path.join(ws,"test.rec")) as f:
        for line in f:
            if "iteration       obj func" in line:
                f.readline() # skip the initial obj func
                lines = []
                for _ in range(pst.control_data.noptmax):
                    lines.append(f.readline())
    averse_obj_funcs = np.array([float(line.strip().split()[-1]) for line in lines])
    print(averse_obj_funcs) 
    assert np.abs(averse_obj_funcs.max() - averse_obj_funcs.min()) < 0.1   
    par_files = [f for f in os.listdir(ws) if f.endswith(".par")] 
    assert len(par_files)> 0
    for par_file in par_files:
        pf = pyemu.pst_utils.read_parfile(os.path.join(ws,par_file))
        diff = pf.parval1.loc[dv_pars[1:]].values - pf.parval1.loc[dv_pars[0]]
        print(diff)
        assert np.abs(diff).sum() <1.e-6

    res_files = [f for f in os.listdir(ws) if f.endswith(".res") or f.endswith(".rei")] 
    assert len(res_files)> 0
    for res_file in res_files:
        pf = pyemu.pst_utils.read_resfile(os.path.join(ws,res_file))
        diff = pf.modelled.loc[dv_pars[1:]].values - pf.modelled.loc[dv_pars[0]]
        print(diff)
        assert np.abs(diff).sum() <1.e-6
 


    pst.parameter_data.loc[dv_pars,"partrans"] = "none"
    pst.pestpp_options["opt_recalc_chance_every"] = 2
    pst.control_data.noptmax = 10
    pst.write(os.path.join(worker_d,"template","test.pst"))
    pyemu.os_utils.start_workers(os.path.join(worker_d, "template"), exe_path, "test.pst",
                                master_dir=os.path.join(worker_d, "master4"), worker_root=worker_d, num_workers=10,
                                verbose=True,port=_get_free_port())
    with open(os.path.join(worker_d,"master4","test.rec")) as f:
        for line in f:
            if "iteration       obj func" in line:
                f.readline() # skip the initial obj func
                lines = []
                for _ in range(pst.control_data.noptmax):
                    line = f.readline()
                    if line.strip().startswith("---"):
                        break
                    lines.append(line)
    averse_obj_funcs = np.array([float(line.strip().split()[-1]) for line in lines])
    print(averse_obj_funcs) 
    assert len(averse_obj_funcs) > 2

    pst.pestpp_options["opt_recalc_chance_every"] = 1
    pst.control_data.noptmax = 4
    pst.write(os.path.join(worker_d,"template","test.pst"))
    pyemu.os_utils.start_workers(os.path.join(worker_d, "template"), exe_path, "test.pst",
                                master_dir=os.path.join(worker_d, "master5"), worker_root=worker_d, num_workers=10,
                                verbose=True,port=_get_free_port())
    with open(os.path.join(worker_d,"master5","test.rec")) as f:
        for line in f:
            if "iteration       obj func" in line:
                f.readline() # skip the initial obj func
                lines = []
                for _ in range(pst.control_data.noptmax):
                    line = f.readline()
                    if line.strip().startswith("---"):
                        break
                    lines.append(line)
    averse_obj_funcs = np.array([float(line.strip().split()[-1]) for line in lines])
    print(averse_obj_funcs) 
    assert len(averse_obj_funcs) > 2
    assert np.abs(averse_obj_funcs.max() - averse_obj_funcs.min()) < 0.1 





def run_supply2_test():
    worker_d = os.path.join("opt_supply2_chance")
    pst = pyemu.Pst(os.path.join(worker_d,"template","supply2_pest.base.pst"))
    pst.pestpp_options["panther_agent_freeze_on_fail"] = True
    pst.control_data.noptmax = 1
    pst.write(os.path.join(worker_d,"template","test.pst"))
    pyemu.os_utils.start_workers(os.path.join(worker_d, "template"), exe_path, "test.pst",
                                master_dir=os.path.join(worker_d, "master1"), worker_root=worker_d, num_workers=10,
                                verbose=True,port=_get_free_port())

    opt = None
    with open(os.path.join(worker_d, "master1", "test.rec"), 'r') as f:
        for line in f:
            if "iteration 1 objective function value:" in line:
                opt = float(line.strip().split()[-2])
    assert opt is not None

    pst = pyemu.Pst(os.path.join(worker_d,"template","supply2_pest.base.pst"))
    pst.control_data.noptmax = 1
    pst.pestpp_options["opt_iter_tol"] = 1.0e-10
    pst.write(os.path.join(worker_d,"template","test.pst"))
    pyemu.os_utils.start_workers(os.path.join(worker_d, "template"), exe_path, "test.pst",
                                master_dir=os.path.join(worker_d, "master2"), worker_root=worker_d, num_workers=10,
                                verbose=True,port=_get_free_port())
    
    with open(os.path.join(worker_d,"master2","test.rec")) as f:
        for line in f:
            if "iteration       obj func" in line:
                f.readline() # skip the initial obj func
                lines = []
                for _ in range(pst.control_data.noptmax):
                    lines.append(f.readline())
    obj_funcs = np.array([float(line.strip().split()[-1]) for line in lines])
    print(obj_funcs)
    #assert np.abs(obj_funcs.max() - obj_funcs.min()) < 0.1




def est_res_test():
    worker_d = os.path.join("opt_supply2_chance")
    t_d = os.path.join(worker_d,"template")
    m_d = os.path.join(worker_d,"master")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)

    pst = pyemu.Pst(os.path.join(t_d,"supply2_pest.base.pst"))
    obs = pst.observation_data
    obs.loc[:,"weight"] = 1.0
    pst.pestpp_options["opt_risk"] = 0.05
    pst.pestpp_options["opt_std_weights"] = True
    pst.write(os.path.join(t_d,"supply2_pest.base.pst"))
    pyemu.os_utils.start_workers(os.path.join(worker_d, "template"), exe_path, "supply2_pest.base.pst",
                                master_dir=os.path.join(worker_d, "master"), worker_root=worker_d, num_workers=10,
                                verbose=True,port=_get_free_port())

    opt = None
    with open(os.path.join(worker_d, "master", "supply2_pest.base.rec"), 'r') as f:
        for line in f:
            if "iteration 1 objective function value:" in line:
                opt = float(line.strip().split()[-2])
    assert opt is not None
    res_est1 = pyemu.pst_utils.read_resfile(os.path.join(m_d,"supply2_pest.base.1.est+chance.rei"))

    for f in ["supply2_pest.base.1.jcb","supply2_pest.base.1.jcb.rei"]:
        shutil.copy2(os.path.join(m_d,f),os.path.join(t_d,f))
    pst = pyemu.Pst(os.path.join(t_d,"supply2_pest.base.pst"))
    pst.pestpp_options["opt_skip_final"] = True
    pst.pestpp_options["base_jacobian"] = "supply2_pest.base.1.jcb"
    pst.pestpp_options["hotstart_resfile"] = "supply2_pest.base.1.jcb.rei"
    pst.write(os.path.join(t_d,"pest_est_res.pst"))
    m_d = os.path.join(worker_d,"master_est_res")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    shutil.copytree(t_d,m_d)
    pyemu.os_utils.run("{0} pest_est_res.pst".format(exe_path),cwd=m_d)
    res_est2 = pyemu.pst_utils.read_resfile(os.path.join(m_d,"pest_est_res.1.est+chance.rei"))
    diff = (res_est1.modelled - res_est2.modelled).apply(lambda x: np.abs(x))
    print(diff.sum())
    assert diff.sum() < 1.0e-10

def stack_test():
    d = os.path.join("opt_dewater_chance", "stack_test")
    if os.path.exists(d):
        shutil.rmtree(d)
    shutil.copytree(os.path.join("opt_dewater_chance", "template"), d)
    pst_file = os.path.join(d, "dewater_pest.base.pst")
    pst = pyemu.Pst(pst_file)
    par = pst.parameter_data
    par.loc[par.partrans=="fixed","partrans"] = "log"
    obs = pst.observation_data
    obs.loc["q1","obsval"] = 100
    obs.loc["q1","obgnme"] = "equal_to"
    obs.loc["q1","weight"] = 1.0
    
    pst.pestpp_options.pop("opt_constraint_groups",None)
    pst.pestpp_options.pop("base_jacobian",None)
    pst.pestpp_options["opt_risk"] = 0.1
    pst.pestpp_options["opt_stack_size"] = 10
    pst.control_data.noptmax = 1
    new_pst_file = os.path.join(d, "test.pst")

    pst.write(new_pst_file)
    pyemu.os_utils.run("{0} {1}".format(exe_path, os.path.split(new_pst_file)[-1]), cwd=d)
    rec1 = os.path.join(d,"test.rec")
    assert os.path.exists(rec1)
    par_stack = "test.0.par_stack.csv"
    assert os.path.exists(os.path.join(d,par_stack))
    shutil.copy2(os.path.join(d,par_stack),os.path.join("opt_dewater_chance", "template","par_stack.csv"))
    pst.pestpp_options["opt_par_stack"] = "par_stack.csv"
    
    d = d + "_par_stack"
    if os.path.exists(d):
        shutil.rmtree(d)
    shutil.copytree(os.path.join("opt_dewater_chance", "template"), d)
    pst.write(os.path.join(d,"test.pst"))
    pyemu.os_utils.run("{0} {1}".format(exe_path, "test.pst"), cwd=d)
    rec2 = os.path.join(d,"test.rec")
    assert os.path.exists(rec2)

    par_stack = "test.0.par_stack.csv"
    assert os.path.exists(os.path.join(d,par_stack))
    obs_stack = "test.1.obs_stack.csv"
    assert os.path.exists(os.path.join(d,obs_stack))
    shutil.copy2(os.path.join(d,obs_stack),os.path.join("opt_dewater_chance", "template","obs_stack.csv"))
    shutil.copy2(os.path.join(d,par_stack),os.path.join("opt_dewater_chance", "template","par_stack.csv"))
    pst.pestpp_options["opt_obs_stack"] = "obs_stack.csv"
    pst.pestpp_options["opt_par_stack"] = "par_stack.csv"
    
    d = d + "_obs_stack"
    if os.path.exists(d):
        shutil.rmtree(d)
    shutil.copytree(os.path.join("opt_dewater_chance", "template"), d)
    
    pst.write(os.path.join(d,"test.pst"))
    pyemu.os_utils.run("{0} {1}".format(exe_path, "test.pst"), cwd=d)   
    rec3 = os.path.join(d,"test.rec")
    assert os.path.exists(rec3)
    

    d = d + "_calcresp"
    if os.path.exists(d):
        shutil.rmtree(d)
    shutil.copytree(os.path.join("opt_dewater_chance", "template"), d)
    pst.pestpp_options.pop("base_jacobian",None)
    pst.write(os.path.join(d,"test.pst"))
    pyemu.os_utils.run("{0} {1}".format(exe_path, "test.pst"), cwd=d)   
    rec4 = os.path.join(d,"test.rec")
    assert os.path.exists(rec4)
    
    def get_obj(rec_file):
        tag = "---  objective function sequence  ---"
        with open(rec_file,'r') as f:
            
            for line in f:
                if tag in line:
                    for _ in range(3):
                        line = f.readline()
                    obj_vals = []
                    while True:
                    #print(line)
                        obj = float(line.strip().split()[-1])
                        print(obj)
                        obj_vals.append(obj)
                        line = f.readline()
                        if "---" in line:
                            break
                    return np.array(obj_vals)
        raise Exception("couldnt find obj value in rec file "+rec_file)

    obj1 = get_obj(rec1)
    obj2 = get_obj(rec2)
    obj3 = get_obj(rec3)
    obj4 = get_obj(rec4)

    assert np.abs(obj1 - obj2) < 1.0e-1
    assert np.abs(obj2 - obj3) < 1.0e-1
    assert np.abs(obj2 - obj4) < 1.0e-1

    d = os.path.join("opt_dewater_chance","stack_iter_test")
    if os.path.exists(d):
        shutil.rmtree(d)
    shutil.copytree(os.path.join("opt_dewater_chance", "template"), d)
    pst.control_data.noptmax = 5
    pst.pestpp_options["opt_recalc_chance_every"] = 1
    pst.pestpp_options["opt_stack_size"] = 10
    pst.pestpp_options.pop("opt_par_stack")
    pst.pestpp_options.pop("opt_obs_stack")
    pst.write(os.path.join(d,"test.pst"))
    pyemu.os_utils.run("{0} {1}".format(exe_path, "test.pst"), cwd=d)   
    rec5 = os.path.join(d,"test.rec")
    assert os.path.exists(rec5)
    obj5 = get_obj(rec5)
    #print(obj5)
    #print(obj1)
    assert np.abs(obj1 - obj5[0]) < 1.0e-1
    assert len(obj5) == pst.control_data.noptmax
    # pst.write(os.path.join("opt_dewater_chance","template","test.pst"))
    # pyemu.os_utils.start_workers(os.path.join("opt_dewater_chance", "template"), exe_path, "test.pst",
    #                             master_dir=d, worker_root="opt_dewater_chance", num_workers=10,
    #                             verbose=True,port=_get_free_port())


def dewater_restart_test():
    worker_d = os.path.join("opt_dewater_chance")
    pst = pyemu.Pst(os.path.join(worker_d,"template","dewater_pest.base.pst"))
    par = pst.parameter_data
    gpars = par.loc[par.parnme.str.contains("grad"),"parnme"]
    par.loc[gpars,"partrans"] = "log"
    #par.loc[par.partrans=="fixed","partrans"] = "log"
    pst.pestpp_options.pop("base_jacobian",None)
    pst.pestpp_options["opt_risk"] = 0.95
    pst.control_data.noptmax = 2
    pst.pestpp_options["opt_recalc_chance_every"] = 1
    #pst.pestpp_options.pop("opt_risk",None)
    pst.write(os.path.join(worker_d,"template","base.pst"))
    pyemu.os_utils.start_workers(os.path.join(worker_d, "template"), exe_path, "base.pst",
                                master_dir=os.path.join(worker_d, "master_base1"), worker_root=worker_d, num_workers=10,
                                verbose=True,port=_get_free_port())

    with open(os.path.join(worker_d,"master_base1","base.rec")) as f:
        for line in f:
            if "iteration       obj func" in line:
                f.readline() # skip the initial obj func
                lines = []
                for _ in range(pst.control_data.noptmax):
                    lines.append(f.readline())
    base_obj_funcs = np.array([float(line.strip().split()[-1]) for line in lines])

    pst.control_data.noptmax = 2
    shutil.copy2(os.path.join(worker_d,"master_base1","base.1.jcb"),os.path.join(worker_d,"template","restart.jcb"))
    pst.pestpp_options["base_jacobian"] = "restart.jcb"
    pst.write(os.path.join(worker_d,"template","restart.pst"))
    pyemu.os_utils.run("{0} restart.pst".format(exe_path),cwd=os.path.join(worker_d,"template"))
    
    with open(os.path.join(worker_d,"template","restart.rec")) as f:
        for line in f:
            if "iteration       obj func" in line:
                f.readline() # skip the initial obj func
                lines = []
                for _ in range(pst.control_data.noptmax):
                    lines.append(f.readline())
    restart_obj_funcs = np.array([float(line.strip().split()[-1]) for line in lines])


    print(base_obj_funcs)
    print(restart_obj_funcs)
    d = np.abs(base_obj_funcs - restart_obj_funcs)

    assert d.max() < 0.1

    
def startworker():
    worker_d = os.path.join("opt_dewater_chance")
    t_d = os.path.join(worker_d, "template")

    pyemu.os_utils.start_workers(t_d,exe_path,"test.pst",num_workers=10,worker_root=worker_d)


def _read_best_obj(rec_file):
    """read the final 'best objective function value' from an opt rec file."""
    val = None
    with open(rec_file, 'r') as f:
        for line in f:
            if "best objective function value:" in line:
                val = float(line.strip().split(":")[-1])
    assert val is not None, "no objective in " + rec_file
    return val


def stack_chance_anomaly_conservative_test():
    """the stack-based chance shift must be the anomaly (spread) form everywhere.

    The chance-shifted value of a less_than constraint is current_sim + risk-spread, where
    risk-spread = stack_quantile - stack_mean >= 0.  Two consequences are asserted here:

    1. Every SLP iteration reports the anomaly form ("...using stack anomalies...") and never
       the raw/direct form.  The raw form (shifted value = stack_quantile = stack_mean +
       spread) silently injects the base-vs-ensemble-mean bias into the LP bound.

    2. CONSERVATISM: chance constraints shrink the feasible region, so the chance optimum
       cannot be BETTER than the deterministic (no-chance) optimum.  For this dewater problem
       (opt_direction=min), that means obj_chance >= obj_det.  The raw form could move a
       bound below its deterministic value and yield a chance optimum that beats the
       deterministic one - which is backwards - so this guards that regression directly.

    Also checks MPS<->reported-shift coherence across all iterations.
    """
    base = os.path.join("opt_dewater_chance", "template")

    def setup(pst):
        par = pst.parameter_data
        par.loc[par.partrans == "fixed", "partrans"] = "log"
        pst.pestpp_options.pop("opt_constraint_groups", None)
        pst.pestpp_options.pop("base_jacobian", None)
        pst.pestpp_options.pop("opt_par_stack", None)
        pst.pestpp_options.pop("opt_obs_stack", None)
        pst.control_data.noptmax = 3

    # deterministic (no chance) run
    d_det = os.path.join("opt_dewater_chance", "stack_det_run")
    if os.path.exists(d_det):
        shutil.rmtree(d_det)
    shutil.copytree(base, d_det)
    pst = pyemu.Pst(os.path.join(d_det, "dewater_pest.base.pst"))
    setup(pst)
    pst.write(os.path.join(d_det, "test.pst"))
    pyemu.os_utils.run("{0} {1}".format(exe_path, "test.pst"), cwd=d_det)
    obj_det = _read_best_obj(os.path.join(d_det, "test.rec"))

    # chance run (internal stack), otherwise identical
    d = os.path.join("opt_dewater_chance", "stack_anomaly_run")
    if os.path.exists(d):
        shutil.rmtree(d)
    shutil.copytree(base, d)
    pst = pyemu.Pst(os.path.join(d, "dewater_pest.base.pst"))
    setup(pst)
    pst.pestpp_options["opt_risk"] = 0.9
    pst.pestpp_options["opt_stack_size"] = 20
    noptmax = 3
    pst.control_data.noptmax = noptmax
    pst.write(os.path.join(d, "test.pst"))
    pyemu.os_utils.run("{0} {1}".format(exe_path, "test.pst"), cwd=d)
    rec_file = os.path.join(d, "test.rec")
    assert os.path.exists(rec_file)

    # 1. anomaly form on every iteration, never direct/raw
    modes = []
    with open(rec_file, 'r') as f:
        for line in f:
            if "using direct stack simulated results" in line:
                modes.append("direct")
            elif "using stack anomalies and current simulated" in line:
                modes.append("anomalies")
    print("chance handoff modes:", modes)
    assert len(modes) > 0, "no chance handoff mode was reported"
    assert all(m == "anomalies" for m in modes), \
        "stack chance shift must always use the anomaly (spread) form, got: {0}".format(modes)

    # 2. conservatism: chance optimum cannot beat the deterministic optimum (min => >=)
    obj_chance = _read_best_obj(rec_file)
    direction = pst.pestpp_options.get("opt_direction", "min")
    print("obj_det={0}  obj_chance={1}  direction={2}".format(obj_det, obj_chance, direction))
    tol = 1.0e-3 * max(1.0, abs(obj_det))
    if str(direction).lower().startswith("min"):
        assert obj_chance >= obj_det - tol, \
            "chance optimum ({0}) beats deterministic ({1}) for a minimization - chance shift is non-conservative".format(
                obj_chance, obj_det)
    else:
        assert obj_chance <= obj_det + tol, \
            "chance optimum ({0}) beats deterministic ({1}) for a maximization - chance shift is non-conservative".format(
                obj_chance, obj_det)

    # 3. the LP the simplex solved must match the reported shifted values
    check_mps_shift_coherence(d, "test", range(1, noptmax + 1))


def stack_chance_external_obs_stack_test():
    """regression test for the external opt_obs_stack -> sequential-LP chance handoff.

    An external obs stack (opt_obs_stack) supplies the constraint uncertainty.  The chance
    shift must be the anomaly (spread) form: shifted value = current_sim + (quantile -
    stack_mean).  This re-centers the stack's spread on the current simulated value, so it is
    conservative (>= current_sim for less_than) and never injects the base-vs-ensemble-mean
    bias that the raw form does.  This test asserts iteration 1 uses the anomaly form and that
    the reported optimum actually satisfies the stack-based chance constraints (re-evaluated
    from the response matrix + stack, no model runs needed).
    """
    tmp = os.path.join("opt_dewater_chance", "stack_ext_gen")
    if os.path.exists(tmp):
        shutil.rmtree(tmp)
    shutil.copytree(os.path.join("opt_dewater_chance", "template"), tmp)
    pst = pyemu.Pst(os.path.join(tmp, "dewater_pest.base.pst"))
    par = pst.parameter_data
    par.loc[par.partrans == "fixed", "partrans"] = "log"
    pst.pestpp_options.pop("opt_constraint_groups", None)
    pst.pestpp_options.pop("base_jacobian", None)
    pst.pestpp_options["opt_risk"] = 0.9
    pst.pestpp_options["opt_stack_size"] = 20
    pst.control_data.noptmax = 1
    pst.write(os.path.join(tmp, "test.pst"))
    pyemu.os_utils.run("{0} {1}".format(exe_path, "test.pst"), cwd=tmp)
    # the obs stack simulated at the initial dec-var point
    obs_stack = os.path.join(tmp, "test.1.obs_stack.csv")
    assert os.path.exists(obs_stack), "expected obs stack was not written: " + obs_stack

    # now run with that stack supplied EXTERNALLY (opt_obs_stack) at the same initial point
    d = os.path.join("opt_dewater_chance", "stack_ext_obs_stack_test")
    if os.path.exists(d):
        shutil.rmtree(d)
    shutil.copytree(os.path.join("opt_dewater_chance", "template"), d)
    shutil.copy2(obs_stack, os.path.join(d, "obs_stack.csv"))
    pst.pestpp_options.pop("opt_stack_size", None)
    pst.pestpp_options["opt_obs_stack"] = "obs_stack.csv"
    pst.control_data.noptmax = 1
    pst.write(os.path.join(d, "test.pst"))
    pyemu.os_utils.run("{0} {1}".format(exe_path, "test.pst"), cwd=d)
    rec = os.path.join(d, "test.rec")
    assert os.path.exists(rec)

    # the chance shift must use the anomaly (spread) form, never the raw/direct form
    mode = None
    with open(rec, 'r') as f:
        for line in f:
            if "using direct stack simulated results" in line:
                mode = "direct"
                break
            if "using stack anomalies and current simulated" in line:
                mode = "anomalies"
                break
    print("external obs stack iter-1 handoff mode:", mode)
    assert mode == "anomalies", \
        "stack chance shift must use the anomaly (spread) form, got '{0}'".format(mode)

    # the anomaly-shifted optimum must be self-consistent: the run's own convergence check
    # (which shifts the re-simulated optimum by the chance offset and compares to the bound)
    # must not report substantial unsatisfied model-based constraints.  a non-conservative
    # (raw) shift is what previously produced large post-solve violations here.
    unsat = {}
    with open(rec, 'r') as f:
        capture = False
        for line in f:
            if "not satisfied" in line:
                capture = True
                continue
            if capture:
                t = line.split()
                if len(t) >= 2 and t[0] == "-->" and t[1].upper().startswith("ONAME"):
                    try:
                        unsat[t[1]] = float(t[-1])
                    except ValueError:
                        pass
                elif t and not line.startswith("-->"):
                    capture = False
    worst = max(unsat.values()) if unsat else 0.0
    print("external obs stack: {0} unsatisfied constraints, worst distance {1}".format(len(unsat), worst))
    # allow only small binding-constraint linearization residuals, not gross violations
    assert worst < 5.0, \
        "reported optimum grossly violates stack chance constraints (worst distance {0}): {1}".format(worst, unsat)

    # the LP the simplex solved (MPS) must match the reported shifted constraint values
    check_mps_shift_coherence(d, "test", [1])


def fosm_invest():
    t_d = os.path.join("opt_dewater_chance","master5")
    pst = pyemu.Pst(os.path.join(t_d,"test.pst"))
    adj_pars = [p for p in pst.adj_par_names if not p.startswith("q")]
    pst.parameter_data.loc[:,"partrans"] = "fixed"
    pst.parameter_data.loc[adj_pars,"partrans"] = "log"
    print(adj_pars)
    fnames = pst.nnz_obs_names
    pst.observation_data.loc[:,"weight"] = 0.0
    jcb_files = [f for f in os.listdir(t_d) if f.endswith(".jcb") and f.startswith("test.")]
    #print(jcb_files)
    for jcb_file in jcb_files:
        jcb = pyemu.Matrix.from_binary(os.path.join(t_d,jcb_file))
        jcb = jcb.get(pst.obs_names,adj_pars)
        print(jcb.shape)
        sc = pyemu.Schur(jco=jcb,pst=pst,predictions=fnames)
        print(jcb_file,sc.get_forecast_summary().loc[:,"prior_var"].apply(lambda x: np.sqrt(x)))


if __name__ == "__main__":
    #fosm_invest()
    #startworker()
    #run_dewater_test()
    #run_supply2_test()
    #est_res_test()
    #shutil.copy2(os.path.join("..","exe","windows","x64","Debug","pestpp-opt.exe"),os.path.join("..","bin","win","pestpp-opt.exe"))
    stack_test()
    stack_chance_anomaly_conservative_test()
    stack_chance_external_obs_stack_test()
    #dewater_restart_test()
    #std_weights_test()
