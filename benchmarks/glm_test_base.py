# TODO: test variance and mean of draws, add chenoliver and test approx and full solution
import os
import shutil
import platform
import numpy as np
import pandas as pd
import platform
import matplotlib.pyplot as plt
import pyemu

bin_path = os.path.join("test_bin")
if "linux" in platform.platform().lower():
    bin_path = os.path.join(bin_path, "linux")
elif "darwin" in platform.platform().lower() or "mac" in platform.platform().lower():
    bin_path = os.path.join(bin_path, "mac")
else:
    bin_path = os.path.join(bin_path, "win")

bin_path = os.path.abspath("test_bin")
os.environ["PATH"] += os.pathsep + bin_path

# case of either appveyor, travis or local
if os.path.exists(os.path.join("pestpp", "bin")):
    bin_path = os.path.join("..", "..", "pestpp", "bin")
else:
    bin_path = os.path.join("..", "..", "..", "..", "pestpp", "bin")

if "windows" in platform.platform().lower():
    exe_path = os.path.join(bin_path, "win", "pestpp-glm.exe")
elif "darwin" in platform.platform().lower() or "mac" in platform.platform().lower():
    exe_path = os.path.join(bin_path, "mac", "pestpp-glm")
else:
    exe_path = os.path.join(bin_path, "linux", "pestpp-glm")

port = 4016


def tenpar_superpar_restart_test():
    model_d = "glm_10par_xsec"
    test_d = os.path.join(model_d, "master_basic_restart_test2")
    template_d = os.path.join(model_d, "template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)
    pst.parameter_data.loc[:, "partrans"] = "log"
    # pst.parameter_data.loc[pst.par_names[:2],"partrans"] = "log"
    pst.control_data.noptmax = -1
    pst_name = os.path.join(template_d, "pest_restart.pst")
    pst.write(pst_name)
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_restart.pst", num_workers=5,
                                 master_dir=test_d, verbose=True, worker_root=model_d,
                                 port=port)
    shutil.copy2(os.path.join(test_d, "pest_restart.jcb"), os.path.join(template_d, "restart.jcb"))
    shutil.copy2(os.path.join(test_d, "pest_restart.rei"), os.path.join(template_d, "restart.rei"))
    pst.pestpp_options["base_jacobian"] = "restart.jcb"
    pst.pestpp_options["hotstart_resfile"] = "restart.rei"
    pst.control_data.noptmax = -1
    pst.pestpp_options["glm_num_reals"] = 10
    pst_name = os.path.join(template_d, "pest_restart1.pst")
    pst.write(pst_name)
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_restart1.pst", num_workers=5,
                                 master_dir=test_d, verbose=True, worker_root=model_d,
                                 port=port)
    pst = pyemu.Pst(os.path.join(test_d, "pest_restart1.pst"))
    print(pst.phi)
    assert pst.phi == 50.0
    assert os.path.exists(os.path.join(test_d, "pest_restart1.post.obsen.csv"))

    pst.control_data.noptmax = 5
    pst.pestpp_options["glm_num_reals"] = 10
    pst_name = os.path.join(template_d, "pest_restart1.pst")
    pst.write(pst_name)
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_restart1.pst", num_workers=5,
                                 master_dir=test_d, verbose=True, worker_root=model_d,
                                 port=port)
    pst = pyemu.Pst(os.path.join(test_d, "pest_restart1.pst"))
    print(pst.phi)
    assert pst.phi < 1.0e-10, pst.phi
    assert os.path.exists(os.path.join(test_d, "pest_restart1.post.obsen.csv"))

    pst.control_data.noptmax = 5
    pst.pestpp_options["n_iter_base"] = -1
    pst.pestpp_options["n_iter_super"] = pst.control_data.noptmax
    pst.pestpp_options["glm_num_reals"] = 10

    pst_name = os.path.join(template_d, "pest_restart1.pst")
    pst.write(pst_name)
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_restart1.pst", num_workers=5,
                                 master_dir=test_d, verbose=True, worker_root=model_d,
                                 port=port)
    pst = pyemu.Pst(os.path.join(test_d, "pest_restart1.pst"))
    print(pst.phi)
    assert pst.phi < 0.05
    assert os.path.exists(os.path.join(test_d, "pest_restart1.post.obsen.csv"))
    par_unc1 = pd.read_csv(os.path.join(test_d, "pest_restart1.par.usum.csv"), index_col=0)

    cov = pyemu.Cov.from_parameter_data(pst)
    cov.to_binary(os.path.join(template_d, "prior.jcb"))
    pst.pestpp_options["parcov"] = "prior.jcb"
    pst.write(pst_name)
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_restart1.pst", num_workers=5,
                                 master_dir=test_d, verbose=True, worker_root=model_d,
                                 port=port)
    assert os.path.exists(os.path.join(test_d, "pest_restart1.post.obsen.csv"))
    par_unc2 = pd.read_csv(os.path.join(test_d, "pest_restart1.par.usum.csv"), index_col=0)

    diff = par_unc1 - par_unc2
    print(diff.apply(np.abs).sum())
    assert diff.apply(np.abs).sum().sum() == 0.0


def tenpar_base_test():
    """tenpar basic test"""

    model_d = "glm_10par_xsec"
    test_d = os.path.join(model_d, "master_basic_test2")
    template_d = os.path.join(model_d, "template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    # shutil.copytree(template_d, test_d)
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)
    pst.parameter_data.loc[:, "partrans"] = "log"
    # pst.parameter_data.loc[:,"parval1"] = pst.parameter_data.parubnd

    pst.pestpp_options = {}
    pst.pestpp_options["glm_num_reals"] = 10
    pst.control_data.noptmax = 3

    # pst.pestpp_options["ies_verbose_level"] = 3
    pst_name = os.path.join(template_d, "pest_basic.pst")
    pst.write(pst_name)
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_basic.pst", num_workers=10,
                                 master_dir=test_d, verbose=True, worker_root=model_d,
                                 port=port)
    pst_master = pyemu.Pst(os.path.join(test_d, "pest_basic.pst"))
    print(pst_master.phi)
    assert pst_master.phi < 1.0e-10, pst_master.phi
    assert os.path.exists(os.path.join(test_d, "pest_basic.post.obsen.csv"))

    pst.parameter_data.loc[pst.par_names[1:], "partrans"] = "fixed"
    pst.write(pst_name)
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_basic.pst", num_workers=5,
                                 master_dir=test_d, verbose=True, worker_root=model_d,
                                 port=port)
    pst_master = pyemu.Pst(os.path.join(test_d, "pest_basic.pst"))
    print(pst_master.phi)
    assert pst_master.phi < 1.0e-10, pst_master.phi
    assert os.path.exists(os.path.join(test_d, "pest_basic.post.obsen.csv"))

    pst.parameter_data.loc[pst.par_names[0], "parval1"] = pst.parameter_data.loc[pst.par_names[0], "parubnd"]
    pst.write(pst_name)
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_basic.pst", num_workers=5,
                                 master_dir=test_d, verbose=True, worker_root=model_d,
                                 port=port)
    pst_master = pyemu.Pst(os.path.join(test_d, "pest_basic.pst"))
    print(pst_master.phi)
    assert pst_master.phi < 1.0e-10, pst_master.phi
    assert os.path.exists(os.path.join(test_d, "pest_basic.post.obsen.csv"))


def freyberg_basic_restart_test():
    model_d = "glm_freyberg"
    test_d = os.path.join(model_d, "master_basic_test2")
    template_d = os.path.join(model_d, "template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    # shutil.copytree(template_d, test_d)
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)
    jco = pyemu.Matrix.from_binary(os.path.join(template_d, "pest_basic_restart.jcb"))
    jco = jco.get(pst.obs_names, pst.par_names)
    jco.to_binary(os.path.join(template_d, "temp.jcb"))
    pst.prior_information = pst.null_prior
    pst.control_data.pestmode = "estimation"
    pst.control_data.noptmax = 1
    pst.pestpp_options["glm_num_reals"] = 10
    pst.pestpp_options["base_jacobian"] = "temp.jcb"
    pst.pestpp_options["n_iter_base"] = -1
    pst.pestpp_options["n_iter_super"] = pst.control_data.noptmax
    pst.write(os.path.join(template_d, "pest_basic.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_basic.pst", num_workers=10,
                                 master_dir=test_d, verbose=True, worker_root=model_d,
                                 port=port)
    assert os.path.exists(os.path.join(test_d, "pest_basic.post.obsen.csv"))
    df_fosm_rei = pyemu.pst_utils.read_resfile(os.path.join(test_d, "pest_basic.fosm_reweight.rei"))
    pst = pyemu.Pst(os.path.join(test_d, "pest_basic.pst"))
    pst.adjust_weights_discrepancy()
    diff = (df_fosm_rei.loc[pst.nnz_obs_names, "weight"] - pst.observation_data.loc[pst.nnz_obs_names, "weight"]).apply(
        np.abs)
    print(df_fosm_rei.loc[pst.nnz_obs_names, "weight"])
    print(pst.observation_data.loc[pst.nnz_obs_names, "weight"])
    print(diff)
    assert diff.max() < 1.0e-5, diff.max()

    pst.control_data.noptmax = 4
    pst.pestpp_options["glm_num_reals"] = 10
    pst.pestpp_options["n_iter_base"] = -1
    pst.pestpp_options["n_iter_super"] = pst.control_data.noptmax
    pst.write(os.path.join(template_d, "pest_basic.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_basic.pst", num_workers=10,
                                 master_dir=test_d, verbose=True, worker_root=model_d,
                                 port=port)
    assert os.path.exists(os.path.join(test_d, "pest_basic.post.obsen.csv"))
    df_fosm_rei = pyemu.pst_utils.read_resfile(os.path.join(test_d, "pest_basic.fosm_reweight.rei"))
    pst = pyemu.Pst(os.path.join(test_d, "pest_basic.pst"))
    pst.adjust_weights_discrepancy()
    diff = (df_fosm_rei.loc[pst.nnz_obs_names, "weight"] - pst.observation_data.loc[pst.nnz_obs_names, "weight"]).apply(
        np.abs)
    print(df_fosm_rei.loc[pst.nnz_obs_names, "weight"])
    print(pst.observation_data.loc[pst.nnz_obs_names, "weight"])
    print(diff)
    assert diff.max() < 1.0e-5, diff.max()


def jac_diff_invest():
    model_d = "glm_10par_xsec"
    test_d = os.path.join(model_d, "master_basic_restart_test2")
    jco = pyemu.Jco.from_binary(os.path.join(test_d, "pest_restart1.jcb")).to_dataframe()
    pst = pyemu.Pst(os.path.join(test_d, "pest_restart1.pst"))
    for o in pst.nnz_obs_names:
        print(jco.loc[o, :])


def new_fmt_load_test():
    pst_d = os.path.join("ctl_read_test", "pst")
    test_d = os.path.join("ctl_read_test", "test2")
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    os.makedirs(test_d)
    pst_files = os.listdir(pst_d)

    def rec_scrape(pst_file):
        rec_file = pst_file.replace(".pst",".rec")
        assert os.path.exists(os.path.join(test_d,rec_file)),rec_file
        with open(os.path.join(test_d,rec_file),'r') as f:
            for line in f:

                if "number of adjustable parameters" in line.lower():
                    raw = line.split('=')
                    npar_adj = int(raw[1])
                elif "number of observations" in line.lower():
                    raw = line.split('=')
                    nobs = int(raw[1])
                elif "number of prior estimates" in line.lower():
                    raw = line.split('=')
                    npi = int(raw[1])
                elif "forecasts:" in line.lower():
                    raw = line.split(":")
                    forecasts = [r.strip().lower() for r in raw[1].split(',') if len(r.strip())]
                    #print(forecasts)

        return npar_adj,nobs,npi,forecasts

    for pst_file in pst_files:
        if ".pst" not in pst_file:
            continue
        print(pst_file)
        pst = pyemu.Pst(os.path.join(pst_d, pst_file))
        pst.pestpp_options["debug_parse_only"] = True
        pst.control_data.noptmax = 0
        if pst.nprior == 0:
            pst.control_data.pestmode = "estimation"

        pst.write(os.path.join(test_d, pst_file))
        pyemu.os_utils.run("{0} {1}".format(exe_path,pst_file), cwd=test_d)
        npar_adj,nobs,npi,forecasts = rec_scrape(pst_file)
        assert npar_adj == pst.npar_adj
        assert nobs == pst.nobs
        assert npi == pst.nprior
        if "predictions" in pst.pestpp_options:
            py_preds = set(pst.pestpp_options["predictions"].split(','))
            forecasts = set(forecasts)
            d = py_preds - forecasts
            print(pst_file,d)
            if len(d) > 0:
                raise Exception(pst_file + ": "+str(d))
        if  "forecasts" in pst.pestpp_options:
            py_preds = set(pst.pestpp_options["forecasts"].split(','))
            forecasts = set(forecasts)
            d = py_preds - forecasts
            print(pst_file,py_preds,forecasts)
            if len(d) > 0:
                raise Exception(pst_file + ": "+str(d))
        pst.write(os.path.join(test_d, "test.pst"), version=2)
        pyemu.os_utils.run("{0} {1}".format(exe_path,pst_file), cwd=test_d)
        npar_adj,nobs,npi,forecasts = rec_scrape(pst_file)
        assert npar_adj == pst.npar_adj
        assert nobs == pst.nobs
        assert npi == pst.nprior
        if "predictions" in pst.pestpp_options:
            py_preds = set(pst.pestpp_options["predictions"].split(','))
            forecasts = set(forecasts)
            d = py_preds - forecasts
            print(pst_file,d)
            if len(d) > 0:
                raise Exception(pst_file + ": "+str(d))
        if  "forecasts" in pst.pestpp_options:
            py_preds = set(pst.pestpp_options["forecasts"].split(','))
            forecasts = set(forecasts)
            d = py_preds - forecasts
            print(pst_file,py_preds,forecasts)
            if len(d) > 0:
                raise Exception(pst_file + ": "+str(d))

def tenpar_hotstart_test():
    model_d = "glm_10par_xsec"
    test_d = os.path.join(model_d, "master_hotstart_test")
    template_d = os.path.join(model_d, "template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    # shutil.copytree(template_d, test_d)
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)
    pst.control_data.pestmode = "estimation"
    pst.control_data.noptmax = -1
    pyemu.helpers.zero_order_tikhonov(pst)
  
    
    pst.write(os.path.join(template_d,"pest_temp.pst"))
    #pyemu.os_utils.run("{0} pest_temp.pst".format(exe_path),cwd=template_d)
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_temp.pst", num_workers=10,
                                 master_dir=test_d, verbose=True, worker_root=model_d,
                                 port=port)
    def super_failed(rec_file):
        with open(rec_file,'r') as f:
            for line in f:
                if "super parameter process" in line.lower():
                    return True
        return False

    pst.control_data.noptmax = 2
    shutil.copy2(os.path.join(test_d,"pest_temp.jcb"),os.path.join(template_d,"pest_temp.jcb"))
    shutil.copy2(os.path.join(test_d,"pest_temp.rei"),os.path.join(template_d,"pest_temp.rei"))
    pst.pestpp_options["base_jacobian"] = "pest_temp.jcb"
    pst.pestpp_options["n_iter_base"] = -1
    pst.pestpp_options["n_iter_super"] = 2
    #pst.pestpp_options["n_iter_super"] = pst.control_data.noptmax
    pst.pestpp_options["glm_num_reals"] = 5
    pst.write(os.path.join(template_d, "pest_hotstart.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_hotstart.pst", num_workers=10,
                                 master_dir=test_d, verbose=True, worker_root=model_d,
                                 port=port)
    assert os.path.exists(os.path.join(test_d,"pest_hotstart.post.obsen.csv"))
    if super_failed(os.path.join(test_d,"pest_hotstart.rec")):
        raise Exception("super failed")
    
    pst.control_data.noptmax = 2
    shutil.copy2(os.path.join(test_d,"pest_temp.jcb"),os.path.join(template_d,"pest_temp.jcb"))
    shutil.copy2(os.path.join(test_d,"pest_temp.rei"),os.path.join(template_d,"pest_temp.rei"))
    pst.pestpp_options["base_jacobian"] = "pest_temp.jcb"
    pst.pestpp_options["n_iter_base"] = 1
    #pst.pestpp_options["n_iter_super"] = pst.control_data.noptmax
    pst.pestpp_options["hotstart_resfile"] = "pest_temp.rei"
    pst.pestpp_options["glm_num_reals"] = 5

    pst.write(os.path.join(template_d, "pest_hotstart.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_hotstart.pst", num_workers=10,
                                 master_dir=test_d, verbose=True, worker_root=model_d,
                                 port=port)
    assert os.path.exists(os.path.join(test_d,"pest_hotstart.post.obsen.csv"))
    if super_failed(os.path.join(test_d,"pest_hotstart.rec")):
        raise Exception("super failed")

    pst.pestpp_options["n_iter_base"] = -1
    pst.pestpp_options["n_iter_super"] = 2
    pst.control_data.noptmax = 2
    pst.write(os.path.join(template_d, "pest_hotstart.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_hotstart.pst", num_workers=10,
                                 master_dir=test_d, verbose=True, worker_root=model_d,
                                 port=port)
    assert os.path.exists(os.path.join(test_d,"pest_hotstart.post.obsen.csv"))
    if super_failed(os.path.join(test_d,"pest_hotstart.rec")):
        raise Exception("super failed")
    
def tenpar_normalform_test():
    model_d = "glm_10par_xsec"
    test_d = os.path.join(model_d, "master_normal_test")
    template_d = os.path.join(model_d, "template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    # shutil.copytree(template_d, test_d)
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)
    with open(os.path.join(template_d,"dum.dat.tpl"),'w') as f:
        f.write("ptf ~\n~   dum1   ~\n")
    pst.add_parameters(os.path.join(template_d,"dum.dat.tpl"),pst_path=".")
    pst.control_data.pestmode = "estimation"
    pst.control_data.noptmax = 5
    pst.pestpp_options["glm_num_reals"] = 5
    pst.pestpp_options["n_iter_super"] = 3
    pst.pestpp_options["n_iter_base"] = 1
    pst.pestpp_options["glm_normal_form"] = "prior"
    pst.pestpp_options["max_n_super"] = 2
    pst.svd_data.maxsing = 2
    pst.write(os.path.join(template_d, "pest_prior.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_prior.pst", num_workers=10,
                                master_dir=test_d, verbose=True, worker_root=model_d,
                                port=port)

    pst.control_data.noptmax = -1
    par = pst.parameter_data
    par.loc["k_03","parval1"] = 0.9
    par.loc["k_03","parlbnd"] = 0.8
    par.loc["k_03","parubnd"] = 1.0
    par.loc["k_05","parubnd"] = 2.55
    pst.write(os.path.join(template_d,"pest_prior.pst"))
    #pyemu.os_utils.run("{0} pest_temp.pst".format(exe_path),cwd=template_d)
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_prior.pst", num_workers=10,
                                 master_dir=test_d, verbose=True, worker_root=model_d,
                                 port=port)
    shutil.copy2(os.path.join(test_d,"pest_prior.jcb"),os.path.join(template_d,"restart.jcb"))
    shutil.copy2(os.path.join(test_d,"pest_prior.rei"),os.path.join(template_d,"restart.rei"))
    pst.pestpp_options["base_jacobian"] = "restart.jcb"
    pst.pestpp_options["hotstart_resfile"] = "restart.rei"
    pst.pestpp_options["n_iter_base"] = -1
    pst.control_data.noptmax = 5
    pst.pestpp_options["n_iter_super"] = pst.control_data.noptmax
    par = pst.parameter_data
    par.loc["k_03","parval1"] = 0.9
    par.loc["k_03","parlbnd"] = 0.8
    par.loc["k_03","parubnd"] = 1.0
    par.loc["k_05","parubnd"] = 2.55
    pst.write(os.path.join(template_d,"pest_prior.pst"))
    #pyemu.os_utils.run("{0} pest_temp.pst".format(exe_path),cwd=template_d)
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_prior.pst", num_workers=10,
                                 master_dir=test_d, verbose=True, worker_root=model_d,
                                 port=port)
    pst.pestpp_options["glm_normal_form"] = "diag"
    pst.control_data.noptmax = 3
    pst.write(os.path.join(template_d, "pest_diag.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_diag.pst", num_workers=10,
                                 master_dir=test_d, verbose=True, worker_root=model_d,
                                 port=port)
    
    pst.pestpp_options["glm_normal_form"] = "ident"
    pst.write(os.path.join(template_d, "pest_diag.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_diag.pst", num_workers=10,
                                 master_dir=test_d, verbose=True, worker_root=model_d,
                                 port=port)


def freyberg_stress_test():
    model_d = "glm_freyberg"
    test_d = os.path.join(model_d, "master_stress")
    template_d = os.path.join(model_d, "template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    # shutil.copytree(template_d, test_d)
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)
    jco = pyemu.Matrix.from_binary(os.path.join(template_d, "pest_basic_restart.jcb"))
    jco = jco.get(pst.obs_names, pst.par_names)
    jco.to_binary(os.path.join(template_d, "temp.jcb"))
    pst.pestpp_options = {}
    par = pst.parameter_data
    par.loc[:,"parubnd"] *= 1.25
    par.loc[:,"parubnd"] *= 0.75
    pst.svd_data.maxsing = 10
    pst.prior_information = pst.null_prior
    pst.control_data.pestmode = "estimation"
    pst.control_data.noptmax = 5
    pst.pestpp_options["max_n_super"] = 10
    pst.pestpp_options["glm_num_reals"] = 10
    pst.pestpp_options["base_jacobian"] = "temp.jcb"
    pst.pestpp_options["n_iter_base"] = 1
    pst.pestpp_options["n_iter_super"] = pst.control_data.noptmax
    #pst.pestpp_options["glm_debug_der_fail"] = True
    pst.pestpp_options["glm_debug_lamb_fail"] = True
    pst.pestpp_options["glm_normal_form"] = "diag"
    pst.pestpp_options["parcov"] = "prior.jcb"
    pst.pestpp_options["glm_accept_mc_phi"] = True
    pst.write(os.path.join(template_d, "pest_stress.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_stress.pst", num_workers=10,
                                 master_dir=test_d, verbose=True, worker_root=model_d,
                                 port=port)
    pst = pyemu.Pst(os.path.join(test_d,"pest_stress.pst"))
    print(pst.phi)
    assert pst.phi < 35
    oe = pd.read_csv(os.path.join(test_d,"pest_stress.post.obsen.csv"),index_col=0)
    assert oe.dropna().shape == (int(pst.pestpp_options["glm_num_reals"]),pst.nobs),oe.dropna().shape
    

def tenpar_xsec_stress_test():
    model_d = "glm_10par_xsec"
    test_d = os.path.join(model_d, "master_stress")
    template_d = os.path.join(model_d, "template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    # shutil.copytree(template_d, test_d)
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)
    pst.pestpp_options = {}
    par = pst.parameter_data
    par.loc[:,"parubnd"] *= 1.25
    par.loc[:,"parubnd"] *= 0.75
    pst.svd_data.maxsing = 10
    pst.prior_information = pst.null_prior
    pst.control_data.pestmode = "estimation"
    pst.control_data.noptmax = 4
    pst.pestpp_options["glm_num_reals"] = 10
    pst.pestpp_options["n_iter_base"] = 1
    pst.pestpp_options["n_iter_super"] = pst.control_data.noptmax
    pst.pestpp_options["glm_debug_der_fail"] = False
    pst.pestpp_options["glm_debug_lamb_fail"] = True
    pst.pestpp_options["glm_normal_form"] = "prior"
    pst.pestpp_options["glm_accept_mc_phi"] = True
    pst.write(os.path.join(template_d, "pest_stress.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_stress.pst", num_workers=10,
                                 master_dir=test_d, verbose=True, worker_root=model_d,
                                 port=port)
    pst = pyemu.Pst(os.path.join(test_d,"pest_stress.pst"))
    print(pst.phi)
    assert pst.phi < .0002
    oe = pd.read_csv(os.path.join(test_d,"pest_stress.post.obsen.csv"),index_col=0)
    assert oe.dropna().shape == (int(pst.pestpp_options["glm_num_reals"]),pst.nobs),oe.dropna().shape

    

def threept_fail_test():
    model_d = "glm_10par_xsec"
    test_d = os.path.join(model_d, "master_3pt_fail")
    template_d = os.path.join(model_d, "template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)

    pst_name = os.path.join(test_d, "pest.pst")
    pst = pyemu.Pst(pst_name)
    pst.pestpp_options = {}
    pst.parameter_groups.loc[:,"forcen"] = "always_3"
    pst.control_data.noptmax = 1
    pst.write(os.path.join(test_d,"pest_fail.pst"))


def tenpar_xsec_high_phi_test():
    model_d = "glm_10par_xsec"
    test_d = os.path.join(model_d, "master_high_phi")
    template_d = os.path.join(model_d, "template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    # shutil.copytree(template_d, test_d)
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)
    pst.pestpp_options = {}
    par = pst.parameter_data
    par.loc[:,"parubnd"] *= 1.25
    par.loc[:,"parubnd"] *= 0.75
    pst.svd_data.maxsing = 10
    pst.prior_information = pst.null_prior
    pst.control_data.pestmode = "estimation"
    pst.control_data.noptmax = 2
    pst.pestpp_options["glm_num_reals"] = 10
    pst.pestpp_options["n_iter_base"] = 1
    pst.pestpp_options["n_iter_super"] = 2
    
    pst.pestpp_options["glm_debug_lamb_fail"] = True
    pst.pestpp_options["glm_debug_high_2nd_iter_phi"] = True
    pst.write(os.path.join(template_d, "pest_test.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_test.pst", num_workers=10,
                                 master_dir=test_d, verbose=True, worker_root=model_d,
                                 port=port)
    pst = pyemu.Pst(os.path.join(test_d,"pest_test.pst"))
    print(pst.phi)
    assert pst.phi < .5
    oe = pd.read_csv(os.path.join(test_d,"pest_test.post.obsen.csv"),index_col=0)
    assert oe.dropna().shape == (int(pst.pestpp_options["glm_num_reals"]),pst.nobs),oe.dropna().shape
    p2df = pyemu.pst_utils.read_parfile(os.path.join(test_d,"pest_test.2.par")).parval1
    pdf = pyemu.pst_utils.read_parfile(os.path.join(test_d,"pest_test.par")).parval1
    d = (p2df - pdf).apply(np.abs).sum()
    print(d)   
    assert d > 0.1 

def tenpar_xsec_stress_test2a():
    model_d = "glm_10par_xsec"
    test_d = os.path.join(model_d, "master_stress2")
    template_d = os.path.join(model_d, "template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    # shutil.copytree(template_d, test_d)
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)
    pst.pestpp_options = {}
    par = pst.parameter_data
    par.loc[:,"parubnd"] *= 1.25
    par.loc[:,"parubnd"] *= 0.75
    pst.svd_data.maxsing = 10
    #pst.prior_information = pst.null_prior
    pyemu.helpers.zero_order_tikhonov(pst)
    #phimaccept", "fracphim", "wfinit"]
    pst.reg_data.phimlim = 1.0
    pst.reg_data.phimaccept = 1.05
    pst.reg_data.fracphim = 0.5
    pst.reg_data.wfinit = 1.5


    pst.control_data.pestmode = "regularization"
    pst.control_data.noptmax = 10
    pst.pestpp_options["glm_num_reals"] = 10
    pst.pestpp_options["n_iter_base"] = 10
    #pst.pestpp_options["n_iter_super"] = pst.control_data.noptmax
    #pst.pestpp_options["glm_debug_der_fail"] = True
    pst.pestpp_options["glm_debug_lamb_fail"] = True
    pst.pestpp_options["glm_normal_form"] = "diag"
    pst.pestpp_options["glm_accept_mc_phi"] = True
    pst.write(os.path.join(template_d, "pest_stress.pst"))
    #return
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_stress.pst", num_workers=10,
                                 master_dir=test_d, verbose=True, worker_root=model_d,
                                 port=port)
    pst = pyemu.Pst(os.path.join(test_d,"pest_stress.pst"))
    print(pst.phi)
    assert np.abs((pst.phi - pst.reg_data.phimlim) / pst.phi) < 0.1,pst.phi
    oe = pd.read_csv(os.path.join(test_d,"pest_stress.post.obsen.csv"),index_col=0)
    assert oe.dropna().shape == (int(pst.pestpp_options["glm_num_reals"]),pst.nobs),oe.dropna().shape


def tenpar_xsec_stress_test_2():
    model_d = "glm_10par_xsec"
    test_d = os.path.join(model_d, "master_stress2")
    template_d = os.path.join(model_d, "template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    # shutil.copytree(template_d, test_d)
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)
    pst.pestpp_options = {}
    par = pst.parameter_data
    par.loc[:,"parubnd"] *= 1.25
    par.loc[:,"parubnd"] *= 0.75
    pst.svd_data.maxsing = 10
    pst.control_data.pestmode = "estimation"
    pst.prior_information = pst.null_prior
    pst.control_data.noptmax = 4
    pst.pestpp_options["glm_num_reals"] = 10
    pst.pestpp_options["n_iter_base"] = 10
    #pst.pestpp_options["n_iter_super"] = pst.control_data.noptmax
    pst.pestpp_options["glm_debug_der_fail"] = True
    pst.pestpp_options["glm_debug_lamb_fail"] = True
    pst.pestpp_options["glm_accept_mc_phi"] = True
    #pst.pestpp_options["glm_debug_high_2nd_iter_phi"] = True
    pst.pestpp_options["glm_iter_mc"] = True
    #pst.reg_data.phimlim = 10
    #pst.reg_data.fracphim = 0.5
    #pst.reg_data.phimaccept = 11
    #pst.reg_data.wfinit = 10
    pst.write(os.path.join(template_d, "pest_stress.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_stress.pst", num_workers=10,
                                 master_dir=test_d, verbose=True, worker_root=model_d,
                                 port=port)
    pst = pyemu.Pst(os.path.join(test_d,"pest_stress.pst"))
    print(pst.phi)
    assert os.path.exists(os.path.join(test_d,"pest_stress.4.par"))
    oe = pd.read_csv(os.path.join(test_d,"pest_stress.post.obsen.csv"),index_col=0)
    assert oe.dropna().shape == (int(pst.pestpp_options["glm_num_reals"]),pst.nobs),oe.dropna().shape


def tenpar_xsec_stress_test_3():
    model_d = "glm_10par_xsec"
    test_d = os.path.join(model_d, "master_stress3")
    template_d = os.path.join(model_d, "template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    # shutil.copytree(template_d, test_d)
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)
    pst.pestpp_options = {}
    par = pst.parameter_data
    par.loc[:,"parubnd"] *= 1.25
    par.loc[:,"parubnd"] *= 0.75
    pst.svd_data.maxsing = 10
    pst.control_data.pestmode = "regularization"
    pst.prior_information = pst.null_prior
    pyemu.helpers.zero_order_tikhonov(pst)
    pst.control_data.noptmax = 4
    pst.pestpp_options["glm_num_reals"] = 10
    pst.pestpp_options["n_iter_base"] = 10
    #pst.pestpp_options["n_iter_super"] = pst.control_data.noptmax
    pst.pestpp_options["glm_debug_der_fail"] = False
    pst.pestpp_options["glm_debug_lamb_fail"] = False
    pst.pestpp_options["glm_accept_mc_phi"] = True
    #pst.pestpp_options["glm_debug_high_2nd_iter_phi"] = True
    pst.pestpp_options["glm_iter_mc"] = True
    #pst.reg_data.phimlim = 10
    #pst.reg_data.fracphim = 0.5
    #pst.reg_data.phimaccept = 11
    #pst.reg_data.wfinit = 10
    pst.write(os.path.join(template_d, "pest_stress.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_stress.pst", num_workers=10,
                                 master_dir=test_d, verbose=True, worker_root=model_d,
                                 port=port)
    pst = pyemu.Pst(os.path.join(test_d,"pest_stress.pst"))
    print(pst.phi)
    assert os.path.exists(os.path.join(test_d,"pest_stress.4.par"))
    oe = pd.read_csv(os.path.join(test_d,"pest_stress.post.obsen.csv"),index_col=0)
    assert oe.dropna().shape == (int(pst.pestpp_options["glm_num_reals"]),pst.nobs),oe.dropna().shape

 
    shutil.copy2(os.path.join(test_d,"pest_stress.jcb"),os.path.join(template_d,"restart.jcb"))
    pst.pestpp_options["base_jacobian"] = "restart.jcb"
    pst.control_data.noptmax = 1
    pst.write(os.path.join(template_d, "pest_stress.pst"))
    test_d += "a"
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_stress.pst", num_workers=10,
                                 master_dir=test_d, verbose=True, worker_root=model_d,
                                 port=port)


def tenpar_xsec_stress_test_4():
    model_d = "glm_10par_xsec"
    test_d = os.path.join(model_d, "master_stress4")
    template_d = os.path.join(model_d, "template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    # shutil.copytree(template_d, test_d)
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)
    pst.pestpp_options = {}
    par = pst.parameter_data
    par.loc[:,"parubnd"] *= 1.25
    par.loc[:,"parubnd"] *= 0.75
    pst.svd_data.maxsing = 10
    pst.control_data.pestmode = "regularization"
    pst.prior_information = pst.null_prior
    pyemu.helpers.zero_order_tikhonov(pst)
    pst.control_data.noptmax = -1
    pst.pestpp_options["glm_num_reals"] = 10
    pst.pestpp_options["n_iter_base"] = 10
    #pst.pestpp_options["n_iter_super"] = pst.control_data.noptmax
    pst.pestpp_options["glm_debug_der_fail"] = False
    pst.pestpp_options["glm_debug_lamb_fail"] = False
    pst.pestpp_options["glm_accept_mc_phi"] = False
    #pst.pestpp_options["glm_debug_high_2nd_iter_phi"] = True
    pst.pestpp_options["glm_iter_mc"] = True
    #pst.reg_data.phimlim = 10
    #pst.reg_data.fracphim = 0.5
    #pst.reg_data.phimaccept = 11
    #pst.reg_data.wfinit = 10
    pst.write(os.path.join(template_d, "pest_stress.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_stress.pst", num_workers=10,
                                 master_dir=test_d, verbose=True, worker_root=model_d,
                                 port=port)
    
    
    shutil.copy2(os.path.join(test_d,"pest_stress.jcb"),os.path.join(template_d,"restart.jcb"))
    pst.pestpp_options["base_jacobian"] = "restart.jcb"
    pst.pestpp_options["glm_normal_form"] = "prior"
    pyemu.Cov.from_parameter_data(pst).to_ascii(os.path.join(template_d,"prior.cov"))
    pst.pestpp_options["parcov"] = "prior.cov"
    pst.control_data.pestmode = "estimation"
    pst.control_data.noptmax = 1
    pst.pestpp_options["n_iter_base"] = -1
    pst.pestpp_options["n_iter_super"] = 2
    pst.write(os.path.join(template_d, "pest_stress.pst"))
    test_d += "a"
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_stress.pst", num_workers=10,
                                 master_dir=test_d, verbose=True, worker_root=model_d,
                                 port=port)

def invest():
    model_d = "glm_10par_xsec"
    test_d = os.path.join(model_d, "master_stress2")
    jco = pyemu.Jco.from_binary(os.path.join(test_d,"pest_stress.jcb"))
    print(jco.to_dataframe())


def tenpar_xsec_stress_test_5():
    model_d = "glm_10par_xsec"
    test_d = os.path.join(model_d, "master_stress5")
    template_d = os.path.join(model_d, "template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    # shutil.copytree(template_d, test_d)
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)
    pst.pestpp_options = {}
    par = pst.parameter_data
    #par["partrans"] = "none"
    #par.loc[:,"parubnd"] *= 1.25
    #par.loc[:,"parubnd"] *= 0.75
    
    obs = pst.observation_data
    obs["weight"] = 100
    #np.random.seed(1123441)
    #obs["obsval"] += np.random.normal(0,10,obs.shape[0])
    obs["obsval"] += 100
    pst.svd_data.maxsing = 10

    pst.prior_information = pst.null_prior
    pst.control_data.pestmode = "estimation"
    pst.control_data.noptmax = 4 #hard coded to result check below...
    pst.pestpp_options["glm_num_reals"] = 10
    #pst.pestpp_options["n_iter_base"] = 1
    #pst.pestpp_options["n_iter_super"] = pst.control_data.noptmax
    pst.pestpp_options["glm_debug_der_fail"] = False
    pst.pestpp_options["glm_debug_lamb_fail"] = True
    pst.pestpp_options["glm_normal_form"] = "prior"
    pst.pestpp_options["glm_accept_mc_phi"] = True
    pst.write(os.path.join(template_d, "pest_stress.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_stress.pst", num_workers=10,
                                 master_dir=test_d, verbose=True, worker_root=model_d,
                                 port=port)
    pst = pyemu.Pst(os.path.join(test_d,"pest_stress.pst"))
    print(pst.phi)
    assert os.path.exists(os.path.join(test_d,"pest_stress.4.par"))

    df = pd.read_csv(os.path.join(test_d,"pest_stress.upg.csv"),index_col=0)
    violations = []
    for name,lb,ub in zip(par.parnme,par.parlbnd,par.parubnd):
        vals = df[name]
        lb_out = vals[vals<lb]
        ub_out = vals[vals>ub]
        if lb_out.shape[0]:
            print(name,lb_out)
            violations.append(name)
        if ub_out.shape[0]:
            print(name,ub_out)
            violations.append(name)
    if len(violations) > 0:
        raise Exception("the following pars are out of bounds in upg.csv:"+",".join(violations))

    df = pd.read_csv(os.path.join(test_d,"pest_stress.ipar"),index_col=0)
    print(df.shape)
    assert df.shape[0] == 5
    violations = []
    for name,lb,ub in zip(par.parnme,par.parlbnd,par.parubnd):
        vals = df[name]
        lb_out = vals[vals<lb]
        ub_out = vals[vals>ub]
        if lb_out.shape[0]:
            print(name,lb_out)
            violations.append(name)
        if ub_out.shape[0]:
            print(name,ub_out)
            violations.append(name)
    if len(violations) > 0:
        raise Exception("the following pars are out of bounds in ipar:"+",".join(violations))


def tenpar_xsec_stress_test_5super():
    model_d = "glm_10par_xsec"
    test_d = os.path.join(model_d, "master_stress5")
    template_d = os.path.join(model_d, "template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    # shutil.copytree(template_d, test_d)
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)
    pst.pestpp_options = {}
    par = pst.parameter_data
    #par["partrans"] = "none"
    #par.loc[:,"parubnd"] *= 1.25
    #par.loc[:,"parubnd"] *= 0.75
    
    obs = pst.observation_data
    obs["weight"] = 100
    #np.random.seed(1123441)
    #obs["obsval"] += np.random.normal(0,10,obs.shape[0])
    obs["obsval"] += 100
    pst.svd_data.maxsing = 10

    pst.prior_information = pst.null_prior
    pst.control_data.pestmode = "estimation"
    pst.control_data.noptmax = 4
    pst.pestpp_options["glm_num_reals"] = 10
    pst.pestpp_options["n_iter_base"] = 1
    pst.pestpp_options["n_iter_super"] = pst.control_data.noptmax
    pst.pestpp_options["glm_debug_der_fail"] = False
    pst.pestpp_options["glm_debug_lamb_fail"] = True
    pst.pestpp_options["glm_normal_form"] = "prior"
    pst.pestpp_options["glm_accept_mc_phi"] = True
    pst.write(os.path.join(template_d, "pest_stress.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_stress.pst", num_workers=10,
                                 master_dir=test_d, verbose=True, worker_root=model_d,
                                 port=port)
    pst = pyemu.Pst(os.path.join(test_d,"pest_stress.pst"))
    print(pst.phi)
    assert os.path.exists(os.path.join(test_d,"pest_stress.4.par"))

    df = pd.read_csv(os.path.join(test_d,"pest_stress.upg.csv"),index_col=0)
    violations = []
    for name,lb,ub in zip(par.parnme,par.parlbnd,par.parubnd):
        vals = df[name]
        lb_out = vals[vals<lb]
        ub_out = vals[vals>ub]
        if lb_out.shape[0]:
            print(name,lb_out)
            violations.append(name)
        if ub_out.shape[0]:
            print(name,ub_out)
            violations.append(name)
    if len(violations) > 0:
        raise Exception("the following pars are out of bounds in upg.csv:"+",".join(violations))

    df = pd.read_csv(os.path.join(test_d,"pest_stress.ipar"),index_col=0)
    print(df.shape)
    assert df.shape[0] == 5
    violations = []
    for name,lb,ub in zip(par.parnme,par.parlbnd,par.parubnd):
        vals = df[name]
        lb_out = vals[vals<lb]
        ub_out = vals[vals>ub]
        if lb_out.shape[0]:
            print(name,lb_out)
            violations.append(name)
        if ub_out.shape[0]:
            print(name,ub_out)
            violations.append(name)
    if len(violations) > 0:
        raise Exception("the following pars are out of bounds in ipar:"+",".join(violations))

if __name__ == "__main__":
    #tenpar_base_test()
    #tenpar_superpar_restart_test()
    #freyberg_basic_restart_test()
    # jac_diff_invest()
    new_fmt_load_test()
    #tenpar_hotstart_test()
    #tenpar_normalform_test()
    #freyberg_stress_test()
    #shutil.copy2(os.path.join("..","exe","windows","x64","Debug","pestpp-glm.exe"),os.path.join("..","bin","win","pestpp-glm.exe"))
    #tenpar_xsec_stress_test_2()
    #tenpar_xsec_stress_test_3()
    #tenpar_xsec_stress_test_5super()
    
    #tenpar_xsec_stress_test_5()

    

        
    #invest()
    #tenpar_xsec_high_phi_test()
    
    #new_fmt_load_test()
    #threept_fail_test()
