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
exe_path = os.path.join(bin_path, "pestpp-ies" + exe)


noptmax = 4
num_reals = 20
port = 4021

def basic_test(model_d="ies_10par_xsec"):
    pyemu.Ensemble.reseed()
    base_d = os.path.join(model_d, "template")
    new_d = os.path.join(model_d, "test_template")
    if os.path.exists(new_d):
        shutil.rmtree(new_d)
    shutil.copytree(base_d, new_d)
    print(platform.platform().lower())
    local=True
    if "linux" in platform.platform().lower() and "10par" in model_d:
        #print("travis_prep")
        #prep_for_travis(model_d)
        local=False
    pst = pyemu.Pst(os.path.join(new_d, "pest.pst"))
    print(pst.model_command)
    
    # set first par as fixed
    #pst.parameter_data.loc[pst.par_names[0], "partrans"] = "fixed"

    pst.observation_data.loc[pst.nnz_obs_names,"weight"] = 1.0

    # set noptmax
    pst.control_data.noptmax = noptmax

    # wipe all pestpp options
    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = num_reals
    pst.pestpp_options["lambda_scale_fac"] = [0.5,0.75,1.0]
    pst.pestpp_options["ies_lambda_mults"] = 1.0
    # write a generic 2D cov
    if os.path.exists(os.path.join(new_d,"prior.jcb")):
        cov = pyemu.Cov.from_binary(os.path.join(new_d,"prior.jcb"))
        #cov.to_ascii(os.path.join(new_d,"prior.cov"))
    elif os.path.exists(os.path.join(new_d, "prior.cov")):
        cov = pyemu.Cov.from_ascii(os.path.join(new_d, "prior.cov"))
    else:
        cov = pyemu.Cov.from_parameter_data(pst)
        cov = pyemu.Cov(cov.as_2d, names=cov.row_names)
        #cov.to_ascii(os.path.join(new_d, "prior.cov"))
        cov.to_binary(os.path.join(new_d, "prior.jcb"))

    # draw some ensembles
    idx = [i for i in range(num_reals)]
    idx[-1] = "base"
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst, cov=cov, num_reals=num_reals)
    pe.index = idx
    pe.to_csv(os.path.join(new_d, "par.csv"))
    pe.to_binary(os.path.join(new_d, "par.jcb"))
    pe.to_csv(os.path.join(new_d, "sweep_in.csv"))
    pe.loc[:, pst.adj_par_names].to_csv(os.path.join(new_d, "par_some.csv"))
    pe.iloc[:-3, :].to_csv(os.path.join(new_d, "restart_failed_par.csv"))
    oe = pyemu.ObservationEnsemble.from_gaussian_draw(pst, num_reals=num_reals)
    oe.index = idx
    oe.to_csv(os.path.join(new_d, "obs.csv"))
    oe.iloc[:-3, :].to_csv(os.path.join(new_d, "restart_failed_base_obs.csv"))
    oe.to_binary(os.path.join(new_d, "obs.jcb"))

    pst.control_data.noptmax = 0
    pst.write(os.path.join(new_d, "pest.pst"))
    pyemu.os_utils.run("{0} pest.pst".format(exe_path),cwd=new_d)
    df = pd.read_csv(os.path.join(new_d,"pest.phi.group.csv"))
    assert df.loc[0,"head"] == 0.5,df
    #return
    pst.control_data.noptmax = noptmax
    pst.write(os.path.join(new_d, "pest.pst"))
    

    
    m_d = os.path.join(model_d,"master_pestpp_sen")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pyemu.os_utils.start_workers(new_d, exe_path.replace("-ies","-sen"), "pest.pst", 5, master_dir=m_d,
                           worker_root=model_d,local=local,port=port)
    df = pd.read_csv(os.path.join(m_d, "pest.mio"),index_col=0)

    # run sweep
    m_d = os.path.join(model_d,"master_sweep1")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pyemu.os_utils.start_workers(new_d, exe_path.replace("-ies","-swp"), "pest.pst", 5, master_dir=m_d,
                           worker_root=model_d,local=local,port=port)
    df = pd.read_csv(os.path.join(m_d, "sweep_out.csv"),index_col=0)
    
    m_d = os.path.join(model_d,"master_pestpp-glm")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pyemu.os_utils.start_workers(new_d, exe_path.replace("-ies","-glm"), "pest.pst", 10, master_dir=m_d,
                           worker_root=model_d,local=local,port=port)

    m_d = os.path.join(model_d,"master_pestpp-ies")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pyemu.os_utils.start_workers(new_d, exe_path, "pest.pst", 10, master_dir=m_d,
                           worker_root=model_d,local=local,port=port)



def glm_save_binary_test():
    model_d = "ies_10par_xsec"
    local = True
    if "linux" in platform.platform().lower() and "10par" in model_d:
        # print("travis_prep")
        # prep_for_travis(model_d)
        local = False

    t_d = os.path.join(model_d, "template")
    m_d = os.path.join(model_d, "master_save_binary")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pst = pyemu.Pst(os.path.join(t_d, "pest.pst"))
    pst.pestpp_options = {"glm_num_reals":30,"save_binary":True}
    pst.control_data.noptmax = 1
    pst.write(os.path.join(t_d, "pest_save_binary.pst"))
    pyemu.os_utils.start_workers(t_d, exe_path.replace("-ies", "-glm"), "pest_save_binary.pst", 10, master_dir=m_d,
                                worker_root=model_d, local=local, port=port)

    pe = pyemu.ParameterEnsemble.from_binary(pst=pst,filename=os.path.join(m_d,"pest_save_binary.post.paren.jcb"))
    pe = pyemu.ObservationEnsemble.from_binary(pst=pst,filename=os.path.join(m_d, "pest_save_binary.post.obsen.jcb"))

def sweep_forgive_test():
    model_d = "ies_10par_xsec"
    local=True
    if "linux" in platform.platform().lower() and "10par" in model_d:
        #print("travis_prep")
        #prep_for_travis(model_d)
        local=False
    
    t_d = os.path.join(model_d,"template")
    m_d = os.path.join(model_d,"master_sweep_forgive")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pst = pyemu.Pst(os.path.join(t_d,"pest.pst"))
    pe = pyemu.ParameterEnsemble.from_uniform_draw(pst,num_reals=50)#.loc[:,pst.par_names[:2]]
    pe.loc[:,pst.par_names[2:]] = pst.parameter_data.loc[pst.par_names[2:],"parval1"].values
    pe.to_csv(os.path.join(t_d,"sweep_in.csv"))
    pst.pestpp_options["sweep_forgive"] = True
    pst.write(os.path.join(t_d,"pest_forgive.pst"))
    pyemu.os_utils.start_workers(t_d, exe_path.replace("-ies","-swp"), "pest_forgive.pst", 10, master_dir=m_d,
                           worker_root=model_d,local=local,port=port)
    df1 = pd.read_csv(os.path.join(m_d, "sweep_out.csv"),index_col=0)

    pe = pe.loc[:,pst.par_names[:2]]
    pe.to_csv(os.path.join(t_d,"sweep_in.csv"))
    pst.pestpp_options["sweep_forgive"] = True
    pst.write(os.path.join(t_d,"pest_forgive.pst"))
    pyemu.os_utils.start_workers(t_d, exe_path.replace("-ies","-swp"), "pest_forgive.pst", 10, master_dir=m_d,
                           worker_root=model_d,local=local,port=port)
    df2 = pd.read_csv(os.path.join(m_d, "sweep_out.csv"),index_col=0)
    diff = df1 - df2
    print(diff.max())
    assert diff.max().max() == 0.0


def inv_regul_test():
    model_d = "ies_10par_xsec"
    local=True
    if "linux" in platform.platform().lower() and "10par" in model_d:
        #print("travis_prep")
        #prep_for_travis(model_d)
        local=False
    
    t_d = os.path.join(model_d,"template")
    m_d = os.path.join(model_d,"master_inv_regul")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pst = pyemu.Pst(os.path.join(t_d,"pest.pst"))
    #pyemu.helpers.zero_order_tikhonov(pst)
    #pst.control_data.pestmode = "regularization"
    pst.reg_data.phimlim = 2
    pst.reg_data.phimaccept = 2.2
    pst.control_data.noptmax = 10
    pst.write(os.path.join(t_d,"pest_regul.pst"))
    pyemu.os_utils.start_workers(t_d, exe_path.replace("-ies","-glm"), "pest_regul.pst", 10, master_dir=m_d,
                           worker_root=model_d,local=local,port=port)
    

def tie_by_group_test():
    model_d = "ies_10par_xsec"
    local=True
    if "linux" in platform.platform().lower() and "10par" in model_d:
        #print("travis_prep")
        #prep_for_travis(model_d)
        local=False
    
    t_d = os.path.join(model_d,"template")
    m_d = os.path.join(model_d,"master_tie_by_group")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pst = pyemu.Pst(os.path.join(t_d,"pest.pst")) 
    par = pst.parameter_data
    tied_names = pst.adj_par_names[:3]
    par.loc[tied_names[1:3],"partrans"] = "tied"
    par.loc[tied_names[1:3],"partied"] = tied_names[0]
    par.loc[tied_names[1:3],"parval1"] = 1.0
    par.loc[tied_names[1:3],"parubnd"] = par.loc[tied_names[1:3],"parval1"] * 1.0001
    par.loc[tied_names[1:3],"parlbnd"] = par.loc[tied_names[1:3],"parval1"] * 0.9999
    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options["ies_lambda_mults"] = 1.0
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["tie_by_group"] = True
    pst.pestpp_options["ies_include_base"] = False
    pst.pestpp_options["enforce_tied_bounds"] = True
    pst.control_data.noptmax = 1


    pst.write(os.path.join(t_d,"pest_tied.pst"))
    pyemu.os_utils.start_workers(t_d, exe_path.replace("-ies","-sen"), "pest_tied.pst", 5, master_dir=m_d,
                           worker_root=model_d,local=local,port=port)
    df = pd.read_csv(os.path.join(m_d,"pest_tied.sen.par.csv"),index_col=0)
    df.columns = df.columns.str.lower()
    print(df.loc[:,tied_names[1:3]])
    print(df.loc[:,tied_names[1:3]].std(axis=1))
    print(df.loc[:,tied_names[1:3]].std(axis=1).apply(np.abs).max())
    assert df.loc[:,tied_names[1:3]].std(axis=1).apply(np.abs).max() < 1.0e-8
    for real in df.index:
        too_low = df.loc[real,df.loc[real,par.parnme] < par.parlbnd]
        assert too_low.shape[0] == 0, "sen,{0},{1}".format(real,too_low)
        too_high = df.loc[real, df.loc[real, par.parnme] > par.parubnd]
        assert too_high.shape[0] == 0, "sen,{0},{1}".format(real,too_high)
    
    #pst.write(os.path.join(t_d,"pest_tied.pst"))
    pyemu.os_utils.start_workers(t_d, exe_path.replace("-ies","-glm"), "pest_tied.pst", 5, master_dir=m_d,
                           worker_root=model_d,local=local,port=port)
    jco = pyemu.Jco.from_binary(os.path.join(m_d,"pest_tied.jcb"))
    assert jco.shape[1] == 2,jco.shape
    par_df = pyemu.pst_utils.read_parfile(os.path.join(m_d,"pest_tied.par"))
    print(par_df)
    too_low = par.loc[par_df.parval1 < par.parlbnd,"parnme"]
    assert too_low.shape[0] == 0,too_low
    too_high = par.loc[par_df.parval1 > par.parubnd, "parnme"]
    assert too_high.shape[0] == 0, too_high

    
    pst.control_data.noptmax = 1
    pst.write(os.path.join(t_d, "pest_tied.pst"))

    pyemu.os_utils.start_workers(t_d, exe_path, "pest_tied.pst", 10, master_dir=m_d,
                           worker_root=model_d,local=local,port=port)
    
    for nopt in range(pst.control_data.noptmax+1):
        df = pd.read_csv(os.path.join(m_d,"pest_tied.{0}.par.csv".format(nopt)),index_col=0)
        df.columns = df.columns.str.lower()
        print(df.loc[:,tied_names[1:3]])
        print(df.loc[:,tied_names[1:3]].std(axis=1))
        print(df.loc[:,tied_names[1:3]].std(axis=1).apply(np.abs).max())
        assert df.loc[:,tied_names[1:3]].std(axis=1).apply(np.abs).max() < 1.0e-8
        for real in df.index:
            too_low = df.loc[real,df.loc[real,par.parnme] < par.parlbnd]
            assert too_low.shape[0] == 0, "ies,{0},{1},{2}".format(nopt,real,too_low)
            too_high = df.loc[real, df.loc[real, par.parnme] > par.parubnd]
            assert too_high.shape[0] == 0, "ies,{0},{1},{2}".format(nopt,real,too_high)


    par.loc[tied_names[1:3],"parval1"] = par.loc[tied_names[0],"parval1"]
    print(par.parval1)
    par.loc[tied_names[1:3], "parubnd"] = par.loc[tied_names[1:3], "parval1"] * 1.5
    par.loc[tied_names[1:3], "parlbnd"] = par.loc[tied_names[1:3], "parval1"] * 0.5
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options["ies_lambda_mults"] = 1.0
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["tie_by_group"] = True


    pst.control_data.noptmax = 3
    pst.write(os.path.join(t_d, "pest_tied.pst"))

    pyemu.os_utils.start_workers(t_d, exe_path, "pest_tied.pst", 10, master_dir=m_d,
                                worker_root=model_d, local=local, port=port)
    df = pd.read_csv(os.path.join(m_d, "pest_tied.{0}.par.csv".format(pst.control_data.noptmax)), index_col=0)
    df.columns = df.columns.str.lower()
    print(df.loc[:, tied_names].std(axis=1).apply(np.abs).max())
    assert df.loc[:, tied_names].std(axis=1).apply(np.abs).max() < 1.0e-8
    for real in df.index:
        too_low = df.loc[real,df.loc[real,par.parnme] < par.parlbnd]
        assert too_low.shape[0] == 0, "ies,{0},{1}".format(real,too_low)
        too_high = df.loc[real, df.loc[real, par.parnme] > par.parubnd]
        assert too_high.shape[0] == 0, "ies,{0},{1}".format(real,too_high)
    
    df.to_csv(os.path.join(t_d,"sweep_in.csv"))
    pyemu.os_utils.start_workers(t_d, exe_path.replace("-ies","-swp"), "pest_tied.pst", 5, master_dir=m_d,
                           worker_root=model_d,local=local,port=port)
    pst.control_data.noptmax = 3
    pst.pestpp_options["enforce_tied_bounds"] = False
    pst.write(os.path.join(t_d, "pest_tied.pst"))
    pyemu.os_utils.start_workers(t_d, exe_path.replace("-ies","-glm"), "pest_tied.pst", 5, master_dir=m_d,
                           worker_root=model_d,local=local,port=port)
    jco = pyemu.Jco.from_binary(os.path.join(m_d,"pest_tied.jcb"))
    assert jco.shape[1] == 2,jco.shape



def unc_file_test():
    model_d = "ies_10par_xsec"
    local=True
    if "linux" in platform.platform().lower() and "10par" in model_d:
        #print("travis_prep")
        #prep_for_travis(model_d)
        local=False
    
    t_d = os.path.join(model_d,"template")
    m_d = os.path.join(model_d,"master_uncfile")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    shutil.copytree(t_d,m_d)
    pst = pyemu.Pst(os.path.join(m_d,"pest.pst"))
    cov = pyemu.Cov.from_parameter_data(pst)
    cov.to_uncfile(os.path.join(m_d,"pest.unc"),covmat_file=os.path.join(m_d,"cov.mat"),var_mult=2.0,include_path=False)
    pst.pestpp_options = {}
    pst.pestpp_options["parcov"] = "pest.unc"
    pst.pestpp_options["ies_num_reals"] = 10000
    pst.pestpp_options["ies_verbose_level"] = 3
    pst.control_data.noptmax = -2
    pst.write(os.path.join(m_d,"pest_unc.pst"))
    pyemu.os_utils.run("{0} {1}".format(exe_path,"pest_unc.pst"),cwd=m_d)

    pe_1 = pd.read_csv(os.path.join(m_d,"pest_unc.0.par.csv"),index_col=0).apply(np.log10)

    cov = pyemu.Cov.from_parameter_data(pst)
    cov *= 2.0
    cov.to_uncfile(os.path.join(m_d,"pest.unc"),covmat_file=os.path.join(m_d,"cov.mat"),var_mult=1.0,include_path=False)
    pst.pestpp_options = {}
    pst.pestpp_options["parcov"] = "cov.mat"
    pst.pestpp_options["ies_num_reals"] = 10000
    pst.pestpp_options["ies_verbose_level"] = 3 
    pst.control_data.noptmax = -2
    pst.write(os.path.join(m_d,"pest_unc.pst"))
    pyemu.os_utils.run("{0} {1}".format(exe_path,"pest_unc.pst"),cwd=m_d)
    pe_2 = pd.read_csv(os.path.join(m_d,"pest_unc.0.par.csv"),index_col=0).apply(np.log10)
    diff = pe_1 - pe_2
    print(pe_1.std(ddof=0)**2)
    print(pe_2.std(ddof=0)**2)
    print(diff.sum())
    assert diff.sum().max() < 1.0e-10

def parchglim_test():
    model_d = "ies_10par_xsec"
    local=True
    if "linux" in platform.platform().lower() and "10par" in model_d:
        #print("travis_prep")
        #prep_for_travis(model_d)
        local=False
    
    t_d = os.path.join(model_d,"template")
    m_d = os.path.join(model_d,"master_parchglim")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    shutil.copytree(t_d,m_d)
    pst = pyemu.Pst(os.path.join(m_d,"pest.pst"))
    fpm = 1.05
    pst.control_data.facparmax = fpm
    par = pst.parameter_data
    par.loc[pst.par_names[1:],"partrans"] = "fixed"
    par.loc[pst.par_names[0],"partrans"] = "log"
    par.loc[pst.par_names[0],"parchglim"] = "factor"
    par.loc[pst.par_names[0],"parval1"] = 1.0
    
    pst.control_data.noptmax = 1
    pst.pestpp_options["lambdas"] = 1.0
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.write(os.path.join(m_d,"pest_parchglim.pst"))
    pyemu.os_utils.run("{0} pest_parchglim.pst".format(exe_path.replace("-ies","-glm")),cwd=m_d)
    p_df = pyemu.pst_utils.read_parfile(os.path.join(m_d,"pest_parchglim.par"))
    assert p_df.loc["stage","parval1"] == fpm

    rpm = 0.1
    par.loc[pst.par_names[0],"parchglim"] = "relative"
    pst.control_data.relparmax = rpm
    pst.write(os.path.join(m_d,"pest_parchglim.pst"))
    pyemu.os_utils.run("{0} pest_parchglim.pst".format(exe_path.replace("-ies","-glm")),cwd=m_d)
    p_df = pyemu.pst_utils.read_parfile(os.path.join(m_d,"pest_parchglim.par"))
    print(par)
    print(p_df)
    assert p_df.loc["stage","parval1"] == par.loc["stage","parval1"] + (rpm * par.loc["stage","parval1"])


    par.loc[pst.par_names[0],"partrans"] = "none"
    par.loc[pst.par_names[0],"parlbnd"] = -10.0
    par.loc[pst.par_names[0],"parubnd"] = 0.0   
    par.loc[pst.par_names[0],"parchglim"] = "factor"
    par.loc[pst.par_names[0],"parval1"] = -1.0
    pst.write(os.path.join(m_d,"pest_parchglim.pst"))
    pyemu.os_utils.run("{0} pest_parchglim.pst".format(exe_path.replace("-ies","-glm")),cwd=m_d)
    p_df = pyemu.pst_utils.read_parfile(os.path.join(m_d,"pest_parchglim.par"))
    print(p_df)
    print(p_df.loc["stage","parval1"],par.loc["stage","parval1"] + np.abs(par.loc["stage","parval1"] * (fpm-1)))
    assert p_df.loc["stage","parval1"] <= par.loc["stage","parval1"] + np.abs(par.loc["stage","parval1"] * (fpm-1))

    rpm = 1.1
    par.loc[pst.par_names[0],"partrans"] = "none"
    par.loc[pst.par_names[0],"parlbnd"] = -10.0
    par.loc[pst.par_names[0],"parubnd"] = 10.0   
    par.loc[pst.par_names[0],"parchglim"] = "relative"
    par.loc[pst.par_names[0],"parval1"] = -1.0
    pst.control_data.relparmax = rpm
    pst.write(os.path.join(m_d,"pest_parchglim.pst"))
    pyemu.os_utils.run("{0} pest_parchglim.pst".format(exe_path.replace("-ies","-glm")),cwd=m_d)
    p_df = pyemu.pst_utils.read_parfile(os.path.join(m_d,"pest_parchglim.par"))
    print(p_df)
    print(p_df.loc["stage","parval1"],par.loc["stage","parval1"] + rpm)
    assert np.abs(p_df.loc["stage","parval1"] - (par.loc["stage","parval1"] + rpm)) < 1.0e-6


    par.loc[pst.par_names[1:],"partrans"] = "log"
    par.loc[pst.par_names[1:],"parchglim"] = "factor"
    pst.control_data.facparmax = 5.0
    
    pst.write(os.path.join(m_d,"pest_parchglim.pst"))
    pyemu.os_utils.run("{0} pest_parchglim.pst".format(exe_path.replace("-ies","-glm")),cwd=m_d)
    p_df = pyemu.pst_utils.read_parfile(os.path.join(m_d,"pest_parchglim.par"))
    print(p_df)
    print(p_df.loc["stage","parval1"],par.loc["stage","parval1"] + rpm)
    d = np.abs(p_df.loc["stage","parval1"] - (par.loc["stage","parval1"] + rpm))
    assert d < 1.0e-6,d

    rpm = 1.1
    par.loc[pst.par_names[1:],"partrans"] = "fixed"
    par.loc[pst.par_names[1:],"parchglim"] = "factor"
    par.loc[pst.par_names[0],"partrans"] = "none"
    par.loc[pst.par_names[0],"parlbnd"] = -10.0
    par.loc[pst.par_names[0],"parubnd"] = 10.0   
    par.loc[pst.par_names[0],"parchglim"] = "relative"
    par.loc[pst.par_names[0],"parval1"] = 0.0
    pst.control_data.relparmax = rpm
    pst.write(os.path.join(m_d,"pest_parchglim.pst"))
    try:

        pyemu.os_utils.run("{0} pest_parchglim.pst".format(exe_path.replace("-ies","-glm")),cwd=m_d)
    except:
        pass
    else:
        raise Exception("should have failed")
    

    rpm = 100
    fpm = 100
    par.loc[pst.par_names[1:],"partrans"] = "fixed"
    par.loc[pst.par_names[1:],"parchglim"] = "factor"
    par.loc[pst.par_names[0],"partrans"] = "none"
    par.loc[pst.par_names[0],"parlbnd"] = 0.9
    par.loc[pst.par_names[0],"parubnd"] = 1.1   
    par.loc[pst.par_names[0],"parchglim"] = "relative"
    par.loc[pst.par_names[0],"parval1"] = 1.0
    pst.control_data.relparmax = rpm
    pst.control_data.facparmax = fpm
    
    pst.write(os.path.join(m_d,"pest_parchglim.pst"))
    pyemu.os_utils.run("{0} pest_parchglim.pst".format(exe_path.replace("-ies","-glm")),cwd=m_d)
    p_df = pyemu.pst_utils.read_parfile(os.path.join(m_d,"pest_parchglim.par"))
    print(p_df)
    assert p_df.loc["stage","parval1"] == par.loc["stage","parubnd"]

    
def sen_plusplus_test():
    model_d = "ies_10par_xsec"
    local=True
    if "linux" in platform.platform().lower() and "10par" in model_d:
        #print("travis_prep")
        #prep_for_travis(model_d)
        local=False
    
    t_d = os.path.join(model_d,"template")
    m_d = os.path.join(model_d,"master_sen_plusplus")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pst = pyemu.Pst(os.path.join(t_d,"pest.pst"))
    pst.pestpp_options = {}
    pst.pestpp_options["gsa_method"] = "morris"
    pst.pestpp_options["gsa_sobol_samples"] = 50
    pst.pestpp_options["gsa_sobol_par_dist"] = "unif"
    pst.pestpp_options["gsa_morris_r"] = 4
    pst.pestpp_options["gsa_morris_p"] = 5
    pst.pestpp_options["gsa_morris_delta"] = 2
    pst.write(os.path.join(t_d,"pest_sen.pst"))
    pyemu.os_utils.start_workers(t_d, exe_path.replace("-ies","-sen"), "pest_sen.pst", 5, master_dir=m_d,
                           worker_root=model_d,local=local,port=port)

def secondary_marker_test():
    t_d = os.path.join("secondary_marker_test","template")
    tpl_file = os.path.join(t_d,"par.dat.tpl")

    with open(tpl_file,'w') as f:
        f.write("ptf ~\n")
        f.write("~ p1    ~\n")

    tpl_file = "par.dat.tpl"
    b_d = os.getcwd()
    os.chdir(t_d)
    try:

        ins_files = [f for f in os.listdir(".") if f.endswith(".ins")]
        with open("forward_run.py",'w') as f:
            f.write("import shutil\n")
            for ins_file in ins_files:
                out_file = ins_file.replace(".ins","")
                f.write("shutil.copy2('{0}','{1}')\n".format(out_file+"_bak",out_file))

        for ins_file in ins_files:

            shutil.copy2(out_file+"_bak",out_file)
            pst = pyemu.Pst.from_io_files(tpl_file,tpl_file.replace(".tpl",""),
                ins_file,ins_file.replace(".ins",""))
            pst.control_data.noptmax = 0
            pst.pestpp_options["additional_ins_delimiters"] = "|"
            pst.model_command = "python forward_run.py"
            pst.write(os.path.join("test.pst"))

            pyemu.os_utils.run("{0} test.pst".format(exe_path))
            pst = pyemu.Pst("test.pst")
            assert pst.res is not None
            d = pst.res.loc[pst.obs_names,"modelled"] - pst.observation_data.loc[pst.obs_names,"obsval"]
            l2_d = (d.apply(np.abs)**2).sum()
            
    except Exception as e:
       os.chdir(b_d)
       raise Exception(e)
    os.chdir(b_d)

def sen_basic_test():
    local = True
    model_d = "sen_invest"
    t_d = os.path.join(model_d, "template")
    if os.path.exists(t_d):
        shutil.rmtree(t_d)
    os.makedirs(t_d)
    if "linux" in platform.platform().lower() and "10par" in model_d:
        # print("travis_prep")
        # prep_for_travis(model_d)
        local = False
    par_names = ["p1","p2"]
    obs_names = ["p1","p2","p1+p2","p1*p2","p1^p2","const"]

    tpl_file = os.path.join(t_d,"in.dat.tpl")
    with open(tpl_file,'w') as f:
        f.write("ptf ~\n")
        for par_name in par_names:
            f.write("{0}  ~     {0}      ~\n".format(par_name))
    ins_file = os.path.join(t_d,"out.dat.ins")
    with open(ins_file,'w') as f:
        f.write("pif ~\n")
        for obs_name in obs_names:
            f.write("l1 w !{0}!\n".format(obs_name))

    with open(os.path.join(t_d,"forward_run.py"),'w') as f:
        f.write("import pandas as pd\n")
        f.write("df = pd.read_csv('in.dat',index_col=0,delim_whitespace=True,names=['name','value'])\n")
        f.write("df.loc['p1+p2','value'] = df.loc['p1','value'] + df.loc['p2','value']\n")
        f.write("df.loc['p1*p2','value'] = df.loc['p1','value'] * df.loc['p2','value']\n")
        f.write("df.loc['p1^p2','value'] = df.loc['p1','value'] * df.loc['p2','value']\n")
        f.write("df.loc['const','value'] = 1.0\n")
        f.write("df.to_csv('out.dat',sep=' ',header=False)\n")

    with open(os.path.join(t_d,"in.dat"),'w') as f:
        f.write("p1 1.0\n")
        f.write("p2 1.0\n")
        f.write("p3 1.0\n")
    pyemu.os_utils.run("python forward_run.py",cwd=t_d)

    pst = pyemu.Pst.from_io_files(tpl_files=tpl_file,in_files=tpl_file.replace(".tpl",""),
                                  ins_files=ins_file,out_files=ins_file.replace(".ins",""),pst_path=".")
    pst.model_command = "python forward_run.py"
    pst.control_data.noptmax = 0
    pst.parameter_data.loc[:,"partrans"] = "log"
    pst.parameter_data.loc[:,"parchglim"] = "relative"
    pst.parameter_data.loc[:,"parubnd"] = 10.0
    pst.parameter_data.loc[:,"parlbnd"] = .1
    pst.parameter_data.loc[:,"parval1"] = 1.0

    msn_file = os.path.join(t_d,"pest.msn")
    mio_file = os.path.join(t_d, "pest.mio")

    obs = pst.observation_data
    obs.loc[:, "weight"] = 0.0
    obs.loc["const", "weight"] = 1.0
    pst.write(os.path.join(t_d, "pest.pst"))
    pyemu.os_utils.run("{0} pest.pst".format(exe_path.replace("-ies", "-sen")), cwd=t_d)
    df = pd.read_csv(msn_file, index_col=0)
    df.columns = df.columns.map(str.lower)
    df.columns = df.columns.map(str.strip)
    df.index = df.index.map(str.lower)
    print(df)
    assert df.sen_mean_abs.sum() == 0.0
    assert df.sen_std_dev.sum() == 0.0
    df = pd.read_csv(mio_file, index_col=0)
    df.columns = df.columns.map(str.lower)

    df.loc[:, "parameter_name"] = df.parameter_name.apply(str.lower)
    df.index = df.index.map(str.lower)
    print(df)
    assert df.loc[df.parameter_name == "p2", :].loc["p1", "sen_mean_abs"] == 0
    assert df.loc[df.parameter_name == "p1", :].loc["p2", "sen_mean_abs"] == 0

    pst.pestpp_options["gsa_method"] = "sobol"
    pst.pestpp_options["gsa_sobol_samples"] = 5
    pst.write(os.path.join(t_d, "pest.pst"))
    #pyemu.os_utils.run("{0} pest.pst".format(exe_path.replace("-ies", "-sen")), cwd=t_d)
    m_d = os.path.join(model_d,"master_sobol")
    pyemu.os_utils.start_workers(t_d, exe_path.replace("-ies", "-sen"), "pest.pst", 5, master_dir=m_d,
                                 worker_root=model_d, local=local, port=port)
    si_vals = pd.read_csv(os.path.join(m_d,"pest.sobol.si.csv"),index_col=0)
    sti_vals = pd.read_csv(os.path.join(m_d,"pest.sobol.sti.csv"),index_col=0)
    v_d = os.path.join("sen_invest","verf")
    si_verf_vals = pd.read_csv(os.path.join(v_d, "si.csv"), index_col=0)
    sti_verf_vals = pd.read_csv(os.path.join(v_d, "sti.csv"), index_col=0)
    d_si = (si_vals.loc[pst.obs_names,:] - si_verf_vals.loc[pst.obs_names,:]).apply(np.abs)
    print(d_si.max())
    assert d_si.max().max() < .001
    d_sti = (sti_vals.loc[pst.obs_names, :] - sti_verf_vals.loc[pst.obs_names, :]).apply(np.abs)
    print(d_sti.max())
    assert d_sti.max().max() < .001


def salib_verf():
    import pyemu
    from SALib.sample import saltelli
    from SALib.analyze import sobol
    m_d = os.path.join("sen_invest","master_sobol")
    v_d = os.path.join("sen_invest","verf")
    if os.path.exists(v_d):
        shutil.rmtree(v_d)
    os.makedirs(v_d)
    pst = pyemu.Pst(os.path.join(m_d,"pest.pst"))
    pst.add_transform_columns()
    bounds = [[l,u] for l,u in zip(pst.parameter_data.parlbnd_trans,pst.parameter_data.parubnd_trans)]
    problem = {"num_vars":pst.npar_adj,"names":pst.par_names,"bounds":bounds}
    test = saltelli.sample(problem,100,calc_second_order=False)
    out_df = pd.read_csv(os.path.join(m_d,"pest.sobol.obs.csv"),index_col=0)
    reorder_df = out_df.copy()
    idx = [0,3,1,2]
    for i in range(4):
        s = i*5
        e = s + 5
        chunk = out_df.iloc[s:e,:].copy()
        reorder_df.iloc[idx[i]::4,:] = chunk.values
        print(chunk.p1,reorder_df.iloc[idx[i]::4,:].p1)
        pass
    si_vals = pd.DataFrame(columns=pst.par_names,index=pst.obs_names)
    sti_vals = pd.DataFrame(columns=pst.par_names, index=pst.obs_names)

    for obs_name in pst.obs_names:
        #if obs_name != "p1":
        #    continue
        si = sobol.analyze(problem,reorder_df.loc[:,obs_name].values,calc_second_order=False,num_resamples=5)
        print(obs_name,si)
        si_vals.loc[obs_name,:] = si["S1"]
        sti_vals.loc[obs_name, :] = si["ST"]
    si_vals.to_csv(os.path.join(v_d,"si.csv"))
    sti_vals.to_csv(os.path.join(v_d, "sti.csv"))

    # in_df = pd.read_csv(os.path.join(m_d,"pest.sobol.par.csv"),index_col=0)
    # import matplotlib.pyplot as plt
    # fig, ax = plt.subplots(1,1)
    # print(test)
    # test = test ** 10
    # print(test)
    # ax.scatter(test[:,0],test[:,1],marker='.',color='g')
    # ax.scatter(in_df.iloc[:,0],in_df.iloc[:,1],marker='.',color='r')
    #
    #
    # plt.show()


def tplins1_test():
    model_d = "tplins_test_1"
    t_d = os.path.join(model_d, "test")
    if os.path.exists(t_d):
        shutil.rmtree(t_d)
    shutil.copytree(os.path.join(model_d,"template"),t_d)
    pst = pyemu.Pst(os.path.join(t_d,"pest.pst"))
    ins_file = os.path.join(t_d,"AOC_obs.txt.ins")
    pst.add_observations(ins_file,ins_file.replace(".ins",""),pst_path=".")
    
    pyemu.os_utils.run("{0} pest.pst".format(exe_path.replace("-ies","-glm")),cwd=t_d)
    obf_df = pd.read_csv(os.path.join(t_d,"out1.dat.obf"),delim_whitespace=True,header=None,names=["obsnme","obsval"])
    obf_df.index = obf_df.obsnme
    pst = pyemu.Pst(os.path.join(t_d,"pest.pst"))
    res_df = pst.res
    
    d = (obf_df.obsval - res_df.modelled).apply(np.abs)
    #print(d)
    print(d.max())
    assert d.max() < 1.0e-5, d

    jco = pyemu.Jco.from_binary(os.path.join(t_d,"pest.jcb")).to_dataframe().apply(np.abs)
    assert jco.sum().sum() == 0, jco.sum()

    # check the input file - the last two number should be the same
    arr = np.loadtxt(os.path.join(t_d,"hk_Layer_1.ref"))
    assert arr[-2] == arr[-1]

    lines_tpl = open(os.path.join(t_d,"hk_Layer_1.ref.tpl"),'r').readlines()
    lines_in = open(os.path.join(t_d,"hk_Layer_1.ref"),'r').readlines()
    assert len(lines_tpl) - 1 == len(lines_in)

    pst = pyemu.Pst(os.path.join(t_d,"pest.pst"))
    pst.pestpp_options["fill_tpl_zeros"] = True
    pst.write(os.path.join(t_d,"pest_fill.pst"))
    pyemu.os_utils.run("{0} pest_fill.pst".format(exe_path.replace("-ies","-glm")),cwd=t_d)
    obf_df = pd.read_csv(os.path.join(t_d,"out1.dat.obf"),delim_whitespace=True,header=None,names=["obsnme","obsval"])
    obf_df.index = obf_df.obsnme
    pst = pyemu.Pst(os.path.join(t_d,"pest_fill.pst"))
    res_df = pst.res
    
    d = (obf_df.obsval - res_df.modelled).apply(np.abs)
    #print(d)
    print(d.max())
    assert d.max() < 1.0e-5, d

    jco = pyemu.Jco.from_binary(os.path.join(t_d,"pest_fill.jcb")).to_dataframe().apply(np.abs)
    assert jco.sum().sum() == 0, jco.sum()

    # check the input file - the last two number should be the same
    arr = np.loadtxt(os.path.join(t_d,"hk_Layer_1.ref"))
    assert arr[-2] == arr[-1]

    lines_tpl = open(os.path.join(t_d,"hk_Layer_1.ref.tpl"),'r').readlines()
    lines_in = open(os.path.join(t_d,"hk_Layer_1.ref"),'r').readlines()
    assert len(lines_tpl) - 1 == len(lines_in)

    pst = pyemu.Pst(os.path.join(t_d, "pest.pst"))
    dum_obs = ['h01_03', 'h01_07']
    pst.observation_data.drop(index=dum_obs, inplace=True)
    pst.model_output_data = pd.DataFrame({"pest_file":"out1dum.dat.ins",
                                          "model_file":'out1.dat'},index=["out1dum.dat.ins"])
    #pst.instruction_files = ['out1dum.dat.ins']
    pst.write(os.path.join(t_d, "pest_dum.pst"))
    pyemu.os_utils.run("{0} pest_dum.pst".format(exe_path.replace("-ies", "-glm")), cwd=t_d)
    obf_df = pd.read_csv(os.path.join(t_d, "out1.dat.obf"), delim_whitespace=True, header=None,
                         names=["obsnme", "obsval"])
    obf_df.index = obf_df.obsnme
    pst = pyemu.Pst(os.path.join(t_d, "pest_dum.pst"))
    res_df = pst.res

    d = (obf_df.obsval - res_df.modelled).apply(np.abs)
    # print(d)
    print(d.max())
    assert d.max() < 1.0e-5, d

    jco = pyemu.Jco.from_binary(os.path.join(t_d, "pest_dum.jcb")).to_dataframe().apply(np.abs)
    assert jco.sum().sum() == 0, jco.sum()



def ext_stdcol_test():
    model_d = "ies_10par_xsec"
    local=True
    if "linux" in platform.platform().lower() and "10par" in model_d:
        #print("travis_prep")
        #prep_for_travis(model_d)
        local=False
    
    t_d = os.path.join(model_d,"template")
    m_d = os.path.join(model_d,"master_ext_stdcol")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    shutil.copytree(t_d,m_d)
    

    pst = pyemu.Pst(os.path.join(m_d,"pest.pst"))  
    obs = pst.observation_data
    obs.loc[pst.nnz_obs_names,"standard_deviation"] = 1/obs.loc[pst.nnz_obs_names,"weight"]
    pst.add_transform_columns()
    par = pst.parameter_data
    par.loc[pst.adj_par_names,"standard_deviation"] = (par.loc[pst.adj_par_names,"parubnd_trans"] - par.loc[pst.adj_par_names,"parlbnd_trans"]) / 4.0
    #par.loc[pst.adj_par_names[0],"mean"] = par.loc[pst.adj_par_names[0],"parubnd"]
    pst.pestpp_options["ies_num_reals"] = 10
    pst.control_data.noptmax = -1
    pst.write(os.path.join(m_d,"pest_base.pst"))
    pyemu.os_utils.run("{0} pest_base.pst".format(exe_path),cwd=m_d)

    pst.write(os.path.join(m_d,"pest_ext_stdcol.pst"),version=2)
    pyemu.os_utils.run("{0} pest_ext_stdcol.pst".format(exe_path),cwd=m_d)
    df1 = pd.read_csv(os.path.join(m_d,"pest_base.phi.meas.csv"),index_col=0)
    df2 = pd.read_csv(os.path.join(m_d,"pest_ext_stdcol.phi.meas.csv"),index_col=0)

    d = (df1 - df2).apply(np.abs)
    print(d.max())
    assert d.max().max() < 1.0e-6,d.max().max()

    pst.pestpp_options["ies_num_reals"] = 100000
    pst.control_data.noptmax = -2
    obs = pst.observation_data
    obs.loc[pst.nnz_obs_names,"standard_deviation"] = 7.5
    pst.write(os.path.join(m_d,"pest_ext_stdcol.pst"),version=2)
    pyemu.os_utils.run("{0} pest_ext_stdcol.pst".format(exe_path),cwd=m_d)
    df = pd.read_csv(os.path.join(m_d,"pest_ext_stdcol.obs+noise.csv"),index_col=0).loc[:,pst.nnz_obs_names]
    d = (df.std() - obs.loc[pst.nnz_obs_names,"standard_deviation"]).apply(np.abs)
    print(d)
    assert d.max() < 0.1,d.max()
    
    obs = pst.observation_data
    obs.loc[pst.nnz_obs_names,"upper_bound"] = obs.loc[pst.nnz_obs_names,"obsval"] * 1.1
    obs.loc[pst.nnz_obs_names,"lower_bound"] = obs.loc[pst.nnz_obs_names,"obsval"] * 0.9
    par = pst.parameter_data
    par.loc[pst.adj_par_names[0],"mean"] = par.loc[pst.adj_par_names[0],"parubnd"]
    pst.write(os.path.join(m_d,"pest_ext_stdcol.pst"),version=2)
    pyemu.os_utils.run("{0} pest_ext_stdcol.pst".format(exe_path),cwd=m_d)
    df = pd.read_csv(os.path.join(m_d,"pest_ext_stdcol.obs+noise.csv"),index_col=0).loc[:,pst.nnz_obs_names]
    mn = df.min()
    mx = df.max()
    dmn = mn - obs.loc[pst.nnz_obs_names,"obsval"] * 0.9
    print(obs.loc[pst.nnz_obs_names,"obsval"] * 0.9)
    print(mn)  
    print(dmn)
    dmx = mx - obs.loc[pst.nnz_obs_names,"obsval"] * 1.1
    print(obs.loc[pst.nnz_obs_names,"obsval"] * 1.1)
    print(mx)
    print(dmx)

    dmn = dmn.apply(np.abs)
    dmx = dmx.apply(np.abs)

    assert dmn.max() < 1.0e-6,dmn
    assert dmx.max() < 1.0e-6,dmx


def mf6_v5_ies_test():
    model_d = "mf6_freyberg"
    local=True
    if "linux" in platform.platform().lower() and "10par" in model_d:
        #print("travis_prep")
        #prep_for_travis(model_d)
        local=False
    
    t_d = os.path.join(model_d,"template")
    m_d = os.path.join(model_d,"master_ies")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pst = pyemu.Pst(os.path.join(t_d,"freyberg6_run_ies.pst"))
    pst.control_data.noptmax = 0
    pst.write(os.path.join(t_d,"freyberg6_run_ies.pst"))
    pyemu.os_utils.run("pestpp-ies freyberg6_run_ies.pst",cwd=t_d)

    pst.control_data.noptmax = 2
    pst.write(os.path.join(t_d,"freyberg6_run_ies.pst"))
    pyemu.os_utils.start_workers(t_d, "pestpp-ies", "freyberg6_run_ies.pst", num_workers=15,
                                master_dir=m_d,worker_root=model_d,port=port)

    
    oe_file = os.path.join(m_d,"freyberg6_run_ies.{0}.obs.csv".format(pst.control_data.noptmax))
    assert os.path.exists(oe_file)
    pe_file = oe_file.replace(".obs.",".par.")
    assert os.path.exists(pe_file)
    pcs_file = oe_file.replace(".obs.",".pcs.")
    assert os.path.exists(pcs_file)
    df = pd.read_csv(pcs_file,index_col=0)
    pst_pargp = set(list(pst.parameter_data.pargp.unique()))
    df_pargp = set(df.index.to_list())
    d = pst_pargp.symmetric_difference(df_pargp)
    print(d)
    assert len(d) == 0,d


def mf6_v5_sen_test():
    model_d = "mf6_freyberg"
    local=True
    if "linux" in platform.platform().lower() and "10par" in model_d:
        #print("travis_prep")
        #prep_for_travis(model_d)
        local=False
    
    t_d = os.path.join(model_d,"template")
    m_d = os.path.join(model_d,"master_sen")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pst = pyemu.Pst(os.path.join(t_d,"freyberg6_run_sen.pst"))
    m_d = os.path.join(model_d,"master_sen")
    pyemu.os_utils.start_workers(t_d, "pestpp-sen", "freyberg6_run_sen.pst",
                                 num_workers=15, master_dir=m_d, worker_root=model_d,
                                 port=port)

    pst = pyemu.Pst(os.path.join(m_d,"freyberg6_run_sen.pst"))
    mio_file = os.path.join(m_d,"freyberg6_run_sen.mio")
    assert os.path.exists(mio_file)
    df = pd.read_csv(mio_file)
    assert df.shape[0] > 1
    msn_file = mio_file.replace(".mio",".msn")
    assert os.path.exists(msn_file)

def mf6_v5_opt_stack_test():
    model_d = "mf6_freyberg"
    local=True
    if "linux" in platform.platform().lower() and "10par" in model_d:
        #print("travis_prep")
        #prep_for_travis(model_d)
        local=False
    
    t_d = os.path.join(model_d,"template")
    m_d = os.path.join(model_d,"master_opt_stack")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pst = pyemu.Pst(os.path.join(t_d,"freyberg6_run_opt.pst"))
    m_d = os.path.join(model_d,"master_opt_stack")
    pyemu.os_utils.start_workers(t_d, "pestpp-opt", "freyberg6_run_opt.pst", 
                                 num_workers=15, master_dir=m_d,worker_root=model_d,
                                 port=port)

    assert os.path.exists(os.path.join(m_d,"freyberg6_run_opt.1.sim+chance.rei"))
    assert os.path.exists(os.path.join(m_d,"freyberg6_run_opt.1.obs_stack.csv"))


def mf6_v5_glm_test():
    model_d = "mf6_freyberg"
    local=True
    if "linux" in platform.platform().lower() and "10par" in model_d:
        #print("travis_prep")
        #prep_for_travis(model_d)
        local=False
    
    t_d = os.path.join(model_d,"template")
    m_d = os.path.join(model_d,"master_glm")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pst = pyemu.Pst(os.path.join(t_d,"freyberg6_run_glm.pst"))
    m_d = os.path.join(model_d,"master_glm")
    pyemu.os_utils.start_workers(t_d, "pestpp-glm", "freyberg6_run_glm.pst", 
                                 num_workers=15, master_dir=m_d,worker_root=model_d,
                                 port=port)

    oe_file = os.path.join(m_d,"freyberg6_run_glm.post.obsen.csv")
    assert os.path.exists(oe_file)
    oe = pd.read_csv(oe_file)
    assert oe.shape[0] == pst.pestpp_options["glm_num_reals"],"{0},{1}".\
        format(oe.shape[0],pst.pestpp_options["glm_num_reals"])


def cmdline_test():
    model_d = "mf6_freyberg"
    local=True
    if "linux" in platform.platform().lower() and "10par" in model_d:
        #print("travis_prep")
        #prep_for_travis(model_d)
        local=False
    
    t_d = os.path.join(model_d,"template")
    pst_name = "freyberg6_run_glm.pst"
    pst = pyemu.Pst(os.path.join(t_d,"freyberg6_run_glm.pst"))
    pst.pestpp_options["debug_parse_only"] = True
    pst_name = "CmdLine_test.pst" #camel case on purpose for linux testing
    pst.write(os.path.join(t_d,pst_name))
    pyemu.os_utils.run("{0} {1}".format(exe_path,pst_name),cwd=t_d)
    pyemu.os_utils.run("{0} {1} /h :4004".format(exe_path,pst_name),cwd=t_d)
    pyemu.os_utils.run("{0} {1} /r /h :4004".format(exe_path.replace("-ies","-glm"),pst_name),cwd=t_d) 
    pyemu.os_utils.run("{0} {1} /r ".format(exe_path.replace("-ies","-glm"),pst_name),cwd=t_d) 
    
    try:
        pyemu.os_utils.run("{0} {1} \\h :4004".format(exe_path,pst_name),cwd=t_d) 
       
    except:
        pass
    else:
        raise Exception("should have failed")
    
    try:
        pyemu.os_utils.run("{0} {1} :4004".format(exe_path,pst_name),cwd=t_d) 
       
    except:
        pass
    else:
        raise Exception("should have failed")

    try:
        pyemu.os_utils.run("{0} {1} /h 4004".format(exe_path,pst_name),cwd=t_d) 
       
    except:
        pass
    else:
        raise Exception("should have failed")
    

def da_prep_4_mf6_freyberg_seq():
    t_d = os.path.join("mf6_freyberg","template_seq")
    if os.path.exists(t_d):
        shutil.rmtree(t_d)
    shutil.copytree(os.path.join("mf6_freyberg","template"),t_d)
    for f in os.listdir(t_d):
        for tag in ["ies","opt","glm"]:
            if tag in f.lower():
                os.remove(os.path.join(t_d,f))

    #first modify the tdis
    with open(os.path.join(t_d,"freyberg6.tdis"),'w') as f:
        f.write("BEGIN Options\n  TIME_UNITS  days\nEND Options\n")
        f.write("BEGIN Dimensions\n  NPER  1\nEND Dimensions\n")
        f.write("BEGIN PERIODDATA\n31.00000000  1       1.00000000\nEND PERIODDATA\n")
    #make sure it runs
    pyemu.os_utils.run("mf6",cwd=t_d)

    # write a tdis template file - could possibly keep all 25 stress periods to
    # simulate a 2-year-ahead forecast...

    with open(os.path.join(t_d,"freyberg6.tdis.tpl"),'w') as f:
        f.write("ptf  ~\n")
        f.write("BEGIN Options\n  TIME_UNITS  days\nEND Options\n")
        f.write("BEGIN Dimensions\n  NPER  1\nEND Dimensions\n")
        f.write("BEGIN PERIODDATA\n~  perlen  ~  1       1.00000000\nEND PERIODDATA\n")
    new_tpl,new_in = [os.path.join(t_d,"freyberg6.tdis.tpl")],[os.path.join(t_d,"freyberg6.tdis")]
    new_tpl_cycle = [-1]

    # mod the sto to make sp 1 transient - or should this be a pre-processor so that cycle 0 is 
    # ss and the rest are transient?

    # split out the head, sfr and list instruction files into multiple instruction file
    #eventually, want to move to da par and obs cycle tables for heads and gage_1 obs
    lines = open(os.path.join(t_d,"heads.csv.ins"),'r').readlines()[2:]
    new_ins,new_out,new_ins_cycle = [],[],[]
    #print(lines)
    for icycle, line in enumerate(lines):
        ins_name = os.path.join(t_d,"heads_{0}.csv.ins".format(icycle))
        with open(ins_name,'w') as f:
            f.write("pif ~\nl1\n")
            f.write(line)
        new_ins.append(ins_name)
        new_out.append(os.path.join(t_d,"heads.csv"))
        new_ins_cycle.append(icycle)
    remove_ins = ["heads.csv.ins"]

    lines = open(os.path.join(t_d,"sfr.csv.ins"),'r').readlines()[2:]
    #print(lines)
    for icycle, line in enumerate(lines):
        ins_name = os.path.join(t_d,"sfr_{0}.csv.ins".format(icycle))
        with open(ins_name,'w') as f:
            f.write("pif ~\nl1\n")
            f.write(line)
        new_ins.append(ins_name)
        new_out.append(os.path.join(t_d,"sfr.csv"))
        new_ins_cycle.append(icycle)
    remove_ins.append("sfr.csv.ins")

    lines = open(os.path.join(t_d,"freyberg6.lst.ins"),'r').readlines()[1:]
    icycle = 0
    tag_line = lines[0]
    for s in range(0,len(lines),13):
        ins_name = os.path.join(t_d,"freyberg6_{0}.lst.ins".format(icycle))
        with open(os.path.join(t_d,"freyberg6_{0}.lst.ins".format(icycle)),'w') as f:
            f.write("pif ~\n")
            f.write(tag_line)
            for line in lines[s+1:s+13]:
                f.write(line)
        new_ins.append(ins_name)
        new_out.append(os.path.join(t_d,"freyberg6.lst"))
        new_ins_cycle.append(icycle)
        icycle += 1
    remove_ins.append("freyberg6.lst.ins")

    # modify the ic file
    k = 0
    with open(os.path.join(t_d,"freyberg6.ic"),'r') as f:
        while True:
            line = f.readline()
            if line == "":
                break
            if line.lower().strip().startswith("internal"):
                arr_lines = []
                while True:
                    line = f.readline()
                    if line == "":
                        raise Exception
                    if line.lower().strip().startswith("end"):
                        break
                    if line.lower().strip().startswith("internal"):
                        arr = np.array(arr_lines,dtype=np.float)
                        np.savetxt(os.path.join(t_d,"heads_{0}.dat_in".format(k)),arr,fmt="%15.6E")
                        k += 1
                        arr_lines = []
                    else:
                        arr_lines.append(line.strip().split())

        arr = np.array(arr_lines, dtype=np.float)
        np.savetxt(os.path.join(t_d, "heads_{0}.dat_in".format(k)), arr, fmt="%15.6E")
    with open(os.path.join(t_d,"freyberg6.ic"),'w') as f:
        f.write("begin griddata\nstrt layered\n")
        for k in range(3):
            f.write("open/close 'heads_{0}.dat_in' FACTOR 1.0\n".format(k))
        f.write("end griddata\n")


    # write a python script to extract final heads and save to files
    with open(os.path.join(t_d,"forward_run.py"),'w') as f:
        f.write("import numpy as np\nimport flopy\nimport pyemu\n")
        f.write("pyemu.os_utils.run('mf6')\n")
        f.write("hds = flopy.utils.HeadFile('freyberg6_freyberg.hds')\n")
        f.write("arr = hds.get_data()\n")
        f.write("for k,a in enumerate(arr):\n")
        f.write("    np.savetxt('heads_'+str(k)+'.dat',a,fmt='%15.6E')\n")

    # dont run it so we can harvest the ic values in the arrays for setting the parval1 values
    pyemu.os_utils.run("python forward_run.py",cwd=t_d)

    # now write ins and tpl file for these
    ic_parvals = {}
    for k in range(3):
        fname = os.path.join(t_d,"heads_{0}.dat".format(k))
        assert os.path.exists(fname),fname
        arr = np.loadtxt(fname)
        fname_ins = fname + "_out.ins"
        fname_tpl = fname + "_in.tpl"
        in_arr = np.loadtxt(fname_tpl.replace(".tpl",""))
        ft = open(fname_tpl,'w')
        with open(fname_ins,'w') as f:
            f.write("pif ~\n")
            ft.write("ptf ~\n")
            for i in range(arr.shape[0]):
                f.write("l1 ")
                for j in range(arr.shape[1]):
                    if np.abs(arr[i,j]) > 100 or np.abs(in_arr[i,j]) > 100:
                        f.write(" !dum! ")
                        ft.write(" 40 ")
                    else:
                        oname = "head_{0:02d}_{1:03d}_{2:03d}".format(k,i,j)
                        f.write(" !{0}! ".format(oname))
                        ft.write(" ~  {0} ~ ".format(oname))
                        ic_parvals[oname] = in_arr[i,j]
                f.write("\n")
                ft.write("\n")
        ft.close()
        new_tpl.append(fname_tpl)
        new_in.append(fname_tpl.replace(".tpl",""))
        new_tpl_cycle.append(-1)
        new_ins.append(fname_ins)
        new_out.append(os.path.join(t_d,"heads_{0}.dat".format(k)))
        new_ins_cycle.append(-1)

        i = pyemu.pst_utils.InstructionFile(fname_ins)
        df = i.read_output_file(fname)
        #print(df)

    # split out the wel and rch tpl files into cycle files
    lines = []
    with open(os.path.join(t_d,"freyberg6.wel.tpl"),'r') as f:
        for i in range(19):
            lines.append(f.readline())
    print(lines)

    for icycle in range(25):
        tpl_file = os.path.join(t_d,"freyberg6.wel_{0}.tpl".format(icycle))
        with open(tpl_file,'w') as f:
            for line in lines:
                new_line = line.replace("_0","_{0}".format(icycle))
                f.write(new_line)

        new_tpl.append(tpl_file)
        new_in.append(os.path.join(t_d,"freyberg6.wel"))
        new_tpl_cycle.append(icycle)
    remove_tpl = ["freyberg6.wel.tpl"]

    lines = []
    with open(os.path.join(t_d, "freyberg6.rch.tpl"), 'r') as f:
        for i in range(11):
            lines.append(f.readline())
    print(lines)

    for icycle in range(25):
        tpl_file = os.path.join(t_d,"freyberg6.rch_{0}.tpl".format(icycle))
        with open(tpl_file,'w') as f:
            for line in lines:
                new_line = line.replace("_0","_{0}".format(icycle))
                f.write(new_line)

        new_tpl.append(tpl_file)
        new_in.append(os.path.join(t_d,"freyberg6.rch"))
        new_tpl_cycle.append(icycle)
    remove_tpl.append('freyberg6.rch.tpl')

    # now for the fun part: modify the pst
    pst = pyemu.Pst(os.path.join(t_d,"freyberg6_run.pst"))

    print(pst.npar_adj,pst.nnz_obs)

    # swap out obs info
    dropped_dfs = []
    for ins in remove_ins:
        dropped_dfs.append(pst.drop_observations(os.path.join(t_d, ins), '.'))
    for insf, outf, cy in zip(new_ins, new_out, new_ins_cycle):
        df = pst.add_observations(insf,outf, pst_path=".")
        pst.observation_data.loc[df.obsnme, "cycle"] = cy
        pst.model_output_data.loc[os.path.join(".",os.path.split(insf)[1]),"cycle"] = cy
    pst.observation_data.loc[:,"weight"] = 0.0
    for df in dropped_dfs:
        for c in ["obsval","weight"]:
            pst.observation_data.loc[df.obsnme, c] = df.loc[:, c]

    # swap out par info
    dropped_dfs = []
    for tpl in remove_tpl:
        dropped_dfs.append(pst.drop_parameters(os.path.join(t_d,tpl),'.'))
    pst.parameter_data.loc[:,"cycle"] = -1
    pst.model_input_data.loc[:,"cycle"] = -1
    for tplf, inf, cy in zip(new_tpl,new_in,new_tpl_cycle):
        df = pst.add_parameters(tplf,inf,pst_path=".")
        pst.parameter_data.loc[df.parnme,"cycle"] = cy
        pst.model_input_data.loc[os.path.join(".",os.path.split(tplf)[1]),"cycle"] = cy
    for df in dropped_dfs:
        for c in ["parval1","parubnd","parlbnd","pargp"]:
            pst.parameter_data.loc[df.parnme,c] = df.loc[:,c]

    #set the state parameter info
    for p,v in ic_parvals.items():
        pst.parameter_data.loc[p,"parval1"] = v
        pst.parameter_data.loc[p, "parlbnd"] = v * 0.5
        pst.parameter_data.loc[p, "parubnd"] = v * 1.5
        pst.parameter_data.loc[p,"pargp"] = "head_state"
        pst.parameter_data.loc[p,"partrans"] = "none"

    pst.control_data.noptmax = 2
    pst.model_command = "python forward_run.py"
    pst.pestpp_options.pop("ies_par_en")
    pst.parameter_data.loc["perlen","partrans"] = "fixed"
    pst.pestpp_options["ies_num_reals"] = 5
    pst.pestpp_options["da_num_reals"] = 5
    pst.write(os.path.join(t_d,"freyberg6_run_da1.pst"),version=2)
    return pst

def da_mf6_freyberg_test_1():
    pst = da_prep_4_mf6_freyberg_seq()
    test_d = "mf6_freyberg"
    t_d = os.path.join(test_d,"template_seq")
    print(exe_path.replace("ies","da"))
    pst.control_data.noptmax = 0
    pst.pestpp_options["ies_verbose_level"] = 4
    pst.write(os.path.join(t_d, "freyberg6_run_da1.pst"), version=2)
    pyemu.os_utils.run("{0} freyberg6_run_da1.pst".format(exe_path.replace("ies","da")),cwd=t_d)
    #pst.control_data.noptmax = 2
    # pst.write(os.path.join(t_d, "freyberg6_run_da1.pst"), version=2)
    # pyemu.os_utils.start_workers(t_d,exe_path.replace("ies","da"),"freyberg6_run_da1.pst",
    #                              num_workers=5,worker_root=test_d,port=port,
    #                              master_dir=os.path.join(test_d,"master_da_1"),verbose=True)


def da_prep_4_mf6_freyberg_seq_tbl():
    test_d = "mf6_freyberg"
    t_d = os.path.join(test_d,"template_seq")
    assert os.path.exists(t_d)
    pst = pyemu.Pst(os.path.join(t_d,"freyberg6_run_da1.pst"))
    pdf = pd.DataFrame({"perlen":31},index=np.arange(25))
    pdf.T.to_csv(os.path.join(t_d,"par_cycle_tbl.csv"))
    pst.pestpp_options["da_parameter_cycle_table"] = "par_cycle_tbl.csv"

    # mod sfr_out
    sfr_ins_file = os.path.join(t_d, "sfr.csv.ins")
    with open(sfr_ins_file, 'w') as f:
        f.write("pif ~\n")
        f.write("l1\n")
        f.write("l1 ~,~ !headwater!  ~,~ !tailwater!  ~,~ !gage_1!\n")

    # and lst budget
    lines = open(os.path.join(t_d, "freyberg6_0.lst.ins"), 'r').readlines()
    lst_ins_file = os.path.join(t_d, "freyberg6.lst.ins")
    with open(lst_ins_file, 'w') as f:
        for line in lines:
            f.write(line.replace("_20151231",""))

    obs = pst.observation_data.copy()
    tr_obs = obs.loc[obs.obsnme.str.startswith("trgw"),:].copy()
    tr_obs.loc[tr_obs.obsnme,"datetime"] = pd.to_datetime(tr_obs.obsnme.apply(lambda x: x.split('_')[-1]))
    tr_obs.loc[tr_obs.obsnme,"k"] = tr_obs.obsnme.apply(lambda x: np.int(x.split('_')[1]))
    tr_obs.loc[tr_obs.obsnme, "i"] = tr_obs.obsnme.apply(lambda x: np.int(x.split('_')[2]))
    tr_obs.loc[tr_obs.obsnme, "j"] = tr_obs.obsnme.apply(lambda x: np.int(x.split('_')[3]))
    tr_obs.loc[tr_obs.obsnme,"obgnme"] = tr_obs.obsnme.apply(lambda x: "_".join(x.split("_")[:-1]))

    head_obs = obs.loc[obs.obsnme.str.startswith("head_"),:].copy()
    head_obs.loc[head_obs.obsnme, "k"] = head_obs.obsnme.apply(lambda x: np.int(x.split('_')[1]))
    head_obs.loc[head_obs.obsnme, "i"] = head_obs.obsnme.apply(lambda x: np.int(x.split('_')[2]))
    head_obs.loc[head_obs.obsnme, "j"] = head_obs.obsnme.apply(lambda x: np.int(x.split('_')[3]))

    #print(pst.nobs)
    for ins_file in pst.model_output_data.pest_file.copy():
        if "heads_" in ins_file and ins_file.endswith("dat_out.ins"):
            continue
        pst.drop_observations(os.path.join(t_d,ins_file),pst_path=".")
    for ins_file in [sfr_ins_file,lst_ins_file]:
        pst.add_observations(ins_file,pst_path=".")

    # work out which heads are obs site
    obs_heads = {}
    odf_names = []
    pst.observation_data.loc[:,"org_obgnme"] = np.NaN
    pst.observation_data.loc[:, "weight"] = 0.0
    for og in tr_obs.obgnme.unique():
        site_obs = tr_obs.loc[tr_obs.obgnme==og,:]
        site_obs.sort_values(by="datetime",inplace=True)
        head_name = "head_{0:02d}_{1:03d}_{2:03d}".format(site_obs.k[0],site_obs.i[0],site_obs.j[0])
        for i,oname in enumerate(site_obs.obsnme):
            obs_heads[oname] = (head_name,i)
        # assign max weight in the control file since some have zero weight and
        # we are covering weights in the weight table
        #print(og,site_obs.weight.max())
        pst.observation_data.loc[head_name,"weight"] = site_obs.weight.max()
        pst.observation_data.loc[head_name,"org_obgnme"] = og
        odf_names.append(head_name)

    odf_names.append("gage_1")

    odf = pd.DataFrame(columns=odf_names,index=np.arange(25))
    wdf = pd.DataFrame(columns=odf_names,index=np.arange(25))
    for tr_name,(head_name,icycle) in obs_heads.items():
        odf.loc[icycle,head_name] = obs.loc[tr_name,"obsval"]
        wdf.loc[icycle, head_name] = obs.loc[tr_name, "weight"]

    g_obs = obs.loc[obs.obsnme.str.startswith("gage_1"),:].copy()
    #give these obs the max weight since some have zero weight
    pst.observation_data.loc["gage_1", "weight"] = g_obs.weight.max()
    g_obs.sort_index(inplace=True)
    for i,name in enumerate(g_obs.obsnme):
        odf.loc[i,"gage_1"] = g_obs.loc[name,"obsval"]
        wdf.loc[i, "gage_1"] = g_obs.loc[name, "weight"]

    # now drop any entries that have zero weight across all cycles
    print(pst.nnz_obs_names)
    odf = odf.loc[:,pst.nnz_obs_names]
    wdf = wdf.loc[:, pst.nnz_obs_names]

    odf.T.to_csv(os.path.join(t_d,"obs_cycle_tbl.csv"))
    pst.pestpp_options["da_observation_cycle_table"] = "obs_cycle_tbl.csv"
    wdf.T.to_csv(os.path.join(t_d, "weight_cycle_tbl.csv"))
    pst.pestpp_options["da_weight_cycle_table"] = "weight_cycle_tbl.csv"

    pst.observation_data.loc[:,"cycle"] = -1
    pst.model_output_data.loc[:,"cycle"] = -1

    pst.write(os.path.join(t_d,"freyberg6_run_da2.pst"),version=2)
    return pst

def da_mf6_freyberg_test_2():
    pst = da_prep_4_mf6_freyberg_seq_tbl()
    test_d = "mf6_freyberg"
    t_d = os.path.join(test_d, "template_seq")
    print(exe_path.replace("ies", "da"))
    pst.control_data.noptmax = 0
    pst.pestpp_options["ies_verbose_level"] = 4
    pst.write(os.path.join(t_d, "freyberg6_run_da2.pst"), version=2)
    pyemu.os_utils.run("{0} freyberg6_run_da2.pst".format(exe_path.replace("ies","da")),cwd=t_d)
    return
    pst.control_data.noptmax = 2
    pst.write(os.path.join(t_d, "freyberg6_run_da2.pst"), version=2)
    pyemu.os_utils.start_workers(t_d, exe_path.replace("ies", "da"), "freyberg6_run_da2.pst",
                                num_workers=5, worker_root=test_d, port=port,
                                master_dir=os.path.join(test_d, "master_da_2"), verbose=True)


def invest():
    t_d = os.path.join("mf6_freyberg","template_seq")
    pst1 = pyemu.Pst(os.path.join(t_d,"freyberg6_run_da1.pst"))
    pst2 = pyemu.Pst(os.path.join(t_d,"freyberg6_run_da2.pst"))
    print(pst1.observation_data.loc[pst1.nnz_obs_names,["obsval","weight"]])
    odf = pd.read_csv(os.path.join(t_d,pst2.pestpp_options["da_observation_cycle_table"]))
    wdf = pd.read_csv(os.path.join(t_d, pst2.pestpp_options["da_weight_cycle_table"]))
    print(odf)
    print(wdf)


if __name__ == "__main__":
    
    #glm_long_name_test()
    #sen_plusplus_test()
    #parchglim_test()
    #unc_file_test()
    #secondary_marker_test()
    #basic_test("ies_10par_xsec")
    #glm_save_binary_test()
    #sweep_forgive_test()
    #inv_regul_test()
    #tie_by_group_test()
    #sen_basic_test()
    #salib_verf()
    #tplins1_test()
    #ext_stdcol_test()

    # parallel_consist_test()
    # ext_stdcol_test()
    #da_prep_4_freyberg_batch()
    #da_prep_4_mf6_freyberg_seq()
    shutil.copy2(os.path.join("..","exe","windows","x64","Debug","pestpp-da.exe"),os.path.join("..","bin","pestpp-da.exe"))
    da_mf6_freyberg_test_1()

    #da_prep_4_mf6_freyberg_seq_tbl()
    da_mf6_freyberg_test_2()
    #mf6_v5_ies_test()
    #mf6_v5_sen_test()
    #mf6_v5_opt_stack_test()
    #mf6_v5_glm_test()
    #cmdline_test()
    invest()
