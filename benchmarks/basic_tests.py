import os
import sys
import shutil
import platform
import numpy as np
import pandas as pd
import platform
import pyemu

bin_path = os.path.join("test_bin")
plat = "unknown"
if "linux" in platform.platform().lower():
    bin_path = os.path.join(bin_path,"linux")
    plat = "linux"
elif "darwin" in platform.platform().lower() or "macos" in platform.platform().lower() :
    bin_path = os.path.join(bin_path,"mac")
    plat = "apple"
else:
    bin_path = os.path.join(bin_path,"win")
    plat = "windows"

bin_path = os.path.abspath("test_bin")
os.environ["PATH"] += os.pathsep + bin_path


# bin_path = os.path.join("..","..","..","bin")
# exe = ""
# if "windows" in platform.platform().lower():
#     exe = ".exe"
# exe_path = os.path.join(bin_path, "pestpp-ies" + exe)

# case of either appveyor, travis or local
if os.path.exists(os.path.join("pestpp","bin")):
    bin_path = os.path.join("..","..","pestpp","bin")
else:
    bin_path = os.path.join("..","..","..","..","pestpp","bin")

        
if "windows" in platform.platform().lower():
    exe_path = os.path.join(bin_path, "win", "pestpp-ies.exe")
elif "darwin" in platform.platform().lower() or "macos" in platform.platform().lower() :
    exe_path = os.path.join(bin_path,  "mac", "pestpp-ies")
else:
    exe_path = os.path.join(bin_path, "linux", "pestpp-ies")

noptmax = 4
num_reals = 20
port = 4021


def nonascii_path_test(model_d="ies_10par_xsec"):
    pyemu.Ensemble.reseed()
    base_d = os.path.join(model_d, "template")
    new_d = os.path.join(model_d, "test_template_\u0187")
    if os.path.exists(new_d):
        shutil.rmtree(new_d)
    shutil.copytree(base_d, new_d)
    print(platform.platform().lower())
    pst = pyemu.Pst(os.path.join(new_d, "pest.pst"))
    cmd = pst.model_command[0].split()
    print(cmd)
    cmd = "\"\"{0}\" \"{1}\"\"".format(cmd[0],cmd[1])
    print(cmd)
    pst.model_command.append(cmd)
    cmd = pst.model_command[0].split()
    cmd = "\"\'{0}\' \'{1}\'\"".format(cmd[0],cmd[1])
    pst.model_command.append(cmd)

    tpl_data = pst.model_input_data
    tpl_data["pest_file"] = tpl_data.pest_file.apply(lambda x: "\"{0}\"".format(x))
    tpl_data["model_file"] = tpl_data.model_file.apply(lambda x: "\"{0}\"".format(x))
    
    ins_data = pst.model_output_data
    ins_data["pest_file"] = ins_data.pest_file.apply(lambda x: "\"{0}\"".format(x))
    ins_data["model_file"] = ins_data.model_file.apply(lambda x: "\"{0}\"".format(x))
    
    pst.control_data.noptmax = 1
    pst.observation_data.loc[pst.nnz_obs_names,"weight"] = 1.0
    #pst.pestpp_options["panther_agent_freeze_on_fail"] = True
    pst.pestpp_options["ies_num_reals"] = 5
    pst.write(os.path.join(new_d, "pest.pst"))

    
    m_d = os.path.join(model_d,"master_pestpp")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    shutil.copytree(new_d,m_d)

    worker_root = os.path.join(model_d + "_\u0187")
    if os.path.exists(worker_root):
        shutil.rmtree(worker_root)
    os.makedirs(worker_root)    
    try:
        pyemu.os_utils.start_workers(new_d, exe_path, "pest.pst", 1, master_dir=m_d,
                               worker_root=worker_root,port=port,verbose=True)
    except Exception as e:
        if plat != "windows":
            raise Exception(e)
    else:
        if plat == "windows":
            raise Exception("should have failed")


def basic_test(model_d="ies_10par_xsec"):
    pyemu.Ensemble.reseed()
    base_d = os.path.join(model_d, "template")
    new_d = os.path.join(model_d, "test_template")
    if os.path.exists(new_d):
        shutil.rmtree(new_d)
    shutil.copytree(base_d, new_d)
    print(platform.platform().lower())
    pst = pyemu.Pst(os.path.join(new_d, "pest.pst"))
    cmd = pst.model_command[0].split()
    print(cmd)
    cmd = "\"\"{0}\" \"{1}\"\"".format(cmd[0],cmd[1])
    print(cmd)
    pst.model_command.append(cmd)
    cmd = pst.model_command[0].split()
    cmd = "\"\'{0}\' \'{1}\'\"".format(cmd[0],cmd[1])
    pst.model_command.append(cmd)

    tpl_data = pst.model_input_data
    tpl_data["pest_file"] = tpl_data.pest_file.apply(lambda x: "\"{0}\"".format(x))
    tpl_data["model_file"] = tpl_data.model_file.apply(lambda x: "\"{0}\"".format(x))
    
    ins_data = pst.model_output_data
    ins_data["pest_file"] = ins_data.pest_file.apply(lambda x: "\"{0}\"".format(x))
    ins_data["model_file"] = ins_data.model_file.apply(lambda x: "\"{0}\"".format(x))
    
    pst.control_data.noptmax = 0
    pst.write(os.path.join(new_d, "pest.pst"))
    pyemu.os_utils.run("{0} pest.pst".format(exe_path),cwd=new_d)
    # pst.write(os.path.join(new_d, "pest.pst"),version=2)
    # pyemu.os_utils.run("{0} pest.pst".format(exe_path),cwd=new_d)
    
    tpl_data = pst.model_input_data
    tpl_data["pest_file"] = tpl_data.pest_file.apply(lambda x: "\'{0}\'".format(x))
    tpl_data["model_file"] = tpl_data.model_file.apply(lambda x: "\'{0}\'".format(x))
    
    ins_data = pst.model_output_data
    ins_data["pest_file"] = ins_data.pest_file.apply(lambda x: "\'{0}\'".format(x))
    ins_data["model_file"] = ins_data.model_file.apply(lambda x: "\'{0}\'".format(x))

    pst.control_data.noptmax = 0
    pst.write(os.path.join(new_d, "pest.pst"))
    pyemu.os_utils.run("{0} pest.pst".format(exe_path),cwd=new_d)
    # pst.write(os.path.join(new_d, "pest.pst"),version=2)
    # pyemu.os_utils.run("{0} pest.pst".format(exe_path),cwd=new_d)

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
                           worker_root=model_d,port=port,verbose=True)
    #pyemu.os_utils.run("{0} pest.pst".format(exe_path.replace("-ies","-sen")),cwd=new_d)
    df = pd.read_csv(os.path.join(m_d, "pest.mio"),index_col=0)

    # run sweep
    m_d = os.path.join(model_d,"master_sweep1")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pyemu.os_utils.start_workers(new_d, exe_path.replace("-ies","-swp"), "pest.pst", 5, master_dir=m_d,
                           worker_root=model_d,port=port,verbose=True)
    df = pd.read_csv(os.path.join(m_d, "sweep_out.csv"),index_col=0)
    
    m_d = os.path.join(model_d,"master_pestpp-glm")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pyemu.os_utils.start_workers(new_d, exe_path.replace("-ies","-glm"), "pest.pst", 10, master_dir=m_d,
                           worker_root=model_d,port=port,verbose=True)
    #pyemu.os_utils.run("{0} pest.pst".format(exe_path.replace("-ies","-glm")),cwd=new_d)
    df = pd.read_csv(os.path.join(m_d, "pest.par.usum.csv"),index_col=0)

    m_d = os.path.join(model_d,"master_pestpp-ies")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pyemu.os_utils.start_workers(new_d, exe_path, "pest.pst", 10, master_dir=m_d,
                           worker_root=model_d,port=port,verbose=True)



def glm_save_binary_test():
    model_d = "ies_10par_xsec"

    t_d = os.path.join(model_d, "template")
    m_d = os.path.join(model_d, "master_save_binary")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pst = pyemu.Pst(os.path.join(t_d, "pest.pst"))
    pst.pestpp_options = {"glm_num_reals":30,"save_binary":True}
    pst.control_data.noptmax = 1
    pst.write(os.path.join(t_d, "pest_save_binary.pst"))
    pyemu.os_utils.start_workers(t_d, exe_path.replace("-ies", "-glm"), "pest_save_binary.pst", 10, master_dir=m_d,
                                worker_root=model_d, port=port)

    pe = pyemu.ParameterEnsemble.from_binary(pst=pst,filename=os.path.join(m_d,"pest_save_binary.post.paren.jcb"))
    pe = pyemu.ObservationEnsemble.from_binary(pst=pst,filename=os.path.join(m_d, "pest_save_binary.post.obsen.jcb"))

def sweep_forgive_test():
    model_d = "ies_10par_xsec"
    t_d = os.path.join(model_d,"template")
    m_d = os.path.join(model_d,"master_sweep_forgive")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pst = pyemu.Pst(os.path.join(t_d,"pest.pst"))
    pe = pyemu.ParameterEnsemble.from_uniform_draw(pst,num_reals=50)#.loc[:,pst.par_names[:2]]
    pe.loc[:,pst.par_names[2:]] = pst.parameter_data.loc[pst.par_names[2:],"parval1"].values
    pe.to_csv(os.path.join(t_d,"sweep_in.csv"))
    pst.pestpp_options["sweep_forgive"] = True
    pst.control_data.noptmax = -1
    pst.pestpp_options["ies_num_reals"] = 5
    pst.write(os.path.join(t_d,"pest_forgive.pst"))


    pyemu.os_utils.start_workers(t_d, exe_path.replace("-ies","-swp"), "pest_forgive.pst", 10, master_dir=m_d,
                           worker_root=model_d,port=port)
    df1 = pd.read_csv(os.path.join(m_d, "sweep_out.csv"),index_col=0)

    pe = pe.loc[:,pst.par_names[:2]]
    pe.to_csv(os.path.join(t_d,"sweep_in.csv"))
    pst.pestpp_options["sweep_forgive"] = True
    pst.write(os.path.join(t_d,"pest_forgive.pst"))
    pyemu.os_utils.start_workers(t_d, exe_path.replace("-ies","-swp"), "pest_forgive.pst", 10, master_dir=m_d,
                           worker_root=model_d,port=port)
    df2 = pd.read_csv(os.path.join(m_d, "sweep_out.csv"),index_col=0)
    diff = df1 - df2
    print(diff.max())
    assert diff.max().max() == 0.0




def inv_regul_test():
    model_d = "ies_10par_xsec"
    
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
                           worker_root=model_d,port=port)
    

def tie_by_group_test():
    model_d = "ies_10par_xsec"
    
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
    m_d = os.path.join(model_d,"master_tie_by_group_sen")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pyemu.os_utils.start_workers(t_d, exe_path.replace("-ies","-sen"), "pest_tied.pst", 5, master_dir=m_d,
                           worker_root=model_d,port=port)
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
    m_d = os.path.join(model_d,"master_tie_by_group_glm")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pyemu.os_utils.start_workers(t_d, exe_path.replace("-ies","-glm"), "pest_tied.pst", 5, master_dir=m_d,
                           worker_root=model_d,port=port)
    jco = pyemu.Jco.from_binary(os.path.join(m_d,"pest_tied.jcb"))
    assert jco.shape[1] == 2,jco.shape
    par_df = pyemu.pst_utils.read_parfile(os.path.join(m_d,"pest_tied.par"))
    print(par_df)
    print(tied_names)
    too_low = par.loc[par_df.parval1 < par.parlbnd,"parnme"]
    assert too_low.shape[0] == 0,too_low
    too_high = par.loc[par_df.parval1 > par.parubnd, "parnme"]
    assert too_high.shape[0] == 0, too_high

    
    m_d = os.path.join(model_d,"master_tie_by_group_ies")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pst.control_data.noptmax = 1
    pst.write(os.path.join(t_d, "pest_tied.pst"))

    pyemu.os_utils.start_workers(t_d, exe_path, "pest_tied.pst", 10, master_dir=m_d,
                           worker_root=model_d,port=port)
    
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
                                worker_root=model_d, port=port)
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
    m_d = os.path.join(model_d,"master_tie_by_group_swp")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pyemu.os_utils.start_workers(t_d, exe_path.replace("-ies","-swp"), "pest_tied.pst", 5, master_dir=m_d,
                           worker_root=model_d,port=port)
    pst.control_data.noptmax = 3
    pst.pestpp_options["enforce_tied_bounds"] = False
    pst.write(os.path.join(t_d, "pest_tied.pst"))
    pyemu.os_utils.start_workers(t_d, exe_path.replace("-ies","-glm"), "pest_tied.pst", 5, master_dir=m_d,
                           worker_root=model_d,port=port)
    jco = pyemu.Jco.from_binary(os.path.join(m_d,"pest_tied.jcb"))
    assert jco.shape[1] == 2,jco.shape



def unc_file_test():
    model_d = "ies_10par_xsec"
    
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

    cov.to_uncfile(os.path.join(m_d, "pest.unc"), covmat_file=None)
    pst.control_data.noptmax = -2
    pst.pestpp_options["ies_num_reals"] = 100000
    pst.pestpp_options["ies_enforce_bounds"] = False
    pst.write(os.path.join(m_d, "pest_unc.pst"))
    pyemu.os_utils.run("{0} {1}".format(exe_path, "pest_unc.pst"), cwd=m_d)
    pe_3 = pd.read_csv(os.path.join(m_d, "pest_unc.0.par.csv"), index_col=0).apply(np.log10)
    print(pe_3.std(ddof=0))
    pe_std = pe_3.std(ddof=0)
    for r,v  in zip(cov.row_names,cov.x):
        d = np.abs(pe_std.loc[r] - np.sqrt(v))

        print(r,v,np.sqrt(v),d)
        assert d < 0.01
    pst.control_data.noptmax = -1
    pst.write(os.path.join(m_d, "pest_unc.pst"))
    pyemu.os_utils.run("{0} {1}".format(exe_path.replace("-ies","-glm"), "pest_unc.pst"), cwd=m_d)
    fosm_df = pd.read_csv(os.path.join(m_d,"pest_unc.par.usum.csv"),index_col=0)
    cov_df = cov.to_dataframe()
    for pname,prior_std in zip(fosm_df.index,fosm_df.prior_stdev):
        d = np.abs(prior_std - np.sqrt(cov_df.loc[pname,pname]))
        print(pname,d)
        assert d < 1.0e-4


def parchglim_test():
    model_d = "ies_10par_xsec"
    
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
    print(p_df.loc["stage","parval1"],fpm)
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
                           worker_root=model_d,port=port)

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
    model_d = "sen_invest"
    t_d = os.path.join(model_d, "template")
    if os.path.exists(t_d):
        shutil.rmtree(t_d)
    os.makedirs(t_d)
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
        f.write(r"df = pd.read_csv('in.dat',index_col=0,sep='\s+',names=['name','value'])")
        f.write("\n")
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
                                 worker_root=model_d, port=port)
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
    pst.parameter_data.loc["k_10","parval1"] = 12345
    pst.parameter_data.loc["k_10","parubnd"] = 200000
    pst.pestpp_options["tpl_force_decimal"] = True
    pst.control_data.noptmax = 0
    pst.write(os.path.join(t_d,"pest.pst"))
    pyemu.os_utils.run("{0} pest.pst".format(exe_path),cwd=t_d)
    with open(os.path.join(t_d,"hk_Layer_1.ref"),'r') as f:
        for line in f:
            if "e" in line.lower():
                raise Exception(line)
    pst.pestpp_options.pop("tpl_force_decimal")
    pst.control_data.noptmax = -1
    pst.parameter_data.loc["k_10","parval1"] = 120
    pst.parameter_data.loc["k_10","parubnd"] = 200
    pst.write(os.path.join(t_d,"pest.pst"))
    

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

    pst.control_data.noptmax = 0
    pst.write(os.path.join(t_d,"pest.pst"))
    pyemu.os_utils.run("{0} pest.pst".format(exe_path.replace("-ies","-glm")),cwd=t_d)

    # check the input file - the last two number should be the same
    arr = np.loadtxt(os.path.join(t_d,"hk_Layer_1.ref"))
    assert arr[-2] == arr[-1],arr[-2] - arr[-1]

    lines_tpl = open(os.path.join(t_d,"hk_Layer_1.ref.tpl"),'r').readlines()
    lines_in = open(os.path.join(t_d,"hk_Layer_1.ref"),'r').readlines()
    assert len(lines_tpl) - 1 == len(lines_in)

    pst = pyemu.Pst(os.path.join(t_d,"pest.pst"))
    pst.control_data.noptmax = -1
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

    pst.control_data.noptmax = 0
    pst.write(os.path.join(t_d,"pest.pst"))
    pyemu.os_utils.run("{0} pest.pst".format(exe_path.replace("-ies","-glm")),cwd=t_d)

    # check the input file - the last two number should be the same
    arr = np.loadtxt(os.path.join(t_d,"hk_Layer_1.ref"))
    assert arr[-2] == arr[-1]

    lines_tpl = open(os.path.join(t_d,"hk_Layer_1.ref.tpl"),'r').readlines()
    lines_in = open(os.path.join(t_d,"hk_Layer_1.ref"),'r').readlines()
    assert len(lines_tpl) - 1 == len(lines_in)

    pst = pyemu.Pst(os.path.join(os.path.join(model_d,"template"), "pest.pst"))
    #pst.drop_observations(os.path.join(t_d,"AOC_obs.txt.ins"),pst_path=".")
    dum_obs = ['h01_03', 'h01_07', 'dummy_obs']
    pst.observation_data.drop(index=dum_obs, inplace=True)
    pst.model_output_data = pd.DataFrame({"pest_file":"out1dum.dat.ins",
                                          "model_file":'out1.dat'},index=["out1dum.dat.ins"])
    #pst.instruction_files = ['out1dum.dat.ins']
    pst.control_data.noptmax = -1
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
    pst.pestpp_options["ies_no_noise"] = False
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

    t_d = os.path.join(model_d,"template")
    m_d = os.path.join(model_d,"master_ies_glm_loc")
    #if os.path.exists(m_d):
    #    shutil.rmtree(m_d)
    pst = pyemu.Pst(os.path.join(t_d,"freyberg6_run_ies.pst"))
    pst.control_data.noptmax = 0
    pst.write(os.path.join(t_d,"freyberg6_run_ies.pst"))
    pyemu.os_utils.run("{0} freyberg6_run_ies.pst".format(exe_path),cwd=t_d)

    pst.control_data.noptmax = 3
    par = pst.parameter_data

    eff_lb = (par.parlbnd + (np.abs(par.parlbnd.values)*.01)).to_dict()
    eff_ub = (par.parubnd - (np.abs(par.parlbnd.values)*.01)).to_dict()
    log_idx = par.partrans.apply(lambda x: x=="log").to_dict()
    for p,log in log_idx.items():
        if log:
            lb = np.log10(par.loc[p,"parlbnd"])
            eff_lb[p] = (lb + (np.abs(lb)*.01))
            ub = np.log10(par.loc[p,"parubnd"])
            eff_ub[p] = (ub - (np.abs(ub)*.01))

    pargp_map = par.groupby(par.pargp).groups
    print(pargp_map)

    


    m_d = os.path.join(model_d, "master_ies_glm_noloc_standard")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pst = pyemu.Pst(os.path.join(t_d, "freyberg6_run_ies.pst"))
    pst.pestpp_options.pop("ies_localizer",None)
    pst.pestpp_options.pop("ies_autoadaloc",None)
    pst.pestpp_options["ies_bad_phi_sigma"] = 2.5
    pst.pestpp_options["ies_num_reals"] = 30
    pst.pestpp_options["ensemble_output_precision"] = 40
    pst.control_data.noptmax = 3
    pst_name = "freyberg6_run_ies_glm_noloc_standard.pst"
    pst.write(os.path.join(t_d, pst_name))
    pyemu.os_utils.start_workers(t_d, exe_path, pst_name, num_workers=15,
                                 master_dir=m_d, worker_root=model_d, port=port)
    
    


    phidf = pd.read_csv(os.path.join(m_d,pst_name.replace(".pst",".phi.actual.csv")))
    assert phidf.shape[0] == pst.control_data.noptmax + 1
    for i in range(1,pst.control_data.noptmax+1):
        pcs = pd.read_csv(os.path.join(m_d,pst_name.replace(".pst",".{0}.pcs.csv".format(i))),index_col=0)
        #print(pcs)
        pe = pd.read_csv(os.path.join(m_d,pst_name.replace(".pst",".{0}.par.csv".format(i))),index_col=0)
        print(pe.shape)
        #print(pe)
        groups = pcs.index.values.copy()
        groups.sort()
        for group in groups:
            pnames = pargp_map[group].values
            lb_count,ub_count = 0,0
            for pname in pnames:
                lb,ub = eff_lb[pname],eff_ub[pname]
                v = pe.loc[:,pname].values.copy()
                if log_idx[pname]:
                    v = np.log10(v)
                low = np.zeros_like(v,dtype=int)
                low[v < lb] = 1
                high = np.zeros_like(v,dtype=int)
                high[v > ub] = 1
                lb_count += low.sum()
                ub_count += high.sum()
            print(i,group,len(pnames),lb_count,pcs.loc[group,"num_at_near_lbound"],ub_count,pcs.loc[group,"num_at_near_ubound"])
            assert lb_count == pcs.loc[group,"num_at_near_lbound"]
            assert ub_count == pcs.loc[group,"num_at_near_ubound"]
    
    pst.pestpp_options["ies_run_realname"] = "base"
    pst.pestpp_options["ies_par_en"] = "{0}.{1}.par.csv".format(pst_name.replace(".pst",""),3)
    pst.control_data.noptmax = -2
    pst.write(os.path.join(m_d, "test.pst"))
    pyemu.os_utils.run("{0} test.pst".format(exe_path),cwd=m_d)

    pst.pestpp_options.pop("ies_run_realname")
    pst.pestpp_options.pop("ies_par_en")

    pst.write(os.path.join(t_d,"freyberg6_run_ies_glm_loc.pst"))

    m_d = os.path.join(model_d, "master_ies_glm_covloc")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pst.pestpp_options["ies_loc_type"] = "cov"
    pst.pestpp_options["ies_lambda_mults"] = [1.0]
    pst.pestpp_options["lambda_scale_fac"] = [1.0]
    pst.control_data.noptmax = 2
    #pst.pestpp_options.pop("ies_localizer",None)
    pst.write(os.path.join(t_d, "freyberg6_run_ies_glm_covloc.pst"))
    pyemu.os_utils.start_workers(t_d, exe_path, "freyberg6_run_ies_glm_covloc.pst", num_workers=15,
                                 master_dir=m_d, worker_root=model_d, port=port)

    m_d = os.path.join(model_d, "master_ies_glm_noloc")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pst = pyemu.Pst(os.path.join(t_d, "freyberg6_run_ies.pst"))
    pst.pestpp_options.pop("ies_localizer",None)
    pst.pestpp_options.pop("ies_autoadaloc",None)
    pst.pestpp_options["ies_lambda_mults"] = [1.0]
    pst.pestpp_options["lambda_scale_fac"] = [1.0]
    pst.control_data.noptmax = 2
    pst.write(os.path.join(t_d, "freyberg6_run_ies_glm_noloc.pst"))
    pyemu.os_utils.start_workers(t_d, exe_path, "freyberg6_run_ies_glm_noloc.pst", num_workers=15,
                                 master_dir=m_d, worker_root=model_d, port=port)

    m_d = os.path.join(model_d, "master_ies_mda_loc")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pst = pyemu.Pst(os.path.join(t_d, "freyberg6_run_ies.pst"))
    pst.pestpp_options["ies_lambda_mults"] = [1.0]
    pst.pestpp_options["lambda_scale_fac"] = [1.0]
    pst.control_data.noptmax = 2
    pst.pestpp_options["ies_use_mda"] = True
    pst.write(os.path.join(t_d, "freyberg6_run_ies_mda_loc.pst"))
    pyemu.os_utils.start_workers(t_d, exe_path, "freyberg6_run_ies_mda_loc.pst", num_workers=15,
                                 master_dir=m_d, worker_root=model_d, port=port)
    
    # m_d = os.path.join(model_d, "master_ies_mda_covloc")
    # if os.path.exists(m_d):
    #     shutil.rmtree(m_d)
    # pst = pyemu.Pst(os.path.join(t_d, "freyberg6_run_ies.pst"))
    # pst.control_data.noptmax = 3
    # pst.pestpp_options["ies_use_mda"] = True
    # pst.pestpp_options["ies_loc_type"] = "cov"
    # pst.write(os.path.join(t_d, "freyberg6_run_ies_mda_covloc.pst"))
    # pyemu.os_utils.start_workers(t_d, exe_path, "freyberg6_run_ies_mda_covloc.pst", num_workers=15,
    #                              master_dir=m_d, worker_root=model_d, port=port)
    
    m_d = os.path.join(model_d, "master_ies_mda_noloc")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pst = pyemu.Pst(os.path.join(t_d, "freyberg6_run_ies.pst"))
    pst.control_data.noptmax = 2
    pst.pestpp_options["ies_lambda_mults"] = [1.0]
    pst.pestpp_options["lambda_scale_fac"] = [1.0]
    pst.pestpp_options["ies_use_mda"] = True
    pst.pestpp_options.pop("ies_localizer", None)
    pst.pestpp_options.pop("ies_autoadaloc", None)
    pst.write(os.path.join(t_d, "freyberg6_run_ies_mda_noloc.pst"))
    pyemu.os_utils.start_workers(t_d, exe_path, "freyberg6_run_ies_mda_noloc.pst", num_workers=15,
                                 master_dir=m_d, worker_root=model_d, port=port)

    m_d = os.path.join(model_d, "master_ies_glm_loc_mm")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pst = pyemu.Pst(os.path.join(t_d, "freyberg6_run_ies.pst"))
    pst.control_data.noptmax = 2
    pst.pestpp_options["ies_lambda_mults"] = [1.0]
    pst.pestpp_options["lambda_scale_fac"] = [1.0]
    pst.pestpp_options["ies_num_threads"] = 1
    pst.pestpp_options["ies_use_mda"] = False
    pst.pestpp_options.pop("ies_localizer", None)
    pst.pestpp_options.pop("ies_autoadaloc", None)
    pst.pestpp_options["ies_multimodal_alpha"] = 0.1
    pst.write(os.path.join(t_d, "freyberg6_run_ies_glm_loc_mm.pst"))
    pyemu.os_utils.start_workers(t_d, exe_path, "freyberg6_run_ies_glm_loc_mm.pst", num_workers=15,
                                 master_dir=m_d, worker_root=model_d, port=port)

    m_d = os.path.join(model_d, "master_ies_glm_noloc_mm")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pst = pyemu.Pst(os.path.join(t_d, "freyberg6_run_ies.pst"))
    pst.control_data.noptmax = 2
    pst.pestpp_options["ies_use_mda"] = False
    pst.pestpp_options["ies_lambda_mults"] = [1.0]
    pst.pestpp_options["lambda_scale_fac"] = [1.0]
    pst.pestpp_options.pop("ies_localizer", None)
    pst.pestpp_options.pop("ies_autoadaloc", None)
    pst.pestpp_options["ies_multimodal_alpha"] = 0.25
    pst.write(os.path.join(t_d, "freyberg6_run_ies_glm_noloc_mm.pst"))
    pyemu.os_utils.start_workers(t_d, exe_path, "freyberg6_run_ies_glm_noloc_mm.pst", num_workers=15,
                                 master_dir=m_d, worker_root=model_d, port=port)




def mf6_v5_sen_test():

    model_d = "mf6_freyberg"

    t_d = os.path.join(model_d,"template")
    m_d = os.path.join(model_d,"master_sen")
    if os.path.exists(m_d):
       shutil.rmtree(m_d)
    pst = pyemu.Pst(os.path.join(t_d,"freyberg6_run_sen.pst"))
    pst.pestpp_options["gsa_morris_p"] = 4
    pst.pestpp_options["gsa_morris_r"] = 4
    pst.pestpp_options["panther_transfer_on_finish"] = ["freyberg6_freyberg.cbc","freyberg6.lst","ies_prior.jcb"]
    pst.write(os.path.join(t_d,"freyberg6_run_sen_trn.pst"))
    m_d = os.path.join(model_d,"master_sen")
    pyemu.os_utils.start_workers(t_d, exe_path.replace("-ies","-sen"), "freyberg6_run_sen_trn.pst",
                                 num_workers=50, worker_root=model_d,
                                 port=4004,verbose=True,master_dir=m_d)

    pst = pyemu.Pst(os.path.join(m_d,"freyberg6_run_sen_trn.pst"))
    mio_file = os.path.join(m_d,"freyberg6_run_sen_trn.mio")
    assert os.path.exists(mio_file),mio_file
    df = pd.read_csv(mio_file)
    assert df.shape[0] > 1
    msn_file = mio_file.replace(".mio",".msn")
    assert os.path.exists(msn_file),msn_file
    msngrp_file = msn_file.replace(".msn",".group.msn")
    assert os.path.exists(msngrp_file),msngrp_file

    jcb_files = [f for f in os.listdir(m_d) if f.lower().startswith("ftx_") and f.lower().endswith(".jcb")]
    print(len(jcb_files))
    assert len(jcb_files) == 52
    for jcb_file in jcb_files:
        j = pyemu.Jco.from_binary(os.path.join(m_d,jcb_file))

    lst_files = [f for f in os.listdir(m_d) if f.lower().startswith("ftx_") and f.lower().endswith(".lst")]
    print(len(lst_files))
    assert len(lst_files) == 52
    
    cbc_files = [f for f in os.listdir(m_d) if f.lower().startswith("ftx") and f.lower().endswith(".cbc")]
    print(len(cbc_files))
    assert len(cbc_files) == 52
    



def mf6_v5_opt_stack_test():
    model_d = "mf6_freyberg"
    
    t_d = os.path.join(model_d,"template")
    m_d = os.path.join(model_d,"master_opt_stack")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pst = pyemu.Pst(os.path.join(t_d,"freyberg6_run_opt.pst"))
    m_d = os.path.join(model_d,"master_opt_stack")
    pyemu.os_utils.start_workers(t_d, exe_path.replace("-ies","-opt"), "freyberg6_run_opt.pst", 
                                 num_workers=15, master_dir=m_d,worker_root=model_d,
                                 port=port)

    assert os.path.exists(os.path.join(m_d,"freyberg6_run_opt.1.sim+chance.rei"))
    assert os.path.exists(os.path.join(m_d,"freyberg6_run_opt.1.obs_stack.csv"))


def mf6_v5_glm_test():
    model_d = "mf6_freyberg"
    
    
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

    
    t_d = os.path.join(model_d,"template")
    pst_name = "freyberg6_run_glm.pst"
    pst = pyemu.Pst(os.path.join(t_d,"freyberg6_run_glm.pst"))
    pst.pestpp_options = {}
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
    
def fr_fail_test():
    model_d = "ies_10par_xsec"
    base_d = os.path.join(model_d, "template")
    new_d = os.path.join(model_d, "test_template")
    if os.path.exists(new_d):
        shutil.rmtree(new_d)
    shutil.copytree(base_d, new_d)
    print(platform.platform().lower())
    pst = pyemu.Pst(os.path.join(new_d, "pest.pst"))
    with open(os.path.join(new_d,"run.py"),'w') as f:
        f.write("import pyemu\npyemu.os_utils.run('mfnwt 10par_xsec.nam')\nprint(junk)\n")
    pst.model_command = "python run.py"
    oe_file = os.path.join(new_d, "pest.0.obs.csv")
    if os.path.exists(oe_file):
        os.remove(oe_file)
    pst.control_data.noptmax = 1
    pst.pestpp_options["panther_transfer_on_fail"] = "10par_xsec.list"
    pst.pestpp_options["ies_num_reals"] = 10
    #pst.pestpp_options["panther_agent_freeze_on_fail"] = True
    pst.write(os.path.join(new_d, "pest.pst"))
    try:
        pyemu.os_utils.run("{0} pest.pst".format(exe_path),cwd=new_d)
    except:
        pass
    else:
        raise Exception("should have failed")

    assert not os.path.exists(oe_file)
    m_d = os.path.join(model_d,"fr_fail_master")
    try:
        pyemu.os_utils.start_workers(new_d,exe_path,"pest.pst",num_workers=5,worker_root=model_d,master_dir=m_d)
    except:
        pass
    else:
        raise Exception("should have failed")
    oe_file = os.path.join(m_d, "pest.0.obs.csv")
    assert not os.path.exists(oe_file)

    trx_files = [f for f in os.listdir(m_d) if f.endswith(".list")]
    print(trx_files)
    assert len(trx_files) == 11,len(trx_files)





def sen_grp_test():
    
    model_d = "ies_10par_xsec"

    t_d = os.path.join(model_d, "template")
    m_d = os.path.join(model_d, "master_sen_group")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    
    pyemu.os_utils.start_workers(t_d, exe_path.replace("-ies", "-sen"), "pest.pst", 10, master_dir=m_d,
                                worker_root=model_d, port=port)

    msn_file = os.path.join(m_d,"pest.msn")
    msndf = pd.read_csv(msn_file)

    grp_file = msn_file.replace(".msn",".group.msn")
    grpdf = pd.read_csv(grp_file)
    assert msndf.shape[0] == grpdf.shape[0]
    for col in ["sen_mean","sen_mean_abs","sen_std_dev"]:
        diff = np.abs(msndf.loc[:,col].sum() - grpdf.loc[:,col].sum())
        print(col,diff)
        assert diff < 1.0e-6


def agnostic_path_test():
    model_d = "ies_10par_xsec"

    t_d = os.path.join(model_d, "template")
    m_d = os.path.join(model_d, "test_path")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    shutil.copytree(t_d,m_d)

    pst = pyemu.Pst(os.path.join(m_d, "pest.pst"))
    pst.parameter_data.loc[pst.adj_par_names,"parval1"] = np.random.random(pst.npar_adj)
    pst.control_data.noptmax = 0
    pst.write(os.path.join(m_d,"pest.pst"))
    pyemu.os_utils.run("{0} pest.pst".format(exe_path),cwd=m_d)

    arr1 = np.loadtxt(os.path.join(m_d,"hk_Layer_1.ref"))
    arr2 = np.loadtxt(os.path.join(m_d,"nested","really","deep","hk_Layer_1.ref"))
    d = np.abs(arr1-arr2).sum()
    print(d)
    assert d == 0,d

def fr_timeout_test():
    model_d = "ies_10par_xsec"
    base_d = os.path.join(model_d, "template")
    new_d = os.path.join(model_d, "test_template")
    if os.path.exists(new_d):
        shutil.rmtree(new_d)
    shutil.copytree(base_d, new_d)
    print(platform.platform().lower())
    pst = pyemu.Pst(os.path.join(new_d, "pest.pst"))
    with open(os.path.join(new_d,"run.py"),'w') as f:
        f.write("import os\nimport time\nimport pyemu\npyemu.os_utils.run('mfnwt 10par_xsec.nam')\n")
        f.write("if not os.path.exists('run.info'):\n    exit()\n")
        f.write("lines = open('run.info','r').readlines()\nrnum = int(lines[-1].split()[-1].split(':')[-1])\n")
        f.write("if rnum % 2 == 0:\n    time.sleep(10000000)\n")
    pst.model_command = "python run.py"
    oe_file = os.path.join(new_d, "pest.0.obs.csv")
    if os.path.exists(oe_file):
        os.remove(oe_file)
    pst.control_data.noptmax = -1
    pst.pestpp_options["overdue_giveup_fac"] = 1.0e+10
    pst.pestpp_options["overdue_giveup_minutes"] = 0.25
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options["ies_include_base"] = False
    pst.pestpp_options["max_run_fail"] = 1

    #pst.pestpp_options["panther_transfer_on_fail"] = "10par_xsec.list"
    pst.pestpp_options["panther_agent_freeze_on_fail"] = False
    pst.write(os.path.join(new_d, "pest.pst"))

    m_d = os.path.join(model_d,"fr_timeout_master")
    pyemu.os_utils.start_workers(new_d,exe_path,"pest.pst",num_workers=5,worker_root=model_d,master_dir=m_d)
    oe_file = os.path.join(m_d, "pest.0.obs.csv")
    assert os.path.exists(oe_file)
    oe = pd.read_csv(oe_file,index_col=0)
    print(oe.shape)
    assert oe.shape[0] == 5,oe.shape

    with open(os.path.join(new_d,"run.py"),'w') as f:
        f.write("import os\nimport time\nimport pyemu\npyemu.os_utils.run('mfnwt 10par_xsec.nam')\n")
        f.write("if not os.path.exists('run.info'):\n    exit()\n")
        f.write("lines = open('run.info','r').readlines()\nrnum = int(lines[-1].split()[-1].split(':')[-1])\n")
        f.write("if rnum % 10 == 0:\n    print(junk)\n")
    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 20 # hard coded to conditional below
    pst.pestpp_options["panther_agent_freeze_on_fail"] = True
    #pst.pestpp_options["overdue_giveup_fac"] = 1.0e+10
    #pst.pestpp_options["overdue_giveup_minutes"] = 0.25
    pst.write(os.path.join(new_d, "pest.pst"))
    pst.control_data.noptmax = 2

    pst.write(os.path.join(new_d, "pest.pst"))
    m_d = os.path.join(model_d,"fr_timeout_master_freeze")
    #num workers hard coded with conditional below
    pyemu.os_utils.start_workers(new_d,exe_path,"pest.pst",num_workers=10,worker_root=model_d,master_dir=m_d)
    #df = pyemu.helpers.parse_rmr_file(os.path.join(m_d,"pest.rmr"))
    #print(df.action.to_list())
    oe = pd.read_csv(os.path.join(m_d,"pest.{0}.obs.csv".format(pst.control_data.noptmax)),index_col=0)
    assert oe.shape[0] == 17 # hard coded to num reals
    with open(os.path.join(m_d,"pest.rmr"),'r') as f:
        for line in f:
            if "timeout" in line.lower():
                raise Exception()
            if line.strip().lower().endswith("agents connected"):
                num = int(line.strip().split()[0])
                print(line.strip())
                assert num == 7 # hard coded above



def ins_missing_e_test():
    import os
    import shutil
    import pyemu
    t_d = os.path.join("tplins_test_1","test_missing_e")
    if os.path.exists(t_d):
        shutil.rmtree(t_d)
    os.makedirs(t_d)
    bd = os.getcwd()
    os.chdir(t_d)
    with open("model.output.bak",'w') as f:
        f.write("12345-123\n")
    pst = pyemu.helpers.pst_from_parnames_obsnames(["p1"],["o1"])
    pst.control_data.noptmax = 0
    with open("forward_run.py",'w') as f:
        f.write("import shutil\n")
        f.write("shutil.copy2('model.output.bak','model.output')\n")
    pst.model_command = "python forward_run.py"
    pst.write("test.pst")
    os.chdir(bd)
    try:
        pyemu.os_utils.run("{0} test.pst".format(exe_path),cwd=t_d)
    except:
        pass
    else:
        raise Exception("should have failed")


def prep_ends():
    model_d = "mf6_freyberg"
    base_d = os.path.join(model_d, "template")
    new_d = os.path.join(model_d, "ends")
    if os.path.exists(new_d):
        shutil.rmtree(new_d)
    os.makedirs(new_d)
    skip = ["pst","csv","log","grb","hds","par","rei","lst","jcb","cov","rec","cbc"]
    files = [f for f in os.listdir(base_d) if f.lower().split('.')[-1] not in skip]
    print(files)
    [shutil.copy2(os.path.join(base_d,f),os.path.join(new_d,f)) for f in files]
    [shutil.copy2(os.path.join(base_d,f),os.path.join(new_d,f)) for f in ["ies_prior.jcb"]]
    
    pyemu.os_utils.run("mf6",cwd=new_d)
    pst = pyemu.Pst(os.path.join(base_d,"freyberg6_run_ies.pst"))
    pst.control_data.noptmax = 0
    pst.pestpp_options = {}
    pst.pestpp_options["ies_par_en"] = "prior.jcb"

    pst.write(os.path.join(new_d,"freyberg6_run_ies.pst"),version=2)
    pyemu.os_utils.run("pestpp-ies freyberg6_run_ies.pst",cwd=new_d)

    build_and_draw_prior(new_d,num_reals=5000)
    pst.control_data.noptmax = -1
    pst.write(os.path.join(new_d,"freyberg6_run_ies.pst"),version=2)
    m_d = os.path.join(model_d,"ends_master")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)

    pyemu.os_utils.start_workers(new_d,"pestpp-ies","freyberg6_run_ies.pst",num_workers=15,worker_root=model_d,master_dir=m_d)


def build_and_draw_prior(t_d="ends",num_reals=500):
    import flopy

    sim = flopy.mf6.MFSimulation.load(sim_ws=t_d)
    m = sim.get_model("freyberg6")
    xgrid = m.modelgrid.xcellcenters
    ygrid = m.modelgrid.ycellcenters
    pst = pyemu.Pst(os.path.join(t_d,"freyberg6_run_ies.pst"))
    par = pst.parameter_data
    static_par = par.loc[par.parnme.apply(lambda x: x[:3] in ["npf","sto"]),:].copy()
    static_par.loc[:, "i"] = static_par.parnme.apply(lambda x: int(x.split('_')[3]))
    static_par.loc[:, "j"] = static_par.parnme.apply(lambda x: int(x.split('_')[4]))
    static_par.loc[:, "x"] = static_par.apply(lambda x: xgrid[x.i,x.j],axis=1)
    static_par.loc[:, "y"] = static_par.apply(lambda x: ygrid[x.i, x.j], axis=1)
    static_par.loc[:,"pargp"] = static_par.parnme.apply(lambda x: "_".join(x.split('_')[:3]))

    wel_par = par.loc[par.parnme.apply(lambda x: x.startswith("wel")),:].copy()
    wel_par.loc[:,"x"] = wel_par.parnme.apply(lambda x: int(x.split('_')[-1]))
    wel_par.loc[:,"y"] = 0.0
    wel_par.loc[:,"pargp"] = wel_par.parnme.apply(lambda x: '_'.join(x.split('_')[:-1]))

    rch_par = par.loc[par.parnme.str.startswith("rch"),:].copy()
    rch_par.loc[:,"x"] = rch_par.parnme.apply(lambda x: int(x.split('_')[-1]))
    rch_par.loc[:,"y"] = 0.0

    spatial_v = pyemu.geostats.ExpVario(contribution=1.0,a=1000.0)
    temporal_v = pyemu.geostats.ExpVario(contribution=1.0,a=3)
    spatial_gs = pyemu.geostats.GeoStruct(variograms=spatial_v)
    temporal_gs = pyemu.geostats.GeoStruct(variograms=temporal_v)

    static_struct_dict = {spatial_gs:[]}
    sgrps = static_par.pargp.unique()
    sgrps.sort()
    for pargp in sgrps:
        static_struct_dict[spatial_gs].append(static_par.loc[static_par.pargp==pargp,["parnme","x","y","i","j"]])
    temporal_struct_dict = {temporal_gs: [rch_par.loc[:, ["parnme", "x", "y"]]]}
    wgrps = wel_par.pargp.unique()
    wgrps.sort()
    for pargp in wgrps:
        temporal_struct_dict[temporal_gs].append(wel_par.loc[wel_par.pargp == pargp, ["parnme", "x", "y"]])

    struct_dict = static_struct_dict
    for k,v in temporal_struct_dict.items():
        struct_dict[k] = v
    print(struct_dict)
    np.random.seed(pyemu.en.SEED)
    pe = pyemu.helpers.geostatistical_draws(pst,struct_dict=struct_dict,num_reals=num_reals)
    pe.to_binary(os.path.join(t_d,"prior.jcb"))




def run():
    model_d = "mf6_freyberg"
    t_d = os.path.join(model_d,"template")
    pst_name = "freyberg6_run_ies_glm_noloc_standard.pst"
    pyemu.os_utils.start_workers(t_d, exe_path, pst_name, num_workers=15,
                                 worker_root=model_d, port=4004)

def sweep_bin_test():

    model_d = "ies_10par_xsec"
    t_d = os.path.join(model_d,"template")
    m_d = os.path.join(model_d,"master_sweep_bin")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pst = pyemu.Pst(os.path.join(t_d,"pest.pst"))
    pe = pyemu.ParameterEnsemble.from_uniform_draw(pst,num_reals=50)#.loc[:,pst.par_names[:2]]

    pe.to_csv(os.path.join(t_d,"sweep_in.csv"))
    pe._df.index = pe.index.map(str)
    print(pe.index)
    pe.to_dense(os.path.join(t_d,"sweep_in.bin"))
    pst.pestpp_options["ies_par_en"] = "sweep_in.csv"
    pst.pestpp_options["sweep_forgive"] = True
    pst.pestpp_options["sweep_parameter_file"] = "sweep_in.bin"
    pst.control_data.noptmax = -1
    pst.pestpp_options.pop("ies_num_reals",None)
    pst.write(os.path.join(t_d,"pest_forgive.pst"))
    pst.pestpp_options["sweep_output_file"] = "sweep_out.bin"
    pst.pestpp_options["sweep_chunk"] = 9
    pst.pestpp_options["ies_include_base"] = False
    pst.write(os.path.join(t_d,"pest_forgive.pst"))
    m_d = os.path.join(model_d,"master_sweep_bin_base")
    pyemu.os_utils.start_workers(t_d, exe_path, "pest_forgive.pst", 10, master_dir=m_d,
                           worker_root=model_d,port=port)
    df1 = pd.read_csv(os.path.join(m_d, "pest_forgive.0.obs.csv"),index_col=0)
    assert df1.shape[0] == pe.shape[0]
    m_d = os.path.join(model_d, "master_sweep_bin")
    pyemu.os_utils.start_workers(t_d, exe_path.replace("-ies", "-swp"), "pest_forgive.pst", 10, master_dir=m_d,
                                 worker_root=model_d, port=port)
    df2 = pyemu.Matrix.from_binary(os.path.join(m_d,"sweep_out.bin")).to_dataframe()
    print(df2)
    print(df1)
    assert df2.shape == df1.shape
    diff = (df1.values - df2.values)
    print(diff)
    print(diff.max())
    print(np.abs(diff).max())
    assert np.abs(diff).max() < 1e-7



    pst.pestpp_options.pop("sweep_output_file")
    pst.pestpp_options["save_dense"] = True
    m_d = os.path.join(model_d,"master_sweep_bin_base2")
    pyemu.os_utils.start_workers(t_d, exe_path, "pest_forgive.pst", 10, master_dir=m_d,
                           worker_root=model_d,port=port)
    df1 = pd.read_csv(os.path.join(m_d, "pest_forgive.0.obs.csv"),index_col=0)
    assert df1.shape[0] == pe.shape[0]
    m_d = os.path.join(model_d, "master_sweep_bin2")
    pyemu.os_utils.start_workers(t_d, exe_path.replace("-ies", "-swp"), "pest_forgive.pst", 10, master_dir=m_d,
                                 worker_root=model_d, port=port)
    df2 = pyemu.Matrix.from_binary(os.path.join(m_d,"sweep_out.bin")).to_dataframe()
    print(df2)
    print(df1)
    assert df2.shape == df1.shape
    diff = (df1.values - df2.values)
    print(diff)
    print(diff.max())
    print(np.abs(diff).max())
    assert np.abs(diff).max() < 1e-7

#def fail_test():
#    raise Exception("fail please")


def tenpar_collapse_invest():
    model_d = "ies_10par_xsec"
    pyemu.Ensemble.reseed()
    base_d = os.path.join(model_d, "template")
    new_d = os.path.join(model_d, "test_template")
    if os.path.exists(new_d):
        shutil.rmtree(new_d)
    shutil.copytree(base_d, new_d)
    print(platform.platform().lower())
    pst = pyemu.Pst(os.path.join(new_d, "pest.pst"))
    print(pst.model_command)
    
    obs = pst.observation_data
    obs.loc[:,"weight"] = 0.0

    nzoname = obs.obsnme.str.startswith("h01").index[3]
    obs.loc[nzoname,"weight"] = 10.0
    
    obs.loc[nzoname,"standard_deviation"] = 0.1

    pst.parameter_data.loc[:,"partrans"] = "log"
    pst.pestpp_options["ies_multimodal_alpha"] = 0.99

    # set noptmax
    num_reals = [10,30,50,100,400,1000]
    pst.control_data.noptmax = 10
    # for num_real in num_reals:

    #     # wipe all pestpp options
    #     pst.pestpp_options = {}
    #     pst.pestpp_options["ies_num_reals"] = num_real
    #     pst.write(os.path.join(new_d, "pest.pst"),version=2)
               
    #     m_d = os.path.join(model_d,"master_ies_base_{0}reals".format(pst.pestpp_options["ies_num_reals"]))
    #     if os.path.exists(m_d):
    #         shutil.rmtree(m_d)
    #     num_workers = 50
    #     if num_real > 500:
    #         num_workers = 200
    #     pyemu.os_utils.start_workers(new_d, exe_path, "pest.pst", num_workers, master_dir=m_d,
    #                            worker_root=model_d,port=port,verbose=True)

    # pst.pestpp_options = {}
    # pst.pestpp_options["ies_num_reals"] = 100000
    # pst.control_data.noptmax = -1
    # pst.write(os.path.join(new_d, "pest.pst"))
           
    # m_d = os.path.join(model_d,"master_ies_base_{0}reals".format(pst.pestpp_options["ies_num_reals"]))
    # if os.path.exists(m_d):
    #     shutil.rmtree(m_d)
    # pyemu.os_utils.start_workers(new_d, exe_path, "pest.pst", 200, master_dir=m_d,
    #                        worker_root=model_d,port=port,verbose=True)


    #pst.observation_data.loc[nzoname,"obsval"] += 8
    pst.observation_data.loc[nzoname,"weight"] = 100.0
    pst.observation_data.loc[nzoname,"obsval"] -= 1.5

    # for num_real in num_reals:

    #     # wipe all pestpp options
    #     pst.pestpp_options = {}
    #     pst.pestpp_options["ies_num_reals"] = num_real
    #     pst.write(os.path.join(new_d, "pest.pst"),version=2)
               
    #     m_d = os.path.join(model_d,"master_ies_corrupt_{0}reals".format(pst.pestpp_options["ies_num_reals"]))
    #     if os.path.exists(m_d):
    #         shutil.rmtree(m_d)
    #     num_workers = 50
    #     if num_real > 500:
    #         num_workers = 200
    #     pyemu.os_utils.start_workers(new_d, exe_path, "pest.pst", num_workers, master_dir=m_d,
    #                            worker_root=model_d,port=port,verbose=True)

    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 100000
    pst.control_data.noptmax = -1
    pst.write(os.path.join(new_d, "pest.pst"),version=2)
           
    m_d = os.path.join(model_d,"master_ies_corrupt_{0}reals".format(pst.pestpp_options["ies_num_reals"]))
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pyemu.os_utils.start_workers(new_d, exe_path, "pest.pst", 200, master_dir=m_d,
                           worker_root=model_d,port=port,verbose=True)

def plot_collapse_invest():
    b_d = "ies_10par_xsec"
    m_ds = [os.path.join(b_d,d) for d in os.listdir(b_d) if os.path.isdir(os.path.join(b_d,d)) and d.startswith("master") and d.endswith("reals")]

    pes, oes, phidfs = {},{},{}
    names = []
    name_dict = {}
    min_phi = 1e300
    for m_d in m_ds:
        phidf = pd.read_csv(os.path.join(m_d,"pest.phi.actual.csv"))
        if "corrupt" in m_d:
            min_phi = min(min_phi,phidf.iloc[:,6:].values.min())
        pst = pyemu.Pst(os.path.join(m_d,"pest.pst"))
        p,o = [],[]
        for i in [0,phidf.iteration.max()]:
            oe = pd.read_csv(os.path.join(m_d,"pest.{0}.obs.csv".format(i)),index_col=0)
            oe.index = oe.index.map(str)
            pe = pd.read_csv(os.path.join(m_d,"pest.{0}.par.csv".format(i)),index_col=0)
            pe.index = pe.index.map(str)
            o.append(oe)
            p.append(pe)
        #if phidf.iteration.max() == 0:
            #realphis = phidf.iloc[0,6:]
            #realkeep = realphis.loc[realphis < pst.nnz_obs * 1.1]
            #print(realkeep.shape)
            #print(o[0].index)
            #print(p[0].index)
            
        #    p.append(p[0])#.loc[realkeep.index,:].copy())

        #    o.append(o[0])#.loc[realkeep.index,:].copy())
            
        if phidf.iteration.max() == 0:
            nreals = oe.shape[0]
            name = "{0}, {1} realizations".format(m_d.split("_")[-2],nreals)
        else:

            name = "{0}, {1} realizations".format(m_d.split("_")[-2],m_d.split("_")[-1].replace("reals",""))
            nreals = int(m_d.split("_")[-1].replace("reals",""))
        if nreals not in name_dict:
            name_dict[nreals] = []
        name_dict[nreals].append(name)
        names.append(name)
        pes[name] = p
        oes[name] = o
        phidfs[name] = phidf
    
    #print(min_phi)
    reals = list(name_dict.keys())
    reals.sort()

    # now filter
    #thresh = min_phi * 5
    thresh = 10
    corrupt_thresh = min_phi * 1.1
    names = []
    for nreals in reals:
        m_ds = name_dict[nreals]
        m_ds.sort()
        for m_d in m_ds:
            p = pes[m_d]
            o = oes[m_d]
            phidf = phidfs[m_d]
            #print(m_d,len(p))
            assert len(p) == 2
            assert len(o) == 2
            ppost = p[-1]
            opost = o[-1]
            phipost = phidf.iloc[-1,:].copy()
            phipost = phipost.iloc[6:]
            #print(phipost)
            if "corrupt" in m_d:
                phipost = phipost.loc[phipost<=corrupt_thresh]
            else:    
                phipost = phipost.loc[phipost<=thresh]
            print(m_d,phipost.shape)
            pes[m_d][-1] = ppost.loc[phipost.index.values,:].copy()
            oes[m_d][-1] = opost.loc[phipost.index.values,:].copy()
            names.append(m_d)   
                 

    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    
    with PdfPages("collapse_compare_2ax.pdf") as pdf:
        for par in pst.adj_par_names:
            fig,axes = plt.subplots(2,1,figsize=(8.5,8.5))
            mn,mx = 1e30,-1e30
            for im_d,m_d in enumerate(names):
                ax = axes[0]
                if "corrupt" in m_d:
                    ax = axes[1]

                p = pes[m_d]
                o = oes[m_d]
                #print(m_d,len(p))
                color = plt.cm.jet(np.linspace(0, 1, len(names))) 
                pp = p[-1]
                if "stage" not in par:
                    pp.loc[:,par] = np.log10(pp.loc[:,par].values)
                    pp.loc[:,par].plot(ax=ax,kind="hist",color=color[im_d],alpha=0.5,label=m_d,density=True)
                mn = min(mn,ax.get_xlim()[0])
                mx = max(mn,ax.get_xlim()[1])
            axes[0].set_title("pure",loc="left")
            axes[1].set_title("corrupt",loc="left")
            
            for ax in axes:
                ax.set_xlim(mn,mx)
                ax.set_yticks([])
                ax.grid()
                ax.legend(loc="upper left")

            plt.tight_layout()
            pdf.savefig()
            plt.close(fig)
             

    with PdfPages("collapse_compare.pdf") as pdf:
        for par in pst.adj_par_names:
            fig,axes = plt.subplots(len(names),1,figsize=(8.5,8.5))
            mn,mx = 1e30,-1e30
            for m_d,ax in zip(names,axes):
                p = pes[m_d]
                o = oes[m_d]
                #print(m_d,len(p))

                colors = ["0.5","b"]
                labels = ["prior","posterior"]
                for pp,oo,c,l in zip(p,o,colors,labels):
                    if "stage" not in par:
                        pp.loc[:,par] = np.log10(pp.loc[:,par].values)
                    pp.loc[:,par].plot(ax=ax,kind="hist",color=c,alpha=0.5,label=l,density=True)
                ax.set_title("{0}, {1} {2} posterior realizations".format(par,m_d,pp.shape[0]),loc="left")
                ax.set_yticks([])
                mn = min(mn,ax.get_xlim()[0])
                mx = max(mn,ax.get_xlim()[1])
                print(m_d)
            
            for ax in axes:
                ax.set_xlim(mn,mx)
            plt.tight_layout()
            pdf.savefig()
            plt.close(fig)
             
        for obs in pst.obs_names:
            fig,axes = plt.subplots(len(names),1,figsize=(8.5,8.5))
            mn,mx = 1e30,-1e30
            for m_d,ax in zip(names,axes):
                p = pes[m_d]
                o = oes[m_d]
                colors = ["0.5","b"]
                labels = ["prior","posterior"]
                for pp,oo,c,l in zip(p,o,colors,labels):
                    oo.loc[:,obs].plot(kind="hist",color=c,alpha=0.5,label=l,ax=ax,density=True)
                ax.set_title("{0}, {1}, {2} posterior realizations".format(obs,m_d,oo.shape[0]),loc="left")
                mn = min(mn,ax.get_xlim()[0])
                mx = max(mn,ax.get_xlim()[1])
            for ax in axes:
                ax.set_xlim(mn,mx)
                
            plt.tight_layout()
            pdf.savefig()
            plt.close(fig)



            
    print(m_ds)


def tenpar_uniform_invest():
    model_d = "ies_10par_xsec"
    pyemu.Ensemble.reseed()
    base_d = os.path.join(model_d, "template")
    new_d = os.path.join(model_d, "test_template_geo")
    if os.path.exists(new_d):
        shutil.rmtree(new_d)
    shutil.copytree(base_d, new_d)
    print(platform.platform().lower())
    pst = pyemu.Pst(os.path.join(new_d, "pest.pst"))
    print(pst.model_command)
    obs = pst.observation_data
    obs.loc[:,"weight"] = 0.0
    obs.loc["h01_03","weight"] = 1.0
    par = pst.parameter_data
    par.loc[:,"partrans"] = "log"
    par.loc["stage","partrans"] = "fixed"
    v = pyemu.geostats.ExpVario(contribution=1.0,a=10)
    gs = pyemu.geostats.GeoStruct(variograms=v)
    x = np.cumsum(np.ones(10))
    y = np.ones(10)
    print(pst.adj_par_names)
    cov = gs.covariance_matrix(x,y,names = pst.adj_par_names).to_dataframe()
    cov.loc["stage",:] = 0.0
    cov.loc[:,"stage"] = 0.0
    
    cov.loc["stage","stage"] = 0.1
    #import matplotlib.pyplot as plt
    #plt.imshow(cov.values)
    #plt.show()
    par.loc["stage","partrans"] = "none"

    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst,cov=pyemu.Cov.from_dataframe(cov),num_reals=20)
    pe.enforce()
    pe.to_csv(os.path.join(new_d,"geoprior.csv"))
    pst.pestpp_options["ies_par_en"] = "geoprior.csv"
    pst.control_data.noptmax = 10
    pst.write(os.path.join(new_d,"pest.pst"))
    pyemu.os_utils.run("pestpp-ies pest.pst",cwd=new_d)

    for i,val in enumerate(np.linspace(1,5,pe.shape[0])):
        print(val)
        pe.iloc[i,:5] = val
    new_d = os.path.join(model_d, "test_template_uniform")
    if os.path.exists(new_d):
        shutil.rmtree(new_d)
    shutil.copytree(base_d, new_d)
    pe.to_csv(os.path.join(new_d,"uniprior.csv"))
    pst.pestpp_options["ies_par_en"] = "uniprior.csv"
    pst.control_data.noptmax = 10
    pst.write(os.path.join(new_d,"pest.pst"))
    pyemu.os_utils.run("pestpp-ies pest.pst",cwd=new_d)


def sweep_large_xfer_test():

    model_d = "ies_10par_xsec"
    t_d = os.path.join(model_d,"template")
    m_d = os.path.join(model_d,"master_sweep_xfer")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pst = pyemu.Pst(os.path.join(t_d,"pest.pst"))
    num_reals = 3
    pe = pyemu.ParameterEnsemble.from_uniform_draw(pst,num_reals=num_reals)#.loc[:,pst.par_names[:2]]

    pe.to_csv(os.path.join(t_d,"sweep_in.csv"))
    pe._df.index = pe.index.map(str)
    print(pe.index)
    pe.to_dense(os.path.join(t_d,"sweep_in.bin"))

    dimen = 1000
    cnames = ["col{0}".format(i) for i in range(dimen)]
    rnames = ["row{0}".format(i) for i in range(dimen)]
    vals = np.random.random((dimen,dimen))
    pyemu.Matrix(x=vals,row_names=rnames,col_names=cnames).to_dense(os.path.join(t_d,"matrix.bin"))

    pst.pestpp_options["ies_par_en"] = "sweep_in.csv"
    pst.pestpp_options["sweep_forgive"] = True
    pst.pestpp_options["sweep_parameter_file"] = "sweep_in.bin"
    pst.control_data.noptmax = -1
    pst.pestpp_options.pop("ies_num_reals",None)
    pst.write(os.path.join(t_d,"pest_forgive.pst"))
    pst.pestpp_options["sweep_output_file"] = "sweep_out.bin"
    pst.pestpp_options["sweep_chunk"] = 9
    pst.pestpp_options["ies_include_base"] = False
    pst.pestpp_options["panther_transfer_on_finish"] = "matrix.bin"
    pst.write(os.path.join(t_d,"pest_forgive.pst"))
    m_d = os.path.join(model_d,"master_sweep_bin_base")

    pyemu.os_utils.start_workers(t_d, exe_path.replace("-ies","-swp"), "pest_forgive.pst", 1, master_dir=m_d,
                           worker_root=model_d,port=port)

    for i in range(num_reals):
        fname = os.path.join(m_d,"ftx_{0}.matrix.bin".format(i))
        assert os.path.exists(fname)
        vals2 = pyemu.Matrix.from_binary(fname).x
        diff = np.abs(vals - vals2).sum()
        print(fname,diff)
        assert diff < 1e-10



def large_fake_test():
    root_d = "large_fake_test"
    if os.path.exists(root_d):
        shutil.rmtree(root_d)
    os.makedirs(root_d)
    t_d = os.path.join(root_d,"template")
    os.makedirs(t_d)

    npar = 10000
    nobs = 10000
    nzobs = 1000

    in_name = os.path.join(t_d,"in.dat")
    out_name = os.path.join(t_d,"out.dat")

    f_in = open(in_name,'w')
    f_tpl = open(in_name+".tpl",'w')
    f_tpl.write("ptf ~\n")

    for i in range(npar):
        f_in.write("par{0:06d},{0}\n".format(i))
        f_tpl.write("par{0:06d},~    par{0:06d}     ~\n".format(i))
    f_tpl.close()
    f_in.close()


    f_out = open(out_name,'w')
    f_ins = open(out_name+".ins",'w')
    f_ins.write("pif ~\n")
    for i in range(nobs):
        f_out.write("obs{0:06d},{0}\n".format(i))
        f_ins.write("l1 ~,~ !obs{0:06d}!\n".format(i))
    f_out.close()
    f_ins.close()

    pst = pyemu.Pst.from_io_files(in_name+".tpl",in_name,out_name+".ins",out_name,pst_path='.')
    obs = pst.observation_data
    obs["weight"] = 0.0
    obs.loc[obs.index[:nzobs],"weight"] = 1.0

    par = pst.parameter_data
    par["parlbnd"] = 0
    par["parubnd"] = npar
    par["partrans"] = "none"

    pst.control_data.noptmax = 0
    pst.write(os.path.join(t_d,"test.pst"),version=2)
    #pyemu.utils.helpers.setup_fake_forward_run(pst, new_pst_name, org_cwd='.', bak_suffix='._bak', new_cwd='.')
    pst = pyemu.helpers.setup_fake_forward_run(pst,"fake.pst",org_cwd=t_d,new_cwd=t_d)

    pst.write(os.path.join(t_d,"fake.pst"),version=2)
    
    frun_name = os.path.join(t_d,"fake_forward_run.py")
    lines = open(frun_name,'r').readlines()
    with open(frun_name,'w') as f:
        for line in lines:
            f.write(line)
        f.write("import time\n")
        f.write("time.sleep(10)\n")

    pyemu.os_utils.run("{0} fake.pst".format(exe_path),cwd=t_d)
    pst.control_data.noptmax = -1
    pst.write(os.path.join(t_d,"fake.pst"),version=2)

    m_d = t_d.replace("template","master")
    pyemu.os_utils.start_workers(t_d, exe_path, "fake.pst", 5, master_dir=m_d,
                           worker_root=root_d,port=port,verbose=True)



def mf6_v5_ies_nonpersistent_test():
    model_d = "mf6_freyberg"

    t_d = os.path.join(model_d,"template")
    m_d = os.path.join(model_d,"master_ies_glm_loc")
    #if os.path.exists(m_d):
    #    shutil.rmtree(m_d)
    pst = pyemu.Pst(os.path.join(t_d,"freyberg6_run_ies.pst"))
    pst.control_data.noptmax = 0
    pst.write(os.path.join(t_d,"freyberg6_run_ies.pst"))
    pyemu.os_utils.run("{0} freyberg6_run_ies.pst".format(exe_path),cwd=t_d)

    pst.control_data.noptmax = -1
    par = pst.parameter_data

    eff_lb = (par.parlbnd + (np.abs(par.parlbnd.values)*.01)).to_dict()
    eff_ub = (par.parubnd - (np.abs(par.parlbnd.values)*.01)).to_dict()
    log_idx = par.partrans.apply(lambda x: x=="log").to_dict()
    for p,log in log_idx.items():
        if log:
            lb = np.log10(par.loc[p,"parlbnd"])
            eff_lb[p] = (lb + (np.abs(lb)*.01))
            ub = np.log10(par.loc[p,"parubnd"])
            eff_ub[p] = (ub - (np.abs(ub)*.01))

    pargp_map = par.groupby(par.pargp).groups
    print(pargp_map)


    m_d = None
    m_d = os.path.join(model_d, "master_ies_nonpersist")
    if os.path.exists(m_d):
         shutil.rmtree(m_d)
    pst = pyemu.Pst(os.path.join(t_d, "freyberg6_run_ies.pst"))
    pst.pestpp_options.pop("ies_localizer",None)
    pst.pestpp_options.pop("ies_autoadaloc",None)
    pst.pestpp_options["panther_persistent_workers"] = False
    pst.pestpp_options["ies_bad_phi_sigma"] = 2.5
    pst.pestpp_options["ies_num_reals"] = 100
    pst.pestpp_options["ensemble_output_precision"] = 40
    pst.pestpp_options["panther_master_timeout_milliseconds"] = 1000
    pst.control_data.noptmax = -1
    pst_name = "freyberg6_run_ies_nonpersist.pst"
    pst.write(os.path.join(t_d, pst_name))
    num_workers = 15
    pyemu.os_utils.start_workers(t_d, exe_path, pst_name, num_workers=num_workers,
                                 master_dir=m_d, worker_root=model_d, port=port)

    found = 0
    with open(os.path.join(m_d,pst_name.replace(".pst",".rmr")),'r') as f:
        for line in f:
            if "using non-persistent agents" in line:
                found += 1
    print("found: ",found,"num workers: ",num_workers)
    assert found ==  num_workers


def parse_pst_test():
    t_d = os.path.join("parse_pst_testfiles","kerry")
    pyemu.os_utils.run("{0} control.pst".format(exe_path),cwd=t_d)



if __name__ == "__main__":
    #parse_pst_test()
    #basic_test()
    nonascii_path_test()

    #mf6_v5_ies_nonpersistent_test()
    #large_fake_test()
    #exit()
    #sweep_large_xfer_test()
    #sweep_bin_test()
    #exit()
    # mf6_v5_sen_test()
    #tie_by_group_test()
    #tenpar_uniform_invest()
    #tenpar_collapse_invest()
    #plot_collapse_invest()

    #run()
    #mf6_v5_ies_test()
    #prep_ends()
    #sweep_bin_test()
    # mf6_v5_sen_test()
    #ext_stdcol_test()
    #shutil.copy2(os.path.join("..","exe","windows","x64","Debug","pestpp-glm.exe"),os.path.join("..","bin","win","pestpp-glm.exe"))
    #shutil.copy2(os.path.join("..", "exe", "windows", "x64", "Debug", "pestpp-ies.exe"),
    #             os.path.join("..", "bin", "win", "pestpp-ies.exe"))
    #ins_missing_e_test()
    #basic_test()
    #agnostic_path_test()
    #glm_long_name_test()
    #sen_plusplus_test()
    #parchglim_test()
    #unc_file_test()
    #cmdline_test()
    #secondary_marker_test()
    #basic_test("ies_10par_xsec")
    #glm_save_binary_test()
    #sweep_forgive_test()
    #inv_regul_test()
    #tie_by_group_test()
    #sen_basic_test()
    #salib_verf()
    #tplins1_test()
    #secondary_marker_test()
    #ext_stdcol_test()

    # parallel_consist_test()
    # ext_stdcol_test()
    #sen_grp_test()
    #da_prep_4_freyberg_batch()
    # da_prep_4_mf6_freyberg_seq()
    # basic_test()
    # da_mf6_freyberg_smoother_test()
    # da_mf6_freyberg_test_1()

    # da_prep_4_mf6_freyberg_seq_tbl()
    # da_mf6_freyberg_test_2()
    #shutil.copy2(os.path.join("..","exe","windows","x64","Debug","pestpp-ies.exe"),os.path.join("..","bin","win","pestpp-ies.exe"))
    #tplins1_test()
    
    #fr_timeout_test()
    #mf6_v5_ies_test()
    #mf6_v5_sen_test()

    #shutil.copy2(os.path.join("..","exe","windows","x64","Debug","pestpp-opt.exe"),os.path.join("..","bin","win","pestpp-opt.exe"))
    #mf6_v5_opt_stack_test()
    # mf6_v5_glm_test()
    # mf6_v5_ies_test()
    #cmdline_test()
    #basic_sqp_test()
    #mf6_v5_ies_test()
    #fr_timeout_test()
    #fr_fail_test()
    #tplins1_test()

    #mf6_v5_glm_test()
    #mf6_v5_sen_test()
