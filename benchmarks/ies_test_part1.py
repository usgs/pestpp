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
    bin_path = os.path.join(bin_path,"linux")
elif "darwin" in platform.platform().lower() or "macos" in platform.platform().lower():
    bin_path = os.path.join(bin_path,"mac")
else:
    bin_path = os.path.join(bin_path,"win")

bin_path = os.path.abspath("test_bin")
os.environ["PATH"] += os.pathsep + bin_path


# case of either appveyor, travis or local
print(platform.platform().lower())
if os.path.exists(os.path.join("pestpp","bin")):
    bin_path = os.path.join("..","..","pestpp","bin")
else:
    bin_path = os.path.join("..","..","..","..","pestpp","bin")
    
if "windows" in platform.platform().lower():
    exe_path = os.path.join(bin_path, "win", "pestpp-ies.exe")
elif "darwin" in platform.platform().lower() or "macos" in platform.platform().lower():
    exe_path = os.path.join(bin_path,  "mac", "pestpp-ies")
else:
    exe_path = os.path.join(bin_path, "linux", "pestpp-ies")





noptmax = 3

compare_files = ["pest.phi.actual.csv", "pest.phi.meas.csv", "pest.phi.regul.csv",
                 "pest.{0}.par.csv".format(noptmax), "pest.{0}.obs.csv".format(noptmax),
                 "pest.{0}.par.csv".format(0), "pest.obs+noise.csv"]
diff_tol = 1.0e-6
port = 4016
num_reals = 10

def tenpar_narrow_range_test():
    """tenpar narrow test"""
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_narrow_test")
    template_d = os.path.join(model_d, "test_template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pst_name = os.path.join(test_d, "pest.pst")
    pst = pyemu.Pst(pst_name)
    par = pst.parameter_data
    #par.loc[:, "partrans"] = "fixed"
    par.loc[:, "parubnd"] = 1.0e+10 #par.parval1 * 1.0001
    par.loc[:, "parlbnd"] = 1.0e-10 #par.parval1 * 0.9999
    #par.loc[pst.par_names[:2], "partrans"] = "none"
    #par.loc[pst.par_names[0],"pargp"] = "stage"

    x = np.zeros((pst.npar_adj, pst.npar_adj)) + 1.0e-11
    for i in range(pst.npar_adj):
        x[i, i] = 5.0e-10
    cov = pyemu.Cov(x, names=pst.adj_par_names)
    cov.to_ascii(os.path.join(test_d, "prior.cov"))
    num_reals = 100000
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst, cov, num_reals=num_reals)
    #pe.enforce()
    #pe.to_csv(os.path.join(test_d,"pyemu_draws.csv"))

    pst.control_data.noptmax = -2

    pst.observation_data.loc[pst.nnz_obs_names,"weight"] *= 1.5

    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = num_reals
    pst.pestpp_options["parcov_filename"] = "prior.cov"
    pst.pestpp_options["ies_enforce_bounds"] = False
    pst.pestpp_options["ies_include_base"] = False
    pst.pestpp_options["ies_group_draws"] = False
    pst.pestpp_options["ies_verbose_level"] = 1
    pst.write(pst_name)
    pyemu.helpers.run(exe_path + " pest.pst", cwd=test_d)

    df = pd.read_csv(os.path.join(test_d, "pest.0.par.csv"), index_col=0)
    df.columns = [c.lower() for c in df.columns]
    p1, p2 = pst.adj_par_names[:2]
    v1,v2 = pe.loc[:,p1].var(),df.loc[:,p1].var()
    diff = np.abs(100 * ((v1 - v2) / v1))
    print(v1,v2,diff)
    assert diff < 2.0

    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = num_reals
    pst.pestpp_options["parcov_filename"] = "prior.cov"
    pst.pestpp_options["ies_enforce_bounds"] = False
    pst.pestpp_options["ies_include_base"] = False
    pst.pestpp_options["ies_group_draws"] = True
    pst.pestpp_options["ies_verbose_level"] = 1
    pst.write(pst_name)
    pyemu.helpers.run(exe_path + " pest.pst", cwd=test_d)

    df = pd.read_csv(os.path.join(test_d, "pest.0.par.csv"), index_col=0)
    df.columns = [c.lower() for c in df.columns]
    p1, p2 = pst.adj_par_names[:2]
    v1, v2 = pe.loc[:, p1].var(), df.loc[:, p1].var()
    diff = np.abs(100 * ((v1 - v2) / v1))
    print(v1, v2, diff)
    assert diff < 2.0

def tenpar_full_cov_test():
    """tenpar full cov test"""
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_full_cov_test")
    template_d = os.path.join(model_d, "test_template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d,test_d)
    pst_name = os.path.join(test_d,"pest.pst")
    pst = pyemu.Pst(pst_name)
    par = pst.parameter_data
    par.loc[:,"partrans"] = "fixed"
    par.loc[:,"parubnd"] = 1.0e+10
    par.loc[:,"parlbnd"] = -1.0e+10
    par.loc[pst.par_names[:2],"partrans"] = "none"
    par.loc[:,"parchglim"] = "relative"
    x = np.zeros((pst.npar_adj,pst.npar_adj)) + 0.25
    for i in range(pst.npar_adj):
        x[i,i] = 1.0
    cov = pyemu.Cov(x,names=pst.adj_par_names)
    cov.to_ascii(os.path.join(test_d,"prior.cov"))
    num_reals = 10000
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst,cov,num_reals=num_reals)
    pe.enforce()

    pst.control_data.noptmax = -2
    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = num_reals
    pst.pestpp_options["parcov_filename"] = "prior.cov"
    pst.pestpp_options["ies_verbose_level"] = 1
    pst.pestpp_options["ies_num_threads"] = 4
    pst.write(pst_name)
    pyemu.helpers.run(exe_path+" pest.pst",cwd=test_d)

    df = pd.read_csv(os.path.join(test_d,"pest.0.par.csv"),index_col=0)
    df.columns = [c.lower() for c in df.columns]

    p1,p2 = pst.adj_par_names
    pe_corr = pe.corr().loc[p1,p2]
    df_corr = df.corr().loc[p1,p2]
    diff = np.abs((pe_corr - df_corr)/pe_corr)
    assert diff < 0.25,"{0},{1},{2}".format(pe_corr,df_corr,diff)

    par.loc[pst.adj_par_names,"partrans"] = "log"
    par.loc[:,"parlbnd"] = 1.0e-10
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst, cov, num_reals=num_reals)
    pe.enforce()
    pst.write(pst_name)
    pyemu.helpers.run(exe_path + " pest.pst", cwd=test_d)
    df = pd.read_csv(os.path.join(test_d, "pest.0.par.csv"), index_col=0)
    df.columns = [c.lower() for c in df.columns]

    p1, p2 = pst.adj_par_names
    pe_corr = pe.apply(lambda x: np.log10(x)).corr().loc[p1, p2]
    df_corr = df.apply(lambda x: np.log10(x)).corr().loc[p1, p2]
    diff = np.abs((pe_corr - df_corr) / pe_corr)
    assert diff < 0.25, "{0},{1},{2}".format(pe_corr, df_corr, diff)

    pst.control_data.noptmax = 1
    pst.pestpp_options["ies_num_reals"] = 10
    pst.write(pst_name)
    pyemu.helpers.run(exe_path+" pest.pst",cwd=test_d)


def tenpar_subset_test():
    """test that using subset gets the same results in the single lambda case"""
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_subset_test")
    template_d = os.path.join(model_d, "test_template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    # shutil.copytree(base_d,test_d)
    pst = pyemu.Pst(os.path.join(template_d, "pest.pst"))
    pst.control_data.noptmax = 3

    # first without subset
    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options["ies_lambda_mults"] = 1.0
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["ies_accept_phi_fac"] = 10000.0

    pst.pestpp_options["ies_subset_size"] = 21
    pst.write(os.path.join(template_d, "pest.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest.pst", num_workers=15,
                               worker_root=model_d, master_dir=test_d,port=port)
    df_base = pd.read_csv(os.path.join(test_d, "pest.phi.meas.csv"),index_col=0)

    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options["ies_lambda_mults"] = 1.0
    pst.pestpp_options["ies_subset_size"] = 5
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["ies_accept_phi_fac"] = 100.0

    pst.write(os.path.join(template_d, "pest.pst"))
    test_d = os.path.join(model_d, "master_subset_test1")
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    pyemu.os_utils.start_workers(template_d, exe_path, "pest.pst", num_workers=15,
                               worker_root=model_d, master_dir=test_d,port=port)
    df_sub = pd.read_csv(os.path.join(test_d, "pest.phi.meas.csv"),index_col=0)
    diff = (df_sub - df_base).apply(lambda x: np.abs(x))
    diff = diff.iloc[:,6:]
    print(diff.max())
    print(df_sub.iloc[-1,:])
    print(df_base.iloc[-1,:])
    assert diff.max().max() == 0.0

def test_freyberg_full_cov():
    """freyberg full cov test"""
    model_d = "ies_freyberg"
    test_d = os.path.join(model_d, "master_draw_test")
    template_d = os.path.join(model_d, "test_template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d,test_d)
    pst = pyemu.Pst(os.path.join(test_d, "pest.pst"))
    pst.parameter_data.loc[:,"partrans"] = "log"
    
    pst.control_data.noptmax = -2
    pst.pestpp_options = {}
    num_reals = 5000

    #diagonal cov
    #pst.pestpp_options["parcov_filename"] = "prior.jcb"
    pst.pestpp_options["ies_num_reals"] = num_reals
    pst.pestpp_options["ies_include_base"] = "false"

    pst.write(os.path.join(test_d, "pest.pst"))
    #cov = pyemu.Cov.from_binary(os.path.join(test_d, "prior.jcb"))
    cov = pyemu.Cov.from_parameter_data(pst)

    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst, cov, num_reals)
    pe.to_csv(os.path.join(test_d, "pyemu_pe.csv"))

    # pyemu.os_utils.start_workers(template_d, exe_path, "pest.pst", num_workers=10,
    #                            worker_root=model_d, master_dir=test_d)
    # pyemu.helpers.run(exe_path + " pest.pst", cwd=test_d)
    # print("loading df")
    # df = pd.read_csv(os.path.join(test_d, "pest.0.par.csv"), index_col=0).apply(np.log10)
    # df.columns = [c.lower() for c in df.columns]
    # pe = pe.apply(np.log10)
    # pe_corr = pe.corr()
    # df_corr = df.corr()

    # diff_tol = 0.05

    # for c in df.columns:
    #     if c not in pe.columns:
    #         continue

    #     m1, m2 = pe.loc[:, c].mean(), df.loc[:, c].mean()
    #     s1, s2 = pe.loc[:, c].std(), df.loc[:, c].std()
    #     mdiff = np.abs((m1 - m2))
    #     sdiff = np.abs((s1 - s2))
    #     print(c, mdiff, sdiff)
    #     assert mdiff < diff_tol, "mean fail {0}:{1},{2},{3}".format(c, m1, m2, mdiff)
    #     assert sdiff < diff_tol, "std fail {0}:{1},{2},{3}".format(c, s1, s2, sdiff)

    # # look for bias
    # diff = df - pe
    # assert diff.mean().mean() < 0.01


    #full cov
    pst.pestpp_options["parcov_filename"] = "prior.jcb"
    pst.pestpp_options["ies_num_reals"] = num_reals
    pst.pestpp_options["ies_include_base"] = "false"
    pst.pestpp_options["ies_group_draws"] = 'true'
    pst.parameter_data.loc[pst.par_names[0],"pargp"] = "test"
    pst.write(os.path.join(test_d,"pest.pst"))
    cov = pyemu.Cov.from_binary(os.path.join(test_d,"prior.jcb"))

    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst,cov,num_reals)
    pe.to_csv(os.path.join(test_d,"pyemu_pe.csv"))

    # pyemu.os_utils.start_workers(template_d, exe_path, "pest.pst", num_workers=10,
    #                            worker_root=model_d, master_dir=test_d)
    pyemu.helpers.run(exe_path+" pest.pst",cwd=test_d)
    df = pd.read_csv(os.path.join(test_d, "pest.0.par.csv"), index_col=0).apply(lambda x: np.log10(x))
    df.columns = [c.lower() for c in df.columns]
    pe = pe.apply(lambda x: np.log10(x))
    pe_corr = pe.corr()
    df_corr = df.corr()

    for i,p1 in enumerate(pst.adj_par_names):
        for p2 in pst.adj_par_names[i+1:]:
            c1 = pe_corr.loc[p1,p2]
            c2 = df_corr.loc[p1,p2]
            #print(p1,p2,c1,c2)

    diff_tol = 0.05

    for c in df.columns:
        if c not in pe.columns:
            continue

        m1, m2 = pe.loc[:,c].mean(), df.loc[:,c].mean()
        s1,s2 =  pe.loc[:,c].std(), df.loc[:,c].std()
        mdiff = np.abs((m1 - m2))
        sdiff = np.abs((s1 - s2))
        #print(c,mdiff,sdiff)
        assert mdiff < diff_tol,"mean fail {0}:{1},{2},{3}".format(c,m1,m2,mdiff)
        assert sdiff < diff_tol,"std fail {0}:{1},{2},{3}".format(c,s1,s2,sdiff)

    #look for bias
    diff = df - pe
    assert diff.mean().mean() < 0.01


def test_freyberg_full_cov_reorder():
    """freyberg full cov reorder test"""
    model_d = "ies_freyberg"
    test_d = os.path.join(model_d, "master_draw_test")
    template_d = os.path.join(model_d, "test_template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d,test_d)
    pst = pyemu.Pst(os.path.join(test_d, "pest.pst"))
    par = pst.parameter_data
    par.loc[:,"partrans"] = "log"
    par.loc[pst.par_names[:5],"pargp"] = pst.par_groups[-1]

    pst.control_data.noptmax = -2
    pst.pestpp_options = {}
    num_reals = 5000

    #diagonal cov
    #pst.pestpp_options["parcov_filename"] = "prior.jcb"
    pst.pestpp_options["ies_num_reals"] = num_reals
    pst.pestpp_options["ies_num_threads"] = 2
    pst.pestpp_options["ies_include_base"] = "false"

    pst.write(os.path.join(test_d, "pest.pst"))
    #cov = pyemu.Cov.from_binary(os.path.join(test_d, "prior.jcb"))
    cov = pyemu.Cov.from_parameter_data(pst)

    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst, cov, num_reals)
    pe.to_csv(os.path.join(test_d, "pyemu_pe.csv"))

    # pyemu.os_utils.start_workers(template_d, exe_path, "pest.pst", num_workers=10,
    #                            worker_root=model_d, master_dir=test_d)
    # pyemu.helpers.run(exe_path + " pest.pst", cwd=test_d)
    # print("loading df")
    # df = pd.read_csv(os.path.join(test_d, "pest.0.par.csv"), index_col=0).apply(np.log10)
    # df.columns = [c.lower() for c in df.columns]
    # pe = pe.apply(np.log10)
    # pe_corr = pe.corr()
    # df_corr = df.corr()

    # diff_tol = 0.05

    # for c in df.columns:
    #     if c not in pe.columns:
    #         continue

    #     m1, m2 = pe.loc[:, c].mean(), df.loc[:, c].mean()
    #     s1, s2 = pe.loc[:, c].std(), df.loc[:, c].std()
    #     mdiff = np.abs((m1 - m2))
    #     sdiff = np.abs((s1 - s2))
    #     print(c, mdiff, sdiff)
    #     assert mdiff < diff_tol, "mean fail {0}:{1},{2},{3}".format(c, m1, m2, mdiff)
    #     assert sdiff < diff_tol, "std fail {0}:{1},{2},{3}".format(c, s1, s2, sdiff)

    # # look for bias
    # diff = df - pe
    # assert diff.mean().mean() < 0.01


    #full cov
    pst.pestpp_options["parcov_filename"] = "prior.jcb"
    pst.pestpp_options["ies_num_reals"] = num_reals
    pst.pestpp_options["ies_include_base"] = "false"
    pst.pestpp_options["ies_group_draws"] = 'true'
    pst.pestpp_options["ies_enforce_bounds"] = False
    pst.parameter_data.loc[pst.par_names[0],"pargp"] = "test"

    pst.write(os.path.join(test_d,"pest.pst"))
    cov = pyemu.Cov.from_binary(os.path.join(test_d,"prior.jcb"))

    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst,cov,num_reals)
    pe.to_csv(os.path.join(test_d,"pyemu_pe.csv"))

    # pyemu.os_utils.start_workers(template_d, exe_path, "pest.pst", num_workers=10,
    #                            worker_root=model_d, master_dir=test_d)
    pyemu.helpers.run(exe_path+" pest.pst",cwd=test_d)
    df = pd.read_csv(os.path.join(test_d, "pest.0.par.csv"), index_col=0).apply(lambda x: np.log10(x))
    df.columns = [c.lower() for c in df.columns]
    pe = pe._df.apply(lambda x: np.log10(x))
    #print('corr')
    #pe_corr = pe.corr()
    #df_corr = df.corr()

    # for i,p1 in enumerate(pst.adj_par_names):
    #     for p2 in pst.adj_par_names[i+1:]:
    #         c1 = pe_corr.loc[p1,p2]
    #         c2 = df_corr.loc[p1,p2]
    #         #print(p1,p2,c1,c2)

    diff_tol = 0.05
    print('compare')
    for c in df.columns[::2]:
        if c not in pe.columns:
            continue

        m1, m2 = pe.loc[:,c].mean(), df.loc[:,c].mean()
        s1,s2 =  pe.loc[:,c].std(), df.loc[:,c].std()
        mdiff = np.abs((m1 - m2))
        sdiff = np.abs((s1 - s2))
        #print(c,mdiff,sdiff)
        assert mdiff < diff_tol,"mean fail {0}:{1},{2},{3}".format(c,m1,m2,mdiff)
        assert sdiff < diff_tol,"std fail {0}:{1},{2},{3}".format(c,s1,s2,sdiff)

    #look for bias
    diff = df - pe
    assert diff.mean().mean() < 0.01


def test_freyberg_full_cov_reorder_run():
    """freyberg full cov reorder run test"""
    model_d = "ies_freyberg"
    test_d = os.path.join(model_d, "master_draw_test")
    template_d = os.path.join(model_d, "test_template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d,test_d)
    pst = pyemu.Pst(os.path.join(test_d, "pest.pst"))
    par = pst.parameter_data
    par.loc[:,"partrans"] = "log"
    #par.loc[pst.par_names[:5],"pargp"] = pst.par_groups[-1]

    pst.control_data.noptmax = 1
    pst.pestpp_options = {}
    num_reals = 30

    #diagonal cov
    #pst.pestpp_options["parcov_filename"] = "prior.jcb"

    pst.pestpp_options["parcov_filename"] = "prior.jcb"
    pst.pestpp_options["ies_num_reals"] = num_reals
    pst.pestpp_options["ies_include_base"] = "true"
    pst.pestpp_options["ies_group_draws"] = 'true'
    pst.pestpp_options["ies_lambda_mults"] = [0.9,1.1]
    pst.pestpp_options["lambda_scale_fac"] = 1.0
   # pst.parameter_data.loc[pst.par_names[0],"pargp"] = "test"
    
    pst.write(os.path.join(template_d, "pest.pst"))
    #cov = pyemu.Cov.from_binary(os.path.join(test_d, "prior.jcb"))
    cov = pyemu.Cov.from_parameter_data(pst)

    #pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst, cov, num_reals, use_homegrown=True)
    #pe.to_csv(os.path.join(test_d, "pyemu_pe.csv"))

    pyemu.os_utils.start_workers(template_d, exe_path, "pest.pst", num_workers=25,
                                worker_root=model_d, master_dir=test_d,port=port)
    

def invest():
    d = os.path.join("ies_freyberg","master_draw_test")
    df = pd.read_csv(os.path.join(d,"draws.dat"),sep='\s+',header=None)
    print(df.std().mean(),df.mean().mean())

    df1 = pd.read_csv(os.path.join(d,"pest.0.par.csv"), index_col=0)
    df2 = pd.read_csv(os.path.join(d,"pyemu_pe.csv"), index_col=0)
    df1.columns = [c.lower() for c in df1.columns]

    df1 = df1.apply(lambda x: np.log10(x))
    df2 = df2.apply(lambda x: np.log10(x))

    for p in df1.columns:
        print(p)
        #print(p)
        #print(df1.loc[:,p])
        #print(df2.loc[:,p])
        print(p,df1.loc[:, p].std(), df2.loc[:, p].std())
        #break

def eval_synth():
    model_d = "ies_synth"
    test_d = os.path.join(model_d,"master")
    template_d = os.path.join(model_d,"template")

    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    print("loading pst")
    pst = pyemu.Pst(os.path.join(template_d,"pest.pst"))
    pst.pestpp_options = {}
    pst.pestpp_options["ies_use_approx"] = "false"
    pst.pestpp_options["ies_use_prior_scaling"] = "true"
    pst.pestpp_options["ies_lambda_mults"] = [0.1,1.0]
    pst.pestpp_options["lambda_scale_fac"] = [0.9,1.1]
    pst.pestpp_options["ies_num_reals"] = 30
    pst.pestpp_options["ies_save_binary"] = True
    pst.control_data.noptmax = 2
    print("writing pst")
    pst.write(os.path.join(template_d,"pest.pst"))
    print("starting slaves")
    pyemu.os_utils.start_workers(template_d,exe_path,"pest.pst",num_workers=15,
        master_dir=test_d,worker_root=model_d,port=port)

def test_chenoliver():
    """chen and oliver test"""
    model_d = "ies_chenoliver"
    test_d = os.path.join(model_d,"master1")
    template_d = os.path.join(model_d,"template")

    # build the prior cov matrix
    x = np.zeros((1,1)) + 0.5
    cov = pyemu.Cov(x=x,names=["par"])
    cov.to_ascii(os.path.join(template_d,"prior.cov"))

    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d,test_d)
    pst = pyemu.Pst(os.path.join(template_d,"pest.pst"))
    pst.parameter_data.loc[:,"parlbnd"] = -1.0e+10
    pst.parameter_data.loc[:,"parubnd"] = 1.0+10
    pst.pestpp_options = {}
    pst.pestpp_options["parcov_filename"] = "prior.cov"
    pst.pestpp_options["ies_num_reals"] = 100
    pst.control_data.noptmax = 0
    par = pst.parameter_data
    par.loc[:,"parchglim"] = "relative"
    pst.write(os.path.join(test_d,"pest.pst"))
    pyemu.helpers.run(exe_path+" pest.pst",cwd=test_d)
    
    num_reals = 30
    noptmax = 3
    
    silent_master = False

    shutil.rmtree(test_d)
    
    pst.observation_data.loc[:,"weight"] = 0.25

    pst.pestpp_options = {}
    pst.pestpp_options["parcov_filename"] = "prior.cov"
    pst.pestpp_options["ies_num_reals"] = num_reals
    pst.pestpp_options["ies_lambda_mults"] = [0.01,1.0,10.0]
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["ies_initial_lambda"] = 0.001
    #pst.pestpp_options["ies_subset_size"] = 10
    pst.pestpp_options["ies_use_approx"] = "false"
    pst.pestpp_options["ies_use_prior_scaling"] = "true"
    pst.control_data.noptmax = noptmax

    pst.write(os.path.join(template_d,"pest.pst"))
    

    pyemu.os_utils.start_workers(template_d,exe_path,"pest.pst",num_workers=20,
        master_dir=test_d,worker_root=model_d,port=port,silent_master=silent_master)
    df_full_obs = pd.read_csv(os.path.join(test_d,"pest.{0}.obs.csv".format(noptmax)),index_col=0)
    df_full_par = pd.read_csv(os.path.join(test_d,"pest.{0}.par.csv".format(noptmax)),index_col=0)

    shutil.rmtree(test_d)
    pst.pestpp_options = {}
    pst.pestpp_options["parcov_filename"] = "prior.cov"
    pst.pestpp_options["ies_num_reals"] = num_reals
    #pst.pestpp_options["ies_lambda_mults"] = "0.01,1.0,100.0"
    pst.pestpp_options["ies_initial_lambda"] = 0.001
    #pst.pestpp_options["ies_subset_size"] = 10
    pst.pestpp_options["ies_use_prior_scaling"] = "true"
    pst.control_data.noptmax = noptmax
    pst.write(os.path.join(template_d,"pest.pst"))

    pyemu.os_utils.start_workers(template_d,exe_path,"pest.pst",num_workers=20,
        master_dir=test_d,worker_root=model_d,port=port,silent_master=silent_master)
    df_approx_obs = pd.read_csv(os.path.join(test_d,"pest.{0}.obs.csv".format(noptmax)),index_col=0)
    df_approx_par = pd.read_csv(os.path.join(test_d,"pest.{0}.par.csv".format(noptmax)),index_col=0)

    ax = plt.subplot(211)
    ax2 = plt.subplot(212)
    #ax.plot(df_full.loc[:,"mean"],color='b',label="full")
    #ax.plot(df_approx.loc[:,"mean"],color='g',label="approx")
    df_full_obs.obs.hist(bins=30,color='b',alpha=0.5, ax=ax)
    df_approx_obs.obs.hist(bins=30,color='0.5',alpha=0.5,ax=ax)
    df_full_par.par.hist(bins=30,color='b',alpha=0.5, ax=ax2)
    df_approx_par.par.hist(bins=30,color='0.5',alpha=0.5,ax=ax2)
    ax.set_title("full: {0:15g}, approx: {1:15g}".format(df_full_obs.obs.mean(),df_approx_obs.obs.mean()))
    ax2.set_title("full: {0:15g}, approx: {1:15g}".format(df_full_par.par.mean(),df_approx_par.par.mean()))
    plt.tight_layout()
    plt.savefig(os.path.join(model_d,"full_approx_ov16.png"))
    plt.close("all")
    d = np.abs(df_full_par.par.mean() - 5.8)
    #assert d < 0.05,"{0},{1}".format(d,df_full_par.PAR.mean())
    print("{0},{1}".format(d,df_full_par.par.mean()))
    d = np.abs(df_approx_par.par.mean() - 6.0)
    #assert d < 0.05,"{0},{1}".format(d,df_approx_par.PAR.mean())
    print("{0},{1}".format(d,df_approx_par.par.mean()))
    pst.observation_data.loc[:,"weight"] = 1.0

    shutil.rmtree(test_d)
    pst.pestpp_options = {}
    pst.pestpp_options["parcov_filename"] = "prior.cov"
    pst.pestpp_options["ies_num_reals"] = num_reals
    #pst.pestpp_options["ies_lambda_mults"] = "0.01,1.0,100.0"
    pst.pestpp_options["ies_initial_lambda"] = 0.001
    #pst.pestpp_options["ies_subset_size"] = 10
    pst.pestpp_options["ies_use_approx"] = "false"
    pst.pestpp_options["ies_use_prior_scaling"] = "true"
    pst.control_data.noptmax = noptmax
    pst.write(os.path.join(template_d,"pest.pst"))

    pyemu.os_utils.start_workers(template_d,exe_path,"pest.pst",num_workers=25,
        master_dir=test_d,worker_root=model_d,port=port,silent_master=silent_master)
    df_full_obs = pd.read_csv(os.path.join(test_d,"pest.{0}.obs.csv".format(noptmax)),index_col=0)
    df_full_par = pd.read_csv(os.path.join(test_d,"pest.{0}.par.csv".format(noptmax)),index_col=0)

    shutil.rmtree(test_d)
    pst.pestpp_options = {}
    pst.pestpp_options["parcov_filename"] = "prior.cov"
    pst.pestpp_options["ies_num_reals"] = num_reals
    #pst.pestpp_options["ies_lambda_mults"] = "0.01,1.0,100.0"
    pst.pestpp_options["ies_initial_lambda"] = 0.001
    pst.pestpp_options["ies_use_prior_scaling"] = "true"
    #pst.pestpp_options["ies_subset_size"] = 10
    pst.control_data.noptmax = noptmax
    pst.write(os.path.join(template_d,"pest.pst"))

    pyemu.os_utils.start_workers(template_d,exe_path,"pest.pst",num_workers=25,
        master_dir=test_d,worker_root=model_d,port=port,silent_master=silent_master)
    df_approx_obs = pd.read_csv(os.path.join(test_d,"pest.{0}.obs.csv".format(noptmax)),index_col=0)
    df_approx_par = pd.read_csv(os.path.join(test_d,"pest.{0}.par.csv".format(noptmax)),index_col=0)

    ax = plt.subplot(211)
    ax2 = plt.subplot(212)
    #ax.plot(df_full.loc[:,"mean"],color='b',label="full")
    #ax.plot(df_approx.loc[:,"mean"],color='g',label="approx")
    df_full_obs.obs.hist(bins=30,color='b',alpha=0.5, ax=ax)
    df_approx_obs.obs.hist(bins=30,color='0.5',alpha=0.5,ax=ax)
    df_full_par.par.hist(bins=30,color='b',alpha=0.5, ax=ax2)
    df_approx_par.par.hist(bins=30,color='0.5',alpha=0.5,ax=ax2)
    ax.set_title("full: {0:15g}, approx: {1:15g}".format(df_full_obs.obs.mean(),df_approx_obs.obs.mean()))
    ax2.set_title("full: {0:15g}, approx: {1:15g}".format(df_full_par.par.mean(),df_approx_par.par.mean()))
    plt.tight_layout()
    plt.savefig(os.path.join(model_d,"full_approx_ov1.png"))
    plt.close("all")

    d = np.abs(df_full_par.par.mean() - 5.8)

    print("{0},{1}".format(d,df_full_par.par.mean()))
    d = np.abs(df_approx_par.par.mean() - 6.0)
    print("{0},{1}".format(d,df_approx_par.par.mean()))
    #assert d < 0.05,"{0},{1}".format(d,df_approx_par.PAR.mean())

def eval_kirishima():

    model_d = "ies_kirishima"
    test_d = os.path.join(model_d, "master")
    template_d = os.path.join(model_d, "template")

    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    print("loading pst")
    pst = pyemu.Pst(os.path.join(template_d, "pest.pst"))
    pst.pestpp_options = {}
    #pst.pestpp_options["ies_num_reals"] = 300
    #pst.pestpp_options["ies_use_approx"] = "true"
    #pst.pestpp_options["ies_use_prior_scaling"] = "true"
    #pst.pestpp_options["ies_subset_size"] = 10
    #pst.pestpp_options["ies_lambda_mults"] = [0.1,1.0,10.0]

    #pst.pestpp_options["ies_initial_lambda"] = 1000.0
    #pst.pestpp_options["ies_bad_phi"] = 80000.0
    pst.control_data.noptmax = 10
    print("writing pst")
    pst.write(os.path.join(template_d, "pest.pst"))
    print("starting slaves")
    pyemu.os_utils.start_workers(template_d, exe_path, "pest.pst", num_workers=15, master_dir=test_d,
                               worker_root=model_d,port=port)

def test_freyberg_ineq():
    """freyberg ineq test"""
    model_d = "ies_freyberg"
    test_d = os.path.join(model_d, "master_ineq")
    template_d = os.path.join(model_d, "template")

    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    print("loading pst")
    pst = pyemu.Pst(os.path.join(template_d, "pest.pst"))
    pst.observation_data.loc[pst.nnz_obs_names[3],"obgnme"] = "less_than"
    pst.observation_data.loc[pst.nnz_obs_names[3],"weight"] = 100.0
    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options["ies_subset_size"] = 4
    pst.pestpp_options["ies_lambda_mults"] = [0.1,1.0,10.0]
    pst.control_data.noptmax = 3
    print("writing pst")
    pst.write(os.path.join(template_d, "pest_ineq.pst"))
    print("starting slaves")
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_ineq.pst", num_workers=10, master_dir=test_d,
                               worker_root=model_d,port=port)
    with open(os.path.join(test_d,"pest_ineq.rmr"),'r') as f:
        for line in f:
            print(line)

    obs_csvs = [f for f in os.listdir(test_d) if f.endswith("obs.csv")]
    print(obs_csvs)
    df = pd.read_csv(os.path.join(test_d,obs_csvs[-1]),index_col=0)
    df.columns = df.columns.map(str.lower)
    pst = pyemu.Pst(os.path.join(template_d, "pest_ineq.pst"))

    axes = [plt.subplot(pst.nnz_obs,1,i+1) for i in range(pst.nnz_obs)]
    #df.loc[:,pst.nnz_obs_names].hist(bins=10,axes=axes)
    for i,n in enumerate(pst.nnz_obs_names):
        v = pst.observation_data.loc[n,"obsval"]
        ax = axes[i]
        df.loc[:,n].hist(ax=ax,bins=10)
        ax.plot([v,v],ax.get_ylim(),"k--",lw=2.0)
    #plt.show()


def tenpar_fixed_test2():
    """tenpar fixed test 2"""
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_fixed2")
    template_d = os.path.join(model_d, "template")
    pst = pyemu.Pst(os.path.join(template_d,"pest.pst"))
    pst.control_data.noptmax = 1
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d,test_d)
    pst.parameter_data.loc[:,"partrans"] = "log"
    pst.parameter_data.loc["k_01","partrans"] = "fixed"
    cov = pyemu.Cov.from_parameter_data(pst)
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst,cov=cov,num_reals=10)
    #pe = pe.loc[:,pst.adj_par_names[3:5]]
    pe.loc[:,"stage"] = np.linspace(0.0,1.0,pe.shape[0])
    pe.to_csv(os.path.join(template_d,"par_fixed.csv"))
    fixed_pars = ["stage","k_01"]
    pst.parameter_data.loc[fixed_pars,"partrans"] = "fixed"
    pst.pestpp_options["ies_par_en"] = "par_fixed.csv"
    pst.pestpp_options["ies_lambda_mults"] = 1.0
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["ies_save_binary"] = False
    pst.pestpp_options["ies_include_base"] = False

    pst.write(os.path.join(template_d, "pest_fixed.pst"))

    # pyemu.helpers.run("{0} pest.pst".format(exe_path), cwd=test_d)
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_fixed.pst", num_workers=5, master_dir=test_d,
                                worker_root=model_d, port=port)
    df = pd.read_csv(os.path.join(test_d,"pest_fixed.{0}.par.csv".format(pst.control_data.noptmax)),index_col=0)
    df.columns = df.columns.map(str.lower)
    #df = df.iloc[:-1, :]
    #df.index = pe.index
    print(df.loc[:,"k_01"])
    print(df.loc[:,"stage"] - pe.loc[:,"stage"])
    assert df.loc[:,"k_01"].mean() == pst.parameter_data.loc["k_01","parval1"]
    assert np.abs(df.loc[:,"stage"] - pe.loc[:,"stage"]).max() < 1.0e-5

    pst.control_data.noptmax = -2
    par = pst.parameter_data
    par.loc["stage","parval1"] = 2
    pst.write(os.path.join(template_d,"pest_fixed.pst"))
    pyemu.os_utils.run("{0} pest_fixed.pst".format(exe_path),cwd=template_d)

    # pe = pe.loc[:,pst.adj_par_names[3:]]
    # pe.to_csv(os.path.join(test_d, "par_fixed.csv"))
    # if "win" in platform.platform().lower(): #bc of the stupid popup
    #     return
    # try:
    #     pyemu.os_utils.run("{0} {1}".format(exe_path,"pest_fixed.pst"),cwd=test_d)
    # except:
    #     pass
    # else:
    #     raise Exception()


def tenpar_fixed_test3():
    """tenpar fixed test 2"""
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_fixed2")
    template_d = os.path.join(model_d, "template")
    pst = pyemu.Pst(os.path.join(template_d,"pest.pst"))
    pst.control_data.noptmax = 2
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d,test_d)
    pst.parameter_data.loc[:,"partrans"] = "log"
    
    cov = pyemu.Cov.from_parameter_data(pst)
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst,cov=cov,num_reals=10) 
    pst.parameter_data.loc["k_01","partrans"] = "fixed"
    pe = pe.loc[:,pst.adj_par_names] 
    pe.to_csv(os.path.join(test_d,"par_fixed.csv"))
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options.pop("ies_par_en",None)
    pst.pestpp_options["ies_par_en"] = "par_fixed.csv"

    pst.parameter_data.loc["k_01","partrans"] = "fixed"
    pst.parameter_data.loc["k_01","scale"] = 10.0
    pst.pestpp_options["ies_lambda_mults"] = 1.0
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["ies_save_binary"] = False
    pst.pestpp_options["ies_include_base"] = False

    pst.control_data.noptmax = 0
    pst.write(os.path.join(test_d, "pest_fixed.pst"))
    pyemu.helpers.run("{0} pest_fixed.pst".format(exe_path), cwd=test_d)
    ins_name = os.path.join(test_d,"hk_Layer_1.ref.ins")
    with open(ins_name,'w') as f:
        f.write("pif ~\nl1 ")
        for i in range(9):
            f.write(" !k_{0:02d}! w ".format(i+1))
        f.write(" !k_10! \n")
    df = pst.add_observations(ins_name,ins_name.replace(".ins",""),pst_path=".")
    pst.observation_data.loc[df.obsnme,"weight"] = 0.0
    pst.write(os.path.join(test_d, "pest_fixed.pst"))
    pyemu.helpers.run("{0} pest_fixed.pst".format(exe_path), cwd=test_d)
    
    pst.control_data.noptmax = 4
    pst.write(os.path.join(test_d, "pest_fixed.pst"))
    pyemu.helpers.run("{0} pest_fixed.pst".format(exe_path), cwd=test_d)
    # pyemu.os_utils.start_workers(template_d, exe_path, "pest_fixed.pst", num_workers=5, master_dir=test_d,
    #                             worker_root=model_d, port=port)
    prev = None
    for itr in range(pst.control_data.noptmax):
        pe = pd.read_csv(os.path.join(test_d,"pest_fixed.{0}.par.csv".format(itr)),index_col=0)
        oe = pd.read_csv(os.path.join(test_d,"pest_fixed.{0}.obs.csv".format(itr)),index_col=0)
        
        pd.columns = pe.columns.map(str.lower)
        #df = df.iloc[:-1, :]
        #df.index = pe.index
        print(itr,pe.loc[:,"k_01"].values,oe.loc[:,"k_01"].values)
        d = pe.loc[:,"k_01"].values - oe.loc[:,"k_01"].values
        assert np.abs(d.max()) < 1e-6
        if prev is not None:
            d = prev - oe.loc[:,"k_01"].values
        assert np.abs(d.max()) < 1e-6

        prev = oe.loc[:,"k_01"].values

    test_d = os.path.join(model_d, "master_fixed2")
    template_d = os.path.join(model_d, "template")
    pst = pyemu.Pst(os.path.join(template_d,"pest.pst"))
    pst.control_data.noptmax = 2
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d,test_d)
    pst.parameter_data.loc[:,"partrans"] = "log"
    
    cov = pyemu.Cov.from_parameter_data(pst)
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst,cov=cov,num_reals=10) 
    #pst.parameter_data.loc["k_01","partrans"] = "fixed"
    pe = pe.loc[:,pst.adj_par_names] 
    pe.to_csv(os.path.join(test_d,"par_fixed.csv"))
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options.pop("ies_par_en",None)
    pst.pestpp_options["ies_par_en"] = "par_fixed.csv"

    pst.parameter_data.loc["k_01","partrans"] = "fixed"
    pst.parameter_data.loc["k_01","scale"] = 10.0
    pst.pestpp_options["ies_lambda_mults"] = 1.0
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["ies_save_binary"] = False
    pst.pestpp_options["ies_include_base"] = False

    pst.control_data.noptmax = 0
    pst.write(os.path.join(test_d, "pest_fixed.pst"))
    pyemu.helpers.run("{0} pest_fixed.pst".format(exe_path), cwd=test_d)
    ins_name = os.path.join(test_d,"hk_Layer_1.ref.ins")
    with open(ins_name,'w') as f:
        f.write("pif ~\nl1 ")
        for i in range(9):
            f.write(" !k_{0:02d}! w ".format(i+1))
        f.write(" !k_10! \n")
    df = pst.add_observations(ins_name,ins_name.replace(".ins",""),pst_path=".")
    pst.observation_data.loc[df.obsnme,"weight"] = 0.0
    pst.write(os.path.join(test_d, "pest_fixed.pst"))
    pyemu.helpers.run("{0} pest_fixed.pst".format(exe_path), cwd=test_d)
    
    pst.control_data.noptmax = 4
    pst.write(os.path.join(test_d, "pest_fixed.pst"))
    pyemu.helpers.run("{0} pest_fixed.pst".format(exe_path), cwd=test_d)
    # pyemu.os_utils.start_workers(template_d, exe_path, "pest_fixed.pst", num_workers=5, master_dir=test_d,
    #                             worker_root=model_d, port=port)
    prev = None
    for itr in range(pst.control_data.noptmax):
        pe = pd.read_csv(os.path.join(test_d,"pest_fixed.{0}.par.csv".format(itr)),index_col=0)
        oe = pd.read_csv(os.path.join(test_d,"pest_fixed.{0}.obs.csv".format(itr)),index_col=0)
        
        pd.columns = pe.columns.map(str.lower)
        #df = df.iloc[:-1, :]
        #df.index = pe.index
        print(itr,pe.loc[:,"k_01"].values,oe.loc[:,"k_01"].values)
        d = pe.loc[:,"k_01"].values - oe.loc[:,"k_01"].values
        assert np.abs(d.max()) < 1e-6
        if prev is not None:
            d = prev - oe.loc[:,"k_01"].values
        assert np.abs(d.max()) < 1e-6

        prev = oe.loc[:,"k_01"].values
    

    # pst.control_data.noptmax = -2
    # par = pst.parameter_data
    # par.loc["stage","parval1"] = 2
    # pst.write(os.path.join(template_d,"pest_fixed.pst"))
    # pyemu.os_utils.run("{0} pest_fixed.pst".format(exe_path),cwd=template_d)


def tenpar_fixed_test():
    """tenpar fixed test"""
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_fixed")
    template_d = os.path.join(model_d, "template")
    pst = pyemu.Pst(os.path.join(template_d,"pest.pst"))
    pst.parameter_data.loc[:, "partrans"] = "log"
    pst.control_data.noptmax = 1
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d,test_d)
    cov = pyemu.Cov.from_parameter_data(pst)
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst,cov=cov,num_reals=10)

    pe.loc[:,"stage"] = np.linspace(0.0,1.0,pe.shape[0])
    pe.loc[:,"k_01"] = 5.0
    pe.to_csv(os.path.join(template_d,"par_fixed.csv"))
    pe.to_binary(os.path.join(template_d, "par_fixed.jcb"))
    fixed_pars = ["stage","k_01"]
    pst.parameter_data.loc[fixed_pars,"partrans"] = "fixed"

    def compare():
        csvs = [f for f in os.listdir(test_d) if f.endswith(".par.csv") and "pest_fixed" in f and "mean" not in f and "base" not in f]
        print(csvs)
        dfs = [pd.read_csv(os.path.join(test_d,csv),index_col=0) for csv in csvs]
        for df in dfs:
            df.columns = df.columns.map(str.lower)
            df = df.loc[:,fixed_pars]
            #df = df.iloc[:-1,:]
            #df.index = pe.index
            diff = pe.loc[df.index,fixed_pars] - df
            assert diff.apply(lambda x: np.abs(x)).sum().sum() < 0.01, diff

    pst.pestpp_options["ies_par_en"] = "par_fixed.csv"
    pst.pestpp_options["ies_lambda_mults"] = 1.0
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["ies_save_binary"] = False
    pst.pestpp_options["ies_include_base"] = False
    pst.write(os.path.join(template_d, "pest_fixed.pst"))
    #pyemu.helpers.run("{0} pest.pst".format(exe_path), cwd=test_d)
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_fixed.pst", num_workers=5, master_dir=test_d,
                               worker_root=model_d,port=port)
    compare()
    pe.to_binary(os.path.join(template_d,"par_fixed.jcb"))
    print(pe)
    pst.pestpp_options["ies_par_en"] = "par_fixed.jcb"
    pst.write(os.path.join(template_d, "pest.pst"))
    #pyemu.helpers.run("{0} pest.pst".format(exe_path), cwd=test_d)
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_fixed.pst", num_workers=5, master_dir=test_d,
                               worker_root=model_d,port=port)
    compare()

    pst.pestpp_options["ies_par_en"] = "par_fixed.jcb"
    pst.pestpp_options["ies_save_binary"] = 'true'
    pst.write(os.path.join(template_d, "pest_fixed.pst"))
    #pyemu.helpers.run("{0} pest.pst".format(exe_path), cwd=test_d)
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_fixed.pst", num_workers=5, master_dir=test_d,
                               worker_root=model_d,port=port)
    pe1 = pyemu.ParameterEnsemble.from_binary(pst=pst,filename=os.path.join(test_d,"pest_fixed.0.par.jcb"))
    #pe1.index = pe.index
    diff = pe - pe1
    assert np.abs(diff.apply(lambda x: np.abs(x)).sum().sum()) < 1.0e-6,diff

# jwhite - disable for now until more testing in the weights stuff
# def tenpar_weights_test():
#     """tenpar weights test"""
#     model_d = "ies_10par_xsec"
#     test_d = os.path.join(model_d, "master_weights")
#     template_d = os.path.join(model_d, "template")
#     pst = pyemu.Pst(os.path.join(template_d,"pest.pst"))
    
#     if os.path.exists(test_d):
#         shutil.rmtree(test_d)
#     shutil.copytree(template_d,test_d)

#     dfs = []
    

#     for i in range(3):
#         obs = pst.observation_data.weight.copy()    
#         dfs.append(obs)

#     df = pd.concat(dfs,axis=1).T
#     df.index = np.arange(df.shape[0])
#     df.to_csv(os.path.join(test_d,"weights.csv"))
#     oe = pyemu.ObservationEnsemble.from_id_gaussian_draw(pst=pst,num_reals=df.shape[0])
#     oe.to_csv(os.path.join(test_d,"obs.csv"))
#     pst.control_data.noptmax = -1
#     #pst.pestpp_options["ies_weights_ensemble"] = "weights.csv"
#     pst.pestpp_options["ies_num_reals"] = df.shape[0]
#     pst.pestpp_options["ies_lambda_mults"] = 1.0
#     pst.pestpp_options["lambda_scale_fac"] = 1.0
#     pst.pestpp_options["ies_obs_en"] = "obs.csv"

#     pst.write(os.path.join(test_d,"pest.pst"))
#     pyemu.helpers.run("{0} pest.pst".format(exe_path), cwd=test_d)
#     df_act = pd.read_csv(os.path.join(test_d,"pest.phi.actual.csv"))
#     df_meas = pd.read_csv(os.path.join(test_d,"pest.phi.meas.csv"))

#     pst.pestpp_options["ies_weights_ensemble"] = "weights.csv"
#     pst.write(os.path.join(test_d,"pest.pst"))
#     pyemu.helpers.run("{0} pest.pst".format(exe_path), cwd=test_d)
#     df_act1 = pd.read_csv(os.path.join(test_d,"pest.phi.actual.csv"))
#     df_meas1 = pd.read_csv(os.path.join(test_d,"pest.phi.meas.csv"))
#     print(df_act.loc[0,"mean"],df_act1.loc[0,"mean"])
#     assert df_act.loc[0,"mean"] == df_act1.loc[0,"mean"]

#     print(df_meas.loc[0,"mean"],df_meas1.loc[0,"mean"])
#     assert df_meas.loc[0,"mean"] == df_meas1.loc[0,"mean"]


def tenpar_tight_tol_test():
    """tenpar tight tol test"""
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_tighttol")
    template_d = os.path.join(model_d, "template")
    pst = pyemu.Pst(os.path.join(template_d,"pest.pst"))
    
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d,test_d)

    
    pst.control_data.noptmax = 3
    #pst.pestpp_options["ies_weights_ensemble"] = "weights.csv"
    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options["ies_lambda_mults"] = 0.00000001
    pst.pestpp_options["lambda_scale_fac"] = 10.0
    pst.pestpp_options["ies_initial_lambda"] = 0.000001
    pst.pestpp_options["ies_accept_phi_fac"] = 1.0

    pst.write(os.path.join(test_d,"pest.pst"))
    pyemu.helpers.run("{0} pest.pst".format(exe_path), cwd=test_d)
    

def tenpar_weight_pareto():

    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_weight_pareto")
    template_d = os.path.join(model_d, "template")
    pst = pyemu.Pst(os.path.join(template_d,"pest.pst"))
    
    #if os.path.exists(test_d):
    #   shutil.rmtree(test_d)
    #shutil.copytree(template_d,test_d)

    dfs = []
    obs = pst.observation_data
    obs.loc["h02_08", "obsval"] = 0
    obs.loc["h02_08", "weight"] = 1.0
    dfs = []
    weights = np.linspace(0.0,10.0,50)
    for weight in weights:
        obs = pst.observation_data.weight.copy()
        obs.loc["h02_08"] = weight
        dfs.append(obs)
        # for i in range(nreal_per):
        #     dfs.append(obs)
        #for weight2 in weights:
        #    obs1 = obs.copy()
        #    obs1.loc["h02_08"] = weight2
        #    #print(obs1)
        #    dfs.append(obs1)

    df = pd.concat(dfs,axis=1).T
    df.index = np.arange(df.shape[0])
    df.to_csv(os.path.join(template_d,"weights.csv"))

    obs = pst.observation_data.obsval
    df = pd.concat([obs.copy() for i in range(df.shape[0])],axis=1).T
    df.index = np.arange(df.shape[0])
    df.to_csv(os.path.join(template_d, "obs.csv"))


    pst.control_data.noptmax = 3
    pst.pestpp_options["ies_weights_ensemble"] = "weights.csv"
    pst.pestpp_options["ies_obs_en"] = "obs.csv"
    pst.pestpp_options["ies_num_reals"] = df.shape[0]
    pst.pestpp_options["ies_subset_size"] = df.shape[0]
    pst.pestpp_options["ies_save_binary"]= False
    #pst.pestpp_options["ies_lambda_mults"] = 1.0
    #pst.pestpp_options["lambda_scale_fac"] = 1.0
    obs = pst.observation_data
    #obs.loc[pst.nnz_obs_names,"obsval"] = obs.loc[pst.nnz_obs_names,"obsval"].mean()
    #obs.loc["h01_04","obsval"] = 0.0
    #obs.loc["h01_06","obsval"] = 6.0

    pst.write(os.path.join(template_d,"pest_pareto.pst"))
    pyemu.os_utils.start_workers(template_d,exe_path,"pest_pareto.pst",num_workers=40,
                                worker_root=model_d,master_dir=test_d,port=port)
    obs = pst.observation_data
    df_init = pd.read_csv(os.path.join(test_d,"pest_pareto.0.obs.csv".format(pst.control_data.noptmax)))
    df = pd.read_csv(os.path.join(test_d,"pest_pareto.{0}.obs.csv".format(pst.control_data.noptmax)))
    df_phi = pd.read_csv(os.path.join(test_d,"pest_pareto.phi.meas.csv"))
    df_phi = df_phi.loc[:,[str(v) for v in df_init.index.values]]
    # fig = plt.figure(figsize=(10,10))
    # ax = plt.subplot(221)
    # ax2 = plt.subplot(223)
    # ax3 = plt.subplot(224)
    # ax2.scatter(((df_init.H01_04-obs.loc["h01_04","obsval"])**2),
    #     ((df_init.H01_06-obs.loc["h01_06","obsval"])**2),s=4, color='0.5')
    # ax2.scatter(((df.H01_04-obs.loc["h01_04","obsval"])**2),
    #     ((df.H01_06-obs.loc["h01_06","obsval"])**2),s=4, color='b')
    # ax.scatter(df_phi.iloc[0,:],((df_init.H01_06-obs.loc["h01_06","obsval"])**2),s=4, color='0.5')
    # ax.scatter(df_phi.iloc[-1,:],((df.H01_06-obs.loc["h01_06","obsval"])**2),s=4, color='b')
    # ax3.scatter(df_phi.iloc[0,:],((df_init.H01_04-obs.loc["h01_04","obsval"])**2),s=4, color='0.5')
    # ax3.scatter(df_phi.iloc[-1,:],((df.H01_04-obs.loc["h01_04","obsval"])**2),s=4, color='b')

    # ax2.scatter(df_init.H01_04,df_init.H01_06, s=4, color='0.5')
    # ax2.scatter(df.H01_04,df.H01_06, s=4, color='b')
    # ax.scatter(df_phi.iloc[0, :], df_init.H01_06, s=4, color='0.5')
    # ax.scatter(df_phi.iloc[-1, :], df.H01_06, s=4, color='b')
    # ax3.scatter(df_phi.iloc[0, :], df_init.H01_04, s=4, color='0.5')
    # ax3.scatter(df_phi.iloc[-1, :], df.H01_04, s=4, color='b')

    # fig = plt.figure(figsize=(10, 10))
    # ax = plt.subplot(221)
    # ax.scatter(df_phi.iloc[0, :], df_init.H02_08, s=4, color='0.5')
    # ax.scatter(df_phi.iloc[-1, :], df.H02_08, s=4, color='b')
    # ax.scatter(df_phi.iloc[0, :], (df_init.H02_08-obs.obsval["h02_08"])**2, s=4, color='0.5')
    # ax.scatter(df_phi.iloc[-1, :], (df.H02_08-obs.obsval["h02_08"])**2, s=4, color='b')

    #plt.show()


def rosenbrock_function():
    par_df = pd.read_csv("par.dat",sep='\s+',index_col=0)
    tot = 0.0
    for i in range(par_df.shape[0]-1):
        tot += 100.0*(par_df.iloc[i+1] - par_df.iloc[i]**2)**2 + (1 - par_df.iloc[i])**2
    with open("obs.dat",'w') as f:
        f.write("obs {0:20.8E}".format(tot))

def setup_rosenbrock():

    npar = 2
    test_d = "ies_rosenbrock"
    if not os.path.exists(test_d):
        os.mkdir(test_d)

    template_d = os.path.join(test_d,"template")
    if os.path.exists(template_d):
        shutil.rmtree(template_d)
    os.mkdir(template_d)

    with open(os.path.join(template_d,"par.dat.tpl"),'w') as f:
        f.write("ptf ~\n")
        for i in range(npar):
            f.write("par{0:04d}  ~   par{0:04d}   ~\n")
    with open(os.path.join(template_d,"obs.dat.ins"),'w') as f:
        f.write("pif ~\n")
        f.write("l1 w !obs1!\n")

    bd = os.getcwd()


def tenpar_localizer_test1():
    """tenpar local 1"""
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_localizer_test11")
    template_d = os.path.join(model_d, "test_template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)

    #mat = pyemu.Matrix.from_names(pst.nnz_obs_names,pst.adj_par_names).to_dataframe()
    mat = pyemu.Matrix.from_names(pst.nnz_obs_names,pst.adj_par_names).to_dataframe()
    mat.loc[:,:] = 1.0
    #mat.iloc[0,:] = 1
    mat = pyemu.Matrix.from_dataframe(mat)
    mat.to_ascii(os.path.join(template_d,"localizer.mat"))

    cov = pyemu.Cov.from_parameter_data(pst)
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst, cov=cov, num_reals=10)
    pe.enforce()
    pe.to_csv(os.path.join(template_d, "par_local.csv"))

    oe = pyemu.ObservationEnsemble.from_gaussian_draw(pst,num_reals=10)
    oe.to_csv(os.path.join(template_d,"obs_local.csv"))

    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options["ies_localizer"] = "localizer.mat"
    pst.pestpp_options["ies_lambda_mults"] = 1.0
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["ies_subset_size"] = 11
    pst.pestpp_options["ies_par_en"] = "par_local.csv"
    pst.pestpp_options["ies_obs_en"] = "obs_local.csv"
    pst.pestpp_options["ies_verbose_level"] = 1
    pst.control_data.noptmax = 3

    #pst.pestpp_options["ies_verbose_level"] = 3
    pst_name = os.path.join(template_d,"pest_local.pst")
    pst.write(pst_name)
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_local.pst", num_workers=10,
                                   master_dir=test_d, verbose=True, worker_root=model_d,
                                   port=port)
    phi_df1 = pd.read_csv(os.path.join(test_d,"pest_local.phi.meas.csv"))

    pst.pestpp_options.pop("ies_localizer")
    pst.write(os.path.join(template_d,"pest_base.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_base.pst", num_workers=10,
                                   master_dir=test_d+"_base", verbose=True, worker_root=model_d,
                                   port=port)
    phi_df2 = pd.read_csv(os.path.join(test_d+"_base","pest_base.phi.meas.csv"))
    diff = phi_df1 - phi_df2
    print(diff.max().max())

    plt.plot(phi_df1.total_runs,phi_df1.loc[:,"mean"], label="local")
    plt.plot(phi_df2.total_runs,phi_df2.loc[:,"mean"], label="full")
    plt.legend()
    plt.savefig(os.path.join(test_d,"local_test.pdf"))
    assert diff.max().max() == 0

def tenpar_localizer_test2():
    """tenpar local 2"""
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_localizer_test2")
    template_d = os.path.join(model_d, "test_template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    #shutil.copytree(template_d, test_d)
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)
    pst.parameter_data.loc[:,"partrans"] = "log"
    cov = pyemu.Cov.from_parameter_data(pst)
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst,cov=cov,num_reals=10)
    pe.enforce()
    pe.to_csv(os.path.join(template_d,"par_local.csv"))

    cov = pyemu.Cov.from_observation_data(pst)
    oe = pyemu.ObservationEnsemble.from_gaussian_draw(pst=pst, num_reals=10)
    oe.to_csv(os.path.join(template_d, "obs_local.csv"))

    #mat = pyemu.Matrix.from_names(pst.nnz_obs_names,pst.adj_par_names).to_dataframe()
    use_pars = pst.adj_par_names[::2]
    pst.parameter_data.loc[:, "pargp"] = "fixed"
    pst.parameter_data.loc[use_pars, "pargp"] = "adj"

    #mat = pyemu.Matrix.from_names(["head"],use_pars).to_dataframe()
    #mat = pyemu.Matrix.from_names(["head"],["adj","fixed"]).to_dataframe()
    mat = pyemu.Matrix.from_names(pst.nnz_obs_names,pst.adj_par_names).to_dataframe()
    mat.loc[:,:] = 0.0
    mat.loc[:,use_pars] = 1.0
    #mat.iloc[0,:] = 1
    mat = pyemu.Matrix.from_dataframe(mat)
    mat.to_ascii(os.path.join(template_d,"localizer.mat"))

    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options["ies_localizer"] = "localizer.mat"
    pst.pestpp_options["ies_lambda_mults"] = 1.0
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["ies_subset_size"] = 11
    pst.pestpp_options["ies_par_en"] = "par_local.csv"
    pst.pestpp_options["ies_obs_en"] = "obs_local.csv"
    pst.pestpp_options["ies_include_base"] = False
    pst.pestpp_options["ies_verbose_level"] = 3

    pst.control_data.noptmax = 2
    
    #pst.pestpp_options["ies_verbose_level"] = 3
    pst_name = os.path.join(template_d,"pest_local.pst")
    pst.write(pst_name)

    pyemu.os_utils.start_workers(template_d, exe_path, "pest_local.pst", num_workers=11,
                                       master_dir=test_d, verbose=True, worker_root=model_d,
                                       port=port)
    phi_df1 = pd.read_csv(os.path.join(test_d,"pest_local.phi.actual.csv"))
    par_df1 = pd.read_csv(os.path.join(test_d,"pest_local.{0}.par.csv".format(pst.control_data.noptmax)),index_col=0)
    par_df1.columns = par_df1.columns.str.lower()

    pst.pestpp_options.pop("ies_localizer")
    pst.parameter_data.loc[:,"partrans"] = "fixed"

    pst.parameter_data.loc[use_pars,"partrans"] = "log"
    pst.write(os.path.join(template_d,"pest_base.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_base.pst", num_workers=11,
                                   master_dir=test_d+"_base", verbose=True, worker_root=model_d,
                                   port=port)
    phi_df2 = pd.read_csv(os.path.join(test_d+"_base","pest_base.phi.actual.csv"))
    par_df2 = pd.read_csv(os.path.join(test_d+"_base","pest_base.{0}.par.csv".format(pst.control_data.noptmax)),index_col=0)
    par_df2.columns = par_df2.columns.str.lower()
    # plt.close("all")
    # plt.plot(phi_df1.total_runs,phi_df1.loc[:,"mean"], label="local")
    # plt.plot(phi_df2.total_runs,phi_df2.loc[:,"mean"], label="full")
    # plt.legend()
    # plt.savefig(os.path.join(test_d,"local_test.pdf"))

    for par in par_df1.columns:
        diff = np.abs(par_df1.loc[:,par] - par_df2.loc[:,par])
        print(diff)
        print(diff.sum())
        assert diff.sum() < 1.0e-3, diff.sum()
        if par in use_pars:
            continue
        diff = np.abs(par_df1.loc[:,par] - pe.loc[:,par])
        print(diff.sum())
        assert diff.sum() < 1.0e-3, diff.sum()
    diff = (phi_df1 - phi_df2).apply(lambda x: np.abs(x))
    print(diff.max().max())
    assert diff.max().max() <  1.0e-3, diff.max().max()

    # for i in range(100):
    #     t = test_d + "_temp"
    #
    #     pyemu.os_utils.start_workers(template_d, exe_path, "pest_local.pst", num_workers=11,
    #                                    master_dir=t, verbose=True, worker_root=model_d,
    #                                    port=port)
    #     phi_df1 = pd.read_csv(os.path.join(t,"pest_local.phi.actual.csv"))
    #     par_df1 = pd.read_csv(os.path.join(t,"pest_local.{0}.par.csv".format(pst.control_data.noptmax)),index_col=0)
    #     par_df1.columns = par_df1.columns.str.lower()
    #     for par in par_df1.columns:
    #         diff = np.abs(par_df1.loc[:,par] - par_df2.loc[:,par])
    #         print(diff)
    #         print(diff.sum())
    #         assert diff.sum() < 1.0e-3, diff.sum()
    #         if par in use_pars:
    #             continue
    #         diff = np.abs(par_df1.loc[:,par] - pe.loc[:,par])
    #         print(diff.sum())
    #         assert diff.sum() < 1.0e-3, diff.sum()
    #     diff = (phi_df1 - phi_df2).apply(np.abs)
    #     print(diff.max().max())
    #     assert diff.max().max() <  1.0e-3, diff.max().max()
    #     #dfs.append(par_df1)
    #


def prep_for_travis(model_d):
    assert os.path.exists(os.path.join(model_d,"test_template"))
    s_d = "test_bin"
    d_d = os.path.join(model_d,"test_template")
    # if os.path.exists(d_d):
    #     shutil.rmtree(d_d)
    # shutil.copytree(s_d,d_d)
    for f in os.listdir(s_d):
        if os.path.exists(os.path.join(d_d,f)):
            os.remove(os.path.join(d_d,f))
        shutil.copy2(os.path.join(s_d,f),os.path.join(d_d,f))
    pst_file = os.path.join(model_d,"test_template","pest.pst")
    pst = pyemu.Pst(pst_file)
    pst.model_command = "./" + pst.model_command[0]
    print(pst.model_command)
    pst.write(pst_file)

def tenpar_incr_num_reals_test():
    """tenpar incr num reals test"""
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_incr_num_reals2")
    template_d = os.path.join(model_d, "template")
    pst = pyemu.Pst(os.path.join(template_d, "pest.pst"))
    num_reals = 10
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pst.control_data.noptmax = 2
    pst.pestpp_options = {"ies_num_reals":num_reals}
    pst.pestpp_options["ies_lambda_mults"] = 1.0
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["ies_num_reals"] = num_reals
    pst.pestpp_options["ies_include_base"] = True
    pst.write(os.path.join(test_d,"pest_restart.pst"))
    pyemu.os_utils.run("{0} {1}".format(exe_path, "pest_restart.pst"), cwd=test_d)

    pst.pestpp_options["ies_par_en"] = "pest_restart.2.par.csv"
    pst.pestpp_options["ies_obs_en"] = "pest_restart.obs+noise.csv"
    pst.pestpp_options["ies_restart_obs_en"] = "pest_restart.2.obs.csv"
    #pst.pestpp_options["ies_num_reals"] += 1
    pst.control_data.noptmax = 2
    pst.write(os.path.join(test_d,"pest.pst"))
    pyemu.os_utils.run("{0} {1}".format(exe_path, "pest.pst"), cwd=test_d)

    df1 = pd.read_csv(os.path.join(test_d,"pest_restart.2.par.csv"),index_col=0)
    df2 = pd.read_csv(os.path.join(test_d,"pest.0.par.csv"),index_col=0)
    diff = df1.loc["base",:] - df2.loc["base",:]
    print(diff)
    print(diff.max())

def tenpar_subset_how_test():
    """tenpar subet how"""
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_subset_how")
    template_d = os.path.join(model_d, "template")
    pst = pyemu.Pst(os.path.join(template_d, "pest.pst"))
    num_reals = 10
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pst.pestpp_options = {"ies_num_reals":num_reals}
    pst.control_data.noptmax = -1
    pst.write(os.path.join(template_d,"pest_restart.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_restart.pst", num_workers=10,
                                worker_root=model_d, master_dir=test_d, port=port)
    pyemu.os_utils.run("{0} {1}".format(exe_path, "pest_restart.pst"), cwd=test_d)
    pst.pestpp_options["ies_par_en"] = "pest_restart.0.par.csv"
    pst.pestpp_options["ies_lambda_mults"] = 1.0
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["ies_num_reals"] = num_reals
    pst.pestpp_options["ies_restart_obs_en"] = "pest_restart.0.obs.csv"
    pst.pestpp_options["ies_obs_en"] = "pest_restart.obs+noise.csv"
    pst.control_data.noptmax = 1
    means = []
    for how in ["first","last","random","phi_based"]:
        pst.pestpp_options["ies_subset_how"] = how
        pst_file = "pest_{0}.pst".format(how)
        pst.write(os.path.join(test_d,pst_file))
        pyemu.os_utils.run("{0} {1}".format(exe_path,pst_file),cwd=test_d)
        df = pd.read_csv(os.path.join(test_d,pst_file.replace(".pst",".phi.composite.csv")))
        means.append(df.iloc[-1].loc["mean"])

    means = np.array(means)
    print(means.mean(),means)
    assert means.mean() == means[0]

    means = []
    pst.pestpp_options["ies_subset_size"] = 2
    for how in ["first", "last", "random", "phi_based"]:
        pst.pestpp_options["ies_subset_how"] = how
        pst_file = "pest_{0}.pst".format(how)
        pst.write(os.path.join(test_d, pst_file))
        pyemu.os_utils.run("{0} {1}".format(exe_path, pst_file), cwd=test_d)
        df = pd.read_csv(os.path.join(test_d, pst_file.replace(".pst", ".phi.composite.csv")))
        means.append(df.iloc[-1].loc["mean"])

    means = np.array(means)
    print(means.mean(), means)
    assert means.mean() == means[0]

    means = []
    pst.pestpp_options["ies_subset_size"] = num_reals
    for how in ["first", "last", "random", "phi_based"]:
        pst.pestpp_options["ies_subset_how"] = how
        pst_file = "pest_{0}.pst".format(how)
        pst.write(os.path.join(test_d, pst_file))
        pyemu.os_utils.run("{0} {1}".format(exe_path, pst_file), cwd=test_d)
        df = pd.read_csv(os.path.join(test_d, pst_file.replace(".pst", ".phi.composite.csv")))
        means.append(df.iloc[-1].loc["mean"])

    means = np.array(means)
    print(means.mean(), means)
    assert means.mean() == means[0]

    #for how,df in obs_dfs.items():


def tenpar_localizer_test3():
    """tenpar local 3"""
    plt.close("all")
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_localizer_test3")
    template_d = os.path.join(model_d, "test_template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    # shutil.copytree(template_d, test_d)
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)
    pst.parameter_data.loc[:,"partrans"] = "log"
    cov = pyemu.Cov.from_parameter_data(pst)
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst, cov=cov, num_reals=10)
    pe.enforce()
    pe.to_csv(os.path.join(template_d, "par_local.csv"))

    cov = pyemu.Cov.from_observation_data(pst)
    oe = pyemu.ObservationEnsemble.from_gaussian_draw(pst=pst, num_reals=10)
    oe.to_csv(os.path.join(template_d, "obs_local.csv"))


    mat = pyemu.Matrix.from_names(pst.nnz_obs_names,pst.adj_par_names).to_dataframe()
    use_pars = pst.adj_par_names[0:2]


    mat.loc[:, :] = 0.0
    mat.loc[pst.nnz_obs_names, use_pars] = 1.0
    # mat.iloc[0,:] = 1
    mat = pyemu.Matrix.from_dataframe(mat)
    mat.to_ascii(os.path.join(template_d, "localizer.mat"))

    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options["ies_localizer"] = "localizer.mat"
    pst.pestpp_options["ies_lambda_mults"] = 1.0
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["ies_subset_size"] = 11
    pst.pestpp_options["ies_par_en"] = "par_local.csv"
    pst.pestpp_options["ies_include_base"] = False
    pst.pestpp_options["ies_accept_phi_fac"] = 1.1
    pst.pestpp_options["ies_obs_en"] = "obs_local.csv"
    pst.pestpp_options["ies_localize_how"] = "pars"
    pst.control_data.nphistp = 10
    pst.control_data.nphinored = 10
    pst.control_data.noptmax = 3

    # pst.pestpp_options["ies_verbose_level"] = 3
    pst_name = os.path.join(template_d, "pest_local.pst")
    pst.write(pst_name)
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_local.pst", num_workers=11,
                               master_dir=test_d, verbose=True, worker_root=model_d,
                               port=port)
    phi_df1 = pd.read_csv(os.path.join(test_d, "pest_local.phi.composite.csv"))
    par_df1 = pd.read_csv(os.path.join(test_d, "pest_local.{0}.par.csv".format(pst.control_data.noptmax)), index_col=0)
    par_df1.columns = par_df1.columns.str.lower()

    pst.pestpp_options.pop("ies_localizer")
    pst.parameter_data.loc[:, "partrans"] = "fixed"
    pst.parameter_data.loc[use_pars, "partrans"] = "log"
    pst.write(os.path.join(template_d, "pest_base.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_base.pst", num_workers=11,
                               master_dir=test_d + "_base", verbose=True, worker_root=model_d,
                               port=port)
    phi_df2 = pd.read_csv(os.path.join(test_d + "_base", "pest_base.phi.composite.csv"))
    par_df2 = pd.read_csv(os.path.join(test_d + "_base", "pest_base.{0}.par.csv".format(pst.control_data.noptmax)),
                          index_col=0)
    par_df2.columns = par_df2.columns.str.lower()
    plt.plot(phi_df1.total_runs, phi_df1.loc[:, "mean"], label="local")
    plt.plot(phi_df2.total_runs, phi_df2.loc[:, "mean"], label="full")
    plt.legend()
    plt.savefig(os.path.join(test_d, "local_test.pdf"))

    for par in par_df1.columns:
        diff = par_df1.loc[:, par] - par_df2.loc[:, par]
        print(diff)
        print(diff.sum())
        diff = (par_df1.loc[:, par] - pe.loc[:, par]).apply(lambda x: np.abs(x))
        print(diff.sum())
        #assert diff.max() < 1e-4
    diff = np.abs(phi_df1.loc[:,"mean"] - phi_df2.loc[:,"mean"])
    print(diff.max())
    #assert diff.max().max() < 0.5
  

def tenpar_restart_similar_test():
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_similar_base")
    template_d = os.path.join(model_d, "template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pst_name = os.path.join(test_d, "pest.pst")
    pst = pyemu.Pst(pst_name)
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst,pyemu.Cov.from_parameter_data(pst),20)
    pe = pe.loc[:,pst.par_names[::-1]]
    pe.index = ["real_{0}".format(i+1) for i in range(pe.shape[0])]
    pe.enforce()
    pe.to_binary(os.path.join(test_d,"prior.jcb"))
    pe.to_csv('org_test.csv')
    #pst.pestpp_options = {}
    #pst.pestpp_options["ies_par_en"] = "prior.jcb"
    #pe = pyemu.ParameterEnsemble.from_binary(pst,os.path.join(template_d,"ies_prior.jcb"))
    #pe = pe.loc[:,pst.par_names]
    #pe.to_binary(os.path.join(test_d,"test_prior.jcb"))
    pst.pestpp_options["ies_par_en"] = "prior.jcb"
    pst.pestpp_options["ies_num_threads"] = 1
    pst.pestpp_options["ies_num_reals"] = 5
    pst.pestpp_options["ies_lambda_mults"] = 1.0
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["ies_verbose_level"] = 1
    #pst.pestpp_options.pop("ies_localizer")
    pst.pestpp_options["ies_autoadaloc"] = False
    pst.pestpp_options["ies_save_lambda_en"] = True
    pst.control_data.noptmax = 1
    pst.write(os.path.join(test_d,"pest.pst"))
    pyemu.os_utils.run("{0} pest.pst".format(exe_path),cwd=test_d)

    phi1_df = pd.read_csv(os.path.join(test_d,"pest.phi.actual.csv"))

    shutil.copy2(os.path.join(test_d,"pest.0.par.csv"),os.path.join(template_d,"par.csv"))
    shutil.copy2(os.path.join(test_d,"pest.0.obs.csv"),os.path.join(template_d,"obs.csv"))
    shutil.copy2(os.path.join(test_d,"pest.obs+noise.csv"),os.path.join(template_d,"noise.csv"))
    pst.pestpp_options["ies_par_en"] = "par.csv"
    pst.pestpp_options["ies_obs_en"] = "noise.csv"
    pst.pestpp_options["ies_restart_obs_en"] = "obs.csv"
    
    
    test_d = os.path.join(model_d, "master_similar_restart")
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pst.write(os.path.join(test_d,"pest.pst"))
    pyemu.os_utils.run("{0} pest.pst".format(exe_path),cwd=test_d)

    phi2_df = pd.read_csv(os.path.join(test_d, "pest.phi.actual.csv"))

    d = (phi1_df.iloc[:,2:] - phi2_df.iloc[:,2:]).apply(lambda x: np.abs(x))
    print(d)
    print(d.max())
    print(d.max().max())
    assert d.max().max() < 1.0e-6


def tenpar_restart_similar_test2():
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_similar_base2")
    template_d = os.path.join(model_d, "template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pst_name = os.path.join(test_d, "pest.pst")
    pst = pyemu.Pst(pst_name)
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst,pyemu.Cov.from_parameter_data(pst),20)
    pe = pe.loc[:,pst.par_names[::-1]]
    pe.index = ["real_{0}".format(i+1) for i in range(pe.shape[0])]
    pe.enforce()
    pe._df.index = ["real_{0}".format(i+1) for i in range(pe.shape[0])]
    pe.to_csv(os.path.join(test_d,"prior.csv"))
    pe.to_csv('org_test.csv')
    pst.pestpp_options["ies_par_en"] = "prior.csv"
    pst.pestpp_options["ies_num_threads"] = 1
    pst.pestpp_options["ies_num_reals"] = 5
    pst.pestpp_options["ies_lambda_mults"] = 1.0
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["ies_verbose_level"] = 1
    pst.pestpp_options["ies_autoadaloc"] = False
    pst.pestpp_options["ies_save_lambda_en"] = True
    pst.pestpp_options["ies_include_base"] = True
    pst.control_data.noptmax = 1
    pst.write(os.path.join(test_d,"pest.pst"))
    pyemu.os_utils.run("{0} pest.pst".format(exe_path),cwd=test_d)

    noise = pd.read_csv(os.path.join(test_d,"pest.obs+noise.csv"),index_col=0)
    oe = pd.read_csv(os.path.join(test_d,"pest.0.obs.csv"),index_col=0)
    pe = pd.read_csv(os.path.join(test_d,"pest.0.par.csv"),index_col=0)
    print(pe.shape,oe.shape)
    assert oe.shape[0] == pe.shape[0],"{0},{1}".format(oe.shape,pe.shape)
    assert noise.shape[0] == pe.shape[0],"{0},{1}".format(noise.shape,pe.shape)


def tenpar_localizer_pdc_test():
    """tenpar local test with  pdc causing some pars to be completely localized out
    """
    model_d = "ies_10par_xsec"
    template_d = os.path.join(model_d, "test_template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name) 
    #spike one of the obs...
    pst.observation_data.loc[pst.nnz_obs_names[0],"obsval"] += 100
    cov = pyemu.Cov.from_parameter_data(pst)
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst, cov=cov, num_reals=10)
    pe.enforce()

    oe = pyemu.ObservationEnsemble.from_gaussian_draw(pst,num_reals=10)
    final_phis = []
    for out_ext in [".mat",".csv",".jcb"]:       
        test_d = os.path.join(model_d, "master_localizer_pdc_test_{0}".format(out_ext.replace('.',''))) 
        if os.path.exists(test_d):
            shutil.rmtree(test_d)
        #shutil.copytree(template_d, test_d)
        
        pe.to_csv(os.path.join(template_d, "par_local.csv"))
        oe.to_csv(os.path.join(template_d,"obs_local.csv"))

        pst_name = os.path.join(template_d, "pest.pst")
        pst = pyemu.Pst(pst_name) 
        #spike one of the obs...
        pst.observation_data.loc[pst.nnz_obs_names[0],"obsval"] += 100

        

        #mat = pyemu.Matrix.from_names(pst.nnz_obs_names,pst.adj_par_names).to_dataframe()
        mat = pyemu.Matrix.from_names(pst.nnz_obs_names,pst.adj_par_names).to_dataframe()
        mat.loc[:,:] = 0.0
        mat.loc[pst.nnz_obs_names[0],pst.adj_par_names[:2]] = 1.0
        mat.loc[pst.nnz_obs_names[1],pst.adj_par_names[2:]] = 1.0
        
        #mat.iloc[0,:] = 1
        loc_name = "localizer" + out_ext
        
        if out_ext == ".mat":
            mat = pyemu.Matrix.from_dataframe(mat)
            mat.to_ascii(os.path.join(template_d,loc_name))
        elif out_ext == ".csv":
            mat.to_csv(os.path.join(template_d,loc_name))
        elif out_ext == ".jcb":
            mat = pyemu.Matrix.from_dataframe(mat)
            mat.to_binary(os.path.join(template_d,loc_name))

        
        pst.pestpp_options = {}
        pst.pestpp_options["ies_num_reals"] = 10
        pst.pestpp_options["ies_localizer"] = loc_name
        pst.pestpp_options["ies_lambda_mults"] = 1.0
        pst.pestpp_options["lambda_scale_fac"] = 1.0
        pst.pestpp_options["ies_subset_size"] = 11
        pst.pestpp_options["ies_par_en"] = "par_local.csv"
        pst.pestpp_options["ies_obs_en"] = "obs_local.csv"
        pst.pestpp_options["ies_drop_conflicts"] = True
        pst.pestpp_options["ies_verbose_level"] = 4
        pst.pestpp_options["ies_localizer_forgive_missing"] = False
        pst.control_data.noptmax = 2

        #pst.pestpp_options["ies_verbose_level"] = 3
        pst_name = os.path.join(template_d,"pest_local_pdc.pst")
        pst.write(pst_name)
        pyemu.os_utils.start_workers(template_d, exe_path, "pest_local_pdc.pst", num_workers=10,
                                       master_dir=test_d, verbose=True, worker_root=model_d,
                                       port=port)
        phi_df1 = pd.read_csv(os.path.join(test_d,"pest_local_pdc.phi.meas.csv"))
        assert phi_df1.shape[0] == pst.control_data.noptmax+1
        assert phi_df1.loc[phi_df1.index[-1],"mean"] < phi_df1.loc[phi_df1.index[0],"mean"]
        mat = pyemu.Matrix.from_ascii(os.path.join(test_d, "initialized_localizer.mat"))
        assert mat.shape == (pst.nnz_obs-1, pst.npar_adj), str(mat.shape)

        # now with a group-based localizer
        test_d = os.path.join(model_d, "master_localizer_pdc_test2_{0}".format(out_ext.replace('.','')))
        if os.path.exists(test_d):
            shutil.rmtree(test_d)
        par = pst.parameter_data
        par.loc[pst.adj_par_names[:4],"pargp"] = "group1"
        par.loc[pst.adj_par_names[4:],"pargp"] = "group2"

        mat = pyemu.Matrix.from_names(pst.nnz_obs_names,["group1","group2"]).to_dataframe()
        mat.loc[:,:] = 0.0
        mat.loc[pst.nnz_obs_names[0],"group1"] = 1.0
        mat.loc[pst.nnz_obs_names[1],"group2"] = 1.0
        
        if out_ext == ".mat":
            mat = pyemu.Matrix.from_dataframe(mat)
            mat.to_ascii(os.path.join(template_d,loc_name))
        elif out_ext == ".csv":
            mat.to_csv(os.path.join(template_d,loc_name))
        elif out_ext == ".jcb":
            mat = pyemu.Matrix.from_dataframe(mat)
            mat.to_binary(os.path.join(template_d,loc_name))
        pst.write(pst_name)
        pyemu.os_utils.start_workers(template_d, exe_path, "pest_local_pdc.pst", num_workers=10,
                                       master_dir=test_d, verbose=True, worker_root=model_d,
                                       port=port)
        phi_df1 = pd.read_csv(os.path.join(test_d,"pest_local_pdc.phi.meas.csv"))
        assert phi_df1.shape[0] == pst.control_data.noptmax+1
        assert phi_df1.loc[phi_df1.index[-1],"mean"] < phi_df1.loc[phi_df1.index[0],"mean"]
        final_phi1 = phi_df1.loc[phi_df1.index[-1],"mean"]
        mat = pyemu.Matrix.from_ascii(os.path.join(test_d, "initialized_localizer.mat"))
        assert mat.shape == (pst.nnz_obs - 1, 2), str(mat.shape)

        #now restart

        shutil.copy2(os.path.join(test_d,"pest_local_pdc.0.obs.csv"),os.path.join(template_d,"restart_local_obs.csv"))
        shutil.copy2(os.path.join(test_d,"pest_local_pdc.0.par.csv"),os.path.join(template_d,"restart_local_par.csv"))
        shutil.copy2(os.path.join(test_d,"pest_local_pdc.obs+noise.csv"),os.path.join(template_d,"restart_local_noise.csv"))
        test_d = os.path.join(model_d, "master_localizer_pdc_test3_{0}".format(out_ext.replace('.','')))
        if os.path.exists(test_d):
            shutil.rmtree(test_d)
        pst.pestpp_options["ies_par_en"] = "restart_local_par.csv"
        pst.pestpp_options["ies_obs_en"] = "restart_local_noise.csv"
        pst.pestpp_options["ies_restart_obs_en"] = "restart_local_obs.csv"
        pst.write(pst_name)
        pyemu.os_utils.start_workers(template_d, exe_path, "pest_local_pdc.pst", num_workers=10,
                                       master_dir=test_d, verbose=True, worker_root=model_d,
                                       port=port)
        phi_df2 = pd.read_csv(os.path.join(test_d,"pest_local_pdc.phi.meas.csv"))
        assert phi_df2.shape[0] == pst.control_data.noptmax+1
        assert phi_df2.loc[phi_df1.index[-1],"mean"] < phi_df2.loc[phi_df1.index[0],"mean"]
        final_phi2 = phi_df2.loc[phi_df2.index[-1],"mean"]
        assert np.abs(final_phi1 - final_phi2) < 0.001,"{0},{1}".format(final_phi1,final_phi2)
        final_phis.append(phi_df2.loc[phi_df2.index[-1],"mean"])
        mat = pyemu.Matrix.from_ascii(os.path.join(test_d, "initialized_localizer.mat"))
        assert mat.shape == (pst.nnz_obs - 1, 2), str(mat.shape)
    final_phis = np.array(final_phis)
    assert np.abs(final_phis.min() - final_phis.max()) < 0.001, "{0},{1}".format(final_phis.min(),final_phis.max())     


def tenpar_localizer_pdc_forgive_test():
    """tenpar local test with  pdc causing some pars to be completely localized out - with extra rows/cols
    """
    model_d = "ies_10par_xsec"
    template_d = os.path.join(model_d, "test_template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)
    cov = pyemu.Cov.from_parameter_data(pst)
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst, cov=cov, num_reals=10)
    pe.enforce()
    oe = pyemu.ObservationEnsemble.from_gaussian_draw(pst,num_reals=10)
    pe.to_csv(os.path.join(template_d, "par_local.csv"))
    oe.to_csv(os.path.join(template_d,"obs_local.csv"))
    final_phis = []
    for out_ext in [".mat",".csv",".jcb"]:
        test_d = os.path.join(model_d, "master_localizer_pdc_forgive_test_{0}".format(out_ext.replace('.','')))
        if os.path.exists(test_d):
            shutil.rmtree(test_d)
        #shutil.copytree(template_d, test_d)
        pst_name = os.path.join(template_d, "pest.pst")
        pst = pyemu.Pst(pst_name)

    
        #mat = pyemu.Matrix.from_names(pst.nnz_obs_names,pst.adj_par_names).to_dataframe()
        #mat = pyemu.Matrix.from_names(pst.nnz_obs_names,pst.adj_par_names).to_dataframe()
        mat = pyemu.Matrix.from_names(pst.obs_names,pst.par_names).to_dataframe()
        mat.loc[:,:] = 0.0
        mat.loc[pst.nnz_obs_names[0],pst.adj_par_names[:2]] = 1.0
        mat.loc[pst.nnz_obs_names[1],pst.adj_par_names[2:]] = 1.0
        
        #mat.iloc[0,:] = 1
        loc_name = "localizer" + out_ext
        
        if out_ext == ".mat":
            mat = pyemu.Matrix.from_dataframe(mat)
            mat.to_ascii(os.path.join(template_d,loc_name))
        elif out_ext == ".csv":
            mat.to_csv(os.path.join(template_d,loc_name))
        elif out_ext == ".jcb":
            mat = pyemu.Matrix.from_dataframe(mat)
            mat.to_binary(os.path.join(template_d,loc_name))

        

        pst.pestpp_options = {}
        pst.pestpp_options["ies_num_reals"] = 10
        pst.pestpp_options["ies_localizer"] = loc_name
        pst.pestpp_options["ies_lambda_mults"] = 1.0
        pst.pestpp_options["lambda_scale_fac"] = 1.0
        pst.pestpp_options["ies_subset_size"] = 11
        pst.pestpp_options["ies_par_en"] = "par_local.csv"
        pst.pestpp_options["ies_obs_en"] = "obs_local.csv"
        pst.pestpp_options["ies_drop_conflicts"] = True
        pst.pestpp_options["ies_verbose_level"] = 4
        pst.pestpp_options["ies_localizer_forgive_missing"] = True
        pst.control_data.noptmax = 2

        #pst.pestpp_options["ies_verbose_level"] = 3
        pst_name = os.path.join(template_d,"pest_local_pdc.pst")
        pst.write(pst_name)
        pyemu.os_utils.start_workers(template_d, exe_path, "pest_local_pdc.pst", num_workers=10,
                                       master_dir=test_d, verbose=True, worker_root=model_d,
                                       port=port)
        phi_df1 = pd.read_csv(os.path.join(test_d,"pest_local_pdc.phi.meas.csv"))
        assert phi_df1.shape[0] == pst.control_data.noptmax+1
        assert phi_df1.loc[phi_df1.index[-1],"mean"] < phi_df1.loc[phi_df1.index[0],"mean"]
        mat = pyemu.Matrix.from_ascii(os.path.join(test_d, "initialized_localizer.mat"))
        assert mat.shape == (pst.nnz_obs, pst.npar_adj), str(mat.shape)

        # now with a group-based localizer
        test_d = os.path.join(model_d, "master_localizer_pdc_forgive_test2_{0}".format(out_ext.replace('.','')))
        if os.path.exists(test_d):
            shutil.rmtree(test_d)
        par = pst.parameter_data
        par.loc[pst.adj_par_names[:4],"pargp"] = "group1"
        par.loc[pst.adj_par_names[4:],"pargp"] = "group2"

        mat = pyemu.Matrix.from_names(pst.obs_names,["group1","group2"]).to_dataframe()
        mat.loc[:,:] = 0.0
        mat.loc[pst.nnz_obs_names[0],"group1"] = 1.0
        mat.loc[pst.nnz_obs_names[1],"group2"] = 1.0
        
        if out_ext == ".mat":
            mat = pyemu.Matrix.from_dataframe(mat)
            mat.to_ascii(os.path.join(template_d,loc_name))
        elif out_ext == ".csv":
            mat.to_csv(os.path.join(template_d,loc_name))
        elif out_ext == ".jcb":
            mat = pyemu.Matrix.from_dataframe(mat)
            mat.to_binary(os.path.join(template_d,loc_name))
        pst.write(pst_name)
        pyemu.os_utils.start_workers(template_d, exe_path, "pest_local_pdc.pst", num_workers=10,
                                       master_dir=test_d, verbose=True, worker_root=model_d,
                                       port=port)
        phi_df1 = pd.read_csv(os.path.join(test_d,"pest_local_pdc.phi.meas.csv"))
        assert phi_df1.shape[0] == pst.control_data.noptmax+1
        assert phi_df1.loc[phi_df1.index[-1],"mean"] < phi_df1.loc[phi_df1.index[0],"mean"]
        mat = pyemu.Matrix.from_ascii(os.path.join(test_d, "initialized_localizer.mat"))
        assert mat.shape == (pst.nnz_obs, 2), str(mat.shape)

        final_phi1 = phi_df1.loc[phi_df1.index[-1],"mean"]
        #now restart

        shutil.copy2(os.path.join(test_d,"pest_local_pdc.0.obs.csv"),os.path.join(template_d,"restart_local_obs.csv"))
        shutil.copy2(os.path.join(test_d,"pest_local_pdc.0.par.csv"),os.path.join(template_d,"restart_local_par.csv"))
        shutil.copy2(os.path.join(test_d,"pest_local_pdc.obs+noise.csv"),os.path.join(template_d,"restart_local_noise.csv"))
        test_d = os.path.join(model_d, "master_localizer_pdc_forgive_test3_{0}".format(out_ext.replace('.','')))
        if os.path.exists(test_d):
            shutil.rmtree(test_d)
        pst.pestpp_options["ies_par_en"] = "restart_local_par.csv"
        pst.pestpp_options["ies_obs_en"] = "restart_local_noise.csv"
        pst.pestpp_options["ies_restart_obs_en"] = "restart_local_obs.csv"
        pst.write(pst_name)
        pyemu.os_utils.start_workers(template_d, exe_path, "pest_local_pdc.pst", num_workers=10,
                                       master_dir=test_d, verbose=True, worker_root=model_d,
                                       port=port)
        phi_df2 = pd.read_csv(os.path.join(test_d,"pest_local_pdc.phi.meas.csv"))
        assert phi_df2.shape[0] == pst.control_data.noptmax+1
        assert phi_df2.loc[phi_df1.index[-1],"mean"] < phi_df2.loc[phi_df1.index[0],"mean"]
        final_phi2 = phi_df2.loc[phi_df2.index[-1],"mean"]
        assert np.abs(final_phi1 - final_phi2) < 0.001,"{0},{1}".format(final_phi1,final_phi2)
        final_phis.append(phi_df2.loc[phi_df2.index[-1],"mean"])
        mat = pyemu.Matrix.from_ascii(os.path.join(test_d, "initialized_localizer.mat"))
        assert mat.shape == (pst.nnz_obs, 2), str(mat.shape)
    final_phis = np.array(final_phis)
    assert np.abs(final_phis.min() - final_phis.max()) < 0.001, "{0},{1}".format(final_phis.min(),final_phis.max())   
        

def tenpar_localizer_pdc_obsgroup_test():
    model_d = "ies_10par_xsec"
    template_d = os.path.join(model_d, "test_template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))

    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)
    # spike one of the obs...

    obs = pst.observation_data
    obs.loc[:, "weight"] = 1.0
    first = pst.obs_names[:4]
    second = pst.obs_names[4:]
    obs.loc[first, "obgnme"] = "group1"
    obs.loc[second, "obgnme"] = "group2"

    obs.loc[first, "obsval"] += 100
    cov = pyemu.Cov.from_parameter_data(pst)
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst, cov=cov, num_reals=10)
    pe.enforce()

    oe = pyemu.ObservationEnsemble.from_gaussian_draw(pst, num_reals=10)
    final_phis = []
    for out_ext in [".mat", ".csv", ".jcb"]:
        test_d = os.path.join(model_d, "master_localizer_pdc_obsgroup_test_{0}".format(out_ext.replace('.', '')))
        if os.path.exists(test_d):
            shutil.rmtree(test_d)
        # shutil.copytree(template_d, test_d)

        pe.to_csv(os.path.join(template_d, "par_local.csv"))
        oe.to_csv(os.path.join(template_d, "obs_local.csv"))

        pst_name = os.path.join(template_d, "pest.pst")
        pst = pyemu.Pst(pst_name)
        # spike one of the obs...
        obs = pst.observation_data
        obs.loc[:,"weight"] = 1.0
        first = pst.obs_names[:4]
        second = pst.obs_names[4:]
        obs.loc[first,"obgnme"] = "group1"
        obs.loc[second, "obgnme"] = "group2"

        obs.loc[first, "obsval"] += 100

        par = pst.parameter_data
        par.loc[pst.adj_par_names[:4], "pargp"] = "group1"
        par.loc[pst.adj_par_names[4:], "pargp"] = "group2"

        # mat = pyemu.Matrix.from_names(pst.nnz_obs_names,pst.adj_par_names).to_dataframe()
        mat = pyemu.Matrix.from_names(["group1","group2"], ["group1","group2"]).to_dataframe()
        mat.loc[:, :] = 0.0
        mat.loc["group1","group1"] = 1.0
        mat.loc["group2","group2"] = 1.0

        # mat.iloc[0,:] = 1
        loc_name = "localizer" + out_ext

        if out_ext == ".mat":
            mat = pyemu.Matrix.from_dataframe(mat)
            mat.to_ascii(os.path.join(template_d, loc_name))
        elif out_ext == ".csv":
            mat.to_csv(os.path.join(template_d, loc_name))
        elif out_ext == ".jcb":
            mat = pyemu.Matrix.from_dataframe(mat)
            mat.to_binary(os.path.join(template_d, loc_name))

        pst.pestpp_options = {}
        pst.pestpp_options["ies_num_reals"] = 10
        pst.pestpp_options["ies_localizer"] = loc_name
        pst.pestpp_options["ies_lambda_mults"] = 1.0
        pst.pestpp_options["lambda_scale_fac"] = 1.0
        pst.pestpp_options["ies_subset_size"] = 11
        pst.pestpp_options["ies_par_en"] = "par_local.csv"
        pst.pestpp_options["ies_obs_en"] = "obs_local.csv"
        pst.pestpp_options["ies_drop_conflicts"] = True
        pst.pestpp_options["ies_verbose_level"] = 4
        pst.pestpp_options["ies_localizer_forgive_missing"] = False
        pst.control_data.noptmax = 2

        # pst.pestpp_options["ies_verbose_level"] = 3
        pst_name = os.path.join(template_d, "pest_local_pdc.pst")
        pst.write(pst_name)
        pyemu.os_utils.start_workers(template_d, exe_path, "pest_local_pdc.pst", num_workers=10,
                                     master_dir=test_d, verbose=True, worker_root=model_d,
                                     port=port)
        phi_df1 = pd.read_csv(os.path.join(test_d, "pest_local_pdc.phi.meas.csv"))
        assert phi_df1.shape[0] == pst.control_data.noptmax + 1
        assert phi_df1.loc[phi_df1.index[-1], "mean"] < phi_df1.loc[phi_df1.index[0], "mean"]
        mat = pyemu.Matrix.from_ascii(os.path.join(test_d, "initialized_localizer.mat"))
        assert mat.shape == (1, 2)


def tenpar_localizer_pdc_pargroup_forgive_test():
    model_d = "ies_10par_xsec"
    template_d = os.path.join(model_d, "test_template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))

    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)
    # spike one of the obs...

    obs = pst.observation_data
    obs.loc[:, "weight"] = 1.0
    first = pst.obs_names[:4]
    second = pst.obs_names[4:]
    obs.loc[first, "obgnme"] = "group1"
    obs.loc[second, "obgnme"] = "group2"
    #obs.loc[first,"weight"] = 0.0

    par = pst.parameter_data
    par.loc[:,"partrans"] = "log"
    par.loc[pst.par_names[:4],"pargp"] = "group1"
    par.loc[pst.par_names[4:],"pargp"] = "group2"
    loc_pnames = pst.par_names[:4]
    loc_pnames.append("group2")



    cov = pyemu.Cov.from_parameter_data(pst)
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst, cov=cov, num_reals=10)
    pe.enforce()

    par.loc[pst.par_names[:4],"partrans"] = "fixed"

    oe = pyemu.ObservationEnsemble.from_gaussian_draw(pst, num_reals=10)
    final_phis = []
    for out_ext in [".mat", ".csv", ".jcb"]:
        test_d = os.path.join(model_d, "master_localizer_pdc_pargroup_forgive_test_{0}".format(out_ext.replace('.', '')))
        if os.path.exists(test_d):
            shutil.rmtree(test_d)
        # shutil.copytree(template_d, test_d)

        pe.to_csv(os.path.join(template_d, "par_local.csv"))
        oe.to_csv(os.path.join(template_d, "obs_local.csv"))

        pst_name = os.path.join(template_d, "pest.pst")
        pst = pyemu.Pst(pst_name)
        # spike one of the obs...
        obs = pst.observation_data
        obs.loc[:,"weight"] = 1.0
        first = pst.obs_names[:4]
        second = pst.obs_names[4:]
        obs.loc[first,"obgnme"] = "group1"
        obs.loc[second, "obgnme"] = "group2"
        #obs.loc[first, "weight"] = 0.0

        par = pst.parameter_data
        par.loc[pst.adj_par_names[:4], "pargp"] = "group1"
        par.loc[pst.adj_par_names[4:], "pargp"] = "group2"
        par.loc[pst.adj_par_names[:4], "partrans"] = "fixed"

        # mat = pyemu.Matrix.from_names(pst.nnz_obs_names,pst.adj_par_names).to_dataframe()

        mat = pyemu.Matrix.from_names(["group1","group2"], loc_pnames).to_dataframe()
        mat.loc[:, :] = 0.0
        mat.loc["group1",pst.par_names[:4]] = 1.0
        mat.loc["group2","group2"] = 1.0

        # mat.iloc[0,:] = 1
        loc_name = "localizer" + out_ext

        if out_ext == ".mat":
            mat = pyemu.Matrix.from_dataframe(mat)
            mat.to_ascii(os.path.join(template_d, loc_name))
        elif out_ext == ".csv":
            mat.to_csv(os.path.join(template_d, loc_name))
        elif out_ext == ".jcb":
            mat = pyemu.Matrix.from_dataframe(mat)
            mat.to_binary(os.path.join(template_d, loc_name))

        pst.pestpp_options = {}
        pst.pestpp_options["ies_num_reals"] = 10
        pst.pestpp_options["ies_localizer"] = loc_name
        pst.pestpp_options["ies_lambda_mults"] = 1.0
        pst.pestpp_options["lambda_scale_fac"] = 1.0
        pst.pestpp_options["ies_subset_size"] = 11
        pst.pestpp_options["ies_par_en"] = "par_local.csv"
        pst.pestpp_options["ies_obs_en"] = "obs_local.csv"
        pst.pestpp_options["ies_drop_conflicts"] = True
        pst.pestpp_options["ies_verbose_level"] = 4
        pst.pestpp_options["ies_localizer_forgive_missing"] = True
        pst.control_data.noptmax = 2

        # pst.pestpp_options["ies_verbose_level"] = 3
        pst_name = os.path.join(template_d, "pest_local_pdc.pst")
        pst.write(pst_name)
        pyemu.os_utils.start_workers(template_d, exe_path, "pest_local_pdc.pst", num_workers=10,
                                     master_dir=test_d, verbose=True, worker_root=model_d,
                                     port=port)
        phi_df1 = pd.read_csv(os.path.join(test_d, "pest_local_pdc.phi.meas.csv"))
        assert phi_df1.shape[0] == pst.control_data.noptmax + 1
        assert phi_df1.loc[phi_df1.index[-1], "mean"] < phi_df1.loc[phi_df1.index[0], "mean"]

        mat = pyemu.Matrix.from_ascii(os.path.join(test_d,"initialized_localizer.mat"))
        assert mat.shape == (2,1),mat.shape

        mat = pyemu.Matrix.from_names(["group1","group2"],["group1","group2"] ).to_dataframe()
        mat.loc[:, :] = 0.0
        mat.loc["group1","group1"] = 1.0
        mat.loc["group2","group2"] = 1.0

        # mat.iloc[0,:] = 1
        loc_name = "localizer" + out_ext

        if out_ext == ".mat":
            mat = pyemu.Matrix.from_dataframe(mat)
            mat.to_ascii(os.path.join(template_d, loc_name))
        elif out_ext == ".csv":
            mat.to_csv(os.path.join(template_d, loc_name))
        elif out_ext == ".jcb":
            mat = pyemu.Matrix.from_dataframe(mat)
            mat.to_binary(os.path.join(template_d, loc_name))

        pst_name = os.path.join(template_d, "pest_local_pdc.pst")
        pst.write(pst_name)
        pyemu.os_utils.start_workers(template_d, exe_path, "pest_local_pdc.pst", num_workers=10,
                                     master_dir=test_d, verbose=True, worker_root=model_d,
                                     port=port)
        phi_df1 = pd.read_csv(os.path.join(test_d, "pest_local_pdc.phi.meas.csv"))
        assert phi_df1.shape[0] == pst.control_data.noptmax + 1
        assert phi_df1.loc[phi_df1.index[-1], "mean"] < phi_df1.loc[phi_df1.index[0], "mean"]


def tenpar_localizer_incomplete_test():
    model_d = "ies_10par_xsec"
    template_d = os.path.join(model_d, "test_template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))

    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)
    # spike one of the obs...

    obs = pst.observation_data
    obs.loc[:, "weight"] = 1.0
    first = pst.obs_names[:4]
    second = pst.obs_names[4:]
    obs.loc[first, "obgnme"] = "group1"
    obs.loc[second, "obgnme"] = "group2"

    cov = pyemu.Cov.from_parameter_data(pst)
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst, cov=cov, num_reals=10)
    pe.enforce()

    oe = pyemu.ObservationEnsemble.from_gaussian_draw(pst, num_reals=10)
    final_phis = []
    for out_ext in [".mat", ".csv", ".jcb"]:
        test_d = os.path.join(model_d, "master_localizer_inc_test_{0}".format(out_ext.replace('.', '')))
        if os.path.exists(test_d):
            shutil.rmtree(test_d)
        # shutil.copytree(template_d, test_d)

        pe.to_csv(os.path.join(template_d, "par_local.csv"))
        oe.to_csv(os.path.join(template_d, "obs_local.csv"))

        pst_name = os.path.join(template_d, "pest.pst")
        pst = pyemu.Pst(pst_name)
        # spike one of the obs...
        obs = pst.observation_data
        obs.loc[:,"weight"] = 1.0
        first = pst.obs_names[:4]
        second = pst.obs_names[4:]
        obs.loc[first,"obgnme"] = "group1"
        obs.loc[second, "obgnme"] = "group2"

        obs.loc[first, "obsval"] += 100

        par = pst.parameter_data
        par.loc[pst.adj_par_names[:4], "pargp"] = "group1"
        par.loc[pst.adj_par_names[4:], "pargp"] = "group2"

        # mat = pyemu.Matrix.from_names(pst.nnz_obs_names,pst.adj_par_names).to_dataframe()
        pnames = [pst.adj_par_names[0],pst.adj_par_names[1]]
        mat = pyemu.Matrix.from_names(["group1","group2"], pnames).to_dataframe()
        mat.loc[:, :] = 0.0
        mat.loc["group1",pnames[0]] = 1.0
        mat.loc["group2",pnames[1]] = 1.0

        # mat.iloc[0,:] = 1
        loc_name = "localizer" + out_ext

        if out_ext == ".mat":
            mat = pyemu.Matrix.from_dataframe(mat)
            mat.to_ascii(os.path.join(template_d, loc_name))
        elif out_ext == ".csv":
            mat.to_csv(os.path.join(template_d, loc_name))
        elif out_ext == ".jcb":
            mat = pyemu.Matrix.from_dataframe(mat)
            mat.to_binary(os.path.join(template_d, loc_name))

        pst.pestpp_options = {}
        pst.pestpp_options["ies_num_reals"] = 10
        pst.pestpp_options["ies_localizer"] = loc_name
        pst.pestpp_options["ies_lambda_mults"] = 1.0
        pst.pestpp_options["lambda_scale_fac"] = 1.0
        pst.pestpp_options["ies_subset_size"] = 11
        pst.pestpp_options["ies_par_en"] = "par_local.csv"
        pst.pestpp_options["ies_obs_en"] = "obs_local.csv"
        pst.pestpp_options["ies_drop_conflicts"] = False
        pst.pestpp_options["ies_verbose_level"] = 4
        pst.pestpp_options["ies_localizer_forgive_missing"] = True
        pst.pestpp_options["ies_save_lambda_en"] = True
        pst.pestpp_options["ies_autoadaloc"] = True
        pst.control_data.noptmax = 2

        # pst.pestpp_options["ies_verbose_level"] = 3
        pst_name = os.path.join(template_d, "pest_local_inc.pst")
        pst.write(pst_name)
        pyemu.os_utils.start_workers(template_d, exe_path, "pest_local_inc.pst", num_workers=10,
                                     master_dir=test_d, verbose=True, worker_root=model_d,
                                     port=port)
        phi_df1 = pd.read_csv(os.path.join(test_d, "pest_local_inc.phi.meas.csv"))
        assert phi_df1.shape[0] == pst.control_data.noptmax + 1
        assert phi_df1.loc[phi_df1.index[-1], "mean"] < phi_df1.loc[phi_df1.index[0], "mean"]
        mat = pyemu.Matrix.from_ascii(os.path.join(test_d, "initialized_localizer.mat"))
        #assert mat.shape == (1, 2)
        print(mat)
        pr_pe = pd.read_csv(os.path.join(test_d,"pest_local_inc.0.par.csv"),index_col=0)
        lam_pe_file = [f for f in os.listdir(test_d) if "lambda" in f and "scale" in f and f.endswith(".par.csv") and f.startswith("pest_local_inc.2.")]
        print(lam_pe_file)
        assert len(lam_pe_file) == 1
        lam_pe = pd.read_csv(os.path.join(test_d,lam_pe_file[0]),index_col=0)
        other_pnames = [n for n in pst.par_names if n not in pnames]
        d = pr_pe.loc[:,other_pnames].values - lam_pe.loc[:,other_pnames].values
        print(np.abs(d).max())
        assert np.abs(d).max() < 1.0e-5


def tenpar_localizer_incomplete_group_test():
    model_d = "ies_10par_xsec"
    template_d = os.path.join(model_d, "test_template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))

    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)
    # spike one of the obs...

    obs = pst.observation_data
    obs.loc[:, "weight"] = 1.0
    first = pst.obs_names[:4]
    second = pst.obs_names[4:]
    obs.loc[first, "obgnme"] = "group1"
    obs.loc[second, "obgnme"] = "group2"

    cov = pyemu.Cov.from_parameter_data(pst)
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst, cov=cov, num_reals=10)
    pe.enforce()

    oe = pyemu.ObservationEnsemble.from_gaussian_draw(pst, num_reals=10)
    final_phis = []
    for out_ext in [".mat", ".csv", ".jcb"]:
        test_d = os.path.join(model_d, "master_localizer_inc_group_test_{0}".format(out_ext.replace('.', '')))
        if os.path.exists(test_d):
            shutil.rmtree(test_d)
        # shutil.copytree(template_d, test_d)

        pe.to_csv(os.path.join(template_d, "par_local.csv"))
        oe.to_csv(os.path.join(template_d, "obs_local.csv"))

        pst_name = os.path.join(template_d, "pest.pst")
        pst = pyemu.Pst(pst_name)
        # spike one of the obs...
        obs = pst.observation_data
        obs.loc[:,"weight"] = 1.0
        first = pst.obs_names[:4]
        second = pst.obs_names[4:]
        obs.loc[first,"obgnme"] = "group1"
        obs.loc[second, "obgnme"] = "group2"

        obs.loc[first, "obsval"] += 100

        par = pst.parameter_data
        par.loc[pst.adj_par_names[:2], "pargp"] = "group1"
        par.loc[pst.adj_par_names[4:], "pargp"] = "group2"

        # mat = pyemu.Matrix.from_names(pst.nnz_obs_names,pst.adj_par_names).to_dataframe()
        pnames = [pst.adj_par_names[0],pst.adj_par_names[1]]
        mat = pyemu.Matrix.from_names(["group1","group2"], ["group1"]).to_dataframe()
        mat.loc[:, :] = 0.0
        mat.loc["group1","group1"] = 1.0
        
        # mat.iloc[0,:] = 1
        loc_name = "localizer" + out_ext

        if out_ext == ".mat":
            mat = pyemu.Matrix.from_dataframe(mat)
            mat.to_ascii(os.path.join(template_d, loc_name))
        elif out_ext == ".csv":
            mat.to_csv(os.path.join(template_d, loc_name))
        elif out_ext == ".jcb":
            mat = pyemu.Matrix.from_dataframe(mat)
            mat.to_binary(os.path.join(template_d, loc_name))

        pst.pestpp_options = {}
        pst.pestpp_options["ies_num_reals"] = 10
        pst.pestpp_options["ies_localizer"] = loc_name
        pst.pestpp_options["ies_lambda_mults"] = 1.0
        pst.pestpp_options["lambda_scale_fac"] = 1.0
        pst.pestpp_options["ies_subset_size"] = 11
        pst.pestpp_options["ies_par_en"] = "par_local.csv"
        pst.pestpp_options["ies_obs_en"] = "obs_local.csv"
        pst.pestpp_options["ies_drop_conflicts"] = False
        pst.pestpp_options["ies_verbose_level"] = 4
        pst.pestpp_options["ies_localizer_forgive_missing"] = True
        pst.pestpp_options["ies_save_lambda_en"] = True
        pst.pestpp_options["ies_autoadaloc"] = True
        pst.control_data.noptmax = 2

        # pst.pestpp_options["ies_verbose_level"] = 3
        pst_name = os.path.join(template_d, "pest_local_inc.pst")
        pst.write(pst_name)
        pyemu.os_utils.start_workers(template_d, exe_path, "pest_local_inc.pst", num_workers=10,
                                     master_dir=test_d, verbose=True, worker_root=model_d,
                                     port=port)
        phi_df1 = pd.read_csv(os.path.join(test_d, "pest_local_inc.phi.meas.csv"))
        assert phi_df1.shape[0] == pst.control_data.noptmax + 1
        assert phi_df1.loc[phi_df1.index[-1], "mean"] < phi_df1.loc[phi_df1.index[0], "mean"]
        mat = pyemu.Matrix.from_ascii(os.path.join(test_d, "initialized_localizer.mat"))
        #assert mat.shape == (1, 2)
        print(mat)
        pr_pe = pd.read_csv(os.path.join(test_d,"pest_local_inc.0.par.csv"),index_col=0)
        lam_pe_file = [f for f in os.listdir(test_d) if "lambda" in f and "scale" in f and f.endswith(".par.csv") and f.startswith("pest_local_inc.2.")]
        print(lam_pe_file)
        assert len(lam_pe_file) == 1
        lam_pe = pd.read_csv(os.path.join(test_d,lam_pe_file[0]),index_col=0)
        other_pnames = [n for n in pst.par_names if n not in pnames]
        d = pr_pe.loc[:,other_pnames].values - lam_pe.loc[:,other_pnames].values
        print(np.abs(d).max())
        assert np.abs(d).max() < 1.0e-5



def tenpar_ineq_test():
    """tenpar ineq test"""
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_ineq_test")
    template_d = os.path.join(model_d, "test_template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)

    pst = pyemu.Pst(os.path.join(template_d, "pest.pst"))
    obs = pst.observation_data
    nznames = pst.nnz_obs_names
    obs.loc[nznames,"obgnme"] = "greater_than"
    

    pst.control_data.noptmax = 2

    num_reals = 10
    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = num_reals
    #pst.pestpp_options["ies_no_noise"] = True
    
    pst.write(os.path.join(test_d,"pest.pst"),version=2)
 
    pyemu.helpers.run(exe_path + " pest.pst", cwd=test_d)
    phidf1 = pd.read_csv(os.path.join(test_d,"pest.phi.meas.csv"))



    test_d = os.path.join(model_d, "master_ineq_test_gtbnd")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    obs.loc[nznames,"obgnme"] = "obsgrp"
    obs.loc[nznames,"greater_than"] = obs.loc[nznames,"obsval"]
    pst.write(os.path.join(test_d,"pest.pst"),version=2)
    pyemu.helpers.run(exe_path + " pest.pst", cwd=test_d)
    phidf2 = pd.read_csv(os.path.join(test_d,"pest.phi.meas.csv"))

    diff = phidf1 - phidf2
    print(diff)
    print(np.abs(diff.values).max())
    assert np.abs(diff.values).max() < 1e-7

    obs.loc[nznames,"obgnme"] = "less_than"
    obs.loc[nznames,"greater_than"] = np.nan
    test_d = os.path.join(model_d, "master_ineq_test")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    #obs.loc[nznames,"obgnme"] = "obsgrp"
    pst.write(os.path.join(test_d,"pest.pst"),version=2)
    pyemu.helpers.run(exe_path + " pest.pst", cwd=test_d)
    phidf1 = pd.read_csv(os.path.join(test_d,"pest.phi.meas.csv"))

    test_d = os.path.join(model_d, "master_ineq_test_ltbnd")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    obs.loc[nznames,"obgnme"] = "obsgrp"
    obs.loc[nznames,"less_than"] = obs.loc[nznames,"obsval"]
    pst.write(os.path.join(test_d,"pest.pst"),version=2)
    pyemu.helpers.run(exe_path + " pest.pst", cwd=test_d)
    phidf2 = pd.read_csv(os.path.join(test_d,"pest.phi.meas.csv"))
    diff = phidf1 - phidf2
    print(diff)
    print(np.abs(diff.values).max())
    assert np.abs(diff.values).max() < 1e-7

    obs.loc[:,"greater_than"] = np.nan
    obs.loc[:,"less_than"] = np.nan
    obs.loc[:,"weight"] = 1.0
    pst.pestpp_options["ies_no_noise"] = True
    test_d = os.path.join(model_d, "master_ineq_test")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    #obs.loc[nznames,"obgnme"] = "obsgrp"
    pst.write(os.path.join(test_d,"pest.pst"),version=2)
    pyemu.helpers.run(exe_path + " pest.pst", cwd=test_d)
    phidf1 = pd.read_csv(os.path.join(test_d,"pest.phi.meas.csv"))

    test_d = os.path.join(model_d, "master_ineq_test_doublebnd")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    obs.loc[:,"greater_than"] = obs.loc[:,"obsval"]
    obs.loc[:,"less_than"] = obs.loc[:,"obsval"]
    pst.write(os.path.join(test_d,"pest.pst"),version=2)
    pyemu.helpers.run(exe_path + " pest.pst", cwd=test_d)
    phidf2 = pd.read_csv(os.path.join(test_d,"pest.phi.meas.csv"))
    diff = phidf1 - phidf2
    print(diff)
    print(np.abs(diff.values).max())
    assert np.abs(diff.values).max() < 1e-7

    m_d = os.path.join(model_d,"master_ineq_test")
    pyemu.os_utils.start_workers(test_d,exe_path,"pest.pst",master_dir=m_d,num_workers=5,worker_root=model_d)




if __name__ == "__main__":
    
    #tenpar_restart_similar_test2()
    #tenpar_localizer_pdc_test()
    #tenpar_localizer_pdc_obsgroup_test()
    #tenpar_localizer_pdc_pargroup_forgive_test()
    #tenpar_localizer_incomplete_test()
    #tenpar_localizer_incomplete_group_test()
    #tenpar_fixed_test3()
    
    # full list of tests
    #tenpar_subset_test()
    #shutil.copy2(os.path.join("..", "exe", "windows", "x64", "Debug", "pestpp-ies.exe"),
    #             os.path.join("..", "bin", "win","pestpp-ies.exe"))
    #invest()
    tenpar_ineq_test()
    #tenpar_restart_similar_test()
    #tenpar_fixed_test()
    # tenpar_full_cov_test()
    # eval_freyberg_full_cov_reorder()
    #test_freyberg_full_cov_reorder()
    # eval_freyberg_full_cov()
    # tenpar_tight_tol_test()
    #test_chenoliver()
    # tenpar_narrow_range_test()
    #test_freyberg_ineq()
    #tenpar_fixed_test()
    #tenpar_full_cov_test()
    #tenpar_fixed_test2()
    # tenpar_subset_how_test()
    # tenpar_localizer_test1()
    # tenpar_localizer_test2()
    # tenpar_localizer_test3()
    # freyberg_localizer_eval1()
    # freyberg_localizer_eval2()
    # freyberg_localizer_test3()
    # freyberg_dist_local_test()
    # freyberg_local_threads_test()
    #tenpar_restart_binary_test()
    #tenpar_weights_test()
    #tenpar_restart_test()
    #csv_tests()
    # tenpar_rns_test()
    # clues_longnames_test()
    # tenpar_localize_how_test()
    #shutil.copy2(os.path.join("..", "exe", "windows", "x64", "Debug", "pestpp-ies.exe"),
    #             os.path.join("..", "bin", "win", "pestpp-ies.exe"))
    #tenpar_localizer_test3()
    #import pyemu
    #m = pyemu.Matrix.from_binary(os.path.join("ies_10par_xsec","master_fixed","pest_fixed.0.par.jcb"))
    #tenpar_fixed_test()
    #tenpar_incr_num_reals_test()
    #test_freyberg_ineq()
    #freyberg_dist_local_invest()
    #test_freyberg_full_cov_reorder_run()
    #test_freyberg_full_cov()