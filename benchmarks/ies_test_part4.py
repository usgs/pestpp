# TODO: test variance and mean of draws, add chenoliver and test approx and full solution
import os
import shutil
import platform
import numpy as np
import scipy.linalg
import pandas as pd
import platform
import matplotlib.pyplot as plt
import pyemu

tests = """0) 10par_xsec "standard user mode" - draw reals from par-bounds prior and obs noise from weights
0a) 10par_xsec same as 0) but with multple lambda 
1) 10par_xsec start with existing par csv and obs csv - using empirical parcov and obscov
1a) 10par_xsec start with existing par csv and obs csv - using parcov file
2) 10par_xsec start with existing par csv and drawing obs en from weights 
3) 10par_xsec restart with full simulated obs en
3a) 10par_xsec restart with failed runs in simulated obs en
3b) 10par_xsec restart with failed runs and bad phi runs in simulated obs en with multiple lam
4) 10par_xsec reg_factor = 0.5 test
5)  10par_xsec full solution test with standard draw mode
5a) 10par_xsec full solution test with empirical parcov
6) freyberg "standard user mode" - draw reals from par-bounds prior and obs noise from weights
6a) freyberg same as 0) but with multple lambda 
7) freyberg draw par en from full parcov supplied in file
8) freyberg full solution with empirical parcov - supplied par csv, obs csv and restart csv with fails, bad phi,MAP solution, prior scaling, lam mults 
9) synth restart and upgrade 1.1M par problem"""

ies_vars = ["ies_par_en", "ies_obs_en", "ies_restart_obs_en",
            "ies_bad_phi", "parcov_filename", "ies_num_reals",
            "ies_use_approx", "ies_use_prior_scaling", "ies_reg_factor",
            "ies_lambda_mults", "ies_initial_lambda", "ies_include_base", "ies_subset_size"]

# the old path system before moving to separate benchmarks repo
# intel = False
# if "windows" in platform.platform().lower():
#     if intel:
#         exe_path = os.path.join("..", "..", "..", "bin", "iwin", "ipestpp-ies.exe")
#     else:
#         exe_path = os.path.join("..", "..", "..", "bin", "win", "pestpp-ies.exe")
# elif "darwin" in platform.platform().lower():
#     exe_path = os.path.join("..", "..", "..", "bin", "mac", "pestpp-ies")
# else:
#     exe_path = os.path.join("..", "..", "..", "bin", "linux", "pestpp-ies")

bin_path = os.path.join("test_bin")
if "linux" in platform.platform().lower():
    bin_path = os.path.join(bin_path, "linux")
elif "darwin" in platform.platform().lower() or "macos" in platform.platform().lower():
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
    exe_path = os.path.join(bin_path, "win", "pestpp-ies.exe")
elif "darwin" in platform.platform().lower() or "macos" in platform.platform().lower():
    exe_path = os.path.join(bin_path, "mac", "pestpp-ies")
else:
    exe_path = os.path.join(bin_path, "linux", "pestpp-ies")

noptmax = 3

compare_files = ["pest.phi.actual.csv", "pest.phi.meas.csv", "pest.phi.regul.csv",
                 "pest.{0}.par.csv".format(noptmax), "pest.{0}.obs.csv".format(noptmax),
                 "pest.{0}.par.csv".format(0), "pest.obs+noise.csv"]
diff_tol = 1.0e-6
port = 4016
num_reals = 10


def tenpar_xsec_aal_sigma_dist_test():
    """testing what happens with really large sigma dist for aal"""

    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_aal_sigma_dist_test")
    template_d = os.path.join(model_d, "test_template")

    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)

    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 50
    pst.pestpp_options["ies_autoadaloc_sigma_dist"] = 2.0
    pst.pestpp_options["ies_autoadaloc"] = True
    pst.pestpp_options["ies_verbose_level"] = 4
    pst.control_data.noptmax = 4
    pst.write(os.path.join(template_d, "pest_aal_sigma_dist.pst"))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d,test_d)


    # pyemu.os_utils.start_workers(template_d, exe_path, "pest_aal_sigma_dist.pst", num_workers=10,
    #                              master_dir=test_d, verbose=True, worker_root=model_d,
    #                             port=port)
    pyemu.os_utils.run("{0} pest_aal_sigma_dist.pst".format(exe_path),cwd=test_d)


    df = pd.read_csv(os.path.join(test_d, "pest_aal_sigma_dist.1.autoadaloc.csv"))
    df.loc[:, "parnme"] = df.parnme.str.lower()
    df.loc[:, "obsnme"] = df.obsnme.str.lower()
    print(df.iloc[0, :])

    test_d = test_d+"_dummy"
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d,test_d)
    with open(os.path.join(test_d,"dummy.dat.tpl"),'w') as f:
        f.write("ptf ~\n")
        f.write(" ~ dummy1  ~\n ~ dummy2   ~\n~   dummy3    ~\n ~    dummy4   ~\n")
    df = pst.add_parameters(os.path.join(test_d,"dummy.dat.tpl"),os.path.join(test_d,"dummy.dat"),pst_path=".")
    pst.parameter_data.loc[df.parnme,"pargp"] = "dummy"
    pst.pestpp_options["ies_autoadaloc_indicator_pars"] = df.parnme.tolist()
    pst.write(os.path.join(test_d, "pest_aal_sigma_dist_indicator.pst"))
    pyemu.os_utils.run("{0} pest_aal_sigma_dist_indicator.pst".format(exe_path),cwd=test_d)

    phi = pd.read_csv(os.path.join(test_d,"pest_aal_sigma_dist_indicator.phi.actual.csv"))
    pcs = pd.read_csv(os.path.join(test_d,"pest_aal_sigma_dist_indicator.{0}.pcs.csv".format(phi.iteration.max())),index_col=0)
    print(pcs.loc["dummy","mean_change"])
    assert np.abs(pcs.loc["dummy","mean_change"]) < 1e-7
    print(pcs.loc["dummy","std_change"])
    assert np.abs(pcs.loc["dummy","std_change"]) < 1e-7
    


    # fig,axes = plt.subplots(pst.npar_adj,pst.nnz_obs,figsize=(6.5,11))

    # for i,pname in enumerate(pst.adj_par_names):
    #     for j,oname in enumerate(pst.nnz_obs_names):
    #         ddf = df.loc[df.apply(lambda x: x.parnme==pname and x.obsnme==oname,axis=1),:].iloc[0,:]
    #         ax = axes[i,j]
    #         print(ddf.iloc[6:])
    #         ddf.iloc[6:].apply(np.float).hist(ax=ax,bins=10,facecolor='b',alpha=0.5)
    #         ylim = ax.get_ylim()
    #         ylim = [0,ylim[1]*1.5]
    #         ax.plot([ddf.loc["correlation_coeff"],ddf.loc["correlation_coeff"]],ylim,"r",label="estimate")
    #         mn,std = ddf.loc["background_mean"], ddf.loc["background_stdev"]
    #         ax.plot([mn,mn],ylim,"b-", label="bg mean")
    #         for m,c in zip([1,2,3],["--","-.",":"]):
    #             ax.plot([mn+(m*std),mn+(m*std)],ylim,color="b",ls=c,label="{0} bg std".format(m))
    #             ax.plot([mn - (m * std), mn - (m * std)], ylim, color="b",ls=c)
    #         ax.grid()
    #         ax.set_ylim(ylim)
    #         ax.set_xlim(-1,1)
    #         ax.set_yticks([])
    #         kept = bool(ddf.loc["kept"])
    #         ax.set_title("{0} - {1}, kept: {2}".format(pname,oname,kept),loc="left")
    # ax.legend()
    # plt.tight_layout()
    # plt.savefig("aal_10par_2sigma.pdf")
    # plt.show()


def tenpar_xsec_aal_invest():
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_aal_test")
    template_d = os.path.join(model_d, "test_template")

    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)

    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 100
    pst.pestpp_options["ies_lambda_mults"] = 0.0000001
    pst.pestpp_options["lambda_scale_fac"] = 0.00001
    # pst.pestpp_options["ies_autoadaloc"] = True
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst, cov=pyemu.Cov.from_parameter_data(pst), num_reals=100)
    pe.loc[:, pst.adj_par_names[:2]] = pst.parameter_data.loc[pst.adj_par_names[0], "parlbnd"]

    pe.to_csv(os.path.join(template_d, "bound_par.csv"))
    pst.pestpp_options["ies_par_en"] = "bound_par.csv"
    pst.pestpp_options["ies_enforce_bounds"] = False
    pst.pestpp_options["ies_verbose_level"] = 3
    pst.control_data.noptmax = 1
    pst.write(os.path.join(template_d, "pest_aal.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_aal.pst", num_workers=30,
                                 master_dir=test_d, verbose=True, worker_root=model_d,
                                 port=port)

    # pst.pestpp_options = {}
    # pst.control_data.noptmax = -1
    # pst.write(os.path.join(test_d, "pest_aal_jco.pst"))
    # pyemu.os_utils.run("{0} pest_aal_jco.pst".format(exe_path.replace("-ies","-glm")),cwd=test_d)

    # df = pd.read_csv(os.path.join(test_d,"pest_aal."))


def tenpar_base_run_test():
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_base_test")
    template_d = os.path.join(model_d, "test_template")

    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)
    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options["ies_include_base"] = True
    pst.control_data.noptmax = -1

    pst.write(os.path.join(template_d, "pest_base.pst"))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pyemu.os_utils.run("{0} pest_base.pst".format(exe_path), cwd=test_d)

    pst.control_data.noptmax = 0
    pst.write(os.path.join(template_d, "pest_base.pst"))
    pyemu.os_utils.run("{0} pest_base.pst".format(exe_path.replace("-ies", "-glm")), cwd=test_d)

    oe = pd.read_csv(os.path.join(test_d, "pest_base.0.obs.csv"), index_col=0)
    pst = pyemu.Pst(os.path.join(test_d, "pest_base.pst"))
    print(oe.loc["base", :])
    print(pst.res.modelled)
    d = oe.loc["base", :] - pst.res.modelled
    assert d.sum() == 0.0, d

    pst.control_data.noptmax = 0
    pst.observation_data.loc[:, "weight"] = 0.0
    pst.write(os.path.join(test_d, "pest_base.pst"))
    pyemu.os_utils.run("{0} pest_base.pst".format(exe_path), cwd=test_d)


def tenpar_base_par_file_test():
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_parfile1")
    template_d = os.path.join(model_d, "test_template")

    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)
    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options["ies_lambda_mults"] = 1.0
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["ies_include_base"] = True
    pst.control_data.noptmax = 2

    pst.write(os.path.join(template_d, "pest_base.pst"))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pyemu.os_utils.run("{0} pest_base.pst".format(exe_path), cwd=test_d)
    assert os.path.exists(os.path.join(test_d, "pest_base.1.base.par"))
    assert os.path.exists(os.path.join(test_d, "pest_base.2.base.par"))

    pe = pd.read_csv(os.path.join(test_d, "pest_base.1.par.csv"), index_col=0)
    pvals = pyemu.pst_utils.read_parfile(os.path.join(test_d, "pest_base.1.base.par"))
    d = pe.loc["base", pvals.index.values].values - pvals.parval1.values
    print(d)
    assert np.abs(d).max() < 1.0e-5
    pe = pd.read_csv(os.path.join(test_d, "pest_base.2.par.csv"), index_col=0)
    pvals = pyemu.pst_utils.read_parfile(os.path.join(test_d, "pest_base.2.base.par"))
    d = pe.loc["base", pvals.index.values].values - pvals.parval1.values
    print(d)
    assert np.abs(d).max() < 1.0e-5


def tenpar_xsec_combined_autoadaloc_test():
    """testing combined matrix + autoadaloc"""
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_comb_aal_test1")
    template_d = os.path.join(model_d, "test_template")

    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)

    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 50
    pst.control_data.noptmax = -1
    pst.write(os.path.join(template_d, "pest_aal.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_aal.pst", num_workers=10,
                                 master_dir=test_d, verbose=True, worker_root=model_d,
                                 port=port)

    mat = pyemu.Matrix.from_names(pst.nnz_obs_names, pst.adj_par_names).to_dataframe()
    mat.loc[:, :] = 1
    mat.loc[:, pst.adj_par_names[::2]] = 0
    pyemu.Matrix.from_dataframe(mat).to_ascii(os.path.join(template_d, "loc.mat"))

    pst.pestpp_options["ies_localizer"] = "loc.mat"
    pst.pestpp_options["ies_par_en"] = "pest_aal.0.par.csv"
    pst.pestpp_options["ies_obs_en"] = "pest_aal.obs+noise.csv"
    pst.pestpp_options["ies_restart_obs_en"] = "pest_aal.0.obs.csv"
    pst.pestpp_options["ies_autoadaloc"] = True
    pst.pestpp_options["ies_verbose_level"] = 3
    pst.control_data.noptmax = 1

    pe = pyemu.ParameterEnsemble.from_dataframe(df=pd.read_csv(os.path.join(test_d, "pest_aal.0.par.csv"), index_col=0),
                                                pst=pst)

    oe = pyemu.ObservationEnsemble.from_dataframe(
        df=pd.read_csv(os.path.join(test_d, "pest_aal.0.obs.csv"), index_col=0), pst=pst)

    for f in ["pest_aal.0.par.csv", "pest_aal.obs+noise.csv", "pest_aal.0.obs.csv"]:
        shutil.copy2(os.path.join(test_d, f), os.path.join(template_d, f))

    pst.write(os.path.join(template_d, "pest_aal_restart.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_aal_restart.pst", num_workers=10,
                                 master_dir=test_d, verbose=True, worker_root=model_d,
                                 port=port)
    df = pyemu.Matrix.from_ascii(os.path.join(test_d, "pest_aal_restart.1.autoadaloc.tCC.mat")).to_dataframe()
    print(df.loc[:, pst.adj_par_names[::2]].sum())
    pe2 = pd.read_csv(os.path.join(test_d, "pest_aal_restart.0.par.csv"))
    diff = pe - pe2
    print(diff.loc[:, pst.adj_par_names[::2]].sum())
    assert df.loc[:, pst.adj_par_names[::2]].sum().sum() == 0.0
    assert diff.loc[:, pst.adj_par_names[::2]].sum().sum() == 0.0


def freyberg_aal_test():
    import flopy
    model_d = "ies_freyberg"
    test_d = os.path.join(model_d, "master_aal_test")
    template_d = os.path.join(model_d, "template")
    m = flopy.modflow.Modflow.load("freyberg.nam", model_ws=template_d, load_only=[], check=False)
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    # print("loading pst")
    pst = pyemu.Pst(os.path.join(template_d, "pest.pst"))

    par = pst.parameter_data

    par = pst.parameter_data
    par.loc[:, "partrans"] = "fixed"
    par.loc[par.pargp == "hk", "partrans"] = "log"

    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options["ies_subset_size"] = 10
    pst.pestpp_options["ies_num_threads"] = 20
    pst.pestpp_options["ies_lambda_mults"] = [1.0]
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    # pst.pestpp_options["ies_include_base"] = False
    # pst.pestpp_options["ies_par_en"] = "par_local.csv"
    pst.pestpp_options["ies_use_approx"] = False
    pst.pestpp_options["ies_use_prior_scaling"] = True
    # pst.pestpp_options["ies_localizer"] = "localizer.mat"
    pst.pestpp_options["ies_localize_how"] = "par"
    pst.pestpp_options["ies_verbose_level"] = 2
    pst.pestpp_options["ies_save_lambda_en"] = True
    pst.pestpp_options["ies_subset_how"] = "random"
    pst.pestpp_options["ies_accept_phi_fac"] = 1000.0
    pst.pestpp_options["overdue_giveup_fac"] = 10.0
    pst.pestpp_options["ies_autoadaloc"] = True
    pst.pestpp_options["ies_autoadaloc_sigma_dist"] = 1.0
    pst.control_data.noptmax = 1
    pst.write(os.path.join(template_d, "pest_aal.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_aal.pst", num_workers=30, master_dir=test_d,
                                 worker_root=model_d, port=port)


def freyberg_combined_aal_test():
    import flopy
    model_d = "ies_freyberg"
    test_d = os.path.join(model_d, "master_combined_aal_test")
    template_d = os.path.join(model_d, "template")
    m = flopy.modflow.Modflow.load("freyberg.nam", model_ws=template_d, load_only=[], check=False)
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    # print("loading pst")
    pst = pyemu.Pst(os.path.join(template_d, "pest.pst"))

    m = flopy.modflow.Modflow.load("freyberg.nam", model_ws=template_d, load_only=[], check=False)
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    # print("loading pst")

    par = pst.parameter_data
    par.loc[:, "partrans"] = "fixed"
    par.loc[par.pargp == "hk", "partrans"] = "log"

    par_adj = par.loc[pst.adj_par_names, :].copy()
    par_adj.loc[:, "i"] = par_adj.parnme.apply(lambda x: int(x.split('_')[1][1:]))
    par_adj.loc[:, "j"] = par_adj.parnme.apply(lambda x: int(x.split('_')[2][1:]))
    par_adj.loc[:, "x"] = par_adj.apply(lambda x: m.modelgrid.xcellcenters[x.i, x.j], axis=1)
    par_adj.loc[:, "y"] = par_adj.apply(lambda x: m.modelgrid.ycellcenters[x.i, x.j], axis=1)

    pst.observation_data.loc["flx_river_l_19700102", "weight"] = 0.0
    obs_nz = pst.observation_data.loc[pst.nnz_obs_names, :].copy()
    obs_nz.loc[:, "i"] = obs_nz.obsnme.apply(lambda x: int(x[6:8]))
    obs_nz.loc[:, "j"] = obs_nz.obsnme.apply(lambda x: int(x[9:11]))
    obs_nz.loc[:, 'x'] = obs_nz.apply(lambda x: m.modelgrid.xcellcenters[x.i, x.j], axis=1)
    obs_nz.loc[:, 'y'] = obs_nz.apply(lambda x: m.modelgrid.ycellcenters[x.i, x.j], axis=1)

    dfs = []
    v = pyemu.geostats.ExpVario(contribution=1.0, a=1000)
    for name in pst.nnz_obs_names:
        x, y = obs_nz.loc[name, ['x', 'y']].values
        # print(name,x,y)
        p = par_adj.copy()
        # p.loc[:,"dist"] = p.apply(lambda xx: np.sqrt((xx.x - x)**2 + (xx.y - y)**2),axis=1)
        # print(p.dist.max(),p.dist.min())
        cc = v.covariance_points(x, y, p.x, p.y)
        # print(cc.min(),cc.max())
        dfs.append(cc)

    df = pd.concat(dfs, axis=1)
    df.columns = pst.nnz_obs_names

    mat = pyemu.Matrix.from_dataframe(df.T)
    tol = 0.35

    mat.x[mat.x < tol] = 0.0
    mat.to_ascii(os.path.join(template_d, "localizer.mat"))
    df_tol = mat.to_dataframe()
    par_sum = df_tol.sum(axis=0)

    zero_cond_pars = list(par_sum.loc[par_sum == 0.0].index)
    print(zero_cond_pars)

    par = pst.parameter_data

    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options["ies_subset_size"] = 10
    pst.pestpp_options["ies_num_threads"] = 20
    pst.pestpp_options["ies_lambda_mults"] = [1.0]
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    # pst.pestpp_options["ies_include_base"] = False
    # pst.pestpp_options["ies_par_en"] = "par_local.csv"
    pst.pestpp_options["ies_use_approx"] = False
    pst.pestpp_options["ies_use_prior_scaling"] = True
    pst.pestpp_options["ies_localizer"] = "localizer.mat"
    pst.pestpp_options["ies_localize_how"] = "par"
    pst.pestpp_options["ies_verbose_level"] = 2
    pst.pestpp_options["ies_save_lambda_en"] = True
    pst.pestpp_options["ies_subset_how"] = "random"
    pst.pestpp_options["ies_accept_phi_fac"] = 1000.0
    pst.pestpp_options["overdue_giveup_fac"] = 10.0
    pst.pestpp_options["ies_autoadaloc"] = True
    pst.pestpp_options["ies_autoadaloc_sigma_dist"] = 1.0
    pst.control_data.noptmax = 1
    pst.write(os.path.join(template_d, "pest_aal.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_aal.pst", num_workers=30, master_dir=test_d,
                                 worker_root=model_d, port=port)

    pr = pd.read_csv(os.path.join(test_d, "pest_aal.0.par.csv")).loc[:, zero_cond_pars]
    pt = pd.read_csv(os.path.join(test_d, "pest_aal.{0}.par.csv".format(pst.control_data.noptmax))).loc[:,
         zero_cond_pars]

    diff = pr - pt
    print(diff.apply(lambda x: np.abs(x)).max().max())
    assert diff.apply(lambda x: np.abs(x)).max().max() == 0.0


def freyberg_aal_invest():
    import flopy
    model_d = "ies_freyberg"
    test_d = os.path.join(model_d, "master_aal_glm_jco")
    template_d = os.path.join(model_d, "template")
    jco_file = os.path.join(test_d, "pest_aal_jco.jcb")
    if not os.path.exists(jco_file):
        pst = pyemu.Pst(os.path.join(template_d, "pest.pst"))
        pst.control_data.noptmax = -1
        pst.write(os.path.join(template_d, "pest_aal_jco.pst"))
        pyemu.os_utils.start_workers(template_d, exe_path.replace("-ies", "-glm"), "pest_aal_jco.pst", 30,
                                     worker_root=model_d, master_dir=test_d, port=port, verbose=True)
    jco = pyemu.Jco.from_binary(jco_file).to_dataframe()

    test_d = os.path.join(model_d, "master_combined_aal_test")
    m = flopy.modflow.Modflow.load("freyberg.nam", model_ws=template_d, load_only=[], check=False)
    tcc = pyemu.Matrix.from_ascii(os.path.join(test_d, "pest_aal.1.autoadaloc.tCC.mat")).to_dataframe()
    pnames = pd.DataFrame({"name": tcc.columns.values})
    pnames.loc[:, "i"] = pnames.name.apply(lambda x: int(x.split('_')[1][1:]))
    pnames.loc[:, "j"] = pnames.name.apply(lambda x: int(x.split('_')[2][1:]))
    pdict = {n: (i, j) for n, i, j in zip(pnames.name, pnames.i, pnames.j)}

    test_d = os.path.join(model_d, "master_aal_test")
    tcc2 = pyemu.Matrix.from_ascii(os.path.join(test_d, "pest_aal.1.autoadaloc.tCC.mat")).to_dataframe()
    pnames2 = pd.DataFrame({"name": tcc2.columns.values})
    pnames2.loc[:, "i"] = pnames2.name.apply(lambda x: int(x.split('_')[1][1:]))
    pnames2.loc[:, "j"] = pnames2.name.apply(lambda x: int(x.split('_')[2][1:]))
    pdict2 = {n: (i, j) for n, i, j in zip(pnames2.name, pnames2.i, pnames2.j)}

    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(os.path.join(test_d, "compare_sens_cc.pdf")) as pdf:
        for obs in tcc.index:
            i, j = None, None
            if obs.startswith("c00"):
                # continue
                i = int(obs[6:8])
                j = int(obs[9:11])
            tcc_obs = tcc.loc[obs, :].apply(lambda x: np.abs(x))
            jco_obs = jco.loc[obs, :].apply(lambda x: np.abs(x))
            print(tcc_obs)
            arr_cc = np.zeros((m.nrow, m.ncol))
            arr_jco = np.zeros((m.nrow, m.ncol))

            tcc_obs2 = tcc2.loc[obs, :].apply(lambda x: np.abs(x))
            arr_cc2 = np.zeros((m.nrow, m.ncol))

            for n, v, v2 in zip(tcc_obs.index, tcc_obs.values, tcc_obs2.values):
                if not "hk" in n:
                    continue
                arr_cc[pdict[n][0], pdict[n][1]] = v
                arr_jco[pdict[n][0], pdict[n][1]] = jco_obs[n]
                arr_cc2[pdict2[n][0], pdict2[n][1]] = v2
            fig = plt.figure(figsize=(13, 6))

            ax = plt.subplot(131, aspect="equal")
            ax2 = plt.subplot(132, aspect="equal")
            ax3 = plt.subplot(133, aspect="equal")

            arr_cc = arr_cc / arr_cc.max()
            arr_jco = arr_jco / arr_jco.max()
            arr_cc2 = arr_cc2 / arr_cc2.max()
            arr_cc = np.ma.masked_where(arr_cc < 1.0e-6, arr_cc)
            arr_jco = np.ma.masked_where(arr_jco < 1.0e-6, arr_jco)
            c = ax.pcolormesh(m.modelgrid.xcellcenters, m.modelgrid.ycellcenters, arr_cc2, alpha=0.5, vmin=0, vmax=1)
            # plt.colorbar(c,ax=ax)
            c1 = ax2.pcolormesh(m.modelgrid.xcellcenters, m.modelgrid.ycellcenters, arr_cc, alpha=0.5, vmin=0, vmax=1)
            c2 = ax3.pcolormesh(m.modelgrid.xcellcenters, m.modelgrid.ycellcenters, arr_jco, alpha=0.5, vmin=0, vmax=1)

            # plt.colorbar(c2,ax=ax2,fraction=0.046, pad=0.04)
            if i is not None:
                ax.scatter([m.modelgrid.xcellcenters[i, j]], [m.modelgrid.ycellcenters[i, j]], marker='.', s=50)
                ax2.scatter([m.modelgrid.xcellcenters[i, j]], [m.modelgrid.ycellcenters[i, j]], marker='.', s=50)
                ax3.scatter([m.modelgrid.xcellcenters[i, j]], [m.modelgrid.ycellcenters[i, j]], marker='.', s=50)

            ax.set_title("A.) estimated CC".format(obs), fontsize=12, loc="left")
            ax2.set_title("B.) distance loc + estimated CC".format(obs), fontsize=12, loc="left")
            ax3.set_title("C.) normalized JCO row".format(obs), fontsize=12, loc="left")
            for ax in [ax, ax2, ax3]:
                ax.set_xticks([])
                ax.set_yticks([])

            plt.tight_layout()
            pdf.savefig()
            plt.close(fig)


def tenpar_high_phi_test():
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_high_phi_test")
    template_d = os.path.join(model_d, "test_template")

    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)

    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    # shutil.copytree(template_d, test_d)
    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options["ies_lambda_mults"] = 1.0
    pst.pestpp_options["lambda_scale_fac"] = [0.9, 1.0]
    pst.pestpp_options['ies_subset_size'] = 10
    pst.pestpp_options["ies_debug_high_subset_phi"] = True
    pst.pestpp_options["ies_update_by_reals"] = False
    pst.control_data.noptmax = 1
    pst.write(os.path.join(template_d, "pest_high_phi.pst"))
    #pyemu.os_utils.start_workers(template_d, exe_path, "pest_high_phi.pst", num_workers=10,
    #                             master_dir=test_d, verbose=True, worker_root=model_d,
    #                             port=port)
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pyemu.os_utils.run("{0} {1}".format(exe_path,"pest_high_phi.pst"),cwd=test_d)
    phi1 = pd.read_csv(os.path.join(test_d, "pest_high_phi.phi.actual.csv"), index_col=0)
    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options["ies_lambda_mults"] = 1.0
    pst.pestpp_options["lambda_scale_fac"] = [0.9, 1.0]
    pst.pestpp_options['ies_subset_size'] = 10
    pst.pestpp_options["ies_update_by_reals"] = False
    # pst.pestpp_options["ies_debug_high_subset_phi"] =True
    pst.control_data.noptmax = 1
    pst.write(os.path.join(template_d, "pest_high_phi.pst"))
    #pyemu.os_utils.start_workers(template_d, exe_path, "pest_high_phi.pst", num_workers=10,
    #                             master_dir=test_d, verbose=True, worker_root=model_d,
    #                             port=port)
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pyemu.os_utils.run("{0} {1}".format(exe_path,"pest_high_phi.pst"),cwd=test_d)  
    phi2 = pd.read_csv(os.path.join(test_d, "pest_high_phi.phi.actual.csv"), index_col=0)
    diff = phi1 - phi2
    print(diff)
    assert diff.max().max() == 0.0

    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options["ies_lambda_mults"] = 1.0
    pst.pestpp_options["lambda_scale_fac"] = [0.9, 1.0]
    pst.pestpp_options['ies_subset_size'] = 10
    pst.pestpp_options["ies_debug_high_upgrade_phi"] = True
    pst.pestpp_options["ies_update_by_reals"] = False
    pst.control_data.noptmax = 1
    pst.write(os.path.join(template_d, "pest_high_phi.pst"))
    #pyemu.os_utils.start_workers(template_d, exe_path, "pest_high_phi.pst", num_workers=10,
    #                             master_dir=test_d, verbose=True, worker_root=model_d,
    #                             port=port)
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pyemu.os_utils.run("{0} {1}".format(exe_path,"pest_high_phi.pst"),cwd=test_d)
    phi3 = pd.read_csv(os.path.join(test_d, "pest_high_phi.phi.actual.csv"), index_col=0)
    diff = phi3 - phi2
    assert diff.max().max() == 0.0

    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options["ies_lambda_mults"] = [0.5, 1.0]
    pst.pestpp_options["lambda_scale_fac"] = [0.9, 1.0]
    pst.pestpp_options['ies_subset_size'] = 3
    pst.pestpp_options["ies_debug_high_upgrade_phi"] = True
    pst.pestpp_options["ies_debug_fail_subset"] = True
    pst.pestpp_options["ies_debug_fail_remainder"] = True
    pst.pestpp_options["ies_debug_bad_phi"] = True
    pst.control_data.noptmax = 3
    pst.write(os.path.join(template_d, "pest_high_phi.pst"))
    #pyemu.os_utils.start_workers(template_d, exe_path, "pest_high_phi.pst", num_workers=10,
    #                             master_dir=test_d, verbose=True, worker_root=model_d,
    #                             port=port)
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pyemu.os_utils.run("{0} {1}".format(exe_path,"pest_high_phi.pst"),cwd=test_d)
    phi4 = pd.read_csv(os.path.join(test_d, "pest_high_phi.phi.actual.csv"), index_col=0)
    assert os.path.exists(os.path.join(test_d, "pest_high_phi.3.obs.csv"))

    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options["ies_lambda_mults"] = [0.5, 1.0]
    pst.pestpp_options["lambda_scale_fac"] = [0.9, 1.0]
    pst.pestpp_options['ies_subset_size'] = 3
    pst.pestpp_options["ies_debug_high_subset_phi"] = True
    pst.pestpp_options["ies_debug_fail_subset"] = True
    pst.pestpp_options["ies_debug_fail_remainder"] = True
    pst.pestpp_options["ies_debug_bad_phi"] = True
    pst.pestpp_options["ies_center_on"] = "base"
    pst.control_data.noptmax = 3
    pst.write(os.path.join(template_d, "pest_high_phi.pst"))
    #pyemu.os_utils.start_workers(template_d, exe_path, "pest_high_phi.pst", num_workers=10,
    #                             master_dir=test_d, verbose=True, worker_root=model_d,
    #                             port=port)
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pyemu.os_utils.run("{0} {1}".format(exe_path,"pest_high_phi.pst"),cwd=test_d)
    phi5 = pd.read_csv(os.path.join(test_d, "pest_high_phi.phi.actual.csv"), index_col=0)
    assert os.path.exists(os.path.join(test_d, "pest_high_phi.3.obs.csv"))
    


    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options["ies_lambda_mults"] = [0.5, 1.0]
    pst.pestpp_options["lambda_scale_fac"] = [0.9, 1.0]
    pst.pestpp_options['ies_subset_size'] = 3
    #pst.pestpp_options["ies_debug_high_subset_phi"] = True
    #pst.pestpp_options["ies_debug_fail_subset"] = True
    #pst.pestpp_options["ies_debug_fail_remainder"] = True
    #pst.pestpp_options["ies_debug_bad_phi"] = True
    pst.pestpp_options["ies_debug_high_upgrade_phi"] = True
    pst.pestpp_options["ies_center_on"] = "base"
    pst.pestpp_options["ies_no_noise"] = False
    pst.control_data.noptmax = 3
    pst.write(os.path.join(template_d, "pest_high_phi.pst"))
    #pyemu.os_utils.start_workers(template_d, exe_path, "pest_high_phi.pst", num_workers=10,
    #                             master_dir=test_d, verbose=True, worker_root=model_d,
    #                             port=port)
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pyemu.os_utils.run("{0} {1}".format(exe_path,"pest_high_phi.pst"),cwd=test_d)
    phi5 = pd.read_csv(os.path.join(test_d, "pest_high_phi.phi.actual.csv"), index_col=0)
    
    assert os.path.exists(os.path.join(test_d, "pest_high_phi.3.obs.csv"))
    for i in range(2):
        rpe_file = os.path.join(test_d,"pest_high_phi.rejected.{0}.par.csv".format(i+1))
        roe_file = os.path.join(test_d,"pest_high_phi.rejected.{0}.obs.csv".format(i+1))
        assert os.path.exists(rpe_file),rpe_file
        assert os.path.exists(roe_file),roe_file
        roe = pd.read_csv(roe_file,index_col=0)
        rpe = pd.read_csv(rpe_file,index_col=0)
        assert roe.shape[0] == rpe.shape[0]

def freyberg_svd_draws_invest():
    import flopy
    model_d = "ies_freyberg"
    test_d = os.path.join(model_d, "master_svd_draws_test")
    template_d = os.path.join(model_d, "template")
    m = flopy.modflow.Modflow.load("freyberg.nam", model_ws=template_d, load_only=[], check=False)
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    # print("loading pst")
    pst = pyemu.Pst(os.path.join(template_d, "pest.pst"))
    par = pst.parameter_data

    print(par.pargp.value_counts().sort_values())
    pnames = par.loc[par.pargp == "hk", "parnme"].tolist()
    par.loc[par.pargp != "hk", "partrans"] = "fixed"
    cov = pyemu.Cov.from_binary(os.path.join(template_d, "prior.jcb"))
    pcov = cov.get(pnames, pnames)
    ev, ew = np.linalg.eigh(pcov.as_2d)
    u, s, v = np.linalg.svd(pcov.as_2d, full_matrices=True)
    print(np.sqrt(s))
    print(ev)

    s = np.diag(s)
    eproj = np.dot(ew, np.sqrt(np.diag(ev)))
    sproj = np.dot(u, s)
    # print(eproj.shape,sproj.shape,s.shape)
    diff = u - ew

    # pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst,cov=cov,num_reals=10000)
    # pe = pe.loc[:,pnames]
    # pmat = pe.as_pyemu_matrix()
    # pmat =pmat.get(col_names=pnames)
    # ecov = pmat.T * pmat
    # print(ecov.shape)


def freyberg_center_on_test():
    import flopy
    model_d = "ies_freyberg"
    test_d = os.path.join(model_d, "master_center_on1")
    template_d = os.path.join(model_d, "template")
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    # print("loading pst")
    pst = pyemu.Pst(os.path.join(template_d, "pest.pst"))
    pst.pestpp_options = {"ies_num_reals": 5}
    pst.pestpp_options["ies_lambda_mults"] = 1.0
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["ies_subset_size"] = 10
    pst.control_data.nphinored = 20
    pst.control_data.noptmax = 6
    pst.write(os.path.join(template_d, "pest_base.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_base.pst", num_workers=5, master_dir=test_d,
                                 worker_root=model_d, port=port)

    base_phi = pd.read_csv(os.path.join(test_d, "pest_base.phi.actual.csv"), index_col=0)
    shutil.copy2(os.path.join(test_d, "pest_base.0.par.csv"), os.path.join(template_d, "par.csv"))
    shutil.copy2(os.path.join(test_d, "pest_base.0.obs.csv"), os.path.join(template_d, "obs.csv"))
    # shutil.copy2(os.path.join(test_d,"pest_base.obs+noise.csv"),os.path.join(template_d,"noise.csv"))
    pst.pestpp_options["ies_par_en"] = "par.csv"
    # pst.pestpp_options["ies_obs_en"] = "noise.csv"
    pst.pestpp_options["ies_restart_obs_en"] = "obs.csv"

    pst.pestpp_options["ies_center_on"] = "base"
    pst.write(os.path.join(template_d, "pest_center_on.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_center_on.pst", num_workers=5, master_dir=test_d,
                                 worker_root=model_d, port=port)
    center_phi = pd.read_csv(os.path.join(test_d, "pest_center_on.phi.actual.csv"), index_col=0)
    print(base_phi.loc[:, "base"])
    print(center_phi.loc[:, "base"])

    pst.pestpp_options["ies_center_on"] = "_median_"
    pst.write(os.path.join(template_d, "pest_center_on.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_center_on.pst", num_workers=5,
                                 master_dir=test_d + "_median",
                                 worker_root=model_d, port=port)
    center_phi = pd.read_csv(os.path.join(test_d, "pest_center_on.phi.actual.csv"), index_col=0)

    # assert center_phi.loc[pst.control_data.noptmax,"base"] < base_phi.loc[pst.control_data.noptmax,"base"]


def freyberg_pdc_test():
    import flopy
    model_d = "ies_freyberg"
    test_d = os.path.join(model_d, "master_pdc")
    template_d = os.path.join(model_d, "template")
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    # print("loading pst")
    pst = pyemu.Pst(os.path.join(template_d, "pest.pst"))
    pst.observation_data.loc[pst.nnz_obs_names[0], "obsval"] += 20
    pst.pestpp_options = {"ies_num_reals": 5}
    pst.pestpp_options["ies_lambda_mults"] = 1.0
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["ies_subset_size"] = 10
    pst.pestpp_options["ies_drop_conflicts"] = True
    pst.pestpp_options["ies_autoadaloc"] = True
    pst.control_data.nphinored = 20
    pst.control_data.noptmax = -1
    pst.write(os.path.join(template_d, "pest_base.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_base.pst", num_workers=5, master_dir=test_d,
                                 worker_root=model_d, port=port)
    return
    phi_csv = os.path.join(test_d, "pest_base.phi.actual.csv")
    assert os.path.exists(phi_csv), phi_csv
    pdc_phi = pd.read_csv(phi_csv, index_col=0)
    assert pdc_phi.shape[0] == pst.control_data.noptmax + 1
    # scan the rec file for the conflicted obs names
    dropped = []
    with open(os.path.join(test_d, "pest_base.rec"), 'r') as f:
        while True:
            line = f.readline()
            if line == "":
                raise Exception()
            if "...conflicted observations:" in line:
                while True:
                    line = f.readline()
                    if line == "":
                        raise Exception()
                    if line.startswith("...dropping"):
                        break
                    dropped.append(line.strip().lower())
                break
    print(dropped)
    shutil.copy2(os.path.join(test_d, "pest_base.0.par.csv"), os.path.join(template_d, "pdc_par.csv"))
    pst.pestpp_options["ies_par_en"] = "pdc_par.csv"
    shutil.copy2(os.path.join(test_d, "pest_base.obs+noise.csv"), os.path.join(template_d, "pdc_obs.csv"))
    pst.pestpp_options["ies_obs_en"] = "pdc_obs.csv"

    pst.observation_data.loc[dropped, "weight"] = 0.0
    pst.pestpp_options["ies_num_reals"] = 10
    pst.write(os.path.join(template_d, "pest_base.pst"))
    test_d = os.path.join(model_d, "master_pdc_base")
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_base.pst", num_workers=5, master_dir=test_d,
                                 worker_root=model_d, port=port)
    phi_csv = os.path.join(test_d, "pest_base.phi.actual.csv")
    assert os.path.exists(phi_csv), phi_csv
    base_phi = pd.read_csv(phi_csv, index_col=0)
    assert base_phi.shape[0] == pst.control_data.noptmax + 1
    diff = (pdc_phi - base_phi).apply(lambda x: np.abs(x))
    print(diff.max())
    assert diff.max().max() < 0.1, diff.max().max()

    pst.pestpp_options["ies_pdc_sigma_distance"] = 1.0
    pst.write(os.path.join(template_d, "pest_pdc_dist.pst"))
    test_d = os.path.join(model_d, "master_pdc_dist")
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_pdc_dist.pst", num_workers=5, master_dir=test_d,
                                 worker_root=model_d, port=port)

    oe = pd.read_csv(os.path.join(test_d, "pest_pdc_dist.0.obs.csv"), index_col=0)
    oe_base = pd.read_csv(os.path.join(test_d, "pest_pdc_dist.obs+noise.csv"), index_col=0)
    smn, sstd = oe.mean(), oe.std()
    omn, ostd = oe_base.mean(), oe_base.std()
    for name in oe.columns:
        if name not in pst.nnz_obs_names:
            continue
        # print(name,smn[name],sstd[name],omn[name],ostd[name])
    smin = smn - sstd
    smax = smn + sstd
    omin = omn - ostd
    omax = omn + ostd
    conflict = []
    for name, omnn, omx, smnn, smx in zip(oe.columns.values, omin, omax, smin, smax):
        if name not in pst.nnz_obs_names:
            continue
        print(name, smn[name], sstd[name], smnn, smx,
              omn[name], ostd[name], omnn, omx)
        if omx < smnn or omnn > smx:
            conflict.append(name)
    print(conflict)

    pst.pestpp_options["ies_no_noise"] = True
    pst.pestpp_options.pop("ies_obs_en")
    pst.write(os.path.join(template_d, "pest_pdc_dist.pst"))
    test_d = os.path.join(model_d, "master_pdc_dist")
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_pdc_dist.pst", num_workers=5, master_dir=test_d,
                                 worker_root=model_d, port=port)

    oe = pd.read_csv(os.path.join(test_d, "pest_pdc_dist.0.obs.csv"), index_col=0)
    oe_base = pd.read_csv(os.path.join(test_d, "pest_pdc_dist.obs+noise.csv"), index_col=0)
    smn, sstd = oe.mean(), oe.std()
    omn, ostd = oe_base.mean(), oe_base.std()
    for name in oe.columns:
        if name not in pst.nnz_obs_names:
            continue
        # print(name,smn[name],sstd[name],omn[name],ostd[name])
    smin = smn - sstd
    smax = smn + sstd
    omin = omn - ostd
    omax = omn + ostd
    conflict = []
    for name, omnn, omx, smnn, smx in zip(oe.columns.values, omin, omax, smin, smax):
        if name not in pst.nnz_obs_names:
            continue
        print(name, smn[name], sstd[name], smnn, smx,
              omn[name], ostd[name], omnn, omx)
        if omx < smnn or omnn > smx:
            conflict.append(name)
    print(conflict)


def freyberg_rcov_test():
    import flopy
    model_d = "ies_freyberg"
    test_d = os.path.join(model_d, "master_rcov")
    template_d = os.path.join(model_d, "template")
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    # print("loading pst")
    pst = pyemu.Pst(os.path.join(template_d, "pest.pst"))
    pst.observation_data.loc[pst.nnz_obs_names[0], "obsval"] += 20
    pst.pestpp_options = {"ies_num_reals": 8}
    pst.pestpp_options["ies_debug_fail_remainder"] = True
    # pst.pestpp_options["ies_lambda_mults"] = 1.0
    # pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["ies_subset_size"] = 3
    pst.pestpp_options["ies_drop_conflicts"] = True
    pst.pestpp_options["ies_autoadaloc"] = True
    pst.pestpp_options["ies_save_rescov"] = True
    pst.pestpp_options["ies_verbose_level"] = 1
    pst.control_data.nphinored = 20
    pst.control_data.noptmax = 2
    pst.write(os.path.join(template_d, "pest_rescov.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_rescov.pst", num_workers=8, master_dir=test_d,
                                 worker_root=model_d, port=port)
    # check that the shrunk res cov has the same diag as the org res cov
    org_rescov = pyemu.Cov.from_ascii(os.path.join(test_d, "pest_rescov.2.res.cov"))
    shrunk_rescov = pyemu.Cov.from_ascii(os.path.join(test_d, "pest_rescov.2.shrunk_res.cov"))
    diff = np.abs(np.diag(org_rescov.x) - np.diag(shrunk_rescov.x))
    print(diff)
    assert diff.sum() < 1.0e-6, diff.sum()
    shutil.copy2(os.path.join(test_d, "pest_rescov.2.shrunk_res.cov"), os.path.join(template_d, "post_obs.cov"))
    pst.pestpp_options["obscov"] = "post_obs.cov"
    pst.pestpp_options["ies_drop_conflicts"] = False
    
    pst.write(os.path.join(template_d, "pest_bmw.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_bmw.pst", num_workers=8, master_dir=test_d,
                                 worker_root=model_d, port=port)
    org_rescov = pyemu.Cov.from_ascii(os.path.join(test_d, "pest_bmw.2.res.cov"))
    shrunk_rescov = pyemu.Cov.from_ascii(os.path.join(test_d, "pest_bmw.2.shrunk_res.cov"))
    diff = np.abs(np.diag(org_rescov.x) - np.diag(shrunk_rescov.x))
    print(diff)
    assert diff.sum() < 1.0e-6, diff.sum()


def tenpar_align_test():
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_align_test")
    template_d = os.path.join(model_d, "test_template")

    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)

    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d,test_d)
    np.random.seed(1234)
    oe = pyemu.ObservationEnsemble.from_gaussian_draw(pst, num_reals=9)
    oe.loc["base", pst.nnz_obs_names] = pst.observation_data.loc[pst.nnz_obs_names, "obsval"]
    oe.loc[:, "new_index"] = list(oe.index.map(lambda x: str(x)))
    oe.set_index("new_index", inplace=True)
    oe.sort_index(inplace=True, ascending=False)
    oe.to_csv(os.path.join(template_d, "out_of_order_oe.csv"))
    pst.pestpp_options["ies_obs_en"] = "out_of_order_oe.csv"
    pst.pestpp_options["ies_debug_fail_remainder"] = True
    pst.pestpp_options["ies_num_reals"] = 10
    pst.control_data.noptmax = 1
    pst_name = "pest_align.pst"
    pst.write(os.path.join(template_d, pst_name))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d,test_d)
    pyemu.os_utils.run("{0} {1}".format(exe_path,pst_name),cwd=test_d)
    #pyemu.os_utils.start_workers(template_d, exe_path, pst_name, num_workers=8, master_dir=test_d,
    #                             worker_root=model_d, port=port)
    pe_file = os.path.join(test_d, pst_name.replace(".pst", ".1.par.csv"))
    oe_file = os.path.join(test_d, pst_name.replace(".pst", ".1.obs.csv"))
    assert os.path.exists(pe_file), pe_file
    assert os.path.exists(oe_file), oe_file
    pe = pd.read_csv(pe_file)
    oe1 = pd.read_csv(oe_file)
    print(pe.columns)
    print(oe1.columns)
    for i in pe.index:
        pr = pe.loc[i, "real_name"]
        or1 = oe1.loc[i, "real_name"]
        print(pr, or1)
        if (pr != or1):
            raise Exception("real names differ " + pr + "," + or1)
    oe1.index = oe1.pop("real_name")
    oe1.to_csv(os.path.join(template_d, "align_obs_restart.csv"))
    shutil.copy2(os.path.join(test_d, pst_name.replace(".pst", ".0.par.csv")),
                 os.path.join(template_d, "align_par.csv"))
    pst.pestpp_options["ies_restart_obs_en"] = "align_obs_restart.csv"
    pst.pestpp_options["ies_par_en"] = "align_par.csv"
    pst.write(os.path.join(template_d, pst_name))
    # pyemu.os_utils.run("{0} {1}".format(exe_path,pst_name),cwd=test_d)
    #pyemu.os_utils.start_workers(template_d, exe_path, pst_name, num_workers=8, master_dir=test_d,
    #                             worker_root=model_d, port=port)
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d,test_d)
    pyemu.os_utils.run("{0} {1}".format(exe_path,pst_name),cwd=test_d)
    pe_file = os.path.join(test_d, pst_name.replace(".pst", ".1.par.csv"))
    oe_file = os.path.join(test_d, pst_name.replace(".pst", ".1.obs.csv"))
    assert os.path.exists(pe_file), pe_file
    assert os.path.exists(oe_file), oe_file
    pe = pd.read_csv(pe_file)
    oe1 = pd.read_csv(oe_file)
    print(pe.columns)
    print(oe1.columns)
    for i in pe.index:
        pr = pe.loc[i, "real_name"]
        or1 = oe1.loc[i, "real_name"]
        print(pr, or1)
        if (pr != or1):
            raise Exception("real names differ " + pr + "," + or1)


def tenpar_align_test_2():
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_align_test_2")
    template_d = os.path.join(model_d, "test_template")

    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)

    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pst.pestpp_options["ies_num_reals"] = 10
    pst.control_data.noptmax = 2
    pst_name = "pest_align.pst"
    pst.write(os.path.join(test_d, pst_name))
    pyemu.os_utils.run("{0} {1}".format(exe_path, pst_name), cwd=test_d)
    shutil.copy2(os.path.join(test_d, pst_name.replace(".pst", ".0.par.csv")), os.path.join(test_d, "restart.csv"))
    pst.pestpp_options["ies_par_en"] = "restart.csv"
    pst.write(os.path.join(test_d, pst_name))
    pyemu.os_utils.run("{0} {1}".format(exe_path, pst_name), cwd=test_d)

    pe = pd.read_csv(os.path.join(test_d, pst_name.replace(".pst", ".0.par.csv")), index_col=0)
    oe = pd.read_csv(os.path.join(test_d, pst_name.replace(".pst", ".0.obs.csv")), index_col=0)
    assert pe.shape[0] == oe.shape[0], "shape mismatch: {0} vs {1}".format(pe.shape[0], oe.shape[0])
    for i, (p, o) in enumerate(zip(pe.index, oe.index)):
        if p != o:
            raise Exception("misaligned indices in row {0}: {1} vs {2}".format(i, p, o))




def tenpar_upgrade_on_disk_test():
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_upgrade_1")
    template_d = os.path.join(model_d, "test_template")

    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)

    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options["ies_no_noise"] = True
    pst.pestpp_options["ies_lambda_mults"] = [0.1, 1.0]
    pst.pestpp_options["lambda_scale_fac"] = [0.7, 1.0]
    pst.pestpp_options["ies_accept_phi_fac"] = 1000

    par = pst.parameter_data
    par.loc[:,"partrans"] = "log"

    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst,cov=pyemu.Cov.from_parameter_data(pst),num_reals=10)
    
    fixed = pst.par_names[0:2]
    par.loc[fixed,"partrans"] = "fixed"
    adj= pst.adj_par_names
    pe = pe.loc[:,adj]
    print(pe)
    pe.to_binary(os.path.join(template_d,"prior.jcb"))
    
    pst.pestpp_options["ies_par_en"] = "prior.jcb"
    pst.control_data.noptmax = 2
    pst_name = "pest_upgrade.pst"
    pst.write(os.path.join(template_d, pst_name))
    # pyemu.os_utils.run("{0} {1}".format(exe_path, pst_name), cwd=test_d)
    #pyemu.os_utils.start_workers(template_d, exe_path, pst_name, num_workers=8,
    #                             master_dir=test_d, worker_root=model_d, port=port)
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pyemu.os_utils.run("{0} {1}".format(exe_path,pst_name),cwd=test_d)
    phi1 = pd.read_csv(os.path.join(test_d, pst_name.replace(".pst", ".phi.actual.csv")))
    pe1 = pd.read_csv(os.path.join(test_d, pst_name.replace(".pst", ".{0}.par.csv". \
                                                            format(pst.control_data.noptmax))), index_col=0)
    oe1 = pd.read_csv(os.path.join(test_d, pst_name.replace(".pst", ".{0}.obs.csv". \
                                                            format(pst.control_data.noptmax))), index_col=0)

    pst.pestpp_options["ies_upgrades_in_memory"] = False
    pst.write(os.path.join(template_d, pst_name))
    test_d = os.path.join(model_d, "master_upgrade_2")
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)

    #pyemu.os_utils.start_workers(template_d, exe_path, pst_name, num_workers=8,
    #                             master_dir=test_d, worker_root=model_d, port=port)
    pyemu.os_utils.run("{0} {1}".format(exe_path,pst_name),cwd=test_d)
    
    phi2 = pd.read_csv(os.path.join(test_d, pst_name.replace(".pst", ".phi.actual.csv")))
    pe2 = pd.read_csv(os.path.join(test_d, pst_name.replace(".pst", ".{0}.par.csv". \
                                                            format(pst.control_data.noptmax))), index_col=0)
    oe2 = pd.read_csv(os.path.join(test_d, pst_name.replace(".pst", ".{0}.obs.csv". \
                                                            format(pst.control_data.noptmax))), index_col=0)

    d = (phi1 - phi2).apply(lambda x: np.abs(x))
    print(d.max())
    print(d.max().max())
    assert d.max().max() < 1.0e-6

    d = (pe1 - pe2).apply(lambda x: np.abs(x))
    print(d.max())
    print(d.max().max())
    assert d.max().max() < 1.0e-6

    d = (oe1 - oe2).apply(lambda x: np.abs(x))
    print(d.max())
    print(d.max().max())
    assert d.max().max() < 1.0e-6


def multimodal_test():
    noptmax = 2
    num_reals = 100
    # can be "circle" or "h"
    func = "circle"
    model_d = "mm1"
    test_d = os.path.join(model_d, "template")
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    os.makedirs(test_d)
    tpl_file = os.path.join(test_d, "par.dat.tpl")
    with open(tpl_file, 'w') as f:
        f.write("ptf ~\n")
        f.write("par1  ~   par1   ~\n")
        f.write("par2  ~   par2   ~\n")
    ins_file = os.path.join(test_d, "obs.dat.ins")
    with open(ins_file, 'w') as f:
        f.write("pif ~\n")
        f.write("l1 !obs1!\n")
        f.write("l1 !par1!\n")
        f.write("l1 !par2!\n")
        
        
    with open(os.path.join(test_d, "run.py"), 'w') as f:
        f.write("import numpy as np\n")
        f.write("lines =  open('par.dat','r').readlines()\n")
        if func == "circle":
            f.write("p1 = float(lines[0].strip().split()[-1])\n")
            f.write("p2 = float(lines[1].strip().split()[-1])\n")
            f.write("result =  np.sqrt(p1**2 + p2**2)\n")
        else:
            f.write("p1 = float(lines[0].strip().split()[-1])\n")
            f.write("p2 = float(lines[1].strip().split()[-1])\n")
            f.write("result = ((p1**2 + p2 - 11)**2 + (p1 + p2**2 -7)**2)\n")
        f.write("with open('obs.dat','w') as f:\n")
        f.write("    f.write('{0:15.6E}\\n'.format(result))\n")
        f.write("    f.write('{0:15.6E}\\n'.format(p1))\n")
        f.write("    f.write('{0:15.6E}\\n'.format(p2))\n")
        
    pst = pyemu.Pst.from_io_files(tpl_file, tpl_file.replace(".tpl", ""), ins_file, ins_file.replace(".ins", ""),
                                  pst_path=".")
    pst.model_command = "python run.py"
    pst.parameter_data.loc[:, "partrans"] = "none"
    pst.parameter_data.loc[:, "parval1"] = 0
    if func == "circle":
        pst.parameter_data.loc[:, "parubnd"] = 2
        pst.parameter_data.loc[:, "parlbnd"] = -2

        pst.observation_data.loc[:, "obsval"] = 1
    else:
        pst.parameter_data.loc[:, "parubnd"] = 5
        pst.parameter_data.loc[:, "parlbnd"] = -5
        pst.observation_data.loc[:, "obsval"] = 0

    pst.parameter_data.loc[:, "parchglim"] = "relative"

    pst.observation_data.loc[:, "weight"] = 0
    pst.observation_data.loc["obs1", "weight"] = 100.0
    pst.observation_data.loc[:,"obgnme"] = pst.observation_data.obsnme.values
 

    pst.control_data.noptmax = 0
    pst.write(os.path.join(test_d, "mm1.pst"))
    pyemu.os_utils.run("{0} mm1.pst".format(exe_path), cwd=test_d)
    pst.control_data.noptmax = noptmax
    pst.pestpp_options["ies_num_reals"] = num_reals
    #pst.pestpp_options["ies_lambda_mults"] = 1.0
    #pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["ies_subset_size"] = -10
    pst.pestpp_options["ies_multimodal_alpha"] = 0.1
    pst.pestpp_options["ies_verbose_level"] = 3
    pst.pestpp_options["ies_include_base"] = False
    pst.pestpp_options["ies_use_approx"] = True
    pst.pestpp_options["ies_save_lambda_en"] = True
    pst.pestpp_options["ies_num_threads"] = 3
    pst.pestpp_options["ies_debug_bad_phi"] = True
    pst.pestpp_options["ies_debug_fail_remainder"] = True
    pst.pestpp_options["ies_debug_fail_subset"] = True


    #pst.pestpp_options["ies_bad_phi_sigma"] = 1.25
    pst.pestpp_options["ies_use_mda"] = False


    #pst.pestpp_options["ies_no_noise"] = True
    pst.write(os.path.join(test_d, "mm1.pst"))

    if False:
        par = pst.parameter_data
        p1,p2 = par.parnme.iloc[0],par.parnme.iloc[1]
        runs = {p1:[],p2:[]}
        steps = 100
        for i in np.linspace(par.parlbnd.iloc[0],par.parubnd.iloc[0],steps):
            for j in np.linspace(par.parlbnd.iloc[1],par.parubnd.iloc[1],steps):
                runs[p1].append(i)
                runs[p2].append(j)

        df = pd.DataFrame(runs)
        df.to_csv(os.path.join(test_d,"sweep_in.csv"))
        m_d = os.path.join(model_d, "master_sweep_{0}".format(func))
        pyemu.os_utils.start_workers(test_d, exe_path.replace("ies","swp"), "mm1.pst", worker_root=model_d, num_workers=35, master_dir=m_d)
    

    m_d = os.path.join(model_d, "master_mm_{0}_mt".format(func))
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    shutil.copytree(test_d,m_d)
    pyemu.os_utils.run("{0} mm1.pst".format(exe_path),cwd=m_d)
    #pyemu.os_utils.start_workers(test_d, exe_path, "mm1.pst", worker_root=model_d, num_workers=35, master_dir=m_d)
    

    pst.pestpp_options["ies_num_threads"] = 1
    pst.write(os.path.join(test_d, "mm1.pst"))
    m_d = os.path.join(model_d, "master_mm_{0}_single".format(func))
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    shutil.copytree(test_d,m_d)
    pyemu.os_utils.run("{0} mm1.pst".format(exe_path),cwd=m_d)
    
    #pyemu.os_utils.start_workers(test_d, exe_path, "mm1.pst", worker_root=model_d, num_workers=35, master_dir=m_d)
    phidf = pd.read_csv(os.path.join(m_d,"mm1.phi.actual.csv"))
    assert phidf.iteration.max() == pst.control_data.noptmax
    assert phidf.loc[phidf.index[-1],"min"] < 0.1

    pst.pestpp_options["ies_multimodal_alpha"] = 1.0
    pst.write(os.path.join(test_d, "mm1.pst"))
    m_d = os.path.join(model_d, "master_base_{0}".format(func))
    #pyemu.os_utils.start_workers(test_d, exe_path, "mm1.pst", worker_root=model_d, num_workers=35, master_dir=m_d)
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    shutil.copytree(test_d,m_d)
    pyemu.os_utils.run("{0} mm1.pst".format(exe_path),cwd=m_d)
    
    phidf = pd.read_csv(os.path.join(m_d,"mm1.phi.actual.csv"))
    assert phidf.iteration.max() == pst.control_data.noptmax
    
    pst.pestpp_options["ies_multimodal_alpha"] = 0.1
    pst.pestpp_options["ies_num_threads"] = 4
    pst.pestpp_options["ies_include_base"] = True
    pst.pestpp_options["ies_center_on"] = "base"

    pst.write(os.path.join(test_d, "mm1.pst"))
    m_d = os.path.join(model_d, "master_mm_centeron_{0}".format(func))
    #pyemu.os_utils.start_workers(test_d, exe_path, "mm1.pst", worker_root=model_d, num_workers=35, master_dir=m_d)
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    shutil.copytree(test_d,m_d)
    pyemu.os_utils.run("{0} mm1.pst".format(exe_path),cwd=m_d)
    
    phidf = pd.read_csv(os.path.join(m_d,"mm1.phi.actual.csv"))
    assert phidf.iteration.max() == pst.control_data.noptmax
    assert phidf.loc[phidf.index[-1],"min"] < 0.1

    pst.observation_data.loc["par1", "weight"] = 100.0
    pst.observation_data.loc["par1","obsval"] = pst.parameter_data.loc["par1","parubnd"]
    pst.observation_data.loc["obs1", "weight"] = 100.0

    with open(os.path.join(test_d,"phi_facs.dat"),'w') as f:
        f.write("obs1 0.5\n")
        f.write("par1 0.5\n")
    pst.pestpp_options["ies_phi_factor_file"] = "phi_facs.dat"
    pst.write(os.path.join(test_d, "mm1.pst"))
    m_d = os.path.join(model_d, "master_mm_centeron_phifac_{0}".format(func))
    #pyemu.os_utils.start_workers(test_d, exe_path, "mm1.pst", worker_root=model_d, num_workers=35, master_dir=m_d)
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    shutil.copytree(test_d,m_d)
    pyemu.os_utils.run("{0} mm1.pst".format(exe_path),cwd=m_d)
    
    phidf = pd.read_csv(os.path.join(m_d,"mm1.phi.actual.csv"))
    assert phidf.iteration.max() == pst.control_data.noptmax
    
    df = pd.DataFrame(index=np.arange(num_reals),columns=["par1","obs1"])
    df.loc[:,:] = 100
    df.iloc[:int(df.shape[0]/3),0] = 1e-10
    df.iloc[-int(df.shape[0]/3):,1] = 1e-10
    df.to_csv(os.path.join(test_d,"phi_fac.csv"))
    pst.pestpp_options["ies_phi_factor_file"] = "phi_fac.csv"
    pst.pestpp_options["ies_phi_factors_by_real"] = True

    pst.write(os.path.join(test_d, "mm1.pst"))
    m_d = os.path.join(model_d, "master_mm_phifac_byreal_{0}".format(func))
    #pyemu.os_utils.start_workers(test_d, exe_path, "mm1.pst", worker_root=model_d, num_workers=35, master_dir=m_d)
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    shutil.copytree(test_d,m_d)
    pyemu.os_utils.run("{0} mm1.pst".format(exe_path),cwd=m_d)
    
    phidf = pd.read_csv(os.path.join(m_d,"mm1.phi.actual.csv"))
    assert phidf.iteration.max() == pst.control_data.noptmax
    assert phidf.loc[phidf.index[-1],"min"] < 0.1


def plot_mm1_sweep_results():

    import matplotlib.pyplot as plt

    mm_d = os.path.join("mm1", "master_sweep_circle")
    pst = pyemu.Pst(os.path.join(mm_d, "mm1.pst"))
    pe = pd.read_csv(os.path.join(mm_d,"sweep_in.csv"))
    oe = pd.read_csv(os.path.join(mm_d,"sweep_out.csv"))
    pe.loc[:,"obj"] = oe.phi.values
    X = pe.par1.values.reshape(100,100)
    Y = pe.par2.values.reshape(100,100)
    Z = oe.phi.values.reshape(100,100)

    

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    surf = ax.plot_surface(X, Y, Z, cmap="jet",
                       linewidth=0, antialiased=False,alpha=0.75)
    ax.set_xlabel("par1")
    ax.set_ylabel("par2")
    ax.set_zlabel("phi")
    plt.show()




def plot_mm1_results(noptmax=None, func="circle", show_info=False,mm_d = None):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle

    base_d = os.path.join("mm1", "master_base_{0}".format(func))
    if mm_d is None:
        mm_d = os.path.join("mm1", "master_mm_{0}_mt".format(func))
    pst = pyemu.Pst(os.path.join(base_d, "mm1.pst"))
    if noptmax is None:
        noptmax = pst.control_data.noptmax
    org_noptmax = noptmax + 1
    for noptmax in range(org_noptmax):
        fname = os.path.join(base_d, "mm1.{0}.par.csv".format(noptmax))
        for i in range(noptmax):
            if os.path.exists(fname):
                break
            fname = os.path.join(base_d, "mm1.{0}.par.csv".format(noptmax - i))

        pe_pt_base = pd.read_csv(fname)
        oe_pt_base = pyemu.ObservationEnsemble(pst=pst, df=pd.read_csv(fname.replace(".par.", ".obs."), index_col=0))
        base_pv = oe_pt_base.phi_vector

        fname = os.path.join(mm_d, "mm1.{0}.par.csv".format(noptmax))
        for i in range(noptmax):
            if os.path.exists(fname):
                break
            fname = os.path.join(mm_d, "mm1.{0}.par.csv".format(noptmax - i))
        pe_pt_mm = pd.read_csv(fname)
        oe_pt_mm = pyemu.ObservationEnsemble.from_dataframe(pst=pst,
                                                            df=pd.read_csv(fname.replace(".par.", ".obs."), index_col=0))
        mm_pv = oe_pt_mm.phi_vector
        pe_pr = pd.read_csv(os.path.join(mm_d, "mm1.0.par.csv"))
        pe_pt_mm.index = pe_pt_mm.index.map(lambda x: str(int(np.float(x))))
        pe_pr.index = pe_pr.index.map(lambda x: str(int(np.float(x))))
        pe_pr.index = pe_pr.index.map(lambda x: str(x))

        if show_info and noptmax > 0:
            mm_info_fname = [f for f in os.listdir(mm_d) if "mm1.{0}.".format(noptmax) in f and f.endswith(".mm.info.csv")][0]
            print(mm_info_fname)
            mm_df = pd.read_csv(os.path.join(mm_d, mm_info_fname))
            mm_df.index = np.arange(mm_df.shape[0])
            mm_df.loc[:, "pe_real_name"] = mm_df.pe_real_name.apply(lambda x: str(x))
            nei_cols = mm_df.columns[mm_df.columns.map(lambda x: "neighbor" in x)]
            mm_rnames = set(pe_pr.index.tolist())
            df = mm_df.iloc[1, :]
            print(df)
            mm_rname = df.pe_real_name
            neis = df.loc[nei_cols].apply(str)
            neis = neis.loc[neis.apply(lambda x: x in mm_rnames)]
            df = df.loc[neis.index]
            print(neis)
            # print(pe_pt_mm.index)
            pe_pr_nei = pe_pr.loc[neis, :]
            print(mm_rname, pe_pr.loc[mm_rname, :], pe_pt_mm.loc[mm_rname, :])

        fig, axes = plt.subplots(1, 2, figsize=(6, 3))

        axes[1].scatter(pe_pr.par1.values, pe_pr.par2.values, marker=".", color="0.5", alpha=0.5)
        if noptmax != 0:
            axes[1].scatter(pe_pt_mm.par1.values, pe_pt_mm.par2.values, marker=".", c="b", alpha=0.75)
        if show_info and noptmax > 0:
            axes[1].scatter(pe_pr_nei.par1.values, pe_pr_nei.par2.values, marker=".", color="c", s=100)
            #for rname, p1, p2 in zip(pe_pr_nei.index, pe_pr_nei.par1.values, pe_pr_nei.par2.values):
            #    axes[1].text(p1 + 0.025, p2, rname)

            axes[1].scatter([pe_pr.loc[mm_rname, "par1"]], [pe_pr.loc[mm_rname, "par2"]], marker=".", color="r", s=200)
            if noptmax != 0:
                axes[1].scatter([pe_pt_mm.loc[mm_rname, "par1"]], [pe_pt_mm.loc[mm_rname, "par2"]], marker=".", color="m",
                                s=200)

        axes[1].set_title("multimodal upgrade")
        axes[1].set_ylabel("par2")
        axes[1].set_xlabel("par1")

        c = Circle([0,0],1,edgecolor="r",facecolor="none")
        axes[1].add_patch(c)
        
        axes[0].scatter(pe_pr.par1.values, pe_pr.par2.values, marker=".", color="0.5", alpha=0.5)
        if noptmax != 0:
            axes[0].scatter(pe_pt_base.par1.values, pe_pt_base.par2.values, marker=".", c="b", alpha=0.5)

        c = Circle([0,0],1,edgecolor="r",facecolor="none")
        axes[0].add_patch(c)


        axes[0].set_title("unimodal upgrade")
        axes[0].set_ylabel("par2")
        axes[0].set_xlabel("par1")
        if func == "circle":
            fig.suptitle("iteration {0}, model: $obs1 = par1^2 + par2^2$, truth: $obs1 = 1$".format(noptmax))
        else:
            fig.suptitle("iteration {0}: model: Himmelblau's function: $obs1 = (par1^2 + par2 - 11)^2 + (par1 + par2^2 - 7)^2$, truth: $obs1 = 0$".format(noptmax))

            def h(x, y):
                return (x ** 2 + y - 11) ** 2 + (x + y ** 2 - 7) ** 2

            x = np.linspace(-4, 4, 100)
            y = x.copy()

            X, Y = np.meshgrid(x, y)

            H = np.flipud(h(X, Y))
            axes[0].imshow(H,alpha=0.7,extent=(-4,4,-4,4),cmap="jet")
            axes[1].imshow(H, alpha=0.7,extent=(-4,4,-4,4),cmap="jet")
            H = np.flipud(H)
            #axes[0].contour(X,Y,H,color="0.5",linestyles="dashed")

        plt.savefig(os.path.join(mm_d, "compare_{0:02d}.png".format(noptmax)),dpi=1000)
        plt.close(fig)

    pyemu.os_utils.run("ffmpeg -i compare_00.png -vf palettegen=16 -y palette.png",
                       cwd=mm_d)
    cmd = "ffmpeg -i compare_%02d.png -i palette.png -y -filter_complex "
    cmd += "\"fps=10,scale=720:-1:flags=lanczos[x];[x][1:v]paletteuse\" -y  -final_delay 150 compare.gif"
    pyemu.os_utils.run(cmd, cwd=mm_d)


def mm_invest():
    # model_d = "mm1"
    # test_d = os.path.join(model_d,"template")
    # m_d = os.path.join(model_d,"master_mm_circle")
    # pyemu.os_utils.start_workers(test_d, exe_path, "mm1.pst", worker_root=model_d, num_workers=25, master_dir=m_d)

    def h(x,y):
        return (x**2 + y - 11)**2 + (x + y**2 - 7)**2

    x = np.linspace(-4,4,100)
    y = x.copy()

    X,Y = np.meshgrid(x,y)

    H = h(X,Y)

    import matplotlib.pyplot as plt

    fig,ax = plt.subplots(1,1,figsize=(10,10))

    ax.imshow(H)
    plt.show()


def zdt1_ppw_worker(pst_name,host,port,helper=None):
    import pyemu
    ppw = pyemu.os_utils.PyPestWorker(pst_name,host,port,verbose=False)
    pvals = ppw.get_parameters()
    if pvals is None:
        return
    pvals.sort_index(inplace=True)
    #print(pvals)
    while True:

        sim = helper(pdf=pvals.values)
        sim = {n:s for n,s in zip(ppw._pst.obs_names,sim)}
        sim = [sim[o] for o in ppw.obs_names]
        #print(sim,np.array(sim).shape)
        ppw.send_observations(np.array(sim))
        pvals = ppw.get_parameters()
        if pvals is None:
            break
        pvals.sort_index(inplace=True)



def zdt1_weight_test():
    model_d = "zdt1"
    t_d = os.path.join(model_d,"zdt1_template")
    import sys
    sys.path.insert(0,t_d)
    
    from forward_run import helper 
    pst = pyemu.Pst(os.path.join(t_d,"zdt1.pst"))
    num_reals = 200
    np.random.seed(123331)
    pe = pyemu.ParameterEnsemble.from_uniform_draw(pst,num_reals=num_reals)
    oe = pyemu.ObservationEnsemble.from_gaussian_draw(pst,pyemu.Cov.from_observation_data(pst),num_reals=num_reals)
    for oname in pst.obs_names:
        oe.loc[:,oname] = 0.0
        pst.observation_data.loc[oname,"obgnme"] = oname # remove the ineq tag
        pst.observation_data.loc[oname, "obsval"] = 0.0  # remove the ineq tag

    weights = oe.copy()
    wseq = np.linspace(0.05,20,num_reals)
    weights.iloc[:,0] = wseq
    weights.iloc[:,1] = np.flipud(wseq)
    #print(weights)
    pst.pestpp_options = {}
    pst.pestpp_options["ies_par_en"] = "par.csv"
    pe.to_csv(os.path.join(t_d,"par.csv"))
    pst.pestpp_options["ies_obs_en"] = "obs.csv"
    oe.to_csv(os.path.join(t_d, "obs.csv"))
    
    df = pd.DataFrame(index=oe.index,columns=["obj_1","obj_2"])
    #wseq = np.linspace(0.01,0.99,num_reals)
    #df.loc[:,"obj_1"] = wseq
    #df.loc[:,"obj_2"] = np.flipud(wseq)
    
    first_third = np.arange(0,int(num_reals/3))
    middle_third = np.arange(int(first_third[-1]),int((num_reals*2)/3))
    last_third = np.arange(int(middle_third[-1]),num_reals)
    df.iloc[first_third,0] = 0.01
    df.iloc[first_third,1] = 0.99
    df.iloc[middle_third,0] = 0.5
    df.iloc[middle_third,1] = 0.5
    df.iloc[last_third,0] = 0.99
    df.iloc[last_third,1] = 0.01 
    df.to_csv(os.path.join(t_d,"phi_facs.csv"))
    

    pst.control_data.noptmax = 15


    #pst.write(os.path.join(t_d,"zdt1_ies.pst"))
    #m_d = os.path.join(model_d,"zdt1_master1_base")
    #pyemu.os_utils.start_workers(t_d,exe_path,"zdt1_ies.pst",num_workers=30,worker_root=model_d, verbose=True,master_dir=m_d,
    #    port=4200)

    pst.pestpp_options["ies_phi_factor_file"] = "phi_facs.csv"
    pst.pestpp_options["ies_phi_factors_by_real"] = True
    pst.pestpp_options["ies_weights_en"] = "weights.csv"
    pst.pestpp_options["ies_multimodal_alpha"] = .2
    pst.pestpp_options["ies_lambda_mults"] = 1.0
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["panther_agent_freeze_on_fail"] = True
    #pst.pestpp_options["ies_subset_size"] = -20
    weights.to_csv(os.path.join(t_d, "weights.csv"))
    pst.write(os.path.join(t_d, "zdt1_ies.pst"))
    m_d = os.path.join(model_d, "zdt1_master1")
    pyemu.os_utils.start_workers(t_d, exe_path, "zdt1_ies.pst", num_workers=10, worker_root=model_d, verbose=True,
                                 master_dir=m_d,port=4200,ppw_function=zdt1_ppw_worker,ppw_kwargs={"helper":helper})
    oe_file = os.path.join(m_d,"zdt1_ies.{0}.obs.csv".format(pst.control_data.noptmax))
    assert os.path.exists(oe_file)
    oe = pd.read_csv(oe_file,index_col=0)
    assert oe.loc[:,"obj_1"].min() < 0.2
    assert oe.loc[:,"obj_2"].min() < 1.0
    assert oe.loc[:,"obj_1"].max() > 0.5
    assert oe.loc[:,"obj_2"].max() > 4.0

def plot_zdt1_results(noptmax=None):
    import matplotlib.pyplot as plt
    m_d = os.path.join("zdt1","zdt1_master1")
    pst = pyemu.Pst(os.path.join(m_d,"zdt1_ies.pst"))
    if noptmax is None:
        noptmax = pst.control_data.noptmax
    
    for i in range(0,noptmax+1):
        oe_pt_base = pd.read_csv(os.path.join(m_d+"_base","zdt1_ies.{0}.obs.csv".format(i)),index_col=0)
        oe_pt = pd.read_csv(os.path.join(m_d,"zdt1_ies.{0}.obs.csv".format(i)),index_col=0)
        oe_pr = pd.read_csv(os.path.join(m_d,"zdt1_ies.0.obs.csv"),index_col=0)
        fig,ax = plt.subplots(1,1,figsize=(5,5))

        ax.scatter(oe_pr.iloc[:,0],oe_pr.iloc[:,1], marker=".",c="0.5", alpha=0.5,label="prior")
        ax.scatter(oe_pt.iloc[:, 0], oe_pt.iloc[:, 1], marker=".", c="b",label="fancy-sauce posterior")
        ax.scatter(oe_pt_base.iloc[:, 0], oe_pt_base.iloc[:, 1], marker=".", c="m", label="standard IES posterior")
        ax.legend(loc="upper right",fontsize=10)
        #ax.set_title("bi-objective zdt1 optimization benchmark",loc="left",fontsize=10)
        ax.tick_params(axis='both', which='major', labelsize=10)
        ax.set_xlabel("obs group 1",fontsize=12)
        ax.set_ylabel("obs group  2",fontsize=12)
        ax.set_xlim(0,1)
        ax.set_ylim(0,7)
        ax.set_title("iteration: {0}".format(i),loc="left",fontsize=12)
        ax.grid()
        plt.tight_layout()
        plt.savefig(os.path.join(m_d,"compare_{0:03d}.png".format(i)))
        plt.close(fig)
    fps = 3

    pyemu.os_utils.run("ffmpeg -y -i compare_{0:03d}.png -vf palettegen=256 palette.png".format(i),cwd=m_d)
    pyemu.os_utils.run("ffmpeg -r {0} -y -s 1920X1080 -i compare_%03d.png -i palette.png -filter_complex \"scale=720:-1:flags=lanczos[x];[x][1:v]paletteuse\" -final_delay 150 logo.gif".format(fps),
            cwd=m_d)


def tenpar_upgrade_on_disk_test_with_fixed():
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_upgrade_1_fixed")
    template_d = os.path.join(model_d, "test_template")

    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)

    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)

    pst.parameter_data.loc[:,"partrans"] = "log"
    cov = pyemu.Cov.from_parameter_data(pst)
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst,cov=cov,num_reals=10)
    pe.to_csv(os.path.join(template_d,"prior_all.csv"))
    pst.parameter_data.loc[pst.par_names[:-1],"partrans"] = "fixed"

    pst.pestpp_options = {}
    pst.pestpp_options.pop("ies_bad_phi",None)
    pst.pestpp_options.pop("ies_bad_phi_sigma", None)

    pst.pestpp_options["ies_no_noise"] = True
    pst.pestpp_options["ies_lambda_mults"] = [0.1, 1.0]
    pst.pestpp_options["lambda_scale_fac"] = [0.7, 1.0]
    pst.pestpp_options["ies_save_lambda_en"] = True
    #pst.pestpp_options["ies_subset_size"] = 1000
    pst.pestpp_options["ies_par_en"] = "prior_all.csv"
    pst.control_data.noptmax = 2
    pst_name = "pest_upgrade.pst"
    pst.write(os.path.join(template_d, pst_name))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pyemu.os_utils.run("{0} {1}".format(exe_path, pst_name), cwd=test_d)
    
    #pyemu.os_utils.start_workers(template_d, exe_path, pst_name, num_workers=8,
    #                             master_dir=test_d, worker_root=model_d, port=port)
    phi1 = pd.read_csv(os.path.join(test_d, pst_name.replace(".pst", ".phi.actual.csv")))
    pe1 = pd.read_csv(os.path.join(test_d, pst_name.replace(".pst", ".{0}.par.csv". \
                                                            format(pst.control_data.noptmax))), index_col=0)
    oe1 = pd.read_csv(os.path.join(test_d, pst_name.replace(".pst", ".{0}.obs.csv". \
                                                            format(pst.control_data.noptmax))), index_col=0)

    pst.pestpp_options["ies_upgrades_in_memory"] = False
    pst_name = "pest_upgrade_2.pst"
    pst.write(os.path.join(template_d, pst_name))
    test_d = os.path.join(model_d, "master_upgrade_2_fixed")
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pyemu.os_utils.run("{0} {1}".format(exe_path, pst_name), cwd=test_d)
    #pyemu.os_utils.start_workers(template_d, exe_path, pst_name, num_workers=8,
    #                             master_dir=test_d, worker_root=model_d, port=port)

    phi2 = pd.read_csv(os.path.join(test_d, pst_name.replace(".pst", ".phi.actual.csv")))
    pe2 = pd.read_csv(os.path.join(test_d, pst_name.replace(".pst", ".{0}.par.csv". \
                                                            format(pst.control_data.noptmax))), index_col=0)
    oe2 = pd.read_csv(os.path.join(test_d, pst_name.replace(".pst", ".{0}.obs.csv". \
                                                            format(pst.control_data.noptmax))), index_col=0)

    d = (phi1 - phi2).apply(lambda x: np.abs(x))
    print(d.max())
    print(d.max().max())
    assert d.max().max() < 1.0e-6

    d = (pe1 - pe2).apply(lambda x: np.abs(x))
    print(d.max())
    print(d.max().max())
    assert d.max().max() < 1.0e-6

    d = (oe1 - oe2).apply(lambda x: np.abs(x))
    print(d.max())
    print(d.max().max())
    assert d.max().max() < 1.0e-6


def tenpar_upgrade_on_disk_test_with_fixed2():
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_upgrade_1_fixed")
    template_d = os.path.join(model_d, "test_template")

    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)

    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    #shutil.copytree(template_d, test_d)

    #add a mess of fake pars
    fake_tpl = os.path.join(template_d,"fake.dat.tpl")
    with open(fake_tpl,'w') as f:
        f.write("ptf ~\n")
        for i in range(100):
            f.write("fake_{0} ~  fake_{0}   ~\n".format(i))
    pst.add_parameters(fake_tpl,pst_path=".")

    pst.parameter_data.loc[:,"partrans"] = "log"
    cov = pyemu.Cov.from_parameter_data(pst)
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst,cov=cov,num_reals=10)
    pe.to_csv(os.path.join(template_d,"prior_all.csv"))
    pst.parameter_data.loc[pst.par_names[1:],"partrans"] = "fixed"



    pst.pestpp_options["ies_no_noise"] = True
    pst.pestpp_options["ies_lambda_mults"] = [0.1, 1.0]
    pst.pestpp_options["lambda_scale_fac"] = [0.7, 1.0]
    pst.pestpp_options["ies_save_lambda_en"] = True
    pst.pestpp_options["ies_par_en"] = "prior_all.csv"
    pst.control_data.noptmax = 2
    pst_name = "pest_upgrade.pst"
    pst.write(os.path.join(template_d, pst_name))
    # pyemu.os_utils.run("{0} {1}".format(exe_path, pst_name), cwd=test_d)
    #pyemu.os_utils.start_workers(template_d, exe_path, pst_name, num_workers=8,
    #                             master_dir=test_d, worker_root=model_d, port=port)
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pyemu.os_utils.run("{0} {1}".format(exe_path, pst_name), cwd=test_d)
    phi1 = pd.read_csv(os.path.join(test_d, pst_name.replace(".pst", ".phi.actual.csv")))
    pe1 = pd.read_csv(os.path.join(test_d, pst_name.replace(".pst", ".{0}.par.csv". \
                                                            format(pst.control_data.noptmax))), index_col=0)
    oe1 = pd.read_csv(os.path.join(test_d, pst_name.replace(".pst", ".{0}.obs.csv". \
                                                            format(pst.control_data.noptmax))), index_col=0)

    pst.pestpp_options["ies_upgrades_in_memory"] = False
    pst.write(os.path.join(template_d, pst_name))
    test_d = os.path.join(model_d, "master_upgrade_2_fixed")
    #pyemu.os_utils.start_workers(template_d, exe_path, pst_name, num_workers=8,
    #                             master_dir=test_d, worker_root=model_d, port=port)
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pyemu.os_utils.run("{0} {1}".format(exe_path, pst_name), cwd=test_d)

    phi2 = pd.read_csv(os.path.join(test_d, pst_name.replace(".pst", ".phi.actual.csv")))
    pe2 = pd.read_csv(os.path.join(test_d, pst_name.replace(".pst", ".{0}.par.csv". \
                                                            format(pst.control_data.noptmax))), index_col=0)
    oe2 = pd.read_csv(os.path.join(test_d, pst_name.replace(".pst", ".{0}.obs.csv". \
                                                            format(pst.control_data.noptmax))), index_col=0)

    d = (phi1 - phi2).apply(lambda x: np.abs(x))
    print(d.max())
    print(d.max().max())
    assert d.max().max() < 1.0e-6

    d = (pe1 - pe2).apply(lambda x: np.abs(x))
    print(d.max())
    print(d.max().max())
    assert d.max().max() < 1.0e-6

    d = (oe1 - oe2).apply(lambda x: np.abs(x))
    print(d.max())
    print(d.max().max())
    assert d.max().max() < 1.0e-6


def tenpar_upgrade_on_disk_test_weight_ensemble_test():
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_upgrade_weight_1")
    template_d = os.path.join(model_d, "test_template")

    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)

    if os.path.exists(test_d):
        shutil.rmtree(test_d)

    pst.pestpp_options["ies_no_noise"] = True
    pst.pestpp_options["ies_lambda_mults"] = [0.1, 1.0]
    pst.pestpp_options["lambda_scale_fac"] = [0.7, 1.0]
    pst.pestpp_options["ies_save_lambda_en"] = True
    pst.pestpp_options["ies_upgrades_in_memory"] = False
    pst.pestpp_options["ies_debug_fail_remainder"] = True
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options["save_binary"] = True

    pst.control_data.noptmax = -1
    pst_name = "pest_weight.pst"
    pst.write(os.path.join(template_d, pst_name))
    #pyemu.os_utils.start_workers(template_d, exe_path, pst_name, num_workers=8,
    #                             master_dir=test_d, worker_root=model_d, port=port)
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pyemu.os_utils.run("{0} {1}".format(exe_path, pst_name), cwd=test_d)

    phi1 = pd.read_csv(os.path.join(test_d, pst_name.replace(".pst", ".phi.actual.csv")))

    oe1 = pyemu.ObservationEnsemble.from_binary(pst=pst,filename=os.path.join(test_d, pst_name.replace(".pst", ".0.obs.jcb")))
    pv = oe1.phi_vector
    pst = pyemu.Pst(os.path.join(test_d,pst_name))
    weights = {}
    for real,phi in zip(oe1.index.values,pv):
        pst = pyemu.Pst(os.path.join(test_d,pst_name))
        pst.res.loc[pst.obs_names,"modelled"] = oe1._df.loc[real,pst.obs_names].values
        print(real,phi,pst.phi)
        d = {n:0.5 for n in pst.nnz_obs_names}
        pst.adjust_weights(obs_dict=d)
        print("...",pst.phi)
        weights[real] = pst.observation_data.loc[pst.nnz_obs_names,"weight"].values

    wdf = pd.DataFrame(weights,index=pst.nnz_obs_names).T

    #print(wdf)

    wdf = pyemu.ObservationEnsemble(df=wdf,pst=pst)
    wdf.to_binary(os.path.join(template_d,"weights.jcb"))
    pst.pestpp_options["ies_weights_en"] = "weights.jcb"
    pst.pestpp_options["ies_par_en"] = "par1.jcb"
    pst.pestpp_options["ies_obs_en"] = "noise1.jcb"
    pst.pestpp_options["ies_restart_obs_en"] = "obs1.jcb"
    pst.pestpp_options.pop("ies_no_noise",None)
    shutil.copy2(os.path.join(test_d,"pest_weight.0.par.jcb"),os.path.join(template_d,"par1.jcb"))
    shutil.copy2(os.path.join(test_d,"pest_weight.0.obs.jcb"),os.path.join(template_d,"obs1.jcb"))
    #shutil.copy2(os.path.join(test_d,"pest_weight.obs+noise.jcb"),os.path.join(template_d,"noise.jcb"))
    noise = pyemu.ObservationEnsemble.from_binary(pst=pst,filename=os.path.join(test_d,"pest_weight.obs+noise.jcb"))
    noise = noise.loc[oe1.index,:]
    noise.to_binary(os.path.join(template_d,"noise1.jcb"))
    pst.control_data.noptmax = 3
    pst_name = "pest_weight_restart.pst"
    pst.write(os.path.join(template_d,pst_name))
    test_d = os.path.join(model_d, "master_upgrade_weight_2")
    #pyemu.os_utils.start_workers(template_d, exe_path, pst_name, num_workers=8,
    #                             master_dir=test_d, worker_root=model_d, port=port)
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pyemu.os_utils.run("{0} {1}".format(exe_path, pst_name), cwd=test_d)

    oe = pyemu.ObservationEnsemble.from_binary(pst=pst,filename=os.path.join(test_d,"pest_weight_restart.0.obs.jcb"))
    print(oe.shape)
    assert oe.shape[1] == pst.nobs
    pe = pyemu.ParameterEnsemble.from_binary(pst=pst,filename=os.path.join(test_d,"pest_weight_restart.0.par.jcb"))
    pe = pe.iloc[:-1,:]
    pe.to_binary(os.path.join(template_d,"par1.jcb"))
    noise = noise.iloc[:-1,:]
    noise.to_binary(os.path.join(template_d,"noise1.jcb"))
    wdf = wdf.iloc[:-1,:]
    wdf.to_binary(os.path.join(template_d,"weights.jcb"))
    pst.pestpp_options.pop("ies_restart_obs_en")
    pst.write(os.path.join(template_d, pst_name))
    test_d = os.path.join(model_d, "master_upgrade_weight_3")
    #pyemu.os_utils.start_workers(template_d, exe_path, pst_name, num_workers=8,
    #                            master_dir=test_d, worker_root=model_d, port=port)
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pyemu.os_utils.run("{0} {1}".format(exe_path, pst_name), cwd=test_d)
    


def tenpar_extra_binary_vars_test():
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_extra_vars")
    template_d = os.path.join(model_d, "test_template")

    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)

    if os.path.exists(test_d):
        shutil.rmtree(test_d)

    pst.pestpp_options["ies_no_noise"] = False
    pst.pestpp_options["ies_lambda_mults"] = [0.1, 1.0]
    pst.pestpp_options["lambda_scale_fac"] = [0.7, 1.0]
    pst.pestpp_options["ies_save_lambda_en"] = True
    pst.pestpp_options["ies_upgrades_in_memory"] = False
    pst.pestpp_options["ies_debug_fail_remainder"] = True
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options["save_binary"] = True

    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst,cov=pyemu.Cov.from_parameter_data(pst))
    pe.loc[:,"pextra1"] = 1.0
    pe.loc[:,"pextra2"] = 2.0
    pe.to_binary(os.path.join(template_d, "extra_par.jcb"))
    pst.pestpp_options["ies_par_en"] = "extra_par.jcb"

    oe = pyemu.ObservationEnsemble.from_gaussian_draw(pst=pst, cov=pyemu.Cov.from_observation_data(pst))
    oe.loc[:, "oextra1"] = 1.0
    oe.loc[:, "oextra2"] = 2.0
    oe.to_binary(os.path.join(template_d, "extra_noise.jcb"))
    pst.pestpp_options["ies_obs_en"] = "extra_noise.jcb"

    pst.control_data.noptmax = 1
    pst_name = "pest_extra.pst"
    pst.write(os.path.join(template_d,pst_name))
    #pyemu.os_utils.start_workers(template_d, exe_path, pst_name, num_workers=8,
    #                             master_dir=test_d, worker_root=model_d, port=port)
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pyemu.os_utils.run("{0} {1}".format(exe_path, pst_name), cwd=test_d)


def tenpar_adjust_weights_test():
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_adjust_weights")
    template_d = os.path.join(model_d, "test_template")

    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)

    if os.path.exists(test_d):
        shutil.rmtree(test_d)


    obs = pst.observation_data
    obs.loc[pst.obs_names[:4],"obgnme"] = "og1a"
    obs.loc[pst.obs_names[:4],"weight"] = 10 + 3.*(np.random.random(4))
    obs.loc[pst.obs_names[4:8],"obgnme"] = "og1b"
    obs.loc[pst.obs_names[4:8],"weight"] = 5
    obs.loc[pst.obs_names[8:12],"obgnme"] = "og3blahblah"
    obs.loc[pst.obs_names[8:12],"weight"] = 1
    obs.loc[pst.obs_names[12:],"obgnme"] = "og4yadayada"
    obs.loc[pst.obs_names[12:],"weight"] = 1.0

    obs.loc[:,"standard_deviation"] = 0.1
    obs.loc[pst.obs_names[12:],"standard_deviation"] = 1e-11
    obs.loc[pst.obs_names[12:],"obsval"] = 1e-9
    with open(os.path.join(template_d,"phi.csv"),'w') as f:
        f.write("og1,0.333333\n")
        f.write("og3,0.333333\n")
        f.write("og4,0.333333\n")
        

    pst.pestpp_options["ies_no_noise"] = False
    pst.pestpp_options["ies_lambda_mults"] = [0.1, 1.0]
    pst.pestpp_options["lambda_scale_fac"] = [0.7, 1.0]
    pst.pestpp_options["ies_save_lambda_en"] = True
    pst.pestpp_options["ies_upgrades_in_memory"] = False
    pst.pestpp_options["ies_debug_fail_remainder"] = False
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options["save_binary"] = True
    
    
    pst.pestpp_options['ies_verbose_level'] = 4
    pst.pestpp_options["ies_bad_phi_sigma"] = -1.5

    
    pst.control_data.noptmax = 2
    pst.pestpp_options["ies_drop_conflicts"] = False
    pst_name = "pest_adj.pst"
    pst.write(os.path.join(template_d,pst_name),version=2)
    #pyemu.os_utils.start_workers(template_d, exe_path, pst_name, num_workers=8,
    #                             master_dir=test_d, worker_root=model_d, port=port)
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pyemu.os_utils.run("{0} {1}".format(exe_path, pst_name), cwd=test_d)

    noise_file = os.path.join(test_d,"pest_adj.obs+noise.jcb")
    noise = pyemu.ObservationEnsemble.from_binary(pst=pst,filename=noise_file)
    print(noise)
    
    wdf_file = os.path.join(test_d,"pest_adj.weights.jcb")
    wdf = pyemu.ObservationEnsemble.from_binary(pst=pst,filename=wdf_file)
    print(wdf)
    for oname,weight in zip(obs.obsnme,obs.weight):
        print(oname,weight)
        print(wdf.loc[:,oname])
        assert wdf.loc[:,oname].std() < 1.0e-6
        assert np.abs(wdf.loc[:,oname].mean() - weight) < 1.0e-6

    wdf_file = os.path.join(test_d,"pest_adj.adjusted.weights.jcb")
    wdf = pyemu.ObservationEnsemble.from_binary(pst=pst,filename=wdf_file)
    adf_file = os.path.join(test_d,"pest_adj.adjusted.obs_data.csv")
    rei_file = os.path.join(test_d,"pest_adj.0.base.rei")
    pst.set_res(rei_file)
    assert os.path.exists(adf_file)
    adf = pd.read_csv(adf_file,index_col=0)
    for oname,weight in zip(adf.index,adf.weight):
        print(oname,weight)
        print(wdf.loc[:,oname])
        print(pst.res.loc[oname,"weight"])
        assert wdf.loc[:,oname].std() < 1.0e-6
        assert np.abs(wdf.loc[:,oname].mean() - weight) < 1.0e-3 #precision issues...
        assert np.abs(pst.res.loc[oname,"weight"] - weight) < 1.0e-3 #precision issues...


    
    
    pst.pestpp_options["ies_phi_factor_file"] = "phi.csv"
    pst.control_data.noptmax = 2
    pst_name = "pest_adj.pst"
    pst.write(os.path.join(template_d,pst_name),version=2)
    #pyemu.os_utils.start_workers(template_d, exe_path, pst_name, num_workers=8,
    #                             master_dir=test_d, worker_root=model_d, port=port)
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pyemu.os_utils.run("{0} {1}".format(exe_path, pst_name), cwd=test_d)

    wdf_file = os.path.join(test_d,"pest_adj.weights.jcb")
    wdf = pyemu.ObservationEnsemble.from_binary(pst=pst,filename=wdf_file)
    print(wdf)
    adf_file = os.path.join(test_d,"pest_adj.adjusted.obs_data.csv")
    assert os.path.exists(adf_file)
    adf = pd.read_csv(adf_file,index_col=0)
    rei_file = os.path.join(test_d,"pest_adj.0.base.rei")
    pst.set_res(rei_file)

    # for oname,weight in zip(obs.index,obs.weight):
    #     print(oname,weight)
    #     print(wdf.loc[:,oname])
    #     assert wdf.loc[:,oname].std() < 1.0e-6
    #     assert np.abs(wdf.loc[:,oname].mean() - weight) < 1.0e-3 #precision issues...
    
    wdf_file = os.path.join(test_d,"pest_adj.adjusted.weights.jcb")
    wdf = pyemu.ObservationEnsemble.from_binary(pst=pst,filename=wdf_file)
    print(wdf)
    for oname,weight in zip(adf.index,adf.weight):
        print(oname,weight)
        print(wdf.loc[:,oname])
        assert wdf.loc[:,oname].std() < 1.0e-6
        assert np.abs(wdf.loc[:,oname].mean() - weight) < 1.0e-3 #precision issues...
        assert np.abs(pst.res.loc[oname,"weight"] - weight) < 1.0e-3 #precision issues...


    sumfile = os.path.join(test_d,"pest_adj.0.obsgroupadj.summary.csv")
    assert os.path.exists(sumfile),sumfile
    ogdf = pd.read_csv(sumfile)
    phidf = pd.read_csv(os.path.join(test_d,"pest_adj.phi.meas.csv"),index_col=0)
    print(phidf.loc[0,"mean"])
    print(ogdf.adjusted_phi.sum())
    assert np.abs(ogdf.adjusted_phi.sum() - phidf.loc[0,"mean"]) < 1e-3
    assert phidf.loc[0,"mean"] > phidf.loc[2,"mean"]

    with open(os.path.join(template_d,"phi.csv"),'w') as f:
        f.write("og1,-999\n")
        f.write("og3,-999\n")
        f.write("og4,-999\n")

    pst.control_data.noptmax = 2
    pst.pestpp_options["ies_drop_conflicts"] = False
    pst_name = "pest_adj.pst"
    pst.write(os.path.join(template_d,pst_name),version=2)
    #pyemu.os_utils.start_workers(template_d, exe_path, pst_name, num_workers=8,
    #                             master_dir=test_d, worker_root=model_d, port=port)
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pyemu.os_utils.run("{0} {1}".format(exe_path, pst_name), cwd=test_d)

    wdf_file = os.path.join(test_d,"pest_adj.weights.jcb")
    wdf = pyemu.ObservationEnsemble.from_binary(pst=pst,filename=wdf_file)
    print(wdf)
    adf_file = os.path.join(test_d,"pest_adj.adjusted.obs_data.csv")
    assert os.path.exists(adf_file)
    adf = pd.read_csv(adf_file,index_col=0)
    rei_file = os.path.join(test_d,"pest_adj.0.base.rei")
    pst.set_res(rei_file)

    for oname,weight in zip(obs.index,obs.weight):
        print(oname,weight)
        print(wdf.loc[:,oname])
        assert wdf.loc[:,oname].std() < 1.0e-6
        assert np.abs(wdf.loc[:,oname].mean() - weight) < 1.0e-3 #precision issues...
        assert np.abs(pst.res.loc[oname,"weight"] - weight) < 1.0e-3 #precision issues...
    
    wdf_file = os.path.join(test_d,"pest_adj.adjusted.weights.jcb")
    wdf = pyemu.ObservationEnsemble.from_binary(pst=pst,filename=wdf_file)
    print(wdf)
    for oname,weight in zip(adf.index,adf.weight):
        print(oname,weight)
        print(wdf.loc[:,oname])
        assert wdf.loc[:,oname].std() < 1.0e-6
        assert np.abs(wdf.loc[:,oname].mean() - weight) < 1.0e-3 #precision issues...

    pst.set_res(os.path.join(test_d,"pest_adj.0.base.rei"))
    print(pst.phi)
    phidf = pd.read_csv(os.path.join(test_d,"pest_adj.phi.meas.csv"),index_col=0)
    print(phidf.loc[phidf.index[0],"base"])
    d = np.abs(pst.phi - phidf.loc[phidf.index[0],"base"])
    print(d)
    assert d < 0.01

    

def tenpar_adjust_weights_test_by_real():
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_adjust_weights_by_real")
    template_d = os.path.join(model_d, "test_template")

    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)

    if os.path.exists(test_d):
        shutil.rmtree(test_d)


    obs = pst.observation_data
    obs.loc[pst.obs_names[:4],"obgnme"] = "og1a"
    obs.loc[pst.obs_names[:4],"weight"] = 10 + 3.*(np.random.random(4))
    obs.loc[pst.obs_names[4:8],"obgnme"] = "og1b"
    obs.loc[pst.obs_names[4:8],"weight"] = 5
    obs.loc[pst.obs_names[8:12],"obgnme"] = "og3"
    obs.loc[pst.obs_names[8:12],"weight"] = 1
    obs.loc[pst.obs_names[12:],"obgnme"] = "og4"
    obs.loc[pst.obs_names[12:],"weight"] = 0.00001
    obs.loc[:,"standard_deviation"] = 0.1
    #with open(os.path.join(template_d,"phi.csv"),'w') as f:
    #    f.write("og1,0.333333\n")
    #    f.write("og3,0.333333\n")
    #    f.write("og4,0.333333\n")
    gdict = {}
    for g in pst.nnz_obs_groups:
        gdict[g] = obs.loc[obs.obgnme==g,"obsnme"].to_list()

    df = pd.DataFrame(columns=["og1","og3","og4"],index=np.arange(12))
    df.loc[:,:] = 0.33333
    df.loc[np.arange(3),["og1"]] = 0.9
    df.loc[np.arange(3),["og3"]] = 0.05
    df.loc[np.arange(3),["og4"]] = 0.05

    df.loc[np.arange(3,6),["og1"]] = 0.05
    df.loc[np.arange(3,6),["og3"]] = 0.9
    df.loc[np.arange(3,6),["og4"]] = 0.05

    df.loc[np.arange(6,11),["og1"]] = 0.05
    df.loc[np.arange(6,11),["og3"]] = 0.05
    df.loc[np.arange(6,11),["og4"]] = 0.09

    df.to_csv(os.path.join(template_d,"phi.csv"))

    pst.pestpp_options["ies_no_noise"] = False
    pst.pestpp_options["ies_lambda_mults"] = [0.1, 1.0]
    pst.pestpp_options["lambda_scale_fac"] = [0.7, 1.0]
    pst.pestpp_options["ies_save_lambda_en"] = True
    pst.pestpp_options["ies_upgrades_in_memory"] = False
    pst.pestpp_options["ies_debug_fail_remainder"] = True
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options["save_binary"] = True
    pst.pestpp_options["ies_phi_factor_file"] = "phi.csv"
    pst.pestpp_options["ies_phi_factors_by_real"] = True
    pst.pestpp_options["ies_drop_conflicts"] = True
    pst.pestpp_options['ies_verbose_level'] = 4
    pst.pestpp_options["ies_bad_phi_sigma"] = -1.5

    
    pst.control_data.noptmax = 2
    pst_name = "pest_adj.pst"
    pst.write(os.path.join(template_d,pst_name),version=2)
    #pyemu.os_utils.start_workers(template_d, exe_path, pst_name, num_workers=8,
    #                             master_dir=test_d, worker_root=model_d, port=port)
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pyemu.os_utils.run("{0} {1}".format(exe_path, pst_name), cwd=test_d)

    sumfile = os.path.join(test_d,"pest_adj.0.obsgroupadj.summary.csv")
    assert os.path.exists(sumfile),sumfile
    ogdf = pd.read_csv(sumfile)
    phidf = pd.read_csv(os.path.join(test_d,"pest_adj.phi.actual.csv"),index_col=0)
    #print(phidf.loc[0,"mean"])
    #print(ogdf.adjusted_phi.sum())
    #assert np.abs(ogdf.adjusted_phi.sum() - phidf.loc[0,"mean"]) < 1e-3
    #assert phidf.loc[0,"mean"] > phidf.loc[2,"mean"]

    wdf_file = os.path.join(test_d,"pest_adj.weights.jcb")
    wdf = pyemu.ObservationEnsemble.from_binary(pst=pst,filename=wdf_file)
    #print(wdf)
    for oname,weight in zip(obs.index,obs.weight):
        print(oname,weight)
        print(wdf.loc[:,oname])
        assert wdf.loc[:,oname].std() < 1.0e-6
        assert np.abs(wdf.loc[:,oname].mean() - weight) < 1.0e-3 #precision issues...

    adf_file = os.path.join(test_d,"pest_adj.adjusted.obs_data.csv")
    assert os.path.exists(adf_file)
    adf = pd.read_csv(adf_file,index_col=0)
    wdf_file = os.path.join(test_d,"pest_adj.adjusted.weights.jcb")
    wdf = pyemu.ObservationEnsemble.from_binary(pst=pst,filename=wdf_file)
    print(wdf)
    for oname,weight in zip(adf.index,adf.weight):
        
        print(wdf.loc[:,oname])
        print(oname,weight)
        if weight == 0.0:

            assert wdf.loc[:,oname].std() < 1.0e-6
            assert wdf.loc[:,oname].max() == 0.0
            assert np.abs(wdf.loc[:,oname].mean() - weight) < 1.0e-3 #precision issues...

        else:

            assert wdf.loc[:,oname].std() > 0.1
            assert np.abs(wdf.loc[:,oname].mean() - weight) > 1.0e-3 #precision issues...

    shutil.copy2(os.path.join(test_d,"pest_adj.0.obs.jcb"),os.path.join(template_d,"restart.jcb"))
    shutil.copy2(os.path.join(test_d,"pest_adj.0.par.jcb"),os.path.join(template_d,"par.jcb"))
    shutil.copy2(os.path.join(test_d,"pest_adj.obs+noise.jcb"),os.path.join(template_d,"noise.jcb"))
    pst.pestpp_options["ies_par_en"] = "par.jcb"
    pst.pestpp_options["ies_obs_en"] = "noise.jcb"
    pst.pestpp_options["ies_restart_obs_en"] = "restart.jcb"

    pst.control_data.noptmax = 1
    pst_name = "pest_adj.pst"
    pst.write(os.path.join(template_d,pst_name),version=2)
    
    #pyemu.os_utils.start_workers(template_d, exe_path, pst_name, num_workers=8,
    #                             master_dir=test_d, worker_root=model_d, port=port)
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pyemu.os_utils.run("{0} {1}".format(exe_path, pst_name), cwd=test_d)

    sumfile = os.path.join(test_d,"pest_adj.0.obsgroupadj.summary.csv")
    assert os.path.exists(sumfile),sumfile
    ogdf = pd.read_csv(sumfile)
    adf_file = os.path.join(test_d,"pest_adj.adjusted.obs_data.csv")
    assert os.path.exists(adf_file)
    adf = pd.read_csv(adf_file,index_col=0)
    wdf_file = os.path.join(test_d,"pest_adj.adjusted.weights.jcb")
    wdf = pyemu.ObservationEnsemble.from_binary(pst=pst,filename=wdf_file)
    ndf = pyemu.ObservationEnsemble.from_binary(pst=pst,filename=os.path.join(test_d,"pest_adj.obs+noise.jcb"))
    odf = pyemu.ObservationEnsemble.from_binary(pst=pst,filename=os.path.join(test_d,"pest_adj.0.obs.jcb"))
    print(wdf)
    for real in odf.index:
        n = ndf.loc[real,:]
        o = odf.loc[real,:]
        w = wdf.loc[real,:]
        swr = (w * (n-o))**2
        #print(swr)
        for g,onames in gdict.items():
            phi = swr.loc[onames].sum()
            #print(real,g,phi)
            try:
                ophi = ogdf.loc[ogdf.apply(lambda x: str(x.realization).upper() ==str(real).upper() and x.group.upper()==g.upper(),axis=1),"adjusted_phi"].values[0]
            except Exception as e:
                print('error',real,g,e)
                continue
            d = np.abs((phi-ophi)/phi)
            print(real,g,phi,ophi,d)
            assert d < 1e-5
            

    pst.control_data.noptmax = 1
    pst.pestpp_options.pop("ies_obs_en")
    pst.pestpp_options.pop("ies_restart_obs_en")
    pst.pestpp_options.pop("ies_par_en")
    
    pst.pestpp_options["ies_num_reals"] = 20
    pst_name = "pest_adj.pst"
    pst.write(os.path.join(template_d,pst_name),version=2)
    try:
        pyemu.os_utils.start_workers(template_d, exe_path, pst_name, num_workers=8,
                                     master_dir=test_d, worker_root=model_d, port=port)
    except Exception as e:
        pass
    else:
        raise Exception("should have failed")


def tenpar_drop_violations_test():
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_drop_violations")
    template_d = os.path.join(model_d, "test_template")

    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)

    if os.path.exists(test_d):
        shutil.rmtree(test_d)


    obs = pst.observation_data
    obs.loc[:,"drop_violations"] = False

    obs.loc[pst.obs_names[0],"obgnme"] = "less_than"
    obs.loc[pst.obs_names[1],"obgnme"] = "less_than"
    obs.loc[pst.obs_names[0],"weight"] = 1.0
    obs.loc[pst.obs_names[1],"weight"] = 1.0
    obs.loc[pst.obs_names[1],"obgnme"] = "less_than"
    obs.loc[pst.obs_names[1],"drop_violations"] = True
    obs.loc[pst.obs_names[2],"weight"] = 1.0
    obs.loc[pst.obs_names[2],"obgnme"] = "less_than"
    obs.loc[pst.obs_names[2],"drop_violations"] = True
    obs.loc[pst.obs_names[2],"obsval"] = -10000.0


    obs.loc[pst.obs_names[-1],"obgnme"] = "greater_than"
    obs.loc[pst.obs_names[-1],"weight"] = 1.0
    obs.loc[pst.obs_names[-2],"weight"] = 1.0
    obs.loc[pst.obs_names[-2],"obgnme"] = "greater_than"
    obs.loc[pst.obs_names[-2],"drop_violations"] = True
    obs.loc[pst.obs_names[-3],"weight"] = 0.0
    obs.loc[pst.obs_names[-3],"obgnme"] = "greater_than"
    obs.loc[pst.obs_names[-3],"drop_violations"] = True
    obs.loc[pst.obs_names[-3],"obsval"] = 100000



    pst.pestpp_options["ies_no_noise"] = False
    pst.pestpp_options["ies_lambda_mults"] = [0.1, 1.0]
    pst.pestpp_options["lambda_scale_fac"] = [0.7, 1.0]
    pst.pestpp_options["ies_save_lambda_en"] = True
    pst.pestpp_options["ies_upgrades_in_memory"] = False
    pst.pestpp_options["ies_num_reals"] = 20
    pst.pestpp_options["ies_num_threads"] = 4
    pst.pestpp_options["save_binary"] = False
    
    
    pst.pestpp_options['ies_verbose_level'] = 4
    pst.pestpp_options["ies_bad_phi_sigma"] = -1.5
    pst.pestpp_options["ies_autoadaloc"] = True

    
    pst.control_data.noptmax = 4
    pst.pestpp_options["ies_drop_conflicts"] = True
    pst_name = "pest_viol.pst"
    pst.write(os.path.join(template_d,pst_name),version=2)
    try:
        pyemu.os_utils.start_workers(template_d, exe_path, pst_name, num_workers=8,
                                     master_dir=test_d, worker_root=model_d, port=port)
    except Exception as e:
        pass
    else:
        raise Exception("should have failed")

    obs.loc[pst.obs_names[2],"weight"] = 0.0
    pst_name = "pest_viol.pst"
    pst.write(os.path.join(template_d,pst_name),version=2)
    #pyemu.os_utils.start_workers(template_d, exe_path, pst_name, num_workers=8,
    #                             master_dir=test_d, worker_root=model_d, port=port)
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pyemu.os_utils.run("{0} {1}".format(exe_path, pst_name), cwd=test_d)

    for i in range(1,pst.control_data.noptmax+1):
        oe = pd.read_csv(os.path.join(test_d,"pest_viol.{0}.obs.csv".format(i)),index_col=0)
        oe = oe.loc[oe.index.map(lambda x: x != "base"),:]
        print(oe.loc[:,pst.obs_names[1]],obs.loc[pst.obs_names[1],"obsval"])
        assert np.all(oe.loc[:,pst.obs_names[1]]<obs.loc[pst.obs_names[1],"obsval"])
        print(oe.loc[:,pst.obs_names[-2]],obs.loc[pst.obs_names[-2],"obsval"])
        assert np.all(oe.loc[:,pst.obs_names[-2]]>obs.loc[pst.obs_names[-2],"obsval"])

        print(oe.loc[:,pst.obs_names[2]],obs.loc[pst.obs_names[2],"obsval"])
        assert np.any(oe.loc[:,pst.obs_names[2]]>=obs.loc[pst.obs_names[2],"obsval"])
        print(oe.loc[:,pst.obs_names[-3]],obs.loc[pst.obs_names[-3],"obsval"])
        assert np.any(oe.loc[:,pst.obs_names[-3]]<=obs.loc[pst.obs_names[-3],"obsval"])


def tenpar_mean_iter_test():
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_mean_iter")
    template_d = os.path.join(model_d, "test_template")

    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d,test_d)
    pst_name = "pest.pst"
    pst = pyemu.Pst(os.path.join(template_d,pst_name))
    pst.pestpp_options["ies_n_iter_mean"] = 3
    pst.pestpp_options["ies_initial_lambda"] = -100
    pst.control_data.noptmax = 21 # hard coded to test results below
    pst.pestpp_options["ies_num_reals"] = 50
    pst.pestpp_options["ies_debug_fail_remainder"] = True
    pst.pestpp_options["ies_debug_fail_subset"] = True

    
    pst.write(os.path.join(test_d,pst_name))
    pyemu.os_utils.run("{0} {1}".format(exe_path,pst_name),cwd=test_d)
    #pyemu.os_utils.start_workers(template_d, exe_path, pst_name, num_workers=5,
    #                                 master_dir=test_d, worker_root=model_d, port=port)
    test_d1 = test_d + "_base"
    if os.path.exists(test_d1):
        shutil.rmtree(test_d1)
    shutil.copytree(template_d,test_d1)
    
    phidf = pd.read_csv(os.path.join(test_d,"pest.phi.actual.csv"),index_col=0)
    print(phidf.loc[:,"mean"])
    print(phidf.shape)
    assert phidf.shape[0] == 29 #hard coded to noptmax above
    assert phidf.shape[1] == 55 #50 reals + summary stats

    for i in phidf.index.values:
        oe = pd.read_csv(os.path.join(test_d,"pest.{0}.obs.csv".format(i)),index_col=0)
        assert oe.shape[1] == pst.nobs

def tenpar_mean_iter_test_sched():

    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_mean_iter_sched")
    template_d = os.path.join(model_d, "test_template")

    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d,test_d)
    pst_name = "pest.pst"
    pst = pyemu.Pst(os.path.join(template_d,pst_name))
    pst.pestpp_options["ies_n_iter_mean"] = [-1,-3,5,999]
    #pst.pestpp_options["ies_reinflate_factor"] = [1.0,0.9,0.8,0.7]
    
    pst.pestpp_options["ies_initial_lambda"] = -100
    pst.control_data.noptmax = 21 # hard coded to test results below
    pst.pestpp_options["ies_num_reals"] = 50
    pst.pestpp_options["ies_debug_fail_remainder"] = True
    pst.pestpp_options["ies_debug_fail_subset"] = True
    pst.pestpp_options["ies_no_noise"] = False
    

    
    pst.write(os.path.join(test_d,pst_name))
    pyemu.os_utils.run("{0} {1}".format(exe_path,pst_name),cwd=test_d)
    #pyemu.os_utils.start_workers(template_d, exe_path, pst_name, num_workers=5,
    #                                 master_dir=test_d, worker_root=model_d, port=port)
    test_d1 = test_d + "_base"
    if os.path.exists(test_d1):
        shutil.rmtree(test_d1)
    shutil.copytree(template_d,test_d1)
    
    phidf = pd.read_csv(os.path.join(test_d,"pest.phi.actual.csv"),index_col=0)
    print(phidf.loc[:,"mean"])
    print(phidf.shape)
    assert phidf.shape[0] == 25 #hard coded to noptmax above
    assert phidf.shape[1] == 55 #50 reals + summary stats

    for i in phidf.index.values:
        oe = pd.read_csv(os.path.join(test_d,"pest.{0}.obs.csv".format(i)),index_col=0)
        assert oe.shape[1] == pst.nobs

    count = 0
    with open(os.path.join(test_d,"pest.rec"),'r') as f:
        for line in f:
            if "running new mean-shifted prior realizations" in line:
                count += 1
    print(count)
    assert count == 3
    


    test_d = os.path.join(model_d, "master_mean_iter_sched_facs")
    template_d = os.path.join(model_d, "test_template")

    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d,test_d)
    pst_name = "pest.pst"
    pst = pyemu.Pst(os.path.join(template_d,pst_name))
    pst.pestpp_options["ies_n_iter_mean"] = [-1,-3,5,999]
    pst.pestpp_options["ies_reinflate_factor"] = [1.0,0.9,0.8,0.7]
    
    pst.pestpp_options["ies_initial_lambda"] = -100
    pst.control_data.noptmax = 21 # hard coded to test results below
    pst.pestpp_options["ies_num_reals"] = 50
    pst.pestpp_options["ies_debug_fail_remainder"] = True
    pst.pestpp_options["ies_debug_fail_subset"] = True
    pst.pestpp_options["ies_no_noise"] = False

    
    pst.write(os.path.join(test_d,pst_name))
    pyemu.os_utils.run("{0} {1}".format(exe_path,pst_name),cwd=test_d)
    #pyemu.os_utils.start_workers(template_d, exe_path, pst_name, num_workers=5,
    #                                 master_dir=test_d, worker_root=model_d, port=port)
    test_d1 = test_d + "_base"
    if os.path.exists(test_d1):
        shutil.rmtree(test_d1)
    shutil.copytree(template_d,test_d1)
    
    phidf = pd.read_csv(os.path.join(test_d,"pest.phi.actual.csv"),index_col=0)
    print(phidf.loc[:,"mean"])
    print(phidf.shape)
    assert phidf.shape[0] == 25 #hard coded to noptmax above
    assert phidf.shape[1] == 55 #50 reals + summary stats

    for i in phidf.index.values:
        oe = pd.read_csv(os.path.join(test_d,"pest.{0}.obs.csv".format(i)),index_col=0)
        assert oe.shape[1] == pst.nobs

    count = 0
    facs = []
    with open(os.path.join(test_d,"pest.rec"),'r') as f:
        for line in f:
            if "running new mean-shifted prior realizations" in line:
                count += 1
            elif "reinflation factor:" in line:
                facs.append(float(line.split("reinflation factor:")[1]))

    print(count)
    assert count == 3
    print(facs)
    assert len(facs) == 3
    ffacs = pst.pestpp_options["ies_reinflate_factor"][:-1]
    diff = np.array(facs) - ffacs
    print(diff)
    assert diff.sum() == 0.0

    pst.pestpp_options["save_binary"] = True
    pst.pestpp_options["save_dense"]= True
    test_d = os.path.join(model_d, "master_mean_iter_sched_facs_dense")
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d,test_d)
    
    pst.write(os.path.join(test_d,pst_name))
    pyemu.os_utils.run("{0} {1}".format(exe_path,pst_name),cwd=test_d)

    phidf = pd.read_csv(os.path.join(test_d,"pest.phi.actual.csv"),index_col=0)
    print(phidf.loc[:,"mean"])
    print(phidf.shape)
    assert phidf.shape[0] == 25 #hard coded to noptmax above
    assert phidf.shape[1] == 55 #50 reals + summary stats

    pst.pestpp_options["ies_par_en"] = "pest.24.par.bin"
    pst.pestpp_options["ies_restart_obs_en"] = "pest.24.obs.bin"
    pst.pestpp_options["ies_obs_en"] = "pest.obs+noise.bin"
    pst.pestpp_options["ies_weight_en"] = "pest.weights.bin"

    pst.write(os.path.join(test_d,"test.pst"))
    pyemu.os_utils.run("{0} {1}".format(exe_path,"test.pst"),cwd=test_d)

    phidf2 = pd.read_csv(os.path.join(test_d,"test.phi.actual.csv"),index_col=0)
    print(phidf2.loc[:,"mean"])
    print(phidf2.shape)
    assert phidf2.shape[0] == 25 #hard coded to noptmax above
    assert phidf2.shape[1] == 45 #restart with 40 reals + summary stats
    
    pst.pestpp_options.pop("save_dense")

    pst.write(os.path.join(test_d,"test2.pst"))
    pyemu.os_utils.run("{0} {1}".format(exe_path,"test2.pst"),cwd=test_d)

    phidf2 = pd.read_csv(os.path.join(test_d,"test2.phi.actual.csv"),index_col=0)
    print(phidf2.loc[:,"mean"])
    print(phidf2.shape)
    assert phidf2.shape[0] == 25 #hard coded to noptmax above
    assert phidf2.shape[1] == 45 #restart with 40 reals + summary stats
    


def tenpar_mean_iter_sched_phifac_test():

    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_mean_iter_sched_phifac_single")
    template_d = os.path.join(model_d, "test_template")
    
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d,test_d)
    pst_name = "pest.pst"
    pst = pyemu.Pst(os.path.join(template_d,pst_name))
    pst.pestpp_options["ies_n_iter_mean"] = [-1,-3,5,999]
    #pst.pestpp_options["ies_reinflate_factor"] = [1.0,0.9,0.8,0.7]
    
    pst.pestpp_options["ies_initial_lambda"] = -100
    pst.control_data.noptmax = 21 # hard coded to test results below
    pst.pestpp_options["ies_num_reals"] = 50
    pst.pestpp_options["ies_debug_fail_remainder"] = True
    pst.pestpp_options["ies_debug_fail_subset"] = True
    pst.pestpp_options["ies_no_noise"] = False
    onames = pst.obs_names
    obs = pst.observation_data
    obs.loc[:,"weight"] = 1.0
    
    obs.loc[onames[:3],"obgnme"] = "less_than"
    obs.loc[onames[:3],"obsval"] = obs.loc[onames[:3],"obsval"] + 0.1
    obs.loc[onames[3:5],"obgnme"] = "first_group"
    obs.loc[onames[5:10],"obgnme"] = "second_group"
    obs.loc[onames[10:],"obgnme"] = "greater_than"
    obs.loc[onames[10:],"obsval"] = obs.loc[onames[10:],"obsval"] 
    
    facs = {"tag":["less_than","first_group","second_group","greater_than"], "prop":[0.2,0.2,0.1,0.5]}

    df = pd.DataFrame(data=facs)
    df.to_csv(os.path.join(test_d,"phi.csv"),index=False,header=False)
    pst.pestpp_options["ies_phi_factor_file"] = "phi.csv"
    #pst.pestpp_options["debug_parse_only"] = True
    pst.pestpp_options["ies_verbose_level"] = 2
    pst.control_data.noptmax = 20
    
    pst.write(os.path.join(test_d,pst_name),version=2)
    pyemu.os_utils.run("{0} pest.pst".format(exe_path),cwd=test_d)



    test_d = os.path.join(model_d, "master_mean_iter_sched_phifac_multiple")
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d,test_d)

    facs = {"less_than":0.2,"first_group":0.3,"second_group":0.3,"greater_than":0.2}
    df = pd.DataFrame(data=facs,index=np.arange(pst.pestpp_options["ies_num_reals"]))
    df.to_csv(os.path.join(test_d,"phi.csv"))
    pst.pestpp_options["ies_phi_factor_file"] = "phi.csv"
    pst.pestpp_options["ies_phi_factors_by_real"] = True
    pst.pestpp_options["debug_parse_only"] = True
    pst.control_data.noptmax = -2
    
    pst.write(os.path.join(test_d,pst_name),version=2)
    pyemu.os_utils.run("{0} pest.pst".format(exe_path),cwd=test_d)
    
    pst.control_data.noptmax = -1
    pst.write(os.path.join(test_d,pst_name),version=2)
    pyemu.os_utils.run("{0} pest.pst".format(exe_path),cwd=test_d)

    pst.control_data.noptmax = 0
    pst.write(os.path.join(test_d,pst_name),version=2)
    pyemu.os_utils.run("{0} pest.pst".format(exe_path),cwd=test_d)
    
    pst.control_data.noptmax = 21
    pst.write(os.path.join(test_d,pst_name),version=2)
    pyemu.os_utils.run("{0} pest.pst".format(exe_path),cwd=test_d)

    pst.pestpp_options.pop("debug_parse_only")
    
    pst.write(os.path.join(test_d,pst_name))
    pyemu.os_utils.run("{0} {1}".format(exe_path,pst_name),cwd=test_d)
    #pyemu.os_utils.start_workers(template_d, exe_path, pst_name, num_workers=5,
    #                                 master_dir=test_d, worker_root=model_d, port=port)
    test_d1 = test_d + "_base"
    if os.path.exists(test_d1):
        shutil.rmtree(test_d1)
    shutil.copytree(template_d,test_d1)
    
    phidf = pd.read_csv(os.path.join(test_d,"pest.phi.actual.csv"),index_col=0)
    print(phidf.loc[:,"mean"])
    print(phidf.shape)
    assert phidf.shape[0] == 25 #hard coded to noptmax above
    assert phidf.shape[1] == 55 #50 reals + summary stats

    for i in phidf.index.values:
        oe = pd.read_csv(os.path.join(test_d,"pest.{0}.obs.csv".format(i)),index_col=0)
        assert oe.shape[1] == pst.nobs

    count = 0
    adj_count = 0
    with open(os.path.join(test_d,"pest.rec"),'r') as f:
        for line in f:
            if "running new mean-shifted prior realizations" in line:
                count += 1
            elif "adjusting weights using phi factors in file" in line:
                adj_count += 1
    print(count)
    assert count == 3
    assert adj_count == 1
    

    test_d = os.path.join(model_d, "master_mean_iter_sched_facs")
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d,test_d)
    pst_name = "pest.pst"
    pst = pyemu.Pst(os.path.join(template_d,pst_name))
    pst.pestpp_options["ies_n_iter_mean"] = [-1,-3,5,999]
    pst.pestpp_options["ies_reinflate_factor"] = [1.0,0.9,0.8,0.7]
    
    pst.pestpp_options["ies_initial_lambda"] = -100
    pst.control_data.noptmax = 21 # hard coded to test results below
    pst.pestpp_options["ies_num_reals"] = 50
    pst.pestpp_options["ies_debug_fail_remainder"] = True
    pst.pestpp_options["ies_debug_fail_subset"] = True
    pst.pestpp_options["ies_debug_high_upgrade_phi"] = True
    pst.pestpp_options["ies_no_noise"] = False

    onames = pst.obs_names
    obs = pst.observation_data
    obs.loc[:,"weight"] = 1.0
    
    obs.loc[onames[:3],"obgnme"] = "less_than"
    obs.loc[onames[:3],"obsval"] = obs.loc[onames[:3],"obsval"] + 0.1
    obs.loc[onames[3:5],"obgnme"] = "first_group"
    obs.loc[onames[5:10],"obgnme"] = "second_group"
    obs.loc[onames[10:],"obgnme"] = "greater_than"
    obs.loc[onames[10:],"obsval"] = obs.loc[onames[10:],"obsval"] - 1.5
    
    facs = {"less_than":0.2,"first_group":0.3,"second_group":0.3,"greater_than":0.2}
    df = pd.DataFrame(data=facs,index=np.arange(pst.pestpp_options["ies_num_reals"]))
    df.to_csv(os.path.join(test_d,"phi.csv"))
    pst.pestpp_options["ies_phi_factor_file"] = "phi.csv"
    pst.pestpp_options["ies_phi_factors_by_real"] = True
    
    pst.write(os.path.join(test_d,pst_name))
    pyemu.os_utils.run("{0} {1}".format(exe_path,pst_name),cwd=test_d)
    #pyemu.os_utils.start_workers(template_d, exe_path, pst_name, num_workers=5,
    #                                 master_dir=test_d, worker_root=model_d, port=port)
    test_d1 = test_d + "_base"
    if os.path.exists(test_d1):
        shutil.rmtree(test_d1)
    shutil.copytree(template_d,test_d1)
    
    phidf = pd.read_csv(os.path.join(test_d,"pest.phi.actual.csv"),index_col=0)
    print(phidf.loc[:,"mean"])
    print(phidf.shape)
    assert phidf.shape[0] == 25 #hard coded to noptmax above
    assert phidf.shape[1] == 55 #50 reals + summary stats

    for i in phidf.index.values:
        oe = pd.read_csv(os.path.join(test_d,"pest.{0}.obs.csv".format(i)),index_col=0)
        assert oe.shape[1] == pst.nobs

    count = 0
    facs = []
    adj_count = 0
    with open(os.path.join(test_d,"pest.rec"),'r') as f:
        for line in f:
            if "running new mean-shifted prior realizations" in line:
                count += 1
            elif "reinflation factor:" in line:
                facs.append(float(line.split("reinflation factor:")[1]))
            elif "adjusting weights using phi factors in file" in line:
                adj_count += 1

    print(count)
    assert count == 3
    print(facs)
    assert len(facs) == 3
    assert adj_count == 1
    ffacs = pst.pestpp_options["ies_reinflate_factor"][:-1]
    diff = np.array(facs) - ffacs
    print(diff)
    assert diff.sum() == 0.0


def twopar_freyberg_resp_surface_invest():
    model_d = "twopar_freyberg"
    resp_d = os.path.join(model_d, "master_resp_surface")
    template_d = os.path.join(model_d, "template")
    pst_path = os.path.join(template_d,"freyberg.pst")
    pst = pyemu.Pst(pst_path)
    # if os.path.exists(resp_d):
    #    shutil.rmtree(resp_d)
    # shutil.copytree(template_d,resp_d)
    # hk1_vals,rch0_vals = [],[]
    # par = pst.parameter_data
    # for hinc in np.linspace(par.loc["hk1","parlbnd"],par.loc["hk1","parubnd"],40):
    #     for rinc in np.linspace(par.loc["rch0","parlbnd"],par.loc["rch0","parubnd"],40):
    #         hk1_vals.append(hinc)
    #         rch0_vals.append(rinc)
    # df = pd.DataFrame({"hk1":hk1_vals,"rch0":rch0_vals})
    # for pname in pst.par_names:
    #     if pname in ["hk1","rch0"]:
    #         continue
    #     df[pname] = par.loc[pname,"parval1"]
    # df.to_csv(os.path.join(resp_d,"sweep_in.csv"))
    # pyemu.os_utils.run("{0} freyberg.pst".format(exe_path.replace("-ies","-swp")),cwd=resp_d)
    
    #if not os.path.exists(resp_d):
    #    pyemu.os_utils.start_workers(template_d,exe_path.replace("-ies","-swp"),"freyberg.pst",master_dir=resp_d,num_workers=20,worker_root=model_d)
    
    pst.control_data.noptmax = 6
    pst.pestpp_options = {"par_sigma_range":2}
    pst.write(pst_path)
    m_d_base = os.path.join(model_d,"master_base")
    if os.path.exists(m_d_base):
        shutil.rmtree(m_d_base)
    shutil.copytree(template_d,m_d_base)
    pyemu.os_utils.run("{0} freyberg.pst".format(exe_path),cwd=m_d_base)
    #yemu.os_utils.start_workers(template_d,exe_path,"freyberg.pst",master_dir=m_d_base,num_workers=20,worker_root=model_d)
    pst.pestpp_options = {"ies_n_iter_mean":pst.control_data.noptmax-1,"par_sigma_range":2}
    pst.write(pst_path)
    m_d_allmean = os.path.join(model_d,"master_allmean")
    if os.path.exists(m_d_allmean):
        shutil.rmtree(m_d_allmean)
    shutil.copytree(template_d,m_d_allmean)
    pyemu.os_utils.run("{0} freyberg.pst".format(exe_path),cwd=m_d_allmean)
    
    #pyemu.os_utils.start_workers(template_d,exe_path,"freyberg.pst",master_dir=m_d_allmean,num_workers=20,worker_root=model_d)
    pst.pestpp_options = {"ies_n_iter_mean":int(pst.control_data.noptmax/2),"par_sigma_range":2}
    pst.write(pst_path)
    m_d_somemean = os.path.join(model_d,"master_somemean")
    if os.path.exists(m_d_somemean):
        shutil.rmtree(m_d_somemean)
    shutil.copytree(template_d,m_d_somemean)
    pyemu.os_utils.run("{0} freyberg.pst".format(exe_path),cwd=m_d_somemean)
    
    #pyemu.os_utils.start_workers(template_d,exe_path,"freyberg.pst",master_dir=m_d_somemean,num_workers=20,worker_root=model_d)
    
def plot_twopar_resp_results():
    model_d = "twopar_freyberg"
    resp_d = os.path.join(model_d, "master_resp_surface")
    template_d = os.path.join(model_d, "template")
    pst_path = os.path.join(template_d,"freyberg.pst")
    pst = pyemu.Pst(pst_path)
    #m_ds = [os.path.join(model_d,d) for d in os.listdir(model_d) if os.path.isdir(os.path.join(model_d,d)) and d.startswith("master") and "resp" not in d]
    m_ds = ["master_base","master_allmean","master_somemean"]
    labels = ["A) standard iters","B) mean shift iter 6","C) mean shift iter 4"]
    plabels = ["D) phi sequence","E) phi sequence","F) phi sequence"]
    
    m_ds = [os.path.join(model_d,m_d) for m_d in m_ds]
    fig,axes = plt.subplots(2,len(m_ds),figsize=(3*len(m_ds),5.5))
    for im_d,(m_d,ax,label,pax,plabel) in enumerate(zip(m_ds,axes[0,:],labels,axes[1,:],plabels)):
        ax,resp_surf = plot_response_surface(WORKING_DIR=resp_d,ax=ax)
        ax.set_title(label,loc="left")#os.path.split(m_d)[1].replace("master_",""),loc="left")
        phidf = pd.read_csv(os.path.join(m_d,"freyberg.phi.actual.csv"))
        print(phidf)
        iters = phidf.iteration.values
        [pax.plot(iters,vals,lw=0.1,color="0.5",alpha=0.5) for vals in np.log10(phidf.iloc[:,6:].values.T)]
        pax.set_title(plabel,loc="left")
        pax.set_xlabel("iteration")
        pax.set_ylabel("$log_{10} \\phi$")
        pax.set_xlim(0,pst.control_data.noptmax+1)
        iiter = phidf.iteration.max()
        pes = []
        for i in range(iiter+1):
            fname = os.path.join(m_d,"freyberg.{0}.par.csv".format(i))
            if not os.path.exists(fname):
                break
            pe = pd.read_csv(fname,index_col=0)    
            pes.append(pe)
        for ireal,real in enumerate(pes[-1].index):
            xvals  = [pe.loc[real,"hk1"] for pe in pes]
            yvals  = [pe.loc[real,"rch0"] for pe in pes]
            if ireal == 0:
                ax.plot(xvals,yvals,marker=".",c="k",lw=0.4,alpha=0.2,label="realization trajectory")
            else:
                ax.plot(xvals,yvals,marker=".",c="k",lw=0.4,alpha=0.2)

        xvals  = [pe.loc[:,"hk1"].values.mean() for pe in pes]
        yvals  = [pe.loc[:,"rch0"].values.mean() for pe in pes]
        ax.plot(xvals,yvals,marker="^",c="k",lw=1.5,dashes=(1,1),zorder=10,label="mean trajectory")
        ax.set_xlabel("hydraulic conductivity")
        ax.set_ylabel("recharge")
        xvals = pes[-1].loc[:,"hk1"].values
        yvals = pes[-1].loc[:,"rch0"].values
        ax.scatter(xvals,yvals,marker=".",c="b",alpha=0.5,zorder=10,label="posterior")
        if im_d == 0:
            ax.legend(loc="upper left")
        
    plt.tight_layout()
    plt.savefig("mean_iter_compare.pdf")
    print(m_ds)

def plot_response_surface(parnames=['hk1','rch0'], pstfile='freyberg.pst', WORKING_DIR='freyberg_mf6',
                          nanthresh=None,alpha=0.5, label=True, maxresp=None,
                          figsize=(5,5),levels=None, cmap="nipy_spectral",ax=None):
    NUM_STEPS_RESPSURF = 40
    p1,p2 = parnames
    df_in = pd.read_csv(os.path.join(WORKING_DIR, "sweep_in.csv"))
    df_out = pd.read_csv(os.path.join(WORKING_DIR, "sweep_out.csv"))
    resp_surf = np.zeros((NUM_STEPS_RESPSURF, NUM_STEPS_RESPSURF))
    p1_values = df_in[p1].unique()
    p2_values = df_in[p2].unique()
    c = 0
    for i, v1 in enumerate(p1_values):
        for j, v2 in enumerate(p2_values):
            resp_surf[j, i] = df_out.loc[c, "phi"]
            c += 1
    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = plt.subplot(111)
    X, Y = np.meshgrid(p1_values, p2_values)
    if nanthresh is not None:
        resp_surf = np.ma.masked_where(resp_surf > nanthresh, resp_surf)
    if maxresp is None:
        maxresp = np.max(resp_surf)
    if levels is None:
        levels = np.array([0.001, 0.01, 0.02, 0.05, .1, .2, .5])*maxresp
    
    import matplotlib.colors as colors
    vmin=np.min(resp_surf)
    vmax=maxresp
    p = ax.pcolor(X, Y, resp_surf, alpha=alpha, cmap=cmap,# vmin=vmin, vmax=vmax,
                norm=colors.LogNorm(vmin=vmin, vmax=vmax))
    plt.colorbar(p,label="$log_{10} \\phi$")
    c = ax.contour(X, Y, resp_surf,
                   levels=levels,
                   colors='k', alpha=0.5)
    #plt.title('min $\\Phi$ = {0:.2f}'.format(np.nanmin(resp_surf)))
    #if label:
    plt.clabel(c)
    ax.set_xlim(p1_values.min(), p1_values.max())
    ax.set_ylim(p2_values.min(), p2_values.max())
    ax.set_xlabel(p1)
    ax.set_ylabel(p2)
    return ax, resp_surf



def tenpar_noise_invest():
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_noise_invest")
    template_d = os.path.join(model_d, "test_template")

    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)

    # if os.path.exists(test_d):
    #     shutil.rmtree(test_d)
    pst.parameter_data.loc[:,"partrans"] = "log"
    pst.parameter_data.loc[:,"pargp"] = pst.par_names

    obs = pst.observation_data
    # obs.loc[pst.obs_names[:4],"obgnme"] = "og1a"
    # obs.loc[pst.obs_names[:4],"weight"] = 10 + 3.*(np.random.random(4))
    # obs.loc[pst.obs_names[4:8],"obgnme"] = "og1b"
    # obs.loc[pst.obs_names[4:8],"weight"] = 5
    # obs.loc[pst.obs_names[8:12],"obgnme"] = "og3"
    # obs.loc[pst.obs_names[8:12],"weight"] = 1
    # obs.loc[pst.obs_names[12:],"obgnme"] = "og4"
    # obs.loc[pst.obs_names[12:],"weight"] = 0.00001
    # obs.loc[:,"standard_deviation"] = 0.1
    # #with open(os.path.join(template_d,"phi.csv"),'w') as f:
    # #    f.write("og1,0.333333\n")
    # #    f.write("og3,0.333333\n")
    # #    f.write("og4,0.333333\n")
    # gdict = {}
    # for g in pst.nnz_obs_groups:
    #     gdict[g] = obs.loc[obs.obgnme==g,"obsnme"].to_list()

    # df = pd.DataFrame(columns=["og1","og3","og4"],index=np.arange(12))
    # df.loc[:,:] = 0.33333
    # df.loc[np.arange(3),["og1"]] = 0.9
    # df.loc[np.arange(3),["og3"]] = 0.05
    # df.loc[np.arange(3),["og4"]] = 0.05

    # df.loc[np.arange(3,6),["og1"]] = 0.05
    # df.loc[np.arange(3,6),["og3"]] = 0.9
    # df.loc[np.arange(3,6),["og4"]] = 0.05

    # df.loc[np.arange(6,11),["og1"]] = 0.05
    # df.loc[np.arange(6,11),["og3"]] = 0.05
    # df.loc[np.arange(6,11),["og4"]] = 0.09

    # df.to_csv(os.path.join(template_d,"phi.csv"))

    print(pst.nnz_obs_names)
    noise_adds = np.linspace(-1.0,1.0,100)
    noise_arr = np.zeros((len(noise_adds),pst.nnz_obs))
    for j,val in enumerate(obs.loc[pst.nnz_obs_names,"obsval"]):
        noise_arr[:,j] = val + noise_adds


    print(noise_arr)
    df = pd.DataFrame(noise_arr,columns=pst.nnz_obs_names)
    df.to_csv(os.path.join(template_d,"noise.csv"))


    pst.pestpp_options["ies_no_noise"] = False
    pst.pestpp_options["ies_obs_en"] = "noise.csv"

    pst.pestpp_options["ies_num_reals"] = df.shape[0]
    pst.pestpp_options['ies_verbose_level'] = 4
    pst.pestpp_options["ies_bad_phi_sigma"] = -1.2
    #pst.pestpp_options["ies_multimodal_alpha"] = 0.25
    #pst.pestpp_options["ies_n_iter_mean"] = 3
    
    
    
    pst.control_data.noptmax = 12
    pst_name = "pest_adj.pst"
    pst.write(os.path.join(template_d,pst_name),version=2)
    
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d, test_d)
    pyemu.os_utils.run("{0} {1}".format(exe_path, pst_name), cwd=test_d)


    noise = pd.read_csv(os.path.join(test_d,"pest_adj.obs+noise.csv"),index_col=0)
    meas = pd.read_csv(os.path.join(test_d,"pest_adj.phi.meas.csv"))
    actual = pd.read_csv(os.path.join(test_d,"pest_adj.phi.actual.csv"))
    meas.index = meas.iteration
    actual.index = actual.iteration

    obs_dfs = []
    for itr in meas.iteration:
        obs_dfs.append(pd.read_csv(os.path.join(test_d,"pest_adj.{0}.obs.csv".format(itr)),index_col=0))

    last_df = obs_dfs[-1]
    for i,df in enumerate(obs_dfs):
        obs_dfs[i] = df.loc[last_df.index,:]

    fig,axes = plt.subplots(1,3,figsize=(30,10))

    ax = axes[0]
    o1,o2 = pst.nnz_obs_names
    cmap = plt.get_cmap("jet")
    for ireal,real in enumerate(last_df.index):
        c = cmap(float(ireal)/float(last_df.shape[0]))
        
        ax.scatter(noise.loc[real,o1],noise.loc[real,o2],marker='^',s=30,c=c,zorder=10)
        v1 = [odf.loc[real,o1] for odf in obs_dfs]
        v2 = [odf.loc[real,o2] for odf in obs_dfs]
        ax.plot(v1,v2,marker="x",color=c,lw=0.1,alpha=0.25)
        ax.scatter(last_df.loc[real,o1],last_df.loc[real,o2],marker='o',s=30,c=c,zorder=10)
        ax.annotate("{0:3.3f}".format(meas.loc[:,real].values[-1]),
            (last_df.loc[real,o1],last_df.loc[real,o2]),xytext=(0,-4),textcoords="offset pixels",va="top",ha="center")

        ax.annotate("{0:3.3f}".format(actual.loc[:,real].values[-1]),
            (last_df.loc[real,o1],last_df.loc[real,o2]),xytext=(0,4),textcoords="offset pixels",va="bottom",ha="center")
    ax = axes[1]
    ax.hist(np.log10(meas.iloc[0,6:].values),facecolor="0.5",alpha=0.5)
    ax.hist(np.log10(meas.iloc[-1,6:].values),facecolor="b",alpha=0.5)
    
    ax = axes[2]
    ax.hist(np.log10(actual.iloc[0,6:].values),facecolor="0.5",alpha=0.5)
    ax.hist(np.log10(actual.iloc[-1,6:].values),facecolor="b",alpha=0.5)
    

    plt.tight_layout()
    plt.savefig("noise_compare.pdf")
    

    plt.close(fig)

def poly_n_iter_mean_invest(b_d="poly",use_ineq=False,n_iter_mean=3):
    port = 4343
    if os.path.exists(b_d):
        shutil.rmtree(b_d)
    os.makedirs(b_d)
    t_d = os.path.join(b_d,"template")
    if os.path.exists(t_d):
        shutil.rmtree(t_d)
    os.makedirs(t_d)
    pd.DataFrame({"parval1":1.0},index=["par"]).to_csv(os.path.join(t_d,"par.csv"))
    with open(os.path.join(t_d,"par.csv.tpl"),'w') as f:
        f.write("ptf ~\n")
        f.write("parnme,parval1\n")
        f.write("par,~ par   ~\n")

    with open(os.path.join(t_d,"forward_run.py"),'w') as f:
        f.write("import pandas as pd\n")
        f.write("pdf = pd.read_csv('par.csv')\n")
        f.write("pval = float(pdf.iloc[0,1])\n")
        f.write("sim = pval**4 + pval**3 - pval**2 + 1\n")
        f.write("with open('sim.csv','w') as f:\n")
        f.write("    f.write('obsnme,obsval\\n')\n")
        f.write("    f.write('sim,'+str(sim)+'\\n')\n")
    pyemu.os_utils.run("python forward_run.py",cwd=t_d)
    with open(os.path.join(t_d,"sim.csv.ins"),'w') as f:
        f.write("pif ~\n")
        f.write("l1\n")
        f.write("l1 ~,~ !sim!\n")

    pst = pyemu.Pst.from_io_files([os.path.join(t_d,"par.csv.tpl")],
                                  [os.path.join(t_d,"par.csv")],
                                  [os.path.join(t_d,"sim.csv.ins")],
                                  [os.path.join(t_d,"sim.csv")],
                                  pst_path=".")
    pst.model_command = "python forward_run.py"
    par = pst.parameter_data
    par.loc[:,"parval1"] = 1.0
    par.loc[:,"parlbnd"] = 0.0
    par.loc[:,"parubnd"] = 1.8
    par.loc[:,"partrans"] = "none"

    obs = pst.observation_data
    obs.loc[:,"obsval"] = -0.1
    if use_ineq:
        obs.loc[:,"obsval"] = 0.7
        obs.loc[:,"weight"] = 100.0
        obs.loc[:,"obgnme"] = "less_than"

    pst.control_data.noptmax = 0
    pst.control_data.nphinored = 1000
    #pst.pestpp_options["ies_update_by_reals"] = False
    pst.pestpp_options["ies_bad_phi_sigma"] = 1.75
    pst.write(os.path.join(t_d,"pest.pst"))
    pyemu.os_utils.run("{0} pest.pst".format(exe_path),cwd=t_d)
    
    num_reals = 50
    num_workers = 50
    noptmax = 25
    np.random.seed(11233442)
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst,cov=pyemu.Cov.from_parameter_data(pst),num_reals=num_reals)
    pe.enforce()
    pe.to_csv(os.path.join(t_d,"narrow_prior.csv"))

    

    #pst.pestpp_options["ies_initial_lambda"] = 1000
    #pst.pestpp_options["ies_lambda_dec_fac"] = 1.0
    pst.pestpp_options["ies_update_by_reals"] = True
    #pst.pestpp_options["ies_use_mda"] = True
    
    pst.pestpp_options["ies_par_en"] = "narrow_prior.csv"
    pst.control_data.noptmax = noptmax
    pst.write(os.path.join(t_d,"pest.pst"))
    test1_d = os.path.join(b_d,"narrow_base")
    pyemu.os_utils.start_workers(t_d,exe_path,"pest.pst",
                                 num_workers=num_workers,
                                 master_dir=test1_d,worker_root=b_d,
                                 port=port)
    
    par = pst.parameter_data
    par.loc[:,"parval1"] = 0.2
    par.loc[:,"parlbnd"] = -2.0
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst,cov=pyemu.Cov.from_parameter_data(pst),num_reals=num_reals)
    pe.enforce()
    pe.to_csv(os.path.join(t_d,"wide_prior.csv"))

    pst.pestpp_options["ies_par_en"] = "wide_prior.csv"
    pst.write(os.path.join(t_d,"pest.pst"))
    test2_d = os.path.join(b_d,"wide_base")
    pyemu.os_utils.start_workers(t_d,exe_path,"pest.pst",num_workers=num_workers,
                                 master_dir=test2_d,worker_root=b_d,port=port)


    pst.pestpp_options["ies_par_en"] = "narrow_prior.csv"
    pst.pestpp_options["ies_n_iter_mean"] = n_iter_mean
    pst.write(os.path.join(t_d,"pest.pst"))
    test3_d = os.path.join(b_d,"narrow_nim")
    pyemu.os_utils.start_workers(t_d,exe_path,"pest.pst",num_workers=num_workers,
                                 master_dir=test3_d,worker_root=b_d,port=port)
    
    pst.pestpp_options["ies_multimodal_alpha"] = 0.99
    pst.write(os.path.join(t_d,"pest.pst"))
    test8_d = os.path.join(b_d,"narrow_nim_mm99")
    pyemu.os_utils.start_workers(t_d,exe_path,"pest.pst",num_workers=num_workers,
                                 master_dir=test8_d,worker_root=b_d,port=port)
    
    pst.pestpp_options["ies_multimodal_alpha"] = 0.99
    pst.pestpp_options["ies_n_iter_mean"] = 0
    pst.write(os.path.join(t_d,"pest.pst"))
    test4_d = os.path.join(b_d,"narrow_mm99")
    pyemu.os_utils.start_workers(t_d,exe_path,"pest.pst",num_workers=num_workers,
                                 master_dir=test4_d,worker_root=b_d,port=port)
    
    pst.pestpp_options["ies_par_en"] = "wide_prior.csv"
    pst.write(os.path.join(t_d,"pest.pst"))
    test5_d = os.path.join(b_d,"wide_mm99")
    pyemu.os_utils.start_workers(t_d,exe_path,"pest.pst",num_workers=num_workers,
                                 master_dir=test5_d,worker_root=b_d,port=port)
    
    pst.pestpp_options["ies_n_iter_mean"] = n_iter_mean
    pst.write(os.path.join(t_d,"pest.pst"))
    test6_d = os.path.join(b_d,"wide_nim_mm99")
    pyemu.os_utils.start_workers(t_d,exe_path,"pest.pst",num_workers=num_workers,
                                 master_dir=test6_d,worker_root=b_d,port=port)
    
    pst.pestpp_options["ies_par_en"] = "wide_prior.csv"
    pst.pestpp_options["ies_n_iter_mean"] = n_iter_mean
    pst.pestpp_options["ies_multimodal_alpha"] = 1.0
    pst.write(os.path.join(t_d,"pest.pst"))
    test7_d = os.path.join(b_d,"wide_nim")
    pyemu.os_utils.start_workers(t_d,exe_path,"pest.pst",num_workers=num_workers,
                                 master_dir=test7_d,worker_root=b_d,port=port)


def plot_poly(b_d="poly"):
    
    m_ds = [os.path.join(b_d,d) for d in os.listdir(b_d) if os.path.isdir(os.path.join(b_d,d)) and d not in ["template","plot"]]

    def plot_iter(pes,oes,colors,ax,xmin,xmax,ymin,ymax,ineq_level=None):
        x = np.linspace(xmin,xmax,300)

        y = x**4 + x**3 - x**2 + 1
        idx = np.where(y==y.min())
        xmin = x[idx]
        ax.plot(x,y,"k",lw=2.0,alpha=0.5,label="response surface")
        if ineq_level is None:
            ax.scatter(x[idx],y[idx],marker="*",s=70,c='r',zorder=10,label="optimum")
        for oe,pe,color in zip(oes,pes,colors):
            ax.scatter(pe.values,oe.values,marker=".",s=40,c=color,label="realization")
        if ineq_level is not None:
            ineq_idxs = np.where(y<=ineq_level)
            yineq = y[ineq_idxs]
            xineq = x[ineq_idxs]
            
            ineq_level = np.zeros_like(xineq) + ineq_level
            #print(xineq.shape,yineq.shape,ineq_level.shape)
            ax.fill_between(xineq,yineq,ineq_level,facecolor="r",alpha=0.75,label="ineq satisfied")
        ax.set_xlim(x.min(),x.max())
        ax.set_ylim(ymin,ymax)
        ax.legend(loc="upper left")
        ax.set_ylabel("objective function")
        ax.set_xlabel("parameter")
        ax.grid()
        
    results = {}
    m_ds.sort()
    reorder_m_ds = []
    for m_d in m_ds:
        if "narrow_" in os.path.split(m_d)[-1]:
            reorder_m_ds.append(m_d)
    i = 1
    for m_d in m_ds:
        if "wide_" in os.path.split(m_d)[-1]:
            reorder_m_ds.insert(i,m_d)
            i += 2
    #print(m_ds)
    m_ds = reorder_m_ds
    #print(m_ds)
    omin,omax = 1e30,-1e30
    for m_d in m_ds:
        oes,pes = [],[]
        phidf = pd.read_csv(os.path.join(m_d,"pest.phi.actual.csv"))
        pst = pyemu.Pst(os.path.join(m_d,"pest.pst"))
        for i in range(phidf.iteration.max()+1):
            oe = pd.read_csv(os.path.join(m_d,"pest.{0}.obs.csv".format(i)),index_col=0)
            omin = min(oe.values.min(),omin)
            omax = max(oe.values.max(),omax)
            oes.append(oe)
            pes.append(pd.read_csv(os.path.join(m_d,"pest.{0}.par.csv".format(i)),index_col=0))
           
        results[os.path.split(m_d)[-1]] = [pes,oes,phidf]
    
    par = pst.parameter_data
    pmin = par.loc["par","parlbnd"]
    pmax = par.loc["par","parubnd"]
    omin,omax = -0.5,14
    pmin,pmax = -2,2
    plt_d = os.path.join(b_d,"plot")
    if os.path.exists(plt_d):
        shutil.rmtree(plt_d)
    os.makedirs(plt_d)

    ineq_level = None
    obs = pst.observation_data
    if obs.obgnme.iloc[0].startswith("less"):
        ineq_level = obs.obsval.iloc[0]
    def get_tag(_m_d):
        tag = ""
        if "nim" in _m_d:
            tag += "cov re-inflation, "
        if "mm" in _m_d:
            tag += "realization loc"
        return tag
    m_ds = [os.path.split(m_d)[-1] for m_d in m_ds]
    fig,axes = plt.subplots(int(len(m_ds)/2),2,figsize=(6,8))
    axes = axes.flatten()
    mx = -1e30
    for m_d,ax in zip(m_ds,axes):
        phidf = results[m_d][2]
        vals = phidf.iloc[:,6:].values
        navals = vals.copy()
        navals[navals==0] = np.nan
        vals[vals==0] = np.nanmin(navals)
        vals = np.log10(vals)
        itr = phidf.iteration.values
        #print(vals)
        [ax.plot(itr,vals[:,i],"0.5",alpha=0.5,lw=0.5) for i in range(vals.shape[1])]
        tag = get_tag(m_d)
        label = "{1} prior ensemble\n{2}".format(i,m_d.split('_')[0],tag)
        ax.set_title(label,loc="left")
        mx = max(ax.get_ylim()[1],mx)
    for ax in axes:
        ax.set_ylim(0,mx)
        ax.set_xlabel("iteration")
        ax.set_ylabel("$log_{10} \\phi$")
    plt.tight_layout()
    plt.savefig(os.path.join(plt_d,"phi.pdf"))
    plt.close(fig)

    iplt = 0
    for i in range(pst.control_data.noptmax):
        fig,axes = plt.subplots(int(len(m_ds)/2),2,figsize=(6,8))
        axes = axes.flatten()
        for m_d,ax in zip(m_ds,axes):
            #print(i,m_d)
            try:
                pe,oe = results[m_d][0][i],results[m_d][1][i]
            except Exception as e:
                pe = results[m_d][0][results[m_d][2].iteration.max()]
                oe = results[m_d][1][results[m_d][2].iteration.max()]
                
            #    shutil.copy2(os.path.join(plt_d,"fig_{0:03d}.png".format(i-1)),
            #        os.path.join(plt_d,"fig_{0:03d}.png".format(i)))
            #else:
                pass
            plot_iter([pe],[oe],["b"],ax,pmin,pmax,omin,omax,ineq_level=ineq_level)
            tag = get_tag(m_d)
            label = "iteration {0}, {1} prior ensemble\n{2}".format(i,m_d.split('_')[0],tag)
            ax.set_title(label,loc="left")
            # if ineq_level is not None:
            #     xlim = ax.get_xlim()
            #     ax.plot(xlim,[ineq_level,ineq_level],"r--",label="less than")
            #     ax.legend(loc="upper left")
        plt.tight_layout()
        plt.savefig(os.path.join(plt_d,"fig_{0:03d}.png".format(iplt)),dpi=500)
        iplt += 1
        if i == 0:
            for ii in range(5):
                plt.savefig(os.path.join(plt_d,"fig_{0:03d}.png".format(iplt)),dpi=500)
                iplt += 1
        plt.close(fig)
    fps = 1

    pyemu.os_utils.run("ffmpeg -y -i fig_{0:03d}.png -vf palettegen=256 palette.png".format(i),cwd=plt_d)
    pyemu.os_utils.run("ffmpeg -r {0} -y -s 1920X1080 -i fig_%03d.png -i palette.png -filter_complex \"scale=720:-1:flags=lanczos[x];[x][1:v]paletteuse\" -final_delay 150 logo.gif".format(fps),
            cwd=plt_d)


def hosaki_invest(b_d="hosaki",use_ineq=False,n_iter_mean=3,bad_phi_sigma=None):
    port = 4343
    #if os.path.exists(b_d):
    #    shutil.rmtree(b_d)
    #os.makedirs(b_d)
    if not os.path.exists(b_d):
        os.makedirs(b_d)
    t_d = os.path.join(b_d,"template")
    if os.path.exists(t_d):
        shutil.rmtree(t_d)
    os.makedirs(t_d)
    pd.DataFrame({"parval1":[4,2]},index=["par1","par2"]).to_csv(os.path.join(t_d,"par.csv"))
    with open(os.path.join(t_d,"par.csv.tpl"),'w') as f:
        f.write("ptf ~\n")
        f.write("parnme,parval1\n")
        f.write("par1,~   par1   ~\n")
        f.write("par2,~   par2   ~\n")

    with open(os.path.join(t_d,"forward_run.py"),'w') as f:
        f.write("import pandas as pd\n")
        f.write("import numpy as np\n")
        f.write("pdf = pd.read_csv('par.csv')\n")
        f.write("pval1 = float(pdf.iloc[0,1])\n")
        f.write("pval2 = float(pdf.iloc[1,1])\n")
        f.write("term1 = (pval2**2)*np.e**(-pval2)\n")
        f.write("term2 = 1. - (8.*(pval1)) + (7.*(pval1**2))\n")
        f.write("term3 = (-(7./3.)*(pval1**3)) + (0.25*(pval1**4))\n")
        f.write("sim = (term2 + term3) * term1\n")
        f.write("with open('sim.csv','w') as f:\n")
        f.write("    f.write('obsnme,obsval\\n')\n")
        f.write("    f.write('sim,'+str(sim)+'\\n')\n")
    pyemu.os_utils.run("python forward_run.py",cwd=t_d)
    with open(os.path.join(t_d,"sim.csv.ins"),'w') as f:
        f.write("pif ~\n")
        f.write("l1\n")
        f.write("l1 ~,~ !sim!\n")

    pst = pyemu.Pst.from_io_files([os.path.join(t_d,"par.csv.tpl")],
                                  [os.path.join(t_d,"par.csv")],
                                  [os.path.join(t_d,"sim.csv.ins")],
                                  [os.path.join(t_d,"sim.csv")],
                                  pst_path=".")
    pst.model_command = "python forward_run.py"
    par = pst.parameter_data
    par.loc[:,"parval1"] = 1.0
    par.loc[:,"parlbnd"] = 0.0
    par.loc[:,"parubnd"] = 5.0
    par.loc[:,"partrans"] = "none"

    obs = pst.observation_data
    obs.loc[:,"obsval"] = -2.345811576101292
    obs.loc[:,"weight"] = 100
    
    if use_ineq:
        obs.loc[:,"obsval"] = -1.5
        obs.loc[:,"weight"] = 100.0
        obs.loc[:,"obgnme"] = "less_than"

    pst.control_data.noptmax = 0
    pst.control_data.nphinored = 1000
    #pst.pestpp_options["ies_update_by_reals"] = False
    if bad_phi_sigma is not None:
        pst.pestpp_options["ies_bad_phi_sigma"] = 1.75
    pst.write(os.path.join(t_d,"pest.pst"))
    pyemu.os_utils.run("{0} pest.pst".format(exe_path),cwd=t_d)
    
    par.loc[:,"parval1"] = [4,2]
    pst.write(os.path.join(t_d,"pest.pst"))
    pyemu.os_utils.run("{0} pest.pst".format(exe_path),cwd=t_d)
    pst = pyemu.Pst(os.path.join(t_d,"pest.pst"))
    print(pst.phi)
    assert pst.phi < 1e-10
   
    pvals = []
    for p1 in np.linspace(0,5,100):
        for p2 in np.linspace(0,5,100):
            pvals.append([p1,p2])
    pvals = pd.DataFrame(pvals,columns=["par1","par2"])
    pvals.to_csv(os.path.join(t_d,"sweep_in.csv"))
    pst.pestpp_options["ies_par_en"] = "sweep_in.csv"
    pst.control_data.noptmax = -1
    sweep_d = os.path.join(b_d,"sweep")
    pst.pestpp_options["ies_include_base"] = False

    pst.write(os.path.join(t_d,"pest.pst"))
    
    if not os.path.exists(sweep_d):
        num_workers = 200

        pyemu.os_utils.start_workers(t_d,exe_path,"pest.pst",
                                     num_workers=num_workers,
                                     master_dir=sweep_d,worker_root=b_d,
                                     port=port)

    par = pst.parameter_data
    par.loc["par1",["parlbnd","parubnd"]] = [0,2]
    par.loc["par2",["parlbnd","parubnd"]] = [2,5]
    par.loc[:,"parval1"] = [1.,3.5]


    num_reals = 50
    num_workers = 50
    noptmax = 25
    np.random.seed(11233442)
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst,cov=pyemu.Cov.from_parameter_data(pst),num_reals=num_reals)
    pe.enforce()
    pe.to_csv(os.path.join(t_d,"narrow_prior.csv"))

    par = pst.parameter_data
    par.loc["par1",["parlbnd","parubnd"]] = [0,5]
    par.loc["par2",["parlbnd","parubnd"]] = [0,5]
    par.loc[:,"parval1"] = [2.5,2.5]

    #pst.pestpp_options["ies_initial_lambda"] = 1000
    #pst.pestpp_options["ies_lambda_dec_fac"] = 1.0
    pst.pestpp_options["ies_update_by_reals"] = False
    #pst.pestpp_options["ies_use_mda"] = True
    
    pst.pestpp_options["ies_par_en"] = "narrow_prior.csv"

    pst.control_data.noptmax = noptmax
    pst.write(os.path.join(t_d,"pest.pst"))
    test1_d = os.path.join(b_d,"narrow_base")
    pyemu.os_utils.start_workers(t_d,exe_path,"pest.pst",
                                 num_workers=num_workers,
                                 master_dir=test1_d,worker_root=b_d,
                                 port=port)
    
    
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst,cov=pyemu.Cov.from_parameter_data(pst),num_reals=num_reals)
    pe.enforce()
    pe.to_csv(os.path.join(t_d,"wide_prior.csv"))

    pst.pestpp_options["ies_par_en"] = "wide_prior.csv"
    pst.write(os.path.join(t_d,"pest.pst"))
    test2_d = os.path.join(b_d,"wide_base")
    pyemu.os_utils.start_workers(t_d,exe_path,"pest.pst",num_workers=num_workers,
                                 master_dir=test2_d,worker_root=b_d,port=port)


    pst.pestpp_options["ies_par_en"] = "narrow_prior.csv"
    pst.pestpp_options["ies_n_iter_mean"] = n_iter_mean
    pst.write(os.path.join(t_d,"pest.pst"))
    test3_d = os.path.join(b_d,"narrow_nim")
    pyemu.os_utils.start_workers(t_d,exe_path,"pest.pst",num_workers=num_workers,
                                 master_dir=test3_d,worker_root=b_d,port=port)
    
    pst.pestpp_options["ies_multimodal_alpha"] = 0.99
    pst.write(os.path.join(t_d,"pest.pst"))
    test8_d = os.path.join(b_d,"narrow_nim_mm99")
    pyemu.os_utils.start_workers(t_d,exe_path,"pest.pst",num_workers=num_workers,
                                 master_dir=test8_d,worker_root=b_d,port=port)
    
    pst.pestpp_options["ies_multimodal_alpha"] = 0.99
    pst.pestpp_options["ies_n_iter_mean"] = 0
    pst.write(os.path.join(t_d,"pest.pst"))
    test4_d = os.path.join(b_d,"narrow_mm99")
    pyemu.os_utils.start_workers(t_d,exe_path,"pest.pst",num_workers=num_workers,
                                 master_dir=test4_d,worker_root=b_d,port=port)
    
    pst.pestpp_options["ies_par_en"] = "wide_prior.csv"
    pst.write(os.path.join(t_d,"pest.pst"))
    test5_d = os.path.join(b_d,"wide_mm99")
    pyemu.os_utils.start_workers(t_d,exe_path,"pest.pst",num_workers=num_workers,
                                 master_dir=test5_d,worker_root=b_d,port=port)
    
    pst.pestpp_options["ies_n_iter_mean"] = n_iter_mean
    pst.write(os.path.join(t_d,"pest.pst"))
    test6_d = os.path.join(b_d,"wide_nim_mm99")
    pyemu.os_utils.start_workers(t_d,exe_path,"pest.pst",num_workers=num_workers,
                                 master_dir=test6_d,worker_root=b_d,port=port)
    
    pst.pestpp_options["ies_par_en"] = "wide_prior.csv"
    pst.pestpp_options["ies_n_iter_mean"] = n_iter_mean
    pst.pestpp_options["ies_multimodal_alpha"] = 1.0
    pst.write(os.path.join(t_d,"pest.pst"))
    test7_d = os.path.join(b_d,"wide_nim")
    pyemu.os_utils.start_workers(t_d,exe_path,"pest.pst",num_workers=num_workers,
                                 master_dir=test7_d,worker_root=b_d,port=port)


def plot_hosaki(b_d="hosaki",steps=100):
    sweep_d = os.path.join(b_d,"sweep")
    assert os.path.exists(sweep_d)
    sweep_oe = pd.read_csv(os.path.join(sweep_d,"pest.0.obs.csv"),index_col=0)
    sweep_pe = pd.read_csv(os.path.join(sweep_d,"pest.0.par.csv"),index_col=0)
    
    sweep_x = sweep_pe.loc[:,"par1"].values.reshape(steps,steps)
    sweep_y = sweep_pe.loc[:,"par2"].values.reshape(steps,steps)
    sweep_z = sweep_oe.loc[:,"sim"].values.reshape(steps,steps)

    m_ds = [os.path.join(b_d,d) for d in os.listdir(b_d) if os.path.isdir(os.path.join(b_d,d)) and "narrow" in d]
    results = {}
    narrow_m_ds = []
    for m_d in m_ds:
        phidf = pd.read_csv(os.path.join(m_d,"pest.phi.actual.csv"))
        pes = []
        for i in range(phidf.iteration.max()+1):
            pe = pd.read_csv(os.path.join(m_d,"pest.{0}.par.csv".format(i)),index_col=0)
            pes.append(pe)
        results[m_d] = [pes,phidf]
        narrow_m_ds.append(m_d)
        m_d = m_d.replace("narrow","wide")
        print(m_d)
        assert os.path.exists(m_d)
        phidf = pd.read_csv(os.path.join(m_d,"pest.phi.actual.csv"))
        pes = []
        for i in range(phidf.iteration.max()):
            pe = pd.read_csv(os.path.join(m_d,"pest.{0}.par.csv".format(i)),index_col=0)
            pes.append(pe)
        results[m_d] = [pes,phidf]

    narrow_m_ds.sort()

    for m_d in narrow_m_ds:
        plt_d = os.path.join("plot",b_d,os.path.split(m_d)[-1])
        if os.path.exists(plt_d):
            shutil.rmtree(plt_d)
        os.makedirs(plt_d)

        npes,nphidf = results[m_d]
        wpes,wphidf = results[m_d.replace("narrow","wide")]
        itrs = max(nphidf.iteration.max(),wphidf.iteration.max())
        iplt = 0
        for i in range(itrs+1):
            tag = os.path.split(m_d)[-1].split("_")[0]
            if "nim" in m_d:
                tag += ", re-inflation"
            if "mm" in m_d:
                tag += ", realization local"
            tag += "\n iteration {0}".format(i)
            fig,axes = plt.subplots(1,2,figsize=(9,4))
            cb = axes[0].pcolormesh(sweep_x,sweep_y,sweep_z,alpha=0.5,cmap="magma")
            axes[0].contour(sweep_x,sweep_y,sweep_z,levels=[-2,-1],colors='k',linewidths=0.25,linestyles="solid")
            plt.colorbar(cb,ax=axes[0],label="objective function")
            if i >= len(npes):
                pe = npes[-1]
            else:
                pe = npes[i]
            axes[0].scatter(pe.par1.values,pe.par2.values,marker=".",s=10,c="w")
            axes[0].set_title(tag,loc="left")
            cb = axes[1].pcolormesh(sweep_x,sweep_y,sweep_z,alpha=0.5,cmap="magma")
            plt.colorbar(cb,ax=axes[1],label="objective function")
            if i >= len(wpes):
                pe = wpes[-1]
            else:
                pe = wpes[i]
            axes[1].contour(sweep_x,sweep_y,sweep_z,levels=[-2,-1],colors='k',linewidths=0.25,linestyles="solid")   
            axes[1].scatter(pe.par1.values,pe.par2.values,marker=".",s=10,c="w")
            axes[1].set_title(tag.replace("narrow","wide"),loc="left")
            for ax in axes:
                ax.set_yticks([])
                ax.set_xticks([])
                ax.set_aspect("equal")
            
            plt.tight_layout()
            if i == 0:
                for _ in range(4):
                    plt.savefig(os.path.join(plt_d,"fig_{0:03d}.png".format(iplt)),dpi=500)
                    iplt += 1
            plt.savefig(os.path.join(plt_d,"fig_{0:03d}.png".format(iplt)),dpi=500)
            iplt += 1
            plt.close(fig)
        fps = 2
        pyemu.os_utils.run("ffmpeg -y -i fig_{0:03d}.png -vf palettegen=256 palette.png".format(i),cwd=plt_d)
        pyemu.os_utils.run("ffmpeg -r {0} -y -s 1920X1080 -i fig_%03d.png -i palette.png -filter_complex \"scale=720:-1:flags=lanczos[x];[x][1:v]paletteuse\" -final_delay 150 logo.gif".format(fps),
            cwd=plt_d)




def tenpar_fixed_restart_test():

    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_restart_fixed")
    template_d = os.path.join(model_d, "test_template")
    
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d,test_d)
    pst_name = "pest.pst"
    pst = pyemu.Pst(os.path.join(template_d,pst_name))
    pst.pestpp_options["ies_n_iter_mean"] = [-1,-3,5,999]
    #pst.pestpp_options["ies_reinflate_factor"] = [1.0,0.9,0.8,0.7]
    
    pst.pestpp_options["ies_initial_lambda"] = -100
    pst.control_data.noptmax = 21 # hard coded to test results below
    pst.pestpp_options["ies_num_reals"] = 50
    pst.pestpp_options["ies_debug_fail_remainder"] = True
    #pst.pestpp_options["ies_debug_fail_subset"] = True
    pst.pestpp_options["ies_no_noise"] = False

    par = pst.parameter_data
    par.loc[pst.par_names[0],"partrans"] = "fixed"

    onames = pst.obs_names
    obs = pst.observation_data
    obs.loc[:,"weight"] = 1.0
    
    obs.loc[onames[:3],"obgnme"] = "less_than"
    obs.loc[onames[:3],"obsval"] = obs.loc[onames[:3],"obsval"] + 0.1
    obs.loc[onames[3:5],"obgnme"] = "first_group"
    obs.loc[onames[5:10],"obgnme"] = "second_group"
    obs.loc[onames[10:],"obgnme"] = "greater_than"
    obs.loc[onames[10:],"obsval"] = obs.loc[onames[10:],"obsval"] 
    
    facs = {"tag":["less_than","first_group","second_group","greater_than"], "prop":[0.2,0.2,0.1,0.5]}

    df = pd.DataFrame(data=facs)
    df.to_csv(os.path.join(test_d,"phi.csv"),index=False,header=False)
    pst.pestpp_options["ies_phi_factor_file"] = "phi.csv"
    #pst.pestpp_options["debug_parse_only"] = True
    pst.pestpp_options["ies_verbose_level"] = 2
    #pst.pestpp_options["ies_include_base"] = False
    pst.control_data.noptmax = 1
    
    pst.write(os.path.join(test_d,pst_name),version=2)
    pyemu.os_utils.run("{0} pest.pst".format(exe_path),cwd=test_d)

    #shutil.copy2(os.path.join(test_d,"pest.1.par.csv"),os.path.join(test_d,"pest.1.par.csv"))
    #pst.pestpp_options["ies_restart_par_en"] = "pest.1.par.csv"
    #shutil.copy2(os.path.join(test_d,"pest.0.par.csv"),os.path.join(test_d,"pest.0.par.csv"))
    pst.pestpp_options["ies_par_en"] = "pest.1.par.csv"
    #shutil.copy2(os.path.join(test_d,"pest.1.obs.csv"),os.path.join(test_d,"pest.1.obs.csv"))
    #pst.pestpp_options["ies_restart_obs_en"] = "pest.1.obs.csv"
    #shutil.copy2(os.path.join(test_d,"pest.obs+noise.csv"),os.path.join(test_d,"pest.obs+noise.csv"))
    #pst.pestpp_options["ies_obs_en"] = "pest.obs+noise.csv"

    pst.parameter_data.loc[pst.adj_par_names[0],"partrans"] = "fixed"
    
    pst.control_data.noptmax = 1
    pst.write(os.path.join(test_d,"pest1.pst"),version=2)
    pyemu.os_utils.run("{0} pest1.pst".format(exe_path),cwd=test_d)

    oetest1 = pd.read_csv(os.path.join(test_d,"pest.1.obs.csv"),index_col=0)
    oetest1.index = [str(i) for i in oetest1.index]
    oetest2 = pd.read_csv(os.path.join(test_d,"pest1.0.obs.csv"),index_col=0)
    oetest2.index = [str(i) for i in oetest2.index]
    #print(oetest1)
    #print(oetest2)
    diff = np.abs((oetest1 - oetest2).dropna().values)
    #print(diff)
    assert diff.sum() < 1e-6

    it0_base_par = pyemu.pst_utils.read_parfile(os.path.join(test_d,"pest1.0.base.par"))
    pe = pd.read_csv(os.path.join(test_d,"pest.1.par.csv"),index_col=0)
    
    print(pe.loc["base",:])
    diff = (it0_base_par["parval1"] - pe.loc["base",:])
    print(diff)
    assert np.abs(diff.values).max() < 1e-6

    it1_base_par = pyemu.pst_utils.read_parfile(os.path.join(test_d,"pest1.1.base.par"))
    fpars = par.loc[par.partrans=="fixed","parnme"].values
    diff = (it1_base_par.loc[fpars,"parval1"] - pe.loc["base",fpars])
    print(diff)
    assert np.abs(diff.values).max() < 1e-6


    #oe = pd.read_csv(os.path.join(test_d,"restart_obs.csv"),index_col=0)
    oe = pd.read_csv(os.path.join(test_d,"pest1.1.obs.csv"),index_col=0)
    oe.index = [str(i) for i in oe.index]
    r48_ovals = oe.loc["46",pst.obs_names]


    pst.control_data.noptmax = -2
    pst.pestpp_options.pop("ies_restart_obs_en",None)
    pst.pestpp_options.pop("ies_restart_par_en",None)
    pst.pestpp_options.pop("ies_obs_en",None)
    pst.pestpp_options["ies_par_en"] = "pest1.1.par.csv"
    pst.pestpp_options["ies_run_realname"] = "46"
    pst.write(os.path.join(test_d,"pest2.pst"),version=2)
    pyemu.os_utils.run("{0} pest2.pst".format(exe_path),cwd=test_d)
    pst.set_res(os.path.join(test_d,"pest2.46.rei"))

    d1 = r48_ovals - pst.res.modelled.loc[pst.obs_names]
    print(d1)
    assert np.abs(d1.values).sum() < 1.e-7


    

def tenpar_consistency_test():

    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_consist1")
    template_d = os.path.join(model_d, "test_template")
    
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d,test_d)
    pst_name = "pest.pst"
    pst = pyemu.Pst(os.path.join(template_d,pst_name))
    #pst.pestpp_options["ies_n_iter_mean"] = [-1,-3,5,999]
    #pst.pestpp_options["ies_reinflate_factor"] = [1.0,0.9,0.8,0.7]
    
    pst.pestpp_options["ies_initial_lambda"] = -100
    pst.control_data.noptmax = 6 # hard coded to test results below
    pst.pestpp_options["ies_num_reals"] = 50
    pst.pestpp_options["ies_debug_fail_remainder"] = True
    pst.pestpp_options["ies_debug_fail_subset"] = True
    pst.pestpp_options["ies_no_noise"] = False
    pst.pestpp_options["save_binary"] = True
    pst.pestpp_options["ies_save_lambda_en"] = True

    par = pst.parameter_data
    par.loc[pst.par_names[0],"partrans"] = "fixed"

    onames = pst.obs_names
    obs = pst.observation_data
    obs.loc[:,"weight"] = 1.0
    
    obs.loc[onames[:3],"obgnme"] = "less_than"
    obs.loc[onames[:3],"obsval"] = obs.loc[onames[:3],"obsval"] + 0.1
    obs.loc[onames[3:5],"obgnme"] = "first_group"
    obs.loc[onames[5:10],"obgnme"] = "second_group"
    obs.loc[onames[10:],"obgnme"] = "greater_than"
    obs.loc[onames[10:],"obsval"] = obs.loc[onames[10:],"obsval"] 
    
    facs = {"tag":["less_than","first_group","second_group","greater_than"], "prop":[0.2,0.2,0.1,0.5]}

    df = pd.DataFrame(data=facs)
    df.to_csv(os.path.join(test_d,"phi.csv"),index=False,header=False)
    #pst.pestpp_options["ies_phi_factor_file"] = "phi.csv"
    #pst.pestpp_options["debug_parse_only"] = True
    pst.pestpp_options["ies_verbose_level"] = 2
    pst.control_data.noptmax = 6
    
    pst.write(os.path.join(test_d,pst_name),version=2)
    pyemu.os_utils.run("{0} pest.pst".format(exe_path),cwd=test_d)

    new_test_d = os.path.join(model_d, "template_consist2")
    if os.path.exists(new_test_d):
        shutil.rmtree(new_test_d)
    shutil.copytree(test_d,new_test_d)
    m_d = new_test_d.replace("template","master")
    #pyemu.os_utils.run("{0} pest.pst".format(exe_path),cwd=new_test_d)
    pyemu.os_utils.start_workers(new_test_d,exe_path,"pest.pst",worker_root=model_d,num_workers=10,
        master_dir=m_d)
    new_test_d = m_d

    jcb_files = [f for f in os.listdir(test_d) if f.endswith(".jcb")]

    for jcb_file in jcb_files:
        m1 = pyemu.Matrix.from_binary(os.path.join(test_d,jcb_file)).to_dataframe()
        m2 = pyemu.Matrix.from_binary(os.path.join(new_test_d,jcb_file)).to_dataframe()
        m1.loc[:,:] = np.abs(m1.values)
        m2.loc[:,:] = np.abs(m2.values)
        
        diff = m1 - m2
        print(diff.sum,diff.max())
        assert diff.values.sum() == 0.0

def tenpar_uniformdist_invest():

    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_unidist")
    template_d = os.path.join(model_d, "test_template")
    
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(template_d,test_d)
    pst_name = "pest.pst"
    pst = pyemu.Pst(os.path.join(template_d,pst_name))
    #pst.pestpp_options["ies_n_iter_mean"] = [-1,-3,5,999]
    #pst.pestpp_options["ies_reinflate_factor"] = [1.0,0.9,0.8,0.7]
    
    #pst.pestpp_options["ies_initial_lambda"] = -100
    pst.control_data.noptmax = 10
    pst.pestpp_options["ies_num_reals"] = 500
    pst.pestpp_options["ies_bad_phi_sigma"] = 1.5
    
    pst.pestpp_options["ies_no_noise"] = False
    
    par = pst.parameter_data
    par.loc[:,"partrans"] = "none"
    
    draws = []
    for lb,ub in zip(par.parlbnd,par.parubnd):
        vals = np.random.uniform(lb,ub,pst.pestpp_options["ies_num_reals"])
        draws.append(vals)
    pe = pd.DataFrame(data=draws,index=par.parnme.values).T
    pe.to_csv(os.path.join(test_d,"uni_prior.csv"))
    pst.pestpp_options["ies_par_en"] = "uni_prior.csv"

    onames = pst.obs_names
    obs = pst.observation_data
    obs.loc[:,"weight"] = 1.0
    

    pst.pestpp_options["ies_verbose_level"] = 2
    pst.control_data.noptmax = 8
    
    pst.write(os.path.join(test_d,pst_name),version=2)
    pyemu.os_utils.run("{0} pest.pst".format(exe_path),cwd=test_d)

    test_d2 = test_d +"_reinflatemm"
    if os.path.exists(test_d2):
        shutil.rmtree(test_d2)
    shutil.copytree(test_d,test_d2)

    pst.pestpp_options["ies_n_iter_reinflate"] = [-1,-2,999]
    pst.pestpp_options["ies_multimodal_alpha"] = 0.1
    pst.write(os.path.join(test_d2,pst_name),version=2)
    pyemu.os_utils.run("{0} pest.pst".format(exe_path),cwd=test_d2)

    pst.pestpp_options = {"ies_num_reals":pst.pestpp_options["ies_num_reals"],
                          "ies_par_en":"uni_prior.csv"}

    pst.pestpp_options["ies_use_mda"] = True
    test_d3 = test_d +"_mda"
    if os.path.exists(test_d3):
        shutil.rmtree(test_d3)
    shutil.copytree(test_d,test_d3)
    pst.write(os.path.join(test_d3,pst_name),version=2)
    pyemu.os_utils.run("{0} pest.pst".format(exe_path),cwd=test_d3)



def temp_plot():
    from matplotlib.backends.backend_pdf import PdfPages
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_unidist")
    test_d2 = test_d +"_reinflatemm"
    test_d3 = test_d +"_mda"


    m_ds = [test_d,test_d2,test_d3]
    phidfs = []
    pardfs = []
    for m_d in m_ds:
        phidf = pd.read_csv(os.path.join(m_d,"pest.phi.actual.csv"))
        phidfs.append(phidf)
        pardf = pd.read_csv(os.path.join(m_d,"pest.{0}.par.csv".format(phidf.iteration.max())),index_col=0)
        pardfs.append(pardf)

    with PdfPages("compare.pdf") as pdf:
        fig,axes = plt.subplots(3,1,figsize=(8,8))
        for ax,m_d,phidf,color in zip(axes,m_ds,phidfs,["m","c","g"]):
            vals = np.log10(phidf.iloc[:,6:].values)
            itrs = phidf.iteration.values
            print(vals.shape)
            [ax.plot(itrs,vals[:,i],lw=0.05,color=color) for i in range(vals.shape[1])]
            ax.plot(itrs,vals[:,-1],lw=0.1,color=color,label=m_d)
        #ax.legend(loc="upper right")
            ax.set_title("phi "+m_d,loc="left")
        ymax = 0
        for ax in axes:
            ymax = max(ymax,ax.get_ylim()[1])
        for ax in axes:
            ax.set_ylim(0,ymax)
        plt.tight_layout()
        pdf.savefig()
        plt.close(fig)

        fig,ax = plt.subplots(1,1,figsize=(10,6))
        for m_d,phidf,color in zip(m_ds,phidfs,["m","c","g"]):
            vals = np.log10(phidf.iloc[-1,6:].values)
            ax.hist(vals,bins=30,alpha=0.3,facecolor=color,label=m_d)
        ax.legend(loc="upper right")
        ax.set_title("posterior phi")
        plt.tight_layout()
        pdf.savefig()
        plt.close(fig)

        pnames = pardf.columns.values
        for pname in pnames:
            fig,ax = plt.subplots(1,1,figsize=(8,8))
            for m_d,pardf,color in zip(m_ds,pardfs,["m","c","g"]):
                ax.hist(pardf.loc[:,pname].values,alpha=0.3,bins=30,facecolor=color,label=m_d)
            ax.set_title(pname,loc="left")
            ax.legend(loc="upper right")
            plt.tight_layout()
            pdf.savefig()
            plt.close(fig)



def freyberg_mean_invest():
    import flopy
    model_d = "ies_freyberg"
    
    template_d = os.path.join(model_d, "template")
    # print("loading pst")
    pst = pyemu.Pst(os.path.join(template_d, "pest.pst"))
    pst.pestpp_options = {"ies_num_reals": 30}
    
    pst.control_data.noptmax = 3
    pst.observation_data.loc[:,"weight"] = 1.0
    pst.write(os.path.join(template_d, "pest.pst"))

    test_d = os.path.join(model_d, "master_mean_small")
    if os.path.exists(test_d):
       shutil.rmtree(test_d)
    pyemu.os_utils.start_workers(template_d, exe_path, "pest.pst", num_workers=10, master_dir=test_d,
                                worker_root=model_d, port=port)

    pst.pestpp_options["ies_num_reals"] = 300
    pst.write(os.path.join(template_d, "pest.pst"))

    test_d = os.path.join(model_d, "master_mean_bug")
    if os.path.exists(test_d):
       shutil.rmtree(test_d)
    pyemu.os_utils.start_workers(template_d, exe_path, "pest.pst", num_workers=20, master_dir=test_d,
                                worker_root=model_d, port=port)





    

def freyberg_stacked_pe_invest():
    import flopy
    model_d = "ies_freyberg"
    test_d = os.path.join(model_d, "master_stacked_pe")
    template_d = os.path.join(model_d, "template")
    # print("loading pst")
    pst = pyemu.Pst(os.path.join(template_d, "pest.pst"))
    pst.pestpp_options = {"ies_num_reals": 100}
    pst.pestpp_options["ies_save_lambda_en"] = True
    pst.pestpp_options["save_dense"] = True
    pst.pestpp_options["ies_no_noise"] = True
    pst.pestpp_options["ies_initial_lambda"] = 1000.0
    pst.control_data.noptmax = 2
    pst.observation_data.loc[:,"weight"] = 1.0
    pst.write(os.path.join(template_d, "pest.pst"))
    if os.path.exists(test_d):
       shutil.rmtree(test_d)
    pyemu.os_utils.start_workers(template_d, exe_path, "pest.pst", num_workers=10, master_dir=test_d,
                                worker_root=model_d, port=port)


    pelam_files = [f for f in os.listdir(test_d) if "lambda" in f and "par" in f and f.endswith(".bin")]
    assert len(pelam_files) > 0
    pes = [pyemu.ParameterEnsemble.from_binary(pst=None,filename=os.path.join(test_d,f))._df for f in pelam_files]
    pelam = pd.concat(pes)
    pelam.index = ["lambreal{0}".format(real) for real in range(pelam.index.shape[0])]
    
    pst = pyemu.Pst(os.path.join(test_d,"pest.pst"),result_dir=test_d)
    phidf1 = pst.ies.phiactual
    pe = pst.ies.paren
    pe = pe.loc[pe.index.get_level_values(0)<pst.control_data.noptmax,:]
    rnames = ["real{0}-iter{1}".format(rname,itr) for rname,itr in zip(pe.index.get_level_values(1),pe.index.get_level_values(0))]
    pe.index = rnames
    pe = pd.concat([pe,pelam])
    #print(pe)
   
    pyemu.ParameterEnsemble(df=pe,pst=pst).to_dense(os.path.join(template_d,"comb_pe.bin"))
    pst.pestpp_options["ies_par_en"] = "comb_pe.bin"
    pst.pestpp_options.pop("ies_num_reals")
    pst.control_data.noptmax = 1
    pst.write(os.path.join(template_d, "pest.pst"))
    test_d += "comb1"
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    pyemu.os_utils.start_workers(template_d, exe_path, "pest.pst", num_workers=30, master_dir=test_d,
                                 worker_root=model_d, port=port)
    pst = pyemu.Pst(os.path.join(test_d,"pest.pst"),result_dir=test_d)
    phidf2 = pst.ies.phiactual
    p1 = phidf1.loc[phidf1.iteration==phidf1.iteration.max(),:]
    rnames = ["real{0}-iter{1}".format(realname,p1.iteration.values[0]-1) for realname in p1.columns.values[6:]]
    rnames[-1] = "base"
    #print(rnames)
    #print(phidf2.columns.tolist())
    p1 = p1.iloc[:,6:]
    print(p1)
    print(rnames)
    p2 = phidf2.loc[phidf2.iteration==1,rnames]
    print(p2)
    
   
    fig,ax = plt.subplots(1,1)
    #print(p1)
    #print(p2)
    
    #ax.hist(p1.values,alpha=0.5,fc="0.5")
    #ax.hist(p2.values,alpha=0.5,fc="b")
    ax.scatter(p1.values,p2.values)
    #mn = min(ax.get_xlim()[0],ax.get_ylim()[0])
    #mx = max(ax.get_xlim()[1],ax.get_ylim()[1])
    #ax.semilogx()
    #ax.semilogy()
    #ax.set_xlim(mn,mx)
    #ax.set_ylim(mn,mx)
    ax.grid()
    plt.show()


def tenpar_iqr_bad_phi_sigma_test():
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_badphi")
    template_d = os.path.join(model_d, "test_template")

    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)

    if os.path.exists(test_d):
        shutil.rmtree(test_d)


    try:
        pst.pestpp_options.pop("ies_bad_phi")
    except KeyError:
        pass
    try:
        pst.pestpp_options.pop("ies_bad_phi_sigma")
    except KeyError:
        pass
    

    
    pst.control_data.noptmax = 1

    pst_name = "base.pst"
    pst.write(os.path.join(template_d,pst_name),version=2)
    pyemu.os_utils.start_workers(template_d, exe_path, pst_name, num_workers=8,
                                 master_dir=test_d, worker_root=model_d, port=port)
    
    df0 = pd.read_csv(os.path.join(test_d, "base.phi.actual.csv"),index_col=0)

    pst.pestpp_options["ies_bad_phi_sigma"] = -1.5
    pst_name = "iqr.pst"
    pst.write(os.path.join(template_d,pst_name),version=2)
    pyemu.os_utils.start_workers(template_d, exe_path, pst_name, num_workers=8,
                                 master_dir=test_d, worker_root=model_d, port=port)
    # check bad phi is correct 
    df1 = pd.read_csv(os.path.join(test_d, "iqr.phi.actual.csv"))
    vals = df0.iloc[0,5:].values
    vals.sort()
    n = len(vals)
    q1 = vals[int(0.25 * n)]
    q3 = vals[int(0.75 * n)]
    iqr = q3 - q1
    f = pst.pestpp_options["ies_bad_phi_sigma"]
    badphi = q3 + -f*iqr  
    assert df1.iloc[0,5:].max() < badphi


if __name__ == "__main__":
    tenpar_iqr_bad_phi_sigma_test()
    #tenpar_fixed_restart_test()
    #freyberg_stacked_pe_invest()
    #freyberg_mean_invest()
    #exit()
    #tenpar_adjust_weights_test()
    #tenpar_uniformdist_invest()
    #temp_plot()
    #tenpar_fixed_restart_test()
    #tenpar_xsec_aal_sigma_dist_test()
    #tenpar_consistency_test()
    #tenpar_mean_iter_sched_phifac_test()
    #tenpar_mean_iter_test_sched()
    #zdt1_weight_test()
    #tenpar_adjust_weights_test()
    #exit()

    #freyberg_stacked_pe_invest()

    #hosaki_invest(b_d="hosaki_seq")
    #plot_hosaki(b_d="hosaki_seq",n_iter_mean=[-1,-2,-4,999])
    #hosaki_invest(use_ineq=True,b_d="hosaki_ineq_seq",n_iter_mean=[-1,-2,-4,999])
    #plot_hosaki(b_d="hosaki_ineq_seq")
    #hosaki_invest(use_ineq=True,b_d="hosaki_ineq_best_seq",n_iter_mean=[-1,-2,-4,999])
    #plot_hosaki(b_d="hosaki_ineq_best")

    

    # poly_n_iter_mean_invest(b_d="poly_bps")
    # plot_poly(b_d="poly_bps")
    # poly_n_iter_mean_invest("poly_bps_ineq",use_ineq=True)
    # plot_poly(b_d="poly_bps_ineq")
    # poly_n_iter_mean_invest(b_d="poly_bps_minphi",n_iter_mean=-3)
    # plot_poly(b_d="poly_bps_minphi")
    #freyberg_center_on_test()
    #multimodal_test()
    #tenpar_adjust_weights_test()
    #tenpar_noise_invest()
    #shutil.copy2(os.path.join("..","exe","windows","x64","Debug","pestpp-ies.exe"),os.path.join("..","bin","win","pestpp-ies.exe"))
    #twopar_freyberg_resp_surface_invest()
    #plot_twopar_resp_results()
    #tenpar_high_phi_test()
    #zdt1_weight_test()
    #tenpar_mean_iter_test()
    #freyberg_center_on_test()
    #freyberg_rcov_test()
    #tenpar_upgrade_on_disk_test_weight_ensemble_test()
    # tenpar_base_run_test()
    #tenpar_adjust_weights_test()
    #tenpar_drop_violations_test()
    # tenpar_adjust_weights_test_by_real()
    # tenpar_base_par_file_test()
    # tenpar_xsec_autoadaloc_test()
    # tenpar_xsec_combined_autoadaloc_test()
    # tenpar_xsec_aal_sigma_dist_test()
    # tenpar_by_vars_test()
    # tenpar_xsec_aal_invest()
    # temp()
    # tenpar_localize_how_test()
    # clues_longnames_test()
    # freyberg_local_threads_test()
    # freyberg_aal_test()
    #tenpar_high_phi_test()
    #tenpar_adjust_weights_test()
    # freyberg_svd_draws_invest()
    #tenpar_xsec_aal_sigma_dist_test()
    # freyberg_combined_aal_test()
    # freyberg_aal_invest()
    # tenpar_high_phi_test()
    # freyberg_center_on_test()
    #freyberg_pdc_test()
    # freyberg_rcov_tet()
    # freyberg_center_on_test()
    #tenpar_align_test()
    #tenpar_align_test_2()
    # tenpar_covloc_test()
    #tenpar_upgrade_on_disk_test()
    #multimodal_test()
    #mm_invest()
    #plot_mm1_sweep_results()
    #plot_mm1_results(None, func="circle", show_info=True)
    #plot_mm1_results(None, func="circle", show_info=True,mm_d = os.path.join("mm1","master_mm_phifac_byreal_circle"))
    #plot_mm1_results(None, func="circle", show_info=True,mm_d = os.path.join("mm1","master_mm_centeron_phifac_circle"))

    
    
    #mm_invest()
    #zdt1_weight_test()
    #plot_zdt1_results(15)
    
    #tenpar_upgrade_on_disk_test_with_fixed()
    #tenpar_upgrade_on_disk_test_with_fixed2()
    #tenpar_high_phi_test()
    #freyberg_pdc_test()
    #tenpar_upgrade_on_disk_test_weight_ensemble_test()
    #invest()
    #tenpar_extra_binary_vars_test()
    #tenpar_covloc_test()
