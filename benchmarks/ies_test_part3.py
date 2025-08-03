# TODO: test variance and mean of draws, add chenoliver and test approx and full solution
import os
import shutil
import platform
import numpy as np
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
7) freyberg draw par n from full parcov supplied in file
8) freyberg full solution with empirical parcov - supplied par csv, obs csv and restart csv with fails, bad phi,MAP solution, prior scaling, lam mults 
9) synth restart and upgrade 1.1M par problem"""

ies_vars = ["ies_par_en", "ies_obs_en", "ies_restart_obs_en",
            "ies_bad_phi", "parcov_filename", "ies_num_reals",
            "ies_use_approx", "ies_use_prior_scaling", "ies_reg_factor",
            "ies_lambda_mults", "ies_initial_lambda","ies_include_base","ies_subset_size"]


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
    bin_path = os.path.join(bin_path,"linux")
elif "darwin" in platform.platform().lower() or "macos" in platform.platform().lower():
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

def clues_longnames_test():
    """clue long names tests"""
    model_d = "ies_clues"
    test_d = os.path.join(model_d, "master_longnames")
    template_d = os.path.join(model_d, "template")
    pst = pyemu.Pst(os.path.join(template_d, "pest.pst"))
    num_reals = 10
    if os.path.exists(test_d):
       shutil.rmtree(test_d)
    #shutil.copytree(template_d, test_d)
    #cov = pyemu.Cov.from_parameter_data(pst)
    pe = pyemu.ParameterEnsemble.from_uniform_draw(pst=pst, num_reals=num_reals)
    pe.to_csv(os.path.join(template_d,'par.csv'))
    pst.pestpp_options = {"ies_num_reals":num_reals}
    pst.control_data.noptmax = -1
    par = pst.parameter_data
    borked = par.loc[par.apply(lambda x: x.parlbnd >= x.parubnd, axis=1),"parnme"]
    pst.parameter_data.loc[borked,"parubnd"] += 1
    pst.write(os.path.join(template_d,"pest.pst"))
    #pyemu.os_utils.run("{0} {1}".format(exe_path, "pest.pst"), cwd=test_d)
    pyemu.os_utils.start_workers(template_d,exe_path,"pest.pst",5,
                                worker_root=model_d,master_dir=test_d,port=port)
    pdf = pd.read_csv(os.path.join(test_d,"pest.0.par.csv"),index_col=0)
    pdf.columns = pdf.columns.str.lower()
    dset = set(pdf.columns)
    pset = set(pst.par_names)
    d = dset.symmetric_difference(pset)
    assert len(d) == 0,d
    oset = set(pst.obs_names)
    odf = pd.read_csv(os.path.join(test_d,"pest.0.obs.csv"),index_col=0)
    odf.columns = odf.columns.str.lower()
    dset = set(odf.columns)
    d = dset.symmetric_difference(oset)
    assert len(d) == 0, d

    pst.pestpp_options["sweep_parameter_csv_file"] = "par.csv"
    pst.pestpp_options["ies_par_en"] = "par.csv"
    pst.write(os.path.join(template_d,"pest.pst"))

    pyemu.os_utils.start_workers(template_d, exe_path, "pest.pst", 5,
                                worker_root=model_d, master_dir=test_d,port=port)
    pdf = pd.read_csv(os.path.join(test_d, "pest.0.par.csv"), index_col=0)
    pdf.columns = pdf.columns.str.lower()
    dset = set(pdf.columns)
    pset = set(pst.par_names)
    d = dset.symmetric_difference(pset)
    assert len(d) == 0, d
    oset = set(pst.obs_names)
    odf = pd.read_csv(os.path.join(test_d, "pest.0.obs.csv"), index_col=0)
    odf.columns = odf.columns.str.lower()
    dset = set(odf.columns)
    d = dset.symmetric_difference(oset)
    assert len(d) == 0, d
    pyemu.os_utils.start_workers(template_d, exe_path.replace("-ies","-swp"), "pest.pst", 5,
                                worker_root=model_d, master_dir=test_d,port=port)

    odf = pd.read_csv(os.path.join(test_d, "sweep_out.csv"), index_col=0)
    odf.columns = odf.columns.str.lower()
    odf = odf.loc[:,pst.obs_names]
    dset = set(odf.columns)
    oset = set(pst.obs_names)
    d = oset - dset
    assert len(d) == 0, d

    pst.parameter_data.loc[pst.adj_par_names[10:],"partrans"] = "fixed"
    pst.control_data.noptmax = 1
    pst.write(os.path.join(template_d, "pest.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path.replace("-ies", "-glm"), "pest.pst", 5,
                                worker_root=model_d, master_dir=test_d,port=port)


def freyberg_dist_local_test():
    import flopy
    model_d = "ies_freyberg"
    test_d = os.path.join(model_d, "master_dist_local1")
    template_d = os.path.join(model_d, "template")
    pst = pyemu.Pst(os.path.join(template_d, "pest.pst"))
    par = pst.parameter_data
    pt = par.loc[par.partrans=="log","parnme"]
    def compare_par_summary():
        pe_base = pd.read_csv(os.path.join(test_d,"pest_local.0.par.csv"),index_col=0)
        pe_base.columns = pe_base.columns.str.lower()
        pe_base.loc[:,pt] = pe_base.loc[:,pt].apply(lambda x: np.log10(x))
        par_csvs = [f for f in os.listdir(test_d) if ".par.csv" in f and not "0.par.csv" in f]
        base_mean = pe_base.mean(axis=0)
        base_std = pe_base.std(axis=0)
        for par_csv in par_csvs:
            pe = pd.read_csv(os.path.join(test_d,par_csv),index_col=0)
            pe.columns = pe.columns.str.lower()
            pe.loc[:,pt] = pe.loc[:,pt].apply(lambda x: np.log10(x))
            for grp in pst.par_groups:
                pnames = par.loc[par.pargp==grp,"parnme"]
                mean = pe.loc[:,pnames].mean(axis=0)
                std = pe.loc[:,pnames].std(axis=0)
                d = (base_mean.loc[pnames] - mean)
                print(d.loc[d!=0.0])
                mean_diff = np.nanmean(((base_mean.loc[pnames] - mean) / base_mean.loc[pnames]))
                std_diff = np.nanmean(((base_std.loc[pnames] - std) / base_std.loc[pnames]))
                print(par_csv,grp,mean_diff, std_diff)

    

    m = flopy.modflow.Modflow.load("freyberg.nam",model_ws=template_d,load_only=[],check=False)
    if os.path.exists(test_d):
       shutil.rmtree(test_d)
    # print("loading pst")


    par = pst.parameter_data
    par.loc[:,"partrans"] = "fixed"
    par.loc[par.pargp=="hk","partrans"] = "log"

    par_adj = par.loc[pst.adj_par_names,:].copy()
    par_adj.loc[:,"i"] = par_adj.parnme.apply(lambda x: int(x.split('_')[1][1:]))
    par_adj.loc[:,"j"] = par_adj.parnme.apply(lambda x: int(x.split('_')[2][1:]))
    par_adj.loc[:,"x"] = par_adj.apply(lambda x: m.modelgrid.xcellcenters[x.i,x.j],axis=1)
    par_adj.loc[:,"y"] = par_adj.apply(lambda x: m.modelgrid.ycellcenters[x.i,x.j],axis=1)

    pst.observation_data.loc["flx_river_l_19700102","weight"] = 0.0
    obs_nz = pst.observation_data.loc[pst.nnz_obs_names,:].copy()
    obs_nz.loc[:,"i"] = obs_nz.obsnme.apply(lambda x: int(x[6:8]))
    obs_nz.loc[:,"j"] = obs_nz.obsnme.apply(lambda x: int(x[9:11]))
    obs_nz.loc[:,'x'] = obs_nz.apply(lambda x: m.modelgrid.xcellcenters[x.i,x.j],axis=1)
    obs_nz.loc[:,'y'] = obs_nz.apply(lambda x: m.modelgrid.ycellcenters[x.i,x.j],axis=1)

    dfs = []
    v = pyemu.geostats.ExpVario(contribution=1.0, a=1000)
    for name in pst.nnz_obs_names:
        x,y = obs_nz.loc[name,['x','y']].values
        #print(name,x,y)
        p = par_adj.copy()
        #p.loc[:,"dist"] = p.apply(lambda xx: np.sqrt((xx.x - x)**2 + (xx.y - y)**2),axis=1)
        #print(p.dist.max(),p.dist.min())
        cc = v.covariance_points(x,y,p.x,p.y)
        #print(cc.min(),cc.max())
        dfs.append(cc)

    df = pd.concat(dfs,axis=1)
    df.columns = pst.nnz_obs_names

    mat = pyemu.Matrix.from_dataframe(df.T)
    tol = 0.35

    mat.x[mat.x<tol] = 0.0
    mat.to_ascii(os.path.join(template_d,"localizer.mat"))
    df_tol = mat.to_dataframe()
    par_sum = df_tol.sum(axis=0)

    zero_cond_pars = list(par_sum.loc[par_sum==0.0].index)
    print(zero_cond_pars)

    cov = pyemu.Cov.from_parameter_data(pst)
    num_reals = 10
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst,cov=cov,num_reals=num_reals)
    pe.enforce()
    pe.to_csv(os.path.join(template_d,"par_local.csv"))

    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = num_reals
    pst.pestpp_options["ies_subset_size"] = num_reals
    pst.pestpp_options["ies_lambda_mults"] = 1.0
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["ies_include_base"] = False
    pst.pestpp_options["ies_par_en"] = "par_local.csv"
    pst.pestpp_options["ies_localizer"] = "localizer.mat"
    pst.pestpp_options["ies_localize_how"] = "par"
    pst.pestpp_options["ies_verbose_level"] = 1
    pst.pestpp_options["ies_num_threads"] = 3
    pst.pestpp_options["ies_subset_how"] = "random"
    pst.pestpp_options["ies_accept_phi_fac"] = 1000.0
    pst.pestpp_options["overdue_giveup_fac"] = 1000.0
    pst.control_data.noptmax = 2
    pst.write(os.path.join(template_d, "pest_local.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_local.pst", num_workers=20, master_dir=test_d,
                               worker_root=model_d,port=port)

    par_df_org = pd.read_csv(os.path.join(test_d, "pest_local.0.par.csv"), index_col=0)
    #par_df_org.index = pe.index
    par_df_org.columns = par_df_org.columns.str.lower()
    for i in range(0,pst.control_data.noptmax):
        f = os.path.join(test_d, "pest_local.{0}.par.csv".format(i+1))
        if not os.path.exists(f):
            continue
        # try:
        #     par_df = pd.read_csv(f,index_col=0)
        # except:
        #     break
        # par_df.index = pe.index
        par_df = pd.read_csv(f, index_col=0)
        par_df.index = pe.index
        par_df.columns = par_df.columns.str.lower()
        diff = par_df_org - par_df
        print(f)
        print(diff.loc[:,zero_cond_pars].max(axis=1))
        assert diff.loc[:,zero_cond_pars].max().max() < 1.0e-6, diff.loc[:,zero_cond_pars].max()
    compare_par_summary()


def freyberg_dist_local_invest():
    import flopy
    model_d = "ies_freyberg"
    test_d = os.path.join(model_d, "master_dist_local2")
    template_d = os.path.join(model_d, "template")
    m = flopy.modflow.Modflow.load("freyberg.nam", model_ws=template_d, load_only=[], check=False)
    #if os.path.exists(test_d):
    #   shutil.rmtree(test_d)
    # print("loading pst")
    pst = pyemu.Pst(os.path.join(template_d, "pest.pst"))

    obs = pst.observation_data
    pobs = obs.loc[obs.obgnme=="pothead","obsnme"]
    obs.loc[pobs[::2],"weight"] = 2.0
    obs.loc[pobs[::2],"obsval"] += np.random.normal(0.0,4.0,pobs[::2].shape[0])

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
    obs_nz.loc[:, "i"] = obs_nz.obsnme.apply(lambda x: int(x[6:8])-1)
    obs_nz.loc[:, "j"] = obs_nz.obsnme.apply(lambda x: int(x[9:11])-1)
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
    tol = 0.45
    mat.x[mat.x < tol] = 0.0
    #mat.x[mat.x > 0.0] = 1.0
    mat.to_ascii(os.path.join(template_d, "localizer.mat"))
    df_tol = mat.to_dataframe()
    par_sum = df_tol.sum(axis=0)
    par_sum = par_sum / par_sum.max()
    print(par_sum.max(0))

    zero_cond_pars = list(par_sum.loc[par_sum == 0.0].index)
    print(zero_cond_pars)

    arr = np.zeros((m.nrow, m.ncol))
    arr[par_adj.i, par_adj.j] = par_sum.values
    arr[arr == 0] = np.NaN
    # plt.imshow(arr)
    # plt.show()
    # return

    cov = pyemu.Cov.from_parameter_data(pst)
    num_reals = 8
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst, cov=cov, num_reals=num_reals)
    pe.enforce()
    pe.to_csv(os.path.join(template_d, "par_local.csv"))

    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = num_reals
    pst.pestpp_options["ies_subset_size"] = 3
    pst.pestpp_options["ies_lambda_mults"] = [1.0]
    pst.pestpp_options["lambda_scale_fac"] = [1.0]
    pst.pestpp_options["ies_include_base"] = True
    pst.pestpp_options["ies_par_en"] = "par_local.csv"
    pst.pestpp_options["ies_localizer"] = "localizer.mat"
    pst.pestpp_options["ies_verbose_level"] = 1
    pst.pestpp_options["ies_subset_how"] = "random"
    pst.pestpp_options["ies_accept_phi_fac"] = 1000.0
    pst.control_data.nphinored = 20
    #pst.pestpp_options["ies_initial_lambda"] = 0.1
    pst.control_data.noptmax = 20
    pst.write(os.path.join(template_d, "pest_local.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_local.pst", num_workers=20, master_dir=test_d,
                                worker_root=model_d, port=port)

    par_df_org = pd.read_csv(os.path.join(test_d, "pest_local.0.par.csv"), index_col=0)
    # par_df_org.index = pe.index
    par_df_org.columns = par_df_org.columns.str.lower()
    fig = plt.figure(figsize=(10 * pst.control_data.noptmax, 10))
    # axes = [plt.subplot(2,pst.control_data.noptmax+1,i+1) for i in range(pst.control_data.noptmax+1)]
    # axes[0].imshow(arr)
    ax = plt.subplot2grid((1, pst.control_data.noptmax + 1), (0, 0))
    cb = ax.imshow(arr)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_title("localizer function")
    # plt.colorbar(cb)

    diff_arrs = []
    for i in range(0, pst.control_data.noptmax):
        try:
            par_df = pd.read_csv(os.path.join(test_d, "pest_local.{0}.par.csv".format(i + 1)),
                                 index_col=0)
        except:
            break
        # par_df.index = pe.index
        par_df.columns = par_df.columns.str.lower()
        diff = par_df_org - par_df
        # print(diff.loc[:,zero_cond_pars].max(axis=1))
        #assert diff.loc[:, zero_cond_pars].max().max() < 1.0e-6, diff.loc[:, zero_cond_pars].max().max()
        dsum = diff.apply(lambda x: np.abs(x)).sum().loc[par_adj.index]
        d_arr = np.zeros((m.nrow, m.ncol))
        d_arr[par_adj.i, par_adj.j] = dsum
        d_arr[d_arr == 0.0] = np.NaN
        diff_arrs.append(d_arr)

        # #axes[i+1].imshow(d_arr)
        # ax = plt.subplot2grid((2, pst.control_data.noptmax + 2), (0, i+1))
        # ax.imshow(d_arr)
        # s_arr = np.zeros_like(d_arr)
        # s_arr[par_adj.i,par_adj.j] = par_df.loc[:,par_adj.index].std()
        # print(i)
        # ax = plt.subplot2grid((2, pst.control_data.noptmax + 2), (1,i+1))
        # ax.imshow(s_arr)
    # mn,mx = 1.0e+10,-1.0e+10

    mx = max([np.nanmax(a) for a in diff_arrs])
    mn = min([np.nanmin(a) for a in diff_arrs])
    for i, d_arr in enumerate(diff_arrs):
        ax = plt.subplot2grid((1, pst.control_data.noptmax + 1), (0, i + 1))
        ax.imshow(d_arr, vmax=mx, vmin=mn)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_title("HK diff from initial en, iteration {0}".format(i + 1))
    plt.savefig(os.path.join(test_d,"hk_compare.pdf"))

    pst.pestpp_options.pop("ies_localizer")
    pst.write(os.path.join(template_d, "pest_base.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_base.pst", num_workers=20, master_dir=test_d+"_base",
                                worker_root=model_d, port=port)

    phi_loc = pd.read_csv(os.path.join(test_d,"pest_local.phi.meas.csv"),index_col=0)
    phi_all = pd.read_csv(os.path.join(test_d+"_base","pest_base.phi.meas.csv"),index_col=0)

    fig = plt.figure()
    ax = plt.subplot(111)
    phi_loc.loc[:,"mean"].plot(ax=ax,color='g', label="local")
    phi_all.loc[:,"mean"].plot(ax=ax,color='b', label="full")
    ax.set_ylabel("phi")
    ax.legend()
    plt.savefig(os.path.join(test_d,"phi.pdf"))


def tenpar_localize_how_test():
    """tenpar local 1"""
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_localize_how_test")
    template_d = os.path.join(model_d, "test_template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    # shutil.copytree(template_d, test_d)
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)

    # mat = pyemu.Matrix.from_names(pst.nnz_obs_names,pst.adj_par_names).to_dataframe()
    mat = pyemu.Matrix.from_names(["head"], ["k"]).to_dataframe()
    mat.loc[:, :] = 1.0
    # mat.iloc[0,:] = 1
    mat = pyemu.Matrix.from_dataframe(mat)
    mat.to_ascii(os.path.join(template_d, "localizer.mat"))

    cov = pyemu.Cov.from_parameter_data(pst)
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst, cov=cov, num_reals=10)
    pe.enforce()
    pe.to_csv(os.path.join(template_d, "par_local.csv"))

    oe = pyemu.ObservationEnsemble.from_gaussian_draw(pst, num_reals=10)
    oe.to_csv(os.path.join(template_d, "obs_local.csv"))

    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options["ies_localizer"] = "localizer.mat"
    pst.pestpp_options["ies_lambda_mults"] = 1.0
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["ies_subset_size"] = 11
    pst.pestpp_options["ies_par_en"] = "par_local.csv"
    pst.pestpp_options["ies_obs_en"] = "obs_local.csv"
    #pst.pestpp_options["ies_localize_how"] = "o"
    #pst.pestpp_options["ies_verbose_level"] = 3
    #pst.control_data.noptmax = 3

    # pst.pestpp_options["ies_verbose_level"] = 3
    #pst_name = os.path.join(template_d, "pest_local_o.pst")
    #pst.write(pst_name)
    #pyemu.os_utils.start_workers(template_d, exe_path, "pest_local_o.pst", num_workers=10,
    #                            master_dir=test_d+"_o", verbose=True, worker_root=model_d,
    #                            port=port)
    #phi_df1 = pd.read_csv(os.path.join(test_d+"_o", "pest_local_o.phi.meas.csv"))

    pst.pestpp_options["ies_localize_how"] = "p"
    pst.write(os.path.join(template_d, "pest_local_p.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_local_p.pst", num_workers=10,
                                master_dir=test_d + "_p", verbose=True, worker_root=model_d,
                                port=port)
    phi_df2 = pd.read_csv(os.path.join(test_d + "_p", "pest_local_p.phi.meas.csv"))
    #diff = phi_df1 - phi_df2
    #print(diff.max().max())

    # plt.plot(phi_df1.total_runs, phi_df1.loc[:, "mean"], label="local")
    # plt.plot(phi_df2.total_runs, phi_df2.loc[:, "mean"], label="full")
    # plt.legend()
    # plt.savefig(os.path.join(test_d+"_p", "local_test.pdf"))
    #assert diff.max().max() == 0

    mat = pyemu.Matrix.from_names(pst.nnz_obs_names,pst.adj_par_names).to_dataframe()
    mat.loc[:, :] = 1.0
    # mat.iloc[0,:] = 1
    mat = pyemu.Matrix.from_dataframe(mat)
    mat.to_ascii(os.path.join(template_d, "localizer.mat"))

    pst.pestpp_options["ies_localize_how"] = "p"
    pst.control_data.noptmax = 2
    pst.write(os.path.join(template_d, "pest_local_p.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_local_p.pst", num_workers=10,
                                master_dir=test_d + "_p", verbose=True, worker_root=model_d,
                                port=port)
    mat = pyemu.Matrix.from_names(pst.nnz_obs_names,pst.adj_par_names).to_dataframe()
    mat.loc[:, :] = 1.0
    mat.iloc[0,1:] = 0
    mat.iloc[1,:] = np.arange(mat.shape[1]) / mat.shape[1]
    mat.iloc[:2,:] = 1
    mat.iloc[:,-1] = 0
    print(mat)
    mat = pyemu.Matrix.from_dataframe(mat)
    mat.to_ascii(os.path.join(template_d, "localizer.mat"))
    pst.control_data.noptmax = 2
    par = pst.parameter_data

    pst.write(os.path.join(template_d, "pest_local_p.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_local_p.pst", num_workers=10,
                                master_dir=test_d + "_p", verbose=True, worker_root=model_d,
                                port=port)
    
    par.loc[pst.adj_par_names[:2],"pargp"] = "test_group"
    adj_groups = par.loc[pst.adj_par_names,"pargp"].unique()
    mat = pyemu.Matrix.from_names(pst.nnz_obs_names,adj_groups).to_dataframe()
    mat.loc[:, :] = 0.1
    mat.iloc[:,0] = 0.0
    mat = pyemu.Matrix.from_dataframe(mat)
    mat.to_ascii(os.path.join(template_d, "localizer.mat"))
    pst.write(os.path.join(template_d, "pest_local_p.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_local_p.pst", num_workers=10,
                                master_dir=test_d + "_p", verbose=True, worker_root=model_d,
                                port=port)

    #pst.pestpp_options["ies_localize_how"] = "o"
    #pst.control_data.noptmax = 2
    #pst.write(os.path.join(template_d, "pest_local_p.pst"))
    #pyemu.os_utils.start_workers(template_d, exe_path, "pest_local_p.pst", num_workers=10,
    #                            master_dir=test_d + "_p", verbose=True, worker_root=model_d,
    #                            port=port)

def freyberg_local_threads_test():
    import flopy
    model_d = "ies_freyberg"
    test_d = os.path.join(model_d, "master_local_threads")
    template_d = os.path.join(model_d, "template")
    m = flopy.modflow.Modflow.load("freyberg.nam",model_ws=template_d,load_only=[],check=False)
    if os.path.exists(test_d):
       shutil.rmtree(test_d)
    # print("loading pst")
    pst = pyemu.Pst(os.path.join(template_d, "pest.pst"))

    par = pst.parameter_data
    par.loc[:,"partrans"] = "fixed"
    par.loc[par.pargp=="hk","partrans"] = "log"
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst,
        cov=pyemu.Cov.from_parameter_data(pst),num_reals=num_reals)
    pe.to_csv(os.path.join(template_d,"pars1.csv"))
    
    pyemu.ObservationEnsemble.from_gaussian_draw(pst=pst,num_reals=num_reals).to_csv(os.path.join(template_d,"obs1.csv"))

    par_adj = par.loc[pst.adj_par_names,:].copy()
    par_adj.loc[:,"i"] = par_adj.parnme.apply(lambda x: int(x.split('_')[1][1:]))
    par_adj.loc[:,"j"] = par_adj.parnme.apply(lambda x: int(x.split('_')[2][1:]))
    par_adj.loc[:,"x"] = par_adj.apply(lambda x: m.modelgrid.xcellcenters[x.i,x.j],axis=1)
    par_adj.loc[:,"y"] = par_adj.apply(lambda x: m.modelgrid.ycellcenters[x.i,x.j],axis=1)

    pst.observation_data.loc["flx_river_l_19700102","weight"] = 0.0
    obs_nz = pst.observation_data.loc[pst.nnz_obs_names,:].copy()
    obs_nz.loc[:,"i"] = obs_nz.obsnme.apply(lambda x: int(x[6:8]))
    obs_nz.loc[:,"j"] = obs_nz.obsnme.apply(lambda x: int(x[9:11]))
    obs_nz.loc[:,'x'] = obs_nz.apply(lambda x: m.modelgrid.xcellcenters[x.i,x.j],axis=1)
    obs_nz.loc[:,'y'] = obs_nz.apply(lambda x: m.modelgrid.ycellcenters[x.i,x.j],axis=1)

    dfs = []
    v = pyemu.geostats.ExpVario(contribution=1.0, a=1000)
    for name in pst.nnz_obs_names:
        x,y = obs_nz.loc[name,['x','y']].values
        #print(name,x,y)
        p = par_adj.copy()
        #p.loc[:,"dist"] = p.apply(lambda xx: np.sqrt((xx.x - x)**2 + (xx.y - y)**2),axis=1)
        #print(p.dist.max(),p.dist.min())
        cc = v.covariance_points(x,y,p.x,p.y)
        #print(cc.min(),cc.max())
        dfs.append(cc)

    df = pd.concat(dfs,axis=1)
    df.columns = pst.nnz_obs_names

    mat = pyemu.Matrix.from_dataframe(df.T)
    tol = 0.35

    mat.x[mat.x<tol] = 0.0
    mat.to_ascii(os.path.join(template_d,"localizer.mat"))
    df_tol = mat.to_dataframe()
    par_sum = df_tol.sum(axis=0)

    zero_cond_pars = list(par_sum.loc[par_sum==0.0].index)
    print(zero_cond_pars)

    #cov = pyemu.Cov.from_parameter_data(pst)
    #num_reals = 10
    #pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst,cov=cov,num_reals=num_reals)
    #pe.enforce()
    #pe.to_csv(os.path.join(template_d,"par_local.csv"))

    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = num_reals
    pst.pestpp_options["ies_subset_size"] = num_reals
    pst.pestpp_options["ies_lambda_mults"] = [0.1,1.0,5.0]
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    #pst.pestpp_options["ies_include_base"] = False
    #pst.pestpp_options["ies_par_en"] = "par_local.csv"
    pst.pestpp_options["ies_use_approx"] = False
    pst.pestpp_options["ies_use_prior_scaling"] = True
    pst.pestpp_options["ies_localizer"] = "localizer.mat"
    pst.pestpp_options["ies_localize_how"] = "par"
    pst.pestpp_options["ies_verbose_level"] = 3
    pst.pestpp_options["ies_save_lambda_en"] = True
    pst.pestpp_options["ies_subset_how"] = "random"
    pst.pestpp_options["ies_accept_phi_fac"] = 1000.0
    pst.pestpp_options["overdue_giveup_fac"] = 10.0
    pst.pestpp_options["ies_par_en"] = "pars1.csv"
    pst.pestpp_options["ies_obs_en"] = "obs1.csv"

    pst.control_data.noptmax = 3
    pst.write(os.path.join(template_d, "pest_local.pst"))
    d = test_d+"_base"
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_local.pst", num_workers=20, master_dir=d,
                               worker_root=model_d,port=port, verbose=True)

    def get_results(dd):

        phi_df = pd.read_csv(os.path.join(dd,"pest_local.phi.actual.csv"),index_col=0)
        upgrade_dfs = {}
        for f in os.listdir(dd):
            if f.endswith(".scale.par.csv"):
                df = pd.read_csv(os.path.join(dd,f),index_col=0)
                upgrade_dfs[f] = df
        return phi_df,upgrade_dfs

    base_phi,base_dfs = get_results(d)

    pst.pestpp_options["ies_num_threads"] = 1

    pst.write(os.path.join(template_d, "pest_local.pst"))
    d = test_d + "_1thread"

    pyemu.os_utils.start_workers(template_d, exe_path, "pest_local.pst", num_workers=20, master_dir=d,
                                worker_root=model_d, port=port)
    # phi,dfs = get_results(d)
    # phi_diff = (base_phi - phi).apply(np.abs)
    # assert phi_diff.max().max() == 0.0, phi_diff.max()
    # for f,df in base_dfs.items():
    #     assert f in dfs,f
    #     diff = (df - dfs[f]).apply(np.abs)
    #     assert diff.max().max() == 0.0,diff.max()

    pst.pestpp_options["ies_num_threads"] = 10
    pst.write(os.path.join(template_d, "pest_local.pst"))
    d = test_d + "_10thread"
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_local.pst", num_workers=20, master_dir=d,
                                worker_root=model_d, port=port)
    
    phi_1thread, dfs_1thread = get_results(test_d+"_1thread")
    phi_10thread, dfs_10thread = get_results(test_d+"_10thread")
    for c in phi_1thread.columns:
        print(c,base_phi.loc[2,c],phi_1thread.loc[2,c],phi_10thread.loc[2,c])

    for d in [test_d + "_1thread",test_d+"_10thread"]:
        phi, dfs = get_results(d)
        phi_diff = (base_phi - phi).apply(lambda x: np.abs(x))
        assert phi_diff.max().max() == 0.0, phi_diff.max()
        for f, df in base_dfs.items():
            assert f in dfs, f
            diff = (df - dfs[f]).apply(lambda x: np.abs(x))
            assert diff.max().max() == 0.0, diff.max()


def tenpar_tied_test():
    """tenpar tied test"""
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_tied_test")
    template_d = os.path.join(model_d, "test_template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    # shutil.copytree(template_d, test_d)
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)
    print(pst.adj_par_names)
    tnames = ["k_03","k_04"]

    par = pst.parameter_data
    par.loc[tnames,"partrans"] = "tied"
    par.loc[tnames,"partied"] = "k_02"
    par.loc["k_02","parval1"] = 2.5

    
    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options["ies_lambda_mults"] = [1.0,2.0]
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    #pst.pestpp_options["ies_subset_size"] = 11
    pst.pestpp_options["ies_save_binary"] = True
    
    pst.control_data.noptmax = 0
    pst_name = os.path.join(template_d, "pest_tied.pst")
    pst.write(pst_name)
    pyemu.os_utils.run("{0} {1}".format(exe_path,"pest_tied.pst"),cwd=template_d)
    df = pyemu.pst_utils.read_parfile(os.path.join(template_d,"pest_tied.base.par"))
    assert df.loc[tnames[0],"parval1"] == df.loc["k_02","parval1"]
    assert df.loc[tnames[1], "parval1"] == df.loc["k_02", "parval1"]

    #pst.pestpp_options["ies_verbose_level"] = 3
    pst.control_data.noptmax = 1

    # pst.pestpp_options["ies_verbose_level"] = 3
    pst_name = os.path.join(template_d, "pest_tied.pst")
    pst.write(pst_name)
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_tied.pst", num_workers=5,
                                master_dir=test_d, verbose=True, worker_root=model_d,
                                port=port)
    
    pe = pyemu.ParameterEnsemble.from_binary(filename=os.path.join(test_d,"pest_tied.0.par.jcb"),pst=pst)
    assert pe.shape[1] == pst.npar,"{0} vs {1}".format(pe.shape[1],pst.npar)
    assert list(pe.columns) == pst.par_names
    d = pe._df.loc[:,tnames[0]] - pe._df.loc[:,"k_02"].values
    print(d)
    assert np.abs(d).max() < 0.0001
    d = pe._df.loc[:,tnames[1]] - pe._df.loc[:,"k_02"].values
    print(d)
    assert np.abs(d).max() < 0.0001

    pe = pyemu.ParameterEnsemble.from_binary(filename=os.path.join(test_d,"pest_tied.1.par.jcb"),pst=pst)
    assert pe.shape[1] == pst.npar,"{0} vs {1}".format(pe.shape[1],pst.npar)
    assert list(pe.columns) == pst.par_names
    d = pe._df.loc[:,tnames[0]] - pe._df.loc[:,"k_02"].values
    print(d)
    assert np.abs(d).max() < 0.0001
    d = pe._df.loc[:,tnames[1]] - pe._df.loc[:,"k_02"].values
    print(d)
    assert np.abs(d).max() < 0.0001

    pst.pestpp_options["ies_save_binary"] = False
    pst.control_data.noptmax = 1
    # pst.pestpp_options["ies_verbose_level"] = 3
    pst_name = os.path.join(template_d, "pest_tied.pst")
    pst.write(pst_name)
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_tied.pst", num_workers=5,
                                master_dir=test_d, verbose=True, worker_root=model_d,
                                port=port)
    pe = pyemu.ParameterEnsemble.from_csv(filename=os.path.join(test_d, "pest_tied.0.par.csv"), pst=pst)
    assert pe.shape[1] == pst.npar, "{0} vs {1}".format(pe.shape[1], pst.npar)
    assert list(pe.columns) == pst.par_names
    d = pe._df.loc[:, tnames[0]] - pe._df.loc[:, "k_02"].values
    print(d)
    assert np.abs(d).max() < 0.0001
    d = pe._df.loc[:, tnames[1]] - pe._df.loc[:, "k_02"].values
    print(d)
    assert np.abs(d).max() < 0.0001

    pe = pyemu.ParameterEnsemble.from_csv(filename=os.path.join(test_d, "pest_tied.1.par.csv"), pst=pst)
    assert pe.shape[1] == pst.npar, "{0} vs {1}".format(pe.shape[1], pst.npar)
    assert list(pe.columns) == pst.par_names
    d = pe._df.loc[:, tnames[0]] - pe._df.loc[:, "k_02"].values
    print(d)
    assert np.abs(d).max() < 0.0001
    d = pe._df.loc[:, tnames[1]] - pe._df.loc[:, "k_02"].values
    print(d)
    assert np.abs(d).max() < 0.0001
    





def tenpar_by_vars_test():
    """tenpar vars test"""
    model_d = "ies_10par_xsec"
    test_d = os.path.join(model_d, "master_vars_test")
    template_d = os.path.join(model_d, "test_template")
    if not os.path.exists(template_d):
        raise Exception("template_d {0} not found".format(template_d))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    # shutil.copytree(template_d, test_d)
    pst_name = os.path.join(template_d, "pest.pst")
    pst = pyemu.Pst(pst_name)
    tnames = ["k_02","k_03"]
    fnames = ["k_04"]
    par = pst.parameter_data
    par.loc[tnames,"partrans"] = "tied"
    par.loc[tnames,"partied"] = "k_01"
    par.loc[fnames,"partrans"] = "fixed"
    par.loc["k_01","partrans"] = "log"
    
    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = 10
    pst.pestpp_options["ies_lambda_mults"] = [0.1,1.0,2.0]
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    #pst.pestpp_options["ies_subset_size"] = 11
    pst.pestpp_options["ies_save_binary"] = False
    pst.pestpp_options["ies_csv_by_reals"] = True
    pst.pestpp_options["ensemble_output_precision"] = 12
    #pst.pestpp_options["ies_verbose_level"] = 3
    pst.control_data.noptmax = 2

    # pst.pestpp_options["ies_verbose_level"] = 3
    pst_name = os.path.join(template_d, "pest_vars.pst")
    pst.write(pst_name)
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_vars.pst", num_workers=10,
                                master_dir=test_d, verbose=True, worker_root=model_d,
                                port=port)
    
    pe1 = pd.read_csv(os.path.join(test_d,"pest_vars.2.par.csv"),index_col=0)
    
    pst.pestpp_options["ies_csv_by_reals"] = False
    pst.write(pst_name)
    
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_vars.pst", num_workers=10,
                                master_dir=test_d, verbose=True, worker_root=model_d,
                                port=port)
    
    pe2 = pd.read_csv(os.path.join(test_d,"pest_vars.2.par.csv"),index_col=0).T
    pe2.index = pe1.index
    diff = pe2 - pe1
    print(diff.apply(lambda x: np.abs(x)).max())
    assert diff.apply(lambda x: np.abs(x)).max().max() == 0.0

    shutil.copy2(os.path.join(test_d,"pest_vars.0.par.csv"),os.path.join(template_d,"restart_by_vars.par.csv"))
    shutil.copy2(os.path.join(test_d, "pest_vars.obs+noise.csv"), os.path.join(template_d, "restart_by_vars.obs.csv"))
    pst.pestpp_options['ies_par_en'] = "restart_by_vars.par.csv"
    pst.pestpp_options['ies_obs_en'] = "restart_by_vars.obs.csv"
    df = pd.read_csv(os.path.join(template_d,"restart_by_vars.par.csv"),index_col=0)
    fval = 2.50000012345
    df.loc[fnames,:] = fval
    df.to_csv(os.path.join(template_d,"restart_by_vars.par.csv"))
    pst.write(pst_name)
    test_d = os.path.join(model_d, "master_vars_test1")
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_vars.pst", num_workers=10,
                                master_dir=test_d, verbose=True, worker_root=model_d,
                                port=port)
    pe2 = pd.read_csv(os.path.join(test_d, "pest_vars.2.par.csv"), index_col=0).T
    pe2.index = pe1.index
    diff = pe2 - pe1
    print(diff.apply(lambda x: np.abs(x)).max())
    assert diff.apply(lambda x: np.abs(x)).max().max() < 1.0e-3
    with pd.option_context('display.precision', 12):
        print(pe2.loc[:,fnames])
        d = 1000000. * np.abs(pe2.loc[:,fnames] - fval) / fval
        print(d)
    assert d.max().max() == 0.0,d.max()

def tenpar_xsec_autoadaloc_test():
    """testing the CC calculations made by autoadaloc"""
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
    pst.pestpp_options["ies_num_reals"] = 50
    pst.control_data.noptmax = -1
    pst.write(os.path.join(template_d,"pest_aal.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_aal.pst", num_workers=10,
                               master_dir=test_d, verbose=True, worker_root=model_d,
                               port=port)
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

    for f in ["pest_aal.0.par.csv","pest_aal.obs+noise.csv","pest_aal.0.obs.csv"]:
        shutil.copy2(os.path.join(test_d,f),os.path.join(template_d,f))



    def compare(case):
        df_aal = pd.read_csv(os.path.join(test_d,"pest_aal_restart.1.autoadaloc.csv"))
        df_aal.loc[:,"parnme"] = df_aal.parnme.str.lower()
        df_aal.loc[:, "obsnme"] = df_aal.obsnme.str.lower()
        df_aal.loc[:,"names"] = df_aal.apply(lambda x: (x.obsnme,x.parnme),axis=1)

        #pname = pst.adj_par_names[0]
        #oname = pst.nnz_obs_names[0]
        pe.transform()
        for pname in pst.adj_par_names:
            for oname in pst.nnz_obs_names:
                df = pd.DataFrame(data={pname:pe.loc[:,pname].values,oname:oe.loc[:,oname].values})
                df_shift = pd.DataFrame(data={pname: pe.loc[:, pname].values, oname: np.roll(oe.loc[:, oname].values,1)})
                cc = df.corr().loc[oname,pname]
                bg_cc_0 = df_shift.corr().loc[oname,pname]
                aal_cc,aal_bg = df_aal.loc[df_aal.names==(oname,pname),["correlation_coeff","0"]].values[0,:]
                print(case,pname,oname,cc,aal_cc,bg_cc_0,aal_bg)
                assert np.abs(cc - aal_cc) < 1.0e-4
                #assert np.abs(bg_cc_0 - aal_bg) < 1.0e-4
        mat = pyemu.Matrix.from_ascii(os.path.join(test_d,"pest_aal_restart.1.autoadaloc.tCC.mat"))

    pst.write(os.path.join(template_d, "pest_aal_restart.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_aal_restart.pst", num_workers=10,
                                master_dir=test_d, verbose=True, worker_root=model_d,
                                port=port)
    compare(case="single thread")

    pst.pestpp_options["ies_num_threads"] = 10
    pst.write(os.path.join(template_d, "pest_aal_restart.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_aal_restart.pst", num_workers=10,
                                master_dir=test_d, verbose=True, worker_root=model_d,
                                port=port)
    compare(case="10 thread")

    pst.pestpp_options["ies_debug_fail_subset"] = True
    pst.pestpp_options["ies_debug_fail_remainder"] = True
    pst.pestpp_options["ies_debug_bad_phi"] = True
    pst.write(os.path.join(template_d, "pest_aal_restart.pst"))
    pyemu.os_utils.start_workers(template_d, exe_path, "pest_aal_restart.pst", num_workers=10,
                                master_dir=test_d, verbose=True, worker_root=model_d,
                                port=port)




if __name__ == "__main__":
    #tenpar_base_run_test()
    #tenpar_xsec_autoadaloc_test()
    #tenpar_xsec_combined_autoadaloc_test()
    #tenpar_xsec_aal_sigma_dist_test()
    #shutil.copy2(os.path.join("..", "exe", "windows", "x64", "Debug", "pestpp-ies.exe"),
    #             os.path.join("..", "bin", "win","pestpp-ies.exe"))
    tenpar_by_vars_test()
    #tenpar_xsec_aal_invest()
    #temp()
    #tenpar_localize_how_test()
    #clues_longnames_test()
    #freyberg_local_threads_test()
    #tenpar_tied_test()
    #freyberg_dist_local_test()
