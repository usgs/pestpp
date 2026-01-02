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



def basic_sqp_test():
    model_d = "mf6_freyberg"
    local=True
    if "linux" in platform.platform().lower() and "10par" in model_d:
        #print("travis_prep")
        #prep_for_travis(model_d)
        local=False
    
    t_d = os.path.join(model_d,"template")

    m_d = os.path.join(model_d,"master_sqp1")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pst = pyemu.Pst(os.path.join(t_d,"freyberg6_run_opt.pst"))
    pst.pestpp_options.pop("ies_par_en")
    pst.control_data.noptmax = 0
    pst.write(os.path.join(t_d,"freyberg6_run_sqp.pst"))
    pyemu.os_utils.run("{0} freyberg6_run_sqp.pst".format(exe_path),cwd=t_d)

    assert os.path.exists(os.path.join(t_d,"freyberg6_run_sqp.base.par"))
    assert os.path.exists(os.path.join(t_d,"freyberg6_run_sqp.base.rei"))

    pst.pestpp_options["sqp_num_reals"] = 10

    pst.control_data.noptmax = -1
    pst.write(os.path.join(t_d,"freyberg6_run_sqp.pst"))
    pyemu.os_utils.start_workers(t_d, exe_path, "freyberg6_run_sqp.pst", 
                                 num_workers=5, master_dir=m_d,worker_root=model_d,
                                 port=port)

    assert os.path.exists(os.path.join(m_d,"freyberg6_run_sqp.0.par.csv"))
    df = pd.read_csv(os.path.join(m_d,"freyberg6_run_sqp.0.par.csv"),index_col=0)
    assert df.shape == (pst.pestpp_options["sqp_num_reals"],pst.npar),str(df.shape)
    assert os.path.exists(os.path.join(m_d,"freyberg6_run_sqp.0.obs.csv"))
    df = pd.read_csv(os.path.join(m_d,"freyberg6_run_sqp.0.obs.csv"),index_col=0)
    assert df.shape == (pst.pestpp_options["sqp_num_reals"],pst.nobs),str(df.shape)


def rosenbrock_setup(version,initial_decvars=1.6,constraints=False,constraint_exp="one_linear"):
    model_d = "rosenbrock"
    t_d = os.path.join(model_d, "template")

    if version == "2par":
        if constraints is True:
            if "two_linear" in constraint_exp:
                raise Exception
                w_d = os.path.join(os.path.join("rosenbrock", "2par_two_linear_constraints"))
            elif "one_linear" in constraint_exp:
                w_d = os.path.join(os.path.join("rosenbrock", "2par_one_linear_constraint"))
        else:
            raise Exception
            w_d = os.path.join(os.path.join("rosenbrock", "2par_unconstrained"))
    elif version == "high_dim":
        if constraints is not None:
            raise Exception
        else:
            raise Exception
            w_d = os.path.join(os.path.join("rosenbrock", "high_dim"))

    if os.path.exists(w_d):
        shutil.rmtree(w_d)
    shutil.copytree(os.path.join(t_d), w_d)
    os.chdir(w_d)

    in_file = os.path.join("par.dat")
    tpl_file = in_file+".tpl"
    out_file = [os.path.join("obs.dat")]
    ins_file = [out_file[0]+".ins"]

    if constraints is True:
        out_file.append(os.path.join("constraints.dat"))
        ins_file.append(out_file[1] + ".ins")

    pst = pyemu.helpers.pst_from_io_files(tpl_file,in_file,ins_file,out_file)

    par = pst.parameter_data
    par.loc[:, "partrans"] = "none"
    if version == "2par" and type(initial_decvars) is not float:
        par.loc[par.parnme[0], "parval1"] = initial_decvars[0]
        par.loc[par.parnme[1], "parval1"] = initial_decvars[1]
    else:
        par.loc[:, "parval1"] = initial_decvars
    par.loc[:, "parubnd"] = 2.2
    par.loc[:, "parlbnd"] = -2.2
    # TODO: repeat with log transform
    par.loc[:, "parchglim"] = "relative"

    obs = pst.observation_data
    obs.loc["obs", "obsval"] = 0.0
    obs.loc["obs", "obgnme"] = "obj_fn"

    pst.pestpp_options["opt_obj_func"] = "obs"

    if constraints is True:
        if "two_linear" in constraint_exp:
            raise Exception
            obs.loc["constraint_0", "obgnme"] = "l_constraint"  #"g_constraint"  # inherit from pestpp_options
            obs.loc["constraint_0", "obsval"] = -1.5  #10.0  # inherit from pestpp_options
            obs.loc["constraint_1", "obgnme"] = "l_constraint"  # inherit from pestpp_options
            obs.loc["constraint_1", "obsval"] = -0.2  #2.0  # inherit from pestpp_options
        elif "one_linear" in constraint_exp:
            obs.loc["constraint", "obgnme"] = "l_constraint"  #"g_constraint"  # inherit sign from pestpp_options
            obs.loc["constraint", "obsval"] = 2.0  #10.0  # inherit from pestpp_options
    obs.loc[:, "weight"] = 1.0

    pst.control_data.noptmax = 0
    if version == "2par":
        if constraints is True:
            if "two_linear" in constraint_exp:
                raise Exception
                pst.model_command = ["python rosenbrock_2par_two_linear_constraints.py"]
                pst.write(os.path.join("rosenbrock_2par_two_linear_constraints.pst"))
            elif "one_linear" in constraint_exp:
                pst.model_command = ["python rosenbrock_2par_constrained.py"]
                pst.write(os.path.join("rosenbrock_2par_constrained.pst"))
        else:
            pst.model_command = ["python rosenbrock_2par.py"]
            pst.write(os.path.join("rosenbrock_2par.pst"))
    elif version == "high_dim":
        if constraints:
            raise Exception
        else:
            pst.model_command = ["python rosenbrock_high_dim.py"]
            pst.write(os.path.join("rosenbrock_high_dim.pst"))

    os.chdir(os.path.join("..", ".."))

def rosenbrock_multiple_update(version,nit=10,draw_mult=3e-5,en_size=20,finite_diff_grad=False,
                               constraints=None,biobj_weight=1.0,biobj_transf=True,
                               hess_self_scaling=True,hess_update=True,damped=True,
                               cma=False,derinc=0.001,alg="BFGS",memory=5,strong_Wolfe=False,
                               rank_one=False,learning_rate=0.5,
                               mu_prop=0.25,use_dist_mean_for_delta=False,mu_learning_prop=0.5,
                               working_set=None,constraint_exp="one_linear",
                               qp_solve_method="null_space",reduced_hessian=False): #filter_thresh=1e-2

    if version == "2par":
        if constraints is True:
            if "two_linear" in constraint_exp:
                raise Exception
                w_d = os.path.join(os.path.join("rosenbrock", "2par_two_linear_constraints"))
                ext = version + "_two_linear_constraints"
            elif "one_linear" in constraint_exp:
                w_d = os.path.join(os.path.join("rosenbrock", "2par_one_linear_constraint"))
                ext = version + "_constrained"
        else:
            raise Exception
            w_d = os.path.join(os.path.join("rosenbrock", "2par"))
            ext = version
    elif version == "high_dim":
        if constraints:
            raise Exception
        else:
            raise Exception
            w_d = os.path.join(os.path.join("rosenbrock", "high_dim"))
    os.chdir(w_d)

    [os.remove(x) for x in os.listdir() if x.endswith("obsensemble.0000.csv")]
    [os.remove(x) for x in os.listdir() if x.endswith("parensemble.0000.csv")]
    [os.remove(x) for x in os.listdir() if (x.startswith("rosenbrock_{}.pst.".format(ext))) and (x.endswith(".csv"))]
    [os.remove(x) for x in os.listdir() if (x.startswith("filter.") and "csv" in x)]
    [os.remove(x) for x in os.listdir() if ("_per_" in x) and ("_alpha_" in x)]
    [os.remove(x) for x in os.listdir() if x == "hess_progress.csv"]

    pcf = "rosenbrock_{}.pst".format(ext)
    pst = pyemu.Pst(pcf)

    # forward checks
    pst.control_data.noptmax = 0
    pst.write(os.path.join(pcf.split(".")[0] + "_run_sqp.pst"))
    pyemu.os_utils.run("{0} {1}".format(exe_path, pcf.split(".")[0] + "_run_sqp.pst"))

    # esqp = pyemu.EnsembleSQP(pst="rosenbrock_{}.pst".format(ext))#,num_slaves=10)
    # esqp.initialize(num_reals=en_size,draw_mult=draw_mult,constraints=constraints,finite_diff_grad=finite_diff_grad,
    #               working_set=working_set)

    # initialize with noptmax = -1, i.e. calc grad
    pst.control_data.noptmax = -1
    #pst.pestpp_options["sqp_use_ensemble_grad"] = True  # with ensembles # omitted as controlled by num_reals arg
    pst.pestpp_options["sqp_num_reals"] = 10
    pst.pestpp_options["par_sigma_range"] = 20
    pst.write(os.path.join(pcf.split(".")[0] + "_run_sqp.pst"))
    pyemu.os_utils.run("{0} {1}".format(exe_path, pcf.split(".")[0] + "_run_sqp.pst"))

    # and with finite differences
    #pst.pestpp_options["sqp_use_ensemble_grad"] = False
    #pst.write(os.path.join(pcf.split(".")[0] + "_run_sqp.pst"))
    #pyemu.os_utils.run("{0} {1}".format(exe_path, pcf.split(".")[0] + "_run_sqp.pst"))

    # do update
    pst.control_data.noptmax = nit
    # TODO: add ++ args in update below
    pst.write(os.path.join(pcf.split(".")[0] + "_run_sqp.pst"))
    pyemu.os_utils.run("{0} {1}".format(exe_path, pcf.split(".")[0] + "_run_sqp.pst"))
    # step_mults = list(np.logspace(-5, 0, 12))  #list(np.linspace(0.1,1.0,10))
    #for it in range(nit):
     #   esqp.update(step_mult=step_mults,constraints=constraints,biobj_weight=biobj_weight,
      #              hess_self_scaling=hess_self_scaling,hess_update=hess_update,damped=damped,
       #             finite_diff_grad=finite_diff_grad,derinc=derinc,alg=alg,memory=memory,strong_Wolfe=strong_Wolfe,
        #            cma=cma,rank_one=rank_one,learning_rate=learning_rate,mu_prop=mu_prop,
         #           use_dist_mean_for_delta=use_dist_mean_for_delta,mu_learning_prop=mu_learning_prop,
          #          qp_solve_method=qp_solve_method,reduced_hessian=reduced_hessian)

    os.chdir(os.path.join("..", ".."))


def rosenbrock_single_linear_constraint(nit):
    constraints, constraint_exp = True, "one_linear"
    if "one_linear" in constraint_exp:
        yy = -1.0  #-2.0  #0.5
        idv = [(2 - yy) / -2.25, yy]  #[(10 - yy) / 6, yy]  #[1.8, yy]
        working_set = ['constraint']
    elif "two_linear" in constraint_exp:
        yy = 1.0
        idv = [(4 - yy) / 1.5, yy]  #[1.75, 1.10]
        working_set = ['constraint_1']  #[]

    rosenbrock_setup(version="2par", constraints=constraints, initial_decvars=idv, constraint_exp=constraint_exp)
    rosenbrock_multiple_update(version="2par", constraints=constraints, finite_diff_grad=True, nit=nit,
                               working_set=working_set, hess_update=False, hess_self_scaling=False,
                               constraint_exp=constraint_exp, qp_solve_method="direct",
                               reduced_hessian=False, damped=False)  # biobj_weight=5.0,alg="LBFGS",damped=False)
    #filter_plot(problem="2par", constraints=True, log_phi=True)
    #plot_2par_rosen(finite_diff_grad=True, constraints=constraints, constraint_exp=constraint_exp,
     #               label="surf_20_no_H_new.pdf")

def dewater_basic_test():
    model_d = "dewater"
    local=True
    if "linux" in platform.platform().lower() and "10par" in model_d:
        #print("travis_prep")
        #prep_for_travis(model_d)
        local=False
    
    t_d = os.path.join(model_d,"template")

    case = "dewater_pest.base"
    pst = pyemu.Pst(os.path.join(t_d,case+".pst"))
    par = pst.parameter_data
    dv_pars = par.loc[par.pargp == "q", "parnme"].tolist()[:3]
    #pst.add_pi_equation(dv_pars,"eq2",1000,obs_group="less_than")
    pst.pestpp_options = {}
    pst.pestpp_options["opt_dec_var_groups"] = "q"
    pst.control_data.noptmax = 0
    case = "test"
    pst.write(os.path.join(t_d,case+".pst"))
    pyemu.os_utils.run("{0} {1}.pst".format(exe_path,case),cwd=t_d)


    assert os.path.exists(os.path.join(t_d,case+".base.par"))
    assert os.path.exists(os.path.join(t_d,case+".base.rei"))

    m_d = os.path.join(model_d,"master_basic")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)

    pst.control_data.noptmax = -1  
    pst.write(os.path.join(t_d,case+".pst"))
    #pyemu.os_utils.run("{0} {1}.pst".format(exe_path,case),cwd=t_d)
    m_d = os.path.join(model_d, "master1")
    pyemu.os_utils.start_workers(t_d, exe_path, case + ".pst", num_workers=10, worker_root=model_d,
                                 master_dir=m_d)

    assert os.path.exists(os.path.join(m_d,case+".base.par"))
    assert os.path.exists(os.path.join(m_d,case+".base.rei"))
    assert os.path.exists(os.path.join(m_d,case+".0.jcb"))

    shutil.copy(os.path.join(m_d,case+".0.jcb"),os.path.join(t_d,"restart.jcb"))
    pst.pestpp_options["base_jacobian"] = "restart.jcb"
    pst.write(os.path.join(t_d,case+".pst"))
    #pyemu.os_utils.run("{0} {1}.pst".format(exe_path,case),cwd=t_d)
    m_d = os.path.join(model_d, "master1")
    pyemu.os_utils.start_workers(t_d, exe_path, case + ".pst", num_workers=10, worker_root=model_d,
                                 master_dir=m_d)

    assert os.path.exists(os.path.join(m_d,case+".base.par"))
    assert os.path.exists(os.path.join(m_d,case+".base.rei"))
    assert os.path.exists(os.path.join(m_d,case+".0.jcb"))

    pst.pestpp_options["base_jacobian"] = "dewater_pest.full.jcb"
    pst.write(os.path.join(t_d,case+".pst"))
    #pyemu.os_utils.run("{0} {1}.pst".format(exe_path,case),cwd=t_d)
    m_d = os.path.join(model_d, "master1")
    pyemu.os_utils.start_workers(t_d, exe_path, case + ".pst", num_workers=10, worker_root=model_d,
                                 master_dir=m_d)

    assert os.path.exists(os.path.join(m_d,case+".base.par"))
    assert os.path.exists(os.path.join(m_d,case+".base.rei"))
    assert os.path.exists(os.path.join(m_d,case+".0.jcb"))

    shutil.copy(os.path.join(m_d,case+".0.jcb"),os.path.join(t_d,"restart.jcb"))
    pst.pestpp_options["base_jacobian"] = "restart.jcb"
    pst.control_data.noptmax = 3
    pst.write(os.path.join(t_d,case+".pst"))
    #pyemu.os_utils.run("{0} {1}.pst".format(exe_path,case),cwd=t_d)
    m_d = os.path.join(model_d, "master2")
    pyemu.os_utils.start_workers(t_d, exe_path, case + ".pst", num_workers=10, worker_root=model_d,
                                 master_dir=m_d)
    assert os.path.exists(os.path.join(m_d, case + ".base.par"))
    assert os.path.exists(os.path.join(m_d, case + ".base.rei"))
    assert os.path.exists(os.path.join(m_d, case + ".0.jcb"))
    assert os.path.exists(os.path.join(m_d, case + ".1.jcb"))
    assert os.path.exists(os.path.join(m_d, case + ".2.jcb"))
    
    pst.pestpp_options["sqp_num_reals"] = 50
    pst.control_data.noptmax = 3
    pst.write(os.path.join(t_d, case + ".pst"))
    # pyemu.os_utils.run("{0} {1}.pst".format(exe_path,case),cwd=t_d)
    m_d = os.path.join(model_d, "master2_enopt")
    pyemu.os_utils.start_workers(t_d, exe_path, case + ".pst", num_workers=20, worker_root=model_d,
                                 master_dir=m_d)



def dewater_slp_opt_test():
    model_d = "dewater"
    local = True
    if "linux" in platform.platform().lower() and "10par" in model_d:
        # print("travis_prep")
        # prep_for_travis(model_d)
        local = False

    t_d = os.path.join(model_d, "template")

    case = "dewater_pest.base"
    pst = pyemu.Pst(os.path.join(t_d, case + ".pst"))
    par = pst.parameter_data
    dv_pars = par.loc[par.pargp == "q", "parnme"].tolist()[:3]
    
    pst.add_pi_equation(par_names=dv_pars, pilbl="eq3", rhs=1000, obs_group="less_than")
    pst.pestpp_options = {}
    pst.pestpp_options["opt_dec_var_groups"] = "q"
    pst.control_data.noptmax = 1
    print(pst.prior_information)
    pst.write(os.path.join(t_d, "test_opt.pst"))
    pyemu.os_utils.run("{0} {1}.pst".format(exe_path.replace("-sqp","-opt"), "test_opt.pst"), cwd=t_d)


    pst.parrep(os.path.join(t_d,"test_opt.par"))
    pst.pestpp_options["hotstart_resfile"] = "test_opt.1.sim.rei"
    pst.pestpp_options["base_jacobian"] = "test_opt.1.jcb"
    pst.control_data.noptmax = 1
    pst.write(os.path.join(t_d,"test_sqp.pst"))
    pyemu.os_utils.run("{0} {1}.pst".format(exe_path, "test_sqp.pst"), cwd=t_d)


def rosenc_test():
    sys.path.append(os.path.join("..","benchmarks"))
    import opt_test_suite_helper as helper

    t_d = helper.setup_problem("rosenc")
    model_d = "mou_tests"
    pst = pyemu.Pst(os.path.join(t_d,"rosenc.pst"))
    print(pst.pestpp_options)

    pst.pestpp_options["opt_objective_function"] = "obj_1"
    pst.observation_data.loc["obj_1", "obgnme"] = "obj"
    # pst.parameter_data.loc[:,"standard_deviation"] = np.nan
    # pst.parameter_data.loc[["dv_0","dv_1"], "standard_deviation"] = 0.1
    # pst.parameter_data.loc["dv_0","parval1"] = -0.052
    # pst.parameter_data.loc["dv_1","parval1"] = -0.1
    pst.pestpp_options["par_sigma_range"] = 25
    pst.parameter_data.loc["dv_0", "parval1"] = -3
    pst.parameter_data.loc["dv_1", "parval1"] = -3
    pst.control_data.noptmax = 10
    pst.pestpp_options["sqp_num_reals"] = 5
    pst.pestpp_options["opt_direction"] = "min"
    pst.write(os.path.join(t_d, "rosenc.pst"),version=2)
    m_d = os.path.join(model_d, "master_rosenc_enopt")
    if not "sqp_num_reals" in pst.pestpp_options:
        m_d = os.path.join(model_d, "master_rosenc_opt")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pyemu.os_utils.start_workers(t_d,exe_path,"rosenc.pst",worker_root=model_d,num_workers=10,master_dir=m_d)
    return m_d


# def rosen_test():
#     sys.path.append(os.path.join("..","benchmarks"))
#     import opt_test_suite_helper as helper

#     t_d = helper.setup_problem("rosen")
#     model_d = "mou_tests"
#     pst = pyemu.Pst(os.path.join(t_d,"rosen.pst"))
#     print(pst.pestpp_options)

#     pst.pestpp_options["opt_objective_function"] = "obj_1"
#     pst.observation_data.loc["obj_1","obgnme"] = "obj"
#     #pst.parameter_data.loc[:,"standard_deviation"] = np.nan
#     #pst.parameter_data.loc[["dv_0","dv_1"], "standard_deviation"] = 0.1
#     #pst.parameter_data.loc["dv_0","parval1"] = -0.052
#     #pst.parameter_data.loc["dv_1","parval1"] = -0.1
#     pst.pestpp_options["par_sigma_range"] = 25
#     pst.parameter_data.loc["dv_0", "parval1"] = -3
#     pst.parameter_data.loc["dv_1", "parval1"] = -3
#     pst.control_data.noptmax = 10
#     pst.pestpp_options["sqp_num_reals"] = 5
#     pst.pestpp_options["opt_direction"] = "min"
#     pst.write(os.path.join(t_d, "rosen.pst"),version=2)

#     m_d = os.path.join(model_d, "master_rosen_enopt")
#     if not "sqp_num_reals" in pst.pestpp_options:
#         m_d = os.path.join(model_d, "master_rosen_opt")

#     if os.path.exists(m_d):
#         shutil.rmtree(m_d)
#     pyemu.os_utils.start_workers(t_d,exe_path,"rosen.pst",worker_root=model_d,num_workers=10,master_dir=m_d)
#     return m_d


def plot_rosen(m_d):
    import matplotlib.pyplot as plt
    case = os.path.split(m_d)[-1].split("_")[1]
    par_files = [f for f in os.listdir(m_d) if f.endswith(".base.par") and len(f.split(".")) == 4 and f.startswith(case)]
    par_iter = [int(f.split('.')[1]) for f in par_files]
    par_dict = {i:f for i,f in zip(par_iter,par_files)}
    par_iter.sort()

    pst = pyemu.Pst(os.path.join(m_d,case+".pst"))

    pnames = pst.par_names
    ubnd = pst.parameter_data.parubnd.to_dict()
    lbnd = pst.parameter_data.parlbnd.to_dict()

    sys.path.append(m_d)
    import forward_run

    if case == "rosenc":
        fxn = forward_run.rosenc
    elif case == "rosen":
        fxn = forward_run.rosen

    vals = []
    allvals = []
    isconst = False
    if "rosenc" in m_d:
        isconst = True
    p1_vals = np.linspace(lbnd[pnames[0]],ubnd[pnames[0]],500)
    p2_vals = np.linspace(lbnd[pnames[1]],ubnd[pnames[1]],500)
    for p1 in p1_vals:
        for p2 in p2_vals:
            objs,cnst = fxn([p1,p2])
            if isconst and max(cnst) > 0:
                vals.append(np.NaN)
            else:
                vals.append(objs[0])
            allvals.append(objs[0])
        #print(p1)
    vals = np.array(vals)
    vals = vals.reshape(p1_vals.shape[0],p2_vals.shape[0])
    vals = np.log10(vals)
    vals = vals.transpose()
    allvals = np.array(allvals)
    allvals = allvals.reshape(p1_vals.shape[0], p2_vals.shape[0])
    allvals = np.log10(allvals)
    allvals = allvals.transpose()
    P1,P2 = np.meshgrid(p1_vals,p2_vals)


    pp1,pp2 = [],[]
    ppp1, ppp2 = [], []

    for i in par_iter:
        par = pyemu.pst_utils.read_parfile(os.path.join(m_d,par_dict[i]))

        print(i)
        fig = plt.figure(figsize=(6,6))
        gs = fig.add_gridspec(5,5)
        ax = fig.add_subplot(gs[:-1,1:])
        #ax.imshow(vals,extent=(lbnd[pnames[0]],ubnd[pnames[0]],lbnd[pnames[1]],ubnd[pnames[1]]),cmap="jet",alpha=0.1)
        cb = ax.pcolormesh(P1,P2,vals,cmap="jet",alpha=0.5)
        ax.contour(P1,P2,allvals,levels=5,colors="k",linewidths=0.3)
        #plt.colorbar(cb,ax=ax,orientation="vertical",label="log $\phi$")
        pe_file = os.path.join(m_d,par_dict[i].replace(".base.par",".par.csv"))
        if os.path.exists(pe_file):
            pe = pd.read_csv(os.path.join(m_d, par_dict[i].replace(".base.par",".par.csv")), index_col=0)
            ax.scatter(pe.loc[:,pnames[0]].values,pe.loc[:,pnames[1]].values,c="m",marker="*",s=30,label="ensemble")
            ax.scatter(pe.loc[:, pnames[0]].values.mean(), pe.loc[:, pnames[1]].values.mean(),marker="^",c="r", s=50,label="ensemble mean")


            ax1 = fig.add_subplot(gs[-1, 1:])
            # ax1 = fig.add_subplot(gs[:-1,0])
            # bins = np.linspace(lbnd[pnames[0]],ubnd[pnames[0]],10)
            bins = np.linspace(pe.loc[:, pnames[0]].min(), pe.loc[:, pnames[0]].max(), 10)
            ax1.hist(pe.loc[:, pnames[0]].values, bins=bins, alpha=0.5, facecolor="0.5")  # ,orientation='horizontal')
            ax1.set_xlim(lbnd[pnames[0]], ubnd[pnames[0]])
            ax1.set_yticks([])
            ax1.set_xlabel("dv1")

            ax2 = fig.add_subplot(gs[:-1, 0])
            bins = np.linspace(pe.loc[:, pnames[1]].min(), pe.loc[:, pnames[1]].max(), 10)
            ax2.hist(pe.loc[:, pnames[1]].values, bins=bins, alpha=0.5, facecolor="0.5", orientation='horizontal')
            ax2.set_ylim(lbnd[pnames[1]], ubnd[pnames[1]])
            ax2.set_xticks([])
            ax2.set_ylabel("dv2")
            ax.set_xticklabels([])
            ax.set_yticklabels([])
        else:
            ax.set_xlabel("dv1")
            ax.set_ylabel("dv2")

        ax.scatter(par.loc[pnames[0],"parval1"], par.loc[pnames[1],"parval1"], marker="o", c="b", s=30,label="current solution")
        if i+1 in par_dict:
            can_file = os.path.join(m_d,par_dict[i+1].replace(".base.par",".dv_candidates.csv"))
            if os.path.exists(can_file):
                can = pd.read_csv(can_file,index_col=0)
                ax.scatter(can.loc[:,pnames[0]],can.loc[:,pnames[1]],marker="+",c='0.5',s=10,label="candidate solution",alpha=0.5)

        #pp1.append(pe.loc[:, pnames[0]].values.mean())
        #pp2.append(pe.loc[:, pnames[1]].values.mean())
        ppp1.append(par.loc[pnames[0],"parval1"])
        ppp2.append(par.loc[pnames[1],"parval1"])
        #ax.plot(pp1,pp2,"r-",lw=0.5,)
        ax.plot(ppp1, ppp2, "b-", lw=0.5, )
        ax.set_xlim(lbnd[pnames[0]],ubnd[pnames[0]])
        ax.set_ylim(lbnd[pnames[1]], ubnd[pnames[1]])
        ax.grid()
        ax.legend(loc="upper left")




        ax.set_title("iteration {0}".format(i),loc="left")
        last = os.path.join(m_d,"iter_{0:03d}.png".format(i))
        plt.savefig(last)
        plt.close(fig)

    prefix = "iter_%03d.png"
    last = os.path.split(last)[-1]
    #pyemu.os_utils.run("ffmpeg -i "+prefix+" -r 1 -y out.gif",cwd=m_d)
    #pyemu.os_utils.run("ffmpeg -i {0} -vf palettegen=256 -y palette.png".format(last),cwd=m_d)
    pyemu.os_utils.run("ffmpeg -f image2 -framerate 1 "
                       "-i iter_%03d.png -vf scale=480X480  -y out.gif",cwd=m_d)


if __name__ == "__main__":

    #if not os.path.exists(os.path.join("..","bin")):
    #    os.mkdir(os.path.join("..","bin"))
    #shutil.copy2(os.path.join("..","exe","windows","x64","Debug","pestpp-sqp.exe"),os.path.join("..","bin","pestpp-sqp.exe"))
    basic_sqp_test()
    #rosenbrock_single_linear_constraint(nit=1)
    #dewater_basic_test()
    #dewater_slp_opt_test()
    #rosenc_test()
    #m_d = rosenc_test()
    #m_d = os.path.join("mou_tests","master_rosenc_enopt")
    #plot_rosen(m_d)
