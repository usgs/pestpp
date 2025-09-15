import os
import sys
import shutil
import platform
import numpy as np
import pandas as pd
import platform
import pyemu
sys.path.append(os.path.join("..","benchmarks"))
import opt_test_suite_helper as mou_suite_helper

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
exe_path = os.path.join(bin_path, "pestpp-mou" + exe)


noptmax = 4
num_reals = 20
port = 4021
test_root = "mou_tests"



def test_sorting_fake_problem():
    test_d = os.path.join(test_root,"sorting_test")
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    os.makedirs(test_d)

    with open(os.path.join(test_d,"par.tpl"),'w') as f:
        f.write("ptf ~\n ~   par1    ~\n")
        f.write("~   par2    ~\n")
        
    with open(os.path.join(test_d,"obs.ins"),'w') as f:
        f.write("pif ~\n")
        for i in range(6): # the number of objs in the test
            f.write("l1 !obj{0}!\n".format(i))
    pst = pyemu.Pst.from_io_files(os.path.join(test_d,"par.tpl"),"par.dat",os.path.join(test_d,"obs.ins"),"obs.dat",pst_path=".")
    obs = pst.observation_data
    obs.loc[:,"obgnme"] = "less_than_obj"
    obs.loc[:,"obsval"] = 0.0
    obs.loc[:,"weight"] = 1.0
    
    par = pst.parameter_data
    par.loc[:,"partrans"] = "none"
    par.loc[:,"parval1"] = 1.0
    par.loc[:,"parubnd"] = 1.5
    par.loc[:,"parlbnd"] = 0.5

    pst.control_data.noptmax = -1
    np.random.seed(111)

    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst,num_reals=50)
    oe = pyemu.ObservationEnsemble.from_gaussian_draw(pst=pst,num_reals=50)

    pe.to_csv(os.path.join(test_d,"par.csv"))
    oe.to_csv(os.path.join(test_d,"obs.csv"))
    pst.pestpp_options["mou_dv_population_file"] = "par.csv"
    pst.pestpp_options["mou_obs_population_restart_file"] = "obs.csv"
    pst.write(os.path.join(test_d,"test.pst"))
    pyemu.os_utils.run("{0} test.pst".format(exe_path),cwd=test_d)

    cov = pyemu.Cov.from_parameter_data(pst).to_2d()
    cov.x[0,1] = 0.0001
    cov.x[1,0] = 0.0001


    pyemu.helpers.first_order_pearson_tikhonov(pst=pst,cov=cov,abs_drop_tol=0.0)
    print(pst.prior_information)
    #pst.prior_information = pst.prior_information.loc[["pcc_3"],:]
    pi = pst.prior_information
    pi.loc["pcc_1","equation"] = pi.loc["pcc_1","equation"].replace("= 0.0","= 1.0").replace(" - "," + ")
    pi.loc[:,"obgnme"] = "less_than_pi"
    pst.write(os.path.join(test_d,"test.pst"))
    pyemu.os_utils.run("{0} test.pst".format(exe_path),cwd=test_d)



def test_risk_obj():
    t_d = mou_suite_helper.setup_problem("zdt1",True,True)
    df = pd.read_csv(os.path.join(t_d,"prior.csv"),index_col=0)
    df.loc[:,"_risk_"] = 0.95
    print(df.columns)
    df.to_csv(os.path.join(t_d,"prior.csv"))
    pst = pyemu.Pst(os.path.join(t_d,"zdt1.pst"))
    pst.pestpp_options["mou_dv_population_file"] = "prior.csv"
    pst.pestpp_options["opt_chance_points"] = "single"
    pst.pestpp_options["opt_recalc_chance_every"] = 1000
    pst.pestpp_options["opt_stack_size"] = 10
    pst.pestpp_options["opt_par_stack"] = "prior.csv"
    pst.pestpp_options["mou_generator"] = "de"
    pst.control_data.noptmax = -1
    pst.write(os.path.join(t_d,"zdt1.pst"))
    m1 = os.path.join("mou_tests","zdt1_test_master_riskobj")
    pyemu.os_utils.start_workers(t_d,exe_path,"zdt1.pst",35,worker_root="mou_tests",
                                 master_dir=m1,verbose=True,port=port)

    t_d = mou_suite_helper.setup_problem("zdt1", True, False)
    df.pop("_risk_")
    df.to_csv(os.path.join(t_d, "prior.csv"))
    pst = pyemu.Pst(os.path.join(t_d, "zdt1.pst"))
    pst.pestpp_options["mou_dv_population_file"] = "prior.csv"
    pst.pestpp_options["opt_chance_points"] = "single"
    pst.pestpp_options["opt_recalc_chance_every"] = 1000
    pst.pestpp_options["opt_risk"] = 0.95
    pst.pestpp_options["opt_stack_size"] = 10
    pst.pestpp_options["opt_par_stack"] = "prior.csv"
    pst.pestpp_options["mou_generator"] = "de"
    pst.control_data.noptmax = -1
    pst.write(os.path.join(t_d, "zdt1.pst"))
    m2 = os.path.join("mou_tests","zdt1_test_master")
    pyemu.os_utils.start_workers(t_d, exe_path, "zdt1.pst", 35,
                                 worker_root="mou_tests",
                                 master_dir=m2,
                                 verbose=True,
                                 port=port)

    test_files = ["zdt1.0.obs_pop.csv","zdt1.0.obs_stack.csv","zdt1.0.chance.obs_pop.csv"]
    for test_file in test_files:
        df1 = pd.read_csv(os.path.join(m1,test_file),index_col=0)
        df2 = pd.read_csv(os.path.join(m2, test_file), index_col=0)
        d = (df1 - df2).apply(lambda x: np.abs(x))
        print(d.max().max())


def test_restart_single():
    t_d = mou_suite_helper.setup_problem("zdt1", True, True)
    pst = pyemu.Pst(os.path.join(t_d, "zdt1.pst"))
    pst.pestpp_options["opt_chance_points"] = "single"
    pst.pestpp_options["opt_recalc_chance_every"] = 1000
    pst.pestpp_options["opt_stack_size"] = 10
    pst.pestpp_options["mou_population_size"] = 10
    pst.pestpp_options["opt_par_stack"] = "prior.csv"
    pst.pestpp_options["mou_generator"] = "de"
    pst.control_data.noptmax = -1
    pst.write(os.path.join(t_d, "zdt1.pst"))
    m1 = os.path.join("mou_tests", "zdt1_test_master_restart1")
    pyemu.os_utils.start_workers(t_d, exe_path, "zdt1.pst", 35, worker_root="mou_tests",
                                 master_dir=m1, verbose=True,port=port)

    shutil.copy2(os.path.join(m1,'zdt1.0.par_stack.csv'),os.path.join(t_d,"par_stack.csv"))
    shutil.copy2(os.path.join(m1, 'zdt1.0.obs_stack.csv'), os.path.join(t_d, "obs_stack.csv"))


    pst.pestpp_options["opt_par_stack"] = "par_stack.csv"
    pst.pestpp_options["opt_obs_stack"] = "obs_stack.csv"
    pst.control_data.noptmax = 3
    pst.pestpp_options["opt_recalc_chance_every"] = 2
    pst.write(os.path.join(t_d, "zdt1.pst"))
    m2 = os.path.join("mou_tests", "zdt1_test_master_restart2")
    pyemu.os_utils.start_workers(t_d, exe_path, "zdt1.pst", 35, worker_root="mou_tests",
                                 master_dir=m2, verbose=True,port=port)

    chance_file = "zdt1.0.chance.obs_pop.csv"
    d1 = pd.read_csv(os.path.join(m1,chance_file),index_col=0)
    d2 = pd.read_csv(os.path.join(m2, chance_file), index_col=0)
    d = (d1-d2).apply(lambda x: np.abs(x))
    print(d.max().max())
    assert d.max().max() < 0.01

def test_restart_all():
    t_d = mou_suite_helper.setup_problem("zdt1", True, True)
    pst = pyemu.Pst(os.path.join(t_d, "zdt1.pst"))
    pst.pestpp_options["opt_chance_points"] = "all"
    pst.pestpp_options["opt_recalc_chance_every"] = 1000
    pst.pestpp_options["opt_stack_size"] = 5
    pst.pestpp_options["mou_population_size"] = 10
    pst.pestpp_options["opt_par_stack"] = "prior.csv"
    pst.pestpp_options["mou_generator"] = "de"
    pst.control_data.noptmax = -1
    pst.write(os.path.join(t_d, "zdt1.pst"))
    m1 = os.path.join("mou_tests", "zdt1_test_master_restart1")
    pyemu.os_utils.start_workers(t_d, exe_path, "zdt1.pst", 35, worker_root="mou_tests",
                                 master_dir=m1, verbose=True,port=port)

    shutil.copy2(os.path.join(m1, 'zdt1.0.nested.par_stack.csv'), os.path.join(t_d, "par_stack.csv"))
    shutil.copy2(os.path.join(m1, 'zdt1.0.nested.obs_stack.csv'), os.path.join(t_d, "obs_stack.csv"))

    pst.pestpp_options["opt_par_stack"] = "par_stack.csv"
    pst.pestpp_options["opt_obs_stack"] = "obs_stack.csv"
    pst.control_data.noptmax = 3
    pst.pestpp_options["opt_recalc_chance_every"] = 2
    pst.write(os.path.join(t_d, "zdt1.pst"))
    m2 = os.path.join("mou_tests", "zdt1_test_master_restart2")
    pyemu.os_utils.start_workers(t_d, exe_path, "zdt1.pst", 35, worker_root="mou_tests",
                                 master_dir=m2, verbose=True,port=port)

    chance_file = "zdt1.0.chance.obs_pop.csv"
    d1 = pd.read_csv(os.path.join(m1, chance_file), index_col=0)
    d2 = pd.read_csv(os.path.join(m2, chance_file), index_col=0)
    d = (d1 - d2).apply(lambda x: np.abs(x))
    print(d.max().max())
    assert d.max().max() < 0.01

    t_d = mou_suite_helper.setup_problem("zdt1", True, True)
    pst = pyemu.Pst(os.path.join(t_d, "zdt1.pst"))
    pst.pestpp_options["opt_chance_points"] = "all"
    pst.pestpp_options["opt_recalc_chance_every"] = 1000
    pst.pestpp_options["opt_stack_size"] = 5
    pst.pestpp_options["mou_population_size"] = 10
    pst.pestpp_options["opt_par_stack"] = "prior.csv"
    pst.pestpp_options["mou_generator"] = "pso"
    pst.control_data.noptmax = -1
    pst.write(os.path.join(t_d, "zdt1.pst"))
    m1 = os.path.join("mou_tests", "zdt1_test_master_restart1")
    pyemu.os_utils.start_workers(t_d, exe_path, "zdt1.pst", 35, worker_root="mou_tests",
                                 master_dir=m1, verbose=True,port=port)

    shutil.copy2(os.path.join(m1, 'zdt1.0.nested.par_stack.csv'), os.path.join(t_d, "par_stack.csv"))
    shutil.copy2(os.path.join(m1, 'zdt1.0.nested.obs_stack.csv'), os.path.join(t_d, "obs_stack.csv"))

    pst.pestpp_options["opt_par_stack"] = "par_stack.csv"
    pst.pestpp_options["opt_obs_stack"] = "obs_stack.csv"
    pst.control_data.noptmax = 3
    pst.pestpp_options["opt_recalc_chance_every"] = 2
    pst.write(os.path.join(t_d, "zdt1.pst"))
    m2 = os.path.join("mou_tests", "zdt1_test_master_restart2")
    pyemu.os_utils.start_workers(t_d, exe_path, "zdt1.pst", 35, worker_root="mou_tests",
                                 master_dir=m2, verbose=True,port=port)

    chance_file = "zdt1.0.chance.obs_pop.csv"
    d1 = pd.read_csv(os.path.join(m1, chance_file), index_col=0)
    d2 = pd.read_csv(os.path.join(m2, chance_file), index_col=0)
    d = (d1 - d2).apply(lambda x: np.abs(x))
    print(d.max().max())
    assert d.max().max() < 0.01

def invest_risk_obj():
    t_d = mou_suite_helper.setup_problem("zdt1",True,True)
    pst = pyemu.Pst(os.path.join(t_d,"zdt1.pst"))
    pst.pestpp_options["opt_chance_points"] = "single"
    pst.pestpp_options["opt_recalc_chance_every"] = 1000
    pst.pestpp_options["opt_stack_size"] = 50
    pst.pestpp_options["mou_generator"] = "de"
    pst.control_data.noptmax = 300
    pst.write(os.path.join(t_d,"zdt1.pst"))
    m1 = os.path.join("mou_tests","zdt1_test_master_riskobj_full")
    pyemu.os_utils.start_workers(t_d,exe_path,"zdt1.pst",35,worker_root="mou_tests",
                                 master_dir=m1,verbose=True,port=port)
    plot_results(os.path.join("mou_tests", "zdt1_test_master_riskobj_full"))

def chance_consistency_test():
    t_d = mou_suite_helper.setup_problem("constr", additive_chance=True, risk_obj=False)
    pst = pyemu.Pst(os.path.join(t_d, "constr.pst"))
    pst.pestpp_options["opt_chance_points"] = "all"
    pst.pestpp_options["opt_recalc_chance_every"] = 5
    pst.pestpp_options["opt_stack_size"] = 10
    pst.pestpp_options["mou_generator"] = "de"
    pst.pestpp_options["mou_population_size"] = 10
    pst.pestpp_options["opt_risk"] = 0.95
    pst.control_data.noptmax = -1
    pst.write(os.path.join(t_d, "constr.pst"))
    m1 = os.path.join("mou_tests", "constr_test_master_fail_3")
    pyemu.os_utils.start_workers(t_d, exe_path, "constr.pst", 35, worker_root="mou_tests",
                                 master_dir=m1, verbose=True,port=port)


    
    tag = "chance-shifted initial population objective function summary"
    with open(os.path.join(m1,"constr.rec"),'r') as f:
        while True:
            line = f.readline()
            if line == "":
                raise Exception()
            if tag in line:
                line = f.readline()
                mname = line.strip().split()[2].lower()
                for _ in range(3):
                    f.readline()
                o1 = float(f.readline().strip().split()[-1])
                o2 = float(f.readline().strip().split()[-1])
                break
        

    print(mname,o1,o2)
    df = pd.read_csv(os.path.join(m1,"constr.0.chance.obs_pop.csv"),index_col=0)
    d1 = np.abs(df.loc[mname,:].iloc[0] - o1)
    d2 = np.abs(df.loc[mname,:].iloc[1] - o2)
    print(mname,o1,df.loc[mname,:].iloc[0],d1,o2,df.loc[mname,:].iloc[1],d2)
    assert d1 < 1.e-5
    assert d2 < 1.e-5
    
def fail_test():
    t_d = mou_suite_helper.setup_problem("zdt1", additive_chance=True, risk_obj=False)
    pst = pyemu.Pst(os.path.join(t_d, "zdt1.pst"))
    pst.pestpp_options["opt_chance_points"] = "all"
    pst.pestpp_options["opt_recalc_chance_every"] = 5
    pst.pestpp_options["opt_stack_size"] = 10
    pst.pestpp_options["mou_generator"] = "de"
    pst.pestpp_options["mou_population_size"] = 30
    pst.pestpp_options["ies_debug_fail_remainder"] = True
    pst.pestpp_options["opt_risk"] = 0.95
    pst.control_data.noptmax = 4
    pst.write(os.path.join(t_d, "zdt1.pst"))
    m1 = os.path.join("mou_tests", "zdt1_test_master_fail_1")
    pyemu.os_utils.start_workers(t_d, exe_path, "zdt1.pst", 35, worker_root="mou_tests",
                                 master_dir=m1, verbose=True,port=port)

    dp_file = os.path.join(m1,"zdt1.0.nested.par_stack.csv")
    dp = pd.read_csv(dp_file,index_col=0)
    op_file = os.path.join(m1,"zdt1.0.nested.obs_stack.csv")
    op = pd.read_csv(dp_file,index_col=0)
    print(dp.shape,op.shape)
    assert dp.shape[0] == op.shape[0]
    
    shutil.copy2(dp_file,os.path.join(t_d,"par_stack.csv"))
    shutil.copy2(op_file,os.path.join(t_d,"obs_stack.csv"))
    pst.pestpp_options["opt_par_stack"] = "par_stack.csv"
    pst.pestpp_options["opt_obs_stack"] = "obs_stack.csv"
    pst.control_data.noptmax = 2
    pst.write(os.path.join(t_d, "zdt1.pst"))
    m1 = os.path.join("mou_tests", "zdt1_test_master_fail_1")
    pyemu.os_utils.start_workers(t_d, exe_path, "zdt1.pst", 35, worker_root="mou_tests",
                                 master_dir=m1, verbose=True,port=port)


    t_d = mou_suite_helper.setup_problem("zdt1", additive_chance=True, risk_obj=False)
    pst = pyemu.Pst(os.path.join(t_d, "zdt1.pst"))
    pst.pestpp_options["opt_chance_points"] = "all"
    pst.pestpp_options["opt_recalc_chance_every"] = 5
    pst.pestpp_options["opt_stack_size"] = 10
    pst.pestpp_options["mou_generator"] = "de"
    pst.pestpp_options["mou_population_size"] = 30
    pst.pestpp_options["ies_debug_fail_subset"] = True
    pst.pestpp_options["opt_risk"] = 0.95
    pst.control_data.noptmax = 4
    pst.write(os.path.join(t_d, "zdt1.pst"))
    m1 = os.path.join("mou_tests", "zdt1_test_master_fail_1")
    pyemu.os_utils.start_workers(t_d, exe_path, "zdt1.pst", 35, worker_root="mou_tests",
                                 master_dir=m1, verbose=True)


    dp_file = os.path.join(m1,"zdt1.0.nested.par_stack.csv")
    dp = pd.read_csv(dp_file,index_col=0)
    op_file = os.path.join(m1,"zdt1.0.nested.obs_stack.csv")
    op = pd.read_csv(dp_file,index_col=0)
    print(dp.shape,op.shape)
    assert dp.shape[0] == op.shape[0]
    
    shutil.copy2(dp_file,os.path.join(t_d,"par_stack.csv"))
    shutil.copy2(op_file,os.path.join(t_d,"obs_stack.csv"))
    pst.pestpp_options["opt_par_stack"] = "par_stack.csv"
    pst.pestpp_options["opt_obs_stack"] = "obs_stack.csv"
    pst.control_data.noptmax = 2
    pst.write(os.path.join(t_d, "zdt1.pst"))
    m1 = os.path.join("mou_tests", "zdt1_test_master_fail_1")
    pyemu.os_utils.start_workers(t_d, exe_path, "zdt1.pst", 35, worker_root="mou_tests",
                                 master_dir=m1, verbose=True,port=port)

    t_d = mou_suite_helper.setup_problem("zdt1", additive_chance=True,risk_obj=False)
    pst = pyemu.Pst(os.path.join(t_d, "zdt1.pst"))
    pst.pestpp_options["opt_chance_points"] = "single"
    pst.pestpp_options["opt_recalc_chance_every"] = 1000
    pst.pestpp_options["opt_stack_size"] = 10
    pst.pestpp_options["mou_generator"] = "de"
    pst.pestpp_options["mou_population_size"] = 30
    pst.pestpp_options["ies_debug_fail_subset"] = True
    pst.pestpp_options["opt_risk"] = 0.95
    pst.control_data.noptmax = 4
    pst.write(os.path.join(t_d, "zdt1.pst"))
    m1 = os.path.join("mou_tests", "zdt1_test_master_fail_1")
    pyemu.os_utils.start_workers(t_d, exe_path, "zdt1.pst", 35, worker_root="mou_tests",
                                 master_dir=m1, verbose=True,port=port)


def pso_fail_test():
    t_d = mou_suite_helper.setup_problem("zdt1", additive_chance=True, risk_obj=False)
    pst = pyemu.Pst(os.path.join(t_d, "zdt1.pst"))
    pst.pestpp_options["opt_chance_points"] = "all"
    pst.pestpp_options["opt_recalc_chance_every"] = 5
    pst.pestpp_options["opt_stack_size"] = 10
    pst.pestpp_options["mou_generator"] = "pso"
    pst.pestpp_options["mou_population_size"] = 30
    pst.pestpp_options["ies_debug_fail_remainder"] = True
    pst.pestpp_options["opt_risk"] = 0.95
    pst.control_data.noptmax = 4
    pst.write(os.path.join(t_d, "zdt1.pst"))
    m1 = os.path.join("mou_tests", "zdt1_test_master_fail_1")
    pyemu.os_utils.start_workers(t_d, exe_path, "zdt1.pst", 35, worker_root="mou_tests",
                                 master_dir=m1, verbose=True,port=port)

    dp_file = os.path.join(m1,"zdt1.0.nested.par_stack.csv")
    dp = pd.read_csv(dp_file,index_col=0)
    op_file = os.path.join(m1,"zdt1.0.nested.obs_stack.csv")
    op = pd.read_csv(dp_file,index_col=0)
    print(dp.shape,op.shape)
    assert dp.shape[0] == op.shape[0]
    
    shutil.copy2(dp_file,os.path.join(t_d,"par_stack.csv"))
    shutil.copy2(op_file,os.path.join(t_d,"obs_stack.csv"))
    pst.pestpp_options["opt_par_stack"] = "par_stack.csv"
    pst.pestpp_options["opt_obs_stack"] = "obs_stack.csv"
    pst.control_data.noptmax = 2
    pst.write(os.path.join(t_d, "zdt1.pst"))
    m1 = os.path.join("mou_tests", "zdt1_test_master_fail_1")
    pyemu.os_utils.start_workers(t_d, exe_path, "zdt1.pst", 35, worker_root="mou_tests",
                                 master_dir=m1, verbose=True,port=port)


    t_d = mou_suite_helper.setup_problem("zdt1", additive_chance=True, risk_obj=False)
    pst = pyemu.Pst(os.path.join(t_d, "zdt1.pst"))
    pst.pestpp_options["opt_chance_points"] = "all"
    pst.pestpp_options["opt_recalc_chance_every"] = 5
    pst.pestpp_options["opt_stack_size"] = 10
    pst.pestpp_options["mou_generator"] = "pso"
    pst.pestpp_options["mou_population_size"] = 30
    pst.pestpp_options["ies_debug_fail_subset"] = True
    pst.pestpp_options["opt_risk"] = 0.95
    pst.control_data.noptmax = 4
    pst.write(os.path.join(t_d, "zdt1.pst"))
    m1 = os.path.join("mou_tests", "zdt1_test_master_fail_1")
    pyemu.os_utils.start_workers(t_d, exe_path, "zdt1.pst", 35, worker_root="mou_tests",
                                 master_dir=m1, verbose=True)


    dp_file = os.path.join(m1,"zdt1.0.nested.par_stack.csv")
    dp = pd.read_csv(dp_file,index_col=0)
    op_file = os.path.join(m1,"zdt1.0.nested.obs_stack.csv")
    op = pd.read_csv(dp_file,index_col=0)
    print(dp.shape,op.shape)
    assert dp.shape[0] == op.shape[0]
    
    shutil.copy2(dp_file,os.path.join(t_d,"par_stack.csv"))
    shutil.copy2(op_file,os.path.join(t_d,"obs_stack.csv"))
    pst.pestpp_options["opt_par_stack"] = "par_stack.csv"
    pst.pestpp_options["opt_obs_stack"] = "obs_stack.csv"
    pst.control_data.noptmax = 2
    pst.write(os.path.join(t_d, "zdt1.pst"))
    m1 = os.path.join("mou_tests", "zdt1_test_master_fail_1")
    pyemu.os_utils.start_workers(t_d, exe_path, "zdt1.pst", 35, worker_root="mou_tests",
                                 master_dir=m1, verbose=True,port=port)

    t_d = mou_suite_helper.setup_problem("zdt1", additive_chance=True,risk_obj=False)
    pst = pyemu.Pst(os.path.join(t_d, "zdt1.pst"))
    pst.pestpp_options["opt_chance_points"] = "single"
    pst.pestpp_options["opt_recalc_chance_every"] = 1000
    pst.pestpp_options["opt_stack_size"] = 10
    pst.pestpp_options["mou_generator"] = "pso"
    pst.pestpp_options["mou_population_size"] = 30
    pst.pestpp_options["ies_debug_fail_subset"] = True
    pst.pestpp_options["opt_risk"] = 0.95
    pst.control_data.noptmax = 4
    pst.write(os.path.join(t_d, "zdt1.pst"))
    m1 = os.path.join("mou_tests", "zdt1_test_master_fail_1")
    pyemu.os_utils.start_workers(t_d, exe_path, "zdt1.pst", 35, worker_root="mou_tests",
                                 master_dir=m1, verbose=True,port=port)



def invest_5():
    cases = ["kur","sch","srn","tkn","constr","zdt2","zdt3","zdt4","zdt6"]
    noptmax = 100
    for case in cases:
        #t_d = mou_suite_helper.setup_problem("fon")
        mou_suite_helper.run_problem(case,noptmax=noptmax)
        # mou_suite_helper.run_problem_chance(case, noptmax=noptmax,risk_obj=True,chance_points="all",
        #                                     recalc=10000)
        mou_suite_helper.run_problem_chance(case, noptmax=noptmax, risk_obj=True, chance_points="single",
                                            recalc=10000)
        mou_suite_helper.run_problem_chance(case, noptmax=noptmax, risk_obj=False, chance_points="all",
                                            recalc=10000)
        mou_suite_helper.run_problem_chance(case, noptmax=noptmax, risk_obj=False, chance_points="single",
                                            recalc=10000)

def invest_2():
    cases = ["zdt1","zdt2","zdt3","zdt4","zdt6"]
    noptmax = 100
    for case in cases:
        #t_d = mou_suite_helper.setup_problem("fon")

        mou_suite_helper.run_problem(case, noptmax=noptmax, generator="sbx")
        mou_suite_helper.run_problem(case, noptmax=noptmax, generator="de")

        # mou_suite_helper.run_problem_chance(case, noptmax=noptmax,risk_obj=True,chance_points="all",
        #                                     recalc=10000)

def invest_3():
    # cases = ["zdt2", "zdt3", "zdt4", "zdt6", "sch", "srn", "tkn", "constr"]
    # for case in cases:
    #     mou_suite_helper.run_problem(case, generator="sbx", env="spea")
    #     mou_suite_helper.run_problem(case, generator="sbx", env="nsga")
    #     mou_suite_helper.run_problem(case, generator="de,sbx,pm", env="nsga")
    #     mou_suite_helper.run_problem(case, generator="de,sbx,pm", env="spea")
    cases = ["tkn", "constr","zdt1", "zdt3"]
    for case in cases:
        mou_suite_helper.run_problem(case, generator="de", env="spea",self_adaptive=True)
        mou_suite_helper.run_problem(case, generator="de", env="nsga",self_adaptive=True)
        # mou_suite_helper.run_problem(case, generator="pm", env="spea",self_adaptive=True)
        # mou_suite_helper.run_problem(case, generator="pm", env="nsga",self_adaptive=True)
        # mou_suite_helper.run_problem(case, generator="sbx", env="spea",self_adaptive=True)
        # mou_suite_helper.run_problem(case, generator="sbx", env="nsga",self_adaptive=True)
        # mou_suite_helper.run_problem(case, generator="de,sbx,pm", env="nsga",self_adaptive=True)
        # mou_suite_helper.run_problem(case, generator="de,sbx,pm", env="spea",self_adaptive=True)
        # mou_suite_helper.run_problem(case, generator="de", env="spea")
        # mou_suite_helper.run_problem(case, generator="de", env="nsga")
        # mou_suite_helper.run_problem(case, generator="pm", env="spea")
        # mou_suite_helper.run_problem(case, generator="pm", env="nsga")
        # mou_suite_helper.run_problem(case, generator="sbx", env="spea")
        # mou_suite_helper.run_problem(case, generator="sbx", env="nsga")
        # mou_suite_helper.run_problem(case, generator="de,sbx,pm", env="nsga")
        # mou_suite_helper.run_problem(case, generator="de,sbx,pm", env="spea")

def all_infeas_test():
    t_d = mou_suite_helper.setup_problem("tkn")
    pst = pyemu.Pst(os.path.join(t_d,"tkn.pst"))
    obs = pst.observation_data
    print(obs)
    obs.loc["const_1","obsval"] = -1e10
    pst.pestpp_options["mou_population_size"] = 15
    pst.pestpp_options["mou_generator"] = "de"
    pst.pestpp_options["mou_env_selector"] = "spea"
    pst.control_data.noptmax = 2
    pst.write(os.path.join(t_d,"tkn.pst"))
    m1 = os.path.join("mou_tests", "test_master_all_infeas")
    pyemu.os_utils.start_workers(t_d, exe_path, "tkn.pst", 15, worker_root="mou_tests",
                                 master_dir=m1, verbose=True,port=port)


    out_file = os.path.join(m1,"tkn.obs_pop.csv".format(pst.control_data.noptmax))
    assert os.path.exists(out_file)
    df = pd.read_csv(out_file)
    assert df.shape[0] == pst.pestpp_options["mou_population_size"]

    pst.pestpp_options["mou_generator"] = "sbx"
    pst.pestpp_options["mou_env_selector"] = "spea"
    pst.write(os.path.join(t_d,"tkn.pst"))
    m1 = os.path.join("mou_tests", "test_master_all_infeas")
    pyemu.os_utils.start_workers(t_d, exe_path, "tkn.pst", 15, worker_root="mou_tests",
                                 master_dir=m1, verbose=True,port=port)


    out_file = os.path.join(m1,"tkn.obs_pop.csv".format(pst.control_data.noptmax))
    assert os.path.exists(out_file)
    df = pd.read_csv(out_file)
    assert df.shape[0] == pst.pestpp_options["mou_population_size"]

    pst.pestpp_options["mou_generator"] = "de"
    pst.pestpp_options["mou_env_selector"] = "nsga"
    pst.write(os.path.join(t_d,"tkn.pst"))
    m1 = os.path.join("mou_tests", "test_master_all_infeas")
    pyemu.os_utils.start_workers(t_d, exe_path, "tkn.pst", 15, worker_root="mou_tests",
                                 master_dir=m1, verbose=True,port=port)


    out_file = os.path.join(m1,"tkn.obs_pop.csv".format(pst.control_data.noptmax))
    assert os.path.exists(out_file)
    df = pd.read_csv(out_file)
    assert df.shape[0] == pst.pestpp_options["mou_population_size"]

    pst.pestpp_options["mou_generator"] = "pso"
    pst.pestpp_options["mou_env_selector"] = "nsga"
    pst.write(os.path.join(t_d,"tkn.pst"))
    m1 = os.path.join("mou_tests", "test_master_all_infeas")
    pyemu.os_utils.start_workers(t_d, exe_path, "tkn.pst", 15, worker_root="mou_tests",
                                 master_dir=m1, verbose=True,port=port)


    out_file = os.path.join(m1,"tkn.obs_pop.csv".format(pst.control_data.noptmax))
    assert os.path.exists(out_file)
    df = pd.read_csv(out_file)
    assert df.shape[0] == pst.pestpp_options["mou_population_size"]

def invest_4():
    #mou_suite_helper.run_problem(test_case="zdt1", pop_size=100, noptmax=100, generator="de", env="nsga", self_adaptive=False)
    #mou_suite_helper.run_problem(test_case="zdt1", pop_size=100, noptmax=100, generator="de", env="nsga",
    #                             self_adaptive=True)
    mou_suite_helper.plot_results(os.path.join("mou_tests","zdt1_master_generator=de_env=nsga_popsize=100_risk=0.5_riskobj=False_adaptive=False"))
    mou_suite_helper.plot_results(os.path.join("mou_tests",
                                               "zdt1_master_generator=de_env=nsga_popsize=100_risk=0.5_riskobj=False_adaptive=True"))


def restart_dv_test():

    t_d = mou_suite_helper.setup_problem("tkn")
    pst = pyemu.Pst(os.path.join(t_d, "tkn.pst"))
    obs = pst.observation_data
    print(obs)
    obs.loc["const_1", "obsval"] = -1e10
    pst.pestpp_options["mou_population_size"] = 15
    pst.pestpp_options["mou_generator"] = "de"
    pst.pestpp_options["mou_env_selector"] = "spea"
    pst.control_data.noptmax = 2
    pst.write(os.path.join(t_d, "tkn.pst"))
    pyemu.os_utils.run("{0} tkn.pst".format(exe_path), cwd=t_d)
    shutil.copy2(os.path.join(t_d,"tkn.dv_pop.csv".format(pst.control_data.noptmax)),os.path.join(t_d,"restart.csv"))
    pst.pestpp_options["mou_dv_population_file"] = "restart.csv"
    pst.control_data.noptmax = 1
    pst.write(os.path.join(t_d, "tkn.pst"))
    pyemu.os_utils.run("{0} tkn.pst".format(exe_path), cwd=t_d)
    df = pd.read_csv(os.path.join(t_d,"tkn.dv_pop.csv".format(pst.control_data.noptmax)))
    gen_num = df.real_name.apply(lambda x: int(x.split("=")[1].split('_')[0]))
    print(gen_num.max())
    assert gen_num.max() == 3

def chance_all_binary_test():

    t_d = mou_suite_helper.setup_problem("constr", additive_chance=True, risk_obj=False)
    pst = pyemu.Pst(os.path.join(t_d, "constr.pst"))
    pst.pestpp_options["opt_chance_points"] = "all"
    pst.pestpp_options["opt_recalc_chance_every"] = 5
    pst.pestpp_options["opt_stack_size"] = 4
    pst.pestpp_options["mou_generator"] = "de"
    pst.pestpp_options["mou_population_size"] = 5
    pst.pestpp_options["opt_risk"] = 0.95
    pst.pestpp_options["save_binary"] = True
    par = pst.parameter_data
    par.loc["dv_1","partrans"] = "fixed"
    par.loc["obj1_add_par","partrans"] = "fixed"
    pst.control_data.noptmax = 1
    pst.write(os.path.join(t_d, "constr.pst"))
    m1 = os.path.join("mou_tests", "constr_test_master_chance_binary")
    pyemu.os_utils.start_workers(t_d, exe_path, "constr.pst", 35, worker_root="mou_tests",
                                 master_dir=m1, verbose=True,port=port)
    pe = pyemu.ParameterEnsemble.from_binary(pst=pst,filename=os.path.join(m1,"constr.0.nested.par_stack.jcb"))
    oe = pyemu.ObservationEnsemble.from_binary(pst=pst, filename=os.path.join(m1, "constr.0.nested.obs_stack.jcb"))
    dva = pyemu.ParameterEnsemble.from_binary(pst=pst,filename=os.path.join(m1,"constr.archive.dv_pop.jcb"))
    dv = pyemu.ParameterEnsemble.from_binary(pst=pst, filename=os.path.join(m1, "constr.dv_pop.jcb"))
    opa = pyemu.ObservationEnsemble.from_binary(pst=pst, filename=os.path.join(m1, "constr.archive.obs_pop.jcb"))
    op = pyemu.ObservationEnsemble.from_binary(pst=pst, filename=os.path.join(m1, "constr.obs_pop.jcb"))

    shutil.copy(os.path.join(m1,"constr.0.nested.par_stack.jcb"),os.path.join(t_d,"par_stack.jcb"))
    shutil.copy(os.path.join(m1, "constr.0.nested.obs_stack.jcb"), os.path.join(t_d, "obs_stack.jcb"))
    shutil.copy(os.path.join(m1, "constr.0.dv_pop.jcb"), os.path.join(t_d, "dv_pop.jcb"))
    pst.pestpp_options["opt_par_stack"] = "par_stack.jcb"
    pst.write(os.path.join(t_d, "constr.pst"))
    m2 = os.path.join("mou_tests", "constr_test_master_chance_binary_restart")
    pyemu.os_utils.start_workers(t_d, exe_path, "constr.pst", 35, worker_root="mou_tests",
                                 master_dir=m2, verbose=True, port=port)
    pe = pyemu.ParameterEnsemble.from_binary(pst=pst, filename=os.path.join(m2, "constr.0.nested.par_stack.jcb"))
    oe = pyemu.ObservationEnsemble.from_binary(pst=pst, filename=os.path.join(m2, "constr.0.nested.obs_stack.jcb"))
    dva = pyemu.ParameterEnsemble.from_binary(pst=pst,filename=os.path.join(m2,"constr.archive.dv_pop.jcb"))
    dv = pyemu.ParameterEnsemble.from_binary(pst=pst, filename=os.path.join(m2, "constr.dv_pop.jcb"))
    opa = pyemu.ObservationEnsemble.from_binary(pst=pst, filename=os.path.join(m2, "constr.archive.obs_pop.jcb"))
    op = pyemu.ObservationEnsemble.from_binary(pst=pst, filename=os.path.join(m2, "constr.obs_pop.jcb"))

    shutil.copy(os.path.join(m1, "constr.0.nested.par_stack.jcb"), os.path.join(t_d, "par_stack.jcb"))
    shutil.copy(os.path.join(m1, "constr.0.nested.obs_stack.jcb"), os.path.join(t_d, "obs_stack.jcb"))
    shutil.copy(os.path.join(m1, "constr.0.dv_pop.jcb"), os.path.join(t_d, "dv_pop.jcb"))
    pst.pestpp_options["opt_par_stack"] = "par_stack.jcb"
    pst.write(os.path.join(t_d, "constr.pst"))
    m2 = os.path.join("mou_tests", "constr_test_master_chance_binary_restart")
    pyemu.os_utils.start_workers(t_d, exe_path, "constr.pst", 35, worker_root="mou_tests",
                                 master_dir=m2, verbose=True, port=port)
    pe = pyemu.ParameterEnsemble.from_binary(pst=pst, filename=os.path.join(m2, "constr.0.nested.par_stack.jcb"))
    oe = pyemu.ObservationEnsemble.from_binary(pst=pst, filename=os.path.join(m2, "constr.0.nested.obs_stack.jcb"))
    dva = pyemu.ParameterEnsemble.from_binary(pst=pst, filename=os.path.join(m2, "constr.archive.dv_pop.jcb"))
    dv = pyemu.ParameterEnsemble.from_binary(pst=pst, filename=os.path.join(m2, "constr.dv_pop.jcb"))
    opa = pyemu.ObservationEnsemble.from_binary(pst=pst, filename=os.path.join(m2, "constr.archive.obs_pop.jcb"))
    op = pyemu.ObservationEnsemble.from_binary(pst=pst, filename=os.path.join(m2, "constr.obs_pop.jcb"))



def risk_demo(case="zdt1",noptmax=100,std_weight=0.05,mou_gen="de",pop_size=100,num_workers=30):

    obj_names = ["obj_1"]
    if "zdt" in case:
        obj_names.append("obj_2")
    constr_bnd = std_weight * 2

    # spec'd risk = 0.95, stack, reuse
    t_d = mou_suite_helper.setup_problem(case, additive_chance=True, risk_obj=False)
    pst = pyemu.Pst(os.path.join(t_d, case + ".pst"))
    par = pst.parameter_data
    add_par = par.loc[par.parnme.str.contains("_add_"), "parnme"]
    if case == "constr":
        par.loc[add_par, "parubnd"] = constr_bnd
        par.loc[add_par, "parlbnd"] = -constr_bnd

    pst.pestpp_options["opt_chance_points"] = "single"
    pst.pestpp_options["opt_recalc_chance_every"] = 10000
    pst.pestpp_options["opt_stack_size"] = 100
    pst.pestpp_options["mou_generator"] = mou_gen
    pst.pestpp_options["mou_population_size"] = pop_size
    pst.pestpp_options["opt_risk"] = 0.95
    pst.pestpp_options["save_binary"] = True
    pst.observation_data.loc[pst.nnz_obs_names, "weight"] = std_weight
    pst.pestpp_options["opt_std_weights"] = False
    pst.control_data.noptmax = noptmax
    pst.write(os.path.join(t_d, case + ".pst"))
    m1 = os.path.join("mou_tests", case + "_test_master_95_singlept_match_" + mou_gen)
    pyemu.os_utils.start_workers(t_d, exe_path, case + ".pst", num_workers, worker_root="mou_tests",
                                 master_dir=m1, verbose=True, port=port)

    # spec'd risk = 0.05, stack, reuse
    t_d = mou_suite_helper.setup_problem(case, additive_chance=True, risk_obj=False)
    pst = pyemu.Pst(os.path.join(t_d, case + ".pst"))
    par = pst.parameter_data
    add_par = par.loc[par.parnme.str.contains("_add_"), "parnme"]
    if case == "constr":
        par.loc[add_par, "parubnd"] = constr_bnd
        par.loc[add_par, "parlbnd"] = -constr_bnd

    pst.pestpp_options["opt_chance_points"] = "single"
    pst.pestpp_options["opt_recalc_chance_every"] = 10000
    pst.pestpp_options["opt_stack_size"] = 100
    pst.pestpp_options["mou_generator"] = mou_gen
    pst.pestpp_options["mou_population_size"] = pop_size
    pst.pestpp_options["opt_risk"] = 0.05
    pst.pestpp_options["save_binary"] = True
    pst.observation_data.loc[pst.nnz_obs_names, "weight"] = std_weight
    pst.pestpp_options["opt_std_weights"] = False
    pst.control_data.noptmax = noptmax
    pst.write(os.path.join(t_d, case + ".pst"))
    m1 = os.path.join("mou_tests", case + "_test_master_05_singlept_match_" + mou_gen)
    pyemu.os_utils.start_workers(t_d, exe_path, case + ".pst", num_workers, worker_root="mou_tests",
                                 master_dir=m1, verbose=True, port=port)

    # risk obj, stack, all chance pts, reuse
    t_d = mou_suite_helper.setup_problem(case, additive_chance=True, risk_obj=True)
    pst = pyemu.Pst(os.path.join(t_d, case+".pst"))
    par = pst.parameter_data
    add_par = par.loc[par.parnme.str.contains("_add_"),"parnme"]
    if case == "constr":
        par.loc[add_par, "parubnd"] = constr_bnd
        par.loc[add_par, "parlbnd"] = -constr_bnd

    pst.pestpp_options["opt_chance_points"] = "all"
    pst.pestpp_options["opt_recalc_chance_every"] = 10000
    pst.pestpp_options["opt_stack_size"] = 100
    pst.pestpp_options["mou_generator"] = mou_gen
    pst.pestpp_options["mou_population_size"] = pop_size
    pst.pestpp_options["opt_risk"] = 0.95
    pst.pestpp_options["save_binary"] = True
    pst.observation_data.loc[pst.nnz_obs_names,"weight"] = std_weight
    pst.pestpp_options["opt_std_weights"] = False
    pst.control_data.noptmax = noptmax * 3
    pst.write(os.path.join(t_d, case+".pst"))
    m1 = os.path.join("mou_tests", case+"_test_master_riskobj_allpts_more_"+mou_gen)
    pyemu.os_utils.start_workers(t_d, exe_path, case+".pst", num_workers, worker_root="mou_tests",
                                 master_dir=m1, verbose=True, port=port)

    # risk obj, stack, single pt, reuse
    t_d = mou_suite_helper.setup_problem(case, additive_chance=True, risk_obj=True)
    pst = pyemu.Pst(os.path.join(t_d, case+".pst"))
    par = pst.parameter_data
    add_par = par.loc[par.parnme.str.contains("_add_"), "parnme"]
    if case == "constr":
        par.loc[add_par, "parubnd"] = constr_bnd
        par.loc[add_par, "parlbnd"] = -constr_bnd

    pst.pestpp_options["opt_chance_points"] = "single"
    pst.pestpp_options["opt_recalc_chance_every"] = 10000
    pst.pestpp_options["opt_stack_size"] = 100
    pst.pestpp_options["mou_generator"] = mou_gen
    pst.pestpp_options["mou_population_size"] = pop_size
    pst.pestpp_options["opt_risk"] = 0.95
    pst.pestpp_options["save_binary"] = True
    pst.observation_data.loc[pst.nnz_obs_names,"weight"] = std_weight
    pst.pestpp_options["opt_std_weights"] = False
    pst.control_data.noptmax = noptmax * 3
    pst.write(os.path.join(t_d, case+".pst"))
    m1 = os.path.join("mou_tests", case+"_test_master_riskobj_singlept_more_"+mou_gen)
    pyemu.os_utils.start_workers(t_d, exe_path, case+".pst", num_workers, worker_root="mou_tests",
                                 master_dir=m1, verbose=True, port=port)

    # deterministic
    t_d = mou_suite_helper.setup_problem(case, additive_chance=False, risk_obj=False)
    pst = pyemu.Pst(os.path.join(t_d, case+".pst"))
    par = pst.parameter_data
    add_par = par.loc[par.parnme.str.contains("_add_"), "parnme"]
    if case == "constr":
        par.loc[add_par, "parubnd"] = constr_bnd
        par.loc[add_par, "parlbnd"] = -constr_bnd

    pst.pestpp_options["mou_generator"] = mou_gen
    pst.pestpp_options["mou_population_size"] = pop_size
    pst.pestpp_options["opt_risk"] = 0.5
    pst.pestpp_options["save_binary"] = True
    pst.observation_data.loc[pst.nnz_obs_names,"weight"] = std_weight
    pst.control_data.noptmax = noptmax
    pst.write(os.path.join(t_d, case+".pst"))
    m1 = os.path.join("mou_tests", case+"_test_master_deter_"+mou_gen)
    pyemu.os_utils.start_workers(t_d, exe_path, case+".pst", num_workers, worker_root="mou_tests",
                                 master_dir=m1, verbose=True, port=port)


    # spec'd risk = 0.95, std weights
    t_d = mou_suite_helper.setup_problem(case, additive_chance=True, risk_obj=False)
    pst = pyemu.Pst(os.path.join(t_d, case+".pst"))
    par = pst.parameter_data
    add_par = par.loc[par.parnme.str.contains("_add_"), "parnme"]
    if case == "constr":
        par.loc[add_par, "parubnd"] = constr_bnd
        par.loc[add_par, "parlbnd"] = -constr_bnd

    pst.observation_data.loc[pst.nnz_obs_names,"weight"] = std_weight
    pst.pestpp_options["opt_std_weights"] = True
    pst.pestpp_options["mou_generator"] = mou_gen
    pst.pestpp_options["mou_population_size"] = pop_size
    pst.pestpp_options["opt_risk"] = 0.95
    pst.pestpp_options["save_binary"] = True
    pst.control_data.noptmax = noptmax
    pst.write(os.path.join(t_d, case+".pst"))
    m2 = os.path.join("mou_tests", case+"_test_master_95_"+mou_gen)
    pyemu.os_utils.start_workers(t_d, exe_path, case+".pst", num_workers, worker_root="mou_tests",
                                 master_dir=m2, verbose=True, port=port)


    # spec's risk = 0.05, std weights,
    t_d = mou_suite_helper.setup_problem(case, additive_chance=True, risk_obj=False)
    pst = pyemu.Pst(os.path.join(t_d, case+".pst"))
    par = pst.parameter_data
    add_par = par.loc[par.parnme.str.contains("_add_"), "parnme"]
    if case == "constr":
        par.loc[add_par, "parubnd"] = constr_bnd
        par.loc[add_par, "parlbnd"] = -constr_bnd

    pst.observation_data.loc[pst.nnz_obs_names,"weight"] = std_weight
    pst.pestpp_options["opt_std_weights"] = True
    pst.pestpp_options["mou_generator"] = mou_gen
    pst.pestpp_options["mou_population_size"] = pop_size
    pst.pestpp_options["opt_risk"] = 0.05
    pst.pestpp_options["save_binary"] = True
    pst.control_data.noptmax = noptmax
    pst.write(os.path.join(t_d, case+".pst"))
    m3 = os.path.join("mou_tests", case+"_test_master_05_"+mou_gen)
    pyemu.os_utils.start_workers(t_d, exe_path, case+".pst", num_workers, worker_root="mou_tests",
                                 master_dir=m3, verbose=True, port=port)

    # risk obj, std weights
    t_d = mou_suite_helper.setup_problem(case, additive_chance=True, risk_obj=True)
    pst = pyemu.Pst(os.path.join(t_d, case+".pst"))
    par = pst.parameter_data
    add_par = par.loc[par.parnme.str.contains("_add_"), "parnme"]
    if case == "constr":
        par.loc[add_par, "parubnd"] = constr_bnd
        par.loc[add_par, "parlbnd"] = -constr_bnd

    pst.observation_data.loc[pst.nnz_obs_names,"weight"] = std_weight
    pst.pestpp_options["opt_std_weights"] = True
    pst.pestpp_options["mou_generator"] = mou_gen
    pst.pestpp_options["mou_population_size"] = pop_size
    pst.pestpp_options["opt_risk"] = 0.05
    pst.pestpp_options["save_binary"] = True
    pst.control_data.noptmax = noptmax
    pst.write(os.path.join(t_d, case+".pst"))
    m4 = os.path.join("mou_tests", case+"_test_master_riskobj_match_"+mou_gen)
    pyemu.os_utils.start_workers(t_d, exe_path, case+".pst", num_workers, worker_root="mou_tests",
                                 master_dir=m4, verbose=True, port=port)

    # risk obj ,std weights
    t_d = mou_suite_helper.setup_problem(case, additive_chance=True, risk_obj=True)
    pst = pyemu.Pst(os.path.join(t_d, case+".pst"))
    par = pst.parameter_data
    add_par = par.loc[par.parnme.str.contains("_add_"), "parnme"]
    if case == "constr":
        par.loc[add_par, "parubnd"] = constr_bnd
        par.loc[add_par, "parlbnd"] = -constr_bnd

    pst.observation_data.loc[pst.nnz_obs_names,"weight"] = std_weight
    pst.pestpp_options["opt_std_weights"] = True
    pst.pestpp_options["mou_generator"] = mou_gen
    pst.pestpp_options["mou_population_size"] = pop_size
    pst.pestpp_options["opt_risk"] = 0.05
    pst.pestpp_options["save_binary"] = True
    pst.control_data.noptmax = noptmax * 3
    pst.write(os.path.join(t_d, case+".pst"))
    m5 = os.path.join("mou_tests", case+"_test_master_riskobj_more_"+mou_gen)
    pyemu.os_utils.start_workers(t_d, exe_path, case+".pst", num_workers, worker_root="mou_tests",
                                 master_dir=m5, verbose=True, port=port)


    # risk obj, stack, every
    t_d = mou_suite_helper.setup_problem(case, additive_chance=True, risk_obj=True)
    pst = pyemu.Pst(os.path.join(t_d, case + ".pst"))
    par = pst.parameter_data
    add_par = par.loc[par.parnme.str.contains("_add_"), "parnme"]
    if case == "constr":
        par.loc[add_par, "parubnd"] = constr_bnd
        par.loc[add_par, "parlbnd"] = -constr_bnd

    pst.observation_data.loc[pst.nnz_obs_names, "weight"] = std_weight
    pst.pestpp_options["opt_std_weights"] = True
    pst.pestpp_options["mou_generator"] = mou_gen
    pst.pestpp_options["mou_population_size"] = pop_size
    pst.pestpp_options["opt_risk"] = 0.05
    pst.pestpp_options["save_binary"] = True
    pst.pestpp_options["opt_chance_points"] = "all"
    pst.pestpp_options["opt_recalc_chance_every"] = 1
    pst.pestpp_options["opt_stack_size"] = 100
    pst.control_data.noptmax = noptmax * 3
    pst.write(os.path.join(t_d, case + ".pst"))
    m5 = os.path.join("mou_tests", case + "_test_master_allpts_every_riskobj_more_" + mou_gen)
    pyemu.os_utils.start_workers(t_d, exe_path, case + ".pst", num_workers, worker_root="mou_tests",
                                 master_dir=m5, verbose=True, port=port)


def plot_risk_demo_multi(case = "zdt1", mou_gen="de"):
    import matplotlib.pyplot as plt
    m_deter = os.path.join("mou_tests",case+"_test_master_deter_"+mou_gen)
    m_ravr = os.path.join("mou_tests",case+"_test_master_95_"+mou_gen)
    m_rtol = os.path.join("mou_tests", case+"_test_master_05_"+mou_gen)
    m_robj = os.path.join("mou_tests",case+"_test_master_riskobj_match_"+mou_gen)
    m_robjm = os.path.join("mou_tests", case+"_test_master_riskobj_more_"+mou_gen)

    fig, ax = plt.subplots(1,1,figsize=(10,10))
    for d,c in zip([m_deter,m_ravr,m_rtol,m_robjm],['g','b','r',None,None]):

        pst = pyemu.Pst(os.path.join(d,case+".pst"))
        df = pd.read_csv(os.path.join(d,case+".pareto.archive.summary.csv"))
        mxgen = df.generation.max()
        
        print(d,mxgen)
        df = df.loc[df.generation==mxgen,:]
        print(d,mxgen)
        if "riskobj" in d:
            #print(df.head().loc[:,['obj_1',"obj_2",'_risk_']])
            #ax.scatter(df.obj_1.values[:2],df.obj_2.values[:2],marker="+",c=df._risk_[:2],cmap='jet')
            rdf = df#.loc[df._risk_ < 0.05,:]
            ax.scatter(rdf.obj_1,rdf.obj_2,marker="o",c=1 - rdf._risk_.values,cmap='jet')
        else:
            ax.scatter(df.obj_1.values,df.obj_2.values,marker='.',color=c)

    if case == "zdt1":
        x0 = np.linspace(0,1,1000)
        o1,o2 = [],[]
        for xx0 in x0:
            x = np.zeros(30)
            x[0] = xx0
            ret_vals = mou_suite_helper.zdt1(x)
            o1.append(ret_vals[0][0])
            o2.append(ret_vals[0][1])

        ax.plot(o1,o2,"k",label="truth")
    plt.show()


def plot_risk_demo_multi_3pane(case="zdt1",mou_gen="de"):
    import matplotlib.pyplot as plt
    #m_d = os.path.join("mou_tests", case + "_test_master_riskobj_more_"+mou_gen)
    m_d = os.path.join("mou_tests", case + "_test_master_95_" + mou_gen)
    pst = pyemu.Pst(os.path.join(m_d, case + ".pst"))
    obj_names = pst.pestpp_options["mou_objectives"].lower().split(',')
    df = pd.read_csv(os.path.join(m_d, case + ".pareto.archive.summary.csv"))
    mxgen = df.generation.max()
    print(mxgen)
    df = df.loc[df.generation == mxgen, :]
    df = df.loc[df.nsga2_front==1,:]

    #df = df.loc[df.obj_2<1,:]
    
    #df = pyemu.ObservationEnsemble.from_binary(pst=pst,filename=os.path.join(m_d,case+".archive.obs_pop.jcb"))
    #pe = pyemu.ParameterEnsemble.from_binary(pst=pst,filename=os.path.join(m_d,case+".archive.dv_pop.jcb"))
    #print(pe.loc[:,"dv_1"].min())
    #return
    #df.loc[:,"_risk_"] = pe.loc[df.index,"_risk_"].values
    #df = df.loc[df._risk_.apply(lambda x: x > 0.95),:]
    rvals = None
    if "_risk_" in obj_names:
        rvals = df.loc[:,"_risk_"].values.copy()

    fig, axes = plt.subplots(len(obj_names), len(obj_names), figsize=(10, 10))
    for i,o1 in enumerate(obj_names):
        for j in range(i+1):
            o2 = obj_names[j]
            ax = axes[i, j]
            v1 = df.loc[:, o1].values
            v2 = df.loc[:, o2].values
            if i == j:
                ax.hist(v2,bins=20,facecolor="0.5",edgecolor="none",alpha=0.5)
                ax.set_xlabel(o1)
                ax.set_yticks([])
            else:
                if rvals is None:
                    ax.scatter(v2,v1,marker=".",color="0.5",s=20)
                else:
                    ax.scatter(v2,v1,marker=".",c=rvals,s=20)
                ax.set_xlabel(o2)
                ax.set_ylabel(o1)

            if case == "zdt1" and o1 == "obj_2" and o2 == "obj_1":
                x0 = np.linspace(0, 1, 1000)
                ov1, ov2 = [], []
                for xx0 in x0:
                    x = np.zeros(30)
                    x[0] = xx0
                    ret_vals = mou_suite_helper.zdt1(x)
                    ov1.append(ret_vals[0][0])
                    ov2.append(ret_vals[0][1])

                ax.plot(ov1, ov2, "k", label="truth")

        for j in range(i+1,len(obj_names)):
            o2 = obj_names[j]
            ax = axes[i, j]
            v1 = df.loc[:, o1].values
            v2 = df.loc[:, o2].values
            if rvals is None:
                ax.scatter(v2,v1,marker=".",color="0.5",s=20)
            else:
                ax.scatter(v2,v1,marker=".",c=rvals,s=20)
            ax.set_xlabel(o2)
            ax.set_ylabel(o1)


            if case == "zdt1" and o1 == "obj_1" and o2 == "obj_2":
                x0 = np.linspace(0, 1, 1000)
                ov1, ov2 = [], []
                for xx0 in x0:
                    x = np.zeros(30)
                    x[0] = xx0
                    ret_vals = mou_suite_helper.zdt1(x)
                    ov1.append(ret_vals[0][0])
                    ov2.append(ret_vals[0][1])

                ax.plot(ov2, ov1, "k", label="truth")
    plt.show()


def plot_risk_demo_rosen():
    case = "rosenc"
    import matplotlib.pyplot as plt
    mou_gen = "de"
    m_deter = os.path.join("mou_tests",case+"_test_master_deter_"+mou_gen)
    m_ravr = os.path.join("mou_tests",case+"_test_master_95_"+mou_gen)
    m_rtol = os.path.join("mou_tests", case+"_test_master_05_"+mou_gen)
    #m_robj = os.path.join("mou_tests",case+"_test_master_riskobj_match")
    m_robjm = os.path.join("mou_tests", case+"_test_master_riskobj_match_"+mou_gen)
    bins = 30#np.linspace(-5,5,30)
    fig, axes = plt.subplots(1,2,figsize=(8,4))
    #axes = axes.flatten()
    for d,c,label in zip([m_deter,m_ravr,m_rtol,m_robjm],['g','b','r',"c","m"],["risk neutral (risk=0.5)","risk averse (risk=0.95)","risk tolerant (risk=0.05)","risk as an objective"]):

        pst = pyemu.Pst(os.path.join(d,case+".pst"))
        df = pd.read_csv(os.path.join(d,case+".pareto.archive.summary.csv"))
        df2 = pd.read_csv(os.path.join(d[:-2] + "pso",case+".pareto.archive.summary.csv"))
        mxgen = df.generation.max()
        #mxgen = 10
        #print(d,mxgen)
        df = df.loc[df.generation==mxgen,:]
        mxgen2 = df2.generation.max()
        df2 = df2.loc[df2.generation==mxgen2,:]
        print(d,mxgen)
        ax = axes[0]

        if "riskobj" in d:
            #print(df.head().loc[:,['obj_1',"obj_2",'_risk_']])
            #ax.scatter(df.obj_1.values[:2],df.obj_2.values[:2],marker="+",c=df._risk_[:2],cmap='jet')
            rdf = df#.loc[df._risk_ < 0.05,:]
            #axt = plt.twinx(ax)
            ax.scatter(rdf.obj_1,rdf._risk_,marker="o",c=1 - rdf._risk_.values,cmap='jet',label="risk vs objective function\npareto frontier")
            ax.set_ylabel("risk")


        #ax.hist(df.obj_1.values,bins=bins,facecolor="0.5",edgecolor="none",alpha=0.5,density=False)

        #ax.hist(df.obj_1.values, bins=bins, facecolor=0.5, density=False)

        ax.plot([df.obj_1.mean(),df.obj_1.mean()],[0,1],color=c,ls="--",label=label)
        #ax.set_xlim(bins.min(),bins.max())
        ax.set_xlim(-5,5)
        ax.set_title("A) DE",loc="left")
        ax.set_xlabel("objective function value")

        ax = axes[1]

        if "riskobj" in d:
            # print(df.head().loc[:,['obj_1',"obj_2",'_risk_']])
            # ax.scatter(df.obj_1.values[:2],df.obj_2.values[:2],marker="+",c=df._risk_[:2],cmap='jet')
            rdf = df2 # .loc[df._risk_ < 0.05,:]
            # axt = plt.twinx(ax)
            ax.scatter(rdf.obj_1, rdf._risk_, marker="o", c=1 - rdf._risk_.values, cmap='jet',
                       label="risk vs objective function\npareto frontier")
            ax.set_ylabel("risk")

        # ax.hist(df.obj_1.values,bins=bins,facecolor="0.5",edgecolor="none",alpha=0.5,density=False)

        # ax.hist(df.obj_1.values, bins=bins, facecolor=0.5, density=False)

        ax.plot([df2.obj_1.mean(), df2.obj_1.mean()], [0, 1], color=c, ls="--", label=label)
        # ax.set_xlim(bins.min(),bins.max())
        ax.set_xlim(-5, 5)
        ax.set_title("B) PSO", loc="left")
        ax.set_xlabel("objective function value")
        #ax.set_yticks([])

    for ax in axes:
        ax.set_ylim(0,1)
        ax.legend(loc="upper left")

    #plt.show()
    plt.savefig("rosenc_risk_demo.pdf")

def risk_obj_test():
    t_d = mou_suite_helper.setup_problem("constr", additive_chance=True, risk_obj=True)

    pst = pyemu.Pst(os.path.join(t_d, "constr.pst"))
    pst.pestpp_options["opt_std_weights"] = True
    pst.pestpp_options["opt_recalc_chance_every"] = 5
    pst.pestpp_options["mou_generator"] = "de"
    pst.pestpp_options["mou_population_size"] = 5
    pst.pestpp_options["opt_risk"] = 0.05
    pst.pestpp_options["mou_risk_objective"] = True
    pst.observation_data.loc["obj_1","weight"] = 0.5
    pst.observation_data.loc["obj_2", "weight"] = 0.5
    pst.control_data.noptmax = -1
    pst.write(os.path.join(t_d, "constr.pst"))
    m1 = os.path.join("mou_tests", "constr_riskobj_test_master")
    pyemu.os_utils.start_workers(t_d, exe_path, "constr.pst", 5, worker_root="mou_tests",
                                 master_dir=m1, verbose=True, port=port)
    return

    pst = pyemu.Pst(os.path.join(t_d, "constr.pst"))
    pst.pestpp_options["opt_std_weights"] = False
    pst.pestpp_options["opt_chance_points"] = "all"
    pst.pestpp_options["opt_recalc_chance_every"] = 5
    pst.pestpp_options["opt_stack_size"] = 5
    pst.pestpp_options["mou_generator"] = "de"
    pst.pestpp_options["mou_population_size"] = 5
    pst.pestpp_options["opt_risk"] = 0.05
    pst.pestpp_options["mou_risk_objective"] = True
    pst.control_data.noptmax = -1
    pst.write(os.path.join(t_d, "constr.pst"))
    m1 = os.path.join("mou_tests", "constr_riskobj_test_master")
    pyemu.os_utils.start_workers(t_d, exe_path, "constr.pst", 5, worker_root="mou_tests",
                                 master_dir=m1, verbose=True,port=port)
  
    dv = pd.read_csv(os.path.join(m1,"constr.0.dv_pop.csv"),index_col=0)
    dv.loc[:,"_risk_"] = pst.pestpp_options["opt_risk"]
    dv.to_csv(os.path.join(t_d,"dv.csv"))
    pst.pestpp_options["mou_dv_population_file"] = "dv.csv"
    pst.write(os.path.join(t_d, "constr.pst"))
    m1 = os.path.join("mou_tests", "constr_riskobj_test_master")
    pyemu.os_utils.start_workers(t_d, exe_path, "constr.pst", 5, worker_root="mou_tests",
                                 master_dir=m1, verbose=True, port=port)
    pst.pestpp_options["mou_risk_objective"] = False
    pst.write(os.path.join(t_d, "constr.pst"))
    m2 = os.path.join("mou_tests", "constr_riskobj_test_master2")
    pyemu.os_utils.start_workers(t_d, exe_path, "constr.pst", 5, worker_root="mou_tests",
                                 master_dir=m2, verbose=True, port=port)

    op1 = pd.read_csv(os.path.join(m1, "constr.0.obs_pop.csv"), index_col=0)
    op2 = pd.read_csv(os.path.join(m2, "constr.0.obs_pop.csv"), index_col=0)
    d = (op1 - op2).apply(lambda x: np.abs(x))
    print(d.max())
    assert d.max().max() < 1.0e-10,d.max().max()
    op1 = pd.read_csv(os.path.join(m1, "constr.0.chance.obs_pop.csv"), index_col=0)
    op2 = pd.read_csv(os.path.join(m2, "constr.0.chance.obs_pop.csv"), index_col=0)
    d = (op1 - op2).apply(lambda x: np.abs(x))
    print(d.max())
    assert d.max().max() < 1.0e-10, d.max().max()

    shutil.copy2(os.path.join(m2,"constr.0.nested.par_stack.csv"),os.path.join(t_d,"par_stack.csv"))
    shutil.copy2(os.path.join(m2, "constr.0.nested.obs_stack.csv"), os.path.join(t_d, "obs_stack.csv"))
    pst.pestpp_options["opt_par_stack"] = "par_stack.csv"
    pst.pestpp_options["opt_obs_stack"] = "obs_stack.csv"
    dv = pyemu.ParameterEnsemble.from_uniform_draw(pst,5)
    dv.loc[:,"_risk_"] = pst.pestpp_options["opt_risk"]
    pst.pestpp_options["mou_risk_objective"] = True
    dv.to_csv(os.path.join(t_d, "dv.csv"))
    pst.write(os.path.join(t_d, "constr.pst"))
    m1 = os.path.join("mou_tests", "constr_riskobj_test_master")
    pyemu.os_utils.start_workers(t_d, exe_path, "constr.pst", 5, worker_root="mou_tests",
                                 master_dir=m1, verbose=True, port=port)
    pst.pestpp_options["mou_risk_objective"] = False
    pst.write(os.path.join(t_d, "constr.pst"))
    m2 = os.path.join("mou_tests", "constr_riskobj_test_master2")
    pyemu.os_utils.start_workers(t_d, exe_path, "constr.pst", 5, worker_root="mou_tests",
                                 master_dir=m2, verbose=True, port=port)

    op1 = pd.read_csv(os.path.join(m1, "constr.0.obs_pop.csv"), index_col=0)
    op2 = pd.read_csv(os.path.join(m2, "constr.0.obs_pop.csv"), index_col=0)
    d = (op1 - op2).apply(lambda x: np.abs(x))
    print(d.max())
    assert d.max().max() < 1.0e-10, d.max().max()
    op1 = pd.read_csv(os.path.join(m1, "constr.0.chance.obs_pop.csv"), index_col=0)
    op2 = pd.read_csv(os.path.join(m2, "constr.0.chance.obs_pop.csv"), index_col=0)
    d = (op1 - op2).apply(lambda x: np.abs(x))
    print(d.max())
    assert d.max().max() < 1.0e-10, d.max().max()


def basic_pso_test(case="zdt1"):

    t_d = mou_suite_helper.setup_problem(case, additive_chance=True, risk_obj=False)
    pst = pyemu.Pst(os.path.join(t_d, case+".pst"))
    pst.pestpp_options["mou_generator"] = "pso"
    pst.pestpp_options["opt_risk"] = 0.95
    pst.control_data.noptmax = 20
    pst.pestpp_options["mou_population_size"] = 20
    pst.pestpp_options["mou_save_population_every"] = 1
    pst.write(os.path.join(t_d, case+".pst"))
    m_d = os.path.join("mou_tests", case+"_pso_master_risk")
    pyemu.os_utils.start_workers(t_d, exe_path,  case+".pst", 20, worker_root="mou_tests",
                                master_dir=m_d, verbose=True, port=port)
    assert os.path.exists(os.path.join(m_d,"{0}.{1}.fosm.jcb".format(case,pst.control_data.noptmax)))

    for i in range(0,pst.control_data.noptmax+1):
        dv_file = os.path.join(m_d,"{0}.{1}.archive.dv_pop.csv".format(case,i))
        oe_file = os.path.join(m_d,"{0}.{1}.archive.obs_pop.csv".format(case,i))
        assert os.path.exists(dv_file)
        assert os.path.exists(oe_file)
        dv_pso = pd.read_csv(dv_file)
        oe_pso = pd.read_csv(oe_file)
        assert dv_pso.shape[0] == oe_pso.shape[0]

    pst.pestpp_options["mou_generator"] = "de"
    pst.pestpp_options["opt_risk"] = 0.95
    pst.pestpp_options["mou_save_population_every"] = 1
    pst.write(os.path.join(t_d, case+".pst"))
    m_d = os.path.join("mou_tests", case+"_de_master_risk")
    #pyemu.os_utils.start_workers(t_d, exe_path, case+".pst", 20, worker_root="mou_tests",
    #                             master_dir=m_d, verbose=True, port=port)
    # assert os.path.exists(os.path.join(m_d, "{0}.{1}.fosm.jcb".format(case, pst.control_data.noptmax)))

    # for i in range(0,pst.control_data.noptmax+1):
    #     dv_file = os.path.join(m_d,"{0}.{1}.dv_pop.csv".format(case,i))
    #     oe_file = os.path.join(m_d,"{0}.{1}.obs_pop.csv".format(case,i))
    #     assert os.path.exists(dv_file)
    #     assert os.path.exists(oe_file)
    #     dv_de = pd.read_csv(dv_file)
    #     oe_de = pd.read_csv(oe_file)
    #     assert dv_de.shape[0] == oe_de.shape[0]

    method = mou_suite_helper.zdt1
    x0 = np.linspace(0,1,10000)
    o1,o2 = [],[]
    for xx0 in x0:
        x = np.zeros(30)
        x[0] = xx0
        ret_vals = method(x)
        o1.append(ret_vals[0][0])
        o2.append(ret_vals[0][1])
    

    o1 = np.array(o1)
    o2 = np.array(o2)
    print(o1.min(),oe_pso.loc[:,"obj_1"].values.min())
    diff = np.abs(o1.min() - oe_pso.loc[:,"obj_1"].values.min()) 
    #print(diff)
    #assert diff < 1.0e-4
    diff = np.abs(o1.max() - oe_pso.loc[:,"obj_1"].values.max()) 
    #print(diff)
    #assert diff < 1.0e-4

    diff = np.abs(o2.min() - oe_pso.loc[:,"obj_2"].values.min()) 
    print(diff)
    #assert diff < 1.0e-4
    
    opt_1,opt_2 = o1.min(),o2.min()
    opt = np.array([opt_1,opt_2])
    truth = np.array([o1,o2]).transpose()
    #print(truth)
    #print(opt)
    dist = [(opt-t).sum()**2 for t in truth]
    knee_idx = np.argmin(dist)
    #print(knee_idx,dist[knee_idx],truth[knee_idx])
    knee = truth[knee_idx]

    knee_dist = [(knee-sol).sum()**2 for sol in oe_pso.loc[:,["obj_1","obj_2"]].values]
    #print(knee_dist)
    knee_sol_idx = np.argmin(knee_dist)
    
    knee_sol = oe_pso.loc[:,["obj_1","obj_2"]].values[knee_sol_idx]
    #print(knee_sol_idx,knee_sol)

    dist = (knee - knee_sol).sum()**2
    print(dist)
    #assert dist < 0.001

    # import matplotlib.pyplot as plt
    # fig,axes = plt.subplots(1,2,figsize=(10,5))
    # axes[0].plot(o1,o2)
    # axes[0].scatter(oe_pso.loc[:,"obj_1"].values,oe_pso.loc[:,"obj_2"])
    # axes[1].plot(o1,o2)
    # axes[1].scatter(oe_de.loc[:,"obj_1"].values,oe_de.loc[:,"obj_2"])
    # plt.show()



def plot_zdt_risk_demo_compare(case="zdt1"):
    import matplotlib.pyplot as plt
    # m_d = os.path.join("mou_tests", case + "_test_master_riskobj_more_"+mou_gen)
    truth_func = None
    numdv = 30
    if case == "zdt1":
        truth_func = mou_suite_helper.zdt1
    elif case == "zdt3":
        truth_func = mou_suite_helper.zdt3
    elif case == "zdt6":
        truth_func = mou_suite_helper.zdt6
        numdv = 10
    elif case == "constr":
        numdv = 2
    #else:
    #    raise Exception()

    x0 = np.linspace(-1, 1, 1000)
    ov1, ov2 = [], []
    if truth_func is not None:
        for xx0 in x0:
            x = np.zeros(numdv)
            x[0] = xx0
            ret_vals = truth_func(x)
            ov1.append(ret_vals[0][0])
            ov2.append(ret_vals[0][1])
    master_ds = [case+"_test_master_deter",case+"_test_master_05",case+"_test_master_95"]
    fig, axes = plt.subplots(2, 2, figsize=(8, 8))
    axes = axes.flatten()
    obj_names = ["obj_1","obj_2"]
    colors = ["g","r","b"]
    labels = ["risk neutral (risk=0.5)","risk tolerant (risk=0.05)","risk averse (risk=0.95)"]
    for m_d,c,label in zip(master_ds,colors,labels):
        m_d = os.path.join("mou_tests",m_d)
        pst = pyemu.Pst(os.path.join(m_d+"_de", case + ".pst"))

        df_de = pd.read_csv(os.path.join(m_d+"_de", case + ".pareto.archive.summary.csv"))
        mxgen = df_de.generation.max()
        print(mxgen)
        df_de = df_de.loc[df_de.generation == mxgen, :]
        df_de = df_de.loc[df_de.nsga2_front == 1, :]

        df_pso = pd.read_csv(os.path.join(m_d + "_pso", case + ".pareto.archive.summary.csv"))
        mxgen = df_pso.generation.max()
        print(mxgen)
        df_pso = df_pso.loc[df_pso.generation == mxgen, :]
        df_pso = df_pso.loc[df_pso.nsga2_front == 1, :]


        axes[0].scatter(df_de.obj_1.values,df_de.obj_2.values,color=c,label=label)
        axes[1].scatter(df_pso.obj_1.values, df_pso.obj_2.values, color=c, label=label)

        if label == labels[0]:


            axes[0].plot(ov1, ov2, "k", label="truth")
            axes[1].plot(ov1, ov2, "k", label="truth")
            #print(x0)
        axes[0].set_title("A) DE specified risk",loc="left")
        axes[1].set_title("A) PSO specified risk", loc="left")
        axes[0].set_xlabel("objective 1")
        axes[1].set_xlabel("objective 1")
        axes[0].set_ylabel("objective 2")
        axes[1].set_ylabel("objective 2")
        axes[1].legend(loc="upper right")

    m_d = os.path.join("mou_tests",case + "_test_master_riskobj_more")
    df_de = pd.read_csv(os.path.join(m_d + "_de", case + ".pareto.archive.summary.csv"))
    mxgen = df_de.generation.max()
    print(mxgen)
    df_de = df_de.loc[df_de.generation == mxgen, :]
    df_de = df_de.loc[df_de.nsga2_front == 1, :]

    df_pso = pd.read_csv(os.path.join(m_d + "_pso", case + ".pareto.archive.summary.csv"))
    mxgen = df_pso.generation.max()
    print(mxgen)
    df_pso = df_pso.loc[df_pso.generation == mxgen, :]
    df_pso = df_pso.loc[df_pso.nsga2_front == 1, :]

    ax = axes[2]
    ax.scatter(df_de.obj_1,df_de.obj_2,marker='.',c=df_de._risk_)
    ax.plot(ov1, ov2, "k", label="truth")
    ax.set_title("C) DE risk objective", loc="left")
    ax.set_ylim(0,10)

    ax = axes[3]
    ax.scatter(df_pso.obj_1, df_pso.obj_2, marker='.', c=df_pso._risk_)
    ax.plot(ov1, ov2, "k", label="truth")
    ax.set_title("D) PSO risk objective", loc="left")
    ax.set_ylim(0,10)
    plt.show()

def zdt1_invest():
    m_d = os.path.join("mou_tests","zdt1_test_master_riskobj_more_pso")
    df_sum = pd.read_csv(os.path.join(m_d,"zdt1.pareto.archive.summary.csv"))
    df_sum = df_sum.loc[df_sum.generation==df_sum.generation.max(),:]
    df_sum = df_sum.loc[df_sum.nsga2_front==1,:]
    df_sum.loc[:,"member"] = df_sum.member.str.lower()
    pst = pyemu.Pst(os.path.join(m_d,"zdt1.pst"))

    pe = pyemu.ParameterEnsemble.from_binary(pst=pst,filename=os.path.join(m_d,"zdt1.archive.dv_pop.jcb"))
    oe = pyemu.ObservationEnsemble.from_binary(pst=pst,filename=os.path.join(m_d,"zdt1.archive.obs_pop.jcb"))
    print(oe.head())
    print(pe.loc[oe.head().index,"_risk_"].head())

    #print(pe.index.shape)
    #print(df_sum.member.shape)

    #print(set(pe.index.tolist()) - set(df_sum.member.tolist()))

    #df_sum.sort_values(by="obj_2",inplace=True,ascending=False)
    #print(df_sum.loc[:,["member","obj_1","obj_2"]].head())
    mem = "gen=144_member=7221_pso"
    pvals = pe.loc[mem,:]
    #print(pvals)
    ovals,cvals = mou_suite_helper.zdt1(pvals.values[:29])
    #print(ovals)
    #print(oe.loc["gen=150_member=7501_pso",:])

    #pst.parameter_data.loc[pvals.index,"parval1"] = pvals.values
    # pst.control_data.noptmax = 1
    # pst.pestpp_options["opt_risk"] = 0.5
    # pst.pestpp_options.pop("mou_risk_objective")
    # pst.pestpp_options["mou_dv_population_file"] = "zdt1.dv_pop.jcb"
    # pst.pestpp_options["mou_obs_population_restart_file"] = "zdt1.obs_pop.jcb"
    # pst.write(os.path.join(m_d,"test.pst"))
    # return
    #pyemu.os_utils.run("pestpp-mou test.pst",cwd=m_d)

    #check that mem is nondominated
    # all objs are minimize
    objs = ["obj_1","obj_2","_risk_"]
    df_sum.loc[:, "_risk_"] *= -1
    mem_vals = df_sum.loc[df_sum.member==mem,objs].values[0]

    print(mem_vals)
    for idx in df_sum.index[1:]:
        mem2 = df_sum.loc[idx,"member"]
        if mem == mem2:
            continue
        mem2_vals = df_sum.loc[idx,objs].values
        d = mem_vals - mem2_vals
        if np.all(d<0):
            print(mem2,mem_vals,mem2_vals)
            #print(d)
        #break


def water_invest():
    case = "water"
    t_d = mou_suite_helper.setup_problem(case, additive_chance=False, risk_obj=False)
    pst = pyemu.Pst(os.path.join(t_d, case + ".pst"))
    pst.pestpp_options["mou_generator"] = "pso"
    pst.pestpp_options["mou_population_size"] = 100
    pst.pestpp_options["opt_risk"] = 0.5
    pst.control_data.noptmax = 200
    pst.write(os.path.join(t_d, case + ".pst"))
    m_d = os.path.join("mou_tests", case + "_pso_master_risk")
    pyemu.os_utils.start_workers(t_d, exe_path, case + ".pst", 50, worker_root="mou_tests",
                                 master_dir=m_d, verbose=True, port=port)

    pst.pestpp_options["mou_generator"] = "de"
    pst.pestpp_options["opt_risk"] = 0.5
    pst.pestpp_options["mou_population_size"] = 100
    pst.control_data.noptmax = 200
    pst.write(os.path.join(t_d, "zdt1.pst"))
    m_d = os.path.join("mou_tests", case + "_de_master_risk")
    pyemu.os_utils.start_workers(t_d, exe_path, case + ".pst", 50, worker_root="mou_tests",
                                 master_dir=m_d, verbose=True, port=port)


def plot_constr_risk():
    import matplotlib.pyplot as plt
    case = "constr"
    master_ds = [case + "_test_master_deter",
                 case + "_test_master_05", case + "_test_master_95"]
    fig, axes = plt.subplots(2, 2, figsize=(8, 8))
    axes = axes.flatten()
    for i,ax in enumerate(axes):
        label = False
        if i == 1:
            label = True
        get_constr_base_plot(ax,label=label)
        ax.set_xlabel("objective 1 (minimize)")
        ax.set_ylabel("objective 2 (minimize)")
    obj_names = ["obj_1", "obj_2"]
    cmap = plt.get_cmap("viridis")
    colors = [cmap(r) for r in [0.5,0.05,0.95]]
    labels = ["risk neutral (risk=0.5)", "risk tolerant (risk=0.05)", "risk averse (risk=0.95)"]
    for m_d, c, label in zip(master_ds, colors, labels):
        m_d = os.path.join("mou_tests", m_d)
        pst = pyemu.Pst(os.path.join(m_d + "_de", case + ".pst"))
        df_de = pd.read_csv(os.path.join(m_d + "_de", case + ".pareto.archive.summary.csv"))
        mxgen = df_de.generation.max()
        print(mxgen)
        df_de = df_de.loc[df_de.generation == mxgen, :]
        df_de = df_de.loc[df_de.nsga2_front == 1, :]

        df_pso = pd.read_csv(os.path.join(m_d + "_pso", case + ".pareto.archive.summary.csv"))
        mxgen = df_pso.generation.max()
        print(mxgen)
        df_pso = df_pso.loc[df_pso.generation == mxgen, :]
        df_pso = df_pso.loc[df_pso.nsga2_front == 1, :]
        axes[0].scatter(df_de.obj_1.values, df_de.obj_2.values, color=c, s=4, label=label,zorder=10,alpha=0.5)
        axes[1].scatter(df_pso.obj_1.values, df_pso.obj_2.values, color=c, s=4, label=label,zorder=10,alpha=0.5)

        axes[0].set_title("A) DE specified risk", loc="left")
        axes[1].set_title("B) PSO specified risk", loc="left")
        axes[1].legend(loc="upper right",framealpha=1.0)

    m_d = os.path.join("mou_tests", case + "_test_master_riskobj_more")
    df_de = pd.read_csv(os.path.join(m_d + "_de", case + ".pareto.archive.summary.csv"))
    mxgen = df_de.generation.max()
    print(mxgen)
    df_de = df_de.loc[df_de.generation == mxgen, :]
    df_de = df_de.loc[df_de.nsga2_front == 1, :]

    df_pso = pd.read_csv(os.path.join(m_d + "_pso", case + ".pareto.archive.summary.csv"))
    mxgen = df_pso.generation.max()
    print(mxgen)
    df_pso = df_pso.loc[df_pso.generation == mxgen, :]
    df_pso = df_pso.loc[df_pso.nsga2_front == 1, :]

    ax = axes[2]
    ax.scatter(df_de.obj_1, df_de.obj_2, marker='.', c=df_de._risk_,alpha=0.5,s=10)

    ax.set_title("C) DE objective risk", loc="left")
    ax.set_ylim(0, 10)

    ax = axes[3]
    ax.scatter(df_pso.obj_1, df_pso.obj_2, marker='.', c=df_pso._risk_,alpha=0.5,s=10)

    ax.set_title("D) PSO objective risk", loc="left")
    ax.set_ylim(0, 10)
    plt.savefig("constr_results.pdf")

def get_constr_base_plot(ax,label=False,fontsize=10):

    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon

    def f(x1, x2):
        return (1 + x2) / x1

    # first front
    x1 = np.linspace(0.365, 0.665, 100)
    y1 = 6 - 9 * x1
    f1 = f(x1, y1)

    # second front
    x2 = np.linspace(0.1, 1.0, 100)
    f2 = f(x2, 0)

    # top constraint
    y2 = (9 * x1) - 1
    f3 = f(x1, y2)

    # top feas
    x3 = np.linspace(0.67, 1.0, 100)
    f4 = f(x3, 5)

    # right feas
    y4 = np.linspace(0, 5, 100)
    x4 = np.ones_like(y4)
    f5 = f(1, y4)

    # make the polygon

    # fig,ax = plt.subplots(1,1,figsize=(6,6))
    # ax = axes[0]
    if label:
        ax.plot(x1, f1, "k--", label="constraint")
        ax.plot(x2, f2, "r--", label="pareto frontier")
    else:
        ax.plot(x1, f1, "k--")
        ax.plot(x2, f2, "r--")
    ax.plot(x2, f3, "k--")
    # ax.plot(x3,f4,"")
    # ax.plot(x4,f5,"m")
    ax.set_ylim(0, 10)
    ax.set_xlim(0, 1.1)
    # plt.show()

    # get feasible polygon
    x1 = np.linspace(0.39, 0.67, 100)
    y1 = 6 - 9 * x1
    f1 = f(x1, y1)
    x2 = np.linspace(0.39, 0.66, 100)
    x11 = np.linspace(0.0, 1.0, 100)
    y2 = (9 * x1) - 1
    f3 = f(x1, y2)
    x3 = np.linspace(0.665, 1.0, 100)
    f4 = f(x3, 5)
    f2 = f(x3, 0)
    y4 = np.linspace(0, 5, 100)
    x4 = np.ones_like(y4)
    f5 = f(1, y4)
    xvals = list(x3)
    xvals.extend(list(x4))
    xvals.extend(list(np.flipud(x3)))
    xvals.extend(list(np.flipud(x2)))
    xvals.extend(list(x1))
    yvals = list(f2)
    yvals.extend(list(f5))
    yvals.extend(list(np.flipud(f4)))
    yvals.extend(list(np.flipud(f3)))
    yvals.extend(list(f1))
    xy = np.array([xvals, yvals]).transpose()
    if label:
        p = Polygon(xy, facecolor="0.5", alpha=0.5, edgecolor="none",zorder=0, label="feasible region")
    else:
        p = Polygon(xy, facecolor="0.5", alpha=0.5, edgecolor="none", zorder=0)
    ax.add_patch(p)

    ax.set_ylim(0, 10)
    ax.set_xlim(0, 1.1)
    #ax.set_ylabel("objective 1")
    #ax.set_xlabel("objective 2")
    #ax.legend(loc="lower left")


def plot_constr_risk_2():
    import matplotlib.pyplot as plt
    import string
    case = "constr"
    master_ds = [case + "_test_master_riskobj_singlept_more",
                 case + "_test_master_riskobj_allpts_more"]
    labels = ["single chance point","all chance point"]
    fig, axes = plt.subplots(2, 2, figsize=(8, 8))
    ax_count = 0
    for i,ax in enumerate(axes.flatten()):
        label = False;
        if i == 1:
            label = True
        get_constr_base_plot(ax,label=label)
        ax.set_xlabel("objective 1 (minimize)")
        ax.set_ylabel("objective 2 (minimize)")
    for irow,(m_d,lab) in enumerate(zip(master_ds,labels)):
        m_d = os.path.join("mou_tests", m_d)
        pst = pyemu.Pst(os.path.join(m_d + "_de", case + ".pst"))
        df_de = pd.read_csv(os.path.join(m_d + "_de", case + ".pareto.archive.summary.csv"))
        mxgen = df_de.generation.max()
        print(mxgen)
        df_de = df_de.loc[df_de.generation == mxgen, :]
        df_de = df_de.loc[df_de.nsga2_front == 1, :]

        df_pso = pd.read_csv(os.path.join(m_d + "_pso", case + ".pareto.archive.summary.csv"))
        mxgen = df_pso.generation.max()
        print(mxgen)
        df_pso = df_pso.loc[df_pso.generation == mxgen, :]
        df_pso = df_pso.loc[df_pso.nsga2_front == 1, :]
        axes[irow,0].scatter(df_de.obj_1.values, df_de.obj_2.values, c=df_de._risk_, s=4, label=label,zorder=10,alpha=0.5)
        axes[irow,1].scatter(df_pso.obj_1.values, df_pso.obj_2.values, c=df_pso._risk_, s=4, label=label,zorder=10,alpha=0.5)

        axes[irow,0].set_title("{0}) DE {1}".format(string.ascii_uppercase[ax_count],lab), loc="left")
        ax_count += 1
        axes[irow,1].set_title("{0}) PSO {1}".format(string.ascii_uppercase[ax_count],lab), loc="left")
        axes[irow,1].legend(loc="upper right",framealpha=1.0)



    plt.savefig("constr_results_2.pdf")

def plot_constr_risk_3():
    import matplotlib.pyplot as plt
    import string
    case = "constr"
    master_ds = [case + "_test_master_allpts_every_riskobj_more",
                 case + "_test_master_riskobj_allpts_more"]
    labels = ["all chance point, every", "all chance point, reuse"]
    fig, axes = plt.subplots(4, 2, figsize=(8, 11))
    ax_count = 0
    for i,ax in enumerate(axes.flatten()):
        label = False;
        if i == 1:
            label = True
        get_constr_base_plot(ax,label=label)
        ax.set_xlabel("objective 1 (minimize)")
        ax.set_ylabel("objective 2 (minimize)")
    for irow,(m_d,lab) in enumerate(zip(master_ds,labels)):
        m_d = os.path.join("mou_tests", m_d)
        pst = pyemu.Pst(os.path.join(m_d + "_de", case + ".pst"))
        df_de = pd.read_csv(os.path.join(m_d + "_de", case + ".pareto.archive.summary.csv"))
        mxgen = df_de.generation.max()
        print(mxgen)
        df_de = df_de.loc[df_de.generation == mxgen, :]
        df_de = df_de.loc[df_de.nsga2_front == 1, :]

        df_pso = pd.read_csv(os.path.join(m_d + "_pso", case + ".pareto.archive.summary.csv"))
        mxgen = df_pso.generation.max()
        print(mxgen)
        df_pso = df_pso.loc[df_pso.generation == mxgen, :]
        df_pso = df_pso.loc[df_pso.nsga2_front == 1, :]
        axes[irow,0].scatter(df_de.obj_1.values, df_de.obj_2.values, c=df_de._risk_, s=4, label=label,zorder=10,alpha=0.5)
        axes[irow,1].scatter(df_pso.obj_1.values, df_pso.obj_2.values, c=df_pso._risk_, s=4, label=label,zorder=10,alpha=0.5)

        axes[irow,0].set_title("{0}) DE {1}".format(string.ascii_uppercase[ax_count],lab), loc="left")
        ax_count += 1
        axes[irow,1].set_title("{0}) PSO {1}".format(string.ascii_uppercase[ax_count],lab), loc="left")
        axes[irow,1].legend(loc="upper right",framealpha=1.0)


    plt.tight_layout()
    plt.savefig("constr_results_3.pdf")



def plot_constr_risk_pub():
    import string
    import matplotlib.pyplot as plt
    case = "constr"
    fs = 9
    master_ds = [case + "_test_master_deter",
                 case + "_test_master_05", case + "_test_master_95"]
    fig, axes = plt.subplots(4, 2, figsize=(8,11))
    for i,ax in enumerate(axes.flatten()):
        label = False
        if i == 1:
            label = True
        get_constr_base_plot(ax,label=label,fontsize=fs)
        ax.set_xlabel("$f_1$ (minimize)",fontsize=fs)
        ax.set_ylabel("$f_2$ (minimize)",fontsize=fs)
    obj_names = ["obj_1", "obj_2"]
    cmap = plt.get_cmap("viridis")
    colors = [cmap(r) for r in [0.5,0.05,0.95]]
    labels = ["risk neutral\n(risk=0.5)", "risk tolerant\n(risk=0.05)", "risk averse\n(risk=0.95)"]
    
    for m_d, c, label in zip(master_ds, colors, labels):
        m_d = os.path.join("mou_tests", m_d)
        pst = pyemu.Pst(os.path.join(m_d + "_de", case + ".pst"))
        df_de = pd.read_csv(os.path.join(m_d + "_de", case + ".pareto.archive.summary.csv"))
        mxgen = df_de.generation.max()
        print(mxgen)
        df_de = df_de.loc[df_de.generation == mxgen, :]
        df_de = df_de.loc[df_de.nsga2_front == 1, :]

        df_pso = pd.read_csv(os.path.join(m_d + "_pso", case + ".pareto.archive.summary.csv"))
        mxgen = df_pso.generation.max()
        print(mxgen)
        df_pso = df_pso.loc[df_pso.generation == mxgen, :]
        df_pso = df_pso.loc[df_pso.nsga2_front == 1, :]
        axes[0,0].scatter(df_de.obj_1.values, df_de.obj_2.values, color=c, s=4, label=label,zorder=10,alpha=0.5)
        axes[0,1].scatter(df_pso.obj_1.values, df_pso.obj_2.values, color=c, s=4, label=label,zorder=10,alpha=0.5)

        axes[0,0].set_title("A) DE specified risk, specified uncertainty", loc="left",fontsize=fs)
        axes[0,1].set_title("B) PSO specified risk, specified uncertainty", loc="left",fontsize=fs)
        axes[0,1].legend(loc="upper right",framealpha=1.0,fontsize=8)

    master_ds = [case+"_test_master_riskobj_more",
                 case + "_test_master_riskobj_singlept_more",
                 case + "_test_master_riskobj_allpts_more",
                 ]
    labels = ["objective risk, specified uncertainty",
              "objective risk, stack-based uncertainty,\n     single chance point, reused across generations",
              "objective risk, stack-based uncertainty,\n     all chance point, reused across generations"]
    ax_count = 2
    for irow,(m_d,lab) in enumerate(zip(master_ds,labels)):
        irow += 1
        m_d = os.path.join("mou_tests", m_d)
        pst = pyemu.Pst(os.path.join(m_d + "_de", case + ".pst"))
        df_de = pd.read_csv(os.path.join(m_d + "_de", case + ".pareto.archive.summary.csv"))
        mxgen = df_de.generation.max()
        print(mxgen)
        df_de = df_de.loc[df_de.generation == mxgen, :]
        df_de = df_de.loc[df_de.nsga2_front == 1, :]

        df_pso = pd.read_csv(os.path.join(m_d + "_pso", case + ".pareto.archive.summary.csv"))
        mxgen = df_pso.generation.max()
        print(mxgen)
        df_pso = df_pso.loc[df_pso.generation == mxgen, :]
        df_pso = df_pso.loc[df_pso.nsga2_front == 1, :]
        axes[irow,0].scatter(df_de.obj_1.values, df_de.obj_2.values, c=df_de._risk_, s=4, label=label,zorder=10,alpha=0.5)
        axes[irow,1].scatter(df_pso.obj_1.values, df_pso.obj_2.values, c=df_pso._risk_, s=4, label=label,zorder=10,alpha=0.5)

        axes[irow,0].set_title("{0}) DE {1}".format(string.ascii_uppercase[ax_count],lab), loc="left",fontsize=fs)
        ax_count += 1
        axes[irow,1].set_title("{0}) PSO {1}".format(string.ascii_uppercase[ax_count],lab), loc="left",fontsize=fs)
        ax_count += 1
        #axes[irow,1].legend(loc="upper right",framealpha=1.0)


    plt.tight_layout()
    plt.savefig("constr_results_pub.pdf")


def pop_sched_test():
    case = "zdt1"
    t_d = mou_suite_helper.setup_problem(case, additive_chance=True, risk_obj=False)
    pst = pyemu.Pst(os.path.join(t_d, case + ".pst"))
    pst.control_data.noptmax = 3
    pst.write(os.path.join(t_d, case + ".pst"))
    m_d = os.path.join("mou_tests",case+"_master_pop_sched")
    pyemu.os_utils.start_workers(t_d, exe_path, case + ".pst", 50, worker_root="mou_tests",
                                 master_dir=m_d, verbose=True, port=port)


def simplex_invest_1():
    case = "zdt1"
    t_d = mou_suite_helper.setup_problem(case, additive_chance=False, risk_obj=False)
    pst = pyemu.Pst(os.path.join(t_d, case + ".pst"))
    pst.pestpp_options["mou_generator"] = "de,simplex"
    pst.pestpp_options["mou_simplex_reflections"] = 10
    pst.control_data.noptmax = 100

    pst.write(os.path.join(t_d, case + ".pst"))
    m_d = os.path.join("mou_tests",case+"_master_simplex1")

    pyemu.os_utils.start_workers(t_d, exe_path, case + ".pst", 50, worker_root="mou_tests",
                                 master_dir=m_d, verbose=True, port=port)



def zdt1_tied_test():
    t_d = mou_suite_helper.setup_problem("zdt1",False,False)
    #df = pd.read_csv(os.path.join(t_d,"prior.csv"),index_col=0)
    #df.loc[:,"_risk_"] = 0.95
    #print(df.columns)
    #df.to_csv(os.path.join(t_d,"prior.csv"))
    pst = pyemu.Pst(os.path.join(t_d,"zdt1.pst"))
    par = pst.parameter_data
    first = pst.par_names[1]
    others = pst.par_names[2:-2]
    par.loc[others,"partrans"] = "tied"
    par.loc[others,"partied"] = first
    par.loc[others,"pargp"] = "tiedup"
    #par.loc[first,"parval1"] = 1.0
    opars = par.loc[par.parnme.str.startswith("obj"),"parnme"]
    par.loc[opars,"partrans"] = "log"
    par.loc[opars,"parlbnd"] = 0.01
    par.loc[opars,"parubnd"] = 2.
    par.loc[opars,"parval1"] = 0.5

    #pst.pestpp_options["mou_dv_population_file"] = "prior.csv"
    #pst.pestpp_options["opt_chance_points"] = "single"
    pst.pestpp_options["opt_recalc_chance_every"] = 1000
    pst.pestpp_options["opt_dec_var_groups"] = "decvars"
    pst.pestpp_options["mou_save_population_every"] = 1
    #pst.pestpp_options["opt_stack_size"] = 10
    #pst.pestpp_options["opt_par_stack"] = "prior.csv"
    pst.pestpp_options["mou_generator"] = "de"
    pst.pestpp_options["mou_population_size"] = 10
    pst.control_data.noptmax = 0
    pst.write(os.path.join(t_d,"zdt1.pst"))
    pyemu.os_utils.run("{0} {1}".format(exe_path,"zdt1.pst"),cwd=t_d)
    df = pd.read_csv(os.path.join(t_d,"dv.dat"),header=None,names=["dv","val"],sep='\s+')
    df.index = df.dv
    assert (df.loc[others,"val"] - df.loc[first,"val"]).apply(lambda x: np.abs(x)).sum() == 0.0

    m1 = os.path.join("mou_tests","zdt1_tied_test")
    pst.control_data.noptmax = 3
    pst.write(os.path.join(t_d,"zdt1.pst"))
    #pyemu.os_utils.start_workers(t_d,exe_path,"zdt1.pst",10,worker_root="mou_tests",
    #                             master_dir=m1,verbose=True,port=port)
    if os.path.exists(m1):
        shutil.rmtree(m1)
    shutil.copytree(t_d,m1)
    pyemu.os_utils.run("{0} {1}".format(exe_path,"zdt1.pst"),cwd=m1)
    for i in range(pst.control_data.noptmax+1):
        dp = pd.read_csv(os.path.join(m1,"zdt1.{0}.dv_pop.csv").format(i),index_col=0)
        op = pd.read_csv(os.path.join(m1, "zdt1.{0}.obs_pop.csv").format(i), index_col=0)
        assert dp.shape[0] == op.shape[0]
        print(dp.shape[0])
        for ii in range(dp.shape[0]):
            assert dp.index[ii] == op.index[ii]
            ii = dp.index[ii]
            d = np.abs(dp.loc[ii,others].values - dp.loc[ii,first]).sum()
            assert d == 0,d

        dp = pd.read_csv(os.path.join(m1, "zdt1.{0}.archive.dv_pop.csv").format(i), index_col=0)
        op = pd.read_csv(os.path.join(m1, "zdt1.{0}.archive.obs_pop.csv").format(i), index_col=0)
        assert dp.shape[0] == op.shape[0]
        print(dp.shape[0])
        for ii in range(dp.shape[0]):
            assert dp.index[ii] == op.index[ii]
            ii = dp.index[ii]
            d = np.abs(dp.loc[ii, others].values - dp.loc[ii, first]).sum()
            assert d == 0, d

    shutil.copy2(os.path.join(m1,"zdt1.0.dv_pop.csv"),os.path.join(t_d,"dvpop.csv"))
    org_dv = pd.read_csv(os.path.join(t_d,"dvpop.csv"),index_col=0)
    pst.pestpp_options["mou_dv_population_file"] = "dvpop.csv"
    pst.write(os.path.join(t_d,"zdt1_restart.pst"),version=2)
    if os.path.exists(m1):
        shutil.rmtree(m1)
    shutil.copytree(t_d,m1)
    pyemu.os_utils.run("{0} {1}".format(exe_path,"zdt1_restart.pst"),cwd=m1)
    new_dv = pd.read_csv(os.path.join(t_d,"dvpop.csv"),index_col=0)
    diff = np.abs(org_dv.values - new_dv.values)
    print(diff.max())
    assert diff.max() < 1.0e-6

            

    pst.pestpp_options["save_binary"] = True
    pst.write(os.path.join(t_d, "zdt1.pst"))
    m1 = os.path.join("mou_tests", "zdt1_tied_test_bin")
    pyemu.os_utils.start_workers(t_d, exe_path, "zdt1.pst", 10, worker_root="mou_tests",
                                 master_dir=m1, verbose=True, port=port)
    for i in range(pst.control_data.noptmax + 1):
        dp = pyemu.ParameterEnsemble.from_binary(pst=pst, filename=os.path.join(m1, "zdt1.{0}.dv_pop.jcb").format(i))
        op = pyemu.ObservationEnsemble.from_binary(pst=pst, filename=os.path.join(m1, "zdt1.{0}.obs_pop.jcb").format(i))
        assert dp.shape[0] == op.shape[0]
        print(dp.shape[0])
        for ii in range(dp.shape[0]):
            assert dp.index[ii] == op.index[ii]
            ii = dp.index[ii]
            d = np.abs(dp._df.loc[ii, others].values - dp._df.loc[ii, first]).sum()
            assert d == 0, d

        dp = pyemu.ParameterEnsemble.from_binary(pst=pst,
                                                 filename=os.path.join(m1, "zdt1.{0}.archive.dv_pop.jcb").format(i))
        op = pyemu.ObservationEnsemble.from_binary(pst=pst,
                                                   filename=os.path.join(m1, "zdt1.{0}.archive.obs_pop.jcb").format(i))
        assert dp.shape[0] == op.shape[0]
        print(dp.shape[0])
        for ii in range(dp.shape[0]):
            assert dp.index[ii] == op.index[ii]
            ii = dp.index[ii]
            d = np.abs(dp._df.loc[ii, others].values - dp._df.loc[ii, first]).sum()
            #print(dp._df.loc[ii,others].values)
            #print(dp._df.loc[ii, first])
            assert d == 0, d


def zdt1_stack_run_for_animation(mou_gen):
    t_d = mou_suite_helper.setup_problem("zdt1",True,True)
    
    pst = pyemu.Pst(os.path.join(t_d,"zdt1.pst"))
    par = pst.parameter_data
    
    pop_size=10
    pst.pestpp_options["opt_recalc_chance_every"] = 10000
    pst.pestpp_options["mou_save_population_every"] = 1
    pst.pestpp_options["opt_stack_size"] = 100
    #pst.pestpp_options["opt_par_stack"] = "prior.csv"
    pst.pestpp_options["mou_generator"] = mou_gen
    pst.pestpp_options["mou_population_size"] = pop_size
    pst.pestpp_options["mou_verbose_level"] = 4
    pst.pestpp_options["opt_chance_points"] = "single"
    pst.pestpp_options["mou_max_archive_size"] = 100000
    pst.control_data.noptmax = 40
    #pst.write(os.path.join(t_d,"zdt1.pst"))
    #m1 = os.path.join("mou_tests","master_zdt1_test_{0}".format(mou_gen))
    #pyemu.os_utils.start_workers(t_d,exe_path,"zdt1.pst",20,worker_root="mou_tests",
    #                             master_dir=m1,verbose=True,port=port)

    # now robust opt
    t_d = mou_suite_helper.setup_problem("zdt1",True,False)

    pe = pyemu.ParameterEnsemble.from_uniform_draw(pst=pst,num_reals=pop_size)
    pe.to_csv(os.path.join(t_d,"init_pop_and_stack.csv"))
    pst = pyemu.Pst(os.path.join(t_d,"zdt1.pst"))
    par = pst.parameter_data
    par.loc[par.parnme.str.contains("add"),"partrans"] = "fixed"

    pst.pestpp_options["mou_save_population_every"] = 1
    #pst.pestpp_options["opt_par_stack"] = "prior.csv"
    pst.pestpp_options["mou_generator"] = mou_gen
    pst.pestpp_options["mou_population_size"] = pop_size
    pst.pestpp_options["mou_verbose_level"] = 4
    pst.pestpp_options["mou_max_archive_size"] = 100000
    pst.pestpp_options["mou_dv_population_file"] = "init_pop_and_stack.csv"
    pst.pestpp_options["mou_shuffle_fixed_pars"] = True
    pst.control_data.noptmax = 40
    pst.write(os.path.join(t_d,"zdt1.pst"))
    m1 = os.path.join("mou_tests","master_zdt1_test_{0}_ro".format(mou_gen))
    pyemu.os_utils.start_workers(t_d,exe_path,"zdt1.pst",20,worker_root="mou_tests",
                                 master_dir=m1,verbose=True,port=port)


def invest(name="zdt1"):
    t_d = mou_suite_helper.setup_problem(name,False,False)
    pst = pyemu.Pst(os.path.join(t_d,"{0}.pst".format(name)))
    par = pst.parameter_data
    #par_org = par.copy()
    #par.loc[["dv_0","dv_1"],"parval1"] = 9
    #par.loc[["dv_0","dv_1"],"parlbnd"] = 9
    #pst.parameter_data = par_org
    pe = pyemu.ParameterEnsemble.from_uniform_draw(pst=pst,num_reals=30)
    pe.to_csv(os.path.join(t_d,"init.csv"))
    
    pst.pestpp_options["mou_save_population_every"] = 1
    pst.pestpp_options["mou_generator"] = "de"
    pst.pestpp_options["mou_population_size"] = 30
    pst.pestpp_options["mou_verbose_level"] = 4
    pst.pestpp_options["mou_dv_population_file"] = "init.csv"

    #pst.pestpp_options["opt_chance_points"] = "all"
    pst.control_data.noptmax = 30
    pst.write(os.path.join(t_d,"{0}.pst".format(name)))
    m1 = os.path.join("mou_tests","master_{0}_test_{1}".format(name,pst.pestpp_options["mou_generator"]))
    pyemu.os_utils.start_workers(t_d,exe_path,"{0}.pst".format(name),30,worker_root="mou_tests",
                                 master_dir=m1,verbose=True,port=port)


def plot_zdt1(name="zdt1",m_d=None):
    if m_d is None:
        m_d = os.path.join("mou_tests","master_{0}_test".format(name))
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    pst = pyemu.Pst(os.path.join(m_d,"{0}.pst".format(name)))
    df = pd.read_csv(os.path.join(m_d,"{0}.pareto.archive.summary.csv".format(name)))
    gens = df.generation.unique()
    gens.sort()
    plt_d = m_d+"_plots"
    if os.path.exists(plt_d):
        shutil.rmtree(plt_d)
    os.makedirs(plt_d)
    for i,g in enumerate(gens):
        #gdf = df.loc[df.generation==g,:].copy()
        gdf = pd.read_csv(os.path.join(m_d,"{0}.{1}.archive.obs_pop.csv".format(name,g)))
        dvdf = pd.read_csv(os.path.join(m_d,"{0}.{1}.archive.dv_pop.csv".format(name,g)))
        dvdf.index = dvdf.real_name
        fig,ax = plt.subplots(1,1,figsize=(5.5,5))
        c = 'b'
        if "_risk_" in dvdf.columns:
            c = dvdf.loc[gdf.real_name.values,"_risk_"].values
        #print(g,c)
        ax.scatter(gdf.obj_1.values,gdf.obj_2.values,marker=".",s=20,c=c,cmap="jet_r")

        gdf = pd.read_csv(os.path.join(m_d,"{0}.{1}.obs_pop.csv".format(name,g)))
        dvdf = pd.read_csv(os.path.join(m_d,"{0}.{1}.dv_pop.csv".format(name,g)))
        dvdf.index = dvdf.real_name
        ax.scatter(gdf.obj_1.values,gdf.obj_2.values,marker='.',s=20,c='0.5',alpha=0.5)
        #ax.set_xlim(2,-0.75)
        #ax.set_ylim(-1,1)
        ax.set_xlim(50,-0.75)
        ax.set_ylim(-1,50)
        ax.grid()
        cax = fig.colorbar(mpl.cm.ScalarMappable(cmap="jet_r"),ax=ax,orientation="vertical",
            shrink=0.8,location="right",pad=.025)
        cax.set_ticks([0.05,0.5,0.95])
        cax.set_ticklabels(["tolerant","neutral","averse"],rotation=90,fontsize=12,va="center")
        #cax.set_label("increasing reliability",fontsize=12)
        ax.set_title("generation {0:03d}".format(g),loc="left",fontsize=12)
        plt.tight_layout()
        plt.savefig(os.path.join(plt_d,"mou_{0:03d}.png".format(i)),dpi=400)
        plt.close(fig)

    fps = 15
    pyemu.os_utils.run("ffmpeg -i mou_{0:03d}.png -vf palettegen=256 palette.png".format(i),cwd=plt_d)
    pyemu.os_utils.run("ffmpeg -r {0} -y -s 1920X1080 -i mou_%03d.png -i palette.png -filter_complex \"scale=720:-1:flags=lanczos[x];[x][1:v]paletteuse\" logo.gif".format(fps),
            cwd=plt_d)



def stack_invest():
    name = "constroc"
    t_d = mou_suite_helper.setup_problem(name,True,True)
    
    pst = pyemu.Pst(os.path.join(t_d,"{0}.pst".format(name)))
    par = pst.parameter_data
    pst.control_data.noptmax = 50
    pst.pestpp_options["opt_recalc_chance_every"] = pst.control_data.noptmax
    pst.pestpp_options["mou_save_population_every"] = 1
    pst.pestpp_options["opt_stack_size"] = 30
    #pst.pestpp_options["opt_par_stack"] = "prior.csv"
    pst.pestpp_options["mou_generator"] = "pso"
    pst.pestpp_options["mou_population_size"] = 30
    pst.pestpp_options["mou_verbose_level"] = 4
    pst.pestpp_options["opt_chance_points"] = "all"
    pst.pestpp_options["mou_max_archive_size"] = 100000
    pst.pestpp_options["opt_risk"] = 0.95
    
    pst.write(os.path.join(t_d,"{0}.pst".format(name)))
    m1 = os.path.join("mou_tests","master_stack_test")
    pyemu.os_utils.start_workers(t_d,exe_path,"{0}.pst".format(name),20,worker_root="mou_tests",
                                 master_dir=m1,verbose=True,port=port)

    pdf = pd.read_csv(os.path.join(m1,"{0}.{1}.nested.par_stack.csv".\
        format(name,pst.control_data.noptmax)),index_col=0)
    pdf.loc[:,"member"] = pdf.index.map(lambda x: x.split("||")[1])
    umem = pdf.member.unique()
    umem.sort()
    print(umem)
    pdf0 = pdf.loc[pdf.member==umem[-1],:].copy()
    pdf0.index = pdf0.pop("member")
    print(pdf0)
    pdf0.to_csv(os.path.join(t_d,"sweep_in.csv"))
    
    m2 = os.path.join("mou_tests","master_sweep_test")
    # pyemu.os_utils.start_workers(t_d,exe_path.replace("-mou","-swp"),"{0}.pst".format(name),20,worker_root="mou_tests",
    #                              master_dir=m2,verbose=True,port=port)


    odf = pd.read_csv(os.path.join(m1,"{0}.{1}.nested.obs_stack.csv".\
        format(name,pst.control_data.noptmax)),index_col=0)
    odf.loc[:,"member"] = odf.index.map(lambda x: x.split("||")[1])
    odf0 = odf.loc[odf.member == umem[-1],:].copy()
    odf0.index = odf0.pop("member")

    #sdf = pd.read_csv(os.path.join(m2,"sweep_out.csv"),index_col=1)
    #diff = sdf.loc[:,odf0.columns] - odf0
    #print(diff.apply(np.abs).sum())

    arc = pd.read_csv(os.path.join(m1,"{0}.pareto.summary.csv".format(name)))
    arc.loc[:,"memgen"] = arc.member.apply(lambda x: int(x.split('_')[0].split('=')[1]))
    arc = arc.loc[arc.memgen==pst.control_data.noptmax,:]
    umem = arc.member.unique()
    umem.sort()

    chance = pd.read_csv(os.path.join(m1,"constroc.{0}.chance.obs_pop.csv".format(pst.control_data.noptmax)),index_col=0)
    actual = pd.read_csv(os.path.join(m1,"constroc.{0}.obs_pop.csv".format(pst.control_data.noptmax)),index_col=0)

    import matplotlib.pyplot as plt

    from matplotlib.backends.backend_pdf import PdfPages
    obs = pst.observation_data
    with PdfPages("chance.pdf") as pdfpages:
        for mem in umem:

            fig,axes = plt.subplots(obs.shape[0],1,figsize=(5*obs.shape[0],10))
            for oname,ax in zip(obs.obsnme,axes):
                ax.hist(odf.loc[odf.member==mem,:].loc[:,oname],alpha=0.25,facecolor="0.5",label="stack")
                if "const" in oname:
                    ax.plot([obs.loc[oname,"obsval"],obs.loc[oname,"obsval"]],ax.get_ylim(),"r--",lw=2.5,label="rhs")
                ax.plot([chance.loc[mem, oname], chance.loc[mem, oname]], ax.get_ylim(), "b--", lw=2.5,label="chance value")
                ax.plot([actual.loc[mem, oname], actual.loc[mem, oname]], ax.get_ylim(), "k--", lw=2.5,label="sim value")

                risk = pdf.loc[pdf.member==mem,"_risk_"].values[0]
                ax.set_title("output:{0}, group:{1}, member:{2}, risk:{3:3.2f}".format(oname,obs.loc[oname,"obgnme"],mem,float(risk)),loc="left",fontsize=12)
                ax.set_yticks([])
                print(mem, oname)
                ax.set_xlim(min(obs.loc[oname,"obsval"]*.75,odf.loc[:,oname].min()),max(obs.loc[oname,"obsval"]*1.25,odf.loc[:,oname].max()))
                ax.legend(loc="upper left",fontsize=12)
            plt.tight_layout()
            pdfpages.savefig()
            plt.close(fig)



def run():
    pyemu.os_utils.start_workers(os.path.join("mou_tests","constroc_template"),exe_path,"constroc.pst",worker_root="mou_tests",num_workers=20,master_dir=os.path.join("mou_tests","constroc_mimic"))


def zdt1_fixed_scaleoffset_test():
    t_d = mou_suite_helper.setup_problem("zdt1",False,False)
    pst = pyemu.Pst(os.path.join(t_d,"zdt1.pst"))
    pe = pyemu.ParameterEnsemble.from_uniform_draw(pst,num_reals=10)
    pe.to_csv(os.path.join(t_d,"init_pop.csv"))
    par = pst.parameter_data
    first = pst.par_names[1]
    others = pst.par_names[2:5]
    par.loc[others,"partrans"] = "fixed"
    par.loc[others,"scale"] = -1
    
    pst.pestpp_options["opt_recalc_chance_every"] = 1000
    pst.pestpp_options["opt_dec_var_groups"] = "decvars"
    pst.pestpp_options["mou_save_population_every"] = 1
    #pst.pestpp_options["opt_stack_size"] = 10
    #pst.pestpp_options["opt_par_stack"] = "prior.csv"
    pst.pestpp_options["mou_generator"] = "de"
    pst.pestpp_options["mou_population_size"] = 10
    pst.control_data.noptmax = 0
    pst.write(os.path.join(t_d,"zdt1.pst"))
    pyemu.os_utils.run("{0} {1}".format(exe_path,"zdt1.pst"),cwd=t_d)
    df = pd.read_csv(os.path.join(t_d,"dv.dat"),header=None,names=["dv","val"],sep='\s+')
    df.index = df.dv 

    print(df.loc[others,:])
    assert np.all(df.loc[others,"val"].values < 0)
    
    
    pst.pestpp_options["mou_dv_population_file"] = "init_pop.csv"
    pst.control_data.noptmax = -1
    pst.write(os.path.join(t_d,"zdt1.pst"))
    
    pyemu.os_utils.run("{0} {1}".format(exe_path,"zdt1.pst"),cwd=t_d)
    df = pd.read_csv(os.path.join(t_d,"dv.dat"),header=None,names=["dv","val"],sep='\s+')
    df.index = df.dv 

    print(df.loc[others,:])
    assert np.all(df.loc[others,"val"].values <= 0)
    

def zdt1_fixedtied_stack_test():
    t_d = mou_suite_helper.setup_problem("zdt1",True,True)
    #df = pd.read_csv(os.path.join(t_d,"prior.csv"),index_col=0)
    #df.loc[:,"_risk_"] = 0.95
    #print(df.columns)
    #df.to_csv(os.path.join(t_d,"prior.csv"))
    pst = pyemu.Pst(os.path.join(t_d,"zdt1.pst"))
    
    df = pd.read_csv(os.path.join(t_d,"dv.dat"),header=None,names=["dv","val"],sep='\s+') 
    ins_file = os.path.join(t_d,"dv.dat.ins")
    with open(ins_file,'w') as f:
        f.write("pif ~\n")
        for name in df.dv:
            f.write("l1 w !{0}!\n".format(name))
    pst.add_observations(ins_file,ins_file.replace(".ins",""),pst_path=".")

    df = pd.read_csv(os.path.join(t_d,"additive_par.dat"),header=None,names=["dv","val"],sep='\s+') 
    ins_file = os.path.join(t_d,"additive_par.dat.ins")
    with open(ins_file,'w') as f:
        f.write("pif ~\n")
        for name in df.dv:
            f.write("l1 w !{0}!\n".format(name))
    pst.add_observations(ins_file,ins_file.replace(".ins",""),pst_path=".")

    par = pst.parameter_data
    log_pars = pst.par_names[:2]
    par.loc[log_pars,"partrans"] = "log"
    par.loc[log_pars, "parlbnd"] = 0.0001
    first = pst.par_names[1]
    others = pst.par_names[2:-4]
    par.loc[others,"partrans"] = "tied"
    par.loc[others,"partied"] = first
    #par.loc[others,"pargp"] = "tiedup"
    par.loc[others[0],"partrans"] = "fixed"
    par.loc[others[0],"parval1"] = par.loc[others[0],"parubnd"]
    tied = par.loc[par.partrans=="tied","parnme"].values
    
    pe = pyemu.ParameterEnsemble.from_uniform_draw(pst,num_reals=10)
    pe.to_csv(os.path.join(t_d,"init_pop.csv"))


    #par.loc[first,"parval1"] = 1.0

    #pst.pestpp_options["mou_dv_population_file"] = "prior.csv"
    #pst.pestpp_options["opt_chance_points"] = "single"
    pst.pestpp_options["opt_recalc_chance_every"] = 1
    pst.pestpp_options["opt_dec_var_groups"] = "decvars"
    pst.pestpp_options["mou_save_population_every"] = 1
    pst.pestpp_options["opt_stack_size"] = 10
    #pst.pestpp_options["opt_par_stack"] = "prior.csv"
    pst.pestpp_options["mou_generator"] = "de"
    pst.pestpp_options["mou_population_size"] = 10
    pst.pestpp_options["opt_risk"] = 0.95
    pst.pestpp_options["ensemble_output_precision"] = 10
    #pst.pestpp_options["mou_dv_population_file"] = "init_pop.csv"
    pst.control_data.noptmax = 0
    pst.write(os.path.join(t_d,"zdt1.pst"))
    pyemu.os_utils.run("{0} {1}".format(exe_path,"zdt1.pst"),cwd=t_d)
    df = pd.read_csv(os.path.join(t_d,"dv.dat"),header=None,names=["dv","val"],sep='\s+')
    df.index = df.dv
    assert (df.loc[df.dv!=others[0],"val"] - df.loc[first,"val"]).apply(lambda x: np.abs(x)).sum() == 0.0
    df.loc[others[0],"val"] == par.loc[others[0],"parubnd"]
  
    m1 = os.path.join("mou_tests","zdt1_fixedtied_stack_test")
    pst.control_data.noptmax = 1
    pst.write(os.path.join(t_d,"zdt1.pst"))
    pyemu.os_utils.start_workers(t_d,exe_path,"zdt1.pst",10,worker_root="mou_tests",
                                 master_dir=m1,verbose=True,port=port)

    for i in range(pst.control_data.noptmax+1):


        dp = pd.read_csv(os.path.join(m1,"zdt1.{0}.dv_pop.csv").format(i),index_col=0)
        op = pd.read_csv(os.path.join(m1, "zdt1.{0}.obs_pop.csv").format(i), index_col=0)
        assert dp.shape[0] == op.shape[0]
        print(dp.shape[0])
        for ii in range(dp.shape[0]):
            assert dp.index[ii] == op.index[ii]
            ii = dp.index[ii]

            d = np.abs(dp.loc[ii,tied].values - dp.loc[ii,first]).sum()
            assert d == 0,d

        dp = pd.read_csv(os.path.join(m1, "zdt1.{0}.archive.dv_pop.csv").format(i), index_col=0)
        op = pd.read_csv(os.path.join(m1, "zdt1.{0}.archive.obs_pop.csv").format(i), index_col=0)
        assert dp.shape[0] == op.shape[0]
        print(dp.shape[0])
        for ii in range(dp.shape[0]):
            assert dp.index[ii] == op.index[ii]
            ii = dp.index[ii]

            d = np.abs(dp.loc[ii, tied].values - dp.loc[ii, first]).sum()
            assert d == 0, d

        dps = pd.read_csv(os.path.join(m1,"zdt1.{0}.par_stack.csv").format(i),index_col=0)
        ops = pd.read_csv(os.path.join(m1, "zdt1.{0}.obs_stack.csv").format(i), index_col=0)
        for pname in pst.par_names:
            if pname in ops.columns:
                d = (dps.loc[:,pname] - ops.loc[:,pname]).apply(lambda x: np.abs(x)).sum()
                #print(pname,d)
                assert d < 1.0e-7
                if par.loc[pname,"partrans"] == "fixed":
                    assert np.abs(ops.loc[:,pname].values - par.loc[pname,"parval1"]).sum() == 0
                else:    
                     assert np.abs(ops.loc[:,pname].values - par.loc[pname,"parval1"]).min() > 0
                
            for ii in op.index:
                d = np.abs(op.loc[ii,tied].values - op.loc[ii,first]).sum()
            #print(d)
                assert d < 1.0e-6

            #make sure the dv value in the stack is in the dv population - checking for transform issues
            # compare with the archive
            if pname.startswith('dv'):
                stack_ovals = list(np.round(ops.loc[:,pname].values,5))
                pop_ovals = list(np.round(op.loc[:,pname].values,5))
                #print(pop_ovals)
                #print(stack_ovals)
                in_pop = stack_ovals[0] in pop_ovals
                assert in_pop

    pst.pestpp_options["opt_chance_points"] = "all"
    m1 = os.path.join("mou_tests","zdt1_fixedtied_stack_every_test")
    pst.control_data.noptmax = 1
    pst.write(os.path.join(t_d,"zdt1.pst"))
    pyemu.os_utils.start_workers(t_d,exe_path,"zdt1.pst",10,worker_root="mou_tests",
                                 master_dir=m1,verbose=True,port=port)
    for i in range(pst.control_data.noptmax+1):

        dp = pd.read_csv(os.path.join(m1, "zdt1.{0}.archive.dv_pop.csv").format(i), index_col=0)
        op = pd.read_csv(os.path.join(m1, "zdt1.{0}.archive.obs_pop.csv").format(i), index_col=0)
        assert dp.shape[0] == op.shape[0]
        print(dp.shape[0])
        for ii in range(dp.shape[0]):
            assert dp.index[ii] == op.index[ii]
            ii = dp.index[ii]

            d = np.abs(dp.loc[ii, tied].values - dp.loc[ii, first]).sum()
            assert d == 0, d

        dp = pd.read_csv(os.path.join(m1, "zdt1.{0}.dv_pop.csv").format(i), index_col=0)
        op = pd.read_csv(os.path.join(m1, "zdt1.{0}.obs_pop.csv").format(i), index_col=0)
        assert dp.shape[0] == op.shape[0]
        print(dp.shape[0])
        for ii in range(dp.shape[0]):
            assert dp.index[ii] == op.index[ii]
            ii = dp.index[ii]

            d = np.abs(dp.loc[ii, tied].values - dp.loc[ii, first]).sum()
            assert d == 0, d


        dps = pd.read_csv(os.path.join(m1,"zdt1.{0}.nested.par_stack.csv").format(i),index_col=0)
        ops = pd.read_csv(os.path.join(m1, "zdt1.{0}.nested.obs_stack.csv").format(i), index_col=0)
        for pname in pst.par_names:
            if pname in ops.columns:
                d = (dps.loc[:,pname] - ops.loc[:,pname]).apply(lambda x: np.abs(x)).sum()
                #print(pname,d)
                assert d < 1.0e-7
                if par.loc[pname,"partrans"] == "fixed":
                    assert np.abs(ops.loc[:,pname].values - par.loc[pname,"parval1"]).sum() == 0
                else:    
                     assert np.abs(ops.loc[:,pname].values - par.loc[pname,"parval1"]).min() > 0
            for ii in ops.index:
                d = np.abs(ops.loc[ii,tied].values - ops.loc[ii,first]).sum()
                #print(d)
                assert d < 1.0e-6
                # make sure the dv value in the stack is in the dv population - checking for transform issues
                # compare with the non-archive
            if pname.startswith('dv'):
                for real in op.index:
                    rvals = np.round(ops.loc[ops.index.map(lambda x: real in x), pname].values, 5)
                    stack_ovals = list(rvals)
                    pop_oval = op.loc[real, pname]
                    pop_oval = list(np.round([pop_oval], 5))[0]
                    print(pop_oval)
                    print(stack_ovals)
                    print(i,pname,real)
                    in_pop = pop_oval in stack_ovals
                    assert in_pop
    
def multigen_test():
    t_d = mou_suite_helper.setup_problem("tkn")
    pst = pyemu.Pst(os.path.join(t_d,"tkn.pst"))
    obs = pst.observation_data
    #print(obs)
    #obs.loc["const_1","obsval"] = -1
    df = pd.DataFrame(data={"gen":[0,1],"pop_size":[100,20]})
    df.to_csv(os.path.join(t_d,"pop_sched.dat"),index=False,header=False,sep=" ")
    pst.pestpp_options["mou_population_schedule"] = "pop_sched.dat"
    pst.pestpp_options["mou_population_size"] = 15
    pst.pestpp_options["mou_generator"] = "pso"
    pst.pestpp_options["mou_env_selector"] = "nsga"
    pst.pestpp_options["mou_use_multigen_population"] = True
    pst.pestpp_options["mou_save_population_every"] = 1
    pst.pestpp_options["mou_verbose_level"] = 4
    pst.control_data.noptmax = 10
    pst.write(os.path.join(t_d,"tkn.pst"))
    m1 = os.path.join("mou_tests", "test_master_multigen_pso")
    pyemu.os_utils.start_workers(t_d, exe_path, "tkn.pst", 15, worker_root="mou_tests",
                                 master_dir=m1, verbose=True,port=port)

    dp = pd.read_csv(os.path.join(m1,"tkn.0.dv_pop.csv"),index_col=0)
    print(dp.shape)
    assert dp.shape[0] == 100
    op = pd.read_csv(os.path.join(m1,"tkn.0.obs_pop.csv"),index_col=0)
    print(op.shape)
    assert op.shape[0] == 100

    dp = pd.read_csv(os.path.join(m1,"tkn.1.dv_pop.csv"),index_col=0)
    print(dp.shape)
    assert dp.shape[0] == 20
    op = pd.read_csv(os.path.join(m1,"tkn.1.obs_pop.csv"),index_col=0)
    print(op.shape)
    assert op.shape[0] == 20

    dp = pd.read_csv(os.path.join(m1,"tkn.{0}.dv_pop.csv".format(pst.control_data.noptmax)),index_col=0)
    print(dp.shape)
    assert dp.shape[0] == 15
    op = pd.read_csv(os.path.join(m1,"tkn.{0}.obs_pop.csv".format(pst.control_data.noptmax)),index_col=0)
    print(op.shape)
    assert op.shape[0] == 15
    
    
def zdt1_fixed_robust_opt_test():
    t_d = mou_suite_helper.setup_problem("zdt1",True,False)
    pst = pyemu.Pst(os.path.join(t_d,"zdt1.pst"))
    
    df = pd.read_csv(os.path.join(t_d,"dv.dat"),header=None,names=["dv","val"],sep='\s+') 
    ins_file = os.path.join(t_d,"dv.dat.ins")
    with open(ins_file,'w') as f:
        f.write("pif ~\n")
        for name in df.dv:
            f.write("l1 w !{0}!\n".format(name))
    pst.add_observations(ins_file,ins_file.replace(".ins",""),pst_path=".")

    df = pd.read_csv(os.path.join(t_d,"additive_par.dat"),header=None,names=["dv","val"],sep='\s+') 
    ins_file = os.path.join(t_d,"additive_par.dat.ins")
    with open(ins_file,'w') as f:
        f.write("pif ~\n")
        for name in df.dv:
            f.write("l1 w !{0}!\n".format(name))
    pst.add_observations(ins_file,ins_file.replace(".ins",""),pst_path=".")

    par = pst.parameter_data
    pe = pyemu.ParameterEnsemble.from_uniform_draw(pst,num_reals=10)
    pe.to_csv(os.path.join(t_d,"init_pop.csv"))
    par.loc[df.dv.values,"partrans"] = "fixed"
    
    pst.pestpp_options["opt_dec_var_groups"] = "decvars"
    pst.pestpp_options["mou_save_population_every"] = 1
    pst.pestpp_options["mou_generator"] = "pso"
    pst.pestpp_options["mou_population_size"] = 10
    pst.pestpp_options["opt_risk"] = 0.5
    pst.pestpp_options["ensemble_output_precision"] = 10
    pst.pestpp_options["mou_dv_population_file"] = "init_pop.csv"
    pst.control_data.noptmax = 0
    pst.write(os.path.join(t_d,"zdt1.pst"))
    pyemu.os_utils.run("{0} {1}".format(exe_path,"zdt1.pst"),cwd=t_d)
    df = pd.read_csv(os.path.join(t_d,"additive_par.dat"),header=None,names=["dv","val"],sep='\s+')
    assert np.all(par.loc[df.dv.values,"parval1"].values == df.loc[:,"val"].values)
  
    m1 = os.path.join("mou_tests","zdt1_fixed_robust_opt_test")
    pst.control_data.noptmax = 3
    pst.write(os.path.join(t_d,"zdt1.pst"))
    pyemu.os_utils.start_workers(t_d,exe_path,"zdt1.pst",10,worker_root="mou_tests",
                                 master_dir=m1,verbose=True,port=port)
    pe0 = pd.read_csv(os.path.join(m1,"zdt1.0.dv_pop.csv"),index_col=0)
    assert pe0.shape == pe0.shape
    d = np.abs(pe._df.values - pe0.values)
    print(d.max())
    assert d.max() < 1.0e-6

    oe0 = pd.read_csv(os.path.join(m1,"zdt1.0.obs_pop.csv"),index_col=0)
    oe0 = oe0.loc[:,pe.columns]
    assert oe0.shape == oe0.dropna().shape
    d = np.abs(pe._df.values - oe0.values)
    print(d.max())
    assert d.max() < 1.0e-6

    pe1 = pd.read_csv(os.path.join(m1,"zdt1.1.dv_pop.csv"),index_col=0)
    assert pe.shape == pe1.shape
    fpar = par.loc[par.partrans=="fixed","parnme"].values
    pe1f = pe1.loc[:,fpar]
    pef = pe0.loc[:,fpar]
    for f in fpar:
        vals = set(list(pef.loc[:,f].values))
        for v in pe1.loc[:,f].values:
            if v not in vals:
                print(v,vals)
                raise Exception("not found")
    oe1 = pd.read_csv(os.path.join(m1,"zdt1.1.obs_pop.csv"),index_col=0)
    oe1 = oe1.loc[:,pe.columns]
    assert oe1.shape == oe1.dropna().shape 
    for f in fpar:
        vals = set(list(pef.loc[:,f].values))
        for v in oe1.loc[:,f].values:
            if v not in vals:
                print(v,vals)
                raise Exception("not found") 

    pe0 = pe0.iloc[:5,:]
    pe0.to_csv(os.path.join(t_d,"init_pop.csv"))
    oe0 = pd.read_csv(os.path.join(m1,"zdt1.0.obs_pop.csv"),index_col=0)
    oe0 = oe0.iloc[:5,:]
    oe0.to_csv(os.path.join(t_d,"init_obs_pop.csv"))
    pst.pestpp_options["mou_obs_population_restart_file"] = "init_obs_pop.csv"
    pst.control_data.noptmax = 3
    pst.write(os.path.join(t_d,"zdt1.pst"))

    pyemu.os_utils.start_workers(t_d,exe_path,"zdt1.pst",10,worker_root="mou_tests",
                                 master_dir=m1,verbose=True,port=port)


def zdt1_chance_schedule_test():
    t_d = mou_suite_helper.setup_problem("zdt1",True,True)
    pst = pyemu.Pst(os.path.join(t_d,"zdt1.pst"))
    


    pst.pestpp_options["opt_recalc_chance_every"] = 1000
    gens_in = [0,5,8,9,10]
    with open(os.path.join(t_d,"chance_schedule.dat"),'w') as f:  
        for gen in gens_in:
            f.write("{0},True\n".format(gen))
    pst.pestpp_options["opt_chance_schedule"] = "chance_schedule.dat"


    pst.pestpp_options["opt_dec_var_groups"] = "decvars"
    pst.pestpp_options["mou_save_population_every"] = 1
    #pst.pestpp_options["opt_stack_size"] = 10
    #pst.pestpp_options["opt_par_stack"] = "prior.csv"
    pst.pestpp_options["mou_generator"] = "de"
    pst.pestpp_options["mou_population_size"] = 10
    pst.pestpp_options["opt_stack_size"] = 5
    pst.pestpp_options["opt_risk"] = 0.95
    pst.pestpp_options["opt_chance_points"] = "all"
    pst.pestpp_options["opt_recalc_chance_every"] = 10000
    
    pst.control_data.noptmax = 11
    pst.write(os.path.join(t_d,"zdt1.pst"))
    pyemu.os_utils.run("{0} {1}".format(exe_path,"zdt1.pst"),cwd=t_d)
    t_d = os.path.join("mou_tests","zdt1_template")
    with open(os.path.join(t_d,"zdt1.rec"),'r') as f:
        for line in f:
            if "chance runs for generation" in line:
                gen = int(line.strip().split()[-1])
                gens_in.remove(gen)
    if len(gens_in) > 0:
        raise Exception(str(gens_in))


    # chance_files = [f for f in os.listdir(t_d) if "stack_summary." in f]
    # print(chance_files)
    # for f in chance_files:
    #     gen = int(f.split(".")[1])
    #     print(gen)
    #     if gen not in gens_in:
    #         raise Exception(str(gen))
    #     gens_in.remove(gen)
    # if len(gens_in) > 0:
    #     raise Exception(str(gens_in))
    

def basic_empcov_invest(case="hosaki"):

    t_d = mou_suite_helper.setup_problem(case, additive_chance=True, risk_obj=False)
    pst = pyemu.Pst(os.path.join(t_d, case+".pst"))
    pst.pestpp_options["mou_generator"] = "empcov"
    pst.pestpp_options["opt_risk"] = 0.95
    pst.control_data.noptmax = 20
    pst.pestpp_options["mou_population_size"] = 10
    pst.pestpp_options["mou_save_population_every"] = 1
    pst.write(os.path.join(t_d, case+".pst"))
    m_d = os.path.join("mou_tests", case+"_empcov_master_risk")
    pyemu.os_utils.start_workers(t_d, exe_path,  case+".pst", 50, worker_root="mou_tests",
                                master_dir=m_d, verbose=True, port=port)
    assert os.path.exists(os.path.join(m_d,"{0}.{1}.fosm.jcb".format(case,pst.control_data.noptmax)))

    for i in range(0,pst.control_data.noptmax+1):
        dv_file = os.path.join(m_d,"{0}.{1}.dv_pop.csv".format(case,i))
        oe_file = os.path.join(m_d,"{0}.{1}.obs_pop.csv".format(case,i))
        assert os.path.exists(dv_file)
        assert os.path.exists(oe_file)
        dv_pso = pd.read_csv(dv_file)
        oe_pso = pd.read_csv(oe_file)
        assert dv_pso.shape[0] == oe_pso.shape[0]

    # pst.pestpp_options["mou_generator"] = "de"
    # pst.pestpp_options["opt_risk"] = 0.95
    # pst.pestpp_options["mou_save_population_every"] = 1
    # pst.write(os.path.join(t_d, case+".pst"))
    # m_d = os.path.join("mou_tests", case+"_de_master_risk")
    #pyemu.os_utils.start_workers(t_d, exe_path, case+".pst", 20, worker_root="mou_tests",
    #                             master_dir=m_d, verbose=True, port=port)
    # assert os.path.exists(os.path.join(m_d, "{0}.{1}.fosm.jcb".format(case, pst.control_data.noptmax)))

    # for i in range(0,pst.control_data.noptmax+1):
    #     dv_file = os.path.join(m_d,"{0}.{1}.dv_pop.csv".format(case,i))
    #     oe_file = os.path.join(m_d,"{0}.{1}.obs_pop.csv".format(case,i))
    #     assert os.path.exists(dv_file)
    #     assert os.path.exists(oe_file)
    #     dv_de = pd.read_csv(dv_file)
    #     oe_de = pd.read_csv(oe_file)
    #     assert dv_de.shape[0] == oe_de.shape[0]

    method = mou_suite_helper.zdt1
    x0 = np.linspace(0,1,10000)
    o1,o2 = [],[]
    for xx0 in x0:
        x = np.zeros(30)
        x[0] = xx0
        ret_vals = method(x)
        o1.append(ret_vals[0][0])
        o2.append(ret_vals[0][1])
    

    # o1 = np.array(o1)
    # o2 = np.array(o2)
    # diff = np.abs(o1.min() - oe_pso.loc[:,"obj_1"].values.min()) 
    # #print(diff)
    # assert diff < 1.0e-4
    # diff = np.abs(o1.max() - oe_pso.loc[:,"obj_1"].values.max()) 
    # #print(diff)
    # assert diff < 1.0e-4

    # diff = np.abs(o2.min() - oe_pso.loc[:,"obj_2"].values.min()) 
    # print(diff)
    # #assert diff < 1.0e-4
    
    # opt_1,opt_2 = o1.min(),o2.min()
    # opt = np.array([opt_1,opt_2])
    # truth = np.array([o1,o2]).transpose()
    # #print(truth)
    # #print(opt)
    # dist = [(opt-t).sum()**2 for t in truth]
    # knee_idx = np.argmin(dist)
    # #print(knee_idx,dist[knee_idx],truth[knee_idx])
    # knee = truth[knee_idx]

    # knee_dist = [(knee-sol).sum()**2 for sol in oe_pso.loc[:,["obj_1","obj_2"]].values]
    # #print(knee_dist)
    # knee_sol_idx = np.argmin(knee_dist)
    
    # knee_sol = oe_pso.loc[:,["obj_1","obj_2"]].values[knee_sol_idx]
    # #print(knee_sol_idx,knee_sol)

    # dist = (knee - knee_sol).sum()**2
    # print(dist)
    #assert dist < 0.001

    # import matplotlib.pyplot as plt
    # fig,axes = plt.subplots(1,2,figsize=(10,5))
    # axes[0].plot(o1,o2)
    # axes[0].scatter(oe_pso.loc[:,"obj_1"].values,oe_pso.loc[:,"obj_2"])
    # axes[1].plot(o1,o2)
    # axes[1].scatter(oe_de.loc[:,"obj_1"].values,oe_de.loc[:,"obj_2"])
    # plt.show()

def plot_hosaki(m_d):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    name = "hosaki"
    pst = pyemu.Pst(os.path.join(m_d,"{0}.pst".format(name)))
    df = pd.read_csv(os.path.join(m_d,"{0}.pareto.archive.summary.csv".format(name)))
    gens = df.generation.unique()
    gens.sort()
    plt_d = m_d+"_plots"
    if os.path.exists(plt_d):
        shutil.rmtree(plt_d)
    os.makedirs(plt_d)
    for i,g in enumerate(gens):
        #gdf = df.loc[df.generation==g,:].copy()
        gdf = pd.read_csv(os.path.join(m_d,"{0}.{1}.archive.obs_pop.csv".format(name,g)))
        dvdf = pd.read_csv(os.path.join(m_d,"{0}.{1}.archive.dv_pop.csv".format(name,g)))
        dvdf.index = dvdf.real_name
        fig,ax = plt.subplots(1,1,figsize=(5.5,5))
        c = 'b'
        if "_risk_" in dvdf.columns:
            c = dvdf.loc[gdf.real_name.values,"_risk_"].values
        #print(g,c)
        ax.scatter(dvdf.dv_0.values,dvdf.dv_1.values,marker=".",s=20,c=c,cmap="jet_r")

        gdf = pd.read_csv(os.path.join(m_d,"{0}.{1}.obs_pop.csv".format(name,g)))
        dvdf = pd.read_csv(os.path.join(m_d,"{0}.{1}.dv_pop.csv".format(name,g)))
        dvdf.index = dvdf.real_name
        ax.scatter(dvdf.dv_0.values,dvdf.dv_1.values,marker='.',s=20,c='0.5',alpha=0.5)
        #ax.set_xlim(2,-0.75)
        #ax.set_ylim(-1,1)
        ax.set_xlim(0,5)
        ax.set_ylim(0,5)
        ax.grid()
        cax = fig.colorbar(mpl.cm.ScalarMappable(cmap="jet_r"),ax=ax,orientation="vertical",
            shrink=0.8,location="right",pad=.025)
        cax.set_ticks([0.05,0.5,0.95])
        cax.set_ticklabels(["tolerant","neutral","averse"],rotation=90,fontsize=12,va="center")
        #cax.set_label("increasing reliability",fontsize=12)
        ax.set_title("generation {0:03d}".format(g),loc="left",fontsize=12)
        plt.tight_layout()
        plt.savefig(os.path.join(plt_d,"mou_{0:03d}.png".format(i)),dpi=400)
        plt.close(fig)

    fps = 15
    pyemu.os_utils.run("ffmpeg -i mou_{0:03d}.png -vf palettegen=256 palette.png".format(i),cwd=plt_d)
    pyemu.os_utils.run("ffmpeg -r {0} -y -s 1920X1080 -i mou_%03d.png -i palette.png -filter_complex \"scale=720:-1:flags=lanczos[x];[x][1:v]paletteuse\" logo.gif".format(fps),
            cwd=plt_d)


def pi_output_test():
    t_d = mou_suite_helper.setup_problem("constr", additive_chance=True, risk_obj=False)
    pst = pyemu.Pst(os.path.join(t_d, "constr.pst"))
    pst.pestpp_options["opt_chance_points"] = "all"
    #pst.pestpp_options["opt_recalc_chance_every"] = 5
    #pst.pestpp_options["opt_stack_size"] = 10
    pst.pestpp_options["mou_generator"] = "de"
    pst.pestpp_options["mou_population_size"] = 10
    pst.pestpp_options["opt_risk"] = 0.5
    pst.control_data.noptmax = 1
    pst.write(os.path.join(t_d, "constr.pst"))
    m1 = os.path.join("mou_tests", "constr_test_master_pieq")
    pyemu.os_utils.start_workers(t_d, exe_path, "constr.pst", 35, worker_root="mou_tests",
                                 master_dir=m1, verbose=True,port=port)
    assert len([f for f in os.listdir(m1) if "pi_pop" in f]) > 0





if __name__ == "__main__":
    basic_pso_test()
    #test_restart_all()
    #chance_consistency_test()
    #zdt1_chance_schedule_test()
    #gpr_run_riskobj_baselines()
    #gpr_compare_invest()
    
    #zdt1_fixed_robust_opt_test()
    #multigen_test()
    #basic_empcov_invest()
    #plot_hosaki(m_d=os.path.join("mou_tests","hosaki_empcov_master_risk"))
    #basic_pso_test()
    #zdt1_fixedtied_stack_test()
    #zdt1_fixed_scaleoffset_test()
    #shutil.copy2(os.path.join("..","exe","windows","x64","Debug","pestpp-mou.exe"),os.path.join("..","bin","pestpp-mou.exe"))
    #invest()
    #plot()
    #run()
    #stack_invest()
    #zdt1_stack_run_for_animation(mou_gen="pso")
    #plot_zdt1(name="zdt1",m_d=os.path.join("mou_tests","zdt1_empcov_master_risk"))
    #plot_zdt1(name="zdt1",m_d=os.path.join("mou_tests","master_zdt1_test_de_ro"))

    #zdt1_stack_run_for_animation(mou_gen="pso")
    #plot_zdt1(name="zdt1",m_d=os.path.join("mou_tests","master_zdt1_test_pso"))
    #plot_zdt1(name="zdt1",m_d=os.path.join("mou_tests","master_zdt1_test_pso_ro"))

    #zdt1_tied_test()
    #basic_pso_test()

    #shutil.copy2(os.path.join("..", "bin", "win", "pestpp-mou.exe"),
    #             os.path.join("..", "bin", "pestpp-mou.exe"))
    #basic_pso_test()
    #risk_obj_test()
    #invest_2()
    #chance_consistency_test()
    #invest_3()
    # mou_suite_helper.start_workers("zdt1")
    #all_infeas_test()
    #invest_4()
    #restart_dv_test()
    #chance_all_binary_test()
    #invest_5()
    #constr_risk_demo()
    #plot_constr_risk_demo()

    #risk_demo(case='kur',noptmax=300,std_weight=0.01)
    #plot_risk_demo_multi(case='kur')

    #risk_demo(case='zdt1',noptmax=300,std_weight=0.00001)
    #plot_risk_demo_multi(case='zdt1')

    #risk_demo(case="rosenc",std_weight=1.0,noptmax=500)
    #plot_risk_demo_multi()
    #plot_risk_demo_rosen()
    #risk_demo(case='constr', noptmax=150, std_weight=0.05, pop_size=100, num_workers=50, mou_gen="de")
    #risk_demo(case='constr', noptmax=150, std_weight=0.05, pop_size=100, mou_gen="de", num_workers=50)

    #risk_demo(case='constr',noptmax=20,std_weight=0.05,pop_size=100,num_workers=50,mou_gen="simplex,de")
    #risk_demo(case='constr', noptmax=20, std_weight=0.05, pop_size=100,mou_gen="de",num_workers=50)
    #plot_risk_demo_multi_3pane(case='zdt1',mou_gen="de")
    #plot_zdt_risk_demo_compare(case="constr")
    #zdt1_invest()
    #test_sorting_fake_problem()
    #plot_risk_demo_rosen()
    #risk_demo(case="rosenc",std_weight=1.0,mou_gen="pso",pop_size=100,noptmax=30)
    #plot_risk_demo_rosen(mou_gen="pso")
    #risk_demo(case="rosenc", std_weight=1.0, mou_gen="de",pop_size=100,noptmax=30)
    #plot_risk_demo_rosen(mou_gen="de")
    #all_infeas_test()
    #case = "constr"
    #basic_pso_test(case=case)
    #water_invest()
    #mou_suite_helper.plot_results(os.path.join("mou_tests",case+"_pso_master_risk"),sequence=True)
    #mou_suite_helper.plot_results(os.path.join("mou_tests",case+"_de_master_risk"),sequence=True)
    #plot_constr_risk()
    #plot_constr_risk_pub()
    #stack_map_invest()

    #pop_sched_test()
    #simplex_invest_1()
    pi_output_test()
