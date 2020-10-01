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
exe_path = os.path.join(bin_path, "pestpp-mou" + exe)


noptmax = 4
num_reals = 20
port = 4021
test_root = "mou_tests"



# def mf6_v5_glm_test():
#     model_d = "mf6_freyberg"
#     local=True
#     if "linux" in platform.platform().lower() and "10par" in model_d:
#         #print("travis_prep")
#         #prep_for_travis(model_d)
#         local=False
    
#     t_d = os.path.join(model_d,"template")
#     m_d = os.path.join(model_d,"master_glm")
#     if os.path.exists(m_d):
#         shutil.rmtree(m_d)
#     pst = pyemu.Pst(os.path.join(t_d,"freyberg6_run_glm.pst"))
#     m_d = os.path.join(model_d,"master_glm")
#     pyemu.os_utils.start_workers(t_d, "pestpp-glm", "freyberg6_run_glm.pst", 
#                                  num_workers=15, master_dir=m_d,worker_root=model_d,
#                                  port=port)

#     oe_file = os.path.join(m_d,"freyberg6_run_glm.post.obsen.csv")
#     assert os.path.exists(oe_file)
#     oe = pd.read_csv(oe_file)
#     assert oe.shape[0] == pst.pestpp_options["glm_num_reals"],"{0},{1}".\
#         format(oe.shape[0],pst.pestpp_options["glm_num_reals"])


def zdt1(x):
    g = 1 + 9 * np.sum(x[1:]) / (len(x) - 1)
    return x[0], g * (1 - np.sqrt(x[0] / g))

def zdt2(x):
    g = 1 + 9 * np.sum(x[1:]) / (len(x) - 1)
    return x[0], g * (1 - np.power(x[0] / g, 2))

def zdt3(x):
    g = 1 + 9 * np.sum(x[1:]) / (len(x) - 1)
    return x[0], g * (1 - np.sqrt(x[0] / g) - (x[0] / g) *\
           np.sin(10 * np.pi * x[0]))

def zdt4(x):
    g = 1 + 10 * (len(x) - 1) + np.sum(np.power(x[1:], 2) - 10 * np.cos(4 * np.pi * x[1:]))
    return x[0], g * (1 - np.sqrt(x[0] / g))

def zdt6(x):
    f1 = 1 - np.exp(-4 * x[0]) * np.power(np.sin(6 * np.pi * x[0]), 6)
    g = 1 + 9 * np.power(np.sum(x[1:]) / (len(x) - 1), 0.25)
    return f1, g * (1 - np.power(f1 / g, 2))

def constr(x):
    return x[0],(1 + x[1]) / x[0]

def srn(x):
    const1 = np.power(x[0], 2) + np.power(x[1], 2)  # lest than or equal to 225
    const2 =  3 * x[1] - x[0]  # greater than or equal to 10
    f1 = np.power(x[0] - 2, 2) + np.power(x[1] - 1, 2) + 2
    f2 = 9 * x[0] - np.power(x[1] - 1, 2)
    return const1, const2, f1, f2


def zdt_helper(func):
    pdf = pd.read_csv("dv.dat",delim_whitespace=True,index_col=0, header=None, names=["parnme","parval1"])
    obj1,obj2 = func(pdf.values)
    if os.path.exists("additive_par.dat"):
        cdf = pd.read_csv("additive_par.dat", delim_whitespace=True, index_col=0,header=None, names=["parnme","parval1"])
        obj1[0] += cdf.parval1.values[0]
        obj2[0] += cdf.parval1.values[1]

    with open("obj.dat",'w') as f:
        f.write("obj_1 {0}\n".format(float(obj1)))
        f.write("obj_2 {0}\n".format(float(obj2)))

def setup_zdt_problem(name,num_dv,additive_chance=False):
    test_d = os.path.join(test_root,"{0}_template".format(name))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    os.makedirs(test_d)
    
    # write a generic template file for the dec vars
    tpl_file = "dv.dat.tpl".format(name)
    with open(os.path.join(test_d,tpl_file),'w') as f:
        f.write("ptf ~\n")
        for i in range(num_dv):
            f.write("dv_{0} ~ dv_{0}     ~\n".format(i))
    
    additive_chance_tpl_file = None
    if additive_chance:
        additive_chance_tpl_file = "additive_par.dat.tpl"
        with open(os.path.join(test_d,additive_chance_tpl_file),'w') as f:
            f.write("ptf ~\n")
            f.write("obj1_add_par ~   obj1_add_par   ~\n")
            f.write("obj2_add_par ~   obj2_add_par   ~\n")

        with open(os.path.join(test_d,additive_chance_tpl_file.replace(".tpl","")),'w') as f:
            f.write("obj1_add_par 0.0\n")
            f.write("obj2_add_par 0.0\n")
            
    with open(os.path.join(test_d,tpl_file.replace(".tpl","")),'w') as f:
        for i in range(num_dv):
            f.write("dv_{0} 0.5\n".format(i))

    # write a generic ins file for the two objectives
    ins_file = "obj.dat.ins".format(name)
    with open(os.path.join(test_d,ins_file),'w') as f:
        f.write("pif ~\n")
        f.write("l1 w !obj_1!\n")
        f.write("l1 w !obj_2!\n")

    # now scape this python file to get the function lines and
    # the helper lines
    lines = open("mou_tests.py",'r').readlines()
    func_lines = []
    for i in range(len(lines)):
        if lines[i].startswith("def {0}(x):".format(name)):
            func_lines.append(lines[i])
            for ii in range(i+1,len(lines)):
                if lines[ii][0] not in [" ","   ","\n"]:
                    break
                func_lines.append(lines[ii])
            break
    #print(func_lines)
    if len(func_lines) == 0:
        raise Exception()
    helper_lines = []
    for i in range(len(lines)):
        if lines[i].startswith("def zdt_helper(".format(name)):
            helper_lines.append(lines[i])
            for ii in range(i+1,len(lines)):
                if lines[ii][0] not in [" ","   ","\n"]:
                    break
                helper_lines.append(lines[ii])
            break
    #print(helper_lines)
    if len(helper_lines) == 0:
        raise Exception()
    #helper_lines[0] = "def zdt_helper({0}):\n".format(name)

    # write these functions to the forward run script
    with open(os.path.join(test_d,"forward_run.py"),'w') as f:
        f.write("import os\nimport numpy as np\nimport pandas as pd\n")
        for func_line in func_lines:
            f.write(func_line)
        for helper_line in helper_lines:
            f.write(helper_line)

        f.write("if __name__ == '__main__':\n    zdt_helper({0})\n".format(name))

    # make sure it runs
    pyemu.os_utils.run("python forward_run.py",cwd=test_d)
  
    # create the control file
    tpl_file = os.path.join(test_d,tpl_file)
    ins_file = os.path.join(test_d,ins_file)
    pst = pyemu.Pst.from_io_files(tpl_files=tpl_file,in_files=tpl_file.replace(".tpl",""),
                                   ins_files=ins_file,out_files=ins_file.replace(".ins",""),pst_path=".")

    par = pst.parameter_data
    par.loc[:,"parubnd"] = 1.0
    par.loc[:,"parlbnd"] = 0.0
    par.loc[:,"partrans"] = "none"
    par.loc[:,"parval1"] = 0.5
    par.loc[:,"pargp"] = "decvars"
    if additive_chance_tpl_file is not None:
        adf = pst.add_parameters(os.path.join(test_d,additive_chance_tpl_file),pst_path=".")
        print(adf)
        
        par = pst.parameter_data
        par.loc[adf.parnme,"partrans"] = "none"
        par.loc[adf.parnme,"parubnd"] = 1.0
        par.loc[adf.parnme,"parval1"] = 0.0
        par.loc[adf.parnme,"parlbnd"] = -1.0
        par.loc[adf.parnme,"parchglim"] = "relative"
        par.loc[adf.parnme,"pargp"] = "obj_add"
        pst.rectify_pgroups()
        pst.parameter_groups.loc["obj_add","inctyp"] = "absolute"


    obs = pst.observation_data
    obs.loc[:,"weight"] = 1.0
    obs.loc[:,"obgnme"] = "less_than_obj" # all these zdt probs are min-min

    pst.pestpp_options["opt_dec_var_groups"] = "decvars"
    pst.model_command = "python forward_run.py"
    pst.control_data.noptmax = 0
    pst.write(os.path.join(test_d,name+".pst"))
    
    # run mou with noptmax = 0 to make sure we are getting something
    pyemu.os_utils.run("{0} {1}.pst".format(exe_path,name),cwd=test_d)
    pst = pyemu.Pst(os.path.join(test_d,name+".pst"))
    print(pst.phi)
    assert pst.phi < 1.0e-10
    return test_d

def test_zdt1():
    test_case = "zdt1"
    test_d = setup_zdt_problem(test_case,30,additive_chance=False)
    pst = pyemu.Pst(os.path.join(test_d,"{0}.pst".format(test_case)))
    pst.control_data.noptmax = 1
    pst.pestpp_options["mou_population_size"] = 100
    pst.write(os.path.join(test_d,"{0}.pst".format(test_case)))
    #pyemu.os_utils.run("{0} {1}.pst".format(exe_path,test_case),cwd=test_d)
    master_d = test_d.replace("template","master")
    pyemu.os_utils.start_workers(test_d, exe_path, "{0}.pst".format(test_case), 
                                  num_workers=15, master_dir=master_d,worker_root=test_root,
                                  port=port)

    #TODO: need some asserts here
    dv_pop_file = "{0}.0.dv_pop.csv".format(test_case)
    assert os.path.exists(os.path.join(master_d,dv_pop_file)),dv_pop_file
    obs_pop_file = "{0}.0.obs_pop.csv".format(test_case)
    assert os.path.exists(os.path.join(master_d,obs_pop_file)),obs_pop_file
    dv_df = pd.read_csv(os.path.join(master_d,dv_pop_file),index_col=0)
    obs_df = pd.read_csv(os.path.join(master_d,obs_pop_file),index_col=0)
    assert dv_df.shape[0] == pst.pestpp_options["mou_population_size"]
    assert dv_df.shape[1] == pst.npar
    assert obs_df.shape[0] == pst.pestpp_options["mou_population_size"]
    assert obs_df.shape[1] == pst.nobs
    assert dv_df.index.to_list() == obs_df.index.to_list()

    dv_df = dv_df.iloc[1:,:]
    obs_df = obs_df.iloc[:-1,:]
    dv_df.to_csv(os.path.join(test_d,"restart_dv.csv"))
    obs_df.to_csv(os.path.join(test_d,"restart_obs.csv"))

    shutil.copy2(os.path.join(master_d,dv_pop_file),os.path.join(test_d,"restart_dv.csv"))
    shutil.copy2(os.path.join(master_d,obs_pop_file),os.path.join(test_d,"restart_obs.csv"))
    pst.pestpp_options["mou_dv_population_file"] = "restart_dv.csv"
    pst.pestpp_options["mou_obs_population_restart_file"] = "restart_obs.csv"
    pst.control_data.noptmax = 3
    pst.write(os.path.join(test_d,"{0}.pst".format(test_case)))
    pyemu.os_utils.start_workers(test_d, exe_path, "{0}.pst".format(test_case), 
                                      num_workers=15, master_dir=master_d,worker_root=test_root,
                                      port=port)
    dv_pop_file = "{0}.0.dv_pop.csv".format(test_case)
    assert os.path.exists(os.path.join(master_d,dv_pop_file)),dv_pop_file
    obs_pop_file = "{0}.0.obs_pop.csv".format(test_case)
    assert os.path.exists(os.path.join(master_d,obs_pop_file)),obs_pop_file
    dv_df = pd.read_csv(os.path.join(master_d,dv_pop_file),index_col=0)
    obs_df = pd.read_csv(os.path.join(master_d,obs_pop_file),index_col=0)
    assert dv_df.shape[0] == pst.pestpp_options["mou_population_size"]
    assert dv_df.shape[1] == pst.npar
    assert obs_df.shape[0] == pst.pestpp_options["mou_population_size"]
    assert obs_df.shape[1] == pst.nobs
    assert dv_df.index.to_list() == obs_df.index.to_list()


def test_zdt1_chance():
    test_case = "zdt1"
    test_d = setup_zdt_problem(test_case,30,additive_chance=True)
    pst = pyemu.Pst(os.path.join(test_d,"{0}.pst".format(test_case)))
    pst.control_data.noptmax = 1
    pst.pestpp_options["mou_population_size"] = 5    
    pst.pestpp_options["opt_risk"] = 0.95
    pst.pestpp_options["opt_stack_size"] = 3
    pst.pestpp_options["opt_chance_points"] = "all"
    pst.write(os.path.join(test_d,"{0}.pst".format(test_case)))
    #pyemu.os_utils.run("{0} {1}.pst".format(exe_path,test_case),cwd=test_d)
    master_d = test_d.replace("template","master_chance")
    pyemu.os_utils.start_workers(test_d, exe_path, "{0}.pst".format(test_case), 
                                  num_workers=15, master_dir=master_d,worker_root=test_root,
                                  port=port)


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


if __name__ == "__main__":
        
    #zdt1_test()
    # setup_zdt_problem("zdt1",30)
    # setup_zdt_problem("zdt2",30)
    # setup_zdt_problem("zdt3",30)
    # setup_zdt_problem("zdt4",10)
    # setup_zdt_problem("zdt6",10)
    shutil.copy2(os.path.join("..","exe","windows","x64","Debug","pestpp-mou.exe"),os.path.join("..","bin","pestpp-mou.exe"))
    #setup_zdt_problem("zdt1",30, additive_chance=True)
    #test_zdt1()
    #test_zdt1_chance()
    #setup_zdt_problem("zdt1",30, additive_chance=True)
    test_sorting_fake_problem()