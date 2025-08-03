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
port = 4022
test_root = "mou_tests"


def fon(x):
    f1 = 1.0 - np.exp(-1.0 * np.sum(np.power(x - (1.0 / np.sqrt(3.0)), 2)))
    f2 = 1.0 - np.exp(-1.0 * np.sum(np.power(x + (1.0 / np.sqrt(3.0)), 2)))
    return (np.array([f1]),np.array([f2])),[]

def kur(x):
    xsq = np.power(x,2)
    xsq_sum = xsq[:-1] + xsq[1:]
    f1 = np.sum(-10. * np.exp(-0.2 * np.sqrt(xsq_sum)))
    f2 = np.sum(np.power(np.abs(x),0.8) + (5.0 * np.sin(np.power(x,3))))
    return (np.array([f1]), np.array([f2])), []

def tkn(x):
    g1 = -np.power(x[0],2) - np.power(x[1],2) + 1.0 + 0.1 * np.cos(16. * np.arctan(x[0]/x[1]))
    g2 = np.power((x[0] - 0.5),2) + np.power((x[1]-0.5),2)
    return (x[0],x[1]),[g1,g2]

def water(x):
    f1 = (107680.37 * (x[1] + x[2])) + 61704.67
    f2 = 3000.0 * x[0]
    f3 = (305700. * (2289 * x[1])) / (0.06 * np.power(2289,0.65))
    f4 = 250. * 2289 * np.exp((-39.75 * x[1]) + (9.9 * x[2]) + 2.74) 
    f5 = 25.0 * ((1.39/(x[0] * x[1])) + (4940 * x[2]) - 80)
    g1 = (0.00139 / (x[0] * x[1])) + (4.94 * x[2]) - 0.08
    g2 = 0.000306 / (x[0] * x[1]) + (1.082 * x[2]) - 0.0986
    g3 = 12.307 / (x[0] * x[1]) + (49408.24 * x[2]) + 4051.02
    g4 = 2.098 / (x[0] * x[1]) + (8046.33 * x[2]) + 696.71
    g5 = 2.138 / (x[0] * x[1]) + (7883.39 * x[2]) - 705.04
    g6 = 0.417 / (x[0] * x[1]) + (1721.26 * x[2]) - 136.54
    g7 = 0.164 / (x[0] * x[1]) + (631.13 * x[2]) - 54.48

    return (f1,f2,f3,f4,f5),[g1,g2,g3,g4,g5,g6,g7]
def sch(x):
    return (np.power(x,2),np.power(x-2,2)),[]


def zdt1(x):
    g = 1 + 9 * np.sum(x[1:]) / (len(x) - 1)
    return (x[0]*10, g * (1 - np.sqrt(x[0] / g))*10),[]

def zdt2(x):
    g = 1 + 9 * np.sum(x[1:]) / (len(x) - 1)
    return (x[0], g * (1 - np.power(x[0] / g, 2))),[]

def zdt3(x):
    g = 1 + 9 * np.sum(x[1:]) / (len(x) - 1)
    return (x[0], g * (1 - np.sqrt(x[0] / g) - (x[0] / g) *\
           np.sin(10 * np.pi * x[0]))),[]

def zdt4(x):
    g = 1 + 10 * (len(x) - 1) + np.sum(np.power(x[1:], 2) - 10 * np.cos(4 * np.pi * x[1:]))
    return (x[0], g * (1 - np.sqrt(x[0] / g))),[]

def zdt6(x):
    f1 = 1 - np.exp(-4 * x[0]) * np.power(np.sin(6 * np.pi * x[0]), 6)
    g = 1 + 9 * np.power(np.sum(x[1:]) / (len(x) - 1), 0.25)
    return (f1, g * (1 - np.power(f1 / g, 2))),[]

def constr(x):
    return (x[0],(1 + x[1]) / x[0]),[]

def constroc(x):
    return (x[0],(1 + x[1]) / x[0]),[(9. * x[0])+x[1],(9. * x[0]) - x[1]]

def srn(x):
    const1 = np.power(x[0], 2) + np.power(x[1], 2)  # lest than or equal to 225
    const2 =  x[0] - (3 * x[1])  # less than or equal to -10
    f1 = np.power(x[0] - 2, 2) + np.power(x[1] - 1, 2) + 2
    f2 = 9 * x[0] - np.power(x[1] - 1, 2)
    return (f1, f2), [const1, const2]

def rosen(x):
    f1 = np.power(1 - x[0],2) + (100 * np.power(x[1] - np.power(x[0],2),2))
    return (f1,np.array([-999])),[]

def rosenc(x):
    f1 = np.power(1 - x[0],2) + (100 * np.power(x[1] - np.power(x[0],2),2))
    const1 = (x[0] - 1)**3 - x[1] + 1
    const2 = x[0] + x[1] - 2
    return (f1,np.array([-999])),[const1, const2]

def ackley(x):
    t1 = -20. * np.exp(-0.2 * np.sqrt(0.5 * np.power(x[0],2) + np.power(x[1],2)))
    t2 = -1. * np.exp(0.5 * (np.cos(2. * np.pi * x[0]) + np.cos(2.0 * np.pi * x[1]))) + np.e + 20.0
    return (t1 + t2,np.array([-999])),[]

def binh(x):
    f1 = (4. * x[0]**2) + (4. * x[1]**2)
    f2 = (x[0] - 5)**2 + (x[1]-5)**2
    return (f1,f2),[]

def hosaki(x):    
    pval1 = x[0]
    pval2 = x[1]
    term1 = (pval2**2)*np.e**(-pval2)
    term2 = 1. - (8.*(pval1)) + (7.*(pval1**2))
    term3 = (-(7./3.)*(pval1**3)) + (0.25*(pval1**4))
    sim = (term2 + term3) * term1
    return (sim,np.array([-999])),[]

def helper(func):
    pdf = pd.read_csv("dv.dat",sep='\s+',index_col=0, header=None, names=["parnme","parval1"])
    #obj1,obj2 = func(pdf.values)
    objs,constrs = func(pdf.values)
    
    if os.path.exists("additive_par.dat"):
        obj1,obj2 = objs[0],objs[1]
        cdf = pd.read_csv("additive_par.dat", sep='\s+', index_col=0,header=None, names=["parnme","parval1"])
        obj1[0] += cdf.parval1.values[0]
        obj2[0] += cdf.parval1.values[1]
        for i in range(len(constrs)):
            constrs[i] += cdf.parval1.values[i+2]

    with open("obj.dat",'w') as f:
        for i,obj in enumerate(objs):
            f.write("obj_{0} {1}\n".format(i+1,float(obj)))
        #f.write("obj_2 {0}\n".format(float(obj2)))
        for i,constr in enumerate(constrs):
            f.write("constr_{0} {1}\n".format(i+1,float(constr)))
    return objs,constrs

def setup_problem(name,additive_chance=False, risk_obj=False, self_adaptive=False):
    test_d = os.path.join(test_root,"{0}_template".format(name))
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    os.makedirs(test_d)
    
    num_dv = 30
    if name.lower() in ["zdt4","zdt6"]:
        num_dv = 10
    elif name.lower() == "constr" or name.lower() == "constroc":
        num_dv = 2
    elif name.lower() == "srn":
        num_dv = 2
    elif name.lower() == "sch":
        num_dv = 1
    elif name.lower() == "water":
        num_dv = 3
    elif name.lower() in ["rosen","ackley","rosenc"]:
        num_dv = 2
    elif name.lower() == "tkn":
        num_dv = 2
    elif name.lower() == "fon":
        num_dv = 3
    elif name.lower() == "kur":
        num_dv = 100
    elif name.lower() == "binh":
        num_dv = 2
    elif name.lower() == "hosaki":
        num_dv = 2 

    # write a generic template file for the dec vars
    tpl_file = "dv.dat.tpl".format(name)
    with open(os.path.join(test_d,tpl_file),'w') as f:
        f.write("ptf ~\n")
        for i in range(num_dv):
            f.write("dv_{0} ~ dv_{0}     ~\n".format(i))
    

    additive_chance_tpl_file = "additive_par.dat.tpl"
    with open(os.path.join(test_d,additive_chance_tpl_file),'w') as f:
        f.write("ptf ~\n")
        f.write("obj1_add_par ~   obj1_add_par   ~\n")
        f.write("obj2_add_par ~   obj2_add_par   ~\n")
        if name.lower() in ["srn","constr","tkn","rosenc","constroc"]:
            f.write("constr1_add_par ~   constr1_add_par   ~\n")
            f.write("constr2_add_par ~   constr2_add_par   ~\n")
        elif name.lower() == "water":
            [f.write("constr{0}_add_par ~   constr{0}_add_par   ~\n".format(i)) for i in range(7)]


    with open(os.path.join(test_d,additive_chance_tpl_file.replace(".tpl","")),'w') as f:
        f.write("obj1_add_par 0.0\n")
        f.write("obj2_add_par 0.0\n")
        if name.lower() in ["srn","constr","tkn","rosenc","constroc"]:
            f.write("constr1_add_par 0.0\n")
            f.write("constr2_add_par 0.0\n")
        elif name.lower() == "water":
            [f.write("constr{0}_add_par 0.0\n".format(i)) for i in range(7)]

    if risk_obj:
        risk_tpl_file = os.path.join(test_d,"risk.dat.tpl")
        with open(risk_tpl_file, 'w') as f:
            f.write("ptf ~\n")
            f.write("_risk_ ~   _risk_   ~\n")

    if self_adaptive:
        adaptive_tpl_file = os.path.join(test_d,"adaptive.dat.tpl")
        with open(adaptive_tpl_file, 'w') as f:
            f.write("ptf ~\n")
            f.write("_de_f_ ~   _de_f_   ~\n")
            f.write("_cr_ ~   _cr_   ~\n")
            f.write("_mr_ ~   _mr_   ~\n")

    with open(os.path.join(test_d,tpl_file.replace(".tpl","")),'w') as f:
        for i in range(num_dv):
            f.write("dv_{0} 0.5\n".format(i))

    # write a generic ins file for the two objectives
    ins_file = "obj.dat.ins".format(name)
    with open(os.path.join(test_d,ins_file),'w') as f:
        f.write("pif ~\n")
        if name.lower() == "water":
            for i in range(1,6):
                f.write("l1 w !obj_{0}!\n".format(i))
            for i in range(1,8):
                f.write("l1 w !const_{0}!\n".format(i))
        else:
            f.write("l1 w !obj_1!\n")
            if name.lower() not in ["rosen","ackley","rosenc","hosaki"]:
                f.write("l1 w !obj_2!\n")
            else:
                f.write("l1\n")
            if name.lower() in ["srn","tkn","rosenc","constroc"]:
                f.write("l1 w !const_1!\n")
                f.write("l1 w !const_2!\n")
        

    # now scape this python file to get the function lines and
    # the helper lines
    path_to_file = os.path.relpath(__file__)
    print(path_to_file)
    lines = open(path_to_file,'r').readlines()
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
        raise Exception(name + " function not found")
    helper_lines = []
    for i in range(len(lines)):
        if lines[i].startswith("def helper(".format(name)):
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

        f.write("if __name__ == '__main__':\n    helper({0})\n".format(name))

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
    par.loc[:,"parchglim"] = "relative"

    if name.lower() == "zdt4":
        par.loc[pst.par_names[0],"parubnd"] = 1.0
        par.loc[pst.par_names[0],"parlbnd"] = 0.0
        par.loc[pst.par_names[1:],"parubnd"] = 5.0
        par.loc[pst.par_names[1:],"parlbnd"] = -5.0

    elif name.lower() == "srn":
        par.loc[:,"parubnd"] = 20.0
        par.loc[:,"parlbnd"] = -20.0
        par.loc[:,"partrans"] = "none"
        par.loc[:,"parval1"] = 5.0

    elif name.lower() == "tkn":
        par.loc[:,"parubnd"] = np.pi - 1.0e-10
        par.loc[:,"parlbnd"] = 1.0e-10
        par.loc[:,"partrans"] = "none"
        par.loc[:,"parval1"] = 1.0

    elif name.lower() == "fon":
        par.loc[:,"parubnd"] = 4
        par.loc[:,"parlbnd"] = -4
        par.loc[:,"partrans"] = "none"
        par.loc[:,"parval1"] = 1.0

    elif name.lower() == "kur":
        par.loc[:,"parubnd"] = 5
        par.loc[:,"parlbnd"] = -5
        par.loc[:,"partrans"] = "none"
        par.loc[:,"parval1"] = 1.0

    elif name.lower() == "constr" or name.lower() == "constroc":
        par.loc["dv_0","parlbnd"] = 0.1
        par.loc["dv_0","parubnd"] = 1.0
        par.loc["dv_0","parval1"] = 0.5
        par.loc["dv_1","parlbnd"] = 0.0
        par.loc["dv_1","parubnd"] = 5.0
        par.loc["dv_1","parval1"] = 2.5
        if name.lower() == "constr":
            pst.prior_information = pst.null_prior
            pi = pst.prior_information
            pi.loc["const_1","pilbl"] = "const_1"
            pi.loc["const_1","equation"] = "9.0 * dv_0 + 1.0 * dv_1 = 6.0"
            pi.loc["const_1","weight"] = 1.0
            
            pi.loc["const_1","obgnme"] = "greater_than"
            pi.loc["const_2","pilbl"] = "const_2"
            pi.loc["const_2","equation"] = "9.0 * dv_0 - 1.0 * dv_1 = 1.0"
            pi.loc["const_2","weight"] = 1.0
            pi.loc["const_2","obgnme"] = "greater_than"
        
      
    elif name.lower() == "sch":
        par.loc[:,"parubnd"] = 1000.0
        par.loc[:,"parlbnd"] = -1000.0
        par.loc[:,"partrans"] = "none"
        par.loc[:,"parval1"] = 0.0

    elif name.lower() == "water":
        par.loc["dv_0","parlbnd"] = 0.01
        par.loc["dv_0","parubnd"] = 0.45
        par.loc["dv_0","parval1"] = 0.2
        par.loc["dv_1","parlbnd"] = 0.01
        par.loc["dv_1","parubnd"] = 0.1
        par.loc["dv_1","parval1"] = 0.05
        par.loc["dv_2","parlbnd"] = 0.01
        par.loc["dv_2","parubnd"] = 0.1
        par.loc["dv_2","parval1"] = 0.05

    elif name.lower() in ["rosen","ackley","rosenc"]:
        par.loc["dv_0","parlbnd"] = -4.
        par.loc["dv_0","parubnd"] = 4
        par.loc["dv_0","parval1"] = -1.
        par.loc["dv_1","parlbnd"] = -4.
        par.loc["dv_1","parubnd"] = 4
        par.loc["dv_1","parval1"] = -1.

    elif name.lower() == "binh":
        par.loc["dv_0","parlbnd"] = -10.
        par.loc["dv_0","parubnd"] = 10
        par.loc["dv_0","parval1"] = -1.
        par.loc["dv_1","parlbnd"] = -10.
        par.loc["dv_1","parubnd"] = 10
        par.loc["dv_1","parval1"] = -1.

    elif name.lower() == "hosaki":
        par.loc["dv_0","parlbnd"] = 0
        par.loc["dv_0","parubnd"] = 5
        par.loc["dv_0","parval1"] = 2.5
        par.loc["dv_1","parlbnd"] = 0
        par.loc["dv_1","parubnd"] = 5
        par.loc["dv_1","parval1"] = 2.5
    


    adf = pst.add_parameters(os.path.join(test_d,additive_chance_tpl_file),pst_path=".")
    print(adf)

    par = pst.parameter_data

    par.loc[adf.parnme,"parubnd"] = 0.5
    par.loc[adf.parnme,"parval1"] = 0.0
    par.loc[adf.parnme,"parlbnd"] = -0.5
    par.loc[adf.parnme,"parchglim"] = "relative"
    par.loc[adf.parnme,"pargp"] = "obj_add"
    #much less uncertainty in the second obj
    par.loc[adf.parnme[1],"parubnd"] = 0.5
    par.loc[adf.parnme[1],"parlbnd"] = -0.5

    if additive_chance:
        par.loc[adf.parnme, "partrans"] = "none"
        if name == "srn":
            par.loc[adf.parnme, "parubnd"] = 5.0
            par.loc[adf.parnme, "parlbnd"] = -5.0
            par.loc[adf.parnme, "parchglim"] = "relative"
            par.loc[adf.parnme, "pargp"] = "obj_add"
            # much less uncertainty in the second obj
            par.loc[adf.parnme[1], "parubnd"] = 2.5
            par.loc[adf.parnme[1], "parlbnd"] = -2.5
    else:
        par.loc[adf.parnme, "partrans"] = "fixed"


    pst.rectify_pgroups()
    pst.parameter_groups.loc["obj_add","inctyp"] = "absolute"

    if risk_obj:
        rdf = pst.add_parameters(risk_tpl_file,pst_path=".")
        par = pst.parameter_data
        par.loc[rdf.parnme,"partrans"] = "none"
        par.loc[rdf.parnme, "parubnd"] = 1.0
        par.loc[rdf.parnme, "parval1"] = 0.5
        par.loc[rdf.parnme, "parlbnd"] = 0.0
        par.loc[rdf.parnme, "parchglim"] = "relative"
        par.loc[rdf.parnme, "pargp"] = "decvars"
        pst.add_pi_equation(["_risk_"],"_risk_",obs_group="greater_than")

    if self_adaptive:
        addf = pst.add_parameters(adaptive_tpl_file,pst_path=".")
        par = pst.parameter_data
        par.loc[addf.parnme,"partrans"] = "none"
        par.loc[addf.parnme, "parchglim"] = "relative"
        par.loc[addf.parnme, "pargp"] = "decvars"
        par.loc["_de_f_", "parubnd"] = 1.0
        par.loc["_de_f_", "parval1"] = 0.8
        par.loc["_de_f_", "parlbnd"] = 0.5
        par.loc["_cr_", "parubnd"] = 1.0
        par.loc["_cr_", "parval1"] = 0.9
        par.loc["_cr_", "parlbnd"] = 0.8
        par.loc["_mr_", "parubnd"] = (1.0/par.shape[0]) * 1.25
        par.loc["_mr_", "parval1"] = (1.0/par.shape[0])
        par.loc["_mr_", "parlbnd"] = 0.01



    obs = pst.observation_data
    obs.loc[:,"weight"] = 1.0
    obs.loc[:,"obgnme"] = "less_than"
    #obs.loc[["obj_1","obj_2"],"obgnme"] = "less_than_obj" # all these zdt probs are min-min

    if name.lower() == "srn":
        obs.loc["const_1","obsval"] = 225
        obs.loc["const_2","obsval"] = -10
    elif name.lower() == "tkn":
        obs.loc["const_1","obsval"] = 0.0
        obs.loc["const_2","obsval"] = 0.5
    elif name.lower() == "water":
        obs.loc["const_1","obsval"] = 1
        obs.loc["const_2","obsval"] = 1
        obs.loc["const_3","obsval"] = 50000
        obs.loc["const_4","obsval"] = 16000
        obs.loc["const_5","obsval"] = 10000
        obs.loc["const_6","obsval"] = 2000
        obs.loc["const_7","obsval"] = 550
    elif name.lower() == "rosenc":
        obs.loc["const_1", "obsval"] = 0
        obs.loc["const_2", "obsval"] = 0
    elif name.lower() == "constroc":
        obs.loc["const_1", "obsval"] = 6
        obs.loc["const_2", "obsval"] = 1
        obs.loc["const_1","obgnme"] = "greater_than"
        obs.loc["const_2","obgnme"] = "greater_than"


           
    pst.pestpp_options["opt_dec_var_groups"] = "decvars"
    if name.lower() in ["rosen","rosenc","hosaki"]:
        pst.pestpp_options["mou_objectives"] = "obj_1"
    else:
        pst.pestpp_options["mou_objectives"] = "obj_1,obj_2"
    if name.lower() == "water":
        pst.pestpp_options["mou_objectives"] = "obj_1,obj_2,obj_3,obj_4,obj_5"

    if risk_obj:
        pst.pestpp_options["mou_objectives"] += ",_risk_"
        pst.pestpp_options["mou_risk_objective"] = True

    pst.model_command = "python forward_run.py"
    pst.control_data.noptmax = 0
    pst.write(os.path.join(test_d,name+".pst"))
    
    # run mou with noptmax = 0 to make sure we are getting something
    pyemu.os_utils.run("{0} {1}.pst".format(exe_path,name),cwd=test_d)
    pst = pyemu.Pst(os.path.join(test_d,name+".pst"))
    print(pst.phi)
    #if name.lower() in ["zdt1","zdt2","zdt3"]:
    #    assert pst.phi < 1.0e-10

    cov = pyemu.Cov.from_parameter_data(pst)
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst,cov,100)
    pe.enforce()
    pe.to_csv(os.path.join(test_d,"prior.csv"))
    return test_d

def run_problem_chance_external_fixed(test_case="zdt1"):
    assert "zdt" in test_case
    test_d = setup_problem(test_case,additive_chance=True)
    pst = pyemu.Pst(os.path.join(test_d,"{0}.pst".format(test_case)))
    par = pst.parameter_data
    par.loc["dv_9","partrans"] = "fixed"
    par.loc["obj1_add_par","partrans"] = "fixed"
    pst.control_data.noptmax = 2
    pst.pestpp_options["mou_population_size"] = 10
    pst.pestpp_options["panther_echo"] = True
    pst.pestpp_options["mou_generator"] = "de"
    pst.pestpp_options["panther_agent_freeze_on_fail"] = False
    pst.pestpp_options["opt_stack_size"] = 10
    pst.pestpp_options["opt_risk"] = 0.95
    pst.pestpp_options["opt_recalc_chance_every"] = 100
    pst.write(os.path.join(test_d,"{0}.pst".format(test_case)))
    #pyemu.os_utils.run("{0} {1}.pst".format(exe_path,test_case),cwd=test_d)
    master_d = test_d.replace("template","master_chance_external_fixed")
    pyemu.os_utils.start_workers(test_d, exe_path, "{0}.pst".format(test_case),
                                  num_workers=35, master_dir=master_d,worker_root=test_root,
                                  port=port)
    return master_d

def run_problem(test_case="zdt1",pop_size=100,noptmax=100,generator="de",env="nsga", self_adaptive=False):
    test_d = setup_problem(test_case,additive_chance=False,risk_obj=False,self_adaptive=self_adaptive)
    pst = pyemu.Pst(os.path.join(test_d,"{0}.pst".format(test_case)))
    pst.control_data.noptmax = noptmax
    pst.pestpp_options["mou_population_size"] = pop_size
    pst.pestpp_options["panther_echo"] = False
    pst.pestpp_options["mou_generator"] = generator
    pst.pestpp_options["mou_env_selector"] = env
    pst.pestpp_options["panther_agent_freeze_on_fail"] = True
    pst.write(os.path.join(test_d,"{0}.pst".format(test_case)))
    #pyemu.os_utils.run("{0} {1}.pst".format(exe_path,test_case),cwd=test_d)
    master_d = test_d.replace("template","master_generator={0}_env={1}_popsize={2}_risk=0.5_riskobj=False_adaptive={3}".
                              format(generator,env,pop_size,self_adaptive))
    pyemu.os_utils.start_workers(test_d, exe_path, "{0}.pst".format(test_case), 
                                  num_workers=35, master_dir=master_d,worker_root=test_root,
                                  port=port)
    
    #TODO: need some asserts here
    return master_d


def plot_results(master_d, sequence=False):
    plt_dir = os.path.join("mou_tests","test_result_plots")
    if not os.path.exists(plt_dir):
        os.mkdir(plt_dir)
    #master_d = os.path.join("mou_tests","zdt1_master")
    assert os.path.exists(master_d)

    #odf_files = [f for f in os.listdir(master_d) if f.endswith("obs_pop.csv") and "archive" in f]
    #odf_files_all = [f for f in os.listdir(master_d) if f.endswith("obs_pop.csv") and "archive" not in f]
    #odfs = [pd.read_csv(os.path.join(master_d,f),index_col=0) for f in odf_files]
    #odfs_all = [pd.read_csv(os.path.join(master_d,f),index_col=0) for f in odf_files_all]
    case = os.path.split(master_d)[1].split('_')[0]
    df = pd.read_csv(os.path.join(master_d,"{0}.pareto.summary.csv".format(case)))
    df_arc = pd.read_csv(os.path.join(master_d,"{0}.pareto.archive.summary.csv".format(case)))
    print(df_arc.columns)

    import matplotlib.colors as colors
    import matplotlib.pyplot as plt
    cols = ["obj_1","obj_2"]
    if "water" in master_d.lower():
        cols = ["obj_4","obj_1"]
    # cols = odfs[0].columns
    # mins,maxs = {},{}
    # for oname in cols:
    #     mins[oname] = 1.0e+30
    #     maxs[oname] = -1.0e+30
    #     for df in odfs_all:
    #         mins[oname] = min(mins[oname],df.loc[:,oname].min())
    #         maxs[oname] = max(maxs[oname], df.loc[:, oname].max())
    # print(mins,maxs)

    #colors = ["0.5",'m','b','c','g','y','r']
    #df = df.loc[df.generation<80,:]
    gens = df.generation.unique()
    gens.sort()
    print(gens)
    if sequence is False:
        gens = [gens[-2]]
    cmap = plt.get_cmap("jet",lut=len(gens))
    #for i,df in enumerate(odfs_all):
    #for i,gen in enumerate(gens):
    #    ax.scatter(df_arc.loc[df_arc.generation==gen,cols[0]],df_arc.loc[df_arc.generation==gen,cols[1]],marker=".",
    #        c=cmap(i/len(gens)),s=50,alpha=0.25)
    for igen,gen in enumerate(gens):
        fig, ax = plt.subplots(1, 1, figsize=(6, 6))
        if "_risk_" in df_arc.columns:
            rcmap = plt.get_cmap("jet")
            rvals = df_arc._risk_.unique()
            #df = df_arc.loc[df_arc._risk_>0.9,:]
            df = df_arc.copy()

            df = df.loc[df.generation == gen,:]
            print(df.shape)
            cb = ax.scatter(df.loc[:, cols[0]],
                       df.loc[:, cols[1]],
                       marker="+", c=df._risk_.values, s=100, label="final non dom solutions",
                       cmap="jet")

            cb = plt.colorbar(cb)
            cb.set_label("risk")
            ax.set_title("{0}, {1} generation shown, {2} members in archive". \
                         format(case, gen, df.shape[0]))
        else:
            if df_arc.shape[0] > 0:
                ax.scatter(df_arc.loc[df_arc.generation==gen,cols[0]],
                    df_arc.loc[df_arc.generation==gen,cols[1]],
                    marker="+",c='k',s=100,label="final non dom solutions")
            else:
                ax.scatter(df.loc[df.generation == gen, cols[0]],
                           df.loc[df.generation == gen, cols[1]],
                           marker="+", c='k', s=100, label="final non dom solutions")
            ax.set_title("{0}, {1} generation shown, {2} members in archive". \
                         format(case, gen, df_arc.shape[0]))

        possibles = globals().copy()
        possibles.update(locals())
        method = possibles.get(case)

        if "zdt" in master_d.lower():
            x0 = np.linspace(0,1,1000)
            o1,o2 = [],[]
            for xx0 in x0:
                x = np.zeros(30)
                x[0] = xx0
                ret_vals = method(x)
                o1.append(ret_vals[0][0])
                o2.append(ret_vals[0][1])

            ax.plot(o1,o2,"k",label="truth")
            ax.set_xlim(min(df_arc.loc[:,cols[0]].min(),min(o1)),
                max(df_arc.loc[:,cols[0]].max(),max(o1)))
            ax.set_ylim(min(df_arc.loc[:,cols[1]].min(),min(o2)),
                max(df_arc.loc[:,cols[1]].max(),max(o2)))
        else:    
            ax.set_xlim(df_arc.loc[:,cols[0]].min(),df_arc.loc[:,cols[0]].max())
            ax.set_ylim(df_arc.loc[:,cols[1]].min(),df_arc.loc[:,cols[1]].max())

        ax.legend()
        
        ax.set_xlabel("objective 1 (minimize)")
        ax.set_ylabel("objective 2 (minimize)")

        plt.tight_layout()
        if sequence:
            plt.savefig(os.path.join(plt_dir,os.path.split(master_d)[-1]+"_{0:03d}.png".format(igen)))
        else:
            plt.savefig(os.path.join(plt_dir, os.path.split(master_d)[-1] + "_{0:03d}.pdf".format(gen)))
        plt.close("all")
    if sequence:
        prefix = os.path.split(master_d)[-1]+"_%03d.png"
        pyemu.os_utils.run("ffmpeg -r 60 -y -f image2 -s 1920x1080 -i {0} -vcodec libx264 -crf 25  -pix_fmt yuv420p {1}.mp4".format(prefix,os.path.split(master_d)[-1]),
                           cwd=plt_dir)


def run_problem_chance(test_case="zdt1",pop_size=100,noptmax=10,stack_size=30,
                       chance_points="single",recalc=100,risk_obj=False,
                       risk=0.95,generator="de",env="nsga",self_adaptive=False):
    
    test_d = setup_problem(test_case,additive_chance=True, risk_obj=risk_obj, self_adaptive=self_adaptive)
    pst = pyemu.Pst(os.path.join(test_d,"{0}.pst".format(test_case)))
    pst.control_data.noptmax = noptmax
    pst.pestpp_options["mou_population_size"] = pop_size   
    pst.pestpp_options["opt_risk"] = risk
    pst.pestpp_options["opt_stack_size"] = stack_size
    pst.pestpp_options["opt_chance_points"] = chance_points
    pst.pestpp_options["panther_echo"] = False
    pst.pestpp_options["mou_generator"] = generator
    pst.pestpp_options["mou_env_selector"] = env
    pst.pestpp_options["opt_recalc_chance_every"] = recalc
    pst.write(os.path.join(test_d,"{0}.pst".format(test_case)))
    master_d = test_d.replace("template","master_generator={3}_env={4}_risk={0:04.3f}_riskobj={1}_chancepoints={2}_adaptive={5}".\
                              format(risk,risk_obj,chance_points,generator,env,self_adaptive))
    pyemu.os_utils.start_workers(test_d, exe_path, "{0}.pst".format(test_case), 
                                  num_workers=35, master_dir=master_d,worker_root=test_root,
                                  port=port)
    return master_d




def start_workers(case="zdt1"):
    pyemu.os_utils.start_workers(os.path.join("mou_tests","{0}_template".format(case)),
                                 exe_path, "{0}.pst".format(case),
                                  num_workers=35, worker_root="mou_tests",
                                  port=4004)

def run_single_obj_sch_prob(risk_obj):
    test_d = os.path.join("mou_tests","sch_template")
    setup_problem("sch",True,risk_obj=risk_obj)
    pst = pyemu.Pst(os.path.join(test_d,"sch.pst"))
    pst.pestpp_options["mou_objectives"] = "obj_1"
    if risk_obj:
        pst.pestpp_options["mou_objectives"] += ",_risk_"
    pst.observation_data.loc["obj_2","obsval"] = 2.0
    pst.write(os.path.join(test_d,"sch.pst"))
    pst.control_data.noptmax = 100
    pst.pestpp_options["mou_population_size"] = 50
    pst.pestpp_options["panther_echo"] = True
    pst.pestpp_options["mou_generator"] = "de"
    pst.pestpp_options["panther_agent_freeze_on_fail"] = True
    pst.pestpp_options["opt_chance_points"] = "all"
    pst.pestpp_options["opt_recalc_chance_every"] = 100
    pst.pestpp_options["opt_stack_size"] = 50
    pst.write(os.path.join(test_d,"sch.pst"))
    #pyemu.os_utils.run("{0} {1}.pst".format(exe_path,test_case),cwd=test_d)
    master_d = test_d.replace("template","master")
    pyemu.os_utils.start_workers(test_d, exe_path, "sch.pst", 
                                  num_workers=35, master_dir=master_d,worker_root=test_root,
                                  port=port)
    
    #TODO: need some asserts here
    return master_d


def plot_results_single(master_d):
    plt_dir = os.path.join("mou_tests","test_result_plots")
    if not os.path.exists(plt_dir):
        os.mkdir(plt_dir)
    #master_d = os.path.join("mou_tests","zdt1_master")
    assert os.path.exists(master_d)

    odf_files_all = [f for f in os.listdir(master_d) if f.endswith("obs_pop.csv") and "archive" not in f]
    odfs_all = [pd.read_csv(os.path.join(master_d,f),index_col=0) for f in odf_files_all]
    
    pdf_files_all = [f for f in os.listdir(master_d) if f.endswith("dv_pop.csv") and "archive" not in f]
    pdfs_all = [pd.read_csv(os.path.join(master_d,f),index_col=0) for f in pdf_files_all]
    
    import matplotlib.colors as colors
    import matplotlib.pyplot as plt
    
    possibles = globals().copy()
    possibles.update(locals())
    case = os.path.split(master_d)[-1].split("_")[0]
    method = possibles.get(case)

    fig,ax = plt.subplots(1,1,figsize=(6,6))
    cmap = plt.get_cmap("jet",lut=len(odfs_all))
    for i,(df,pdf) in enumerate(zip(odfs_all,pdfs_all)):
       #ax.scatter(pdf.loc[:,"dv_0"],pdf.loc[:,"dv_1"],
       # c=cmap(i/len(odfs_all)), marker=".",s=50,alpha=0.25)
       ax.scatter(pdf.loc[:,"dv_0"],pdf.loc[:,"dv_1"],
        c="k", marker=".",s=20,alpha=i/len(odfs_all))

    x = np.linspace(-2,2,1000)
    y = np.linspace(-2,2,1000)
    X,Y = np.meshgrid(x,y)
    
    z = []
    #for xx,yy in zip(X.flatten(),Y.flatten()):
    for xx in x:
        for yy in y:
            z.append(method((xx,yy))[0])
    Z = np.array(z).reshape((x.shape[0],x.shape[0]))
    ax.imshow(Z,alpha=0.5,extent=(-2,2,-2,2),interpolation="none")
    ax.legend()
    ax.set_xlim(-2,2.0)
    ax.set_ylim(-2,2.0)
    plt.tight_layout()
    plt.savefig(os.path.join(plt_dir,os.path.split(master_d)[-1]+".pdf"))
    plt.close("all")



def risk_compare_plot():
    plt_dir = os.path.join("mou_tests", "test_result_plots")
    if not os.path.exists(plt_dir):
        os.mkdir(plt_dir)
    master_d = os.path.join("mou_tests","zdt1_master_chance_0.950_True")
    assert os.path.exists(master_d)

    case = os.path.split(master_d)[1].split('_')[0]
    df = pd.read_csv(os.path.join(master_d, "{0}.pareto.summary.csv".format(case)))
    df_arc = pd.read_csv(os.path.join(master_d, "{0}.pareto.archive.summary.csv".format(case)))
    print(df_arc.columns)

    import matplotlib.colors as colors
    import matplotlib.pyplot as plt
    cols = ["obj_1", "obj_2"]

    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    # df = df.loc[df.generation<80,:]
    gens = df.generation.unique()
    gens.sort()
    print(gens)
    ax = axes[0]
    cmap = plt.get_cmap("jet", lut=len(gens))
    # for i,df in enumerate(odfs_all):
    # for i,gen in enumerate(gens):
    #    ax.scatter(df_arc.loc[df_arc.generation==gen,cols[0]],df_arc.loc[df_arc.generation==gen,cols[1]],marker=".",
    #        c=cmap(i/len(gens)),s=50,alpha=0.25)
    if "_risk_" in df_arc.columns:
        rcmap = plt.get_cmap("jet")
        rvals = df_arc._risk_.unique()
        # df = df_arc.loc[df_arc._risk_>0.9,:]
        df = df_arc.copy()
        print(df._risk_.values)
        df = df.loc[df.generation == gens[-1], :]
        cb = ax.scatter(df.loc[:, cols[0]],
                        df.loc[:, cols[1]],
                        marker="+", c=df._risk_.values, s=100, label="final non dom solutions",
                        cmap="jet")
        cb.set_label("risk")
        plt.colorbar(cb)
    else:
        ax.scatter(df_arc.loc[df_arc.generation == gens[-1], cols[0]],
                   df_arc.loc[df_arc.generation == gens[-1], cols[1]],
                   marker="+", c='k', s=100, label="final non dom solutions")

    ax.set_title("{0}, {1} generations shown, {2} members in archive, {3} total members". \
                 format(case, len(gens), df_arc.shape[0], df.shape[0]))
    ax.legend()

    ax = axes[1]
    master_ds = [os.path.join("mou_tests",m_d) for m_d in os.listdir("mou_tests")\
                 if "zdt1_master_chance_" in m_d and "False" in m_d]

    print(master_ds)
    for master_d in master_ds:
        r = float(master_d.split('_')[-2])
        case = os.path.split(master_d)[1].split('_')[0]
        df = pd.read_csv(os.path.join(master_d, "{0}.pareto.summary.csv".format(case)))
        df_arc = pd.read_csv(os.path.join(master_d, "{0}.pareto.archive.summary.csv".format(case)))
        gens = df.generation.unique()
        gens.sort()

        ax.scatter(df_arc.loc[df_arc.generation == gens[-1], cols[0]],
                   df_arc.loc[df_arc.generation == gens[-1], cols[1]],
                   marker="+", c=rcmap(r), s=100)

    possibles = globals().copy()
    possibles.update(locals())
    method = possibles.get(case)

    for ax in axes:
        x0 = np.linspace(0, 1, 1000)
        o1, o2 = [], []
        for xx0 in x0:
            x = np.zeros(30)
            x[0] = xx0
            ret_vals = method(x)
            o1.append(ret_vals[0][0])
            o2.append(ret_vals[0][1])

        ax.plot(o1, o2, "k", label="truth")

    # ax.set_xlim(-0.1,1.1)
    # ax.set_ylim(-0.1,7.0)
    plt.tight_layout()
    plt.savefig(os.path.join(plt_dir, "risk.pdf"))
    plt.close("all")

if __name__ == "__main__":
        
    #zdt1_test()
    # setup_zdt_problem("zdt1",30)
    # setup_zdt_problem("zdt2",30)
    # setup_zdt_problem("zdt3",30)
    # setup_zdt_problem("zdt4",10)
    # setup_zdt_problem("zdt6",10)
    #shutil.copy2(os.path.join("..","exe","windows","x64","Debug","pestpp-mou.exe"),os.path.join("..","bin","pestpp-mou.exe"))
    start_workers("zdt1")
    #setup_problem("rosenc")
