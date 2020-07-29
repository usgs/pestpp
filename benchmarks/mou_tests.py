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
    return x[0],return (1 + x[1]) / x[0]

def srn(x):
    const1 = np.power(x[0], 2) + np.power(x[1], 2)  # lest than or equal to 225
    const2 =  3 * x[1] - x[0]  # greater than or equal to 10
    f1 = np.power(x[0] - 2, 2) + np.power(x[1] - 1, 2) + 2
    f2 = 9 * x[0] - np.power(x[1] - 1, 2)
    return const1, const2, f1, f2


def setup_zdt_problem(name,num_dv):
    tpl_file = "dv.dat.tpl".format(name)
    with open(os.path.join("{0}_template".format(name)tpl_file,'w') as f:
        f.write("ptf ~\n")
        for i in range(num_dv):
            f.write("dv_{1} ~ {0}_dv_{1}     ~\n".format(name, i))
    ins_file = "obj.dat.ins".format(name)
    with open(os.path.join("{0}_template".format(name),ins_file,'w') as f:
        f.write("pif ~\n")
        f.write("l1 !obj_1!\n")
        f.write("l2 !obj_2!\n")

    
    

    with open(os.path.join("{0}_template".format(name),"forward_run.py"),'w') as f:


if __name__ == "__main__":
        
    zdt1_test()
