import os
import sys
import shutil
import platform
import numpy as np
import pandas as pd
import platform
import pyemu
sys.path.append(os.path.join("..","benchmarks"))

import opt_test_suite_helper as mou_suite_heler

bin_path = os.path.join("test_bin")
if "linux" in platform.platform().lower():
    bin_path = os.path.join(bin_path,"linux")
elif "darwin" in platform.platform().lower() or "macos" in platform.platform().lower():
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




def test_setup_and_three_iters():
    cases = ["zdt1"]#,"water","constr","zdt2","zdt3","zdt4","zdt6","sch","srn","ackley","rosen","tkn"]
    noptmax = 3
    for case in cases:
        print("\n\n\n\n\n-----------------------------------------")
        print("                 {0}                  ".format(case))
        print("-----------------------------------------\n\n\n\n")

        #t_d = setup_problem(case)
        m_d = mou_suite_heler.run_problem(case,noptmax=noptmax)
        arc_file = os.path.join(m_d,"{0}.pareto.summary.csv".format(case))
        assert os.path.exists(arc_file), arc_file
        arc_df = pd.read_csv(arc_file,index_col=0)
        assert arc_df.shape[0] > 0, case
        return
        if "zdt" not in case:
            continue

        #fosm
        m_d = mou_suite_heler.run_problem_chance(case, noptmax=noptmax, pop_size=10, chance_points="single", recalc=200,
                                 stack_size=0)
        arc_file = os.path.join(m_d, "{0}.pareto.summary.csv".format(case))
        assert os.path.exists(arc_file),arc_file
        arc_df = pd.read_csv(arc_file, index_col=0)
        assert arc_df.shape[0] > 0, case

        # stack with all point and full reuse
        m_d = mou_suite_heler.run_problem_chance(case,noptmax=noptmax,pop_size=10,chance_points="all",recalc=200,
                                 stack_size=10)
        arc_file = os.path.join(m_d,"{0}.pareto.summary.csv".format(case))
        assert os.path.exists(arc_file),arc_file
        arc_df = pd.read_csv(arc_file,index_col=0)
        assert arc_df.shape[0] > 0, case

        # stack with single point and full reuse
        m_d = mou_suite_heler.run_problem_chance(case,noptmax=noptmax,pop_size=10,chance_points="single",recalc=200,
                                 stack_size=10)
        arc_file = os.path.join(m_d,"{0}.pareto.summary.csv".format(case))
        assert os.path.exists(arc_file), arc_file
        arc_df = pd.read_csv(arc_file,index_col=0)
        assert arc_df.shape[0] > 0, case

        # stack with single point, full reuse and risk obj
        m_d = mou_suite_heler.run_problem_chance(case,noptmax=noptmax,pop_size=10,chance_points="single",recalc=200, 
                                 risk_obj=True, stack_size=10)
        arc_file = os.path.join(m_d,"{0}.pareto.summary.csv".format(case))
        assert os.path.exists(arc_file), arc_file
        arc_df = pd.read_csv(arc_file,index_col=0)
        assert arc_df.shape[0] > 0, case

        # stack with single point, recalc every iter
        m_d = mou_suite_heler.run_problem_chance(case,noptmax=noptmax,pop_size=10,chance_points="single",recalc=1,
                                 stack_size=10)
        arc_file = os.path.join(m_d,"{0}.pareto.summary.csv".format(case))
        assert os.path.exists(arc_file), arc_file
        arc_df = pd.read_csv(arc_file,index_col=0)
        assert arc_df .shape[0] > 0, case


        # a restart test with fixed dec vars
        m_d = mou_suite_heler.run_problem_chance_external_fixed(case)
        arc_file = os.path.join(m_d, "{0}.pareto.summary.csv".format(case))
        assert os.path.exists(arc_file), arc_file
        arc_df = pd.read_csv(arc_file, index_col=0)
        assert arc_df.shape[0] > 0, case




if __name__ == "__main__":
    #shutil.copy2(os.path.join("..","exe","windows","x64","Debug","pestpp-mou.exe"),os.path.join("..","bin","pestpp-mou.exe"))
    test_setup_and_three_iters()
    