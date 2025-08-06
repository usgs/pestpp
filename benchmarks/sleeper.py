import numpy as np 
import time
import pyemu

def persistence_big_test():
    model_d = "ies_10par_xsec"
    base_d = os.path.join(model_d, "template")
    new_d = os.path.join(model_d, "timeout_template")
    if os.path.exists(new_d):
        shutil.rmtree(new_d)
    shutil.copytree(base_d, new_d)
    print(platform.platform().lower())
    pst = pyemu.Pst(os.path.join(new_d, "pest.pst"))
    with open(os.path.join(new_d,"run.py"),'w') as f:
        f.write("import os\nimport numpy as np\nimport time\nimport pyemu\npyemu.os_utils.run('mfnwt 10par_xsec.nam')\n")
        f.write("minutes = np.random.uniform(0.1,10)\nstart = time.time()\n")
        f.write("print(f'sleeping for {minutes:0.2f} minutes')")
        f.write("time.sleep(minutes*60)\nprint(f'wake up! Its been {(time.time()-start)/60} minutes')\n")
    pst.model_command = "python run.py"
    oe_file = os.path.join(new_d, "pest.0.obs.csv")
    if os.path.exists(oe_file):
        os.remove(oe_file)
    pst.control_data.noptmax = -1
    pst.pestpp_options["panther_persistent_workers"] = False
    pst.pestpp_options["ies_bad_phi_sigma"] = 2.5
    pst.pestpp_options["ies_num_reals"] = 200
    pst.pestpp_options["ensemble_output_precision"] = 40
    pst.pestpp_options["panther_master_timeout_milliseconds"] = 1000
    pst.pestpp_options["panther_agent_freeze_on_fail"] = False
    pst.write(os.path.join(new_d, "pest.pst"))
