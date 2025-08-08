import platform
import numpy as np 
import pyemu
import os, shutil 
import pathlib as pl

IP_ADDRESS = None # set this to a real one to run

def persistence_big_test():
    """
    test for panther_persistent_workers -> requires running on 
    a large cluster with Master noted in the IP_ADDRESS arg above.
    
    This test is left out of 
    """
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
        f.write("print(f'sleeping for {minutes:0.2f} minutes')\n")
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
    shutil.copy2('../bin/pestpp-ies',os.path.join(new_d,'pestpp-ies'))
    shutil.copy2('./test_bin/linux/mfnwt',os.path.join(new_d,'mfnwt'))
    
    # make condor files for the persistence test 
  
    wkdir = pl.Path("./ies_10par_xsec")
    if (wkdir / 'log').exists():
        shutil.rmtree(wkdir / 'log')
    (wkdir / 'log').mkdir()
    os.system(f"tar --use-compress-program='pigz' -cf {str(wkdir / 'timeout_template.tar')} {str(wkdir / 'timeout_template')}")
    with open(wkdir / 'worker.sh', 'w') as ofp:
        ofp.write("#!/bin/sh\n")
        ofp.write("tar xzf timeout_template.tar\n")
        ofp.write("mkdir -p tmp\n")
        ofp.write("tar -xzf tmp.tar.gz -C tmp\n")
        ofp.write("cd ies_10par_xsec/timeout_template\n")
        ofp.write("export PATH=../tmp/bin:$PATH \n")
        ofp.write(f"./pestpp-ies pest.pst /h {IP_ADDRESS}:4242\n")
    os.system(f"chmod +x {str(wkdir/'worker.sh')}")
    with open(wkdir / 'test.sub', 'w') as ofp:
        ofp.write("notification=Never\n")
        ofp.write("executable = worker.sh\n")
        ofp.write("universe=Vanilla\n")
        ofp.write("should_transfer_files   = YES\n")
        ofp.write("transfer_output_files = ""\n")
        ofp.write("transfer_input_files = timeout_template.tar, tmp.tar.gz\n")
        ofp.write("stream_output = True\n")
        ofp.write("stream_error = True\n")
        ofp.write("output                  = log/$(Cluster)_out.$(Process)\n")
        ofp.write("error                   = log/$(Cluster)_err.$(Process)\n")
        ofp.write("log                     = log/$(Cluster).log\n")
        ofp.write("queue 200\n")

       
                
        
if __name__=="__main__":
    persistence_big_test()
