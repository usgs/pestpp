"""
This script is used to setup the Fonseca-Fleming problem.

Inputs devision variables: x1, x2
Objectives: obj1, obj2
Bounds: x1, x2 in [-4, 4]
"""
import os
import sys
import shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pyemu

template_dir = os.path.join("fonseca_fleming_demo", "FON_template")

def fon_setup():
    # Create a basic PST file with minimal required content
    with open(os.path.join(template_dir, "fon.pst"), "w") as f:
        f.write("pcf\n")
        f.write("* control data\n")
        f.write("restart estimation\n")
        f.write("2 2 1 0 1\n")
        f.write("1 1 single point 1 0 0\n")
        f.write("5.0 2.0 0.3 0.01 10\n")
        f.write("5.0 5.0 0.001\n")
        f.write("0.1\n")
        f.write("30 0.01 3 3 0.01 3\n")
        f.write("1 1 1\n")
        f.write("* parameter groups\n")
        f.write("* parameter data\n")
        f.write("* observation data\n")
        f.write("* model command line\n")
        f.write("* model input/output\n")
        f.write("fon.tpl fon.par\n")
        f.write("fon.ins fon.out\n")
    pst = pyemu.Pst(os.path.join(template_dir, "fon.pst"))

    pst.control_data.noptmax = 20

    par_data = []
    for i in range(2):
        par_name = f"x{i+1}"
        par_data.append(
            {"parnme": par_name, 
            "partrans": "none", 
            "parchglim": "relative",
            "pargp": "decvar", 
            "parval1": 0.0,  
            "parlbnd": -4.0,  
            "parubnd": 4.0, 
            "scale": 1.0, 
            "offset": 0.0}
        )
    pst.add_parameters(par_data)

    # Add observation groups for the two objectives
    pst.observation_groups = ["l_obj"]

    # Add the two objective functions as observations
    # For Fonseca-Fleming, we want to minimize both objectives
    obs_data = []
    for i in range(2):
        obs_name = f"obj{i+1}"
        obs_data.append(
            {"obsnme": obs_name,
            "obsval": 0.0,
            "weight": 1.0,
            "obgnme": "l_obj"}
        )
    pst.add_observations(obs_data)

    # Create a template file for the parameters
    with open(os.path.join(template_dir, "fon.tpl"), "w") as f:
        f.write("ptf !\n")
        f.write("!    x1    !\n")
        f.write("!    x2    !\n")

    # Create an instruction file to read the objectives
    with open(os.path.join(template_dir, "fon.ins"), "w") as f:
        f.write("pif #\n")
        f.write("l1 \n")
        f.write("l1 #,#!obj1!\n")
        f.write("l1 #,#!obj2!\n")

    # Set up the model command
    pst.model_command = f"python forward_run.py"

    # Set up the template and instruction files
    pst.template_files = ["fon.tpl"]
    pst.input_files = ["fon.par"]
    pst.instruction_files = ["fon.ins"]
    pst.output_files = ["fon.out"]

    pst.pestpp_options["opt_dec_var_groups"] = "decvar"
    pst.pestpp_options["mou_objectives"] = "obj1,obj2"
    pst.pestpp_options["panther_agent_freeze_on_fail"] = "True"
    pst.pestpp_options["mou_save_population_every"] = 1
    pst.pestpp_options["mou_population_size"] = 20
    pst.pestpp_options["mou_generator"] = "pso"
    pst.pestpp_options["opt_iter_tol"] = 0.0001
    pst.pestpp_options["mou_max_generations"] = 50

    # Write the PST file
    pst_file = os.path.join(template_dir, "fon.pst")
    pst.write(pst_file)

    print(f"PST file created: {pst_file}")

if __name__ == "__main__":
    fon_setup()
