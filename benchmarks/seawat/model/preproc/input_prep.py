# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 11:10:36 2023

@author: rqmac
"""

import pandas as pd
import os
import glob

outer_dirlist = [os.path.basename(x) for x in glob.glob("./outer_[0-999]*", recursive=True)]
outer_dirlist = sorted(outer_dirlist, key = lambda x: int(x.split("_")[1]))

dv_file = glob.glob(os.path.join(".", outer_dirlist[-1], "*0.dv_pop.csv"), recursive= True)
obs_file = glob.glob(os.path.join(".", outer_dirlist[-1], "*0.obs_pop.csv"), recursive= True)

path_in = os.path.join("template", "model", "input")
if outer_dirlist[-1].endswith("_0"):
    if os.path.exists(os.path.join(path_in, "trainingdata.csv")):
        os.remove( os.path.join(path_in, "trainingdata.csv"))
     
data_x = pd.read_csv(dv_file[0])
data_y = pd.read_csv(obs_file[0])

infill_data = pd.concat([data_x,data_y],axis=1)
infill_data.to_csv(os.path.join(path_in, "trainingdata.csv"), index = False)

sd_init = pd.read_csv(os.path.join("template", "model", "preproc", "sd_aug.csv"))
obs_pop = pd.concat([data_y,sd_init], axis=1)
obs_pop.to_csv(os.path.join(path_in,"mogp_henry.obs_pop.csv"), index=False)

dv_pop = data_x
dv_pop.to_csv(os.path.join(path_in,"mogp_henry.dv_pop.csv"), index=False)