# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 21:15:54 2023

@author: mac732
"""

import pandas as pd
import os
import shutil
import glob
import pyemu
from distutils.dir_util import copy_tree

#update outer repository
outer_dirlist = [os.path.basename(x) for x in glob.glob("./outer_[0-999]*", recursive=True)]
outer_dirlist = sorted(outer_dirlist, key = lambda x: int(x.split("_")[1]))

prev_dv_file = glob.glob(os.path.join(".", outer_dirlist[-2], "*archive.dv_pop.csv"), recursive= True)
prev_dv = pd.read_csv(prev_dv_file[0])
prev_obs_file = glob.glob(os.path.join(".", outer_dirlist[-2], "*archive.obs_pop.csv"), recursive= True)
prev_obs = pd.read_csv(prev_obs_file[0])

curr_dv_file = glob.glob(os.path.join(".", outer_dirlist[-1], "*0.dv_pop.csv"), recursive= True)
curr_dv = pd.read_csv(curr_dv_file[0])
curr_obs_file = glob.glob(os.path.join(".", outer_dirlist[-1], "*0.obs_pop.csv"), recursive= True)
curr_obs = pd.read_csv(curr_obs_file[0])

merged_dv_file = pd.concat([prev_dv, curr_dv], ignore_index= True)
merged_obs_file = pd.concat([prev_obs, curr_obs], ignore_index= True)

merged_dv_file.to_csv(os.path.join(".", "template", "update_repo", "merged.dv_pop.csv"), index = False)
merged_obs_file.to_csv(os.path.join(".", "template", "update_repo", "merged.obs_pop.csv"), index = False)

shutil.copytree(os.path.join("outer_rundir", "template"), "temp")
copy_tree(os.path.join("template", "update_repo"), "temp")

#update training dataset
path_in = os.path.join("template", "model", "input")
infill_data = pd.concat([curr_dv, curr_obs], axis=1)

training_data = pd.read_csv(os.path.join(path_in, "trainingdata.csv"))
training_data.columns = infill_data.columns.values
training_data = pd.concat([training_data, infill_data], ignore_index= True)
training_data.to_csv(os.path.join(path_in, "trainingdata.csv"))

#run pestpp mou in pareto sorting mode
os.chdir("./temp")
os.system("pestpp-mou outer_repo.pst")
os.chdir((".."))

outer_repo_sumlist = glob.glob(os.path.join("temp", "outer_repo.pareto*"), recursive=True)
for name in outer_repo_sumlist:
    shutil.copy(name, outer_dirlist[-1])

#generate restart files for next set of inner iters
pop_size = pd.read_csv(glob.glob(os.path.join("template", "model", "input", "*.dv_pop.csv"))[0]).shape[0]

dv_restart = pd.read_csv(os.path.join("temp", "outer_repo.archive.dv_pop.csv"))
obs_restart = pd.read_csv(os.path.join("temp", "outer_repo.archive.obs_pop.csv"))
shutil.rmtree("temp")

curr_dv = curr_dv[~curr_dv['real_name'].isin(dv_restart['real_name'].values)]
curr_dv_subsample = curr_dv.sample(n=pop_size-dv_restart.shape[0])
curr_obs_subsample = curr_obs[curr_obs['real_name'].isin(curr_dv_subsample['real_name'])]

dv_restart = pd.concat([dv_restart, curr_dv_subsample], ignore_index= True)
dv_restart = dv_restart.sort_values('real_name')
obs_restart = pd.concat([obs_restart, curr_obs_subsample], ignore_index= True)
obs_restart = obs_restart.sort_values('real_name')

sd_init = pd.read_csv(os.path.join("template", "model", "preproc", "sd_aug.csv"))
obs_restart.to_csv(os.path.join(path_in,"mogp_henry.obs_pop.csv"), index=False)
dv_restart.to_csv(os.path.join(path_in,"mogp_henry.dv_pop.csv"), index=False)


