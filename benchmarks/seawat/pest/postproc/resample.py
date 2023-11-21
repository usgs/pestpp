# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 21:42:06 2023

@author: mac732
"""

import pandas as pd
import os
import glob

max_infill = 100

#get inner and outer pareto data
outer_dirlist = [os.path.basename(x) for x in glob.glob("./outer_[0-999]*", recursive=True)]
outer_dirlist = sorted(outer_dirlist, key = lambda x: int(x.split("_")[1]))
outer_pareto_file = glob.glob(os.path.join(outer_dirlist[-1],"*.pareto.archive.summary.csv"), recursive=True)
outer_pareto = pd.read_csv(outer_pareto_file[0])

inner_dirlist = [os.path.basename(x) for x in glob.glob("./mogp_[0-999]*", recursive=True)]
inner_dirlist = sorted(inner_dirlist, key = lambda x: int(x.split("_")[1]))
inner_pareto_file = glob.glob(os.path.join(inner_dirlist[-1], "pest", "*.pareto.archive.summary.csv"), recursive=True)
inner_pareto = pd.read_csv(inner_pareto_file[0])
        
nmax_inner = max(inner_pareto.generation)

inner_pareto = inner_pareto[inner_pareto["generation"] == nmax_inner]
inner_pareto = inner_pareto[~inner_pareto['member'].isin(outer_pareto['member'].values)]

maxppf = max(inner_pareto['nsga2_front'])

#get all dv data
all_dv = pd.DataFrame()
csvfiles = glob.glob(os.path.join(inner_dirlist[-1], "pest",'*[0-999].dv_pop.csv'), recursive=True)
csvfiles = sorted(csvfiles, key = lambda x: int(x.split(".dv")[0].split(".")[1]))
for name in csvfiles:
    inner_dv = pd.read_csv(name)
    all_dv = pd.concat([all_dv, inner_dv], ignore_index= True)

#perform infill selection
n_infill = 0
ppd_front = 1
infill_pop = pd.DataFrame(columns = inner_pareto.columns.values)
while ((n_infill < max_infill) & (ppd_front <= maxppf)):
    infill_sort = inner_pareto[inner_pareto['nsga2_front']==ppd_front]
    infill_sort = infill_sort[~infill_sort['member'].isin(infill_pop['member'])]
    infill_sort.drop_duplicates(subset = 'member', inplace = True)
    
    if infill_sort.shape[0]>0:
        infill_sort = infill_sort.sort_values(by='ehvi', ascending = False)
        infill_add = infill_sort.head(max_infill-n_infill)
        infill_pop = pd.concat([infill_pop, infill_add], ignore_index=True)
        n_infill = infill_pop.shape[0]
    else:
        ppd_front += 1

infill_dv_pop = all_dv[all_dv['real_name'].isin(infill_pop['member'].values)]

#save to outer iter dir
infill_dv_pop.to_csv(os.path.join(os.path.join("outer_rundir","template"),"infill.dv_pop.csv"), index=False)
infill_dv_pop.to_csv(os.path.join(os.path.join("outer_rundir","master"),"infill.dv_pop.csv"), index=False)

