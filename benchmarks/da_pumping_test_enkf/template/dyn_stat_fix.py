import os, sys
import pandas as pd


fn_par = r"D:\Workspace\projects\mississippi\pure_c_pestpp\pestpp\benchmarks\da_pumping_test_enkf\template\test.par_data.csv"
fn_obs = r"D:\Workspace\projects\mississippi\pure_c_pestpp\pestpp\benchmarks\da_pumping_test_enkf\template\test.obs_data.csv"
fn_ensemble = r"D:\Workspace\projects\mississippi\pure_c_pestpp\pestpp\benchmarks\da_pumping_test_enkf\template\kh_ensemble0.csv"
par_df = pd.read_csv(fn_par)
obs_df = pd.read_csv(fn_obs)
ens_df = pd.read_csv(fn_ensemble)

dyn_par_names = par_df.loc[par_df['pargp']=='heads','parnme']
dyn_par_names = dyn_par_names.astype(str) + "_dyn"
par_df.loc[par_df['pargp'] == 'heads', 'parnme'] = dyn_par_names

dyn_obs_names = obs_df.loc[obs_df['obgnme']=='heads', 'obsnme']
old_obs_names = dyn_obs_names.copy()
dyn_obs_names = dyn_obs_names.astype(str) + "_dyn"
obs_df.loc[obs_df['obgnme'] == 'heads', 'state_par_link'] = dyn_obs_names
ens_columns = ens_df.columns
for field in ens_columns:
    if field in old_obs_names.values:
        new_field = field + "_dyn"
        ens_df[new_field] = ens_df[field]
        del(ens_df[field])


par_df.to_csv(r"D:\Workspace\projects\mississippi\pure_c_pestpp\pestpp\benchmarks\da_pumping_test_enkf\template\test.par_data_dyn.csv")
obs_df.to_csv(r"D:\Workspace\projects\mississippi\pure_c_pestpp\pestpp\benchmarks\da_pumping_test_enkf\template\test.obs_data_dyn.csv")
ens_df.to_csv(r"D:\Workspace\projects\mississippi\pure_c_pestpp\pestpp\benchmarks\da_pumping_test_enkf\template\kh_ensemble_dyn.csv")



# template
fn = r"D:\Workspace\projects\mississippi\pure_c_pestpp\pestpp\benchmarks\da_pumping_test_enkf\template\stats_head.tpl"
fid = open(fn, 'r')
content = fid.readlines()
fid.close()

for iline , line in enumerate(content):
    if iline == 0:
        continue
    par_name = content[iline].strip().replace("~", "").strip()
    new_par_name = par_name + "_dyn"
    content[iline] = content[iline].replace(par_name, new_par_name)
fid = open(fn, 'w')
for line in content:
    fid.write(line)
fid.close()
xx = 1

