# %%
"""
# Data Assimilation for One-dimensional GW flow model
Data assimilation combine model predictions and observations. In this notebook, the model will be setup first and then it will be used to generate a set  of synthatic observations.
"""

# %%
import time
start_time = time.time()
import os, sys
import numpy as np
import pandas as pd
import configparser
import pyemu
import matplotlib.pyplot as plt
import inspect
import shutil
import flopy

# %%
"""
## Config Parameters
In this exercise, we will use multiple user configuration parameters. To make things clean, we used the a configuration file to hold all these parameters.

"""

# %%
config = configparser.ConfigParser()
config.read("config.init")
fidw = open("config.init", 'r')
settings = fidw.read()
fidw.close()
print(settings)

# %%
"""
## Setup workspace
Create a template folder that will contain model executable, model input files, and model output files. It is 
called "template" because the worker will be generated from it.
"""

# %%
template_ws = os.path.join(".", config.get('workspace', 'template_folder')) # od_flow os one dim flow model
if not(os.path.isdir(template_ws)):
    os.mkdir(template_ws)
else:
    shutil.rmtree(template_ws)
print(os.path.abspath(template_ws))

# %%
"""
## Model description
A Python function "generate_1d_model" is used to setup a 1D GW flow model as well as misc files that we will be used and showed in this notebook. The following are the main features in the model:- 
* The model consists of one layer, one row, and 100 columns. 
* The model generates a one-row model where a general head boundary (GHB = -10 m) exist at the upstream and another on GHB downstream (GHB = -20). We should interpret the negative number as depth below groundsurface level.
* The model consists of 24 stress periods, and each stress period has 30 time steps.
* Recharge is applied uniformly in space and variably in time at all cells. 
* In the middle cell a groundwater pumping is applied. 

Let us generate the model and look more at it's component. For detailed in formation about model setup, please see "generate_1d_model" in file model_setup. Look at output at the end of the simulation.
"""
def csv_to_tpl(csv_file, name_col, par_col, tpl_file ):
    df = pd.read_csv(csv_file)

    fidw = open(tpl_file, 'w')
    fidw.write("ptf ~\n")

    for i, col in enumerate(df.columns):
        if i == 0:
            csv_header = str(col)
            continue
        csv_header = csv_header + "," + str(col)
    csv_header = csv_header + "\n"
    fidw.write(csv_header)

    line = ""
    for irow, row in df.iterrows():
        row[par_col] =  " ~   {0}    ~".format(row[name_col])
        line = ",".join(row.astype(str).values.tolist()) + "\n"
        fidw.write(line)

    fidw.close()

def csv_to_ins(csv_file, name_col, obs_col, ins_file):
    df = pd.read_csv(csv_file)
    part1 = "l1"

    for col in df.columns:
        if col in [obs_col]:
            break
        part1 = part1 + " ~,~"

    obs_names = df[name_col]
    fidw = open(ins_file, 'w')
    fidw.write("pif ~\n")
    fidw.write("l1\n") # header

    for irow, row in df.iterrows():
        line = part1 + "   !{0}!    ~".format(row[name_col])
        line = line + ",~\n"

        fidw.write(line)

    fidw.close()


def generate_random_ensemble(dx = np.array([1]* 100),
                             dy = np.array([1])
                             , corr_scale = 10, N = 100 ):
    v = pyemu.utils.geostats.ExpVario(a=corr_scale,contribution=1.0)
    gs = pyemu.utils.geostats.GeoStruct(variograms=v,nugget=0.0)

    ss = pyemu.utils.geostats.SpecSim2d(dx,dy,gs)
    arrays = ss.draw_arrays(num_reals=N)
    arrays = arrays.squeeze().T
    return arrays

# %%
#2) create a simple flow model
from model_setup import generate_1d_model
figs = generate_1d_model() # this is a simple script that setup the 1D flow problem. 


# %%
"""
The GW pumping causes a depression cone to develop and grow with time.
"""

# %%
#%matplotlib inline
figs[0]

# %%
figs[1]

# %%
# Copy the generated model to the template folder
print(os.listdir())
shutil.copytree('model_dataset', os.path.join(template_ws, 'model_dataset'))
shutil.rmtree('model_dataset')
shutil.copy('..\param_utils.py', os.path.join(template_ws, 'param_utils.py'))
shutil.copy('..\obs_utils.py', os.path.join(template_ws, 'obs_utils.py'))

# %%
"""
### Simplifying Model I/O
We write a Python warapper that takes one file for each parameter group and generate one output file as output. Let us look at the model input files. and output file is similar to observation data section. This is just a personal preference! you may still deals with the native MODFLOW I/O files.
* input_dis.csv
* hk_1.dat
* iheads_1.dat
"""

# %%


# %%
"""
The model script (model_setup.py) has a function called "forward_model" that will work as our model.
"""

# %%
import model_setup
src_forward_run = inspect.getsource(model_setup.forward_run)
print(src_forward_run)

# %%
# let us assemble the script that works as forward_model
script = """
import os, sys
import shutil
import numpy as np
import pandas as pd
import flopy
import matplotlib.pyplot as plt
import pyemu
import obs_utils, param_utils

{}
if __name__ == "__main__":
     forward_run()
 
""".format(src_forward_run)
print(script)

# %%
# Let us write the model
with open(os.path.join(template_ws, r"forward_model.py"), 'w') as fidw:
    fidw.write(script)

# %%
fnam = os.path.join(template_ws, r'model_dataset\flow_1d.nam')
mf = flopy.modflow.Modflow.load(os.path.basename(fnam), model_ws= os.path.dirname(fnam))
fig = plt.figure(figsize=(15, 5))
ax = fig.add_subplot(1, 1, 1)
xsect = flopy.plot.PlotCrossSection(model=mf, line={'Row': 0})
linecollection = xsect.plot_grid()
t = ax.set_title('Cross-Section - Model Grid')

# %%


# %%
np.random.seed(int(config.get('da_ensemble', 'seed_number' )))
               
corr_scale = 20.0*float(config.get('da_ensemble', 'correlation_length' ))
nx, ny = mf.nrow, mf.ncol # of cells
N = int(config.get('da_ensemble', 'number_of_realizations' )) # number of realizations
#m = int(config.get('da_ensemble', 'number_of_realizations' ) # number of measurements
refence_realization = 4
delc = mf.dis.delc.array * 20 # cells must be square to use geostat simulator!
delr = mf.dis.delr.array

# pyemu allow us to define the spatial structure of grid. This object (sr) will be used later.
sr = pyemu.helpers.SpatialReference(delr=delr, delc=delc )


# %%
# Before running the model. Let us see something else. The trasnient data table.
trans_data = pd.read_csv(os.path.join(template_ws, r".\model_dataset\temporal_param.csv" ))
trans_data

# %%
base_folder = os.getcwd()
os.chdir(template_ws)
cmd = sys.executable
cmd = cmd + " " + "forward_model.py"
os.system(cmd)
os.chdir(base_folder)

# %%
out_df = pd.read_csv(os.path.join(template_ws, r".\model_dataset\output_sim.csv" ))
out_df

# %%
# We can change the start time and end time and run the model. Since our model has 24 stress period,
# let us divide it into 4 cycles, each has 6 stress period
df_input = pd.read_csv(os.path.join(template_ws, r"model_dataset\input_dis.csv" ))
df_input.loc[df_input['parname']=='start_sp', 'parval'] = 0
df_input.loc[df_input['parname']=='end_sp', 'parval'] = 5
df_input.to_csv(os.path.join(template_ws, r"model_dataset\input_dis.csv" ))
base_folder = os.getcwd()
os.chdir(template_ws)
cmd = sys.executable
cmd = cmd + " " + "forward_model.py"
os.system(cmd)
os.chdir(base_folder)

# %%
"""
# PstFrom to setup the problem.
"""

# %%
# Generate a pf object
# new_d =  os.path.join(".", config.get('workspace', 'worker_folder'))
# pf = pyemu.utils.PstFrom(original_d=template_ws, new_d=new_d,
#                  remove_existing=True,
#                  longnames=True, spatial_reference=sr,
#                  zero_based=False,start_datetime="1-1-2021")

# %%
pst_files = {}
par_fn = os.path.join(template_ws,r"model_dataset\input_dis.csv")
df_par1 = pd.read_csv(par_fn)
pst_files[os.path.basename(par_fn)] = os.path.basename(par_fn) + ".tpl"
csv_to_tpl(csv_file=par_fn, par_col='parval', name_col='parname', tpl_file= par_fn+".tpl")
parnames = []
parnames = parnames + df_par1['parname'].values.tolist()

# pyemu.utils.write_list_tpl(filenames=par_fn, dfs = df_par1, name= 'dis',
#                            tpl_filename= par_fn + ".tpl", index_cols = 'parval')
# pf.add_parameters(filenames=par_fn, use_cols=["parval"], par_type='grid', index_cols=['parname'],
#                   par_name_base="dis_", pargp="dis_sps", ult_ubound=50, ult_lbound=1e-3,
#                   par_style='direct', transform='log')

# %%
# add input file that contains parameters
par_fn =  os.path.join(template_ws,r"model_dataset\hk.dat")
arr = np.loadtxt(par_fn)
nms = ["k_{}".format(i) for i in range(len(arr))]
pyemu.utils.simple_tpl_from_pars(parnames=nms, tplfilename=par_fn + ".tpl")
parnames = parnames + nms
# par_name_base = 'k'
# pargp = 'KH'
#
# v = pyemu.utils.geostats.ExpVario(a=corr_scale,contribution=1.0)
# grid_gs = pyemu.utils.geostats.GeoStruct(variograms=v,nugget=0.0)

# pf.add_parameters(filenames=par_fn, par_type="grid",
#                   par_name_base=par_name_base, pargp=pargp,
#                   upper_bound=3., lower_bound=-3, ult_ubound=3.0, ult_lbound=-3,
#                   geostruct=grid_gs, transform='none', par_style='direct')

# %%
# add input file that contains parameters
par_fn =  os.path.join(template_ws,r"model_dataset\iheads.dat")
arr = np.loadtxt(par_fn)
nms = ["h_ini_{}".format(i) for i in range(len(arr))]
pyemu.utils.simple_tpl_from_pars(parnames=nms, tplfilename=par_fn + ".tpl")
parnames = parnames + nms

#
# par_name_base = 'h'
# pargp = 'Hs'
#
# v = pyemu.utils.geostats.ExpVario(a=corr_scale,contribution=1.0)
# grid_gs = pyemu.utils.geostats.GeoStruct(variograms=v,nugget=0.0)
#
# pf.add_parameters(filenames=par_fn, par_type="grid",
#                   par_name_base=par_name_base, pargp=pargp,
#                   upper_bound=0, lower_bound=-50, ult_ubound= 0, ult_lbound=-50,
#                   geostruct=grid_gs, transform='none', par_style='direct')

# %%
output_file = r".\model_dataset\output_sim.csv"
#output_file = os.path.join(template_ws, output_file)
#output_file = r"D:\Workspace\projects\mississippi\manuscripts\pestppda_paper\notebooks\one_dim_flow_enkf\flow1d_worker\model_dataset\output_sim.csv"

fn_out = os.path.join(template_ws, output_file)
df_out = pd.read_csv(fn_out)
ins_file =  fn_out + ".ins"
# odf = pyemu.pst_utils.csv_to_ins_file(csv_filename = os.path.join(template_ws, output_file),
#                                 ins_filename = os.path.join(template_ws, output_file)+".ins",
#                                 only_cols = "simval" , prefix='')
csv_to_ins(csv_file = fn_out, name_col = 'obsnme', obs_col = 'simval',
           ins_file = ins_file)

pst = pyemu.Pst.from_par_obs_names(par_names= parnames,
                                   obs_names= df_out['obsnme'].values)


#
# obs_prefix = 'h' ## Note: using capital letter will cause issues
# obsgp = 'heads'
# hds_df = pf.add_observations(output_file,insfile=ins_file,index_cols="obsnme",
#                     use_cols= ['simval'],prefix=obs_prefix, obsgp = obsgp, ofile_sep = ",")
#
#

# %%
# pf.mod_sys_cmds.append(cmd)
# pst = pf.build_pst(version=2)


ihead = np.loadtxt(os.path.join(template_ws, "model_dataset\iheads.dat"))

pe_k = generate_random_ensemble(dx = np.array([1]* 100),
                             dy = np.array([1])
                             , corr_scale = 20, N = 50 )
pe_h = np.zeros_like(pe_k) + ihead[:,np.newaxis]
pe = np.vstack([pe_k, pe_h])

nms1 = ["k_{}".format(i) for i in range(100)]
nms2 = ["h_ini_{}".format(i) for i in range(100)]

pe = pd.DataFrame(pe.T, columns=nms1+nms2)
pe.to_csv(os.path.join(template_ws, 'HKensemble.csv'))
## Select one realization to represent the unknown truth (the actual reference field). From this truth realization choose observations.
x_true = pe.iloc[refence_realization].values

# %%


# %%
plt.figure()
plt.plot(pe.T.values)
plt.show()

ens =  pe.values.T
y_ref = ens[:, refence_realization]

plt.figure()
plt.plot(ens, color = [0.7,0.7,0.7], zorder=1)
plt.plot(ens[:,refence_realization])
plt.show()

# %%
# this a synthatic problem, so we need observations. The observation is generated from the reference realization
true_parameters = ens[:,refence_realization]
df_input = pd.read_csv(os.path.join(template_ws, r"model_dataset\input_dis.csv" ))
df_input.loc[df_input['parname']=='start_sp', 'parval'] = 0
df_input.loc[df_input['parname']=='end_sp', 'parval'] = 23
df_input.to_csv(os.path.join(template_ws, r"model_dataset\input_dis.csv" ))

hk = np.loadtxt(os.path.join(template_ws, "model_dataset\hk.dat"))

hk = np.power(10,true_parameters[0:100])
np.savetxt(os.path.join(template_ws, "model_dataset\hk.dat"), hk )

base_folder = os.getcwd()
os.chdir(template_ws)
cmd = sys.executable
cmd = cmd + " " + "forward_model.py"
os.system(cmd)
os.chdir(base_folder)
df_truth_outupt = pd.read_csv(os.path.join(template_ws, "model_dataset\output_sim.csv"))
xx = 1
# pyemu.pst_utils.write_to_template(pst.parameter_data.parval1_trans,
#                                   os.path.join(pst_path, tpl_file),
#                                   os.path.join(pst_path, in_file))

actual_obs_for_cycles = df_truth_outupt[df_truth_outupt['obgnme']=='simHO'].copy()
actual_obs_for_cycles["totim"] = actual_obs_for_cycles['obsnme'].str.split("_").str[-1].astype(int)
totims = np.sort(actual_obs_for_cycles["totim"].unique())
totims = totims.reshape(4,6)
for cycle in range(4):
    ttims = totims[cycle]
    actual_obs_for_cycles.loc[actual_obs_for_cycles['totim'].isin(ttims), 'cycle'] = cycle

## use observation from the reference
pst.observation_data['obsval'] = 0.0
pst.observation_data['weight'] = 1000 # this is 1/std
pst.observation_data['cycle'] = -1

obsnames = actual_obs_for_cycles[actual_obs_for_cycles['cycle']==0]['obsnme'].values

obs_cycle = pd.DataFrame(index = obsnames, columns = range(0,4))
weight_cycle = pd.DataFrame(index = obsnames, columns = range(0,4))
par_cycle = pd.DataFrame(index = ['start_sp', 'end_sp' ], columns = range(0,4))
sp_strt = 0
sp_end = 5
for cycle in range(4):
    vals = actual_obs_for_cycles[actual_obs_for_cycles['cycle'] == cycle]['simval'].values
    obs_cycle[cycle] = vals
    weight_cycle[cycle] = 1000.0
    par_cycle[cycle] = [sp_strt, sp_end]
    sp_strt = sp_end + 1
    sp_end = sp_strt + 5

pst.observation_data['cycle'] = -1
pst.observation_data.loc[pst.observation_data['obsnme'].str.contains("hf_"), 'obgnme'] = "Hf"
pst.observation_data.loc[pst.observation_data['obsnme'].str.contains("hf_"), 'weight'] = 0
pst.observation_data.loc[pst.observation_data['obsnme'].str.contains("qghb_"), 'obgnme'] = "GHBQ"
pst.observation_data.loc[pst.observation_data['obsnme'].str.contains("qghb_"), 'weight'] = 0
#pst.observation_data['state_par_link']

pst.parameter_data.loc[pst.parameter_data['parnme'].str.contains("h_ini_"), 'pargp'] = "Hi"
pst.parameter_data.loc[pst.parameter_data['parnme'].str.contains("k_"), 'pargp'] = "KH"
pst.parameter_data.loc[pst.parameter_data['parnme'].str.contains("_sp"), 'pargp'] = "DIS"
dyn_par = pst.parameter_data.loc[pst.parameter_data['pargp'].isin(['Hi']), 'parnme'].values
pst.observation_data.loc[pst.observation_data['obgnme'].isin(["Hf"]), 'state_par_link'] =dyn_par


pst.parameter_data['cycle'] = -1
pst.parameter_data['parchglim'] = 'relative'

pst.model_input_data['pest_file'] = ['input_dis.csv.tpl', 'iheads.dat.tpl', 'hk.dat.tpl']
pst.model_input_data['model_file'] = ['input_dis.csv', 'iheads.dat', 'hk.dat']
pst.model_input_data['cycle'] = -1
pst.model_output_data['cycle'] = -1

pst.model_command = cmd

pst.svd_data.eigthresh =  1e-5
pst.pestpp_options['da_add_base']= False
pst.pestpp_options['da_parameter_ensemble']= 'HKensemble.csv'
#pst.pestpp_options['DA_SUBSET_SIZE'] = 50
pst.pestpp_options['da_num_reals'] = N
pst.pestpp_options['ies_init_lam'] = [1]
pst.pestpp_options['ies_lambda_mults'] = 1
pst.pestpp_options['lambda_scale_fac'] = 1
pst.control_data.noptmax = 1

# %%
#pst.model_command[0] = pst.model_command[0].replace('python', sys.executable)
pst.write("flow1d.pst", version=2)

# %%
"""
## Important check locations of ins and templpe
"""

# %%
shutil.copy2(os.path.join(vs_code_dir,"exe","windows","x64","Debug","pestpp-da.exe"),os.path.join(new_d,"pestpp-da.exe"))

# %%
import subprocess
base_folder = os.getcwd()
os.chdir(new_d)
argv = ["pestpp-da.exe", os.path.basename(pf.pst.filename)]
def run_command(cmd):
    with subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True) as p:
        while True:
            line = p.stdout.readline()
            if not line:
                break
            print(line)    
        exit_code = p.poll()
    return exit_code


run_command(argv)

os.chdir(base_folder)

# %%
posterior_par = pd.read_csv(os.path.join(new_d, template_ws+".global.0.pe.csv"))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(ens, color = [0.7,0.7,0.7], zorder=1, label = 'Prior K')
handles, labels = ax.get_legend_handles_labels()
display = []
display.append(len(labels))

ax.plot(ens[:,refence_realization], 'g', zorder = 3)
handles, labels = ax.get_legend_handles_labels()
display.append(len(labels))

ax.scatter(y_index, y_ref[y_index], color = 'r', zorder = 3)
handles, labels = ax.get_legend_handles_labels()
display.append(display[-1] + len(labels))

del(posterior_par['real_name'])
ax.plot(posterior_par.values.T, color = 'b', label = 'Posterior K')
handles, labels = ax.get_legend_handles_labels()
display.append(display[-1] + len(labels))

ax.legend([handle for i,handle in enumerate(handles) if i in display],
      [label for i,label in enumerate(labels) if i in display], loc = 'best')


# %%
## Try more iterations for 

# %%
pst.pestpp_options['da_num_reals'] = N
pst.pestpp_options['DA_SUBSET_SIZE'] = 10
pst.pestpp_options['ies_init_lam'] = [10]
pst.pestpp_options['ies_lambda_mults'] = [0.1, 1, 10]
pst.pestpp_options['lambda_scale_fac'] = 1
pst.control_data.noptmax = 5 # number of iterations
pst.write(os.path.abspath(pf.pst.filename), version=2)

# %%
import subprocess
base_folder = os.getcwd()
os.chdir(new_d)
argv = ["pestpp-da.exe", os.path.basename(pf.pst.filename)]
def run_command(cmd):
    with subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True) as p:
        while True:
            line = p.stdout.readline()
            if not line:
                break
            print(line)    
        exit_code = p.poll()
    return exit_code


run_command(argv)

os.chdir(base_folder)

# %%
posterior_par = pd.read_csv(os.path.join(new_d, template_ws+".global.0.pe.csv"))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(ens, color = [0.7,0.7,0.7], zorder=1, label = 'Prior K')
handles, labels = ax.get_legend_handles_labels()
display = []
display.append(len(labels))

ax.plot(ens[:,refence_realization], 'g', zorder = 3)
handles, labels = ax.get_legend_handles_labels()
display.append(len(labels))

ax.scatter(y_index, y_ref[y_index], color = 'r', zorder = 3)
handles, labels = ax.get_legend_handles_labels()
display.append(display[-1] + len(labels))

del(posterior_par['real_name'])
ax.plot(posterior_par.values.T, color = 'b', label = 'Posterior K')
handles, labels = ax.get_legend_handles_labels()
display.append(display[-1] + len(labels))

ax.legend([handle for i,handle in enumerate(handles) if i in display],
      [label for i,label in enumerate(labels) if i in display], loc = 'best')

# %%
sys.path.append(r"..")
import rec_util 
import rec_util
ws = r"D:\Workspace\projects\mississippi\manuscripts\pestppda_paper\notebooks\one_dim_flow\new_od_flow_template"
fname = r"od_flow_template.rec"

fname = os.path.join(ws, fname)
rec = rec_util.RecFile(fname= os.path.join(new_d, template_ws+".rec"))



# %%
