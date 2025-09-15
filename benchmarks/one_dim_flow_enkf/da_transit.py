# %%
"""
# Data Assimilation for One-dimensional GW flow model
Data assimilation combines model predictions and observations. In this notebook, the model will be setup first and then it will be used to generate a set  of synthatic observations.
"""

# %%
import time
start_time = time.time()
import os, sys
import numpy as np
import pandas as pd
import configparser
sys.path.insert(0, r"D:\Workspace\Codes\pyemu_myfork\pyemu")
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
# Model description
A Python function "generate_1d_model" is used to setup a 1D GW flow model as well as misc files that we will be used and shown in this notebook. The following are the main features in the model:- 
* The model consists of one layer, one row, and 100 columns. 
* The model generates a one-row model where a general head boundary (GHB = -10 m) exist at the upstream and another on GHB downstream (GHB = -20). We should interpret the negative number as depth below groundsurface level.
* The model consists of 24 stress periods, and each stress period has 30 time steps.
* Recharge is applied uniformly in space and variably in time at all cells. 
* In the middle cell a groundwater pumping is applied. 

Let us generate the model and look more at it's component. For detailed information about model setup, please see "generate_1d_model" in file model_setup. Let Look us at output at the end of the simulation.


"""

# %%
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
We write a Python warapper that takes one file for each parameter group and generate one output file as output. Let us look at the model input files. Output file is similar to observation data section. This is just a personal preference! you may still deals with the native MODFLOW I/O files.
* input_dis.csv
* hk_1.dat
* iheads_1.dat
"""

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
"""
## Problem dimensions
"""

# %%

np.random.seed(int(config.get('da_ensemble', 'seed_number' )))               
corr_scale = 20.0*float(config.get('da_ensemble', 'correlation_length' ))
nx, ny = mf.nrow, mf.ncol # of cells
N = int(config.get('da_ensemble', 'number_of_realizations' )) # number of realizations
#m = int(config.get('da_ensemble', 'number_of_realizations' ) # number of measurements
refence_realization = 4
delc = mf.dis.delc.array * 20 # cells must be square to use geostat simulator!
delr = mf.dis.delr.array



# %%


# %%
"""
Let us run the model and look at output file
"""

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
"""
In sequential Data Assimilation, we need to be able to rerun the model starting from any arbitrary time point and stop at any other time point . The Python wrapper, that works as our model, uses a simple csv input file with two parameters: start time and end time. In this example we are going to divide the 24 stress periods into 4 cycles. 
"""

# %%
# We can change the start time and end time and run the model. Since our model has 24 stress period,
# let us divide it into 4 cycles, each has 6 stress period
df_input = pd.read_csv(os.path.join(template_ws, r"model_dataset\input_dis.csv" ))
df_input.loc[df_input['parname']=='start_sp', 'parval'] = 0
df_input.loc[df_input['parname']=='end_sp', 'parval'] = 5
df_input.to_csv(os.path.join(template_ws, r"model_dataset\input_dis.csv" ))
df_input

# %%
"""
Additionally, we have a master file that holds all transient information that define what pumping (Qw), recharge (Rch), and reference heads at boundary conditions, that are used. In MODFLOW, we also we need to define the stress period length and number of time steps in each stress period. In a more complex model setup, we might need to define a file groundwater pumping; in this case this csv file might contain a pumping file name for each stress period.  
"""

# %%
trans_data = pd.read_csv(os.path.join(template_ws, r".\model_dataset\temporal_param.csv" ))
trans_data

# %%


# %%
# now let us run the model using the new start and end times
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
"""
## Setup the PEST files
### Input/Template files 
We have three files..
"""

# %%
# ---------------------------
#(1) input_dis.csv File
# ---------------------------
input_file_list = [] # a container for all pst files input/output and tpl/ins file
tpl_file_list = []
parnames = []
obsnames = []

par_fn = os.path.join(template_ws,r"model_dataset\input_dis.csv")
df_par1 = pd.read_csv(par_fn)

# add the file to pst files container
input_file_list.append(r"model_dataset\input_dis.csv")
tpl_file_list.append(r"model_dataset\input_dis.csv"+".tpl")

#generate template file by converting csv file into simple tpl file
csv_to_tpl(csv_file=par_fn, par_col='parval', name_col='parname', tpl_file= par_fn+".tpl")

# we will keep adding params to this list
parnames = parnames + df_par1['parname'].values.tolist()

# see the template file
f = open(par_fn+".tpl", 'r')
print(f.read())
f.close()


# %%
# ---------------------------
#(2) hk.dat File: the hydraulic conductivity field
# ---------------------------
par_fn =  os.path.join(template_ws,r"model_dataset\hk.dat")
arr = np.loadtxt(par_fn)

nms = ["k_{}".format(i) for i in range(len(arr))]

#simple temple
pyemu.utils.simple_tpl_from_pars(parnames=nms, tplfilename=par_fn + ".tpl")

# add the file to pst files container
input_file_list.append(r"model_dataset\hk.dat")
tpl_file_list.append(r"model_dataset\hk.dat"+".tpl")

# we will keep adding params to this list
parnames = parnames + nms

f = open(par_fn+".tpl", 'r')
print(f.read())
f.close()

# %%
# ---------------------------
#(3) iheads.dat File: the initial hydraulic head field
# ---------------------------
par_fn =  os.path.join(template_ws,r"model_dataset\iheads.dat")
arr = np.loadtxt(par_fn)

nms = ["h_ini_{}".format(i) for i in range(len(arr))]
pyemu.utils.simple_tpl_from_pars(parnames=nms, tplfilename=par_fn + ".tpl")

# add the file to pst files container
input_file_list.append(r"model_dataset\iheads.dat")
tpl_file_list.append(r"model_dataset\iheads.dat"+".tpl")

# we will keep adding params to this list
parnames = parnames + nms

f = open(par_fn+".tpl", 'r')
print(f.read())
f.close()

# %%
"""
### output/instruction files 
"""

# %%
output_file_list = [] # a container for all pst files input/output and tpl/ins file
ins_file_list = []
output_file = r".\model_dataset\output_sim.csv"

fn_out = os.path.join(template_ws, output_file)
df_out = pd.read_csv(fn_out)
ins_file =  fn_out + ".ins"

# add the file to pst files container
output_file_list.append(r".\model_dataset\output_sim.csv")
ins_file_list.append(r".\model_dataset\output_sim.csv" + ".ins")

obsnames = df_out['obsnme'].values
#pyemu.utils.simple_ins_from_obs(obsnames=obsnames, insfilename= os.path.join(template_ws, r".\model_dataset\output_sim.dat") + ".ins")
csv_to_ins(csv_file = fn_out, name_col = 'obsnme', obs_col = 'simval',
           ins_file = ins_file)




# %%



# %%
"""
## Generate prior parameters/states ensemble
Notice that we have three input files, but only two has adjustable parameters/states: hk.dat and iheads.dat. The file input_dis.csv has data that only needed to restart the model. All adjustable states/parameters need prior ensembles
"""

# %%
# (1) Prior K ensemble
pe_k = generate_random_ensemble(dx = np.array([1]* 100),
                             dy = np.array([1])
                             , corr_scale = 20, N = 50 )
print(pe_k)

# %%
# (2) Prior initial head ensemble is assumed here the same. 
ihead = np.loadtxt(os.path.join(template_ws, "model_dataset\iheads.dat"))
pe_h = np.zeros_like(pe_k) + ihead[:,np.newaxis]
print(pe_h)

# %%
# (3) The augmented forecast matrix can be assembled now...
pe = np.vstack([pe_k, pe_h])

nms1 = ["k_{}".format(i) for i in range(100)]
nms2 = ["h_ini_{}".format(i) for i in range(100)]

pe = pd.DataFrame(pe.T, columns=nms1+nms2)
pe.to_csv(os.path.join(template_ws, 'HKensemble.csv'))
x_true = pe.iloc[refence_realization].values

# %%
# (4) Before moving on, let us look visualize the prior ensemble
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
"""
# Observations
As mesnioned before, DA combine model with observations. This is just a synthatic problem, so we need to generate synthatic observations. The observation is generated by simulating a truth parameters (reference realization), which is a a user defined in the init file.

Let us run the model using the assumed synthatic reality and produce our measurements...
"""

# %%
true_parameters = ens[:,refence_realization]

# (1) We need observations at the 24 stress periods.
df_input = pd.read_csv(os.path.join(template_ws, r"model_dataset\input_dis.csv" ))
df_input.loc[df_input['parname']=='start_sp', 'parval'] = 0
df_input.loc[df_input['parname']=='end_sp', 'parval'] = 23
df_input.to_csv(os.path.join(template_ws, r"model_dataset\input_dis.csv" ))

# (2) write input file using true K parameters. 
hk = true_parameters[0:100] # notice our model take the log as input
np.savetxt(os.path.join(template_ws, "model_dataset\hk.dat"), hk )

# (3) run the model
base_folder = os.getcwd()
os.chdir(template_ws)
cmd = sys.executable
cmd = cmd + " " + "forward_model.py"
os.system(cmd)
os.chdir(base_folder)
df_truth_outupt = pd.read_csv(os.path.join(template_ws, "model_dataset\output_sim.csv"))


# %%
df_truth_outupt

# %%
# (4) Assign each observation to a cycle number. In this case study we divide the 24 stress periods into 4 cycles, each has 6 stress periods.
actual_obs_for_cycles = df_truth_outupt[df_truth_outupt['obgnme']=='simHO'].copy()
actual_obs_for_cycles["totim"] = actual_obs_for_cycles['obsnme'].str.split("_").str[-1].astype(int)
totims = np.sort(actual_obs_for_cycles["totim"].unique())
totims = totims.reshape(4,6)
for cycle in range(4):
    ttims = totims[cycle]
    actual_obs_for_cycles.loc[actual_obs_for_cycles['totim'].isin(ttims), 'cycle'] = cycle
actual_obs_for_cycles

# %%
"""
## Cycle tables
At each cycle we need to define the observations to be assimilated, theirs weights, and non-adjustable parameters to used that are needed to restart the model to simulate each cycle. Cycle tables are a nice option to so

"""

# %%
obsnames1 = actual_obs_for_cycles[actual_obs_for_cycles['cycle']==0]['obsnme'].values

obs_cycle = pd.DataFrame(index = obsnames1, columns = range(0,4))
weight_cycle = pd.DataFrame(index = obsnames1, columns = range(0,4))
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
    
obs_cycle.to_csv(os.path.join(template_ws, 'head_obs_cycle.csv'))
weight_cycle.to_csv(os.path.join(template_ws, 'head_weight_cycle.csv'))
par_cycle.to_csv(os.path.join(template_ws, 'par_dis_cycle.csv'))

# %%
"""
# PESTPP-DA Control Files
"""

# %%
## (1) Generate PESTPP-DA control files files
pst = pyemu.Pst.from_par_obs_names(par_names= parnames,
                                   obs_names= obsnames)


# %%
# let us take a look at the observation data in PST
pst.observation_data

# %%
# (2) All observation was assigned to one observation group. Let us fix that. Also let us assign all observations weight to 1000 for now. 
pst.observation_data['obsval'] = 0.0 # let put 0.0 here because observation will be read each cycle from the cycle table
pst.observation_data['weight'] = 1000 # this is 1/std
pst.observation_data['cycle'] = -1   # this means that observation will be used at all cycles

# why are we having three observation groups? do not we have one observation group, which is the head. Yes we do, but we also need to extract dynamic states and
# forecast parameters. Now we have three groups (1) measured heads, (2) initial heads, and (3) forecast flow at the downstream GHB

# (3) Dynamic States (initial head)
pst.observation_data.loc[pst.observation_data['obsnme'].str.contains("hf_"), 'obgnme'] = "Hf"
pst.observation_data.loc[pst.observation_data['obsnme'].str.contains("hf_"), 'weight'] = 0 # Notice that we assign the weight to be zero

# (4) Forecast states that we are interseted in (simulated outflow at the downstream GHB)
pst.observation_data.loc[pst.observation_data['obsnme'].str.contains("qghb_"), 'obgnme'] = "GHBQ"
pst.observation_data.loc[pst.observation_data['obsnme'].str.contains("qghb_"), 'weight'] = 0
pst.observation_data

# %%
# let us look at parameters. We can see that the parameter datadrame is populated with default values. Let us fix that
pst.parameter_data

# %%
# (5-a) Assign group names for parameters. We have three groups: initial head, conductivity, and dis parameters.
pst.parameter_data.loc[pst.parameter_data['parnme'].str.contains("h_ini_"), 'pargp'] = "Hi"
pst.parameter_data.loc[pst.parameter_data['parnme'].str.contains("k_"), 'pargp'] = "KH"
pst.parameter_data.loc[pst.parameter_data['parnme'].str.contains("_sp"), 'pargp'] = "DIS"
pst.parameter_data

# %%
# (5-b) let us make sure that the none adjustable parameters are fixed. The only adjustable parameters is "K". DIS parameters will only be used when model is restarted
# between cycles. 
pst.parameter_data['partrans'] = 'fixed'
pst.parameter_data['cycle'] = -1
pst.parameter_data['parchglim'] = 'relative'
pst.parameter_data.loc[pst.parameter_data['parnme'].str.contains("k_"), 'partrans'] = "none"
pst.parameter_data.loc[pst.parameter_data['parnme'].str.contains("h_ini_"), 'partrans'] = "none" # check if this is needed??
pst.parameter_data.loc[pst.parameter_data['parnme'].str.contains("h_ini_"), 'parlbnd'] = -100
pst.parameter_data.loc[pst.parameter_data['parnme'].str.contains("h_ini_"), 'parubnd'] = 100

# let us also assig lower and upper bounds for K
pst.parameter_data.loc[pst.parameter_data['parnme'].str.contains("k_"), 'parlbnd'] = -4
pst.parameter_data.loc[pst.parameter_data['parnme'].str.contains("k_"), 'parubnd'] = 4
pst.parameter_data

# %%
# (6) The dynamic states is model output that will be used as initial conditions when the model is restarted. 
#  So we need to be able read it from model output and then write it as input for the next cycle IC. This is the reason we need to add the
# dynamic state to both observation data and parameter data. But we need to link both of them... See below.
dyn_par = pst.parameter_data.loc[pst.parameter_data['pargp'].isin(['Hi']), 'parnme'].values

# now link dynamic state in both parameters and observation dataframes
pst.observation_data.loc[pst.observation_data['obgnme'].isin(["Hf"]), 'state_par_link'] =dyn_par

pst.observation_data

# %%
# (7) Let PESTPP-DA be aware of other files
# (7-a) add cycke tables
pst.pestpp_options['da_observation_cycle_table'] = 'head_obs_cycle.csv'
pst.pestpp_options['da_weight_cycle_table'] = 'head_weight_cycle.csv'
pst.pestpp_options['DA_PARAMETER_CYCLE_TABLE'] =  'par_dis_cycle.csv'


# %%
# (7-b) input-template files
# pst.input_files = input_file_list
# pst.output_files = output_file_list
# pst.template_files = tpl_file_list
# pst.instruction_files = ins_file_list
pst.model_input_data['pest_file'] =  tpl_file_list 
pst.model_input_data['model_file'] = input_file_list
pst.model_input_data['cycle'] = -1

# (7-c) output-ins files
df = pd.DataFrame(columns = pst.model_input_data.columns)
df['pest_file'] = ins_file_list
df['model_file'] = output_file_list
df['cycle'] = -1
pst.model_output_data = df


# (7-d) Prior ensemble 
pst.pestpp_options['da_parameter_ensemble']= 'HKensemble.csv'

# %%

pst.model_command = cmd

pst.svd_data.eigthresh =  1e-5
pst.pestpp_options['da_add_base']= False

#pst.pestpp_options['DA_SUBSET_SIZE'] = 50
pst.pestpp_options['da_num_reals'] = N
pst.pestpp_options['ies_init_lam'] = [1]
pst.pestpp_options['ies_lambda_mults'] = 1
pst.pestpp_options['lambda_scale_fac'] = 1
pst.control_data.noptmax = 1

# %%
"""

"""

# %%


# %%
#pst.model_command[0] = pst.model_command[0].replace('python', sys.executable)
pst_file = config.get('workspace', 'basename') + '.pst'
pst.write(os.path.join(template_ws,pst_file), version=2)

# %%
"""
## Run...
"""

# %%
pst_exe = config.get('da_settings', 'pestpp_da_exe')
shutil.copy2(pst_exe,os.path.join(template_ws,"pestpp-da.exe"))

# %%
import subprocess
base_folder = os.getcwd()
os.chdir(template_ws)
argv = ["pestpp-da.exe", os.path.basename(r'flow1d.pst')]
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
col_to_plot = []
for col in pe.columns:
    if "k_" in col:
        col_to_plot.append(col)
ens = pe[col_to_plot]

for cycle in range(4):
    posterior_par = pd.read_csv(os.path.join(template_ws, "flow1d.global.{}.pe.csv".format(cycle)))
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(ens[col_to_plot].values.T, color = [0.7,0.7,0.7], zorder=1, label = 'Prior K')
    handles, labels = ax.get_legend_handles_labels()
    display = []
    display.append(len(labels))

    ax.plot(ens.values.T[:,refence_realization], 'r', zorder = 3)
    handles, labels = ax.get_legend_handles_labels()
    display.append(len(labels))


    del(posterior_par['real_name'])
    ax.plot(posterior_par[col_to_plot].values.T, color = 'b', label = 'Posterior K')
    handles, labels = ax.get_legend_handles_labels()
    display.append(display[-1] + len(labels))

    ax.legend([handle for i,handle in enumerate(handles) if i in display],
          [label for i,label in enumerate(labels) if i in display], loc = 'best')



# %%
sys.path.append(r"..")
import rec_util
recfile = os.path.join(template_ws, config.get('workspace', 'basename')+ ".rec" )
rec = rec_util.RecFile(fname= recfile)
phis = rec.phi
phis[phis['type'] == 'actual'][['cycle', 'iter_no', 'mean']]

# %%
# we can also see changes in parameters ensembel
plt.figure()
pars =rec.pars
pars = pars[pars['group'] == 'kh']
plt.plot(pars['cycle'], pars['init_CV'], label = 'initial CV' )
plt.plot(pars['cycle'], pars['curr_CV'], label = 'Current CV' )

# plot kh group

x = 1