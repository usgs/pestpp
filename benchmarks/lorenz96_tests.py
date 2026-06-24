import copy
import os, sys
import pandas as pd
import shutil
import platform
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D
#from da_engine import Analysis
import pyemu
import lorenz96_model_setup
import inspect


bin_path = os.path.join("test_bin")
if "linux" in platform.platform().lower():
    bin_path = os.path.join(bin_path,"linux")
elif "darwin" in platform.platform().lower():
    bin_path = os.path.join(bin_path,"mac")
else:
    bin_path = os.path.join(bin_path,"win")

bin_path = os.path.abspath("test_bin")
os.environ["PATH"] += os.pathsep + bin_path


bin_path = os.path.join("..","..","bin")
exe = ""
if "windows" in platform.platform().lower():
    exe = ".exe"
exe_path = os.path.join(bin_path, "pestpp-da" + exe)


port = 4021

dim = 36
F = 8.0

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
    #fidw.write(csv_header)
    df.loc[:,par_col] = df.loc[:,name_col].apply(lambda x: "~  {0}  ~".format(x))
    df.to_csv(fidw, mode="a",index=False)


    # line = ""
    # for irow, row in df.iterrows():
    #     row[par_col] =  " ~   {0}    ~".format(row[name_col])
    #     line = ",".join(row.astype(str).values.tolist()) + "\n"
    #     fidw.write(line)

    fidw.close()

def csv_to_ins(csv_file, name_col, obs_col, ins_file):
    df = pd.read_csv(csv_file)
    part1 = "l1 "
    secondary_marker = " ~,~ "
    for icol, col in enumerate(df.columns):
        if col in [obs_col]:
            break
        part1 = part1 + secondary_marker

    obs_names = df[name_col]
    fidw = open(ins_file, 'w')
    fidw.write("pif ~\n")
    fidw.write("l1\n") # header


    #for irow, row in df.iterrows():
    ncol = len(df.columns)
    for irow, name in enumerate(df.loc[:, name_col].values):

        line = part1 + "   !{0}!    ".format(name)
        if ncol> (icol+1):
            line = line + secondary_marker+"\n"
        else:
            line = line + "\n"

        fidw.write(line)

    fidw.close()

def pst_setup_ES_new(N=50,dim=2,obs_time_frac=2,obs_dimen_frac=2):
    """setup for a smoother run"""
    # prepare template file
    template_ws = "template96_pst_ies"


    if template_ws in os.listdir(os.getcwd()):
        shutil.rmtree(os.path.join(os.getcwd(), template_ws))
    os.mkdir(template_ws)

    # generate reference input truth
    df_input = pd.DataFrame(columns=['name', 'value'])
    state_names = ['x{}_00.000'.format(x) for x in range(dim)]
    names = state_names + ['delt', 't_start', 't_end', 'nobs_t', 'nobs_loc', 'is_random', 'is_seq']
    # prior center / first guess (climatological mean F) and prior std - the parvals and the prior
    # ensemble are centered here so the truth is not baked into the prior
    prior_mean = np.ones(dim) * F
    prior_sigma = 1.0
    # truth: a draw from the prior restricted to within 3 sigma of the prior mean, so it is a
    # plausible and recoverable truth (not out in the tails the inversion can't reach)
    np.random.seed(111)
    initial_states = prior_mean + prior_sigma * np.clip(np.random.normal(0.0, 1.0, dim), -3.0, 3.0)
    nobs_t = 100
    nobs_loc = 8
    delt = 0.1
    t_start, t_end = 0.0, 20
    vals = list(initial_states) + [delt, t_start, t_end, nobs_t, nobs_loc, 0, 0]  # note is_random = 0 meaning that deterministic run; last val = is_seq, default to not seq
    df_input['name'] = names
    df_input['value'] = vals
    df_input.to_csv(os.path.join(template_ws, "lorenz96_in.csv"), index=False)
    # persist the truth initial states - lorenz96_in.csv gets overwritten by parvals during runs,
    # and parval1 is now the prior center (not the truth)
    pd.DataFrame({"name": state_names, "value": initial_states}).to_csv(
        os.path.join(template_ws, "truth_states.csv"), index=False)

    src_forward_run = inspect.getsource(lorenz96_model_setup.forward_run)
    src_rand_evolve = inspect.getsource(lorenz96_model_setup.rand_evolve)
    src_f = inspect.getsource(lorenz96_model_setup.f)

    script = """
import os, sys
import shutil
import numpy as np
import pandas as pd
from scipy.integrate import odeint

template_ws = {}
{}
{}
{}

if __name__ == "__main__":
     forward_run()

    """.format("'.'", src_f, src_rand_evolve, src_forward_run)

    with open(os.path.join(template_ws, r"forward_model.py"), 'w') as fidw:
        fidw.write(script)

    # run the model
    base_folder = os.getcwd()
    sys.path.append(template_ws)
    os.chdir(template_ws)
    cmd = sys.executable
    cmd = cmd + " " + "forward_model.py"
    os.system(cmd)
    os.chdir(base_folder)

    obs_df = pd.read_csv(os.path.join(template_ws, r"lorenz96_out.csv"))
    obs_df.to_csv(os.path.join(template_ws, r"lorenz96_out_truth.csv"))

    # generate pest
    final_states = []
    for obs_ in obs_df['obsname'].values:
        if "sim_" in obs_:
            obs_fin = obs_.replace("sim_", "final_")
            final_states.append(obs_fin)

    par_list = df_input['name'].values.tolist() + final_states
    pst = pyemu.helpers.pst_from_parnames_obsnames(parnames=par_list, obsnames=obs_df['obsname'].values)

    pst.parameter_data = pst.parameter_data.loc[par_list, :]
    pst.parameter_data.loc[state_names,"parval1"] = prior_mean
    pst.observation_data = pst.observation_data.loc[obs_df['obsname'].values, :]

    # generate tpl
    par_infile = os.path.join(template_ws, "lorenz96_in.csv")
    csv_to_tpl(csv_file=par_infile
               , par_col='value', name_col='name',
               tpl_file=par_infile + ".tpl")

    # generate ins
    out_file = os.path.join(template_ws, "lorenz96_out.csv")
    csv_to_ins(csv_file=out_file, name_col='obsname', obs_col='simval',
               ins_file=out_file + ".ins")

    # par
    pst.parameter_data['cycle'] = 0
    pst.parameter_data['parchglim'] = 'relative'
    pst.parameter_data['parlbnd'] = -9999
    pst.parameter_data['parubnd'] = 9999
    pst.parameter_data['partrans'] = "none"
    pst.parameter_data['pargp'] = "x"

    pst.parameter_data.loc[pst.parameter_data['parnme'] == 'delt', 'partrans'] = 'fixed'
    pst.parameter_data.loc[pst.parameter_data['parnme'] == 'delt', 'parval1'] = delt

    pst.parameter_data.loc[pst.parameter_data['parnme'] == 't_start', 'partrans'] = 'fixed'
    pst.parameter_data.loc[pst.parameter_data['parnme'] == 't_start', 'parval1'] = t_start

    pst.parameter_data.loc[pst.parameter_data['parnme'] == 't_end', 'partrans'] = 'fixed'
    pst.parameter_data.loc[pst.parameter_data['parnme'] == 't_end', 'parval1'] = t_end

    pst.parameter_data.loc[pst.parameter_data['parnme'] == 'nobs_t', 'partrans'] = 'fixed'
    pst.parameter_data.loc[pst.parameter_data['parnme'] == 'nobs_t', 'parval1'] = 100.0

    pst.parameter_data.loc[pst.parameter_data['parnme'] == 'nobs_loc', 'partrans'] = 'fixed'
    pst.parameter_data.loc[pst.parameter_data['parnme'] == 'nobs_loc', 'parval1'] = 9.0

    pst.parameter_data.loc[pst.parameter_data['parnme'] == 'is_random', 'partrans'] = 'fixed'
    pst.parameter_data.loc[pst.parameter_data['parnme'] == 'is_random', 'parval1'] = 0

    pst.parameter_data.loc[pst.parameter_data['parnme'] == 'is_seq', 'partrans'] = 'fixed'
    pst.parameter_data.loc[pst.parameter_data['parnme'] == 'is_seq', 'parval1'] = 0


    # link initial to final states
    for p in pst.parameter_data['parnme'].values:
        if "final_" in p:
            pp = p.replace("final_", "ini_")
            pst.parameter_data.loc[pst.parameter_data['parnme'] == p, 'state_par_link'] = pp

    # obs
    obs_df = obs_df.set_index(['obsname'])
    pst.observation_data['obsval'] = obs_df['simval']

    pst.observation_data['weight'] = 20  # this is 1/std
    pst.observation_data['cycle'] = 0

    # observe the first half of the simulated times at every location; the last half is left
    # unobserved (zero weight) to be treated as unseen "predictions" - independent of dimension
    t_obs_max = (t_end - t_start) / obs_time_frac

    def set_weight_from_obs(oo_nm):
        otime = float(oo_nm.split("_")[1])
        loc = int(oo_nm.split("_")[0][1:])
        if otime < t_obs_max:
            
            if dim < obs_dimen_frac:
                return 2.5
            elif loc % obs_dimen_frac == 0:
                return 2.5
            else: 
                return 0.0
        else:
            return 0.0

    pst.observation_data.loc[:,"weight"] = pst.observation_data.obsnme.apply(lambda x: set_weight_from_obs(x))
    pst.observation_data["standard_deviation"] = np.nan
    pst.observation_data.loc[pst.nnz_obs_names,"standard_deviation"] = 1.0/pst.observation_data.loc[pst.nnz_obs_names,"weight"]

    #pst.observation_data[['weight','state_par_link']] = pst.observation_data.apply(lambda x: dyn_name_from_obs(x['weight'], x['obsnme']), axis = 1 )

    # tpl
    fn_par = os.path.basename(par_infile)
    pst.model_input_data['pest_file'] = fn_par + ".tpl"
    pst.model_input_data['model_file'] = fn_par
    #pst.model_input_data.loc[0] = [os.path.basename(dummy_file) + ".tpl", os.path.basename(dummy_file)]
    pst.model_input_data['cycle'] = 0

    # ins
    pst.model_output_data['pest_file'] = [os.path.basename(out_file) + ".ins"]
    pst.model_output_data['model_file'] = [os.path.basename(out_file)]
    pst.model_output_data['cycle'] = 0

    # ens - centered on the prior mean (first guess), NOT the truth
    ens = np.array(prior_mean) + prior_sigma * np.random.randn(N, dim)
    #ens0 = 0.0 * np.random.randn(N, num_dynamic_states * dim * 2)  # TODO: check 0 # 3 variables (x,y,z) for each final and initial
    #ens = np.hstack([ens, ens0])
    ens_names = state_names# + dynamic_states + final_states
    ens = pd.DataFrame(ens, columns=ens_names)
    ens.to_csv(os.path.join(template_ws, "initial_x_ens.csv"))

    pst.pestpp_options['da_parameter_ensemble'] = "initial_x_ens.csv"
    #pst.pestpp_options['da_subset_size'] = N
    pst.pestpp_options['da_num_reals'] = N
    # pst.pestpp_options['ies_init_lam'] = [100]
    #pst.pestpp_options['ies_lambda_mults'] = 1
    #pst.pestpp_options['lambda_scale_fac'] = 1
    pst.control_data.noptmax = 1
    pst.pestpp_options['da_use_mda'] = False
    pst.pestpp_options["ies_include_base"] = False

    pst.pestpp_options['da_use_simulated_states'] = False

    #cmd = sys.executable
    cmd = "python forward_model.py"
    pst.model_command = cmd

    pst_file = 'es.pst'
    pst.control_data.noptmax = 0
    pst.write(os.path.join(template_ws, pst_file), version=2)
    pyemu.os_utils.run("pestpp-ies {0}".format(pst_file),cwd=template_ws)

    pass

def run_ies(t_d="template96_pst_ies",num_workers=5):
    pst = pyemu.Pst(os.path.join(t_d,"es.pst"))
    pst.control_data.noptmax = 5
    pst.write(os.path.join(t_d,"es.pst"))
    pyemu.os_utils.start_workers(t_d,"pestpp-ies","es.pst",num_workers=num_workers,master_dir="master_ies")


def mod_to_seq():
    bt_d = "template96_pst_ies"
    assert os.path.exists(bt_d)
    t_d = "template96_pst_da"
    if os.path.exists(t_d):
        shutil.rmtree(t_d)
    shutil.copytree(bt_d,t_d)
    bt_d = None

    # first run the model once for single time.
    df = pd.read_csv(os.path.join(t_d,"lorenz96_in.csv"),index_col=0)
    df.loc["t_start","value"] = 0.0
    df.loc["t_end", "value"] = df.loc["delt","value"]
    #df.loc["nobs_t", "value"] = 2
    #df.loc["is_seq", "value"] = 1
    df.to_csv(os.path.join(t_d,"lorenz96_in.csv"))

    pyemu.os_utils.run("python forward_model.py",cwd=t_d)
    pst = pyemu.Pst(os.path.join(t_d,"es.pst"))
    par = pst.parameter_data
    par.loc["t_end","parval1"] = par.loc["delt","parval1"]
    par.loc["t_start","parval1"] = 0.0


    # get the original obs data and parse cycle and location from obs name
    org_obs = pst.observation_data.copy()
    org_obs.loc[:,"otime"] = org_obs.obsnme.apply(lambda x: float(x.split("_")[1]))
    otimes = org_obs.otime.unique()
    otimes.sort()
    all_cycles = np.arange(otimes.shape[0],dtype=int)
    time_cycle_dict = {t:c for t,c in zip(otimes,all_cycles)}
    org_obs.loc[:,"cycle"] = org_obs.otime.apply(lambda x: time_cycle_dict[x])
    org_obs.loc[:, "location"] = org_obs.obsnme.apply(lambda x: int(x.split("_")[0][1:]))
    nnz_org_obs = org_obs.loc[org_obs.weight>0,:]
    nnz_locs = set(nnz_org_obs.location.unique().tolist())
    print(org_obs)

    # now reset the observation data to just a single cycle

    out_file = os.path.join(t_d, "lorenz96_out.csv")
    pst.drop_observations(out_file+".ins",pst_path=".")
    csv_to_ins(csv_file=out_file, name_col='obsname', obs_col='simval',
               ins_file=out_file + ".ins")
    pst.add_observations(out_file+".ins",pst_path=".")

    # set the weights in the control file
    obs = pst.observation_data
    obs.loc[:,"location"] = obs.obsnme.apply(lambda x: int(x.split("_")[0][1:]))
    obs.loc[:,"weight"] = 20.0
    #obs.loc[obs.location.apply(lambda x: x in nnz_locs),"weight"] = 20.0
    obs.loc[:,"cycle"] = -1

    # write the obs cycle table
    cycles = list(nnz_org_obs.cycle.unique())
    cycles.sort()
    #nnz_obs = pst.nnz_obs_names
    #nnz_locs = obs.loc[nnz_obs,"location"]
    location = obs.location.values
    #all_cycles = np.arange(0, max(cycles)+1, dtype=int)
    wdf = pd.DataFrame(index=obs.obsnme.values,columns=all_cycles)
    odf = pd.DataFrame(index=obs.obsnme.values, columns=all_cycles)
    for cycle in all_cycles:
        cobs = org_obs.loc[org_obs.cycle == cycle,:].copy()
        cobs.index = cobs.location
        odf.loc[:,cycle] = cobs.loc[location,"obsval"].values
        wdf.loc[:, cycle] = cobs.loc[location, "weight"].values
    odf.to_csv(os.path.join(t_d,"obs_cycle_table.csv"))
    wdf.to_csv(os.path.join(t_d, "weight_cycle_table.csv"))
    pst.pestpp_options["da_observation_cycle_table"] = "obs_cycle_table.csv"
    pst.pestpp_options["da_weight_cycle_table"] = "weight_cycle_table.csv"

    #add est state pars
    par = pst.parameter_data
    par.loc[pst.adj_par_names,"est_parnme"] = par.loc[pst.adj_par_names,"parnme"].apply(lambda x: "est_{0}".format(x))
    with open(os.path.join(t_d,"est_states.dat.tpl"),'w') as f:
        f.write("ptf ~\n")
        for name in par.loc[pst.adj_par_names,"est_parnme"]:
            f.write(" ~    {0}    ~\n".format(name))

    df = pst.add_parameters(os.path.join(t_d,"est_states.dat.tpl"),pst_path=".")
    par = pst.parameter_data
    par.loc[df.parnme,"parval1"] = 1.0
    par.loc[df.parnme, "partrans"] = "none"
    par.loc[df.parnme, "parlbnd"] = -9999
    par.loc[df.parnme, "parubnd"] = 9999
    par.loc[:,"state_par_link"] = np.nan
    par.loc[df.parnme,"state_par_link"] = df.parnme.apply(lambda x: x.replace("est_",""))
    pe_file = os.path.join(t_d, "initial_x_ens.csv")
    pe = pd.read_csv(pe_file)
    pe.loc[:, df.parnme] = pe.loc[:, par.loc[df.parnme, "state_par_link"].values].values
    assert pe.shape == pe.dropna().shape
    pe.to_csv(pe_file)

    # write a par cycle table for the fixed pars - this is only needed so that we can
    # drive pestpp-da thru all cycles, not just cycles with nnz obs
    # cycles = np.arange(0, max(cycles), dtype=int)
    # fpar = par.loc[par.partrans=="fixed",:]
    # df = pd.DataFrame(columns=cycles,index=fpar.parnme.values)
    # for fname,fval in zip(fpar.parnme.values,fpar.parval1.values):
    #     df.loc[fname,:] = fval
    # df.to_csv(os.path.join(t_d,"par_cycle_table.csv"))
    # pst.pestpp_options["da_parameter_cycle_table"] = "par_cycle_table.csv"


    pst.parameter_data.loc[:,"cycle"] = -1
    pst.model_input_data.loc[:,"cycle"] = -1
    pst.model_output_data.loc[:,"cycle"] = -1
    pst.control_data.noptmax = 0

    pst.write(os.path.join(t_d,"enkf.pst"),version=2)
    #pyemu.os_utils.run("pestpp-da enkf.pst",cwd=t_d)
    #nzobs = obs.loc[obs.weight>0,:].copy()
    return t_d


    # pdf = pd.DataFrame(index=["t_start","t_end"],columns=cycles[:-1])
    # for icycle in cycles[1:]:
    #     pdf.loc["t_start",icycle-1] = steps[icycle-1]
    #     pdf.loc["t_end", icycle-1] = steps[icycle]
    # print(pdf)
    # pdf.to_csv(os.path.join(t_d,"par_cycle_table.csv"))
    # pst.pestpp_options["da_par_cycle_table"] = "par_cycle_table.csv"

def run_da(t_d="template96_pst_da",num_workers=5,use_sim_states=True,noptmax=1,m_d=None,worker_root="lorenz"):
    pst = pyemu.Pst(os.path.join(t_d,"enkf.pst"))
    pst.control_data.noptmax = noptmax
    pst.pestpp_options["da_use_simulated_states"] = use_sim_states

    if m_d is None:
        if use_sim_states:
            m_d = os.path.join(worker_root,"master_enkf_simstates")
        else:
            m_d = os.path.join(worker_root,"master_enkf_eststates")
            par = pst.parameter_data

            istates_pars = par.loc[par.parnme.str.startswith("x"),"parnme"]
            par.loc[istates_pars,"partrans"] = "fixed"
    if noptmax == -1:
        m_d += "_openloop"
    pst.write(os.path.join(t_d,"enkf.pst"),version=2)
    pyemu.os_utils.start_workers(t_d,exe_path,"enkf.pst",num_workers=num_workers,master_dir=m_d,worker_root="lorenz")


def invest():
    import sys

    t_d = "template96_pst_da"
    sys.path.append(t_d)
    import forward_model
    # start with the initial values in the control file
    pst = pyemu.Pst(os.path.join(t_d,"enkf.pst"))
    delt = 0.01
    par = pst.parameter_data
    par.loc[par.parnme.str.startswith("x"), "parval1"] = 8.0
    par.loc[par.parnme.str.startswith("x"), "parlbnd"] = -20.0
    par.loc[par.parnme.str.startswith("x"), "parubnd"] = 20.0

    par.loc["x19_00.000", "parval1"] = 8.01
    pst.control_data.noptmax = -1
    pst.pestpp_options.pop("da_parameter_ensemble",None)
    pst.pestpp_options["da_num_reals"] = 5
    pst.write(os.path.join(t_d,"test.pst"),version=2)
    #pyemu.os_utils.run("pestpp-da test.pst",cwd=t_d)
    m_d = "enkf_test_master"
    #pyemu.os_utils.start_workers(t_d,"pestpp-da","test.pst",num_workers=5,worker_root=".",master_dir=m_d)
    da_df_dict = {}
    goe_files = [f for f in os.listdir(m_d) if "test.global" in f and f.endswith("oe.csv")]
    print(goe_files)
    for f in goe_files:
        cycle = int(f.split('.')[2])
        otime = (cycle * delt) + delt
        df = pd.read_csv(os.path.join(m_d,f),index_col=0)
        da_df_dict[otime] = df.T

    otimes = list(da_df_dict.keys())
    otimes.sort()

    df_in = pd.read_csv(os.path.join(t_d,"lorenz96_in.csv"),index_col=0)
    idx = df_in.index.map(lambda x: x.startswith("x"))
    df_in.loc[idx,"value"] = 8.0
    df_in.loc["x19_00.000", "value"] = 8.01

    df_in.loc["t_end","value"] = 20.0
    print(df_in)
    df_in.to_csv(os.path.join(t_d,"lorenz96_in.csv"))
    b_d = os.getcwd()
    os.chdir(t_d)
    forward_model.forward_run()
    os.chdir(b_d)
    out_df = pd.read_csv(os.path.join(t_d,"lorenz96_out.csv"))
    out_df.loc[:,"otime"] = out_df.obsname.apply(lambda x: float(x.split('_')[1]))
    x20 = out_df.loc[out_df.obsname.str.startswith("x19"),:].copy()
    x20.sort_values(by="otime")
    x20.index = x20.otime
    #x20.simval.plot()
    #plt.show()
    #return

    pvals = par.loc[par.parnme.str.startswith("x"),"parval1"].values
    pvals[:] = 8
    pvals[19] = 8.01
    #print(pvals)

    istates,fstates = [],[]
    for itime in range(2000):
        print(itime)
        new_vals = forward_model.rand_evolve(state0=pvals,t=[delt])[0]
        istates.append(pvals)
        fstates.append(new_vals)
        pvals = new_vals
        #print(new_vals)

    df = pd.DataFrame(fstates)
    print(df.shape,x20.shape)
    print(df.values)
    print(x20)
    #x20.index = df.index
    df.index = x20.index

    df_in.loc["t_end", "value"] = 0.01
    print(df_in)

    b_d = os.getcwd()
    os.chdir(t_d)
    out_dfs = {}
    for itime in range(2000):
        df_in.to_csv("lorenz96_in.csv")
        forward_model.forward_run()

        out_df_time = pd.read_csv("lorenz96_out.csv",index_col=0)
        otime = (itime * delt) + delt
        out_dfs[otime] = out_df_time
        df_in.loc[out_df_time.index,"value"] = out_df_time.simval.values
    os.chdir(b_d)

    fig,ax = plt.subplots(1,1)
    df.iloc[:,19].plot(ax=ax,color="r")
    x20.simval.plot(ax=ax,color='b',ls="--")
    for otime,df in out_dfs.items():
        ax.scatter([otime],[df.loc["x19_00.000","simval"]],marker="^",c="0.5",alpha=0.1)
    for otime,df in da_df_dict.items():
        ax.scatter([otime for _ in range(df.shape[1])],df.loc["x19_00.000",:].values,marker="+",c="m")
    plt.show()


def plot():
    bat_m_d = "template96_pst_ies"
    seq_m_d = "master_enkf_eststates_lessobs"

    pst = pyemu.Pst(os.path.join(bat_m_d,"es.pst"))
    obs = pst.observation_data
    obs.loc[:,"dim"] = obs.obsnme.apply(lambda x: int(x.split('_')[0][1:]))
    obs.loc[:, "otime"] = obs.obsnme.apply(lambda x: float(x.split('_')[1]))
    utimes = obs.otime.unique()
    utimes.sort()
    cycles = np.arange(len(utimes),dtype=int)
    time_cycle_dict = {t:c for t,c in zip(utimes,cycles)}
    cycle_time_dict = {c:t for t,c in zip(utimes,cycles)}
    obs.loc[:,"cycle"] = obs.otime.apply(lambda x: time_cycle_dict[x])

    global_oe_files = [f for f in os.listdir(seq_m_d) if "global" in f and "oe" in f]
    print(global_oe_files)
    time_globaloe_dict = {}
    for f in global_oe_files:
        df = pd.read_csv(os.path.join(seq_m_d,f),index_col=0)
        c = int(f.split('.')[2])
        time_globaloe_dict[cycle_time_dict[c]] = df
        break
    dim_dict = {int(c.split("_")[0][1:]):c for c in df.columns}
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(os.path.join(seq_m_d,"results_temp.pdf")) as pdf:
        dims = obs.dim.unique()
        dims.sort()
        for dim in dims:
            fig,ax = plt.subplots(1,1,figsize=(6,3))
            #for otime,df in time_globaloe_dict.items():
            #    ax.scatter([otime for _ in range(df.shape[0])],df.loc[:,dim_dict[dim]],marker='.',color='b')
            dobs = obs.loc[obs.dim==dim,:].copy()
            dobs.sort_values(by="otime",inplace=True)
            ax.plot(dobs.otime,dobs.obsval,"r-",lw=0.5)
            nzdobs = dobs.loc[dobs.weight>0,:].copy()
            ax.scatter(nzdobs.otime,nzdobs.obsval,marker='^',color="r")
            ax.set_title("dimension {0}".format(dim),loc="left")
            pdf.savefig()
            plt.close(fig)


def plot_pr_pt(use_range=True):
    bat_m_d = "template96_pst_ies"
    pr_seq_m_d = "master_enkf_eststates_openloop"
    pt_seq_m_d = "master_enkf_eststates"

    pst = pyemu.Pst(os.path.join(bat_m_d,"es.pst"))
    obs = pst.observation_data
    obs.loc[:,"dim"] = obs.obsnme.apply(lambda x: int(x.split('_')[0][1:]))
    obs.loc[:, "otime"] = obs.obsnme.apply(lambda x: float(x.split('_')[1]))
    utimes = obs.otime.unique()
    utimes.sort()
    cycles = np.arange(len(utimes),dtype=int)
    time_cycle_dict = {t:c for t,c in zip(utimes,cycles)}
    cycle_time_dict = {c:t for t,c in zip(utimes,cycles)}
    obs.loc[:,"cycle"] = obs.otime.apply(lambda x: time_cycle_dict[x])

    pt_global_oe_files = [f for f in os.listdir(pt_seq_m_d) if "global" in f and "oe" in f]
    pt_time_globaloe_dict = {}
    for f in pt_global_oe_files:
        df = pd.read_csv(os.path.join(pt_seq_m_d, f), index_col=0)
        c = int(f.split('.')[2])
        pt_time_globaloe_dict[cycle_time_dict[c]] = df

    pr_global_oe_files = [f for f in os.listdir(pr_seq_m_d) if "global" in f and "oe" in f]
    pr_time_globaloe_dict = {}
    for f in pr_global_oe_files:
        df = pd.read_csv(os.path.join(pr_seq_m_d,f),index_col=0)
        c = int(f.split('.')[2])
        pr_time_globaloe_dict[cycle_time_dict[c]] = df
    dim_dict = {int(c.split("_")[0][1:]):c for c in df.columns}


    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(os.path.join(pt_seq_m_d,"results.pdf")) as pdf:
        dims = obs.dim.unique()
        dims.sort()
        for dim in dims:
            fig,ax = plt.subplots(1,1,figsize=(10,5))
            for otime,df in pr_time_globaloe_dict.items():
                ax.scatter([otime for _ in range(df.shape[0])],df.loc[:,dim_dict[dim]],marker='.',color='0.5',alpha=0.1)
            for otime, df in pt_time_globaloe_dict.items():
                ax.scatter([otime for _ in range(df.shape[0])], df.loc[:, dim_dict[dim]], marker='.', color='b',
                           alpha=0.1)

            dobs = obs.loc[obs.dim==dim,:].copy()
            dobs.sort_values(by="otime",inplace=True)
            ax.plot(dobs.otime,dobs.obsval,"r-",lw=0.5)
            nzdobs = dobs.loc[dobs.weight>0,:].copy()
            ax.scatter(nzdobs.otime,nzdobs.obsval,marker='^',color="r")
            ax.set_title("dimension {0}".format(dim),loc="left")
            #ax.set_xlim(0,5)
            pdf.savefig()
            plt.close(fig)


def lorenz96_basic_test():
    pst_setup_ES_new(N=50)
    t_d = mod_to_seq()
    # invest()
    m_d = "master_da_test"
    if os.path.exists(m_d):
        shutil.rmtree(m_d)

    pst = pyemu.Pst(os.path.join(t_d, "enkf.pst"))
    pst.control_data.noptmax = 1
    pst.pestpp_options["da_use_simulated_states"] = False
    pst.pestpp_options["da_num_reals"] = 10
    par = pst.parameter_data
    istates_pars = par.loc[par.parnme.str.startswith("x"), "parnme"]
    par.loc[istates_pars, "partrans"] = "fixed"

    pst.write(os.path.join(t_d, "enkf.pst"), version=2)
    pyemu.os_utils.start_workers(t_d, exe_path, "enkf.pst", num_workers=10, master_dir=m_d,
                                 worker_root=".",port=port)
    # find how many cycles we are using
    num_cycles = None
    with open(os.path.join(m_d,"enkf.rec")) as f:
        for line in f:
            if line.startswith("...assimilating over"):
                raw = line.strip().split()
                num_cycles = int(raw[2])
    assert num_cycles is not None
    # find the global pe files
    gpe_files = [f for f in os.listdir(m_d) if "global" in f and "pe" in f]
    assert len(gpe_files) == num_cycles + 1,"{0} vs {1}".format(num_cycles,len(gpe_files))

def lorenz96_ies_dim_test(dim_use=36, num_reals=50, noptmax=3, num_workers=10):
    """dynamically-sized lorenz96 pestpp-ies smoother test - set dim_use to play with the number
    of states (the number of observations scales with the dimension)"""
    pst_setup_ES_new(N=num_reals, dim=dim_use)
    t_d = "template96_pst_ies"
    pst = pyemu.Pst(os.path.join(t_d, "es.pst"))
    pst.control_data.noptmax = noptmax
    pst.write(os.path.join(t_d, "es.pst"), version=2)
    m_d = "master_ies_dim{0}".format(dim_use)
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    if num_workers == 0:
        shutil.copytree(t_d, m_d)
        pyemu.os_utils.run("pestpp-ies es.pst", cwd=m_d)
    else:
        pyemu.os_utils.start_workers(t_d, "pestpp-ies", "es.pst", num_workers=num_workers,
                                     master_dir=m_d, worker_root=".", port=port)
    phidf = pd.read_csv(os.path.join(m_d, "es.phi.actual.csv"))
    print("dim", dim_use, "npar_adj", pst.npar_adj, "nnz_obs", pst.nnz_obs,
          "phi mean start/end:", phidf["mean"].iloc[0], phidf["mean"].iloc[-1])
    assert phidf.iteration.max() == noptmax
    assert phidf["mean"].iloc[-1] < phidf["mean"].iloc[0]
    return m_d

def _ies_result_handler(m_d, case):
    """return (ies_handler, sorted_iters) from the pyemu result handler for a pestpp-ies master dir.
    iterations are discovered from the per-iteration obs ensemble files the handler already tracks,
    so the plotting routines never have to open csv files themselves."""
    ies = pyemu.pst.result_handler.Results(m_d, case=case).ies
    iters = sorted({int(os.path.basename(f).split(".")[1]) for f in ies.get_files(".obs.")})
    assert len(iters) > 0, "no obs ensembles found in " + m_d
    return ies, iters


def lorenz96_dim3_dashboard(m_d="master_ies_dim3", case="es", which="zero-weighted"):
    """interactive 3D dashboard for the lorenz test: the truth trajectory and the ies estimated
    ensemble of trajectories in the phase space of 3 state dimensions - with a slider to pick the
    ies iteration, one for the max time, and a radio toggle to choose whether the 3 plotted
    dimensions are zero-weighted (unobserved, pure predictions) or weighted (observed) states.  the
    truth's analysis-window portion is solid black, its forecast portion red dashed; an actual-phi
    histogram panel updates with the iteration slider.  `which` sets the initial toggle selection."""
    from matplotlib.widgets import Slider, RadioButtons

    # pyemu result handler for this master dir + the available ies iterations
    ies, iters = _ies_result_handler(m_d, case)

    # truth (obsval) and obs layout from the pst
    pst = pyemu.Pst(os.path.join(m_d, case + ".pst"))
    od = pst.observation_data.copy()
    od["state"] = od.obsnme.apply(lambda x: x.split("_")[0])
    od["time"] = od.obsnme.apply(lambda x: float(x.split("_")[1]))
    # split the states into zero-weight (no observations anywhere - pure predictions) and weighted
    # (observed at some times); the time split into analysis/forecast windows is global (same for both)
    all_states = sorted(od.state.unique(), key=lambda s: int(s[1:]))
    obs_state_set = set(od.loc[od.weight > 0, "state"].unique())
    unobs_states = [s for s in all_states if s not in obs_state_set]
    wt_states = [s for s in all_states if s in obs_state_set]
    times = sorted(od.time.unique())
    times_arr = np.array(times)
    t_bnd = od.loc[od.weight > 0, "time"].max()  # observed/prediction boundary
    n_obs = sum(1 for t in times if t <= t_bnd)   # number of observed time steps

    # build the per-mode data (first 3 states, truth, preloaded per-iteration trajectories, and the
    # locked axis extents so the box does not rescale/"jitter" within a mode) for each mode that has
    # at least 3 state dimensions available
    def build_mode(state_pool):
        states = state_pool[:3]
        cols = {s: ["{0}_{1:06.3f}".format(s, t) for t in times] for s in states}
        truth = np.column_stack([od.loc[od.state == s].set_index("time").loc[times, "obsval"].values for s in states])
        traj = {}
        for it in iters:
            df = ies.get("obsen", it)
            df.columns = df.columns.str.lower()
            traj[it] = np.stack([df.loc[:, cols[s]].values for s in states], axis=-1)
        allpts = np.concatenate([traj[it].reshape(-1, 3) for it in iters] + [truth], axis=0)
        lo, hi = np.nanmin(allpts, axis=0), np.nanmax(allpts, axis=0)
        pad = 0.05 * (hi - lo)
        return {"states": states, "truth": truth, "traj": traj, "lo": lo - pad, "hi": hi + pad}

    modes = {}
    if len(unobs_states) >= 3:
        modes["zero-weighted"] = build_mode(unobs_states)
    if len(wt_states) >= 3:
        modes["weighted"] = build_mode(wt_states)
    assert len(modes) > 0, "need at least 3 weighted or 3 zero-weighted state dims; got {0} weighted, {1} zero-weighted".format(len(wt_states), len(unobs_states))
    labels = [l for l in ["zero-weighted", "weighted"] if l in modes]
    active = which if which in modes else labels[0]

    # per-realization actual phi by iteration (for the histogram panel, if available)
    phidf = ies.phiactual
    have_phi = phidf is not None
    if have_phi:
        real_cols = list(phidf.columns[6:])  # after iteration,total_runs,mean,std,min,max
        allphi = phidf[real_cols].values.astype(float).ravel()
        allphi = allphi[np.isfinite(allphi) & (allphi > 0)]
        pbins = np.logspace(np.log10(allphi.min()), np.log10(allphi.max()), 30)
        pcount_max = 1
        for it in iters:
            v = phidf.loc[phidf.iteration == it, real_cols].values.astype(float).ravel()
            v = v[np.isfinite(v) & (v > 0)]
            if v.size:
                pcount_max = max(pcount_max, int(np.histogram(v, bins=pbins)[0].max()))

    if have_phi:
        fig = plt.figure(figsize=(10.5, 7.4))
        ax = fig.add_axes([0.02, 0.16, 0.62, 0.80], projection="3d")
        axh = fig.add_axes([0.73, 0.42, 0.24, 0.42])
        rax = fig.add_axes([0.76, 0.12, 0.20, 0.14]) if len(modes) > 1 else None
    else:
        fig = plt.figure(figsize=(8, 7.4))
        ax = fig.add_subplot(111, projection="3d")
        plt.subplots_adjust(bottom=0.16)
        rax = fig.add_axes([0.015, 0.02, 0.15, 0.12]) if len(modes) > 1 else None

    def draw(it, tmax):
        m = modes[active]
        states, truth, traj, lo, hi = m["states"], m["truth"], m["traj"], m["lo"], m["hi"]
        obs_lab = "truth (observed)" if active == "weighted" else "truth (analysis window)"
        elev, azim = ax.elev, ax.azim  # preserve the view across redraws
        ax.clear()
        nshow = max(int(np.searchsorted(times_arr, tmax, side="right")), 2)  # time steps to show
        arr = traj[it]
        for r in range(arr.shape[0]):
            ax.plot(arr[r, :nshow, 0], arr[r, :nshow, 1], arr[r, :nshow, 2], "b-", lw=0.5, alpha=0.25, zorder=2)
        n_obs_show = min(n_obs, nshow)
        ax.plot(truth[:n_obs_show, 0], truth[:n_obs_show, 1], truth[:n_obs_show, 2], "k-", lw=2.5, zorder=3, label=obs_lab)
        if nshow > n_obs:
            ax.plot(truth[n_obs - 1:nshow, 0], truth[n_obs - 1:nshow, 1], truth[n_obs - 1:nshow, 2], "r--", lw=2.5, zorder=3, label="truth (prediction)")
        ax.scatter(truth[0, 0], truth[0, 1], truth[0, 2], color="k", s=40, zorder=4)  # start
        ax.set_xlabel(states[0])
        ax.set_ylabel(states[1])
        ax.set_zlabel(states[2])
        # lock the axes to the active mode's global extent so the box stays fixed across redraws
        ax.set_xlim(lo[0], hi[0])
        ax.set_ylim(lo[1], hi[1])
        ax.set_zlim(lo[2], hi[2])
        ax.set_title("lorenz {0}-{1}-{2} ({3}) phase space, iter {4}, t<={5:.2f}".format(states[0], states[1], states[2], active, it, times_arr[nshow - 1]))
        ax.legend(loc="upper left", fontsize=8)
        ax.view_init(elev=elev, azim=azim)
        fig.canvas.draw_idle()

    def draw_hist(it):
        axh.clear()
        vals = phidf.loc[phidf.iteration == it, real_cols].values.astype(float).ravel()
        vals = vals[np.isfinite(vals) & (vals > 0)]
        if vals.size == 0:
            return
        axh.hist(vals, bins=pbins, color="0.6", edgecolor="k")
        axh.axvline(np.median(vals), color="b", ls="--", lw=1.2)
        axh.set_xscale("log")
        axh.set_xlim(pbins[0], pbins[-1])
        axh.set_ylim(0, pcount_max * 1.05)
        axh.set_xlabel("actual phi (log)", fontsize=8)
        axh.set_ylabel("count", fontsize=8)
        axh.set_title("phi, iter {0} (median {1:.3g})".format(it, np.median(vals)), fontsize=8)
        axh.tick_params(labelsize=7)

    ax_it = plt.axes([0.2, 0.07, 0.5, 0.025])
    ax_tm = plt.axes([0.2, 0.02, 0.5, 0.025])
    s_iter = Slider(ax_it, "iteration", min(iters), max(iters), valinit=iters[0], valstep=1)
    s_tmax = Slider(ax_tm, "max time", float(times_arr[1]), float(times_arr[-1]), valinit=float(times_arr[-1]))

    def redraw(_=None):
        it = min(iters, key=lambda x: abs(x - s_iter.val))
        draw(it, s_tmax.val)
        if have_phi:
            draw_hist(it)
    s_iter.on_changed(redraw)
    s_tmax.on_changed(redraw)

    # radio toggle between zero-weighted (unobserved) and weighted (observed) 3-d state subsets
    radio = None
    if rax is not None:
        radio = RadioButtons(rax, labels, active=labels.index(active))
        rax.set_title("dimensions", fontsize=8)

        def on_mode(label):
            nonlocal active
            active = label
            redraw()
        radio.on_clicked(on_mode)

    redraw()
    fig._sliders = (s_iter, s_tmax, radio)  # keep references so the widgets stay responsive
    plt.show()
    return fig

def lorenz96_crps_heatmap(m_d="master_ies_dim3", case="es", predictive_only=True,
                          noise_floor=None, thresh_factor=1.2):
    """heatmap of ensemble CRPS (continuous ranked probability score): rows = predictive (lead)
    times, columns = ies iterations.  CRPS is averaged over all state locations at each time; lower
    is a better (sharper + better-calibrated) ensemble forecast of the truth.  predictive_only keeps
    only the zero-weight (unobserved/prediction) times.  cells whose CRPS is <= thresh_factor *
    noise_floor (i.e. forecasts at/near the observation-noise level) are outlined with a red box.
    noise_floor defaults to the observation-noise std (mean 1/weight over the weighted obs)."""
    # pyemu result handler for this master dir + the available ies iterations
    ies, iters = _ies_result_handler(m_d, case)

    pst = pyemu.Pst(os.path.join(m_d, case + ".pst"))
    od = pst.observation_data.copy()
    od["state"] = od.obsnme.apply(lambda x: x.split("_")[0])
    od["time"] = od.obsnme.apply(lambda x: float(x.split("_")[1]))
    states = sorted(od.state.unique(), key=lambda s: int(s[1:]))
    times = sorted(od.time.unique())
    t_bnd = od.loc[od.weight > 0, "time"].max()  # last observed time
    pred_times = [t for t in times if t > t_bnd] if predictive_only else list(times)
    assert len(pred_times) > 0, "no predictive (zero-weight) times found"

    cols = {s: {t: "{0}_{1:06.3f}".format(s, t) for t in times} for s in states}
    od_i = od.set_index("obsnme")
    truth = np.array([[od_i.loc[cols[s][t], "obsval"] for s in states] for t in pred_times])  # (npt, nstate)

    def crps_cols(ens, tr):
        # empirical CRPS per column: E|X-y| - 0.5 E|X-X'|.  ens (m,K), tr (K,) -> (K,)
        m = ens.shape[0]
        term1 = np.abs(ens - tr[None, :]).mean(axis=0)
        s = np.sort(ens, axis=0)
        i = np.arange(m)[:, None]
        term2 = (2.0 / (m * m)) * ((2 * i - m + 1) * s).sum(axis=0)  # E|X-X'| via sorted form
        return term1 - 0.5 * term2

    crps = np.full((len(pred_times), len(iters)), np.nan)
    escore = np.full((len(pred_times), len(iters)), np.nan)
    for ci, it in enumerate(iters):
        df = ies.get("obsen", it)
        df.columns = df.columns.str.lower()
        for ri, t in enumerate(pred_times):
            ens = df.loc[:, [cols[s][t] for s in states]].values  # (m, nstate)
            crps[ri, ci] = crps_cols(ens, truth[ri]).mean()  # mean of per-state (marginal) CRPS
            # multivariate energy score: E||X-y|| - 0.5 E||X-X'|| (euclidean norm over the full state vector)
            d1 = np.linalg.norm(ens - truth[ri][None, :], axis=1)
            d2 = np.linalg.norm(ens[:, None, :] - ens[None, :, :], axis=2)
            escore[ri, ci] = d1.mean() - 0.5 * d2.mean()

    # observation-noise floor (CRPS is in obs units) and the "good forecast" threshold
    if noise_floor is None:
        noise_floor = float((1.0 / od.loc[od.weight > 0, "weight"]).mean())
    thresh = thresh_factor * noise_floor

    # per-realization actual phi by iteration (for the phi-vs-iteration panel below the heatmap)
    phidf = ies.phiactual
    real_cols = list(phidf.columns[6:])  # cols after iteration,total_runs,mean,std,min,max
    M = phidf.set_index("iteration")
    phi = np.full((len(real_cols), len(iters)), np.nan)  # (nreal, niter)
    for k, it in enumerate(iters):
        if it in M.index:
            phi[:, k] = pd.to_numeric(M.loc[it, real_cols], errors="coerce").values

    fig = plt.figure(figsize=(max(4.5, 1.6 + 0.6 * len(iters)), min(2 + 0.06 * len(pred_times), 9) * 2 + 2.6))
    # 3 rows (CRPS heatmap, energy-score heatmap, phi-vs-iteration) sharing the iteration x-axis;
    # colorbars in thin right-hand cells so the panels keep equal width and their x-axes line up
    gs = fig.add_gridspec(3, 2, width_ratios=[30, 1], height_ratios=[3, 3, 1.3], hspace=0.12, wspace=0.04,
                          left=0.16, right=0.92, top=0.93, bottom=0.08)
    ax = fig.add_subplot(gs[0, 0])                          # CRPS heatmap
    axe = fig.add_subplot(gs[1, 0], sharex=ax, sharey=ax)   # energy-score heatmap
    axp = fig.add_subplot(gs[2, 0], sharex=ax)              # phi vs iteration
    cax = fig.add_subplot(gs[0, 1])
    cae = fig.add_subplot(gs[1, 1])

    yt = list(range(0, len(pred_times), max(1, len(pred_times) // 20)))
    ylab = ["{0:.2f}".format(pred_times[i]) for i in yt]

    # CRPS heatmap (mean over states = mean MARGINAL CRPS), red box where CRPS <= noise-floor thresh
    im = ax.imshow(crps, aspect="auto", origin="lower", cmap="viridis")
    for ri in range(crps.shape[0]):
        for ci in range(crps.shape[1]):
            if crps[ri, ci] <= thresh:
                ax.add_patch(plt.Rectangle((ci - 0.5, ri - 0.5), 1, 1, fill=False, edgecolor="red", lw=1.5))
    ax.set_yticks(yt)
    ax.set_yticklabels(ylab)
    ax.set_ylabel("predictive time")
    ax.tick_params(labelbottom=False)
    ax.set_title("mean CRPS over states (marginal) - {0}  (red box: <= {1:.2f}x noise floor {2:.3g})".format(
        os.path.basename(m_d), thresh_factor, noise_floor), fontsize=9)
    fig.colorbar(im, cax=cax, label="CRPS")

    # energy-score heatmap (multivariate / joint over all state dimensions)
    ime = axe.imshow(escore, aspect="auto", origin="lower", cmap="magma")
    axe.set_yticks(yt)
    axe.set_yticklabels(ylab)
    axe.set_ylabel("predictive time")
    axe.tick_params(labelbottom=False)
    axe.set_title("energy score (multivariate, joint over states)", fontsize=9)
    fig.colorbar(ime, cax=cae, label="energy score")

    # phi vs iteration: one grey line per realization, blue median; x lines up with the heatmap columns
    xs = list(range(len(iters)))
    for r in range(phi.shape[0]):
        axp.plot(xs, phi[r], "-", color="0.5", lw=0.6, alpha=0.35)
    axp.plot(xs, np.nanmedian(phi, axis=0), "b-", lw=2.0, label="median")
    axp.set_yscale("log")
    axp.set_xticks(range(len(iters)))
    axp.set_xticklabels(iters)
    axp.set_xlim(-0.5, len(iters) - 0.5)
    axp.set_xlabel("ies iteration")
    axp.set_ylabel("actual phi (log)")
    axp.grid(True, which="both", axis="y", lw=0.3, alpha=0.4)
    axp.legend(fontsize=8, loc="upper right")

    out = os.path.join(m_d, "crps_heatmap.pdf")
    fig.savefig(out)
    print("saved", out, " CRPS {0:.3g}-{1:.3g}, energy {2:.3g}-{3:.3g}; CRPS thresh {4:.3g} ({5} cells boxed)".format(
        np.nanmin(crps), np.nanmax(crps), np.nanmin(escore), np.nanmax(escore), thresh, int(np.nansum(crps <= thresh))))
    return fig

def lorenz96_add_correlated_obs_noise(t_d="template96_pst_ies", case="es", corr_range=1.0,
                                      seed=11, plot_locs=4, plot_nreals=25):
    """generate temporally-autocorrelated observation noise with pyemu, build the noisy observation
    ensemble, attach it to the pst (ies_observation_ensemble) and plot the noise realizations.  the
    noise covariance is an exponential variogram in time (range = corr_range, time units) per
    location (locations independent), built with pyemu.geostats and drawn with
    pyemu.ObservationEnsemble.from_gaussian_draw.  experiment with corr_range: ~0 -> white noise,
    large -> smooth/slowly-varying."""
    pst = pyemu.Pst(os.path.join(t_d, case + ".pst"))
    od = pst.observation_data.copy()
    od["location"] = od.obsnme.apply(lambda x: int(x.split("_")[0][1:]))
    od["time"] = od.obsnme.apply(lambda x: float(x.split("_")[1]))
    nz = od.loc[od.weight > 0].copy()
    assert nz.shape[0] > 0, "no nonzero-weight obs to add noise to"

    pe = pd.read_csv(os.path.join(t_d, "initial_x_ens.csv"), index_col=0)
    num_reals = pe.shape[0]

    # build the correlated obs-noise covariance with pyemu geostats (exponential in time per
    # location); locations are made independent via a large offset in the 2nd coordinate
    if corr_range <= 0:
        cov = None  # white noise from the weights
    else:
        var = float(np.mean((1.0 / nz.weight.values) ** 2))  # obs variance (std = 1/weight)
        gs = pyemu.geostats.GeoStruct(variograms=[pyemu.geostats.ExpVario(contribution=var, a=corr_range)])
        cov = gs.covariance_matrix(x=nz["time"].values, y=nz["location"].values * 1.0e6,
                                   names=list(nz.obsnme.values))

    oe = pyemu.ObservationEnsemble.from_gaussian_draw(pst, cov=cov, num_reals=num_reals,
                                                      by_groups=False, rng=np.random.default_rng(seed))
    oe_file = "obs_noise_corr{0}.csv".format(str(corr_range).replace(".", "p"))
    oe.to_csv(os.path.join(t_d, oe_file))
    pst.pestpp_options["ies_observation_ensemble"] = oe_file
    pst.pestpp_options["ies_num_reals"] = num_reals
    pst.write(os.path.join(t_d, case + ".pst"), version=2)
    print("saved", oe_file, "-> ies_observation_ensemble; num_reals", num_reals)

    # visualize the noise (noisy obs - obsval) per location over the observed times
    oedf = pd.read_csv(os.path.join(t_d, oe_file), index_col=0)
    oedf.columns = oedf.columns.str.lower()
    obsval = od.set_index("obsnme")["obsval"]
    locs = sorted(nz["location"].unique())[:plot_locs]
    fig, axes = plt.subplots(len(locs), 1, figsize=(9, 1.8 * len(locs)), sharex=True, squeeze=False)
    for ax, loc in zip(axes[:, 0], locs):
        g = nz.loc[nz["location"] == loc].sort_values("time")
        names = g.obsnme.values
        t = g["time"].values
        nse = oedf.loc[:, names].values# - obsval.loc[names].values[None, :]
        for r in range(min(num_reals, plot_nreals)):
            ax.plot(t, nse[r], lw=0.7, alpha=0.5)
        ax.axhline(0.0, color="k", lw=0.6)
        ax.set_ylabel("x{0}".format(loc))
    axes[-1, 0].set_xlabel("observation time")
    fig.suptitle("autocorrelated obs noise realizations (corr_range = {0})".format(corr_range))
    fig.tight_layout()
    out = os.path.join(t_d, "obs_noise_corr{0}.pdf".format(str(corr_range).replace(".", "p")))
    fig.savefig(out)
    print("saved", out)
    return fig

def lorenz96_ext_runmanager_test(dim_use=40, num_reals=50, noptmax=5,obs_time_frac=10,obs_dimen_frac=10):
    """clean twin experiment run with the pest++ external ('/e') run manager: the truth observations
    are (re)generated with the SAME vectorized fixed-step model used for inversion, so the model is
    self-consistent (no model error).  the prior is centered on the first guess, NOT the truth, so
    this is a proper twin.  pestpp writes the full ensemble of parameter sets to the run-storage file
    (es.rns), the model command is called once to solve the whole ensemble in one pass and pack the
    results back (lorenz96_model_setup.forward_run_ensemble), then pestpp reads them back."""
    import re
    pst_setup_ES_new(N=num_reals, dim=dim_use,obs_dimen_frac=obs_dimen_frac,obs_time_frac=obs_time_frac)
    t_d = "template96_pst_ies"
    # write a self-contained ensemble forward-run script - the single /e model command
    src = inspect.getsource(lorenz96_model_setup.evolve_lorenz) + "\n\n" + \
          inspect.getsource(lorenz96_model_setup.forward_run_ensemble)
    with open(os.path.join(t_d, "forward_ens.py"), "w") as f:
        f.write("import os, sys, glob, re\nimport numpy as np\nimport pyemu\n\n")
        f.write(src)
        f.write("\nif __name__ == '__main__':\n    forward_run_ensemble()\n")

    pst = pyemu.Pst(os.path.join(t_d, "es.pst"))
    # --- clean twin: regenerate the truth observations with the SAME (fixed-step) model ---
    par = pst.parameter_data
    state_cols = sorted([p for p in pst.par_names if re.match(r"x\d+_0*0\.000$", p)],
                        key=lambda p: int(p.split("_")[0][1:]))
    # the truth ICs are persisted in truth_states.csv (parval1 is now the prior center, not the truth)
    truth_in = pd.read_csv(os.path.join(t_d, "truth_states.csv")).set_index("name")["value"]
    X0_truth = truth_in.loc[state_cols].values.astype(float)
    delt = float(par.loc["delt", "parval1"]); t_start = float(par.loc["t_start", "parval1"]); t_end = float(par.loc["t_end", "parval1"])
    times = np.arange(0.0, t_end - t_start, delt)
    traj = lorenz96_model_setup.evolve_lorenz(X0_truth[None, :], len(times))  # (n_out, 1, N)
    snames = [c.split("_")[0] for c in state_cols]
    truthmap = {"{0}_{1:06.3f}".format(sn, otime): traj[k, 0, i]
                for k, otime in enumerate(times) for i, sn in enumerate(snames)}
    pst.observation_data["obsval"] = pst.observation_data.obsnme.map(truthmap)
    assert pst.observation_data.obsval.notna().all(), "some obs not covered by the regenerated truth"
    # --- end twin truth regeneration ---

    pst.model_command = ["python forward_ens.py"]
    pst.control_data.noptmax = noptmax
    pst.pestpp_options.pop("ies_lambda_mults",None)
    pst.pestpp_options.pop("lambda_scale_fac",None)
    #pst.pestpp_options["ies_n_iter_reinflate"] = [-5,-5,999]
    pst.pestpp_options["save_binary"] = True
    pst.pestpp_options["save_dense"] = True
    pst.pestpp_options["ies_bad_phi_sigma"] = 1.5
    pst.pestpp_options["ies_multimodal_alpha"] = 1.0
    #pst.pestpp_options["ies_multimodal_phi_weight"] = 0.5
    #pst.pestpp_options["ies_multimodal_weight_exponent"] = 4
    pst.pestpp_options["ies_num_threads"] = 10
    
    pst.write(os.path.join(t_d, "es.pst"), version=2)
    lorenz96_add_correlated_obs_noise(t_d)
    m_d = "master_ies_ext_dim{0}".format(dim_use)
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    shutil.copytree(t_d, m_d)
    pyemu.os_utils.run("{0} es.pst /e".format(exe_path.replace("-da","-ies")), cwd=m_d)
    phidf = pd.read_csv(os.path.join(m_d, "es.phi.actual.csv"))
    print("ext run mgr (twin): dim", dim_use, "phi mean start/end:", phidf["mean"].iloc[0], phidf["mean"].iloc[-1])
   #assert phidf.iteration.max() == noptmax
    assert phidf["mean"].iloc[-1] < phidf["mean"].iloc[0]
    return m_d

if __name__ == "__main__":
    #lorenz96_basic_test()
    #m_d = lorenz96_ies_dim_test(dim_use=400,num_workers=10)
    m_d = lorenz96_ext_runmanager_test(dim_use=400,num_reals=500,noptmax=10,obs_dimen_frac=5,obs_time_frac=5)
    #lorenz96_add_correlated_obs_noise()
    lorenz96_crps_heatmap(m_d)
    lorenz96_dim3_dashboard(m_d=m_d)
    #pst_setup_ES_new(N=50)
    #forward_run()
    #pst_setup()
    #run_ies(num_workers=50)
    #mod_to_seq()
    #invest()

    #run_da(num_workers=50, use_sim_states=False,noptmax=-1)
    #run_da(num_workers=50, use_sim_states=False, noptmax=1)

    #run_da(num_workers=50, use_sim_states=True,noptmax=-1)

    #run_da(num_workers=50,use_sim_states=False)
    #run_da(num_workers=50,use_sim_states=True)

    #plot()
    #plot_pr_pt()


