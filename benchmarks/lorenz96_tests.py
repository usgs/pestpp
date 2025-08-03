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

def pst_setup_ES_new(N=50):
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
    initial_states = np.ones(dim) * F
    np.random.seed(111)
    initial_states += np.random.normal(0.0,3.0,dim)
    nobs_t = 100
    nobs_loc = 8
    delt = 0.01
    t_start, t_end = 0.0, 0.25
    vals = list(initial_states) + [delt, t_start, t_end, nobs_t, nobs_loc, 0, 0]  # note is_random = 0 meaning that deterministic run; last val = is_seq, default to not seq
    df_input['name'] = names
    df_input['value'] = vals
    df_input.to_csv(os.path.join(template_ws, "lorenz96_in.csv"), index=False)

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
    pst.parameter_data.loc[state_names,"parval1"] = initial_states
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

    # reduce to only 9 locs
    #obs_data['loc'] = obs_data.obsnme.apply(lambda x: x.split("_").replace("x", "") if obs_data.obgnme=="obs" else "none")
    #l = [x for x in range(4, 36 + 1)[::4]]
    #for ll in l:
        #obs_data.loc[(obs_data.obgnme == "obs") & (obs_data['loc'] == str(ll)), "weight"] = 20.0
    #dt_m = nobs_t / ((t_end - t_start) / delt)
    #ind_m = set(list((np.linspace(int(dt_m / delt), int(t_end / delt) - 1, int(nobs_t))).astype(int)))
    ind_loc = set(list(np.linspace(0, len(state_names) - 1, int(nobs_loc)).astype(int)))

    def set_weight_from_obs(oo_nm):
        oloc = int(oo_nm.split("_")[0][1:])
        otime = float(oo_nm.split("_")[1])
        #if "sim_" in oo_nm:
        if oloc in ind_loc:
            #if np.isclose(otime,t_end,rtol=delt/20,atol=delt/20) or np.isclose(otime % (delt * 20),0.0,rtol=delt,atol=delt):
            #    print(otime,otime % (delt * 20))
            return 2.5
        else:
            return 0.0

    pst.observation_data.loc[:,"weight"] = pst.observation_data.obsnme.apply(lambda x: set_weight_from_obs(x))

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

    # ens
    ens = np.array(initial_states) + np.random.randn(N, dim)
    #ens0 = 0.0 * np.random.randn(N, num_dynamic_states * dim * 2)  # TODO: check 0 # 3 variables (x,y,z) for each final and initial
    #ens = np.hstack([ens, ens0])
    ens_names = state_names# + dynamic_states + final_states
    ens = pd.DataFrame(ens, columns=ens_names)
    ens.to_csv(os.path.join(template_ws, "initial_x_ens.csv"))

    pst.pestpp_options['da_parameter_ensemble'] = "initial_x_ens.csv"
    pst.pestpp_options['da_subset_size'] = N
    pst.pestpp_options['da_num_reals'] = N
    # pst.pestpp_options['ies_init_lam'] = [100]
    pst.pestpp_options['ies_lambda_mults'] = 1
    pst.pestpp_options['lambda_scale_fac'] = 1
    pst.control_data.noptmax = 1
    pst.pestpp_options['da_use_mda'] = 'True'
    pst.pestpp_options['da_mda_init_fac'] = 1.0
    pst.pestpp_options['da_mda_dec_fac'] = 1.0

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

if __name__ == "__main__":
    lorenz96_basic_test()
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


