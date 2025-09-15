import os
import shutil
import platform
import string
import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import flopy
import pyemu

bin_path = os.path.join("test_bin")
if "linux" in platform.platform().lower():
    bin_path = os.path.join(bin_path,"linux")
elif "darwin" in platform.platform().lower():
    bin_path = os.path.join(bin_path,"mac")
else:
    bin_path = os.path.join(bin_path,"win")

#TODO:update mf6 bins for transport
bin_path = os.path.abspath(bin_path)
os.environ["PATH"] += os.pathsep + bin_path
mf_exe = os.path.join(bin_path,"mf6")
print(os.environ["PATH"])


bin_path = os.path.join("..","..","..","bin")
exe = ""
if "windows" in platform.platform().lower():
    exe = ".exe"
exe_path = os.path.join(bin_path, "pestpp-mou" + exe)


noptmax = 4
num_reals = 20
port = 4021
test_root = "henry"


org_d = os.path.join("henry","ex-gwt-henry-b")

def run_and_plot_results(cwd):
    pyemu.os_utils.run("mf6",cwd=cwd,verbose=True)
    plt_dir = "plots"
    if not os.path.exists(plt_dir):
        os.mkdir(plt_dir)
    ucn = flopy.utils.HeadFile(os.path.join(cwd,"trans.ucn"),text="concentration")
    hds = flopy.utils.HeadFile(os.path.join(cwd,"flow.hds"))

    wel_df = pd.read_csv(os.path.join(cwd, "flow.wel_stress_period_data_historic.txt"), header=None,
                         names=["l", "r", "c", "flux", "concen"], delim_whitespace=True)
    wel_df = wel_df.loc[wel_df.flux < 0, :]


    with PdfPages(os.path.join(plt_dir,"henry.pdf")) as pdf:
        for i,time in enumerate(ucn.get_times()):
            d = ucn.get_data(totim=time)
            fig,ax = plt.subplots(1,1,figsize=(8,3))
            cb = ax.imshow(d[:,0,:],interpolation="none",vmin=0.0,vmax=35.0)
            levels = [0.5]
            ax.contour(d[:,0,:],levels=levels,colors="w")
            d = hds.get_data(totim=time)
            levels = np.linspace(d.min(),d.mean(),5)
            #ax.contour(d[:,0,:],levels,colors="k",)
            ax.scatter(wel_df.c-1, wel_df.l-1, marker="^", color="k", label="extraction wells")
            ax.set_title("{0}".format(time))
            #plt.savefig(os.path.join(plt_dir,"{0:03d}_henry.pdf".format(i)))
            pdf.savefig()
            plt.close("all")
            print(time,d.min(),d.max())

def prep_model():
    new_d = os.path.join("henry","henry_temp")
    if os.path.exists(new_d):
        shutil.rmtree(new_d)
    shutil.copytree(org_d,new_d)

    #add zero-flux wells to the 3rd wel list file
    

def add_artrch():
    args = {}
    with open("artrch.dat",'r') as f:
        for line in f:
            raw = line.strip().split()
            args[raw[0]] = float(raw[1])

    org_wel_list_file = "flow.wel_stress_period_data_scenario_base.txt"
    new_wel_list_file = "flow.wel_stress_period_data_scenario.txt"
    lines = []
    with open(org_wel_list_file,'r') as f:
        for line in f:
            if len(line.strip()) > 0:
                lines.append(line)
    idist = max(40,int(args["dist"]))
    iwidth = max(1,int(args["width"]))
    iend = max(idist+iwidth,150)
    rate_per_well = args["rate"] / iwidth
    for i in range(iwidth):
        lines.append("1 1 {0} {1} {2}\n".
            format(idist + i,rate_per_well,args["concen"]))
    with open(new_wel_list_file,'w') as f:
        for line in lines:
            f.write(line)

def eval_add_artrch(test_d,dist=80,width=10,rate=1,concen=0.0,write_tpl=False):
    with open(os.path.join(test_d,"artrch.dat"),'w') as f:
        f.write("dist {0}\n".format(dist))
        f.write("width {0}\n".format(width))
        f.write("rate {0}\n".format(rate))
        f.write("concen {0}\n".format(concen))
    tpl_file = "artrch.dat.tpl"
    if write_tpl:
        with open(os.path.join(test_d, tpl_file), 'w') as f:
            f.write("ptf ~\n")
            f.write("dist ~   ar_dist    ~\n")
            f.write("width ~  ar_width    ~\n")
            f.write("rate   ~   ar_rate     ~\n")
            f.write("concen ~   ar_concen    ~\n")
    base_file = os.path.join(test_d, "flow.wel_stress_period_data_scenario_base.txt")
    if not os.path.exists(base_file):
        shutil.copy2(os.path.join(test_d, "flow.wel_stress_period_data_historic.txt"),
                     base_file)
    bd = os.getcwd()
    os.chdir(test_d)
    add_artrch()
    os.chdir(bd)
    if not write_tpl:
        run_and_plot_results(test_d)
    return tpl_file

def eval_head_at_artrch(d):
    cwd = os.getcwd()
    os.chdir(d)
    df = head_at_artrch()
    os.chdir(cwd)
    return df

def head_at_artrch():
    import flopy
    hds = flopy.utils.HeadFile("flow.hds")
    df = pd.read_csv("artrch.dat",header=None,names=["var","value"],index_col=0,delim_whitespace=True)

    start = int(df.loc["dist"])
    stop = start + int(df.loc["width"])
    hist_vals = hds.get_data((19,0))[0,0,start:stop]
    scen_vals = hds.get_data()[0,0,start:stop]
    d = scen_vals - hist_vals
    perd = 100. * (d / hist_vals)
    #print(hist_vals)
    #print(scen_vals)
    #print(d)
    #print(perd)
    with open("ar_heads.csv",'w') as f:
        f.write("qname,time,arhead\n")
        f.write("hist_mean,0,{0:15.6E}\n".format(hist_vals.mean()))
        f.write("scen_mean,1,{0:15.6E}\n".format(scen_vals.mean()))
        f.write("hist_max,0,{0:15.6E}\n".format(hist_vals.max()))
        f.write("scen_max,1,{0:15.6E}\n".format(scen_vals.max()))


def setup_pst():
    old_dir = os.path.join("henry","henry_temp")
    new_dir = os.path.join("henry","henry_template")
    sim = flopy.mf6.MFSimulation.load(sim_ws=old_dir)
    m = sim.get_model("flow")
    # setup a fake model grid
    delx = m.dis.delr.data
    botm = m.dis.botm.data[:,0,0]
    botm = np.insert(botm,0,1.0)
    dely = botm[0:-1] - botm[1:]
    #print(dely.shape)
    mg = flopy.discretization.StructuredGrid(dely,delx)

    v_k = pyemu.geostats.ExpVario(1.0,0.1,bearing=90,anisotropy=100)
    gs_k = pyemu.geostats.GeoStruct(variograms=v_k)

    pot_lim = 0.5 # potable limit, g/l salinity

    pf = pyemu.utils.PstFrom(old_dir,new_dir,spatial_reference=mg,
                             remove_existing=True,zero_based=False)


    # setup pars for k using aniso to represent vk
    pf.add_parameters("flow.npf_k.txt",par_type="grid",upper_bound=864*10.0,lower_bound=864*0.1,
                      par_name_base="k",pargp="k",
                      geostruct=gs_k,par_style="direct")

    # setup pars for porosity
    pf.add_parameters("trans.mst_porosity.txt", par_type="grid", upper_bound=0.4, lower_bound=0.3,
                     par_name_base="pr", pargp="pr",
                     geostruct=gs_k,par_style="direct")



    # copy the current historic well list file to a "base" file that will get modified...
    shutil.copy2(os.path.join(new_dir,"flow.wel_stress_period_data_historic.txt"),
                 os.path.join(new_dir,"flow.wel_stress_period_data_scenario_base.txt"))
    pf.add_parameters("flow.wel_stress_period_data_scenario_base.txt",par_type="grid",par_style="direct",
                      index_cols=[0,1,2],use_cols=3,pargp="wel",par_name_base="wel")

    pf.add_parameters("flow.ghb_stress_period_data_3.txt", par_type="constant", par_style="direct",
                      index_cols=[0, 1, 2], use_cols=3, pargp="stage", par_name_base="stage")

    # setup obs for all concentrations at the end of the 3 periods
    pump_filename,conc_filenames = eval_process_unc(new_dir)
    df = pd.read_csv(os.path.join(new_dir,pump_filename))
    cols = df.columns.to_list()
    pf.add_observations(pump_filename, index_cols=["time"], use_cols=cols[1:],
                        ofile_sep=",")
    ins_file = os.path.join(new_dir,"mean_concen.dat.ins")
    with open(ins_file,'w') as f:
        f.write("pif  ~\n")
        f.write("l1 w !mean_concen_time:0.0!\n")
    pf.add_observations_from_ins(ins_file,pst_path=".")

    for conc_filename in conc_filenames:
        pf.add_observations(conc_filename,index_cols=["time","k","j"],use_cols=["conc"],
                            ofile_sep=",")
        break
    pf.extra_py_imports.append("flopy")
    path_to_file = os.path.relpath(__file__)
    pf.add_py_function(path_to_file,"process_unc()",is_pre_cmd=False)
    pf.mod_sys_cmds.append("mf6")
    pf.add_py_function(path_to_file,"add_artrch()",is_pre_cmd=True)
    pf.tmp_files.append("flow.wel_stress_period_data_scenario.txt")

    # add artificial recharge basin dvs
    tpl_file = eval_add_artrch(new_dir, write_tpl=True)


    eval_head_at_artrch(new_dir)
    pf.add_observations("ar_heads.csv", ofile_sep=",", index_cols=[0,1], use_cols=[2], prefix="arhead")

    pf.add_py_function(path_to_file, "head_at_artrch()", is_pre_cmd=False)

    pf.build_pst()

    df = pf.pst.add_parameters(os.path.join(new_dir, tpl_file), pst_path=".")
    with open(os.path.join(new_dir, "risk.dat.tpl"), 'w') as f:
        f.write("ptf ~\n")
        f.write("_risk_ ~   _risk_    ~\n")
    pf.pst.add_parameters(os.path.join(new_dir, "risk.dat.tpl"), pst_path=".")

    par = pf.pst.parameter_data
    par.loc["_risk_","pargp"] = "dv_pars"
    par.loc["_risk_","parlbnd"] = 0.001
    par.loc["_risk_", "parubnd"] = 0.99
    # start at lower bound for cases without risk obj
    par.loc["_risk_", "parval1"] = 0.001
    par.loc["_risk_", "partrans"] = "none"

    par.loc[df.parnme,"pargp"] = "dv_pars"
    par.loc["ar_concen","parval1"] = 3.5
    par.loc["ar_concen", "parubnd"] = 17.0
    par.loc["ar_concen", "parlbnd"] = pot_lim
    par.loc["ar_concen", "partrans"] = "none"

    #dont let this go to zero
    par.loc["ar_rate", "parval1"] = 5.5
    par.loc["ar_rate", "parubnd"] = 12.0
    par.loc["ar_rate", "parlbnd"] = 0.001
    par.loc["ar_rate", "partrans"] = "none"


    par.loc["ar_dist", "parval1"] = 80
    par.loc["ar_dist", "parubnd"] = 140
    par.loc["ar_dist", "parlbnd"] = 1
    par.loc["ar_dist", "partrans"] = "none"

    par.loc["ar_width", "parval1"] = 10
    par.loc["ar_width", "parubnd"] = 20
    par.loc["ar_width", "parlbnd"] = 2
    par.loc["ar_width", "partrans"] = "fixed"

    wel_par = par.loc[par.pargp=="wel",:]
    wpar = wel_par.loc[wel_par.parval1>0,"parnme"]

    par.loc[wpar, "partrans"] = "log"
    par.loc[wpar,"pargp"] = "wel_rch"
    par.loc[wpar, "parubnd"] = par.loc[wpar,"parval1"] * 1.1
    par.loc[wpar, "parlbnd"] = par.loc[wpar,"parval1"] * 0.9

    wpar = wel_par.loc[wel_par.parval1<0,"parnme"]
    par.loc[wpar, "partrans"] = "none"
    par.loc[wpar, "pargp"] = "dv_pars"
    par.loc[wpar, "parubnd"] = 0.0
    par.loc[wpar, "parlbnd"] = -3.0
    # this one is the objective
    pf.pst.add_pi_equation(wpar.to_list(),pilbl="pump_rate",obs_group="less_than")
    # this one is the constraint
    #tot = par.loc[wpar,"parval1"].sum()
    #pf.pst.add_pi_equation(wpar.to_list(), pilbl="constraint_pump_rate", obs_group="less_than",rhs=tot)

    stage_par = par.loc[par.pargp == "stage", "parnme"].values[0]
    par.loc[stage_par, "partrans"] = "fixed"
    #par.loc[stage_par,"parubnd"] = 1.05
    #par.loc[stage_par, "parlbnd"] = 1.00
    #par.loc[stage_par,"pargp"] = "dv_pars"


    pf.pst.control_data.noptmax = 0

    #pf.pst.add_pi_equation([stage_par], obs_group="greater_than", pilbl=stage_par)
    #pf.pst.add_pi_equation(["ar_width"],obs_group="less_than",pilbl="ar_width")
    pf.pst.add_pi_equation(["ar_rate"], obs_group="less_than",pilbl="ar_rate")
    pf.pst.add_pi_equation(["ar_concen"], obs_group="greater_than",pilbl="ar_concen")
    pf.pst.add_pi_equation(["ar_dist"], obs_group="greater_than", pilbl="ar_dist")

    pf.pst.add_pi_equation(["_risk_"], obs_group="greater_than",pilbl="_risk_")
    #pf.pst.pestpp_options["mou_objectives"] = ["ar_width","ar_rate","ar_concen","pump_rate", "_risk_"]
    pf.pst.pestpp_options["mou_objectives"] = ["_risk_","ar_dist","ar_rate","ar_concen","pump_rate"]

    pf.pst.pestpp_options["opt_dec_var_groups"] = "dv_pars"
    pf.pst.pestpp_options["panther_echo"] = True
    pf.pst.pestpp_options["mou_risk_objective"] = True
    pf.pst.pestpp_options["mou_generator"] = "de"
    pf.pst.pestpp_options["mou_population_size"] = 100

    pf.pst.try_parse_name_metadata()
    obs = pf.pst.observation_data
    mx_time = obs.time.max()
    # use the last output time as the constraint on concen at the pumping wells
    #obs.loc[obs.time==mx_time,"obsval"] = pot_lim
    #obs.loc[obs.time==mx_time, "obgnme"] = "less_than"
    obs.loc["mean_concen_time:0.0","obsval"] = pot_lim
    obs.loc["mean_concen_time:0.0", "obgnme"] = "less_than"

    obs.loc["oname:arhead_otype:lst_usecol:arhead_qname:scen_max_time:1","weight"] = 1.0
    obs.loc["oname:arhead_otype:lst_usecol:arhead_qname:scen_max_time:1", "obsval"] = 1.025
    obs.loc["oname:arhead_otype:lst_usecol:arhead_qname:scen_max_time:1", "obgnme"] = "less_than"

    pf.pst.write(os.path.join(new_dir,"henry.pst"))
    pe = pf.draw(100,use_specsim=True)
    pe.to_binary(os.path.join(new_dir,"prior.jcb"))

    pyemu.os_utils.run("{0} henry.pst".format(exe_path),cwd=new_dir)

def eval_process_unc(test_d):
    bd = os.getcwd()
    os.chdir(test_d)
    ret = process_unc()
    os.chdir(bd)
    return ret

def process_unc():
    ucn = flopy.utils.HeadFile("trans.ucn", text="concentration")
    df = pd.read_csv(os.path.join("flow.wel_stress_period_data_scenario.txt"),
                    delim_whitespace=True, header=None, names=["l", "r", "c", "flux", "concen"])

    df = df.loc[df.flux<=0, :]
    d = ucn.get_alldata()
    # print(d.shape)
    wel_ucn = {}
    k,j = {},{}
    flx_dict = {}
    for l,r,c,flx in zip(df.l,df.r,df.c,df.flux):
        dd = d[:,l-1,r-1,c-1]
        name = "well_{0:03d}_{1:03d}_{2:03d}".format(l,r,c)
        wel_ucn[name] = dd
        k[name] = l-1
        j[name] = c-1
        flx_dict[name] = -1. * flx

    df = pd.DataFrame(wel_ucn,index=np.round(ucn.get_times(),5))
    df.loc[:,"time"] = df.index
    cols = df.columns.to_list()
    cols.sort()
    cols.remove("time")
    cols.insert(0,"time")
    df = df.loc[:,cols]
    last_concen = df.loc[df.time==df.time.max(),:].iloc[0].to_dict()
    tflx,tmas = 0,0
    for n,f in flx_dict.items():
        tflx += f
        tmas += f * last_concen[n]
    mn_concen = tmas / tflx
    with open("mean_concen.dat",'w') as f:
        f.write("mean_concen {0}\n".format(mn_concen))

    pump_file_name = "pump_well_concen.csv"
    df.to_csv(pump_file_name,index=False)

    col = np.zeros((d.shape[1],d.shape[3]))
    lay = np.zeros((d.shape[1],d.shape[3]))
    for k in range(d.shape[1]):
        col[k,:] = np.arange(d.shape[3])
        lay[k,:] = k
    col = col.flatten()
    lay = lay.flatten()

    filenames = []
    # for i,time in enumerate(ucn.get_times()):
    #     dd = d[i,:,:,:].flatten()
    #     time = np.round(time,5)
    #     df = pd.DataFrame({"conc":dd})
    #     df.loc[:,"time"] = time
    #     df.loc[:,"j"] = col
    #     df.loc[:,"k"] = lay
    #     filename = "concen_{0}.csv".format(i)
    #     df.loc[:,["time","k","j","conc"]].to_csv(filename,index=False)
    #     filenames.append(filename)

    # dd = d[-1,:,:,:].flatten()
    # time = np.round(ucn.get_times()[-1],5)
    # df = pd.DataFrame({"conc":dd})
    # df.loc[:,"time"] = time
    # df.loc[:,"j"] = col
    # df.loc[:,"k"] = lay
    # filename = "concen_{0}.csv".format(d.shape[0]-1)
    # df.loc[:,["time","k","j","conc"]].to_csv(filename,index=False)
    # filenames.append(filename)

    return pump_file_name,filenames

def plot_pr_real():
    pst = pyemu.Pst(os.path.join("henry_template","henry.pst"))
    pe = pyemu.ParameterEnsemble.from_binary(pst,os.path.join("henry_template","prior.jcb"))
    pst.parameter_data.loc[:,"parval1"] = pe.loc[pe.index[0],pst.par_names]
    pst.write(os.path.join("henry_template","test.pst"))
    pyemu.os_utils.run("pestpp-ies test.pst",cwd="henry_template")
    pst.write_input_files("henry_template")
    k = np.log10(np.loadtxt(os.path.join("henry_template","org","flow.npf_k.txt")))
    plt.imshow(k)
    plt.show()

def start_workers_for_debug(with_master=True):
    t_d = os.path.join("henry", "henry_template")
    m_d = os.path.join("henry","henry_master_chance_risk_obj")
    if with_master:
        if os.path.exists(m_d):
            shutil.rmtree(m_d)
        shutil.copytree(t_d,m_d)
        pst = pyemu.Pst(os.path.join(m_d,"henry.pst"))
        pst.control_data.noptmax = 200
        pst.pestpp_options["opt_par_stack"] = "prior.jcb"
        pst.pestpp_options["opt_stack_size"] = 30
        pst.pestpp_options["opt_recalc_chance_every"] = 1000
        pst.pestpp_options["mou_population_size"] = 100
        pst.pestpp_options["opt_chance_points"] = "all"
        pst.pestpp_options["opt_risk"] = 0.70

        pst.write(os.path.join(m_d,"henry.pst"))
    pyemu.os_utils.start_workers(t_d, exe_path, "henry.pst",
                                  num_workers=12, worker_root="henry",
                                  port=4004)

def plot_results(m_d,risk_thres=0.0,tag=""):
    fs=8
    # plt_d = "henry_results"
    # if os.path.exists(plt_d):
    #     shutil.rmtree(plt_d)
    # os.mkdir(plt_d)
    df_arc = pd.read_csv(os.path.join(m_d,"henry.pareto.archive.summary.csv"))

    pst = pyemu.Pst(os.path.join(m_d,"henry.pst"))
    par = pst.parameter_data
    par = par.loc[par.pargp == "dv_pars", :].copy()
    rate_pars = par.loc[par.parnme.apply(lambda x: not x.startswith("ar") and not "risk" in x), "parnme"]
    obj_names = pst.pestpp_options["mou_objectives"].split(',')

    bnd_dict = {dv: [par.loc[dv, "parlbnd"], par.loc[dv, "parubnd"]] for dv in obj_names if dv in pst.par_names}
    bnd_dict["pump_rate"] = [par.loc[rate_pars, "parubnd"].sum() * -1., par.loc[rate_pars, "parlbnd"].sum() * -1.]


    if "_risk_" in obj_names and "riskobj" not in m_d:
        obj_names.remove("_risk_")


    dmn = 0
    dmx = 160

    gens = df_arc.generation.unique()
    gens.sort()
    print(gens)
    obj_names_dict = {"ar_width": "basin width (columns)",
                      "ar_rate": "recharge rate ($\\frac{m^3}{d}$)",
                      "pump_rate": "extraction rate ($\\frac{m^3}{d}$)",
                      "ar_concen" : "recharge salinity ($\\frac{g}{l}$)",
                      "_risk_":"risk","stage_inst:0_usecol:3_direct":"coastal stage ($m$)",
                      "ar_dist":"basin edge (column)"}

    cmap = plt.get_cmap("jet")
    obj_names.sort()

    with PdfPages(os.path.join(m_d,os.path.split(m_d)[-1]+tag+".pdf")) as pdf:
        for gen in [gens[-1]]:
            ax_count = 0
            #df_gen = df_arc.loc[df_arc.generation==gen,:].copy()
            df_gen = pd.read_csv(os.path.join(m_d,"henry.archive.dv_pop.csv"))
            df_gen.loc[:,"front"] = df_arc.loc[df_gen.index,"nsga2_front"]

            # flip the rate pars
            df_gen.loc[:, rate_pars] *= -1.
            df_gen.loc[:,"pump_rate"] = df_gen.loc[:, rate_pars].sum(axis=1)

            # only show solutions with some min amount of pumping
            df_gen = df_gen.loc[df_gen.pump_rate > 2.1]
            print(df_gen.loc[:,"ar_rate"])

     

            # only with risk averse
            if "_risk_" in obj_names:
                df_gen = df_gen.loc[df_gen._risk_ > risk_thres,:]

            if df_gen.shape[0] == 0:
                continue
            norm = Normalize(df_gen.loc[:, "ar_concen"].min(), df_gen.loc[:, "ar_concen"].max())

            fig = plt.figure(figsize=(7.5,9))

            gs = GridSpec(len(obj_names)+1,len(obj_names),figure=fig)
            ax = fig.add_subplot(gs[0,:])
            mx = 0
            mn = 1.0e+10
            for d,w,r,c,t in zip(df_gen.ar_dist,df_gen.ar_width,df_gen.ar_rate,
                               df_gen.ar_concen,df_gen.pump_rate,):
                rect = Rectangle((d,0),w,r/t,facecolor=cmap(norm(c)),edgecolor="k",alpha=0.25)
                ax.add_patch(rect)
                #mpt = d + (w/2.)
                mx = max(mx,r/t)
                mn = min(mn,r/t)
                #ax.scatter([mpt],[r],marker="o",s=(100 * (t-df_gen.pump_rate.min())/(df_gen.pump_rate.max()-df_gen.pump_rate.min())),color="k",zorder=10,alpha=0.5)

            ax.set_xlim(dmn,dmx)
            #ax.set_ylim(par.loc["ar_rate","parlbnd"],par.loc["ar_rate","parubnd"])
            ax.set_ylim(mn,mx)
            ax.set_xlabel("model column",fontsize=fs)
            ax.set_ylabel("ratio of recharge\nto extraction",fontsize=fs)
            ax.plot([60,60],ax.get_ylim(),"k--")
            ax.text(65,(mn+mx)/3.7,"extraction\nwells",rotation=90,va='bottom',ha="center",fontsize=fs)
            cb = plt.colorbar(mpl.cm.ScalarMappable(norm=norm,cmap=cmap),ax=ax,alpha=0.25)
            cb.set_label(label=obj_names_dict["ar_concen"],size=fs)
            cb.ax.tick_params(labelsize=fs)
            #cb = plt.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap))
            #cb.set_label("recharge concentration")


            ax.set_title("{0}) optimal artificial recharge basin designs, {2} feasible non-dominated solutions".\
                         format(string.ascii_uppercase[ax_count],gen,df_gen.shape[0]),loc="left",fontsize=fs)
            ax.tick_params(axis='both', which='major', labelsize=fs)
            ax_count += 1
            for i,o1 in enumerate(obj_names):
                for j in range(i+1):
                    ax = fig.add_subplot(gs[i+1,j])
                    ax.set_title("{0})".format(string.ascii_uppercase[ax_count]),loc="left",fontsize=fs)
                    ax_count += 1
                    if i == j:
                        ax.hist(df_gen.loc[:,o1],facecolor="0.5",edgecolor="none",alpha=0.5)
                        ax.set_yticks([])
                        ax.set_xlabel(obj_names_dict[o1],fontsize=fs)
                        #ax.set_xlim(bnd_dict[o1])
                        ax.set_ylabel("increasing probability",fontsize=fs)
                    else:
                        if "_risk_" in obj_names:
                            ax.scatter(df_gen.loc[:, obj_names[j]], df_gen.loc[:, o1], marker=".", c=df_gen._risk_.values)
                        else:
                            ax.scatter(df_gen.loc[:,obj_names[j]],df_gen.loc[:,o1],marker=".",color="0.5")
                        ax.set_xlabel(obj_names_dict[obj_names[j]],fontsize=fs)
                        ax.set_ylabel(obj_names_dict[o1],fontsize=fs)
                    ax.tick_params(axis='both', which='major', labelsize=fs)
                        #ax.set_xlim(bnd_dict[obj_names[j]])
                for j in range(i,len(obj_names)):

                    if i != j:
                        ax = fig.add_subplot(gs[i + 1, j])
                        ax.set_title("{0})".format(string.ascii_uppercase[ax_count]), loc="left",fontsize=fs)
                        ax_count += 1
                        if "_risk_" in obj_names:
                            ax.scatter(df_gen.loc[:, obj_names[j]], df_gen.loc[:, o1], marker=".", c=df_gen._risk_.values)
                        else:
                            ax.scatter(df_gen.loc[:,obj_names[j]],df_gen.loc[:,o1],marker=".",color="0.5")
                        #ax.scatter(df_gen.loc[:,obj_names[j]],df_gen.loc[:,o1],marker=".",color="0.5")
                        #ax.set_xlabel(obj_names[j])
                        #ax.set_ylabel(o1)
                        ax.set_xlabel(obj_names_dict[obj_names[j]],fontsize=fs)
                        ax.set_ylabel(obj_names_dict[o1],fontsize=fs)
                        ax.tick_params(axis='both', which='major', labelsize=fs)



            plt.tight_layout()

            #plt.show()
            pdf.savefig()
            plt.close(fig)
            #if gen > 10:
            #    break

def invest():
    m_d = os.path.join("henry","henry_master_riskobj_all_once")
    pst = pyemu.Pst(os.path.join(m_d,"henry.pst"))
    onames = pst.nnz_obs_names
    df = pd.read_csv(os.path.join(m_d,"henry.0.nested.obs_stack.csv"))
    df.loc[:,"member"] = df.real_name.apply(lambda x: x.split("||")[1])
    df.loc[:, "real"] = df.real_name.apply(lambda x: x.split("||")[0])
    mnames = df.member.unique()
    mnames.sort()
    rnames = df.real.unique()
    rnames.sort()
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(os.path.join(m_d,"stack_summary.pdf")) as pdf:
        for oname in onames:
            fig,ax = plt.subplots(1,1,figsize=(6,6))
            for rname in rnames:
                dfm = df.loc[df.real == rname]
                dfm.loc[:,oname].hist(ax=ax,alpha=0.5,bins=20)
            ax.set_title(oname)
            pdf.savefig()
            plt.close(fig)




def run_mou(risk_obj=False,chance_points="single",risk=0.5,stack_size=100,
            num_workers=12,pop_size=100,tag="",recalc_every=100000,noptmax=100):
    t_d = os.path.join("henry","henry_template")
    pst = pyemu.Pst(os.path.join(t_d,"henry.pst"))
    pst.pestpp_options["opt_par_stack"] = "prior.jcb"
    pst.pestpp_options["mou_risk_objective"] = risk_obj
    pst.pestpp_options["opt_recalc_chance_every"] = recalc_every
    pst.pestpp_options["opt_chance_points"] = chance_points
    pst.pestpp_options["opt_risk"] = risk
    pst.pestpp_options["mou_population_size"] = pop_size
    pst.pestpp_options["opt_stack_size"] = stack_size
    if risk == 0.5:
        objs = pst.pestpp_options["mou_objectives"].split(",")
        if "_risk_" in objs:
            objs.remove("_risk_")
        pst.pestpp_options["mou_objectives"] = objs
    else:
        objs = pst.pestpp_options["mou_objectives"].split(",")
        if "_risk_" not in objs:
            objs.append("_risk_")
        pst.pestpp_options["mou_objectives"] = objs
    pst.control_data.noptmax = noptmax
    pst.write(os.path.join(t_d,"henry.pst"))

    m_d = t_d.replace("template","master")
    if len(tag) > 0:
        m_d += "_" + tag
    pyemu.os_utils.start_workers(t_d,"pestpp-mou","henry.pst",
                                 num_workers=num_workers,master_dir=m_d,
                                 verbose=True, worker_root="henry")



def plot_domain(cwd,include_pred=False):
    fs = 10
    wel_df = pd.read_csv(os.path.join(cwd,"flow.wel_stress_period_data_historic.txt"),header=None,
                         names=["l","r","c","flux","concen"],delim_whitespace=True)
    wel_df = wel_df.loc[wel_df.flux<0,:]
    swel_df = pd.read_csv(os.path.join(cwd,"flow.wel_stress_period_data_scenario.txt"),header=None,
                         names=["l","r","c","flux","concen"],delim_whitespace=True)
    swel_df = swel_df.loc[swel_df.c>66,:]

   

    #pyemu.os_utils.run("mf6", cwd=cwd, verbose=True)
    plt_dir = "plots"
    if not os.path.exists(plt_dir):
        os.mkdir(plt_dir)
    ucn = flopy.utils.HeadFile(os.path.join(cwd, "trans.ucn"), text="concentration")
    hds = flopy.utils.HeadFile(os.path.join(cwd, "flow.hds"))
    print(ucn.get_times())
    d = ucn.get_data(kstpkper=(0,1))
    if include_pred:
        fig, axes = plt.subplots(3, 1, figsize=(8.5, 7.5))
        levels = [0.5,3.5,17.5]
    else:
        fig, axes = plt.subplots(2, 1, figsize=(8.5, 5))
        levels = [0.5]

    ax = axes[0]
    ax.set_title("A) pre-developement salinity concentration",loc="left",fontsize=fs)
    ax.set_ylabel("layer (elevation ($m$))", fontsize=fs)
    #ax.set_ylabel("layer",fontsize=fs)
    #ax.set_xlabel("column",fontsize=fs)
    cb = ax.imshow(d[:, 0, :], interpolation="none", vmin=0.0, vmax=35.0)
    cb = plt.colorbar(cb,ax=ax)
    cb.set_label("salinity ($\\frac{g}{l}$)",fontsize=fs)
    cb.ax.tick_params(axis='both', which='major', labelsize=fs)
    ax.contour(d[:, 0, :], levels=levels, colors="w",label="potable salinity limit")
    ax.scatter(wel_df.c-1,wel_df.l-1,marker="^",color="0.5",label="extraction wells")
    ylim = ax.get_ylim()
    ax.plot([0,0],ylim,"m",lw=5,label="upgradient boundary")
    xmx = ax.get_xlim()[1]
    ax.plot([xmx, xmx], ylim, "r", lw=5, label="coastal boundary")
    #ax.plot([40,55],[0,0],"g",lw=5,label="feasible artifical recharge basin")
    ax.legend(loc="upper left",fontsize=fs)
    ax.tick_params(axis='both', which='major', labelsize=fs)
    #d = hds.get_data(totim=time)
    #levels = np.linspace(d.min(), d.mean(), 5)
    #ax.contour(d[:, 0, :], levels, colors="k", )
    #plt.savefig("henry_domain.pdf")
    ax = axes[1]
    ax.set_title("B) salinity concentration after historic water use", loc="left",fontsize=fs)
    ax.set_ylabel("layer (elevation ($m$))",fontsize=fs)
    ax.set_xlabel("column\n (x distance ($m$))",fontsize=fs)

    d = ucn.get_data(kstpkper=(0, 20))
    cb = ax.imshow(d[:, 0, :], interpolation="none", vmin=0.0, vmax=35.0)
    cb = plt.colorbar(cb, ax=ax)
    cb.set_label("salinity ($\\frac{g}{l}$)",fontsize=fs)
    cb.ax.tick_params(axis='both', which='major', labelsize=fs)

    ax.contour(d[:, 0, :], levels=levels, colors="w", label="potable salinity limit")
    ax.scatter(wel_df.c-1, wel_df.l-1, marker="^", color="0.5", label="extraction wells")
    ylim = ax.get_ylim()
    ax.plot([0, 0], ylim, "m", lw=5, label="upgradient boundary")
    xmx = ax.get_xlim()[1]
    ax.plot([xmx, xmx], ylim, "r", lw=5, label="coastal boundary")
    #ax.plot([40, 55], [0, 0], "g", lw=5, label="feasible artifical recharge basin")
    ax.legend(loc="upper left",fontsize=fs)
    ax.tick_params(axis='both', which='major', labelsize=fs)
    if include_pred:
        ax = axes[2]
        ax.set_title("C) salinity concentration after predictive", loc="left",fontsize=fs)
        ax.set_ylabel("layer (elevation ($m$))",fontsize=fs)
        ax.set_xlabel("column\n (x distance ($m$))",fontsize=fs)

        d = ucn.get_data()
        cb = ax.imshow(d[:, 0, :], interpolation="none", vmin=0.0, vmax=35.0)
        cb = plt.colorbar(cb, ax=ax)
        cb.set_label("salinity ($\\frac{g}{l}$)",fontsize=fs)
        cb.ax.tick_params(axis='both', which='major', labelsize=fs)

        
        ax.contour(d[:, 0, :], levels=levels, colors="w", label="potable salinity limit")
        ax.scatter(wel_df.c-1, wel_df.l-1, marker="^", color="0.5", label="extraction wells")
        ax.scatter(swel_df.c-1, swel_df.l -1, marker="*", color="r", label="recharge basin cells",s=100)
        
        ylim = ax.get_ylim()
        ax.plot([0, 0], ylim, "m", lw=5, label="upgradient boundary")
        xmx = ax.get_xlim()[1]
        ax.plot([xmx, xmx], ylim, "r", lw=5, label="coastal boundary")
        #ax.plot([40, 55], [0, 0], "g", lw=5, label="feasible artifical recharge basin")
        #ax.legend(loc="upper left",fontsize=fs)
        ax.tick_params(axis='both', which='major', labelsize=fs)

    zticks = np.arange(0, 50, 10)
    zdist = zticks * 0.025
    zlabs = ["{0:1.0f} ({1:2.1f})".format(t, d) for t, d in zip(zticks, zdist[::-1])]
    for ax in axes.flatten():
        ax.set_yticks(zticks)
        ax.set_yticklabels(zlabs)
    
    xticks = np.arange(0, 180, 20)
    xdist = xticks * 0.025
    xlabs = ["{0:1.0f}\n({1:2.1f})".format(t, d) for t, d in zip(xticks, xdist)]
    for ax in axes.flatten():
        ax.set_xticks(xticks)
        ax.set_xticklabels(xlabs)
    

    plt.tight_layout()
    plt.savefig("henry_domain.pdf")
    plt.close(fig)


def extract_and_plot_solution():
    m_d = os.path.join("henry","henry_master_riskobj_all_once")
    #m_d = os.path.join("henry","henry_master_95_single_once")
    pst = pyemu.Pst(os.path.join(m_d,"henry.pst"))
    df_arc = pd.read_csv(os.path.join(m_d,"henry.pareto.archive.summary.csv"))
    df_arc = df_arc.loc[df_arc.generation==df_arc.generation.max(),:]
    df_arc.loc[:,"member"] = df_arc.member.str.lower()
    df_arc.sort_values(by="ar_dist",ascending=False,inplace=True)
    df_dv = pd.read_csv(os.path.join(m_d,"henry.archive.dv_pop.csv"))
    #dv_vals = df_dv.loc[df_arc.member.iloc[0],:]
    #df_dv = df_dv.loc[df_dv._risk_>0.95,:]
    df_dv.sort_values(by="_risk_",ascending=False,inplace=True)
    print(df_dv.loc[:,"ar_rate"])
   
    pst.parameter_data.loc[:,"parval1"] = df_dv.loc[df_dv.index[2],pst.par_names]

    pst.control_data.noptmax = 0
    pst.write(os.path.join(m_d,"test.pst"))
    pyemu.os_utils.run("{0} test.pst".format(exe_path),cwd=m_d)
    run_and_plot_results(m_d)

    plot_domain(m_d,include_pred=True)

def simple_henry_test():
    prep_model()
    run_and_plot_results(os.path.join("henry", "henry_temp"))

    plot_domain(os.path.join("henry", "henry_temp"))
    setup_pst()
    run_mou(risk=0.65,tag="65_single_once",num_workers=10,noptmax=10,pop_size=50,stack_size=30)



if __name__ == "__main__":
    simple_henry_test()
    exit()
    #eval_process_unc(os.path.join("henry", "henry_template"))
    #shutil.copy2(os.path.join("..", "bin", "win", "pestpp-mou.exe"), os.path.join("..", "bin", "pestpp-mou.exe"))
    #shutil.copy2(os.path.join("..","exe","windows","x64","Debug","pestpp-mou.exe"),os.path.join("..","bin","pestpp-mou.exe"))
    #prep_model()
    #run_and_plot_results(os.path.join("henry", "henry_temp"))

    #plot_domain(os.path.join("henry", "henry_master_deter_100gen"),include_pred=True)
    extract_and_plot_solution()
    #setup_pst()
    #simple_henry_test()
    #run_mou(risk=0.95,tag="95_single_once",num_workers=40,noptmax=250)
    #run_mou(risk=0.5, tag="deter", num_workers=40, noptmax=100,pop_size=100)
    #run_mou(risk=0.95,risk_obj=True,tag="riskobj_all_once",chance_points="all",
    #        num_workers=40,noptmax=500,pop_size=100)

    #run_mou(risk=0.95,tag="95_single_once",num_workers=40,noptmax=100)
    #run_mou(risk=0.5, tag="deter", num_workers=40, noptmax=100,pop_size=250)
    #run_mou(risk=0.95,risk_obj=True,tag="riskobj_all_once",num_workers=40,noptmax=500,chance_points="all")
    #run_mou(risk=0.95,tag="95_all_once",chance_points="all",num_workers=40,noptmax=400)
    #run_mou(risk=0.95,tag="95_all_100th",chance_points="all",recalc_every=100,num_workers=40,noptmax=500)


    #plot_results(os.path.join("henry","henry_master_deter"))
    #plot_results(os.path.join("henry", "henry_master_95_single_once"))
    plot_results(os.path.join("henry", "henry_master_riskobj_all_once"))
    #plot_results(os.path.join("henry", "henry_master_riskobj_all_once"),risk_thres=0.65,
    #             tag="_riskaverse")

    #xtract_and_plot_solution()
    #invest()