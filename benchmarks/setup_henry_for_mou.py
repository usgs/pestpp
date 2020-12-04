import os
import shutil
import platform
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize
from matplotlib.backends.backend_pdf import PdfPages
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

bin_path = os.path.abspath("test_bin")
os.environ["PATH"] += os.pathsep + bin_path


bin_path = os.path.join("..","..","..","bin")
exe = ""
if "windows" in platform.platform().lower():
    exe = ".exe"
exe_path = os.path.join(bin_path, "pestpp-mou" + exe)


noptmax = 4
num_reals = 20
port = 4021
test_root = "mou_tests"


org_d = os.path.join("mou_tests","ex-gwt-henry-b")

def run_and_plot_results(cwd):
    pyemu.os_utils.run("mf6",cwd=cwd)
    plt_dir = "plots"
    if not os.path.exists(plt_dir):
        os.mkdir(plt_dir)
    ucn = flopy.utils.HeadFile(os.path.join(cwd,"trans.ucn"),text="concentration")
    hds = flopy.utils.HeadFile(os.path.join(cwd,"flow.hds"))
    with PdfPages(os.path.join(plt_dir,"henry.pdf")) as pdf:
        for i,time in enumerate(ucn.get_times()):
            d = ucn.get_data(totim=time)
            fig,ax = plt.subplots(1,1,figsize=(8,3))
            cb = ax.imshow(d[:,0,:],interpolation="none",vmin=0.0,vmax=35.0)
            levels = [35 * f for f in [0.0035, 0.5]]
            ax.contour(d[:,0,:],levels=levels,colors="w")
            d = hds.get_data(totim=time)
            levels = np.linspace(d.min(),d.mean(),5)
            ax.contour(d[:,0,:],levels,colors="k",)

            ax.set_title("{0}".format(time))
            #plt.savefig(os.path.join(plt_dir,"{0:03d}_henry.pdf".format(i)))
            pdf.savefig()
            plt.close("all")
            print(time,d.min(),d.max())

def prep_model():
    new_d = os.path.join("mou_tests","henry_temp")
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

def test_add_artrch(test_d,dist=80,width=10,rate=1,concen=0.0,write_tpl=False):
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

def setup_pst():
    old_dir = os.path.join("mou_tests","henry_temp")
    new_dir = os.path.join("mou_tests","henry_template")
    sim = flopy.mf6.MFSimulation.load(sim_ws=old_dir)
    m = sim.get_model("flow")
    # setup a fake model grid
    delx = m.dis.delr.data
    botm = m.dis.botm.data[:,0,0]
    botm = np.insert(botm,0,1.0)
    dely = botm[0:-1] - botm[1:]
    #print(dely.shape)
    mg = flopy.discretization.StructuredGrid(dely,delx)

    v_k = pyemu.geostats.ExpVario(1.0,0.1)
    gs_k = pyemu.geostats.GeoStruct(variograms=v_k)

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
    pump_filename,conc_filenames = test_process_unc(new_dir)
    df = pd.read_csv(os.path.join(new_dir,pump_filename))
    cols = df.columns.to_list()
    pf.add_observations(pump_filename, index_cols=["time"], use_cols=cols[1:],
                        ofile_sep=",")

    for conc_filename in conc_filenames:
        pf.add_observations(conc_filename,index_cols=["time","k","j"],use_cols=["conc"],
                            ofile_sep=",")
        break
    pf.extra_py_imports.append("flopy")
    pf.add_py_function("setup_henry_for_mou.py","process_unc()",is_pre_cmd=False)
    pf.mod_sys_cmds.append("mf6")
    pf.add_py_function("setup_henry_for_mou.py","add_artrch()",is_pre_cmd=True)
    pf.tmp_files.append("flow.wel_stress_period_data_scenario.txt")
    pf.build_pst()

    # add artificial recharge basin dvs
    tpl_file = test_add_artrch(new_dir, write_tpl=True)
    df = pf.pst.add_parameters(os.path.join(new_dir,tpl_file),pst_path=".")
    with open(os.path.join(new_dir,"risk.dat.tpl"),'w') as f:
        f.write("ptf ~\n")
        f.write("_risk_ ~   _risk_    ~\n")
    pf.pst.add_parameters(os.path.join(new_dir,"risk.dat.tpl"),pst_path=".")

    par = pf.pst.parameter_data
    par.loc["_risk_","pargp"] = "dv_pars"
    par.loc["_risk_","parlbnd"] = 0.001
    par.loc["_risk_", "parubnd"] = 0.999
    par.loc["_risk_", "parval1"] = 0.5
    par.loc["_risk_", "partrans"] = "none"

    par.loc[df.parnme,"pargp"] = "dv_pars"
    par.loc["ar_concen","parval1"] = 0.0
    par.loc["ar_concen", "parubnd"] = 35.0
    par.loc["ar_concen", "parlbnd"] = 0.0
    par.loc["ar_concen", "partrans"] = "none"

    par.loc["ar_rate", "parval1"] = 1.0
    par.loc["ar_rate", "parubnd"] = 2.5
    par.loc["ar_rate", "parlbnd"] = 0.1
    par.loc["ar_rate", "partrans"] = "none"

    par.loc["ar_dist", "parval1"] = 80
    par.loc["ar_dist", "parubnd"] = 120
    par.loc["ar_dist", "parlbnd"] = 40
    par.loc["ar_dist", "partrans"] = "none"

    par.loc["ar_width", "parval1"] = 10
    par.loc["ar_width", "parubnd"] = 20
    par.loc["ar_width", "parlbnd"] = 1
    par.loc["ar_width", "partrans"] = "none"

    wel_par = par.loc[par.pargp=="wel",:]
    wpar = wel_par.loc[wel_par.parval1>0,"parnme"]
    par.loc[wpar, "partrans"] = "none"
    par.loc[wpar,"pargp"] = "wel_rch"
    par.loc[wpar, "parubnd"] = 0.11
    par.loc[wpar, "parlbnd"] = 0.05
    wpar = wel_par.loc[wel_par.parval1<0,"parnme"]
    par.loc[wpar, "partrans"] = "none"
    par.loc[wpar, "pargp"] = "dv_pars"
    par.loc[wpar, "parubnd"] = 0.0
    par.loc[wpar, "parlbnd"] = -2.0
    pf.pst.add_pi_equation(wpar.to_list(),pilbl="pump_rate",obs_group="less_than")

    stage_par = par.loc[par.pargp == "stage", "parnme"]
    par.loc[stage_par, "partrans"] = "fixed"

    pf.pst.control_data.noptmax = 0
    pf.pst.add_pi_equation(["ar_width"],obs_group="less_than",pilbl="ar_width")
    pf.pst.add_pi_equation(["ar_rate"], obs_group="less_than",pilbl="ar_rate")
    pf.pst.add_pi_equation(["ar_concen"], obs_group="greater_than",pilbl="ar_concen")
    pf.pst.add_pi_equation(["_risk_"], obs_group="greater_than",pilbl="_risk_")
    pf.pst.pestpp_options["mou_objectives"] = ["ar_width","ar_rate","ar_concen","pump_rate", "_risk_"]

    pf.pst.pestpp_options["opt_dec_var_groups"] = "dv_pars"
    pf.pst.pestpp_options["panther_echo"] = True
    pf.pst.pestpp_options["mou_risk_objective"] = True
    pf.pst.pestpp_options["mou_generator"] = "de"
    pf.pst.pestpp_options["mou_population_size"] = 100

    pf.pst.try_parse_name_metadata()
    obs = pf.pst.observation_data
    mx_time = obs.time.max()
    # use the last output time as the constraint on concen at the pumping wells
    obs.loc[obs.time==mx_time,"obsval"] = 0.01
    obs.loc[obs.time==mx_time, "obgnme"] = "less_than"

    pf.pst.write(os.path.join(new_dir,"henry.pst"))
    pe = pf.draw(100,use_specsim=True)
    pe.to_binary(os.path.join(new_dir,"prior.jcb"))

    pyemu.os_utils.run("{0} henry.pst".format(exe_path),cwd=new_dir)

def test_process_unc(test_d):
    bd = os.getcwd()
    os.chdir(test_d)
    ret = process_unc()
    os.chdir(bd)
    return ret

def process_unc():
    ucn = flopy.utils.HeadFile("trans.ucn", text="concentration")
    df = pd.read_csv(os.path.join("flow.wel_stress_period_data_historic.txt"),
                    delim_whitespace=True, header=None, names=["l", "r", "c", "flux", "concen"])

    df = df.loc[df.flux<=0, :]
    d = ucn.get_alldata()
    # print(d.shape)
    wel_ucn = {}
    k,j = {},{}
    for l,r,c in zip(df.l,df.r,df.c):
        dd = d[:,l-1,r-1,c-1]
        name = "well_{0:03d}_{1:03d}_{2:03d}".format(l,r,c)
        wel_ucn[name] = dd
        k[name] = l-1
        j[name] = c-1

    df = pd.DataFrame(wel_ucn,index=np.round(ucn.get_times(),5))
    df.loc[:,"time"] = df.index
    cols = df.columns.to_list()
    cols.sort()
    cols.remove("time")
    cols.insert(0,"time")
    df = df.loc[:,cols]
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
    t_d = os.path.join("mou_tests", "henry_template")
    m_d = os.path.join("mou_tests","henry_master_chance_risk_obj")
    if with_master:
        if os.path.exists(m_d):
            shutil.rmtree(m_d)
        shutil.copytree(t_d,m_d)
        pst = pyemu.Pst(os.path.join(m_d,"henry.pst"))
        pst.control_data.noptmax = 100
        pst.pestpp_options["opt_par_stack"] = "prior.jcb"
        pst.pestpp_options["opt_stack_size"] = 20
        pst.pestpp_options["opt_recalc_chance_every"] = 100
        pst.pestpp_options["mou_population_size"] = 100
        pst.pestpp_options["opt_chance_points"] = "all"
        pst.pestpp_options["opt_risk"] = 0.70

        pst.write(os.path.join(m_d,"henry.pst"))
    pyemu.os_utils.start_workers(t_d, exe_path, "henry.pst",
                                  num_workers=15, worker_root="mou_tests",
                                  port=4004)

def plot_results(m_d):
    # plt_d = "henry_results"
    # if os.path.exists(plt_d):
    #     shutil.rmtree(plt_d)
    # os.mkdir(plt_d)
    df_arc = pd.read_csv(os.path.join(m_d,"henry.pareto.archive.summary.csv"))
    pst = pyemu.Pst(os.path.join(m_d,"henry.pst"))
    obj_names = pst.pestpp_options["mou_objectives"].split(',')

    par = pst.parameter_data
    par = par.loc[par.pargp=="dv_pars",:].copy()
    rate_pars = par.loc[par.parnme.apply(lambda x: not x.startswith("ar")),"parnme"]
    print(rate_pars)

    dmn = 60
    dmx = 140

    gens = df_arc.generation.unique()
    gens.sort()
    print(gens)

    cmap = plt.get_cmap("jet")
    norm = Normalize(0.0,35.0)
    with PdfPages(os.path.join(os.path.split(m_d)[-1]+".pdf")) as pdf:
        for gen in [gens[-1]]:
            #df_gen = df_arc.loc[df_arc.generation==gen,:].copy()
            df_gen = pd.read_csv(os.path.join(m_d,"henry.{0}.archive.dv_pop.csv".format(gen)))

            # flip the rate pars
            df_gen.loc[:, rate_pars] *= -1.
            tot = df_gen.loc[:, rate_pars].sum(axis=1)

            # only show solutions with some min amount of pumping
            df_gen = df_gen.loc[tot > 1.0]
            # and only show solutions using concen > some % seawater
            df_gen = df_gen.loc[df_gen.ar_concen >= 17.5]
            df_gen = df_gen.loc[df_gen.ar_rate <= 2.0]

            if df_gen.shape[0] == 0:
                continue
            fig,axes = plt.subplots(4,1,figsize=(15,6))
            ax = axes[0]
            for d,w,r,c,r1,r2,r3 in zip(df_gen.ar_dist,df_gen.ar_width,df_gen.ar_rate,
                               df_gen.ar_concen,df_gen.loc[:,rate_pars[0]],
                                        df_gen.loc[:,rate_pars[1]],
                                        df_gen.loc[:,rate_pars[2]]):
                r = Rectangle((d,0),w,r,facecolor=cmap(c/35.0),edgecolor="none",alpha=0.5)
                ax.add_patch(r)
                mpt = d + (w/2.)
                axes[1].plot([mpt,mpt],[0,r1],color=cmap(c/35.0))
                axes[2].plot([mpt, mpt], [0, r2], color=cmap(c / 35.0))
                axes[3].plot([mpt, mpt], [0, r3], color=cmap(c / 35.0))

            ax.set_xlim(dmn,dmx)
            ax.set_ylim(0,3.5)
            ax.set_xlabel("column")
            ax.set_ylabel("recharge flux rate")
            #cb = plt.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap))
            #cb.set_label("recharge concentration")
            ax.set_title("generation {0}, {1} feasible nondom solutions".format(gen,df_gen.shape[0]))
            for i,ax in enumerate(axes[1:]):
                ax.set_ylim(0,2.0)
                ax.set_xlim(dmn,dmx)
                ax.set_title("pumping well {0}".format(i+1))
                ax.set_ylabel("pumping rate")


            plt.tight_layout()
            #plt.show()
            pdf.savefig()
            plt.close(fig)
            #if gen > 10:
            #    break

def invest():
    m_d = os.path.join("mou_tests","henry_master_chance")
    pst = pyemu.Pst(os.path.join(m_d,"henry.pst"))
    dv = pd.read_csv(os.path.join(m_d,"henry.0.dv_pop.csv"),index_col=0)
    pst.parameter_data.loc[:,"parval1"] = dv.loc[dv.index[0],pst.par_names]
    pst.write_input_files(m_d)
    pyemu.os_utils.run("python forward_run.py",cwd=m_d)

    run_and_plot_results(m_d)

if __name__ == "__main__":
    #shutil.copy2(os.path.join("..", "bin", "win", "pestpp-mou.exe"), os.path.join("..", "bin", "pestpp-mou.exe"))
    shutil.copy2(os.path.join("..","exe","windows","x64","Debug","pestpp-mou.exe"),os.path.join("..","bin","pestpp-mou.exe"))

    #prep_model()
    #run_and_plot_results(os.path.join("mou_tests", "henry_temp"))
    #test_add_artrch("henry_template",write_tpl=False)
    #test_process_unc("henry_temp")
    #setup_pst()
    #run_and_plot_results(os.path.join("mou_tests", "henry_template"))
    start_workers_for_debug(True)
    #plot_pr_real()
    #plot_results(os.path.join("mou_tests","henry_master"))
    #invest()
