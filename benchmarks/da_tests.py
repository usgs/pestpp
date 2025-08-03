import os
import sys
import shutil
import platform
import numpy as np
import pandas as pd
import platform
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
exe_path = os.path.join(bin_path, "pestpp-da" + exe)


noptmax = 4
num_reals = 20
port = 4021


def da_prep_4_mf6_freyberg_seq(sync_state_names=True):
    t_d = os.path.join("mf6_freyberg","template_seq")
    if os.path.exists(t_d):
        shutil.rmtree(t_d)
    shutil.copytree(os.path.join("mf6_freyberg","template"),t_d)
    for f in os.listdir(t_d):
        for tag in ["ies","opt","glm"]:
            if tag in f.lower():
                os.remove(os.path.join(t_d,f))

    #first modify the tdis
    with open(os.path.join(t_d,"freyberg6.tdis"),'w') as f:
        f.write("BEGIN Options\n  TIME_UNITS  days\nEND Options\n")
        f.write("BEGIN Dimensions\n  NPER  1\nEND Dimensions\n")
        f.write("BEGIN PERIODDATA\n31.00000000  1       1.00000000\nEND PERIODDATA\n")
    #make sure it runs
    pyemu.os_utils.run("mf6",cwd=t_d)

    # write a tdis template file - could possibly keep all 25 stress periods to
    # simulate a 2-year-ahead forecast...

    with open(os.path.join(t_d,"freyberg6.tdis.tpl"),'w') as f:
        f.write("ptf  ~\n")
        f.write("BEGIN Options\n  TIME_UNITS  days\nEND Options\n")
        f.write("BEGIN Dimensions\n  NPER  1\nEND Dimensions\n")
        f.write("BEGIN PERIODDATA\n~  perlen  ~  1       1.00000000\nEND PERIODDATA\n")
    new_tpl,new_in = [os.path.join(t_d,"freyberg6.tdis.tpl")],[os.path.join(t_d,"freyberg6.tdis")]
    new_tpl_cycle = [-1]

    # mod the sto to make sp 1 transient - or should this be a pre-processor so that cycle 0 is 
    # ss and the rest are transient?

    # split out the head, sfr and list instruction files into multiple instruction file
    #eventually, want to move to da par and obs cycle tables for heads and gage_1 obs
    lines = open(os.path.join(t_d,"heads.csv.ins"),'r').readlines()[2:]
    new_ins,new_out,new_ins_cycle = [],[],[]
    #print(lines)
    for icycle, line in enumerate(lines):
        ins_name = os.path.join(t_d,"heads_{0}.csv.ins".format(icycle))
        with open(ins_name,'w') as f:
            f.write("pif ~\nl1\n")
            f.write(line)
        new_ins.append(ins_name)
        new_out.append(os.path.join(t_d,"heads.csv"))
        new_ins_cycle.append(icycle)
    remove_ins = ["heads.csv.ins"]

    lines = open(os.path.join(t_d,"sfr.csv.ins"),'r').readlines()[2:]
    #print(lines)
    for icycle, line in enumerate(lines):
        ins_name = os.path.join(t_d,"sfr_{0}.csv.ins".format(icycle))
        with open(ins_name,'w') as f:
            f.write("pif ~\nl1\n")
            f.write(line)
        new_ins.append(ins_name)
        new_out.append(os.path.join(t_d,"sfr.csv"))
        new_ins_cycle.append(icycle)
    remove_ins.append("sfr.csv.ins")

    lines = open(os.path.join(t_d,"freyberg6.lst.ins"),'r').readlines()[1:]
    icycle = 0
    tag_line = lines[0]
    for s in range(0,len(lines),13):
        ins_name = os.path.join(t_d,"freyberg6_{0}.lst.ins".format(icycle))
        with open(os.path.join(t_d,"freyberg6_{0}.lst.ins".format(icycle)),'w') as f:
            f.write("pif ~\n")
            f.write(tag_line)
            for line in lines[s+1:s+13]:
                f.write(line)
        new_ins.append(ins_name)
        new_out.append(os.path.join(t_d,"freyberg6.lst"))
        new_ins_cycle.append(icycle)
        icycle += 1
    remove_ins.append("freyberg6.lst.ins")

    # modify the ic file
    k = 0
    with open(os.path.join(t_d,"freyberg6.ic"),'r') as f:
        while True:
            line = f.readline()
            if line == "":
                break
            if line.lower().strip().startswith("internal"):
                arr_lines = []
                while True:
                    line = f.readline()
                    if line == "":
                        raise Exception
                    if line.lower().strip().startswith("end"):
                        break
                    if line.lower().strip().startswith("internal"):
                        arr = np.array(arr_lines,dtype=float)
                        np.savetxt(os.path.join(t_d,"heads_{0}.dat_in".format(k)),arr,fmt="%15.6E")
                        k += 1
                        arr_lines = []
                    else:
                        arr_lines.append(line.strip().split())

        arr = np.array(arr_lines, dtype=float)
        np.savetxt(os.path.join(t_d, "heads_{0}.dat_in".format(k)), arr, fmt="%15.6E")
    with open(os.path.join(t_d,"freyberg6.ic"),'w') as f:
        f.write("begin griddata\nstrt layered\n")
        for k in range(3):
            f.write("open/close 'heads_{0}.dat_in' FACTOR 1.0\n".format(k))
        f.write("end griddata\n")


    # write a python script to extract final heads and save to files
    with open(os.path.join(t_d,"forward_run.py"),'w') as f:
        f.write("import numpy as np\nimport flopy\nimport pyemu\n")
        f.write("pyemu.os_utils.run('mf6')\n")
        f.write("hds = flopy.utils.HeadFile('freyberg6_freyberg.hds')\n")
        f.write("arr = hds.get_data()\n")
        f.write("for k,a in enumerate(arr):\n")
        f.write("    np.savetxt('heads_'+str(k)+'.dat',a,fmt='%15.6E')\n")

    # dont run it so we can harvest the ic values in the arrays for setting the parval1 values
    pyemu.os_utils.run("python forward_run.py",cwd=t_d)

    # now write ins and tpl file for these
    ic_parvals = {}
    obs_to_par_map = dict()
    for k in range(3):
        fname = os.path.join(t_d,"heads_{0}.dat".format(k))
        assert os.path.exists(fname),fname
        arr = np.loadtxt(fname)
        fname_ins = fname + "_out.ins"
        fname_tpl = fname + "_in.tpl"
        in_arr = np.loadtxt(fname_tpl.replace(".tpl",""))
        ft = open(fname_tpl,'w')
        with open(fname_ins,'w') as f:
            f.write("pif ~\n")
            ft.write("ptf ~\n")
            for i in range(arr.shape[0]):
                f.write("l1 ")
                for j in range(arr.shape[1]):
                    if np.abs(arr[i,j]) > 100 or np.abs(in_arr[i,j]) > 100:
                        f.write(" !dum! ")
                        ft.write(" 40 ")
                    else:
                        oname = "head_{0:02d}_{1:03d}_{2:03d}".format(k,i,j)
                        f.write(" !{0}! ".format(oname))
                        if sync_state_names:
                            ft.write(" ~  {0} ~ ".format(oname))
                            ic_parvals[oname] = in_arr[i, j]
                        else:
                            pname = "p"+oname
                            ft.write(" ~  {0} ~ ".format(pname))
                            obs_to_par_map[oname] = pname
                            ic_parvals[pname] = in_arr[i, j]

                f.write("\n")
                ft.write("\n")
        ft.close()
        new_tpl.append(fname_tpl)
        new_in.append(fname_tpl.replace(".tpl",""))
        new_tpl_cycle.append(-1)
        new_ins.append(fname_ins)
        new_out.append(os.path.join(t_d,"heads_{0}.dat".format(k)))
        new_ins_cycle.append(-1)

        i = pyemu.pst_utils.InstructionFile(fname_ins)
        df = i.read_output_file(fname)
        #print(df)

    # split out the wel and rch tpl files into cycle files
    lines = []
    with open(os.path.join(t_d,"freyberg6.wel.tpl"),'r') as f:
        for i in range(19):
            lines.append(f.readline())
    print(lines)

    for icycle in range(25):
        tpl_file = os.path.join(t_d,"freyberg6.wel_{0}.tpl".format(icycle))
        with open(tpl_file,'w') as f:
            for line in lines:
                new_line = line.replace("_0","_{0}".format(icycle))
                f.write(new_line)

        new_tpl.append(tpl_file)
        new_in.append(os.path.join(t_d,"freyberg6.wel"))
        new_tpl_cycle.append(icycle)
    remove_tpl = ["freyberg6.wel.tpl"]

    lines = []
    with open(os.path.join(t_d, "freyberg6.rch.tpl"), 'r') as f:
        for i in range(11):
            lines.append(f.readline())
    print(lines)

    for icycle in range(25):
        tpl_file = os.path.join(t_d,"freyberg6.rch_{0}.tpl".format(icycle))
        with open(tpl_file,'w') as f:
            for line in lines:
                new_line = line.replace("_0","_{0}".format(icycle))
                f.write(new_line)

        new_tpl.append(tpl_file)
        new_in.append(os.path.join(t_d,"freyberg6.rch"))
        new_tpl_cycle.append(icycle)
    remove_tpl.append('freyberg6.rch.tpl')

    # now for the fun part: modify the pst
    pst = pyemu.Pst(os.path.join(t_d,"freyberg6_run.pst"))

    print(pst.npar_adj,pst.nnz_obs)

    # swap out obs info
    dropped_dfs = []
    for ins in remove_ins:
        dropped_dfs.append(pst.drop_observations(os.path.join(t_d, ins), '.'))
    for insf, outf, cy in zip(new_ins, new_out, new_ins_cycle):
        df = pst.add_observations(insf,outf, pst_path=".")
        pst.observation_data.loc[df.obsnme, "cycle"] = cy
        pst.model_output_data.loc[os.path.join(".",os.path.split(insf)[1]),"cycle"] = cy
    pst.observation_data.loc[:,"weight"] = 0.0
    for df in dropped_dfs:
        for c in ["obsval","weight"]:
            pst.observation_data.loc[df.obsnme, c] = df.loc[:, c]

    # swap out par info
    dropped_dfs = []
    for tpl in remove_tpl:
        dropped_dfs.append(pst.drop_parameters(os.path.join(t_d,tpl),'.'))
    pst.parameter_data.loc[:,"cycle"] = -1
    pst.model_input_data.loc[:,"cycle"] = -1
    for tplf, inf, cy in zip(new_tpl,new_in,new_tpl_cycle):
        df = pst.add_parameters(tplf,inf,pst_path=".")
        pst.parameter_data.loc[df.parnme,"cycle"] = cy
        pst.model_input_data.loc[os.path.join(".",os.path.split(tplf)[1]),"cycle"] = cy
    for df in dropped_dfs:
        for c in ["parval1","parubnd","parlbnd","pargp"]:
            pst.parameter_data.loc[df.parnme,c] = df.loc[:,c]

    #set the state parameter info
    for p,v in ic_parvals.items():
        pst.parameter_data.loc[p,"parval1"] = v
        pst.parameter_data.loc[p, "parlbnd"] = v * 0.5
        pst.parameter_data.loc[p, "parubnd"] = v * 1.5
        pst.parameter_data.loc[p,"pargp"] = "head_state"
        pst.parameter_data.loc[p,"partrans"] = "none"

    pst.control_data.noptmax = 2
    pst.model_command = "python forward_run.py"
    pst.pestpp_options.pop("ies_par_en")
    pst.parameter_data.loc["perlen","partrans"] = "fixed"
    pst.pestpp_options["ies_num_reals"] = 5
    pst.pestpp_options["da_num_reals"] = 5
    if not sync_state_names:
        pst.observation_data.loc[:,"state_par_link"] = np.nan
        obs = pst.observation_data
        obs.loc[:,"state_par_link"] = obs.obsnme.apply(lambda x: obs_to_par_map.get(x,np.nan))
    pst.write(os.path.join(t_d,"freyberg6_run_da1.pst"),version=2)
    return pst

def da_mf6_freyberg_test_1():
    test_d = "mf6_freyberg"
    t_d = os.path.join(test_d,"template_seq")
    
    pst = da_prep_4_mf6_freyberg_seq()
    
    print(exe_path.replace("ies","da"))
    pst.control_data.noptmax = 0
    pst.pestpp_options["ies_verbose_level"] = 4
    pst.pestpp_options["ies_no_noise"] = True
    pst.pestpp_options["da_use_mda"] = False
    pst.pestpp_options["da_lambda_mults"] = 1.0
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.write(os.path.join(t_d, "freyberg6_run_da1.pst"), version=2)
    #pyemu.os_utils.run("{0} freyberg6_run_da1.pst".format(exe_path.replace("ies","da")),cwd=t_d)

    pst.pestpp_options["ies_num_reals"] = 15
    pst.control_data.noptmax = 2
    pst.write(os.path.join(t_d, "freyberg6_run_da1.pst"), version=2)
    pyemu.os_utils.start_workers(t_d,exe_path.replace("ies","da"),"freyberg6_run_da1.pst",
                                 num_workers=15,worker_root=test_d,port=port,
                                 master_dir=os.path.join(test_d,"master_da_1"),verbose=True)
    
def da_prep_4_mf6_freyberg_seq_tbl():
    test_d = "mf6_freyberg"
    t_d = os.path.join(test_d,"template_seq")
    if not os.path.exists(t_d):
        da_prep_4_mf6_freyberg_seq()
    pst = pyemu.Pst(os.path.join(t_d,"freyberg6_run_da1.pst"))
    pdf = pd.DataFrame({"perlen":31},index=np.arange(25))
    pdf.T.to_csv(os.path.join(t_d,"par_cycle_tbl.csv"))
    pst.pestpp_options["da_parameter_cycle_table"] = "par_cycle_tbl.csv"

    # mod sfr_out
    sfr_ins_file = os.path.join(t_d, "sfr.csv.ins")
    with open(sfr_ins_file, 'w') as f:
        f.write("pif ~\n")
        f.write("l1\n")
        f.write("l1 ~,~ !headwater!  ~,~ !tailwater!  ~,~ !gage_1!\n")

    # and lst budget
    lines = open(os.path.join(t_d, "freyberg6_0.lst.ins"), 'r').readlines()
    lst_ins_file = os.path.join(t_d, "freyberg6.lst.ins")
    with open(lst_ins_file, 'w') as f:
        for line in lines:
            f.write(line.replace("_20151231",""))

    obs = pst.observation_data.copy()
    tr_obs = obs.loc[obs.obsnme.str.startswith("trgw"),:].copy()
    tr_obs.loc[tr_obs.obsnme,"datetime"] = pd.to_datetime(tr_obs.obsnme.apply(lambda x: x.split('_')[-1]))
    tr_obs.loc[tr_obs.obsnme,"k"] = tr_obs.obsnme.apply(lambda x: int(x.split('_')[1]))
    tr_obs.loc[tr_obs.obsnme, "i"] = tr_obs.obsnme.apply(lambda x: int(x.split('_')[2]))
    tr_obs.loc[tr_obs.obsnme, "j"] = tr_obs.obsnme.apply(lambda x: int(x.split('_')[3]))
    tr_obs.loc[tr_obs.obsnme,"obgnme"] = tr_obs.obsnme.apply(lambda x: "_".join(x.split("_")[:-1]))

    head_obs = obs.loc[obs.obsnme.str.startswith("head_"),:].copy()
    head_obs.loc[head_obs.obsnme, "k"] = head_obs.obsnme.apply(lambda x: int(x.split('_')[1]))
    head_obs.loc[head_obs.obsnme, "i"] = head_obs.obsnme.apply(lambda x: int(x.split('_')[2]))
    head_obs.loc[head_obs.obsnme, "j"] = head_obs.obsnme.apply(lambda x: int(x.split('_')[3]))

    #print(pst.nobs)
    for ins_file in pst.model_output_data.pest_file.copy():
        if "heads_" in ins_file and ins_file.endswith("dat_out.ins"):
            continue
        pst.drop_observations(os.path.join(t_d,ins_file),pst_path=".")
    for ins_file in [sfr_ins_file,lst_ins_file]:
        pst.add_observations(ins_file,pst_path=".")

    # work out which heads are obs site
    obs_heads = {}
    odf_names = []
    pst.observation_data.loc[:,"org_obgnme"] = np.nan
    pst.observation_data.loc[:, "weight"] = 0.0
    for og in tr_obs.obgnme.unique():
        site_obs = tr_obs.loc[tr_obs.obgnme==og,:]
        site_obs.sort_values(by="datetime",inplace=True)
        head_name = "head_{0:02d}_{1:03d}_{2:03d}".format(int(site_obs.k[0]),int(site_obs.i[0]),int(site_obs.j[0]))
        for i,oname in enumerate(site_obs.obsnme):
            obs_heads[oname] = (head_name,i)
        # assign max weight in the control file since some have zero weight and
        # we are covering weights in the weight table
        #print(og,site_obs.weight.max())
        pst.observation_data.loc[head_name,"weight"] = site_obs.weight.max()
        pst.observation_data.loc[head_name,"org_obgnme"] = og
        odf_names.append(head_name)

    odf_names.append("gage_1")

    odf = pd.DataFrame(columns=odf_names,index=np.arange(25))
    wdf = pd.DataFrame(columns=odf_names,index=np.arange(25))
    for tr_name,(head_name,icycle) in obs_heads.items():
        odf.loc[icycle,head_name] = obs.loc[tr_name,"obsval"]
        wdf.loc[icycle, head_name] = obs.loc[tr_name, "weight"]

    g_obs = obs.loc[obs.obsnme.str.startswith("gage_1"),:].copy()
    #give these obs the max weight since some have zero weight
    pst.observation_data.loc["gage_1", "weight"] = g_obs.weight.max()
    g_obs.sort_index(inplace=True)
    for i,name in enumerate(g_obs.obsnme):
        odf.loc[i,"gage_1"] = g_obs.loc[name,"obsval"]
        wdf.loc[i, "gage_1"] = g_obs.loc[name, "weight"]

    # now drop any entries that have zero weight across all cycles
    print(pst.nnz_obs_names)
    odf = odf.loc[:,pst.nnz_obs_names]
    wdf = wdf.loc[:, pst.nnz_obs_names]

    odf.T.to_csv(os.path.join(t_d,"obs_cycle_tbl.csv"))
    pst.pestpp_options["da_observation_cycle_table"] = "obs_cycle_tbl.csv"
    wdf.T.to_csv(os.path.join(t_d, "weight_cycle_tbl.csv"))
    pst.pestpp_options["da_weight_cycle_table"] = "weight_cycle_tbl.csv"

    pst.observation_data.loc[:,"cycle"] = -1
    pst.model_output_data.loc[:,"cycle"] = -1

    pst.write(os.path.join(t_d,"freyberg6_run_da2.pst"),version=2)
    return pst

def da_mf6_freyberg_test_2():
    pst = da_prep_4_mf6_freyberg_seq_tbl()
    test_d = "mf6_freyberg"
    t_d = os.path.join(test_d, "template_seq")
    print(exe_path.replace("ies", "da"))
    pst.control_data.noptmax = 0
    pst.pestpp_options["ies_verbose_level"] = 4
    pst.pestpp_options["ies_no_noise"] = True
    pst.pestpp_options["da_use_mda"] = False
    pst.pestpp_options["da_lambda_mults"] = 1.0
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.write(os.path.join(t_d, "freyberg6_run_da2.pst"), version=2)
    pyemu.os_utils.run("{0} freyberg6_run_da2.pst".format(exe_path.replace("ies","da")),cwd=t_d)

    pst.control_data.noptmax = -1
    pst.pestpp_options["ies_num_reals"] = 15
    pst.write(os.path.join(t_d, "freyberg6_run_da2.pst"), version=2)
    pyemu.os_utils.start_workers(t_d, exe_path.replace("ies", "da"), "freyberg6_run_da2.pst",
                                num_workers=15, worker_root=test_d, port=port,
                                master_dir=os.path.join(test_d, "master_da_2"), verbose=True)


def invest():
    t_d = os.path.join("mf6_freyberg","template_seq")
    pst1 = pyemu.Pst(os.path.join(t_d,"freyberg6_run_da1.pst"))
    pst2 = pyemu.Pst(os.path.join(t_d,"freyberg6_run_da2.pst"))
    print(pst1.observation_data.loc[pst1.nnz_obs_names,["obsval","weight"]])
    odf = pd.read_csv(os.path.join(t_d,pst2.pestpp_options["da_observation_cycle_table"]))
    wdf = pd.read_csv(os.path.join(t_d, pst2.pestpp_options["da_weight_cycle_table"]))
    print(odf)
    print(wdf)

def da_mf6_freyberg_smoother_test():
    model_d = "mf6_freyberg"
    local=True
    if "linux" in platform.platform().lower() and "10par" in model_d:
        #print("travis_prep")
        #prep_for_travis(model_d)
        local=False
    
    t_d = os.path.join(model_d,"template")

    m_d = os.path.join(model_d,"master_da_smoother")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pst = pyemu.Pst(os.path.join(t_d,"freyberg6_run_ies.pst"))
    pst.control_data.noptmax = 0
    pst.pestpp_options["da_autoadaloc"] = False
    pst.pestpp_options.pop("da_localizer",None)
    pst.pestpp_options["da_num_reals"] = 15
    pst.pestpp_options["da_use_mda"] = False
    pst.pestpp_options["da_lambda_mults"] = 1.0
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    new_dict = {}
    for key,val in pst.pestpp_options.items():
        new_dict[key.replace("ies","da")] = val
    pst.pestpp_options = new_dict

    pst.write(os.path.join(t_d,"freyberg6_run_da.pst"))
    pyemu.os_utils.run("{0} freyberg6_run_da.pst".format(exe_path.replace("-ies","-da")),cwd=t_d)
    
    pst.control_data.noptmax = 2
    
    pst.write(os.path.join(t_d,"freyberg6_run_da.pst"))
    pyemu.os_utils.start_workers(t_d, exe_path.replace("-ies","-da"), "freyberg6_run_da.pst", num_workers=15,
                                master_dir=m_d,worker_root=model_d,port=port)

    
    oe_file = os.path.join(m_d,"freyberg6_run_da.global.0.oe.csv")
    assert os.path.exists(oe_file),oe_file
    pe_file = os.path.join(m_d,"freyberg6_run_da.global.prior.pe.csv")
    assert os.path.exists(pe_file),pe_file
    phi_file = os.path.join(m_d,"freyberg6_run_da.global.phi.actual.csv")
    assert os.path.exists(phi_file),phi_file


def build_lorenz_pst():
    t_d = os.path.join("lorenz96", "template_seq")
    if os.path.exists(t_d):
        shutil.rmtree(t_d)
    shutil.copytree(os.path.join("lorenz96", "template"), t_d)
    os.chdir(os.path.join(t_d))

    # run forward check
    pyemu.os_utils.run("python rk4_lorenz96.py")

    out_csv = "l96_out.csv"
    odf = pyemu.pst_utils.csv_to_ins_file(os.path.join(out_csv))# only_cols=cols)
    ins_file = out_csv + ".ins"

    in_csv = "l96_in.csv"
    F = pd.read_csv(os.path.join(in_csv))['F'][0]
    tpl_file = in_csv + ".tpl"

    pst = pyemu.helpers.pst_from_io_files(tpl_file, in_csv, ins_file, out_csv)

    obs = pst.observation_data
    obs.loc[odf.index, "obsval"] = odf.obsval

    par = pst.parameter_data
    par.loc["forcing", "parval1"] = F
    par.loc["forcing", "parlbnd"] = 0.0
    par.loc["forcing", "parubnd"] = 9.0

    pst.control_data.noptmax = 1
    pst.model_command = "python forward_run.py"

    with open(os.path.join("forward_run.py"), 'w') as f:
        f.write("import rk4_lorenz96\n")
        f.write("rk4_lorenz96.L96()\n")

    pst.write(pst.filename)

    #do pestpp test run here

    os.chdir(os.path.join("..", ".."))


def da_build_mf6_freyberg_seq_localizer_tbl():
    import flopy
    t_d = os.path.join("mf6_freyberg","template_seq")
    assert os.path.exists(t_d)
    sim = flopy.mf6.MFSimulation.load(sim_ws=t_d)
    m = sim.get_model("freyberg6")
    mg = m.modelgrid
    xcel = mg.xcellcenters
    ycel = mg.ycellcenters
    pst = pyemu.Pst(os.path.join(t_d,"freyberg6_run_da2.pst"))

    obs = pst.observation_data
    nzobs = obs.loc[pst.nnz_obs_names,:]
    snzobs = nzobs.loc[nzobs.obsnme.str.contains("head"),:]
    snzobs.loc[:, "j"] = snzobs.obsnme.apply(lambda x: int(x.split("_")[-1]))
    snzobs.loc[:, "i"] = snzobs.obsnme.apply(lambda x: int(x.split("_")[-2]))
    snzobs.loc[:, "k"] = snzobs.obsnme.apply(lambda x: int(x.split("_")[-3]))
    snzobs.loc[:, "x"] = snzobs.apply(lambda x: xcel[x.i, x.j], axis=1)
    snzobs.loc[:, "y"] = snzobs.apply(lambda x: ycel[x.i, x.j], axis=1)

    par = pst.parameter_data
    spar = par.loc[par.pargp.apply(lambda x: "wel" not in x and "rch" not in x and "pargp" not in x),:]
    spar.loc[:,"j"] = spar.parnme.apply(lambda x: int(x.split("_")[-1]))
    spar.loc[:, "i"] = spar.parnme.apply(lambda x: int(x.split("_")[-2]))
    spar.loc[:, "k"] = spar.parnme.apply(lambda x: int(x.split("_")[-3]))
    spar.loc[:,"x"] = spar.apply(lambda x: xcel[x.i,x.j],axis=1)
    spar.loc[:, "y"] = spar.apply(lambda x: ycel[x.i, x.j], axis=1)
    spar_grps = spar.pargp.unique()

    # these are parameter groups that will not be localized
    other_grps = [g for g in pst.adj_par_groups if g not in spar_grps]
    loc_cols = other_grps
    loc_cols.extend(spar.parnme.tolist())
    d_thres = 2500 # 20 model cells
    df = pd.DataFrame(columns=loc_cols,index=pst.nnz_obs_names)
    df.loc[:,:] = 1.0
    # for each obs, find pars that are within d_thres distance
    for oname in snzobs.obsnme:
        ox,oy = snzobs.loc[oname,"x"], snzobs.loc[oname,"y"]
        # calc squared distance from this obs to all spatial pars
        d = spar.apply(lambda x: (x.x - ox)**2 + (x.y-oy)**2,axis=1)
        # ident pars that are too far
        d = d.loc[d>d_thres**2]
        # set pars that are too far to 0.0 in the localizer
        df.loc[oname,d.index] = 0.0
    print(df.sum().sort_values())
    #print(df.shape)
    pyemu.Matrix.from_dataframe(df).to_coo(os.path.join(t_d,"seq_tbl_loc.jcb"))


def da_build_mf6_freyberg_seq_localizer():
    import flopy
    t_d = os.path.join("mf6_freyberg","template_seq")
    assert os.path.exists(t_d)
    sim = flopy.mf6.MFSimulation.load(sim_ws=t_d)
    m = sim.get_model("freyberg6")
    mg = m.modelgrid
    xcel = mg.xcellcenters
    ycel = mg.ycellcenters
    pst = pyemu.Pst(os.path.join(t_d,"freyberg6_run_da1.pst"))

    obs = pst.observation_data
    nzobs = obs.loc[pst.nnz_obs_names,:]
    snzobs = nzobs.loc[nzobs.obsnme.str.contains("trgw"),:]
    snzobs.loc[:, "j"] = snzobs.obsnme.apply(lambda x: int(x.split("_")[-2]))
    snzobs.loc[:, "i"] = snzobs.obsnme.apply(lambda x: int(x.split("_")[-3]))
    snzobs.loc[:, "k"] = snzobs.obsnme.apply(lambda x: int(x.split("_")[-4]))
    snzobs.loc[:, "x"] = snzobs.apply(lambda x: xcel[x.i, x.j], axis=1)
    snzobs.loc[:, "y"] = snzobs.apply(lambda x: ycel[x.i, x.j], axis=1)

    par = pst.parameter_data
    spar = par.loc[par.pargp.apply(lambda x: "wel" not in x and "rch" not in x and "pargp" not in x),:]
    spar.loc[:,"j"] = spar.parnme.apply(lambda x: int(x.split("_")[-1]))
    spar.loc[:, "i"] = spar.parnme.apply(lambda x: int(x.split("_")[-2]))
    spar.loc[:, "k"] = spar.parnme.apply(lambda x: int(x.split("_")[-3]))
    spar.loc[:,"x"] = spar.apply(lambda x: xcel[x.i,x.j],axis=1)
    spar.loc[:, "y"] = spar.apply(lambda x: ycel[x.i, x.j], axis=1)
    spar_grps = spar.pargp.unique()

    # these are parameter groups that will not be localized
    other_grps = [g for g in pst.adj_par_groups if g not in spar_grps]
    loc_cols = other_grps
    loc_cols.extend(spar.parnme.tolist())
    d_thres = 2500 # 20 model cells
    df = pd.DataFrame(columns=loc_cols,index=pst.nnz_obs_names)
    df.loc[:,:] = 1.0
    # for each obs, find pars that are within d_thres distance
    for oname in snzobs.obsnme:
        ox,oy = snzobs.loc[oname,"x"], snzobs.loc[oname,"y"]
        # calc squared distance from this obs to all spatial pars
        d = spar.apply(lambda x: (x.x - ox)**2 + (x.y-oy)**2,axis=1)
        # ident pars that are too far
        d = d.loc[d>d_thres**2]
        # set pars that are too far to 0.0 in the localizer
        df.loc[oname,d.index] = 0.0
    print(df.sum().sort_values())
    #print(df.shape)
    pyemu.Matrix.from_dataframe(df).to_coo(os.path.join(t_d,"seq_loc.jcb"))

def da_mf6_freyberg_test_3():
    pst = da_prep_4_mf6_freyberg_seq(sync_state_names=True)
    test_d = "mf6_freyberg"
    t_d = os.path.join(test_d, "template_seq")
    print(exe_path.replace("ies", "da"))
    pst.control_data.noptmax = 0
    pst.pestpp_options["ies_verbose_level"] = 4
    pst.pestpp_options["ies_no_noise"] = True
    pst.pestpp_options["da_use_mda"] = False
    pst.pestpp_options["da_lambda_mults"] = 1.0
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.write(os.path.join(t_d, "freyberg6_run_da2.pst"), version=2)
    pyemu.os_utils.run("{0} freyberg6_run_da2.pst".format(exe_path.replace("ies","da")),cwd=t_d)

    pst.control_data.noptmax = -1
    pst.pestpp_options["ies_num_reals"] = 15
    pst.write(os.path.join(t_d, "freyberg6_run_da2.pst"), version=2)
    pyemu.os_utils.start_workers(t_d, exe_path.replace("ies", "da"), "freyberg6_run_da2.pst",
                                num_workers=15, worker_root=test_d, port=port,
                                master_dir=os.path.join(test_d, "master_da_2"), verbose=True)


def seq_10par_xsec_state_est_test():
    #todo: add localization to this test!
    test_d = "10par_xsec"
    t_d = os.path.join(test_d, "template")
    pst = pyemu.Pst(os.path.join(t_d,"pest.pst"))
    par = pst.parameter_data
    par.loc[:,"cycle"] = -1
    par.loc[:,"partrans"] = "fixed"
    par.loc[par.parnme.str.contains("strt"),"partrans"] = "log"
    strt_pars = par.loc[par.pargp=="strt","parnme"].tolist()
    obs = pst.observation_data
    obs.loc[obs.obsnme.str.startswith("h01"),"weight"] = 1.0
    obs.loc[:,"state_par_link"] = ""
    obs.loc[obs.obgnme=="head1","state_par_link"] = strt_pars
    obs.loc[:,"cycle"] = -1
    pst.control_data.noptmax = 1

    pst.model_input_data.loc[:,"cycle"] = -1
    pst.model_output_data.loc[:,"cycle"] = -1

    def get_loc(pst):

        loc = pd.DataFrame(index=pst.nnz_obs_names,columns=pst.adj_par_names)
        loc.loc[:,:] = 0.0
        ocells = loc.index.map(lambda x: int(x.split('_')[1]))
        for pname in pst.adj_par_names:
            cstr = pname.split('_')[1]
            cint = int(cstr)
            # strt states can comm with all obs locs
            if "strt" in pname:
                dist = (ocells - cint).map(np.abs)
                loc_vals = 1 / (dist + 0.01)
                loc_vals = loc_vals.values
                loc_vals[loc_vals>1.0] = 1.0
                loc.loc[:,pname] = loc_vals
            # static pars can only comm with obs in the same cell
            else:

                oname = [o for o in pst.nnz_obs_names if cstr in o.split('_')[1] == cstr][0]
                loc.loc[oname, pname] = 1.0
        return loc

    loc = get_loc(pst)
    pyemu.Matrix.from_dataframe(loc).to_ascii(os.path.join(t_d,"loc.mat"))

    cycles = np.arange(0,5)
    odf = pd.DataFrame(index=cycles,columns=pst.nnz_obs_names)
    odf.loc[:,:] = obs.loc[pst.nnz_obs_names,"obsval"].values
    odf.T.to_csv(os.path.join(t_d,"obs_cycle_tbl.csv"))
    wdf = pd.DataFrame(index=cycles,columns=pst.nnz_obs_names)
    wdf.loc[:,:] = 0.0
    wdf.iloc[1,[3,5]] = 0.1
    wdf.iloc[3,:] = 0.1
    wdf.T.to_csv(os.path.join(t_d,"weight_cycle_tbl.csv"))

    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["ies_lambda_mults"] = 1.0
    pst.pestpp_options["da_use_mda"] = True
    pst.pestpp_options["da_observation_cycle_table"] = "obs_cycle_tbl.csv"
    pst.pestpp_options["da_weight_cycle_table"] = "weight_cycle_tbl.csv"
    pst.pestpp_options["ies_num_reals"] = 20

    pst.write(os.path.join(t_d,"pest_seq.pst"),version=2)
    #pyemu.os_utils.run("{0} pest_seq.pst".format(exe_path),cwd=t_d)
    pyemu.os_utils.start_workers(t_d, exe_path.replace("ies", "da"), "pest_seq.pst",
                                 num_workers=pst.pestpp_options["ies_num_reals"], worker_root=test_d, port=port,
                                 master_dir=os.path.join(test_d, "master_se_mda"), verbose=True)

    pst.pestpp_options["ies_localizer"] = "loc.mat"
    pst.write(os.path.join(t_d, "pest_seq.pst"), version=2)
    # pyemu.os_utils.run("{0} pest_seq.pst".format(exe_path),cwd=t_d)
    pyemu.os_utils.start_workers(t_d, exe_path.replace("ies", "da"), "pest_seq.pst",
                                 num_workers=pst.pestpp_options["ies_num_reals"], worker_root=test_d, port=port,
                                 master_dir=os.path.join(test_d, "master_se_mda_local"), verbose=True)

    pst.pestpp_options["ies_autoadaloc"] = True
    pst.write(os.path.join(t_d, "pest_seq.pst"), version=2)
    pyemu.os_utils.start_workers(t_d, exe_path.replace("ies", "da"), "pest_seq.pst",
                                 num_workers=pst.pestpp_options["ies_num_reals"], worker_root=test_d, port=port,
                                 master_dir=os.path.join(test_d, "master_se_mda_local_aad"), verbose=True)

    pst.pestpp_options["ies_loc_type"] = "cov"
    pst.write(os.path.join(t_d, "pest_seq.pst"), version=2)
    #pyemu.os_utils.start_workers(t_d, exe_path.replace("ies", "da"), "pest_seq.pst",
    #                             num_workers=pst.pestpp_options["ies_num_reals"], worker_root=test_d, port=port,
    #                             master_dir=os.path.join(test_d, "master_se_mda_cov_aad"), verbose=True)


    pst.pestpp_options["ies_autoadaloc"] = False
    pst.pestpp_options["ies_loc_type"] = "local"
    pst.pestpp_options["da_use_mda"] = False
    pst.write(os.path.join(t_d, "pest_seq.pst"), version=2)
    pyemu.os_utils.start_workers(t_d, exe_path.replace("ies", "da"), "pest_seq.pst",
                                 num_workers=pst.pestpp_options["ies_num_reals"], worker_root=test_d, port=port,
                                 master_dir=os.path.join(test_d, "master_se_glm_local"), verbose=True)

    pst.pestpp_options["ies_autoadaloc"] = True
    pst.write(os.path.join(t_d, "pest_seq.pst"), version=2)
    pyemu.os_utils.start_workers(t_d, exe_path.replace("ies", "da"), "pest_seq.pst",
                                 num_workers=pst.pestpp_options["ies_num_reals"], worker_root=test_d, port=port,
                                 master_dir=os.path.join(test_d, "master_se_glm_local_aad"), verbose=True)

    pst.pestpp_options["ies_loc_type"] = "cov"
    pst.write(os.path.join(t_d, "pest_seq.pst"), version=2)
    pyemu.os_utils.start_workers(t_d, exe_path.replace("ies", "da"), "pest_seq.pst",
                                 num_workers=pst.pestpp_options["ies_num_reals"], worker_root=test_d, port=port,
                                 master_dir=os.path.join(test_d, "master_se_glm_cov_aad"), verbose=True)

    pst.parameter_data.loc[:, "partrans"] = "log"
    loc = get_loc(pst)
    pyemu.Matrix.from_dataframe(loc).to_ascii(os.path.join(t_d, "loc.mat"))
    pst.pestpp_options["ies_autoadaloc"] = False
    pst.pestpp_options["ies_loc_type"] = "local"
    pst.pestpp_options["da_use_mda"] = True
    pst.pestpp_options["lambda_scale_fac"] = [0.5,1.0]
    pst.pestpp_options["ies_lambda_mults"] = [0.1,1.0,10.0]
    pst.control_data.noptmax = 3
    pst.pestpp_options.pop("ies_localizer")
    pst.write(os.path.join(t_d, "pest_seq.pst"), version=2)
    pyemu.os_utils.start_workers(t_d, exe_path.replace("ies", "da"), "pest_seq.pst",
                                 num_workers=pst.pestpp_options["ies_num_reals"], worker_root=test_d, port=port,
                                 master_dir=os.path.join(test_d, "master_se_glm"), verbose=True)

    pst.pestpp_options["ies_localizer"] = "loc.mat"
    pst.write(os.path.join(t_d, "pest_seq.pst"), version=2)
    pyemu.os_utils.start_workers(t_d, exe_path.replace("ies", "da"), "pest_seq.pst",
                                 num_workers=pst.pestpp_options["ies_num_reals"], worker_root=test_d, port=port,
                                 master_dir=os.path.join(test_d, "master_da_mda_local"), verbose=True)

    pst.pestpp_options["ies_autoadaloc"] = True
    pst.write(os.path.join(t_d, "pest_seq.pst"), version=2)
    pyemu.os_utils.start_workers(t_d, exe_path.replace("ies", "da"), "pest_seq.pst",
                                 num_workers=pst.pestpp_options["ies_num_reals"], worker_root=test_d, port=port,
                                 master_dir=os.path.join(test_d, "master_da_mda_local_aad"), verbose=True)

    pst.pestpp_options["ies_loc_type"] = "cov"
    pst.write(os.path.join(t_d, "pest_seq.pst"), version=2)
    #pyemu.os_utils.start_workers(t_d, exe_path.replace("ies", "da"), "pest_seq.pst",
    #                             num_workers=pst.pestpp_options["ies_num_reals"], worker_root=test_d, port=port,
    #                             master_dir=os.path.join(test_d, "master_da_mda_cov_aad"), verbose=True)

    pst.pestpp_options["ies_autoadaloc"] = False
    pst.pestpp_options["ies_loc_type"] = "local"
    pst.pestpp_options["da_use_mda"] = False
    pst.write(os.path.join(t_d, "pest_seq.pst"), version=2)
    pyemu.os_utils.start_workers(t_d, exe_path.replace("ies", "da"), "pest_seq.pst",
                                 num_workers=pst.pestpp_options["ies_num_reals"], worker_root=test_d, port=port,
                                 master_dir=os.path.join(test_d, "master_da_glm_local"), verbose=True)

    pst.pestpp_options["ies_autoadaloc"] = True
    pst.write(os.path.join(t_d, "pest_seq.pst"), version=2)
    pyemu.os_utils.start_workers(t_d, exe_path.replace("ies", "da"), "pest_seq.pst",
                                 num_workers=pst.pestpp_options["ies_num_reals"], worker_root=test_d, port=port,
                                 master_dir=os.path.join(test_d, "master_da_glm_local_aad"), verbose=True)

    pst.pestpp_options["ies_loc_type"] = "cov"
    pst.write(os.path.join(t_d, "pest_seq.pst"), version=2)
    pyemu.os_utils.start_workers(t_d, exe_path.replace("ies", "da"), "pest_seq.pst",
                                 num_workers=pst.pestpp_options["ies_num_reals"], worker_root=test_d, port=port,
                                 master_dir=os.path.join(test_d, "master_da_glm_cov_aad"), verbose=True)


def seq_10par_xsec_fixed_test():
    
    test_d = "10par_xsec"
    t_d = os.path.join(test_d, "template")
    pst = pyemu.Pst(os.path.join(t_d,"pest.pst"))
    par = pst.parameter_data
    par.loc[:,"cycle"] = -1
    
    par.loc[par.parnme.str.contains("strt"),"partrans"] = "log"
    par.loc[["cnhd_01","strt_02","strt_03"],"partrans"] = "fixed"
    strt_pars = par.loc[par.pargp=="strt","parnme"].tolist()
    obs = pst.observation_data
    obs.loc[obs.obsnme.str.startswith("h01"),"weight"] = 1.0
    obs.loc[:,"state_par_link"] = ""
    obs.loc[obs.obgnme=="head1","state_par_link"] = strt_pars
    obs.loc[:,"cycle"] = -1
    
    pst.control_data.noptmax = 1

    pst.model_input_data.loc[:,"cycle"] = -1
    pst.model_output_data.loc[:,"cycle"] = -1

    def get_loc(pst):

        loc = pd.DataFrame(index=pst.nnz_obs_names,columns=pst.adj_par_names)
        loc.loc[:,:] = 0.0
        ocells = loc.index.map(lambda x: int(x.split('_')[1]))
        for pname in pst.adj_par_names:
            cstr = pname.split('_')[1]
            cint = int(cstr)
            # strt states can comm with all obs locs
            if "strt" in pname:
                dist = (ocells - cint).map(np.abs)
                loc_vals = 1 / (dist + 0.01)
                loc_vals = loc_vals.values
                loc_vals[loc_vals>1.0] = 1.0
                loc.loc[:,pname] = loc_vals
            # static pars can only comm with obs in the same cell
            else:

                oname = [o for o in pst.nnz_obs_names if cstr in o.split('_')[1] == cstr][0]
                loc.loc[oname, pname] = 1.0
        return loc

    loc = get_loc(pst)
    pyemu.Matrix.from_dataframe(loc).to_ascii(os.path.join(t_d,"loc.mat"))

    mx_cycle = 5
    cycles = np.arange(0,mx_cycle)
    odf = pd.DataFrame(index=cycles,columns=pst.nnz_obs_names)
    odf.loc[:,:] = obs.loc[pst.nnz_obs_names,"obsval"].values
    odf.T.to_csv(os.path.join(t_d,"obs_cycle_tbl.csv"))
    wdf = pd.DataFrame(index=cycles,columns=pst.nnz_obs_names)
    wdf.loc[:,:] = 0.0
    wdf.iloc[1,[3,5]] = 1.0
    wdf.iloc[3,:] = 1.0
    wdf.T.to_csv(os.path.join(t_d,"weight_cycle_tbl.csv"))

    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["da_lambda_mults"] = 1.0
    pst.pestpp_options["da_use_mda"] = True
    pst.pestpp_options["da_observation_cycle_table"] = "obs_cycle_tbl.csv"
    pst.pestpp_options["da_weight_cycle_table"] = "weight_cycle_tbl.csv"
    pst.pestpp_options["da_num_reals"] = 10
    pst.pestpp_options["da_localizer"] = "loc.mat"

    pst.write(os.path.join(t_d,"pest_seq.pst"),version=2)
    m_d = os.path.join(test_d, "master_da_fixed")
    pyemu.os_utils.start_workers(t_d, exe_path.replace("ies", "da"), "pest_seq.pst",
                                num_workers=pst.pestpp_options["da_num_reals"], worker_root=test_d, port=port,
                                master_dir=m_d, verbose=True)

    pr_pe = pd.read_csv(os.path.join(m_d,"pest_seq.global.prior.pe.csv"),index_col=0)
    print("testing")
    for cycle in cycles[1:]:
        oe = pd.read_csv(os.path.join(m_d,"pest_seq.global.{0}.oe.csv".format(cycle-1)),index_col=0)
        pe = pd.read_csv(os.path.join(m_d,"pest_seq.{0}.0.par.csv".format(cycle)),index_col=0)
        print(os.path.join(m_d,"pest_seq.global.{0}.oe.csv".format(cycle-1)))
        print(os.path.join(m_d, "pest_seq.{0}.0.par.csv".format(cycle)))
        print(oe)
        print(pe)
        #check that that adj dyn states are being updated and used...
        d = np.abs(oe.loc[:,"h01_02"].values - pe.loc[:,"strt_02"].values)
        print(d)
        assert d.max() < 1.0e-6,d.max()
        d = np.abs(oe.loc[:, "h01_03"].values - pe.loc[:, "strt_03"].values)
        print(d)
        assert d.max() < 1.0e-6, d.max()
        # check that the fixed par (non state) is not being adjusted
        d =  np.abs(pr_pe.loc[:,"cnhd_01"].values - pe.loc[:, "cnhd_01"].values)
        print(d)
        assert d.max() < 1.0e-6, d.max()

    # check that the static pars are being localized before cycle 3 (as per weight table)
    # only static pars in cells 4 and 6 should be allowed to change in cycle 1
    pe_cycle2 = pd.read_csv(os.path.join(m_d,"pest_seq.global.2.pe.csv"),index_col=0)
    fpar_names = par.loc[par.parnme.apply(lambda x: "strt" not in x and "_04"  not in x and "_06" not in x),"parnme"]

    d = (pe_cycle2.loc[:,fpar_names] - pr_pe.loc[:,fpar_names]).apply(lambda x: np.abs(x))
    print(d)
    print(d.max())
    print(d.max().max())
    assert d.max().max() < 1.0e-6


def seq_10par_diff_obspar_cycle_test():
    test_d = "10par_xsec"
    t_d = os.path.join(test_d, "template_diff_parobs")
    tpl_files = [f for f in os.listdir(t_d) if f.lower().endswith(".tpl")]
    in_files = [f.replace(".tpl","") if "cycle" not in f else ".".join([x for x in f.split(".")[:-1] if "cycle" not in x]) for f in tpl_files]

    ins_files = [f for f in os.listdir(t_d) if f.lower().endswith(".ins")]
    out_files = [f.replace(".ins","") if "cycle" not in f else ".".join([x for x in f.split(".")[:-1] if "cycle" not in x]) for f in ins_files]


    tpl_files = [os.path.join(t_d,x) for x in tpl_files]
    in_files = [os.path.join(t_d, x) for x in in_files]
    ins_files = [os.path.join(t_d, x) for x in ins_files]
    out_files = [os.path.join(t_d, x) for x in out_files]

    print(tpl_files)
    print(in_files)
    print(ins_files)
    print(out_files)

    pst = pyemu.Pst.from_io_files(tpl_files,in_files,ins_files,out_files,pst_path=".")
    pst.model_command = "mfnwt 10par_xsec.nam"
    out = pst.model_output_data
    out.loc[out.pest_file.str.contains("cycle"),"cycle"] = out.loc[out.pest_file.str.contains("cycle"),"pest_file"].apply(lambda x: int(x.split('.')[-2].split('_')[-1]))
    print(out)

    out = pst.model_input_data
    out.loc[:,"cycle"] = -1
    out.loc[out.pest_file.str.contains("cycle"), "cycle"] = out.loc[
        out.pest_file.str.contains("cycle"), "pest_file"].apply(lambda x: int(x.split('.')[-2].split('_')[-1]))
    print(out)

    par = pst.parameter_data
    spars = par.loc[par.parnme.str.startswith("ss"),"parnme"]
    assert spars.shape[0] > 0
    print("----->before")
    print(par.loc[spars,["parval1","parubnd","parlbnd"]])
    par.loc[spars,"partrans"] = "fixed"
    par = pst.parameter_data
    par["parubnd"] = par.parval1.values * 1.5
    par["parlbnd"] = par.parval1.values * 0.5

    par.loc[:,"cycle"] = -1
    par.loc[par.parnme.str.startswith("cnhd"),"cycle"] = par.loc[par.parnme.str.startswith("cnhd"),"parnme"].apply(lambda x: int(x.split('_')[-1]))
    print(par)
    print("----->after")
    print(par.loc[spars,["parval1","parubnd","parlbnd"]])
    
    
    obs = pst.observation_data
    obs.loc[:,"cycle"] = obs.obsnme.apply(lambda x: int(x.split('_')[-1]))

    obs.loc[:,"state_par_link"] = ""
    obs.loc[:,"cellid"] = obs.obsnme.apply(lambda x: int(x.split('_')[1]))
    obs.loc[:, "kper"] = obs.obsnme.apply(lambda x: int(x.split('_')[0][1:]))
    spar = par.loc[par.parnme.str.contains("str"),:].copy()
    spar.loc[:,"cellid"] = spar.parnme.apply(lambda x: int(x.split('_')[-1]))
    sobs = obs.loc[obs.kper==1,:]
    sobs = sobs.loc[sobs.cellid != 1,:]
    print(sobs.cellid)
    print(spar.cellid)
    obs.loc[sobs.obsnme,"state_par_link"] = sobs.cellid.apply(lambda x: spar.loc[spar.cellid==x,"parnme"].values[0])
    print(obs)
    pst.control_data.noptmax = 0
    pst.write(os.path.join(t_d, "pest_seq.pst"), version=2)
    t_d1 = t_d.replace("template","test1")
    if os.path.exists(t_d1):
        shutil.rmtree(t_d1)
    shutil.copytree(t_d,t_d1)
    pyemu.os_utils.run("{0} pest_seq.pst".format(exe_path), cwd=t_d1)

    c0_pe = pd.read_csv(os.path.join(t_d1,"pest_seq.global.0.pe.csv"),index_col=0)
    d = np.abs(spar.loc[:,"parval1"].values - c0_pe.loc[c0_pe.index[0],spar.parnme].values)
    print(d)
    print(d.sum())
    assert d.sum() == 0

    c1_pe = pd.read_csv(os.path.join(t_d1,"pest_seq.global.1.pe.csv"),index_col=0)
    res = pyemu.pst_utils.read_resfile(os.path.join(t_d1,"pest_seq.0.base.rei"))
    sres = res.loc[sobs.loc[obs.cycle==0,"obsnme"],:]
    sres = sres.loc[obs.loc[sres.index,"state_par_link"].apply(lambda x: len(x) > 0),:]
    smod = sres.modelled.copy()
    smod.index = smod.index.map(lambda x: obs.loc[x,"state_par_link"])
    print(smod)
    d = np.abs(smod.values - c1_pe.loc[:,smod.index].values)
    print(d)
    print(d.sum())
    assert d.sum() == 0

    pst.control_data.noptmax = -1
    pst.pestpp_options["da_num_reals"] = 10
    pst.pestpp_options["da_lambda_mults"] = [0.1,1.0]
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.write(os.path.join(t_d, "pest_seq.pst"), version=2)
    m_d1 = t_d.replace("template", "master1")
    if os.path.exists(m_d1):
        shutil.rmtree(m_d1)
    #pyemu.os_utils.start_workers(t_d,exe_path,"pest_seq.pst",5,worker_root=test_d,
    #                             master_dir=m_d1)
    shutil.copytree(t_d,m_d1)
    # bd = os.getcwd()
    # os.chdir(m_d1)
    # os.system("{0} pest_seq.pst".format(exe_path))
    # os.chdir(bd)
 
    pyemu.os_utils.run("{0} pest_seq.pst".format(exe_path),cwd=m_d1)
    c0_pe = pd.read_csv(os.path.join(t_d1, "pest_seq.global.0.pe.csv"), index_col=0)
    d = np.abs(spar.loc[:, "parval1"].values - c0_pe.loc["base", spar.parnme].values)
    print(d)
    print(d.sum())
    assert d.sum() == 0

    c1_pe = pd.read_csv(os.path.join(t_d1, "pest_seq.global.1.pe.csv"), index_col=0)
    res = pyemu.pst_utils.read_resfile(os.path.join(t_d1, "pest_seq.0.base.rei"))
    sres = res.loc[sobs.loc[obs.cycle == 0, "obsnme"], :]
    sres = sres.loc[obs.loc[sres.index, "state_par_link"].apply(lambda x: len(x) > 0), :]
    smod = sres.modelled.copy()
    smod.index = smod.index.map(lambda x: obs.loc[x, "state_par_link"])
    print(smod)
    d = np.abs(smod.values - c1_pe.loc["base", smod.index].values)
    print(d)
    print(d.sum())
    assert d.sum() == 0

    c1_oe = pd.read_csv(os.path.join(t_d1, "pest_seq.global.1.oe.csv"), index_col=0)
    d = (c1_pe.loc[:,"cnhd_01_01"] - c1_oe.loc[:,"h01_01_01"]).apply(lambda x: np.abs(x))
    print(d)
    print(d.max())
    assert d.max() == 0
    d = (c1_pe.loc[:, "cnhd_01_00"] - c1_oe.loc[:, "h01_01_00"]).apply(lambda x: np.abs(x))
    print(d)
    print(d.max())
    assert d.max() == 0

def seq_10par_xsec_hotstart_test():
    
    test_d = "10par_xsec"
    t_d = os.path.join(test_d, "template")
    pst = pyemu.Pst(os.path.join(t_d,"pest.pst"))
    par = pst.parameter_data
    par.loc[:,"cycle"] = -1
    
    par.loc[par.parnme.str.contains("strt"),"partrans"] = "log"
    par.loc[["cnhd_01","strt_02","strt_03"],"partrans"] = "fixed"
    strt_pars = par.loc[par.pargp=="strt","parnme"].tolist()
    obs = pst.observation_data
    obs.loc[obs.obsnme.str.startswith("h01"),"weight"] = 1.0
    obs.loc[:,"state_par_link"] = ""
    obs.loc[obs.obgnme=="head1","state_par_link"] = strt_pars
    obs.loc[:,"cycle"] = -1
    
    pst.control_data.noptmax = 1

    pst.model_input_data.loc[:,"cycle"] = -1
    pst.model_output_data.loc[:,"cycle"] = -1

    def get_loc(pst):

        loc = pd.DataFrame(index=pst.nnz_obs_names,columns=pst.adj_par_names)
        loc.loc[:,:] = 0.0
        ocells = loc.index.map(lambda x: int(x.split('_')[1]))
        for pname in pst.adj_par_names:
            cstr = pname.split('_')[1]
            cint = int(cstr)
            # strt states can comm with all obs locs
            if "strt" in pname:
                dist = (ocells - cint).map(np.abs)
                loc_vals = 1 / (dist + 0.01)
                loc_vals = loc_vals.values
                loc_vals[loc_vals>1.0] = 1.0
                loc.loc[:,pname] = loc_vals
            # static pars can only comm with obs in the same cell
            else:

                oname = [o for o in pst.nnz_obs_names if cstr in o.split('_')[1] == cstr][0]
                loc.loc[oname, pname] = 1.0
        return loc

    loc = get_loc(pst)
    pyemu.Matrix.from_dataframe(loc).to_ascii(os.path.join(t_d,"loc.mat"))

    mx_cycle = 5
    cycles = np.arange(0,mx_cycle)
    odf = pd.DataFrame(index=cycles,columns=pst.nnz_obs_names)
    odf.loc[:,:] = obs.loc[pst.nnz_obs_names,"obsval"].values
    odf.T.to_csv(os.path.join(t_d,"obs_cycle_tbl.csv"))
    wdf = pd.DataFrame(index=cycles,columns=pst.nnz_obs_names)
    wdf.loc[:,:] = 0.0
    wdf.iloc[4,:] = 1.0
    wdf.T.to_csv(os.path.join(t_d,"weight_cycle_tbl.csv"))

    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["da_lambda_mults"] = 1.0
    pst.pestpp_options["da_use_mda"] = True
    pst.pestpp_options["da_observation_cycle_table"] = "obs_cycle_tbl.csv"
    pst.pestpp_options["da_weight_cycle_table"] = "weight_cycle_tbl.csv"
    pst.pestpp_options["da_num_reals"] = 10
    pst.pestpp_options["da_localizer"] = "loc.mat"
    pst.pestpp_options["da_hotstart_cycle"] = 3
    pst.pestpp_options["da_subset_how"] = "first"

    pst.write(os.path.join(t_d,"pest_seq.pst"),version=2)
    m_d = os.path.join(test_d, "master_da_hotstart_base")
    pyemu.os_utils.start_workers(t_d, exe_path.replace("ies", "da"), "pest_seq.pst",
                                num_workers=pst.pestpp_options["da_num_reals"], worker_root=test_d, port=port,
                                master_dir=m_d, verbose=True)

    pr_pe = pd.read_csv(os.path.join(m_d,"pest_seq.global.prior.pe.csv"),index_col=0)
    print("testing")
    for cycle in cycles[pst.pestpp_options["da_hotstart_cycle"]+1:]:
        oe = pd.read_csv(os.path.join(m_d,"pest_seq.global.{0}.oe.csv".format(cycle-1)),index_col=0)
        pe = pd.read_csv(os.path.join(m_d,"pest_seq.{0}.0.par.csv".format(cycle)),index_col=0)
        print(os.path.join(m_d,"pest_seq.global.{0}.oe.csv".format(cycle-1)))
        print(os.path.join(m_d, "pest_seq.{0}.0.par.csv".format(cycle)))
        print(oe)
        print(pe)
        #check that that adj dyn states are being updated and used...
        d = np.abs(oe.loc[:,"h01_02"].values - pe.loc[:,"strt_02"].values)
        print(d)
        assert d.max() < 1.0e-6,d.max()
        d = np.abs(oe.loc[:, "h01_03"].values - pe.loc[:, "strt_03"].values)
        print(d)
        assert d.max() < 1.0e-6, d.max()
        # check that the fixed par (non state) is not being adjusted
        d =  np.abs(pr_pe.loc[:,"cnhd_01"].values - pe.loc[:, "cnhd_01"].values)
        print(d)
        assert d.max() < 1.0e-6, d.max()


    pe_cycle2 = pd.read_csv(os.path.join(m_d,"pest_seq.global.3.pe.csv"),index_col=0)
    fpar_names = par.loc[par.parnme.apply(lambda x: "strt" not in x),"parnme"]

    d = (pe_cycle2.loc[:,fpar_names] - pr_pe.loc[:,fpar_names]).apply(lambda x: np.abs(x))
    print(d)
    print(d.max())
    print(d.max().max())
    assert d.max().max() < 1.0e-6
    phi1 = pd.read_csv(os.path.join(m_d, "pest_seq.global.phi.actual.csv"))

    pe_file = os.path.join(m_d,"pest_seq.global.3.pe.csv")
    shutil.copy2(pe_file,os.path.join(t_d,"hs_pe.csv"))
    shutil.copy2(pe_file.replace(".pe.",".oe."),os.path.join(t_d,"hs_oe.csv"))
    shutil.copy2(os.path.join(m_d,"pest_seq.3.obs+noise.csv"),os.path.join(t_d,"hs_noise.csv"));
    pst.pestpp_options["da_par_en"] = "hs_pe.csv"
    pst.pestpp_options["da_restart_obs_en"] = "hs_oe.csv"
    pst.pestpp_options["da_obs_en"] = "hs_noise.csv"

    pst.write(os.path.join(t_d,"pest_seq.pst"),version=2)
    m_d = os.path.join(test_d, "master_da_hotstart_restart")
    pyemu.os_utils.start_workers(t_d, exe_path.replace("ies", "da"), "pest_seq.pst",
                                num_workers=pst.pestpp_options["da_num_reals"], worker_root=test_d, port=port,
                                master_dir=m_d, verbose=True)

    phi2 = pd.read_csv(os.path.join(m_d,"pest_seq.global.phi.actual.csv"))
    d = (phi1-phi2).apply(lambda x: np.abs(x))
    print(d.max())
    print(d.max().max())
    assert d.max().max() < 0.01


def seq_10par_cycle_parse_test():
    
    test_d = "10par_xsec"
    t_d = os.path.join(test_d, "template")
    pst = pyemu.Pst(os.path.join(t_d,"pest.pst"))

    new_ins_file = os.path.join(t_d, "firstcycle.hds.ins")
    org_ins_file = os.path.join(t_d, "10par_xsec.hds.ins")
    lines = open(org_ins_file,'r').readlines()
    with open(new_ins_file, 'w') as f:
        for line in lines:
            f.write(line.replace("h","firstcycle"))

    df = pst.add_observations(new_ins_file,"10par_xsec.hds",pst_path=".")
    obs = pst.observation_data
    obs.loc[obs.obsnme.str.startswith("firstcycle"),"cycle"] = 0
    obs.loc[obs.obsnme.str.startswith("firstcycle"), "weight"] = 1



    tpl_file = os.path.join(t_d,"every_other_cycle.dat.tpl")
    with open(tpl_file,'w') as f:
        f.write("ptf ~\n")
        f.write("every_other_cycle  ~  every_other_cycle  ~\n")
    pst.add_parameters(tpl_file,pst_path=".")
    par = pst.parameter_data
    

    tpl_file = os.path.join(t_d,"once_and_a_while.dat.tpl")
    with open(tpl_file,'w') as f:
        f.write("ptf ~\n")
        f.write("every_other_cycle  ~  once_and_a_while  ~\n")
    pst.add_parameters(tpl_file,pst_path=".")

    df = pst.add_parameters(os.path.join(t_d,"hk_Layer_1_even.ref.tpl"),os.path.join(t_d,"hk_Layer_1.ref"),pst_path=".")
    
    par = pst.parameter_data
    par.loc[:,"cycle"] = -1
    par.loc[par.parnme.str.startswith("k_"),"cycle"] = "1:8:2"
    par.loc[df.parnme,"cycle"] = "::2"
    par.loc[df.parnme,"parval1"] = 2.5
    par.loc[df.parnme,"parlbnd"] = .25
    par.loc[df.parnme,"parubnd"] = 25.

    par.loc["every_other_cycle","cycle"] = "001:005:002"
    par.loc["every_other_cycle","parval1"] = 1.0
    par.loc["every_other_cycle","parubnd"] = 10
    par.loc["every_other_cycle","parlbnd"] = 0.1

    par.loc["once_and_a_while","cycle"] = "003::004"
    par.loc["once_and_a_while","parval1"] = 1.0
    par.loc["once_and_a_while","parubnd"] = 10
    par.loc["once_and_a_while","parlbnd"] = 0.1


    par.loc[par.parnme.str.contains("strt"),"partrans"] = "log"
    par.loc[["cnhd_01","strt_02","strt_03"],"partrans"] = "fixed"
    strt_pars = par.loc[par.pargp=="strt","parnme"].tolist()
    obs = pst.observation_data
    obs.loc[obs.obsnme.str.startswith("h01"),"weight"] = 1.0
    obs.loc[:,"state_par_link"] = ""
    obs.loc[obs.obgnme=="head1","state_par_link"] = strt_pars
    obs.loc[~obs.obsnme.str.startswith("firstcycle"),"cycle"] = -1
    
    pst.control_data.noptmax = 3

    mid = pst.model_input_data
    mid.loc[:,"cycle"] = -1
    mid.loc[mid.model_file.str.contains("every"),"cycle"] = par.loc["every_other_cycle","cycle"]
    mid.loc[mid.model_file.str.contains("once"),"cycle"] = par.loc["once_and_a_while","cycle"]
    mid.loc[mid.pest_file.str.contains("_even"),"cycle"] = par.loc[df.parnme,"cycle"].values[0]
    mid.loc[mid.pest_file.str.contains("hk_Layer_1.ref"),"cycle"] = par.loc[par.parnme.str.startswith("k_"),"cycle"].values[0]



    pst.model_output_data.loc[:,"cycle"] = -1
    mod = pst.model_output_data
    mod.loc[mod.pest_file.str.contains("firstcycle"), "cycle"] = 0

    def get_loc(pst):

        loc = pd.DataFrame(index=pst.nnz_obs_names,columns=pst.adj_par_names)
        loc.loc[:,:] = 0.0
        ocells = loc.index.map(lambda x: int(x.split('_')[1]))
        for pname in pst.adj_par_names:
            if pname in ["every_other_cycle","once_and_a_while"]:
                loc.loc[:,pname] = 0.0
                continue
            cstr = pname.split('_')[1]
            cint = int(cstr)
            # strt states can comm with all obs locs
            if "strt" in pname:
                dist = (ocells - cint).map(np.abs)
                loc_vals = 1 / (dist + 0.01)
                loc_vals = loc_vals.values
                loc_vals[loc_vals>1.0] = 1.0
                loc.loc[:,pname] = loc_vals
            # static pars can only comm with obs in the same cell
            else:

                oname = [o for o in pst.nnz_obs_names if cstr in o.split('_')[1] == cstr][0]
                loc.loc[oname, pname] = 1.0
        return loc

    loc = get_loc(pst)
    pyemu.Matrix.from_dataframe(loc).to_ascii(os.path.join(t_d,"loc.mat"))

    mx_cycle = 10
    cycles = np.arange(0,mx_cycle)
    tbl_onames = [n for n in pst.nnz_obs_names if not n.startswith("firstcycle")]
    odf = pd.DataFrame(index=cycles,columns=tbl_onames)
    odf.loc[:,:] = obs.loc[tbl_onames,"obsval"].values
    odf.T.to_csv(os.path.join(t_d,"obs_cycle_tbl.csv"))
    wdf = pd.DataFrame(index=cycles,columns=tbl_onames)
    wdf.loc[:,:] = 0.0
    wdf.iloc[1,[3,5]] = 1.0
    wdf.iloc[3,:] = 1.0
    wdf.T.to_csv(os.path.join(t_d,"weight_cycle_tbl.csv"))

    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["da_lambda_mults"] = 1.0
    pst.pestpp_options["da_use_mda"] = True
    pst.pestpp_options["da_observation_cycle_table"] = "obs_cycle_tbl.csv"
    pst.pestpp_options["da_weight_cycle_table"] = "weight_cycle_tbl.csv"
    pst.pestpp_options["da_num_reals"] = 10
    pst.pestpp_options["da_localizer"] = "loc.mat"

    pst.write(os.path.join(t_d,"pest_seq.pst"),version=2)
    m_d = os.path.join(test_d, "master_da_cycle_parse")

    pyemu.os_utils.start_workers(t_d, exe_path.replace("ies", "da"), "pest_seq.pst",
                                num_workers=pst.pestpp_options["da_num_reals"], worker_root=test_d, port=port,
                                master_dir=m_d, verbose=True)

    final_oe = os.path.join(m_d,"pest_seq.global.{0}.oe.csv".format(mx_cycle-1))
    assert os.path.exists(final_oe),final_oe


def compare_mf6_freyberg():
    def mod_tdis_sto(org_t_d,t_d):
        tdis_file = "freyberg6.tdis"
        lines = open(os.path.join(org_t_d,tdis_file),'r').readlines()

        with open(os.path.join(t_d,tdis_file),'w') as f:
            iline = 0
            while True:
                line = lines[iline]
                if "begin period" in line.lower():
                    lines[iline+1] = "100000  1   1.0\n"
                    print(lines[iline+1])
                print(line)
                f.write(line)
                iline += 1
                if iline >= len(lines):
                    break
        sto_file = "freyberg6.sto"
        lines = open(os.path.join(org_t_d,sto_file),'r').readlines()

        with open(os.path.join(t_d,sto_file),'w') as f:

            for line in lines:
                f.write(line)
                if line.lower().startswith("end griddata"):
                    break
            f.write("\nbegin period 1\n  transient\nend period 1\n")

    # prep that prior ensemble for da
    da_test_d = "mf6_freyberg"
    da_t_d = os.path.join(da_test_d, "template_seq_native")
    org_da_t_d = da_t_d
    da_t_d = da_t_d + "_compare"
    if os.path.exists(da_t_d):
        shutil.rmtree(da_t_d)
    shutil.copytree(org_da_t_d,da_t_d)
    da_pst = pyemu.Pst(os.path.join(da_t_d,"freyberg6_run_da2.pst"))
    par = da_pst.parameter_data
    par.loc[par.parnme.str.contains("welflx"),"scale"] = -1.0
    par.loc["perlen","parval1"] = 100000
    #par.loc[~par.parnme.str.contains("head"),"partrans"] = "fixed"
      
    ies_test_d = "mf6_freyberg"
    ies_t_d = os.path.join(ies_test_d, "template")
    org_ies_t_d = ies_t_d
    ies_t_d = ies_t_d + "_compare"
    if os.path.exists(ies_t_d):
        shutil.rmtree(ies_t_d)
    shutil.copytree(org_ies_t_d,ies_t_d)
    ies_pst = pyemu.Pst(os.path.join(ies_t_d,"freyberg6_run_ies.pst"))

    mod_tdis_sto(org_ies_t_d,ies_t_d)
    
    pyemu.os_utils.run("mf6",cwd=ies_t_d)
    pyemu.os_utils.run("mf6",cwd=da_t_d)

    ies_pst.pestpp_options["ies_par_en"] = "prior.jcb"
    ies_pe = pyemu.ParameterEnsemble.from_binary(pst=ies_pst,filename=os.path.join(ies_t_d,"prior.jcb"))
    da_pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=da_pst,
        cov=pyemu.Cov.from_parameter_data(da_pst),num_reals=ies_pe.shape[0])
    da_pe.index = ies_pe.index
    d = set(da_pe.columns.tolist()).symmetric_difference(set(ies_pe.columns.tolist()))
    print(d)
    da_pe.loc[:,ies_pe.columns] = ies_pe.values
    da_pe.to_binary(os.path.join(da_t_d,"da_prior.jcb"))
    da_pst.pestpp_options["ies_par_en"] = "da_prior.jcb"

    # turn on all of the org freyberg model obs locations for the first year
    ies_obs = ies_pst.observation_data
    # ies_obs.loc[:,"weight"] = 0.0
    tr_obs = ies_obs.loc[ies_obs.obsnme.apply(lambda x: x.startswith("trgw") and "2017" not in x),:].copy()
    # ies_obs.loc[tr_obs.obsnme,"weight"] = 6
    gage_obs = ies_obs.loc[ies_obs.obsnme.apply(lambda x: x.startswith("gage")),:].copy()
    gage_obs.sort_index(inplace=True)

    mults = np.linspace(0.5,0.5,gage_obs.shape[0])
    new_vals = gage_obs.obsval.values * mults

    # import matplotlib.pyplot as plt
    # plt.plot(gage_obs.obsval.values)
    # plt.plot(new_vals)
    # plt.show()
    # return
    ies_obs.loc[:,"org_obsval"] = np.nan
    ies_obs.loc[gage_obs.obsnme,"org_obsval"] = gage_obs.obsval.values
    ies_obs.loc[gage_obs.obsnme,"obsval"] = new_vals

    #ies_obs.loc[[o for o in ies_pst.nnz_obs_names if not "gage" in o],"weight"] = 50
    #ies_obs.loc[[o for o in ies_pst.nnz_obs_names if "gage" in o],"weight"] = 1.0 / (ies_obs.loc[[o for o in ies_pst.nnz_obs_names if "gage" in o],"obsval"] * 0.15)

    # use this to check the phi contribs
    # ies_pst.control_data.noptmax = 0
    # ies_pst.write(os.path.join(ies_t_d,"reweight.pst"))
    # pyemu.os_utils.run("{0} reweight.pst".format(exe_path.replace("-da","-ies")),cwd=ies_t_d)
    # ies_pst.set_res(os.path.join(ies_t_d,"reweight.base.rei"))
    # print({n:p for n,p in ies_pst.phi_components.items() if n in ies_pst.nnz_obs_groups})
    # ies_pst.plot(kind="phi_pie")
    # import matplotlib.pyplot as plt
    # plt.show()
    # return

    tr_obs.loc[:,"k"] = tr_obs.obsnme.apply(lambda x: int(x.split('_')[1]))
    tr_obs.loc[:, "i"] = tr_obs.obsnme.apply(lambda x: int(x.split('_')[2]))
    tr_obs.loc[:, "j"] = tr_obs.obsnme.apply(lambda x: int(x.split('_')[3]))
    tr_obs.loc[:, "kij"] = tr_obs.apply(lambda x: (x.k,x.i,x.j),axis=1)
    ukij = set(tr_obs.kij.unique().tolist())
    da_obs = da_pst.observation_data
    hd_obs = da_obs.loc[da_obs.obsnme.str.startswith("head_"),:].copy()
    print(hd_obs.obsnme)
    hd_obs.loc[:,"k"]= hd_obs.obsnme.apply(lambda x: int(x.split('_')[1]))
    hd_obs.loc[:, "i"] = hd_obs.obsnme.apply(lambda x: int(x.split('_')[2]))
    hd_obs.loc[:, "j"] = hd_obs.obsnme.apply(lambda x: int(x.split('_')[3]))
    hd_obs.loc[:,"kij"] = hd_obs.apply(lambda x: (x.k,x.i,x.j),axis=1)
    hd_obs.loc[:,"org_obgnme"] = hd_obs.apply(lambda x: "trgw_{0}_{1}_{2}".format(x.k,x.i,x.j),axis=1)

    include = hd_obs.loc[hd_obs.kij.apply(lambda x: x in ukij),"obsnme"].tolist()
    #include.append("gage_1")
    otbl = pd.DataFrame(columns=include,index=np.arange(25)).T
    wtbl = pd.DataFrame(columns=include,index=np.arange(25)).T
    for uog in hd_obs.loc[include,"org_obgnme"].unique():
        oname = hd_obs.loc[hd_obs.org_obgnme==uog,"obsnme"].values[0]
        og_obs = ies_obs.loc[ies_obs.obgnme==uog,:]
        og_obs.sort_index(inplace=True)
        otbl.loc[oname,:] = og_obs.obsval.values
        wtbl.loc[oname,:] = og_obs.weight.values

    otbl.loc["gage_1",:] = ies_obs.loc[ies_obs.obsnme.str.contains("gage"),"obsval"].values
    wtbl.loc["gage_1",:] = ies_obs.loc[ies_obs.obsnme.str.contains("gage"), "weight"].values
    otbl.to_csv(os.path.join(da_t_d,"obs_cycle_tbl.csv"))
    wtbl.to_csv(os.path.join(da_t_d, "weight_cycle_tbl.csv"))
    da_pst.observation_data.loc[otbl.index.values,"weight"] = 1.0
    ies_pst.pestpp_options.pop("ies_num_reals",None)
    da_pst.pestpp_options.pop("da_num_reals",None)
    ies_pst.pestpp_options.pop("da_num_reals",None)
    da_pst.pestpp_options.pop("ies_num_reals",None)
    ies_pst.pestpp_options["ies_no_noise"] = False
    da_pst.pestpp_options["ies_no_noise"] = False
    da_pst.pestpp_options["ies_verbose_level"] = 1
    ies_pst.pestpp_options["ies_verbose_level"] = 1
    ies_pst.pestpp_options.pop("ies_localizer",None)
    da_pst.pestpp_options.pop("ies_localizer", None)
    ies_pst.pestpp_options["ies_autoadaloc"] = False
    da_pst.pestpp_options["ies_autoadaloc"] = False
    ies_pst.pestpp_options["ies_save_lambda_en"] = False
    da_pst.pestpp_options["ies_save_lambda_en"] = False
    ies_pst.pestpp_options["ies_drop_conflicts"] = True
    da_pst.pestpp_options["ies_drop_conflicts"] = True
    da_pst.pestpp_options["ies_num_reals"] = 50
    ies_pst.pestpp_options["ies_num_reals"] = 50

    #da_pst.pestpp_options["da_stop_cycle"] = 1

    ies_pst.control_data.noptmax = 3
    da_pst.control_data.noptmax = 3

    # # run da
    da_pst.pestpp_options["ies_use_mda"] = True
    da_pst.write(os.path.join(da_t_d, "freyberg6_run_da.pst"), version=2)
    da_m_d_mda = os.path.join(da_test_d, "master_da_mda")
    pyemu.os_utils.start_workers(da_t_d, exe_path.replace("-ies", "-da"), "freyberg6_run_da.pst",
                                 num_workers=20, worker_root=da_test_d, port=port,
                                 master_dir=da_m_d_mda, verbose=True)
    da_pst.pestpp_options["ies_use_mda"] = False
    da_pst.write(os.path.join(da_t_d,"freyberg6_run_da.pst"),version=2)
    da_m_d_glm = os.path.join(da_test_d, "master_da_glm")
    pyemu.os_utils.start_workers(da_t_d, exe_path.replace("-ies","-da"), "freyberg6_run_da.pst",
                                num_workers=20, worker_root=da_test_d, port=port,
                                master_dir=da_m_d_glm, verbose=True)

    # run ies    
    ies_pst.pestpp_options["ies_use_mda"] = False
    ies_pst.write(os.path.join(ies_t_d,"freyberg6_run_ies.pst"),version=2)
    ies_m_d_glm = os.path.join(ies_test_d, "master_ies_glm")
    pyemu.os_utils.start_workers(ies_t_d, exe_path.replace("-da","-ies"), "freyberg6_run_ies.pst",
                                num_workers=20, worker_root=ies_test_d, port=port,
                                master_dir=ies_m_d_glm, verbose=True)

    ies_pst.pestpp_options["ies_use_mda"] = True
    ies_pst.write(os.path.join(ies_t_d,"freyberg6_run_ies.pst"),version=2)
    ies_m_d_mda = os.path.join(ies_test_d, "master_ies_mda")
    pyemu.os_utils.start_workers(ies_t_d, exe_path.replace("-da","-ies"), "freyberg6_run_ies.pst",
                                num_workers=20, worker_root=ies_test_d, port=port,
                                master_dir=ies_m_d_mda, verbose=True)


def plot_compare(solution="ies",noptmax=1):
    import matplotlib.pyplot as plt
    da_m_d = os.path.join("mf6_freyberg", "master_da_{0}".format(solution))
    ies_m_d = os.path.join("mf6_freyberg", "master_ies_{0}".format(solution))
    ies_case = "freyberg6_run_ies"
    da_case = "freyberg6_run_da"

    ies_pst = pyemu.Pst(os.path.join(ies_m_d,ies_case+".pst"))
    ies_obs = ies_pst.observation_data#.loc[ies_pst.nnz_obs_names,:]
    ies_obs = ies_obs.loc[ies_obs.obgnme.apply(lambda x: x in ies_pst.nnz_obs_groups),:]
    print(ies_obs)
    ies_obs.loc[:,"datetime"] = pd.to_datetime(ies_obs.obsnme.apply(lambda x: x.split('_')[-1]),format='%Y%m%d')

    ies_pr_oe = pd.read_csv(os.path.join(ies_m_d,ies_case+".0.obs.csv"))
    ies_pt_oe = pd.read_csv(os.path.join(ies_m_d,ies_case+".{0}.obs.csv".format(noptmax)))

    da_pst = pyemu.Pst(os.path.join(da_m_d,da_case+".pst"))
    da_obs = da_pst.observation_data.loc[da_pst.nnz_obs_names,:].copy()
    da_obs.loc[da_obs.obsnme.str.contains("gage"),"org_obgnme"] = "gage"
    print(da_obs)

    #da_pr_oe = pd.read_csv(os.path.join(da_m_d,da_case+".0.obs.csv"))
    num_cycles = 25
    da_pr_dict = {}
    da_pt_dict = {}
    for cycle in range(num_cycles):
        print(cycle)
        da_pr_oe = pd.read_csv(os.path.join(da_m_d,da_case+".{0}.0.obs.csv".format(cycle)))
        da_pr_dict[cycle] = da_pr_oe
        pt_file = os.path.join(da_m_d,da_case+".{0}.{1}.obs.csv".format(cycle,noptmax))
        if (os.path.exists(pt_file)):
            da_pt_oe = pd.read_csv(pt_file)
            da_pt_dict[cycle] = da_pt_oe
        else:
            print("missing posterior",cycle)
    ies_og_uvals = ies_obs.obgnme.unique()
    print(ies_og_uvals)
    plt_d = "compare_plots_{0}_{1}".format(solution,noptmax)
    if not os.path.exists(plt_d):
        os.mkdir(plt_d)
    ies_og_uvals.sort()


    for og in ies_og_uvals:
        ies_obs_og = ies_obs.loc[ies_obs.obgnme == og, :].copy()
        ies_obs_og.sort_values(by="datetime", inplace=True)
        da_obs_og = da_obs.loc[da_obs.org_obgnme == og, :]
        fig_count = 0
        for cycle in range(num_cycles):

            dts = ies_obs_og.datetime.values

            def make_plot(axes):
                ax = axes[0]
                ax.set_title("smoother formulation, observation location: " + og, loc="left")
                [ax.plot(dts,ies_pr_oe.loc[idx,ies_obs_og.obsnme],"0.5",alpha=0.5,lw=0.25) for idx in ies_pr_oe.index]
                [ax.plot(dts,ies_pt_oe.loc[idx,ies_obs_og.obsnme],"b",alpha=0.5,lw=0.5) for idx in ies_pt_oe.index]
                ax.plot(dts, ies_pr_oe.loc[ies_pr_oe.index[0], ies_obs_og.obsnme], "0.5", alpha=0.5,
                        lw=0.1,label="prior real")
                ax.plot(dts, ies_pt_oe.loc[ies_pt_oe.index[0], ies_obs_og.obsnme], "b", alpha=0.5,
                        lw=0.1,label="post real")
                ax.plot(dts,ies_obs_og.obsval,"r",label="truth")
                ies_obs_nz = ies_obs_og.loc[ies_obs_og.weight>0,:]

                ax.scatter(ies_obs_nz.datetime.values, ies_obs_nz.obsval, marker="^", color="r",s=50,
                           zorder=10,label="obs")
                ax.legend(loc="upper right")
                ax = axes[1]

                ax.set_title("filter formulation (monthly cycles), observation location: " + da_obs_og.obsnme.values[0],
                             loc="left")
                ax.plot(dts, ies_obs_og.obsval, "r", label="truth")
                ax.scatter(ies_obs_nz.datetime.values, ies_obs_nz.obsval, marker="^", color="r",
                           s=50, zorder=10, label="obs")

                post_labeled = False
                for ccycle in range(cycle):
                    da_pr_oe = da_pr_dict[ccycle]
                    label = None
                    if ccycle == 0:
                        label = "prior real"
                    ax.scatter([dts[ccycle] for _ in range(da_pr_oe.shape[0])],
                               da_pr_oe.loc[:, da_obs_og.obsnme[0]].values,
                               marker=".", color="0.5", label=label)

                    if ccycle in da_pt_dict:
                        da_pt_oe = da_pt_dict[ccycle]
                        label = None
                        if not post_labeled:
                            label = "post real"
                            post_labeled = True
                        ax.scatter([dts[ccycle] for _ in range(da_pt_oe.shape[0])],
                                   da_pt_oe.loc[:, da_obs_og.obsnme[0]].values,
                                   marker=".", color="b", label=label)
                    ax.set_ylim(axes[0].get_ylim())
                    ax.legend(loc="upper right")
                #ax.set_title("observation location: "+og,loc="left")
                #ax.legend(loc="upper right")

            # prior
            fig, axes = plt.subplots(2, 1, figsize=(8, 8))
            make_plot(axes)
            da_pr_oe = da_pr_dict[cycle]
            ax = axes[1]
            ax.scatter([dts[cycle] for _ in range(da_pr_oe.shape[0])],
                       da_pr_oe.loc[:, da_obs_og.obsnme[0]].values,
                       marker=".", color="0.5")
            ax.set_ylim(axes[0].get_ylim())
            plt.tight_layout()
            plt.savefig(os.path.join(plt_d, "compare_{0}_{1}_{2:03d}.png".format(og, solution, fig_count)))
            plt.close(fig)
            fig_count += 1

            # posterior
            if cycle in da_pt_dict:
                da_pt_oe = da_pt_dict[cycle]
                fig, axes = plt.subplots(2, 1, figsize=(8, 8))
                make_plot(axes)
                da_pr_oe = da_pr_dict[cycle]
                ax.scatter([dts[cycle] for _ in range(da_pr_oe.shape[0])],
                           da_pr_oe.loc[:, da_obs_og.obsnme[0]].values,
                           marker=".", color="0.5")

                ax.scatter([dts[cycle] for _ in range(da_pt_oe.shape[0])],
                           da_pt_oe.loc[:, da_obs_og.obsnme[0]].values,
                           marker=".", color="b")
                ax.set_ylim(axes[0].get_ylim())
                plt.tight_layout()
                plt.savefig(os.path.join(plt_d, "compare_{0}_{1}_{2:03d}.png".format(og, solution, fig_count)))
                plt.close(fig)




            # posterior


    for og in ies_og_uvals:
        pyemu.os_utils.run("ffmpeg -i compare_{0}_{1}_024.png -vf palettegen=16 -y palette.png".format(og,solution),cwd=plt_d)
        cmd = "ffmpeg -i compare_{0}_{1}_%03d.png -i palette.png -y -filter_complex ".format(og,solution)
        cmd += "\"fps=80,scale=720:-1:flags=lanczos[x];[x][1:v]paletteuse\" -y  -final_delay 150 {0}_{1}.gif".format(og,solution)
        pyemu.os_utils.run(cmd,cwd=plt_d)
        #os.system("ffmpeg -r 1 -i {0} -loop 1 -final_delay 100 -y -vf fps=25 {1}.mp4".
        #      format(os.path.join(plt_d,"compare_{0}_{1}_%03d.png".format(og,solution)),
        #             os.path.join(plt_d,solution+"_"+og)))


def da_pareto_demo():

    test_d = "mf6_freyberg"
    t_d = os.path.join(test_d, "template")
    pst = pyemu.Pst(os.path.join(t_d,"freyberg6_run_ies.pst"))
    pst.control_data.noptmax = 0
    #pst.write(os.path.join(t_d,"freyberg6_run_ies.pst"))
    #pyemu.os_utils.run("{0} freyberg6_run_ies.pst".format(exe_path.replace("-da","-ies")),cwd=t_d)

    obs = pst.observation_data
    mx_weight = .05
    mn_weight = 0.05
    mn_foreval = -500
    mx_foreval = 500
    fore_name = "headwater_20171031"
    obs.loc[fore_name,"weight"] = mx_weight
    obs.loc[fore_name,"obsval"] = mx_foreval

    # to find a good max forecast weight
    # import matplotlib.pyplot as plt
    # pst.plot(kind="phi_pie")
    # plt.show()
    # return

    ncycle = 10
    wtbl = pd.DataFrame(index = pst.nnz_obs_names,columns=np.arange(ncycle))
    otbl = pd.DataFrame(index=pst.nnz_obs_names, columns=np.arange(ncycle))
    for col in wtbl.columns:
        wtbl.loc[pst.nnz_obs_names,col] = obs.loc[pst.nnz_obs_names,"weight"].values
        otbl.loc[pst.nnz_obs_names, col] = obs.loc[pst.nnz_obs_names, "obsval"].values

    fore_wghts = np.linspace(mn_weight,mx_weight,ncycle)
    fore_vals = np.linspace(mn_foreval,mx_foreval,ncycle)
    wtbl.loc[fore_name,:] = fore_wghts
    otbl.loc[fore_name,:] = fore_vals

    wtbl.to_csv(os.path.join(t_d,"pareto_weight_table.csv"))
    otbl.to_csv(os.path.join(t_d, "pareto_obs_table.csv"))

    pst.pestpp_options["da_weight_cycle_table"] = "pareto_weight_table.csv"
    pst.pestpp_options["da_observation_cycle_table"] = "pareto_obs_table.csv"

    pst.control_data.noptmax = 3
    pst.pestpp_options["da_verbose_level"] = 3
    #pst.pestpp_options["ies_num_reals"] = 20
    #pst.pestpp_options["ies_lambda_mults"] = 1.0
    #pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options.pop("ies_localizer",None)
    pst.pestpp_options["ies_autoadaloc"] = False
    pst.pestpp_options["ies_no_noise"] = True
    pst.write(os.path.join(t_d,"freyberg6_run_ies_pareto.pst"))
    m_d = os.path.join(test_d,"master_ies_pareto")
    pyemu.os_utils.start_workers(t_d, exe_path, "freyberg6_run_ies_pareto.pst",
                                 num_workers=25, worker_root=test_d, port=port,
                                 master_dir=m_d, verbose=True)


def plot_da_pareto_demo():
    m_d = os.path.join("mf6_freyberg","master_ies_pareto")
    case = "freyberg6_run_ies_pareto"
    pst = pyemu.Pst(os.path.join(m_d,case+".pst"))
    global_oe_files = [f for f in os.listdir(m_d) if case in f and ".global." in f and ".oe." in f]
    print(global_oe_files)
    global_oe_files = {int(f.split('.')[2]):os.path.join(m_d,f) for f in global_oe_files}
    cycles = list(global_oe_files.keys())
    cycles.sort()
    obs = pst.observation_data
    # turn off the forecast weight for phi calcs
    fname = "headwater_20171031"
    obs.loc[fname,"weight"] = 0.0
    import matplotlib.pyplot as plt
    #fig,ax = plt.subplots(1,1,figsize=(10,10))
    fig = plt.figure(constrained_layout=True,figsize=(8,8))
    gs = fig.add_gridspec(4,4)
    ax = fig.add_subplot(gs[:-1,1:])
    ax_phi = fig.add_subplot(gs[:-1,0])
    ax_fore = fig.add_subplot(gs[-1,1:])
    cmap = plt.get_cmap("jet")
    bins = 10
    for cycle in cycles:
        oe = pyemu.ObservationEnsemble.from_csv(pst=pst,filename=global_oe_files[cycle])
        fvec = oe.loc[:,fname]
        pvec = oe.phi_vector
        c = cmap(cycle/max(cycles))
        ax.scatter(fvec,pvec,marker=".",color=c,label="cycle {0}".format(cycle))
        ax_fore.hist(fvec,bins=bins,facecolor=c,edgecolor="none",alpha=0.5)
        ax_phi.hist(pvec, bins=bins, facecolor=c, edgecolor="none", alpha=0.5,orientation="horizontal")
    ax.set_title("B) ensemble-based simulated surface-water/groundwater exchange (SGE) vs $\phi$",loc="left")
    ax_phi.set_title("A) marginal $\phi$ PDFs",loc="left")
    ax_fore.set_title("C) marginal SGE PDFs",loc="left")
    ax.legend(loc="upper left")
    ax_phi.set_ylim(ax.get_ylim())
    ax_phi.set_xticks([])
    ax_fore.set_xlim(ax.get_xlim())
    ax_fore.set_yticks([])
    ax_fore.set_xlabel("sw-gw flux during low-flow period (negative is from gw to sw)".format(fname))
    ax_phi.set_ylabel("sum of squared-weighted residuals ($\phi$)")
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.grid()
    #ax.set_title("A) A demo of direct hypothesis testing for the freyberg model\n"+
    #             "   For each cycle, the value and weight of the forecast was varied\n"+
    #             "   to force the model reduce the gw contrib to sw", loc="left")
    plt.savefig(case+".pdf")
    plt.close(fig)


def seq_10par_xsec_double_state_test():
    """a test of estimating parameters and final states - initial states fixed"""
    test_d = "10par_xsec"
    t_d_org = os.path.join(test_d, "template")
    t_d = os.path.join(test_d, "template_double")
    if os.path.exists(t_d):
        shutil.rmtree(t_d)
    shutil.copytree(t_d_org,t_d)
    init_state_tpl_file = os.path.join(t_d,"strt_Layer_1.ref.tpl")
    lines = open(init_state_tpl_file,'r').readlines()
    final_state_tpl_file = init_state_tpl_file.replace("strt","final_states")
    with open(final_state_tpl_file,'w') as f:
        for line in lines:
            f.write(line.replace("_","final_"))

    pst = pyemu.Pst(os.path.join(t_d, "pest.pst"))
    df = pst.add_parameters(final_state_tpl_file,pst_path=".")
    df.sort_values(by="parnme",inplace=True)
    par = pst.parameter_data
    init_state_pars = df.parnme.str.replace("final_","_").tolist()[1:]
    final_state_pars = [i.replace("_","final_") for i in init_state_pars]
    par.loc[final_state_pars,"parval1"] = par.loc[init_state_pars,"parval1"].values
    par.loc[final_state_pars, "parubnd"] = par.loc[init_state_pars, "parubnd"].values
    par.loc[final_state_pars, "parlbnd"] = par.loc[init_state_pars, "parlbnd"].values
    par.loc[final_state_pars, "pargp"] = "final_states"


    par.loc[:, "cycle"] = -1
    par.loc[:,"state_par_link"] = np.nan
    par.loc[final_state_pars,"state_par_link"] = init_state_pars

    par.loc[par.parnme.str.contains("strt"), "partrans"] = "log"
    par.loc["cnhd_01", "partrans"] = "fixed"
    par.loc[init_state_pars,"partrans"] = "fixed"
    strt_pars = par.loc[par.pargp == "strt", "parnme"].tolist()
    obs = pst.observation_data
    obs.loc[obs.obsnme.str.startswith("h01"), "weight"] = 1.0
    obs.loc[:, "state_par_link"] = ""
    obs.loc[obs.obgnme == "head1", "state_par_link"] = strt_pars
    ostates = obs.loc[obs.obgnme == "head1","obsnme"]
    obs.loc[:, "cycle"] = -1

    pst.control_data.noptmax = 1

    pst.model_input_data.loc[:, "cycle"] = -1
    pst.model_output_data.loc[:, "cycle"] = -1

    def get_loc(pst):

        loc = pd.DataFrame(index=pst.nnz_obs_names, columns=pst.adj_par_names)
        loc.loc[:, :] = 0.0
        ocells = loc.index.map(lambda x: int(x.split('_')[1]))
        for pname in pst.adj_par_names:
            cstr = pname.split('_')[1]
            cint = int(cstr)
            # strt states can comm with all obs locs
            if "strt" in pname:
                dist = (ocells - cint).map(np.abs)
                loc_vals = 1 / (dist + 0.01)
                loc_vals = loc_vals.values
                loc_vals[loc_vals > 1.0] = 1.0
                loc.loc[:, pname] = loc_vals
            # static pars can only comm with obs in the same cell
            else:

                oname = [o for o in pst.nnz_obs_names if cstr in o.split('_')[1] == cstr][0]
                loc.loc[oname, pname] = 1.0
        return loc

    loc = get_loc(pst)
    pyemu.Matrix.from_dataframe(loc).to_ascii(os.path.join(t_d, "loc.mat"))

    # dont change these settings - they are hard coded below in the testing
    mx_cycle = 5
    cycles = np.arange(0, mx_cycle)
    odf = pd.DataFrame(index=cycles, columns=pst.nnz_obs_names)
    odf.loc[:, :] = obs.loc[pst.nnz_obs_names, "obsval"].values
    odf.T.to_csv(os.path.join(t_d, "obs_cycle_tbl.csv"))
    wdf = pd.DataFrame(index=cycles, columns=pst.nnz_obs_names)
    wdf.loc[:, :] = 0.0
    wdf.iloc[1, [3, 5]] = 1.0
    wdf.iloc[3, :] = 1.0
    wdf.T.to_csv(os.path.join(t_d, "weight_cycle_tbl.csv"))

    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["da_lambda_mults"] = 1.0
    pst.pestpp_options["da_use_mda"] = True
    pst.pestpp_options["da_observation_cycle_table"] = "obs_cycle_tbl.csv"
    pst.pestpp_options["da_weight_cycle_table"] = "weight_cycle_tbl.csv"
    pst.pestpp_options["da_num_reals"] = 10
    pst.pestpp_options["da_localizer"] = "loc.mat"
    pst.pestpp_options["da_use_simulated_states"]  = False

    pst.write(os.path.join(t_d, "pest_seq.pst"), version=2)
    m_d = os.path.join(test_d, "master_da_double")
    pyemu.os_utils.start_workers(t_d, exe_path.replace("ies", "da"), "pest_seq.pst",
                                 num_workers=pst.pestpp_options["da_num_reals"], worker_root=test_d, port=port,
                                 master_dir=m_d, verbose=True)

    # first check the global pe's after cycles with and without assimilation
    pe0 = pd.read_csv(os.path.join(m_d,"pest_seq.global.0.pe.csv"),index_col=0)
    pe1 = pd.read_csv(os.path.join(m_d,"pest_seq.global.1.pe.csv"),index_col=0)
    d = np.abs(pe0.loc[:,final_state_pars].values - pe1.loc[:,init_state_pars].values).sum().sum()
    print(d)
    assert d < 1.0e-3

    pe2 = pd.read_csv(os.path.join(m_d, "pest_seq.global.2.pe.csv"), index_col=0)
    # these should be the same since we are not using simulated states (even tho we iterated during cycle 1)
    d = np.abs(pe2.loc[:, init_state_pars].values - pe1.loc[:, final_state_pars].values).sum().sum()
    assert d < 1.0e-3

    pe3 = pd.read_csv(os.path.join(m_d,"pest_seq.global.3.pe.csv"),index_col=0)
    # these should be diff since we assimilated during cycle 3
    d = np.abs(pe3.loc[:, final_state_pars].values - pe3.loc[:, init_state_pars].values).sum().sum()
    print(d)
    assert d > 1.0e-3

    pe4 = pd.read_csv(os.path.join(m_d,"pest_seq.global.4.pe.csv"),index_col=0)
    d = np.abs(pe3.loc[:, final_state_pars].values - pe4.loc[:, init_state_pars].values).sum().sum()
    print(d)
    assert d < 1.0e-3


    oe0 = pd.read_csv(os.path.join(m_d,"pest_seq.global.0.oe.csv"),index_col=0)
    # the final state pars should equal the obs states
    d = np.abs(oe0.loc[:,ostates].values - pe0.loc[:,final_state_pars].values).sum().sum()
    print(d)
    assert d < 1.0e-3
    # the init state pars should be diff than the obs states
    d = np.abs(oe0.loc[:, ostates].values - pe0.loc[:, init_state_pars].values).sum().sum()
    assert d > 1.0e-3

    oe00 = pd.read_csv(os.path.join(m_d, "pest_seq.0.0.obs.csv"), index_col=0)
    # these should be the same since we transferred from oe to pe after the initialization
    d = np.abs(oe00.loc[:, ostates].values - pe0.loc[:, final_state_pars].values).sum().sum()
    assert d < 1.0e-3


def seq_10par_xsec_double_state_test_2():
    """a test of estimating only final states"""
    test_d = "10par_xsec"
    t_d_org = os.path.join(test_d, "template")
    t_d = os.path.join(test_d, "template_double")
    if os.path.exists(t_d):
        shutil.rmtree(t_d)
    shutil.copytree(t_d_org,t_d)
    init_state_tpl_file = os.path.join(t_d,"strt_Layer_1.ref.tpl")
    lines = open(init_state_tpl_file,'r').readlines()
    final_state_tpl_file = init_state_tpl_file.replace("strt","final_states")
    with open(final_state_tpl_file,'w') as f:
        for line in lines:
            f.write(line.replace("_","final_"))

    pst = pyemu.Pst(os.path.join(t_d, "pest.pst"))
    df = pst.add_parameters(final_state_tpl_file,pst_path=".")
    df.sort_values(by="parnme",inplace=True)
    par = pst.parameter_data
    init_state_pars = df.parnme.str.replace("final_","_").tolist()[1:]
    final_state_pars = [i.replace("_","final_") for i in init_state_pars]
    par.loc[final_state_pars,"parval1"] = par.loc[init_state_pars,"parval1"].values
    par.loc[final_state_pars, "parubnd"] = par.loc[init_state_pars, "parubnd"].values
    par.loc[final_state_pars, "parlbnd"] = par.loc[init_state_pars, "parlbnd"].values
    par.loc[final_state_pars, "pargp"] = "final_states"



    par.loc[:, "cycle"] = -1
    par.loc[:,"state_par_link"] = np.nan
    par.loc[final_state_pars,"state_par_link"] = init_state_pars

    par.loc[par.parnme.apply(lambda x: x not in final_state_pars),"partrans"] = "fixed"
    strt_pars = par.loc[par.pargp == "strt", "parnme"].tolist()
    obs = pst.observation_data
    obs.loc[obs.obsnme.str.startswith("h01"), "weight"] = 1.0
    obs.loc[:, "state_par_link"] = ""
    obs.loc[obs.obgnme == "head1", "state_par_link"] = strt_pars
    ostates = obs.loc[obs.obgnme == "head1", "obsnme"]
    obs.loc[:, "cycle"] = -1

    pst.control_data.noptmax = 1

    pst.model_input_data.loc[:, "cycle"] = -1
    pst.model_output_data.loc[:, "cycle"] = -1

    def get_loc(pst):

        loc = pd.DataFrame(index=pst.nnz_obs_names, columns=pst.adj_par_names)
        loc.loc[:, :] = 0.0
        ocells = loc.index.map(lambda x: int(x.split('_')[1]))
        for pname in pst.adj_par_names:
            cstr = pname.split('_')[1]
            cint = int(cstr)
            # strt states can comm with all obs locs
            if "strt" in pname:
                dist = (ocells - cint).map(np.abs)
                loc_vals = 1 / (dist + 0.01)
                loc_vals = loc_vals.values
                loc_vals[loc_vals > 1.0] = 1.0
                loc.loc[:, pname] = loc_vals
            # static pars can only comm with obs in the same cell
            else:

                oname = [o for o in pst.nnz_obs_names if cstr in o.split('_')[1] == cstr][0]
                loc.loc[oname, pname] = 1.0
        return loc

    loc = get_loc(pst)
    pyemu.Matrix.from_dataframe(loc).to_ascii(os.path.join(t_d, "loc.mat"))

    # dont change these settings - they are hard coded below in the testing
    mx_cycle = 5
    cycles = np.arange(0, mx_cycle)
    odf = pd.DataFrame(index=cycles, columns=pst.nnz_obs_names)
    odf.loc[:, :] = obs.loc[pst.nnz_obs_names, "obsval"].values
    odf.T.to_csv(os.path.join(t_d, "obs_cycle_tbl.csv"))
    wdf = pd.DataFrame(index=cycles, columns=pst.nnz_obs_names)
    wdf.loc[:, :] = 0.0
    wdf.iloc[1, [3, 5]] = 1.0
    wdf.iloc[3, :] = 1.0
    wdf.T.to_csv(os.path.join(t_d, "weight_cycle_tbl.csv"))

    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["da_lambda_mults"] = 1.0
    pst.pestpp_options["da_use_mda"] = True
    pst.pestpp_options["da_observation_cycle_table"] = "obs_cycle_tbl.csv"
    pst.pestpp_options["da_weight_cycle_table"] = "weight_cycle_tbl.csv"
    pst.pestpp_options["da_num_reals"] = 10
    pst.pestpp_options["da_localizer"] = "loc.mat"
    pst.pestpp_options["da_use_simulated_states"]  = False

    pst.write(os.path.join(t_d, "pest_seq.pst"), version=2)
    m_d = os.path.join(test_d, "master_da_double")
    pyemu.os_utils.start_workers(t_d, exe_path.replace("ies", "da"), "pest_seq.pst",
                                 num_workers=pst.pestpp_options["da_num_reals"], worker_root=test_d, port=port,
                                 master_dir=m_d, verbose=True)

    # first check the global pe's after cycles with and without assimilation
    pe2 = pd.read_csv(os.path.join(m_d,"pest_seq.global.2.pe.csv"),index_col=0)
    pe1 = pd.read_csv(os.path.join(m_d,"pest_seq.global.1.pe.csv"),index_col=0)
    # cycle 1 final states should equal cycle 2 initial states since we are not using simulated states
    d = np.abs(pe1.loc[:,final_state_pars].values - pe2.loc[:,init_state_pars].values).sum().sum()
    print(d)
    assert d < 1.0e-3

    # cycle 1 final states should be diff than cycle 2 final states - even tho we didnt assimilate in
    # cycle 2, the model was advanced forward and the sim states were transfer to final par states after
    # initialization
    d = np.abs(pe2.loc[:, final_state_pars].values - pe1.loc[:, final_state_pars].values).sum().sum()
    assert d > 1.0e-3
    pe3 = pd.read_csv(os.path.join(m_d,"pest_seq.global.3.pe.csv"),index_col=0)

    # should be equal since not using simulated states
    d = np.abs(pe3.loc[:, init_state_pars].values - pe2.loc[:, final_state_pars].values).sum().sum()
    print(d)
    assert d < 1.0e-3

    pe4 = pd.read_csv(os.path.join(m_d,"pest_seq.global.4.pe.csv"),index_col=0)
    # this is an interesting check: by 5 cycles, the model should have
    # converged on the equilibrium values so init should equal final

    d = np.abs(pe4.loc[:, final_state_pars].values - pe4.loc[:, init_state_pars].values).sum().sum()
    print(d)
    assert d < 1.0e-3
    # and sim states should equal init/final states
    oe4 = pd.read_csv(os.path.join(m_d, "pest_seq.global.4.oe.csv"), index_col=0)
    # the final state pars should equal the obs states
    d = np.abs(oe4.loc[:, ostates].values - pe4.loc[:, final_state_pars].values).sum().sum()
    print(d)
    assert d < 1.0e-3



def seq_10par_xsec_double_state_test_3():
    """a test of carrying but not using a subset of final states"""
    test_d = "10par_xsec"
    t_d_org = os.path.join(test_d, "template")
    t_d = os.path.join(test_d, "template_double")
    if os.path.exists(t_d):
        shutil.rmtree(t_d)
    shutil.copytree(t_d_org,t_d)
    init_state_tpl_file = os.path.join(t_d,"strt_Layer_1.ref.tpl")
    lines = open(init_state_tpl_file,'r').readlines()
    final_state_tpl_file = init_state_tpl_file.replace("strt","final_states")
    with open(final_state_tpl_file,'w') as f:
        for line in lines:
            f.write(line.replace("_","final_"))

    pst = pyemu.Pst(os.path.join(t_d, "pest.pst"))
    df = pst.add_parameters(final_state_tpl_file,pst_path=".")
    df.sort_values(by="parnme",inplace=True)
    par = pst.parameter_data
    init_state_pars = df.parnme.str.replace("final_","_").tolist()[1:]
    final_state_pars = [i.replace("_","final_") for i in init_state_pars]
    par.loc[final_state_pars,"parval1"] = par.loc[init_state_pars,"parval1"].values
    par.loc[final_state_pars, "parubnd"] = par.loc[init_state_pars, "parubnd"].values
    par.loc[final_state_pars, "parlbnd"] = par.loc[init_state_pars, "parlbnd"].values
    par.loc[final_state_pars, "pargp"] = "final_states"


    par.loc[:, "cycle"] = -1
    par.loc[:,"state_par_link"] = np.nan
    par.loc[final_state_pars,"state_par_link"] = init_state_pars

    strt_pars = par.loc[par.pargp == "strt", "parnme"].tolist()
    obs = pst.observation_data
    obs.loc[obs.obsnme.str.startswith("h01"), "weight"] = 1.0
    obs.loc[:, "state_par_link"] = ""
    obs.loc[obs.obgnme == "head1", "state_par_link"] = strt_pars
    ostates = obs.loc[obs.obgnme == "head1", "obsnme"]
    obs.loc[:, "cycle"] = -1

    pst.control_data.noptmax = 1

    pst.model_input_data.loc[:, "cycle"] = -1
    pst.model_output_data.loc[:, "cycle"] = -1

    def get_loc(pst):

        loc = pd.DataFrame(index=pst.nnz_obs_names, columns=pst.adj_par_names)
        loc.loc[:, :] = 0.0
        ocells = loc.index.map(lambda x: int(x.split('_')[1]))
        for pname in pst.adj_par_names:
            cstr = pname.split('_')[1]
            cint = int(cstr)
            # strt states can comm with all obs locs
            if "strt" in pname:
                dist = (ocells - cint).map(np.abs)
                loc_vals = 1 / (dist + 0.01)
                loc_vals = loc_vals.values
                loc_vals[loc_vals > 1.0] = 1.0
                loc.loc[:, pname] = loc_vals
            # static pars can only comm with obs in the same cell
            else:

                oname = [o for o in pst.nnz_obs_names if cstr in o.split('_')[1] == cstr][0]
                loc.loc[oname, pname] = 1.0
        return loc

    loc = get_loc(pst)
    pyemu.Matrix.from_dataframe(loc).to_ascii(os.path.join(t_d, "loc.mat"))

    # dont change these settings - they are hard coded below in the testing
    mx_cycle = 5
    cycles = np.arange(0, mx_cycle)
    odf = pd.DataFrame(index=cycles, columns=pst.nnz_obs_names)
    odf.loc[:, :] = obs.loc[pst.nnz_obs_names, "obsval"].values
    odf.T.to_csv(os.path.join(t_d, "obs_cycle_tbl.csv"))
    wdf = pd.DataFrame(index=cycles, columns=pst.nnz_obs_names)
    wdf.loc[:, :] = 0.0
    wdf.iloc[1, [3, 5]] = 1.0
    wdf.iloc[3, :] = 1.0
    wdf.T.to_csv(os.path.join(t_d, "weight_cycle_tbl.csv"))

    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["da_lambda_mults"] = 1.0
    pst.pestpp_options["da_use_mda"] = True
    pst.pestpp_options["da_observation_cycle_table"] = "obs_cycle_tbl.csv"
    pst.pestpp_options["da_weight_cycle_table"] = "weight_cycle_tbl.csv"
    pst.pestpp_options["da_num_reals"] = 10
    pst.pestpp_options["da_localizer"] = "loc.mat"
    pst.pestpp_options["da_use_simulated_states"] = True

    pst.write(os.path.join(t_d, "pest_seq.pst"), version=2)
    m_d = os.path.join(test_d, "master_da_double")
    pyemu.os_utils.start_workers(t_d, exe_path.replace("ies", "da"), "pest_seq.pst",
                                 num_workers=pst.pestpp_options["da_num_reals"], worker_root=test_d, port=port,
                                 master_dir=m_d, verbose=True)

    # first check the global pe's after cycles with and without assimilation
    pe2 = pd.read_csv(os.path.join(m_d,"pest_seq.global.2.pe.csv"),index_col=0)
    pe1 = pd.read_csv(os.path.join(m_d,"pest_seq.global.1.pe.csv"),index_col=0)

    pst.pestpp_options["da_noptmax_schedule"] = "sched.dat"
    with open(os.path.join(t_d,"sched.dat"),'w') as f:
        f.write("1 3\n")
        f.write("3 2\n")

    pst.write(os.path.join(t_d, "pest_seq.pst"), version=2)
    m_d = os.path.join(test_d, "master_da_double_sched")
    pyemu.os_utils.start_workers(t_d, exe_path.replace("ies", "da"), "pest_seq.pst",
                                 num_workers=pst.pestpp_options["da_num_reals"], worker_root=test_d, port=port,
                                 master_dir=m_d, verbose=True)

def pump_test_2():
    test_d = "pump_test_2"
    t_d = os.path.join(test_d, "template")
    m_d = os.path.join(test_d,"master_da")
    pst = pyemu.Pst(os.path.join(t_d,"es_pmp.pst"))
    pst.pestpp_options["da_lambda_mults"] = 1.0
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.control_data.noptmax = 1
    pyemu.os_utils.start_workers(t_d,exe_path,"es_pmp.pst",num_workers=10,worker_root=test_d,
        master_dir=m_d)
    


def seq_10par_diff_state_cycle_test():

    test_d = "10par_xsec"
    t_d_org = os.path.join(test_d, "template")
    t_d = os.path.join(test_d, "template_double")
    if os.path.exists(t_d):
        shutil.rmtree(t_d)
    shutil.copytree(t_d_org,t_d)
    org_init_state_tpl_file = os.path.join(t_d,"strt_Layer_1.ref.tpl")
    lines = open(org_init_state_tpl_file,'r').readlines()
    pst = pyemu.Pst(os.path.join(t_d, "pest.pst"))

    tpl_df = pst.model_input_data
    tpl_df.index = tpl_df.index.map(lambda x: x.replace(".\\","").replace("./",""))
    for col in ["pest_file","model_file"]:
        tpl_df.loc[:,col] = tpl_df.loc[:,col].apply(lambda x: x.replace(".\\","").replace("./",""))


    pst.drop_parameters(org_init_state_tpl_file,pst_path=".")
    pst.parameter_data.loc[:, "cycle"] = -1
    par = pst.parameter_data

    par.loc[:,"state_par_link"] = np.nan
    pst.model_input_data.loc[:, "cycle"] = -1
    pst.model_output_data.loc[:, "cycle"] = -1
    mx_cycle = 5
    org_hds_ins_file = os.path.join(t_d,"10par_xsec.hds.ins")
    ins_lines = open(org_hds_ins_file,'r').readlines()
    pst.observation_data.loc[:,"cycle"] = -1
    pst.observation_data.loc[:,"state_par_link"] = ""

    pst.drop_observations(org_hds_ins_file,pst_path=".")
    pst.model_input_data.loc[:, "cycle"] = -1
    pst.model_output_data.loc[:, "cycle"] = -1

    for cycle in range(mx_cycle):
        init_state_tpl_file = org_init_state_tpl_file.replace("strt_","cycle{0}_strt".format(cycle))
        final_state_tpl_file = init_state_tpl_file.replace("strt","final_states")
        with open(init_state_tpl_file,'w') as f:
            for line in lines[:2]:
                f.write(line.replace("_","cycle{0}_".format(cycle)))
        with open(final_state_tpl_file,'w') as f:
            for line in lines[:2]:
                f.write(line.replace("_","cycle{0}_final_").format(cycle))

        df_init = pst.add_parameters(init_state_tpl_file,pst_path=".",in_file=os.path.join(t_d,"strt_Layer_1.ref"))
        df = pst.add_parameters(final_state_tpl_file,pst_path=".")
        df.sort_values(by="parnme",inplace=True)
        par = pst.parameter_data
        #print(par.parnme)
        init_state_pars = df.parnme.str.replace("final_","").tolist()
        final_state_pars = [i.replace("_","_final_") for i in init_state_pars]
        par.loc[init_state_pars,"parval1"] = 2
        par.loc[init_state_pars, "parlbnd"] = 0.1
        par.loc[init_state_pars, "parubnd"] = 10.0

        par.loc[final_state_pars,"parval1"] = par.loc[init_state_pars,"parval1"].values
        par.loc[final_state_pars, "parubnd"] = par.loc[init_state_pars, "parubnd"].values
        par.loc[final_state_pars, "parlbnd"] = par.loc[init_state_pars, "parlbnd"].values
        par.loc[final_state_pars, "pargp"] = "final_states"
        par.loc[init_state_pars,"cycle"] = cycle
        par.loc[final_state_pars,"cycle"] = cycle
        par.loc[final_state_pars,"state_par_link"] = init_state_pars
        pst.model_input_data.loc[os.path.join(".",os.path.split(init_state_tpl_file)[-1]),"cycle"] = cycle
        pst.model_input_data.loc[os.path.join(".",os.path.split(final_state_tpl_file)[-1]),"cycle"] = cycle

        cycle_ins_file = org_hds_ins_file.replace(".hds","_cycle{0}.hds".format(cycle))
        with open(cycle_ins_file,'w') as f:
            for line in ins_lines:
                if line.strip() == "":
                    break
                f.write(line.replace("!h0","!cycle{0}_h0").format(cycle))
        df_ins = pst.add_observations(cycle_ins_file,out_file=os.path.join(t_d,"10par_xsec.hds"),pst_path=".")
        pst.observation_data.loc[df_ins.obsnme,"cycle"] = cycle
        pst.observation_data.loc[df_ins.obsnme,"weight"] = 0
        pst.model_output_data.loc[os.path.join(".",os.path.split(cycle_ins_file)[-1]),"cycle"] = cycle
        #print(df_init.parnme.iloc[:df_ins.shape[0]].values)
        pst.observation_data.loc[df_ins.obsnme[:10].values,"state_par_link"] = df_init.parnme.iloc[:df_ins.shape[0]].sort_values().values
    #print(pst.observation_data.state_par_link)


    #print(pst.model_input_data)
    #print(pst.model_output_data)
    #return
    #strt_pars = par.loc[par.pargp == "strt", "parnme"].tolist()
    #obs = pst.observation_data
    #obs.loc[obs.obsnme.str.startswith("h01"), "weight"] = 1.0
    #obs.loc[:, "state_par_link"] = ""
    #obs.loc[obs.obgnme == "head1", "state_par_link"] = strt_pars
    #ostates = obs.loc[obs.obgnme == "head1", "obsnme"]
    #obs.loc[:, "cycle"] = -1

    pst.control_data.noptmax = 1

    
    # dont change these settings - they are hard coded below in the testing

    # cycles = np.arange(0, mx_cycle)
    # odf = pd.DataFrame(index=cycles, columns=pst.nnz_obs_names)
    # odf.loc[:, :] = obs.loc[pst.nnz_obs_names, "obsval"].values
    # odf.T.to_csv(os.path.join(t_d, "obs_cycle_tbl.csv"))
    # wdf = pd.DataFrame(index=cycles, columns=pst.nnz_obs_names)
    # wdf.loc[:, :] = 0.0
    # wdf.iloc[1, [3, 5]] = 1.0
    # wdf.iloc[3, :] = 1.0
    # wdf.T.to_csv(os.path.join(t_d, "weight_cycle_tbl.csv"))

    obs = pst.observation_data
    obs.loc[obs.obsnme.apply(lambda x: ("cycle1" in x or "cycle2" in x ) and ("h01_02" in x or "h01_04" in x)),"weight"] = 1.0

    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["da_lambda_mults"] = 1.0
    pst.pestpp_options["da_use_mda"] = True
    #pst.pestpp_options["da_observation_cycle_table"] = "obs_cycle_tbl.csv"
    #pst.pestpp_options["da_weight_cycle_table"] = "weight_cycle_tbl.csv"
    pst.pestpp_options["da_num_reals"] = 10
    pst.pestpp_options["da_use_simulated_states"] = True

    pst.write(os.path.join(t_d, "pest_seq.pst"), version=2)
    m_d = os.path.join(test_d, "master_da_double")
    pyemu.os_utils.start_workers(t_d, exe_path.replace("ies", "da"), "pest_seq.pst",
                                 num_workers=pst.pestpp_options["da_num_reals"], worker_root=test_d, port=port,
                                 master_dir=m_d, verbose=True)

    # first check the global pe's after cycles with and without assimilation
    #pe2 = pd.read_csv(os.path.join(m_d,"pest_seq.global.2.pe.csv"),index_col=0)
    #pe1 = pd.read_csv(os.path.join(m_d,"pest_seq.global.1.pe.csv"),index_col=0)



def seq_10par_xsec_ineq_test():
    """a test of estimating parameters and final states - initial states fixed"""
    test_d = "10par_xsec"
    t_d_org = os.path.join(test_d, "template")
    t_d = os.path.join(test_d, "template_double")
    if os.path.exists(t_d):
        shutil.rmtree(t_d)
    shutil.copytree(t_d_org,t_d)
    init_state_tpl_file = os.path.join(t_d,"strt_Layer_1.ref.tpl")
    lines = open(init_state_tpl_file,'r').readlines()
    final_state_tpl_file = init_state_tpl_file.replace("strt","final_states")
    with open(final_state_tpl_file,'w') as f:
        for line in lines:
            f.write(line.replace("_","final_"))

    pst = pyemu.Pst(os.path.join(t_d, "pest.pst"))
    df = pst.add_parameters(final_state_tpl_file,pst_path=".")
    df.sort_values(by="parnme",inplace=True)
    par = pst.parameter_data
    init_state_pars = df.parnme.str.replace("final_","_").tolist()[1:]
    final_state_pars = [i.replace("_","final_") for i in init_state_pars]
    par.loc[final_state_pars,"parval1"] = par.loc[init_state_pars,"parval1"].values
    par.loc[final_state_pars, "parubnd"] = par.loc[init_state_pars, "parubnd"].values
    par.loc[final_state_pars, "parlbnd"] = par.loc[init_state_pars, "parlbnd"].values
    par.loc[final_state_pars, "pargp"] = "final_states"


    par.loc[:, "cycle"] = -1
    par.loc[:,"state_par_link"] = np.nan
    par.loc[final_state_pars,"state_par_link"] = init_state_pars

    par.loc[par.parnme.str.contains("strt"), "partrans"] = "log"
    par.loc["cnhd_01", "partrans"] = "fixed"
    par.loc[init_state_pars,"partrans"] = "fixed"
    strt_pars = par.loc[par.pargp == "strt", "parnme"].tolist()
    obs = pst.observation_data
    obs.loc[obs.obsnme.str.startswith("h01"), "weight"] = 1.0
    obs.loc[:, "state_par_link"] = ""
    obs.loc[obs.obgnme == "head1", "state_par_link"] = strt_pars
    ostates = obs.loc[obs.obgnme == "head1","obsnme"]
    obs.loc[:, "cycle"] = -1

    pst.control_data.noptmax = 1

    pst.model_input_data.loc[:, "cycle"] = -1
    pst.model_output_data.loc[:, "cycle"] = -1

    def get_loc(pst):

        loc = pd.DataFrame(index=pst.nnz_obs_names, columns=pst.adj_par_names)
        loc.loc[:, :] = 0.0
        ocells = loc.index.map(lambda x: int(x.split('_')[1]))
        for pname in pst.adj_par_names:
            cstr = pname.split('_')[1]
            cint = int(cstr)
            # strt states can comm with all obs locs
            if "strt" in pname:
                dist = (ocells - cint).map(np.abs)
                loc_vals = 1 / (dist + 0.01)
                loc_vals = loc_vals.values
                loc_vals[loc_vals > 1.0] = 1.0
                loc.loc[:, pname] = loc_vals
            # static pars can only comm with obs in the same cell
            else:

                oname = [o for o in pst.nnz_obs_names if cstr in o.split('_')[1] == cstr][0]
                loc.loc[oname, pname] = 1.0
        return loc

    loc = get_loc(pst)
    pyemu.Matrix.from_dataframe(loc).to_ascii(os.path.join(t_d, "loc.mat"))

    # dont change these settings - they are hard coded below in the testing
    mx_cycle = 5
    cycles = np.arange(0, mx_cycle)
    odf = pd.DataFrame(index=cycles, columns=pst.nnz_obs_names)
    odf.loc[:, :] = obs.loc[pst.nnz_obs_names, "obsval"].values

    wdf = pd.DataFrame(index=cycles, columns=pst.nnz_obs_names)
    wdf.loc[:, :] = 0.0
    wdf.iloc[1, [3, 5]] = 1.0
    wdf.iloc[3, :] = 1.0
    wdf.iloc[4, :] = 1.0
    # this should cause a zero phi for cycle 1
    odf.iloc[1, [3, 5]] = 10000

    # this should cause a conflict for cycle 3
    odf.iloc[3, :] = -10000

    odf.T.to_csv(os.path.join(t_d, "obs_cycle_tbl.csv"))
    wdf.T.to_csv(os.path.join(t_d, "weight_cycle_tbl.csv"))
    obs = pst.observation_data
    obs.loc[pst.nnz_obs_names,"obgnme"] = "less_than"
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["da_lambda_mults"] = 1.0
    pst.pestpp_options["da_use_mda"] = True
    pst.pestpp_options["da_observation_cycle_table"] = "obs_cycle_tbl.csv"
    pst.pestpp_options["da_weight_cycle_table"] = "weight_cycle_tbl.csv"
    pst.pestpp_options["da_num_reals"] = 10
    pst.pestpp_options["da_localizer"] = "loc.mat"
    pst.pestpp_options["da_use_simulated_states"] = False
    pst.pestpp_options["da_drop_conflicts"] = True

    pst.write(os.path.join(t_d, "pest_seq.pst"), version=2)
    m_d = os.path.join(test_d, "master_da_double_ineq")
    pyemu.os_utils.start_workers(t_d, exe_path.replace("ies", "da"), "pest_seq.pst",
                                 num_workers=pst.pestpp_options["da_num_reals"], worker_root=test_d, port=port,
                                 master_dir=m_d, verbose=True)

    phi_df = pd.read_csv(os.path.join(m_d,"pest_seq.global.phi.actual.csv"))
    cycle1_sum = phi_df.loc[phi_df.cycle==1,:].iloc[:,1:].sum().sum()
    print(cycle1_sum)
    assert cycle1_sum == 0.0

    cycle3_sum = phi_df.loc[phi_df.cycle == 3, :].iloc[:, 1:].sum().sum()
    print(cycle3_sum)
    assert cycle3_sum == 0.0

    cycle4_iter0_sum = phi_df.loc[phi_df.apply(lambda x: x.cycle == 4 and x.iteration==0,axis=1), :].iloc[:, 1:].sum().sum()
    cycle4_iter1_sum = phi_df.loc[phi_df.apply(lambda x: x.cycle == 4 and x.iteration==1,axis=1), :].iloc[:, 1:].sum().sum()
    assert cycle4_iter0_sum > 0
    assert cycle4_iter1_sum > 0
    assert cycle4_iter1_sum < cycle4_iter0_sum

    rec_lines = open(os.path.join(m_d,"pest_seq.rec"),'r').readlines()
    found_1,found_2 = False,False
    found_3 = False

    for line in rec_lines:
        #print(line.strip())
        if "initial actual phi mean too low but only inequality obs are being used" in line:
            found_1 = True
        if "current mean actual and/or measurement phi is too low for solution" in line:
            found_2 = True
        if "all non-zero weighted observations in conflict state, continuing to next cycle" in line:
            found_3 = True
        if found_1 and found_2 and found_3:
            break
    assert found_1 is True
    assert found_2 is True
    assert found_3 is True
    rec_lines = None



    # first check the global pe's after cycles with and without assimilation
    pe0 = pd.read_csv(os.path.join(m_d,"pest_seq.global.0.pe.csv"),index_col=0)
    pe1 = pd.read_csv(os.path.join(m_d,"pest_seq.global.1.pe.csv"),index_col=0)
    d = np.abs(pe0.loc[:,final_state_pars].values - pe1.loc[:,init_state_pars].values).sum().sum()
    print(d)
    assert d < 1.0e-3

    pe2 = pd.read_csv(os.path.join(m_d, "pest_seq.global.2.pe.csv"), index_col=0)
    # these should be the same since we are not using simulated states (even tho we iterated during cycle 1)
    d = np.abs(pe2.loc[:, init_state_pars].values - pe1.loc[:, final_state_pars].values).sum().sum()
    assert d < 1.0e-3

    pe3 = pd.read_csv(os.path.join(m_d,"pest_seq.global.3.pe.csv"),index_col=0)
    # these should be diff since we assimilated during cycle 3
    d = np.abs(pe3.loc[:, final_state_pars].values - pe3.loc[:, init_state_pars].values).sum().sum()
    print(d)
    assert d > 1.0e-3

    pe4 = pd.read_csv(os.path.join(m_d,"pest_seq.global.4.pe.csv"),index_col=0)
    d = np.abs(pe3.loc[:, final_state_pars].values - pe4.loc[:, init_state_pars].values).sum().sum()
    print(d)
    assert d < 1.0e-3


    oe0 = pd.read_csv(os.path.join(m_d,"pest_seq.global.0.oe.csv"),index_col=0)
    # the final state pars should equal the obs states
    d = np.abs(oe0.loc[:,ostates].values - pe0.loc[:,final_state_pars].values).sum().sum()
    print(d)
    assert d < 1.0e-3
    # the init state pars should be diff than the obs states
    d = np.abs(oe0.loc[:, ostates].values - pe0.loc[:, init_state_pars].values).sum().sum()
    assert d > 1.0e-3

    oe00 = pd.read_csv(os.path.join(m_d, "pest_seq.0.0.obs.csv"), index_col=0)
    # these should be the same since we transferred from oe to pe after the initialization
    d = np.abs(oe00.loc[:, ostates].values - pe0.loc[:, final_state_pars].values).sum().sum()
    assert d < 1.0e-3


def seq_10par_xsec_double_state_test_with_fail():
    """a test of carrying but not using a subset of final states"""
    test_d = "10par_xsec"
    t_d_org = os.path.join(test_d, "template")
    t_d = os.path.join(test_d, "template_double")
    if os.path.exists(t_d):
        shutil.rmtree(t_d)
    shutil.copytree(t_d_org,t_d)
    init_state_tpl_file = os.path.join(t_d,"strt_Layer_1.ref.tpl")
    lines = open(init_state_tpl_file,'r').readlines()
    final_state_tpl_file = init_state_tpl_file.replace("strt","final_states")
    with open(final_state_tpl_file,'w') as f:
        for line in lines:
            f.write(line.replace("_","final_"))

    pst = pyemu.Pst(os.path.join(t_d, "pest.pst"))
    df = pst.add_parameters(final_state_tpl_file,pst_path=".")
    df.sort_values(by="parnme",inplace=True)
    par = pst.parameter_data
    init_state_pars = df.parnme.str.replace("final_","_").tolist()[1:]
    final_state_pars = [i.replace("_","final_") for i in init_state_pars]
    par.loc[final_state_pars,"parval1"] = par.loc[init_state_pars,"parval1"].values
    par.loc[final_state_pars, "parubnd"] = par.loc[init_state_pars, "parubnd"].values
    par.loc[final_state_pars, "parlbnd"] = par.loc[init_state_pars, "parlbnd"].values
    par.loc[final_state_pars, "pargp"] = "final_states"


    par.loc[:, "cycle"] = -1
    par.loc[:,"state_par_link"] = np.nan
    par.loc[final_state_pars,"state_par_link"] = init_state_pars

    strt_pars = par.loc[par.pargp == "strt", "parnme"].tolist()
    obs = pst.observation_data
    obs.loc[obs.obsnme.str.startswith("h01"), "weight"] = 1.0
    obs.loc[:, "state_par_link"] = ""
    obs.loc[obs.obgnme == "head1", "state_par_link"] = strt_pars
    ostates = obs.loc[obs.obgnme == "head1", "obsnme"]
    obs.loc[:, "cycle"] = -1

    pst.control_data.noptmax = 1

    pst.model_input_data.loc[:, "cycle"] = -1
    pst.model_output_data.loc[:, "cycle"] = -1

    def get_loc(pst):

        loc = pd.DataFrame(index=pst.nnz_obs_names, columns=pst.adj_par_names)
        loc.loc[:, :] = 0.0
        ocells = loc.index.map(lambda x: int(x.split('_')[1]))
        for pname in pst.adj_par_names:
            cstr = pname.split('_')[1]
            cint = int(cstr)
            # strt states can comm with all obs locs
            if "strt" in pname:
                dist = (ocells - cint).map(np.abs)
                loc_vals = 1 / (dist + 0.01)
                loc_vals = loc_vals.values
                loc_vals[loc_vals > 1.0] = 1.0
                loc.loc[:, pname] = loc_vals
            # static pars can only comm with obs in the same cell
            else:

                oname = [o for o in pst.nnz_obs_names if cstr in o.split('_')[1] == cstr][0]
                loc.loc[oname, pname] = 1.0
        return loc

    loc = get_loc(pst)
    pyemu.Matrix.from_dataframe(loc).to_ascii(os.path.join(t_d, "loc.mat"))

    # dont change these settings - they are hard coded below in the testing
    mx_cycle = 5
    cycles = np.arange(0, mx_cycle)
    odf = pd.DataFrame(index=cycles, columns=pst.nnz_obs_names)
    odf.loc[:, :] = obs.loc[pst.nnz_obs_names, "obsval"].values
    odf.T.to_csv(os.path.join(t_d, "obs_cycle_tbl.csv"))
    wdf = pd.DataFrame(index=cycles, columns=pst.nnz_obs_names)
    wdf.loc[:, :] = 0.0
    wdf.iloc[1, [3, 5]] = 1.0
    wdf.iloc[3, :] = 1.0
    wdf.T.to_csv(os.path.join(t_d, "weight_cycle_tbl.csv"))

    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["da_lambda_mults"] = 1.0
    pst.pestpp_options["da_use_mda"] = True
    pst.pestpp_options["da_observation_cycle_table"] = "obs_cycle_tbl.csv"
    pst.pestpp_options["da_weight_cycle_table"] = "weight_cycle_tbl.csv"
    pst.pestpp_options["da_num_reals"] = 10
    pst.pestpp_options["da_localizer"] = "loc.mat"
    pst.pestpp_options["da_use_simulated_states"] = True
    pst.pestpp_options["ies_debug_fail_subset"] = True
    pst.pestpp_options["ies_debug_fail_remainder"] = True
    pst.pestpp_options["ies_update_by_reals"] = False

    pst.write(os.path.join(t_d, "pest_seq.pst"), version=2)
    m_d = os.path.join(test_d, "master_da_double_fail")
    pyemu.os_utils.start_workers(t_d, exe_path.replace("ies", "da"), "pest_seq.pst",
                                 num_workers=pst.pestpp_options["da_num_reals"], worker_root=test_d, port=port,
                                 master_dir=m_d, verbose=True)

    assert os.path.exists(os.path.join(m_d,"pest_seq.global.4.pe.csv"))
    assert os.path.exists(os.path.join(m_d, "pest_seq.global.4.oe.csv"))
    assert os.path.exists(os.path.join(m_d, "pest_seq.1.0.base.rei"))
    assert not os.path.exists(os.path.join(m_d, "pest_seq.1.1.base.rei"))
    assert not os.path.exists(os.path.join(m_d, "pest_seq.2.0.base.rei"))
    assert not os.path.exists(os.path.join(m_d, "pest_seq.3.1.base.rei"))
    assert not os.path.exists(os.path.join(m_d, "pest_seq.4.0.base.rei"))

    files = [f for f in os.listdir(m_d) if f.lower().endswith(".csv") and ("pest_seq.4" in f or "global.4" in f)]
    for f in files:

        df = pd.read_csv(os.path.join(m_d,f),index_col=0)
        print(f,set(list(df.index.values)))
        assert "base" not in set(list(df.index.values))
        assert df.shape[0] == 8


if __name__ == "__main__":
    
    seq_10par_xsec_double_state_test_with_fail()
    #shutil.copy2(os.path.join("..","exe","windows","x64","Debug","pestpp-da.exe"),os.path.join("..","bin","pestpp-da.exe"))
    #seq_10par_cycle_parse_test()
    #seq_10par_xsec_hotstart_test()
    #seq_10par_diff_obspar_cycle_test()
    #da_mf6_freyberg_test_1()
    #da_mf6_freyberg_test_2()
    #da_mf6_freyberg_smoother_test()
    #da_prep_4_mf6_freyberg_seq_tbl()
    #da_build_mf6_freyberg_seq_localizer_tbl()
    #da_build_mf6_freyberg_seq_localizer()
    #da_prep_4_mf6_freyberg_seq(sync_state_names=False)
    #da_mf6_freyberg_test_3()
    #seq_10par_xsec_state_est_test()
    #seq_10par_xsec_fixed_test()
    #compare_mf6_freyberg()
    #plot_compare("glm",noptmax=3)
    #plot_compare("mda",noptmax=3)
    #da_pareto_demo()
    #plot_da_pareto_demo()
    #seq_10par_xsec_double_state_test()
    #seq_10par_xsec_double_state_test_2()
    #seq_10par_xsec_state_est_test()
    #seq_10par_xsec_double_state_test_3()
    #pump_test_2()
    #seq_10par_diff_state_cycle_test()
    #seq_10par_diff_obspar_cycle_test()
    #seq_10par_xsec_ineq_test()
    #seq_10par_xsec_double_state_test_with_fail()
