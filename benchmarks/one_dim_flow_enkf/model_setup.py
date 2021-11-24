import os, sys
import shutil
import numpy as np
import pandas as pd
import flopy
import matplotlib.pyplot as plt
import pyemu
from matplotlib.pyplot import cm

sys.path.append(r"..")
import obs_utils, param_utils


def mf3d_to_arrays(arr3, fname):
    array = arr3
    nlayer = arr3.shape[0]
    basename = os.path.basename(fname)
    dirname = os.path.dirname(fname)
    prefix, ext = os.path.splitext(basename)
    for k in range(nlayer):
        f = prefix + "_" + str(k + 1) + ext
        np.savetxt(os.path.join(dirname, f), array[k])
    return array


def generate_1d_model():
    modelname = "flow_1d"
    model_ws = r"model_dataset"

    if model_ws in os.listdir(os.getcwd()):
        shutil.rmtree(os.path.join(os.getcwd(), model_ws))
    os.mkdir(model_ws)

    exe = r"..\bin\win\mfnwt.exe"
    shutil.copy(src=exe, dst=os.path.join(model_ws, os.path.basename(exe)))

    mf = flopy.modflow.Modflow(modelname, model_ws=model_ws,
                               exe_name=os.path.abspath(exe),
                               version='mfnwt')
    # --- Dis file
    Lx = 2000.0;
    Ly = 1.0
    ztop = 0.0;
    zbot = -50.0
    nlay = 1;
    nrow = 1;
    ncol = 100
    delr = Lx / ncol
    delc = Ly / nrow
    delv = (ztop - zbot) / nlay
    botm = np.linspace(ztop, zbot, nlay + 1)

    nper = 24  # two years
    perlen = nper * [30]
    nstp = nper * [30]
    steady = [True] + (nper - 1) * [False]

    dis = flopy.modflow.ModflowDis(
        mf, nlay, nrow, ncol, delr=delr, delc=delc, top=ztop, botm=botm[1:],
        nper=nper, perlen=perlen, nstp=nstp, steady=steady
    )

    df_dis = pd.DataFrame(columns=['parname', 'parval'])
    df_dis['parname'] = ['start_sp', 'end_sp']
    df_dis['parval'] = [0, 23]
    df_dis.to_csv(os.path.join(model_ws, "input_dis.csv"), index=False)

    temporal_param = pd.DataFrame(columns=[])
    temporal_param['stress_period'] = range(mf.dis.nper)
    temporal_param['perlen'] = mf.dis.perlen.array
    temporal_param['nstp'] = mf.dis.nstp.array

    # --- bas file
    ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
    strt = np.ones((nlay, nrow, ncol), dtype=np.float32)
    strt[:, :, :] = -10.0
    bas = flopy.modflow.ModflowBas(mf, ibound=ibound, strt=strt)

    # --- upw
    laytyp = np.ones(nlay)
    laywet = np.zeros(nlay)
    hk = np.zeros_like(bas.strt.array) + 1.6
    sy = 0.18
    ss = 1e-06
    flopy.modflow.mfupw.ModflowUpw(mf, laytyp=laytyp, layavg=0, chani=1.0, layvka=0,
                                   laywet=laywet, hdry=-1e+30, iphdry=0, hk=hk, hani=1.0,
                                   vka=hk, ss=ss, sy=sy, vkcb=0.0, noparcheck=False, ipakcb=55)

    # --- GHB
    ghb_stress_per = {}
    w = 2.0 * np.pi * 20
    t = np.linspace(0, mf.dis.nper)
    sinPart = 2 * np.sin(w * t)
    for sp in range(mf.dis.nper):
        ghb_data = []  # [lay, row, col,head,cond]
        ghb_data.append([0, 0, 0, -10, 0.1])
        ghb_data.append([0, 0, ncol - 1, -20, 0.1])
        ghb_stress_per[sp] = ghb_data
        temporal_param.loc[temporal_param['stress_period'] == sp, 'GHB_UP'] = -10
        temporal_param.loc[temporal_param['stress_period'] == sp, 'GHB_DN'] = -20

    ghbs = flopy.modflow.mfghb.ModflowGhb(mf, ipakcb=55, stress_period_data=ghb_stress_per,
                                          dtype=None, no_print=False,
                                          options=None, extension='ghb')

    # --- well
    wel_data = dict()
    for sp in range(mf.dis.nper):
        flow = -3.5
        if sp == 0:
            flow = 0
        wel_data[sp] = [[0, 0, 50, flow]]
        wel = flopy.modflow.mfwel.ModflowWel(mf, stress_period_data=wel_data)
        temporal_param.loc[temporal_param['stress_period'] == sp, 'Qw'] = flow

    # rech_data
    rech_data = dict()
    for sp in range(mf.dis.nper):
        # if sp == 0:
        #     continue
        rech_data[sp] = np.zeros_like(rech_data) + 0.0002 + 1e-3 * sinPart[sp]
        rch = flopy.modflow.ModflowRch(mf, nrchop=3, rech=rech_data)
        temporal_param.loc[temporal_param['stress_period'] == sp, 'Rch'] = 0.0002 + 1e-3 * sinPart[sp]

    spd = dict()
    for isp, sp in enumerate(mf.dis.nstp.array):
        spd[(isp, sp - 1)] = ["print head", "print budget", "save head", "save budget"]
    oc = flopy.modflow.ModflowOc(mf, stress_period_data=spd, compact=True)

    nwt = flopy.modflow.mfnwt.ModflowNwt.load(r".\misc_files\solver_options.nwt", mf)

    mf.write_input()

    success, buff = mf.run_model()
    if not success:
        raise Exception("MODFLOW did not terminate normally.")

    hds = flopy.utils.HeadFile(os.path.join(model_ws, modelname + ".hds"))
    totims = hds.get_times()
    color = cm.rainbow(np.linspace(0, 1, len(totims)))
    fig1 = plt.figure()
    for ic, totim in enumerate(totims):
        wl = hds.get_data(totim=totim)
        plt.plot(wl[0][0], c=color[ic])
    plt.xlabel("Distance")
    plt.ylabel("Water level")

    fig2 = plt.figure()
    cbc = flopy.utils.CellBudgetFile(os.path.join(model_ws, modelname + ".cbc"))
    plt.plot(cbc.get_ts(idx=(0, 0, ncol - 1), text=' HEAD DEP BOUNDS')[:, 1])
    plt.xlabel("Time")
    plt.ylabel("Downstream Flow at GHB")

    # generate initial head to be steady-state for k = 1.0
    wl = hds.get_data(totim=30.0)
    np.savetxt(os.path.join(model_ws, "iheads.dat"), wl[0][0])

    # let use log of k
    np.savetxt(os.path.join(model_ws, "hk.dat"), np.log10(mf.upw.hk.array[0][0]))
    #mf3d_to_arrays(wl, os.path.join(model_ws, "iheads.dat"))
    #mf3d_to_arrays(mf.upw.hk.array, os.path.join(model_ws, "hk.dat"))
    # parnme = [];
    # parval1 = [];
    # ii, jj, kk = [], [], []
    #
    # for k in range(nlay):
    #     for i in range(nrow):
    #         for j in range(ncol):
    #             nm = "h_{}_{}_{}".format(k, i, j)
    #             parnme.append(nm)
    #             parval1.append(wl[k, i, j])
    #             ii.append(i)
    #             jj.append(j)
    #             kk.append(k)
    #
    # df = param_utils.add_param(df=df, parnme=parnme,
    #                            parval1=parval1,
    #                            gpname=['ihead'], i=ii,
    #                            j=jj, k=kk)
    # parnme = [];
    # parval1 = [];
    # ii, jj, kk = [], [], []
    #
    # for k in range(nlay):
    #     for i in range(nrow):
    #         for j in range(ncol):
    #             nm = "k_{}_{}_{}".format(k, i, j)
    #             parnme.append(nm)
    #             parval1.append(mf.upw.hk.array[k, i, j])
    #             ii.append(i)
    #             jj.append(j)
    #             kk.append(k)
    #
    # df = param_utils.add_param(df=df, parnme=parnme,
    #                            parval1=parval1,
    #                            gpname=['k'], i=ii,
    #                            j=jj, k=kk)

    # df.to_csv(os.path.join(model_ws, "input_param.csv"), index=False)
    temporal_param.to_csv(os.path.join(model_ws, "temporal_param.csv"), index=False)
    return [fig1, fig2]


def forward_run():
    modelname = "flow_1d"
    model_ws = r"model_dataset"

    # always start by removing old output files
    try:
        os.remove(os.path.join(model_ws, "output_sim.csv"))
    except:
        pass

    input_dis = pd.read_csv(os.path.join(model_ws, "input_dis.csv"))
    temporal_df = pd.read_csv(os.path.join(model_ws, "temporal_param.csv"))

    mf = flopy.modflow.Modflow.load(f=modelname + ".nam", model_ws=model_ws)

    # ------ update dis
    start_sp = input_dis.loc[input_dis['parname'] ==  'start_sp', 'parval'].values[0]
    end_sp = input_dis.loc[input_dis['parname'] ==  'end_sp', 'parval'].values[0]
    nlay = mf.dis.nlay
    nrow = mf.dis.nrow
    ncol = mf.dis.ncol
    delr = mf.dis.delr.array
    delc = mf.dis.delc.array
    ztop = mf.dis.top.array
    botm = mf.dis.botm.array

    time_mask = (temporal_df['stress_period'] >= start_sp) & \
                (temporal_df['stress_period'] <= end_sp)
    curr_tempora_df = temporal_df[time_mask]
    curr_sps = np.sort(curr_tempora_df['stress_period'].unique())
    nper = len(curr_sps)

    perlen = []
    nstp = []
    steady = []
    for sp in curr_sps:
        perlen.append(curr_tempora_df.loc[curr_tempora_df['stress_period'] == sp, 'perlen'].values[0])
        nstp.append(curr_tempora_df.loc[curr_tempora_df['stress_period'] == sp, 'nstp'].values[0])
        steady.append(False)
    mf.remove_package('DIS')
    dis = flopy.modflow.ModflowDis(
        mf, nlay, nrow, ncol, delr=delr, delc=delc, top=ztop, botm=botm,
        nper=nper, perlen=perlen, nstp=nstp, steady=steady)

    # ------ update bas
    initial_head = np.loadtxt(os.path.join(model_ws, r'iheads.dat'))
    strt = np.zeros_like(mf.bas6.strt.array)
    strt[0,0,:] = initial_head
    mf.bas6.strt = strt

    # ------ update upw
    hk2d = np.loadtxt(os.path.join(model_ws, r'hk.dat'))
    hk2d = np.power(10.0, hk2d) # k is in log space
    hk3d = np.zeros_like(mf.upw.hk.array)
    hk3d[0, 0, :] = hk2d
    mf.upw.hk = hk3d

    # ------ update ghb
    ghb_data = {}
    for isp, sp in enumerate(curr_sps):
        ghbb = []
        upval = curr_tempora_df.loc[curr_tempora_df['stress_period'] == sp, 'GHB_UP'].values[0]
        dnval = curr_tempora_df.loc[curr_tempora_df['stress_period'] == sp, 'GHB_DN'].values[0]
        ghbb.append([0, 0, 0, upval, 0.1])
        ghbb.append([0, 0, ncol - 1, dnval, 0.1])
        ghb_data[isp] = ghbb
    mf.ghb.stress_period_data = ghb_data

    # ------ update wel
    wel_data = dict()
    for isp, sp in enumerate(curr_sps):
        welQ = curr_tempora_df.loc[curr_tempora_df['stress_period'] == sp, 'Qw'].values[0]
        wel_data[isp] = [[0, 0, 50, welQ]]
    mf.wel.stress_period_data = wel_data

    # ------ update rech
    rech_data = dict()
    for isp, sp in enumerate(curr_sps):
        rch_val = curr_tempora_df.loc[curr_tempora_df['stress_period'] == sp, 'Rch'].values[0]
        rech_data[isp] = np.zeros_like(rech_data) + rch_val
    mf.remove_package('RCH')
    rch = flopy.modflow.ModflowRch(mf, nrchop=3, rech=rech_data)

    # ------ update OC
    spd = dict()
    for isp, sp in enumerate(mf.dis.nstp.array):
        spd[(isp, sp - 1)] = ["print head", "print budget", "save head", "save budget"]
    mf.remove_package('OC')
    oc = flopy.modflow.ModflowOc(mf, stress_period_data=spd, compact=True)

    mf.write_input()
    exe = os.path.join(model_ws, r"mfnwt.exe")
    mf.exe_name = os.path.abspath(exe)
    success, buff = mf.run_model()
    if not success:
        raise Exception("MODFLOW did not terminate normally.")

    hds = flopy.utils.HeadFile(os.path.join(model_ws, modelname + ".hds"))
    cbc = flopy.utils.CellBudgetFile(os.path.join(model_ws, modelname + ".cbc"))
    totims = hds.get_times()
    wl = hds.get_data(totim=np.max(totims))

    dyn_stat = []
    stat_value = []
    for k in range(nlay):
        for i in range(nrow):
            for j in range(ncol):
                nm = "hf_{}_{}_{}".format(k, i, j)
                dyn_stat.append(nm)
                stat_value.append(wl[k, i, j])

    sim_obs = pd.DataFrame(columns=obs_utils.get_header())
    sim_obs = obs_utils.add_obs(sim_obs, obsnams=dyn_stat, simval=stat_value,
                                obsval=np.zeros_like(stat_value),
                                obsgnme='simH')

    # Assume we have an obs at the end of each stress period
    obs_idx = list(range(10, ncol - 10, 10))
    ho_base_name = ["ho_{}".format(i) for i in obs_idx]
    h_obs = []
    h_nms = []
    for totim in totims:
        wl = hds.get_data(totim=np.max(totim))
        h_obs = h_obs + wl[0][0][obs_idx].tolist()
        nm = ["{}_{}".format(i, int(totim)) for i in ho_base_name]
        h_nms = h_nms + nm

    sim_obs = obs_utils.add_obs(sim_obs, obsnams=h_nms, simval=h_obs,
                                obsval=np.zeros_like(h_obs),
                                obsgnme='simHO')
    q_ghb = cbc.get_ts(idx=(0, 0, ncol - 1), text=' HEAD DEP BOUNDS')[:, 1]
    ghb_nm = ["qghb_{}".format(int(i)) for i in totims]
    sim_obs = obs_utils.add_obs(sim_obs, obsnams=ghb_nm, simval=q_ghb,
                                obsval=np.zeros_like(q_ghb),
                                obsgnme='ghbQ')

    sim_obs.to_csv(os.path.join(model_ws, "output_sim.csv"), index=False)
    np.savetxt(os.path.join(model_ws,"output_sim.dat"), sim_obs['simval'].values)


def forward_run3():
    modelname = "flow_1d"
    model_ws = r"model_dataset"

    # always start by removing old output files
    try:
        os.remove(os.path.join(model_ws, "output_sim.csv"))
    except:
        pass

    input_dis = pd.read_csv(os.path.join(model_ws, "input_dis.csv"))
    temporal_df = pd.read_csv(os.path.join(model_ws, "temporal_param.csv"))

    mf = flopy.modflow.Modflow.load(f=modelname + ".nam", model_ws=model_ws)

    # ------ update dis
    start_sp = input_dis.loc[input_dis['parname'] ==  'start_sp', 'parval'].values[0]
    end_sp = input_dis.loc[input_dis['parname'] ==  'end_sp', 'parval'].values[0]
    nlay = mf.dis.nlay
    nrow = mf.dis.nrow
    ncol = mf.dis.ncol
    delr = mf.dis.delr.array
    delc = mf.dis.delc.array
    ztop = mf.dis.top.array
    botm = mf.dis.botm.array

    time_mask = (temporal_df['stress_period'] >= start_sp) & \
                (temporal_df['stress_period'] <= end_sp)
    curr_tempora_df = temporal_df[time_mask]
    curr_sps = np.sort(curr_tempora_df['stress_period'].unique())
    nper = len(curr_sps)

    perlen = []
    nstp = []
    steady = []
    for sp in curr_sps:
        perlen.append(curr_tempora_df.loc[curr_tempora_df['stress_period'] == sp, 'perlen'].values[0])
        nstp.append(curr_tempora_df.loc[curr_tempora_df['stress_period'] == sp, 'nstp'].values[0])
        steady.append(False)
    mf.remove_package('DIS')
    dis = flopy.modflow.ModflowDis(
        mf, nlay, nrow, ncol, delr=delr, delc=delc, top=ztop, botm=botm,
        nper=nper, perlen=perlen, nstp=nstp, steady=steady)

    # ------ update bas
    initial_head = np.loadtxt(os.path.join(model_ws, r'iheads_1.dat'))
    strt = np.zeros_like(mf.bas6.strt.array)
    strt[0,0,:] = initial_head
    mf.bas6.strt = strt

    # ------ update upw
    hk2d = np.loadtxt(os.path.join(model_ws, r'hk_1.dat'))

    hk3d = np.zeros_like(mf.upw.hk.array)
    hk3d[0, 0, :] = hk2d
    mf.upw.hk = hk3d

    # ------ update ghb
    ghb_data = {}
    for isp, sp in enumerate(curr_sps):
        ghbb = []
        upval = curr_tempora_df.loc[curr_tempora_df['stress_period'] == sp, 'GHB_UP'].values[0]
        dnval = curr_tempora_df.loc[curr_tempora_df['stress_period'] == sp, 'GHB_DN'].values[0]
        ghbb.append([0, 0, 0, upval, 0.1])
        ghbb.append([0, 0, ncol - 1, dnval, 0.1])
        ghb_data[isp] = ghbb
    mf.ghb.stress_period_data = ghb_data

    # ------ update wel
    wel_data = dict()
    for isp, sp in enumerate(curr_sps):
        welQ = curr_tempora_df.loc[curr_tempora_df['stress_period'] == sp, 'Qw'].values[0]
        wel_data[isp] = [[0, 0, 50, welQ]]
    mf.wel.stress_period_data = wel_data

    # ------ update rech
    rech_data = dict()
    for isp, sp in enumerate(curr_sps):
        rch_val = curr_tempora_df.loc[curr_tempora_df['stress_period'] == sp, 'Rch'].values[0]
        rech_data[isp] = np.zeros_like(rech_data) + rch_val
    mf.remove_package('RCH')
    rch = flopy.modflow.ModflowRch(mf, nrchop=3, rech=rech_data)

    # ------ update OC
    spd = dict()
    for isp, sp in enumerate(mf.dis.nstp.array):
        spd[(isp, sp - 1)] = ["print head", "print budget", "save head", "save budget"]
    mf.remove_package('OC')
    oc = flopy.modflow.ModflowOc(mf, stress_period_data=spd, compact=True)

    mf.write_input()
    exe = os.path.join(model_ws, r"mfnwt.exe")
    mf.exe_name = os.path.abspath(exe)
    success, buff = mf.run_model()
    if not success:
        raise Exception("MODFLOW did not terminate normally.")

    hds = flopy.utils.HeadFile(os.path.join(model_ws, modelname + ".hds"))
    cbc = flopy.utils.CellBudgetFile(os.path.join(model_ws, modelname + ".cbc"))
    totims = hds.get_times()
    wl = hds.get_data(totim=np.max(totims))

    dyn_stat = []
    stat_value = []
    for k in range(nlay):
        for i in range(nrow):
            for j in range(ncol):
                nm = "hf_{}_{}_{}".format(k, i, j)
                dyn_stat.append(nm)
                stat_value.append(wl[k, i, j])

    sim_obs = pd.DataFrame(columns=obs_utils.get_header())
    sim_obs = obs_utils.add_obs(sim_obs, obsnams=dyn_stat, simval=stat_value,
                                obsval=np.zeros_like(stat_value),
                                obsgnme='simH')

    # Assume we have an obs at the end of each stress period
    obs_idx = list(range(10, ncol - 10, 10))
    ho_base_name = ["ho_{}".format(i) for i in obs_idx]
    h_obs = []
    h_nms = []
    for totim in totims:
        wl = hds.get_data(totim=np.max(totim))
        h_obs = h_obs + wl[0][0][obs_idx].tolist()
        nm = ["{}_{}".format(i, int(totim)) for i in ho_base_name]
        h_nms = h_nms + nm

    sim_obs = obs_utils.add_obs(sim_obs, obsnams=h_nms, simval=h_obs,
                                obsval=np.zeros_like(h_obs),
                                obsgnme='simHO')
    q_ghb = cbc.get_ts(idx=(0, 0, ncol - 1), text=' HEAD DEP BOUNDS')[:, 1]
    ghb_nm = ["qghb_{}".format(int(i)) for i in totims]
    sim_obs = obs_utils.add_obs(sim_obs, obsnams=ghb_nm, simval=q_ghb,
                                obsval=np.zeros_like(q_ghb),
                                obsgnme='ghbQ')

    sim_obs.to_csv(os.path.join(model_ws, "output_sim.csv"), index=False)
    np.savetxt("output_sim.dat", sim_obs['simval'].values)

def monte_carlo():
    pass


if __name__ == "__main__":
    figs = generate_1d_model()
    x = 1
    forward_run()
    x = 1
