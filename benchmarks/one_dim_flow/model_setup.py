import os, sys
import shutil
import numpy as np
import flopy
import matplotlib.pyplot as plt

def generate_1d_model():
    modelname = "flow_1d"
    model_ws = r"model_dataset"

    if model_ws in os.listdir(os.getcwd()):
        shutil.rmtree(os.path.join(os.getcwd(), model_ws))
    os.mkdir(model_ws)

    exe = r"..\bin\win\mfnwt.exe"
    shutil.copy(src=exe, dst=os.path.join(model_ws, os.path.basename(exe)))

    mf = flopy.modflow.Modflow(modelname, model_ws=model_ws,
                               exe_name = os.path.abspath(exe),
                               version='mfnwt')

    Lx = 100.0
    Ly = 1.0
    ztop = 0.0
    zbot = -50.0
    nlay = 1
    nrow = 1
    ncol = 100
    delr = Lx / ncol
    delc = Ly / nrow
    delv = (ztop - zbot) / nlay
    botm = np.linspace(ztop, zbot, nlay + 1)

    dis = flopy.modflow.ModflowDis(
        mf, nlay, nrow, ncol, delr=delr, delc=delc, top=ztop, botm=botm[1:]
    )

    ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
    #ibound[:, :, 0] = -1
    #ibound[:, :, -1] = -1
    strt = np.ones((nlay, nrow, ncol), dtype=np.float32)
    strt[:, :, :] = -10.0
    bas = flopy.modflow.ModflowBas(mf, ibound=ibound, strt=strt)

    laytyp = np.ones(nlay)
    laywet = np.zeros(nlay)
    hk = np.zeros_like(bas.strt.array) + 10.0
    flopy.modflow.mfupw.ModflowUpw(mf, laytyp=laytyp, layavg=0, chani=1.0, layvka=0,
                                   laywet=laywet, hdry=-1e+30, iphdry=0, hk=hk, hani=1.0,
                                   vka=hk, vkcb=0.0, noparcheck=False, ipakcb = 55)

    ghb_data = [] #[lay, row, col,head,cond]
    ghb_data.append([0,0,0,-10,15])
    ghb_data.append([0,0,ncol-1, -15,15])
    ghb_stress_per = {}
    ghb_stress_per[0] = ghb_data
    ghbs = flopy.modflow.mfghb.ModflowGhb(mf, ipakcb=55, stress_period_data=ghb_stress_per,
                                          dtype=None,  no_print=False,
                                          options=None, extension='ghb')

    spd = {(0, 0): ["print head", "print budget", "save head", "save budget"]}
    oc = flopy.modflow.ModflowOc(mf, stress_period_data=spd, compact=True)

    nwt = flopy.modflow.mfnwt.ModflowNwt.load(r".\misc_files\solver_options.nwt",mf)

    mf.write_input()

    success, buff = mf.run_model()
    if not success:
        raise Exception("MODFLOW did not terminate normally.")

    hds = flopy.utils.HeadFile(os.path.join(model_ws,modelname+".hds"))
    wl = hds.get_data(totim = 1.0)

    plt.plot(wl[0][0])

if __name__ == "__main__":
    generate_1d_model()