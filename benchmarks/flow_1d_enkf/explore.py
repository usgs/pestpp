import os, sys
import flopy
import pandas as pd
import numpy as np
import pyemu


Lx = 1000.
Ly = 1.
ztop = 50.
zbot = 0.
nlay = 1
nrow = 1
ncol = 100
delr = Lx / ncol
delc = Ly / nrow
delv = (ztop - zbot) / nlay
botm = np.linspace(ztop, zbot, nlay + 1)
hk = 1.
vka = 1.
sy = 0.1
ss = 1.e-4
laytyp = 1

# Variables for the BAS package
# Note that changes from the previous tutorial!
ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
strt = 10. * np.ones((nlay, nrow, ncol), dtype=np.float32)

# Time step parameters
nper = 3
perlen = [1, 100, 100]
nstp =   [1, 100, 100]
steady = [True, False, False]

# Flopy objects
modelname = 'tutorial2'
mf = flopy.modflow.Modflow(modelname, exe_name='mf2005')
dis = flopy.modflow.ModflowDis(mf, nlay, nrow, ncol, delr=delr, delc=delc,
                               top=ztop, botm=botm[1:],
                               nper=nper, perlen=perlen, nstp=nstp, steady=steady)
bas = flopy.modflow.ModflowBas(mf, ibound=ibound, strt=strt)
upw = flopy.modflow.ModflowUpw(mf, hk=hk, vka=vka, sy=sy, ss=ss, laytyp=laytyp, ipakcb=53)
pcg = flopy.modflow.ModflowPcg(mf)

# boundary
condA = 50*1*50;
condB = 50*1*50;
stress_period_data = {0: [
    [0, 0, 0, 45, condA],
    [0, 0, 99, 1, condB],
    ]}
xx = 1