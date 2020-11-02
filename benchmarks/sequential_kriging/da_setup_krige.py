import os, sys
import numpy as np
import EnsKF_Tools as KF
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import pyemu
import pandas as pd
## Generate random fields
n = 100
N = 200
m = 4
obs_offset = 10
obs_location = np.linspace(obs_offset, n-obs_offset, m).astype(int)

ref_real = 0 # will be used to generate the obs kriging
x = np.linspace(0,1,n)
x = np.stack((x, np.zeros_like(x)))
distMat = cdist(x.T, x.T, 'euclidean')
covMat = KF.cov(correlation_scale =0.3, stand_dev= 1.0,distMat = distMat)
K = KF.svd_random_generator(covM = covMat, seed = 6543, nreal = N)

obs = K[obs_location, ref_real]

H = K[obs_location, :] # model predictions
Ka = KF.EnsKF_evenson(H = H, K = K, d = obs, err_perc = 100, thr_perc = 5)

## Generate Pst files

# generate parnames and obs names
parnames = ["par_"+str(i) for i in range(n)]
obsnames = ["obs_"+str(i) for i in range(m)]

# get par and obs names to build pst object
par_df = pd.DataFrame(K.T, columns= parnames)
obs_df = pd.DataFrame(obs.reshape(1,m), columns= obsnames)

# generate pst object from parnames and obsnames
pst = pyemu.Pst.from_par_obs_names(par_names= parnames,
                                   obs_names= obsnames)

pyemu.utils.simple_tpl_from_pars(parnames, tplfilename = r"tplkrige.tpl")
pyemu.utils.simple_ins_from_obs2(obsnames, insfilename = r"inskrige.ins")

pst.parameter_data['partrans'] = 'none'
pst.control_data.noptmax = 1
# assign obs values and weights
pst.observation_data['obsval'] = obs
pst.observation_data['obsnme'] = obsnames
pst.observation_data['weight'] = 100
pst.model_command = sys.executable + " " + "forward_run.py"

xx = 1




