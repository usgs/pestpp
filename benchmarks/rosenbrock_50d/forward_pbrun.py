import os
import pandas as pd
import numpy as np

def rosenbrock(x):
    x = np.array(x)
    # Rosenbrock function formula: sum_{i=1}^{d-1} [100(x_{i+1} - x_i^2)^2 + (1 - x_i)^2]
    return np.sum(100.0 * (x[1:] - x[:-1]**2)**2 + (1 - x[:-1])**2)


def helper(pvals=None):
    if pvals is None:
        x = pd.read_csv("dv.dat").values.reshape(-1).tolist()
    else:
        pvals_ordered = {pval: pvals[pval] for pval in sorted(pvals.index, key=lambda x: int(x[1:]))}
        x = np.array(list(pvals_ordered.values()))
    sim = {"func": rosenbrock(x), "func_sd": 0, "func_var": 0, "ei": 0, "cluster_diffct": 0}
    with open('output.dat','w') as f:
        f.write('obsnme,obsval\n')
        f.write('func,'+str(sim["func"])+'\n')
        f.write('func_sd,'+str(sim["func_sd"])+'\n')
        f.write('func_var,'+str(sim["func_var"])+'\n')
        f.write('ei,'+str(sim["ei"])+'\n')
        f.write('cluster_diffct,'+str(sim["cluster_diffct"])+'\n')
    return sim

def ppw_worker(pst_name,host,port):
    import pyemu
    ppw = pyemu.os_utils.PyPestWorker(pst_name,host,port,verbose=False)
    pvals = ppw.get_parameters()
    if pvals is None:
        return

    obs = ppw._pst.observation_data.copy()
    obs = obs.loc[ppw.obs_names,"obsval"]

    while True:

        sim = helper(pvals=pvals)

        obs.update(sim)
        
        ppw.send_observations(obs.values)
        pvals = ppw.get_parameters()
        if pvals is None:
            break


if __name__ == "__main__":
    helper()
