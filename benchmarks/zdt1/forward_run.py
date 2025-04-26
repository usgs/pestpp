import os
import pandas as pd
import numpy as np

def zdt1(x):
    m = len(x)
    f1 = x[0]
    g = 1 + (9 / (m - 1)) * np.sum(x[1:m])
    h = 1 - np.sqrt(f1 / g)
    f2 = g * h
    return f1, f2 


def helper(pvals=None):
    if pvals is None:
        x = pd.read_csv("dv.dat").values.reshape(-1).tolist()
    else:
        pvals_ordered = {pval: pvals[pval] for pval in sorted(pvals.index, key=lambda x: int(x[1:]))}
        x = np.array(list(pvals_ordered.values()))
    sim = {"obj_1": zdt1(x)[0], "obj_2": zdt1(x)[1]}
    with open('output.dat','w') as f:
        f.write('obsnme,obsval\n')
        f.write('obj_1,'+str(sim["obj_1"])+'\n')
        f.write('obj_2,'+str(sim["obj_2"])+'\n')
        
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
