import sys
import math
import os
import pandas as pd
import numpy as np

def fonseca_fleming(x):
    x = np.array(x)
    n = len(x)
    # Fonseca-Fleming objective functions
    # f1(x) = 1 - exp(-sum((xi - 1/sqrt(n))^2))
    # f2(x) = 1 - exp(-sum((xi + 1/sqrt(n))^2))
    term1 = np.sum((x - 1/np.sqrt(n))**2)
    term2 = np.sum((x + 1/np.sqrt(n))**2)
    
    obj1 = 1 - np.exp(-term1)
    obj2 = 1 - np.exp(-term2)
    
    return obj1, obj2

def helper(pvals=None):
    if pvals is None:
        x = pd.read_csv("dv.dat").values.reshape(-1).tolist()
    else:
        pvals_ordered = {pval: pvals[pval] for pval in sorted(pvals.index, key=lambda x: int(x[1:]))}
        x = np.array(list(pvals_ordered.values()))
    
    obj1, obj2 = fonseca_fleming(x)
    sim = {"obj1": obj1, "obj2": obj2,}
    
    with open('output.dat','w') as f:
        f.write('obsnme,obsval\n')
        f.write('obj1,'+str(sim["obj1"])+'\n')
        f.write('obj2,'+str(sim["obj2"])+'\n')
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
