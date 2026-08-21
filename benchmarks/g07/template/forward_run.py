import os
import pandas as pd
import numpy as np

def evaluate(x):
    x = np.array(x)
    obj = x[0]**2 + x[1]**2 + x[0]*x[1] - 14*x[0] - 16*x[1] + (x[2] - 10)**2 + \
        4*(x[3] - 5)**2 + (x[4] - 3)**2 + 2*(x[5] - 1)**2 + 5*x[6]**2 + \
        7*(x[7] - 11)**2 + 2*(x[8] - 10)**2 + (x[9] - 7)**2 + 45
    
    # G07 constraints
    c1 = -105 + 4*x[0] + 5*x[1] - 3*x[6] + 9*x[7]
    c2 = 10*x[0] - 8*x[1] - 17*x[6] + 2*x[7]
    c3 = -8*x[0] + 2*x[1] + 5*x[8] - 2*x[9] -12
    c4 = 3*(x[0]-2)**2 + 4*(x[1]-3)**2 + 2*x[2]**2 -7*x[3] -120
    c5 = 5*x[0]**2 + 8*x[1] + (x[2]-6)**2 - 2*x[3] - 40
    c6 = x[0]**2 + 2*(x[1]-2)**2 - 2*x[0]*x[1] + 14*x[4] - 6*x[5]
    c7 = 0.5*(x[0] - 8)**2 + 2*(x[1]-4)**2 + 3*x[4]**2 - x[5] - 30
    c8 = -3*x[0] + 6*x[1] + 12*(x[8] - 8)**2 - 7*x[9]
    
    return np.array([obj, c1, c2, c3, c4, c5, c6, c7, c8])

def helper(pvals=None):
    if pvals is None:
        x = pd.read_csv("par.dat").values.reshape(-1).tolist()
    else:
        pvals_ordered = {pval: pvals[pval] for pval in sorted(pvals.index, key=lambda x: int(x[1:]))}
        x = np.array(list(pvals_ordered.values()))
    
    results = evaluate(x)
    obj = results[0]
    constraints = results[1:]
    
    sim = {"obj": obj}
    for i, c in enumerate(constraints):
        sim[f"constraint{i+1}"] = c
    
    with open('obs.dat','w') as f:
        f.write('obsnme,obsval\n')
        f.write('obj,'+str(sim["obj"])+'\n')
        for i in range(len(constraints)):
            f.write(f'g{i+1},{sim[f"constraint{i+1}"]}\n')

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
