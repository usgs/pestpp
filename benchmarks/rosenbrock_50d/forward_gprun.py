import pandas as pd
import numpy as np
import glob
import os
import laGPy as gpr
from scipy.stats import norm
import pyemu

temp_dir = "."
#outer_dirs = sorted(glob.glob("outer_*"), key=lambda x: int(os.path.basename(x).split("_")[1]))
pst = pyemu.Pst(glob.glob(os.path.join(temp_dir, "*.pst"))[0])

def emulate(pvals = None):
    if pvals is None:
        decvar = pd.read_csv(os.path.join(temp_dir, "dv.dat")).values.transpose()
    else:
        pvals_ordered = {pval: pvals[pval] for pval in sorted(pvals.index, key=lambda x: int(x[1:]))}
        decvar = np.array(list(pvals_ordered.values())).transpose()
    gp_list = sorted(glob.glob(os.path.join(temp_dir, "*.gp")), 
                     key=lambda x: int(os.path.splitext(os.path.basename(x))[0].split('_')[-1]))
    #preds = [gpr.loadGP(fname=gp).predict_lite(decvar) for gp in gp_list]  
    
    X = pd.read_csv(os.path.join(temp_dir, "gp_0.dv_training.csv")).drop(columns=['real_name']).values
    y = pd.read_csv(os.path.join(temp_dir, "gp_0.obs_training.csv")).drop(columns=['real_name'])['func'].values
    pred = gpr.laGP(Xref=decvar, start=10, end=60, X=X, Z=y)
    sim = {
        'func': pred["mean"].item(),
        'func_sd': np.sqrt(pred["s2"].item()),
        'func_var': pred["s2"].item(),
        'func_var_sd': 0,
        'ei': 0
    }

    #compute EI
    curr_opt = pd.read_csv(os.path.join(temp_dir, "curr_opt.csv"))['func'].values
    s = np.sqrt(sim['func_var'])
    if np.isclose(s, 0):
        ei = 0
    else:
        z = (curr_opt - sim['func']) / s
        ei = (curr_opt - sim['func']) * norm.cdf(z) + s * norm.pdf(z)
    sim['ei'] = ei.item()

    # #cluster analysis
    # training_dv = pd.read_csv(os.path.join(temp_dir, "gp_0.dv_training.csv"))
    # cluster_info = pd.read_csv(os.path.join(temp_dir, "cluster_info.csv"))
    # cluster_summary = pd.read_csv(os.path.join(temp_dir, "cluster_summary.csv"))
    # n_clusters = cluster_summary.shape[0]

    # from sklearn.metrics import pairwise_distances
    # training_dv_values = training_dv.drop(columns=['real_name']).values
    # distances = pairwise_distances(decvar.reshape(1, -1), training_dv_values)
    # nearest_index = np.argmin(distances)

    # # Get the nearest member's real name and corresponding values
    # nearest_member = training_dv.iloc[nearest_index]
    # nearest_real_name = nearest_member['real_name']
    # nearest_cluster_size = cluster_info.loc[cluster_info['real_name'] == nearest_real_name, 'avg_cluster_size'].values[0]
  
    # cluster_diffct = (1.5 * training_dv.shape[0] / n_clusters) - nearest_cluster_size
    # sim['cluster_diffct'] = cluster_diffct

    with open('output.dat','w') as f:
        f.write('obsnme,obsval\n')
        f.write('func_sd,'+str(sim["func_sd"])+'\n')
        f.write('func,'+str(sim["func"])+'\n')
        f.write('func_var,'+str(sim["func_var"])+'\n')
        f.write('func_var_sd,'+str(sim["func_var_sd"])+'\n')
        f.write('ei,'+str(sim["ei"])+'\n')
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
        sim = emulate(pvals=pvals)
        obs.update(sim)
        ppw.send_observations(obs.values)
        pvals = ppw.get_parameters()
        if pvals is None:
            break

if __name__ == "__main__":
    emulate()
