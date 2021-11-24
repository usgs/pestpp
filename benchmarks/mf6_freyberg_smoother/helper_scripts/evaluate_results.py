import os, sys
import pandas as pd
import matplotlib.pyplot as plt
import pyemu
pst = pyemu.Pst(r"..\template\freyberg6_ES.pst")
par = pd.read_csv(r"..\template\cncnc.csv")


ens_par = pd.read_csv(r"../template/prior_par.csv")

mean_par = par.mean(axis = 0)
par_truth = pd.read_csv(r"..\template\truth.par_data.csv")
par_truth = par_truth.set_index(['parnme'])
par_truth['Estimated'] = pd.DataFrame(mean_par)[0]


grps = ['npf_k33_0', 'npf_k33_1', 'npf_k33_2', 'npf_k_0', 'npf_k_1','npf_k_2', 'rch', 'sto_ss_0', 'sto_ss_1', 'sto_ss_2', 'sto_sy_0']
rug_kws={"color": "r", "alpha":0.3, "linewidth": 2, "height":0.4 }
import seaborn as sns
import numpy as np
from statsmodels.distributions.empirical_distribution import ECDF
import scipy.stats
def compute_cdf(x):
    xx = ECDF(x)
    return xx.x, xx.y

def plot_pdf(x):
    hist = np.histogram(x, bins=100)
    hist_dist = scipy.stats.rv_histogram(hist)
    #X = np.linspace(min(x), max(x), 100)
    plt.hist(hist_dist.pdf(x), density=True, bins=100)


for grp in grps:
    val = par_truth[par_truth['pargp'] == grp]
    prior = ens_par[val.index.values.tolist()].values.flatten()
    plt.figure()

    if 1:
        plt.figure()
        plt.title(grp)
        x, y = compute_cdf(np.log10(prior))
        plt.plot(x, y, label = 'prior')
        x, y = compute_cdf(np.log10(val['Estimated']))
        plt.plot(x, y, label='posterior')
        plt.legend()

    if 1:
        plt.figure()
        plt.title(grp)
        sns.distplot(np.log10(prior), hist=False, kde=True, bins=100)
        sns.distplot(np.log10(val['Estimated']), hist=False, rug=True, kde=True, bins=100)
        sns.distplot(np.log10(val['parval1']), hist=False, rug=True, rug_kws = rug_kws, kde=True, bins=1)
    xx = 1

xx = 1