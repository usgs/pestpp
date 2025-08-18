import os
import sys
import shutil
import platform
import numpy as np
import pandas as pd
import platform
import pyemu
sys.path.append(os.path.join("..","benchmarks"))
import opt_test_suite_helper as mou_suite_helper

bin_path = os.path.join("test_bin")
if "linux" in platform.platform().lower():
    bin_path = os.path.join(bin_path,"linux")
elif "darwin" in platform.platform().lower():
    bin_path = os.path.join(bin_path,"mac")
else:
    bin_path = os.path.join(bin_path,"win")

bin_path = os.path.abspath("test_bin")
os.environ["PATH"] += os.pathsep + bin_path


bin_path = os.path.join("..","..","..","bin")
exe = ""
if "windows" in platform.platform().lower():
    exe = ".exe"
exe_path = os.path.join(bin_path, "pestpp-mou" + exe)

port = 4055
num_workers = 10



def ppd_part1_test():
    test_d = "ppd_fitness_test"
    parsum1 = ppd_fitness(ppd_beta = 0.5, master_dir = "master1",test_d=test_d)
    parsum2 = ppd_fitness(ppd_beta = 0.7, master_dir = "master2",test_d=test_d)

    pcloud1 = parsum1.loc[(parsum1['nsga2_front']==1) & (parsum1['generation']==0)]
    pcloud2 = parsum2.loc[(parsum2['nsga2_front']==1) & (parsum2['generation']==0)]

    #starting from the same deterministic population...cloud should be the same at generation 0
    assert pcloud1.shape[0] == pcloud2.shape[0]

    #minimium artificial variance value check
    sd_syn = (pcloud1['obj_1'].max() - pcloud1['obj_1'].min())/pcloud1.shape[0]
    print(sd_syn)
    print(pcloud1['obj_1_sd_syn'].max())
    print(pcloud1['obj_1_sd_syn'].min())
    assert abs(pcloud1['obj_1_sd_syn'].max() - sd_syn) < 1e-6
    assert abs(pcloud1['obj_1_sd_syn'].min() - sd_syn) < 1e-6

    pcloud1 = parsum1.loc[(parsum1['nsga2_front']==1) & (parsum1['generation']==1)]
    pcloud2 = parsum2.loc[(parsum2['nsga2_front']==1) & (parsum2['generation']==1)]
    #pareto cloud with tighter beta should be smaller/narrower
    assert pcloud1.shape[0] < pcloud2.shape[0]

    max_nn = pcloud2['nn_count'].max()
    min_nn = pcloud2['nn_count'].min()
    pcloud2 = pcloud2.loc[pcloud2['nsga2_crowding_distance']!=1]
    pcloud2.sort_values(by='obj_1', inplace=True)
    fitness = 1 - (pcloud2['nn_count'] - min_nn)/(max_nn - min_nn + 1)
    #crowding distance fitness value check
    assert np.all(abs(fitness.values - pcloud2['nsga2_crowding_distance'].values) <= 1e-6)

    pcloud1 = parsum1.loc[parsum1['generation']==1]
    pcloud2 = parsum2.loc[parsum2['generation']==1]
    #there should be more fronts when beta is tighter
    assert pcloud1['nsga2_front'].max() > pcloud2['nsga2_front'].max()


def ppd_fitness(ppd_beta = 0.5, master_dir = "master",test_d="ppd_fitness_test_sklearn"):

    t_d = os.path.join(test_d, "template")
    m_d = None
    if master_dir is not None:
        m_d = os.path.join(test_d, master_dir)

    if m_d is not None and os.path.exists(m_d):
        shutil.rmtree(m_d)
    pst = pyemu.Pst(os.path.join(t_d, "pest.pst"))
    obs = pst.observation_data
    obs["link_to"] = np.nan
    obs.loc["obj_1","link_to"] = "obj_1_stdev"
    obs.loc["obj_2","link_to"] = "obj_2_stdev"
    
    pst.control_data.noptmax = 1

    pst.pestpp_options['mou_generator'] = 'pso'
    pst.pestpp_options['mou_env_selector'] = 'NSGA_PPD'
    pst.pestpp_options['mou_ppd_beta'] = ppd_beta
    pst.pestpp_options['mou_dv_population_file'] = 'fitness_test.dv_pop.csv'
    pst.pestpp_options['mou_obs_population_restart_file'] = 'fitness_test.obs_pop.csv'

    pst.write(os.path.join(t_d, "pest_ppd_fitness.pst"),version=2)
    
    sys.path.insert(0, t_d)
    from forward_gprun import ppw_worker as ppw_function 
    pyemu.os_utils.start_workers(t_d, exe_path, "pest_ppd_fitness.pst", master_dir=m_d,
                                worker_root='.', port=port, num_workers=num_workers, 
                                ppw_function = ppw_function)
    sys.path.remove(t_d)

    return pd.read_csv(os.path.join(m_d, "pest_ppd_fitness.pareto.summary.csv"))


def ppd_part2():
    test_d = "ppd_fitness_test_sklearn"
    parsum1 = ppd_fitness(ppd_beta = 0.5, master_dir = "master1",test_d=test_d)
    parsum2 = ppd_fitness(ppd_beta = 0.8, master_dir = "master2",test_d=test_d)

    pcloud1 = parsum1.loc[(parsum1['nsga2_front']==1) & (parsum1['generation']==0)]
    pcloud2 = parsum2.loc[(parsum2['nsga2_front']==1) & (parsum2['generation']==0)]

    #starting from the same deterministic population...cloud should be the same at generation 0
    assert pcloud1.shape[0] == pcloud2.shape[0]

    #minimium artificial variance value check
    sd_syn = (pcloud1['obj_1'].max() - pcloud1['obj_1'].min())/pcloud1.shape[0]
    assert abs(pcloud1['obj_1_sd_syn'].max() - sd_syn) < 1e-6
    assert abs(pcloud1['obj_1_sd_syn'].min() - sd_syn) < 1e-6

    sd_syn = (pcloud1['obj_2'].max() - pcloud1['obj_2'].min())/pcloud1.shape[0]
    assert abs(pcloud1['obj_2_sd_syn'].max() - sd_syn) < 1e-6
    assert abs(pcloud1['obj_2_sd_syn'].min() - sd_syn) < 1e-6

    sd_syn = (pcloud2['obj_1'].max() - pcloud2['obj_1'].min())/pcloud2.shape[0]
    assert abs(pcloud2['obj_1_sd_syn'].max() - sd_syn) < 1e-6
    assert abs(pcloud2['obj_1_sd_syn'].min() - sd_syn) < 1e-6

    sd_syn = (pcloud2['obj_2'].max() - pcloud2['obj_2'].min())/pcloud2.shape[0]
    assert abs(pcloud2['obj_2_sd_syn'].max() - sd_syn) < 1e-6
    assert abs(pcloud2['obj_2_sd_syn'].min() - sd_syn) < 1e-6

    pcloud1 = parsum1.loc[(parsum1['nsga2_front']==1) & (parsum1['generation']==1)]
    pcloud2 = parsum2.loc[(parsum2['nsga2_front']==1) & (parsum2['generation']==1)]
    #pareto cloud with tighter beta should be smaller/narrower
    assert pcloud1.shape[0] < pcloud2.shape[0]

    max_nn = pcloud2['nn_count'].max()
    min_nn = pcloud2['nn_count'].min()
    pcloud2 = pcloud2.loc[pcloud2['nsga2_crowding_distance']!=1]
    pcloud2.sort_values(by='obj_1', inplace=True)
    fitness = 1 - (pcloud2['nn_count'] - min_nn)/(max_nn - min_nn + 1)
    #crowding distance fitness value check
    assert np.all(abs(fitness.values - pcloud2['nsga2_crowding_distance'].values) <= 1e-6)

    pcloud1 = parsum1.loc[parsum1['generation']==1]
    pcloud2 = parsum2.loc[parsum2['generation']==1]
    #there should be more fronts when beta is tighter
    assert pcloud1['nsga2_front'].max() > pcloud2['nsga2_front'].max()



def ppd_part3():
    test_d = "ppd_fitness_test_obslink"
    parsum1 = ppd_fitness(ppd_beta = 0.5, master_dir = "master1",test_d=test_d)
    parsum2 = ppd_fitness(ppd_beta = 0.8, master_dir = "master2",test_d=test_d)

    pcloud1 = parsum1.loc[(parsum1['nsga2_front']==1) & (parsum1['generation']==0)]
    pcloud2 = parsum2.loc[(parsum2['nsga2_front']==1) & (parsum2['generation']==0)]

    #starting from the same deterministic population...cloud should be the same at generation 0
    assert pcloud1.shape[0] == pcloud2.shape[0]

    #minimium artificial variance value check
    sd_syn = (pcloud1['obj_1'].max() - pcloud1['obj_1'].min())/pcloud1.shape[0]
    assert abs(pcloud1['obj_1_stdev_syn'].max() - sd_syn) < 1e-6
    assert abs(pcloud1['obj_1_stdev_syn'].min() - sd_syn) < 1e-6

    sd_syn = (pcloud1['obj_2'].max() - pcloud1['obj_2'].min())/pcloud1.shape[0]
    assert abs(pcloud1['obj_2_stdev_syn'].max() - sd_syn) < 1e-6
    assert abs(pcloud1['obj_2_stdev_syn'].min() - sd_syn) < 1e-6

    sd_syn = (pcloud2['obj_1'].max() - pcloud2['obj_1'].min())/pcloud2.shape[0]
    assert abs(pcloud2['obj_1_stdev_syn'].max() - sd_syn) < 1e-6
    assert abs(pcloud2['obj_1_stdev_syn'].min() - sd_syn) < 1e-6

    sd_syn = (pcloud2['obj_2'].max() - pcloud2['obj_2'].min())/pcloud2.shape[0]
    assert abs(pcloud2['obj_2_stdev_syn'].max() - sd_syn) < 1e-6
    assert abs(pcloud2['obj_2_stdev_syn'].min() - sd_syn) < 1e-6

    pcloud1 = parsum1.loc[(parsum1['nsga2_front']==1) & (parsum1['generation']==1)]
    pcloud2 = parsum2.loc[(parsum2['nsga2_front']==1) & (parsum2['generation']==1)]
    #pareto cloud with tighter beta should be smaller/narrower
    assert pcloud1.shape[0] < pcloud2.shape[0]

    max_nn = pcloud2['nn_count'].max()
    min_nn = pcloud2['nn_count'].min()
    pcloud2 = pcloud2.loc[pcloud2['nsga2_crowding_distance']!=1]
    pcloud2.sort_values(by='obj_1', inplace=True)
    fitness = 1 - (pcloud2['nn_count'] - min_nn)/(max_nn - min_nn + 1)
    #crowding distance fitness value check
    assert np.all(abs(fitness.values - pcloud2['nsga2_crowding_distance'].values) <= 1e-6)

    pcloud1 = parsum1.loc[parsum1['generation']==1]
    pcloud2 = parsum2.loc[parsum2['generation']==1]
    #there should be more fronts when beta is tighter
    assert pcloud1['nsga2_front'].max() > pcloud2['nsga2_front'].max()

def ppd_part_2_3_test():
    ppd_part2()
    #ppd_part3()

    # arc2 = pd.read_csv(os.path.join("ppd_fitness_test_sklearn","master1","pest_ppd_fitness.pareto.archive.summary.csv"))
    # arc3 = pd.read_csv(os.path.join("ppd_fitness_test_obslink","master1","pest_ppd_fitness.pareto.archive.summary.csv"))
    # assert arc2.shape == arc3.shape
    # for col in ["obj_1","obj_2"]:
    #     diff = np.abs(arc2[col] - arc3[col])
    #     print(diff.max())
    #     assert diff.max() < 1e-6
    # diff = np.abs(arc2["obj_2_sd"].values - arc3["obj_2_stdev"].values)
    # print(diff.max())
    # assert diff.max() < 1e-6
    

if __name__ == "__main__":
    #test_d = "ppd_fitness_test_sklearn"
    #parsum1 = ppd_fitness(ppd_beta = 0.5, master_dir = None,test_d=test_d)
    #ppd_part2()
    ppd_part1_test()
    ppd_part_2_3_test()


