import os, sys
import pandas as pd
import shutil
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
#from da_engine import Analysis

template_ws = "template96_pst_enkf"


def f(state, t, error=[]):
    # The Lorenz 1996 model (Lorenz E., 1996. ``Predictability: a problem partly solved'')
    F = 8.0  # forcing term
    n = len(state)
    d = np.zeros(n)
    x = state
    for i in range(n):
        d[i] = (x[(i + 1) % n] - x[i - 2]) * x[i - 1] - x[i] + F
    return d

#def RK4():  # homegrown integrator not used here (see other scripts); odeint instead

def rand_evolve(state0, F=8.0, t=1, dt=0.01, is_random=0):
    
    #state0_truth = odeint(f, state0, t)

    _stat_0 = state0
    all_states = []
    #all_states.append(_stat_0)
    for t_ in t:  # move one time step
        states1 = odeint(f, _stat_0, np.array([0, dt]))
        _stat_0 = states1[1]
        if is_random == 0:
            model_noise = np.zeros(len(_stat_0))
        else:
            model_noise = np.random.randn(len(_stat_0)) * 0.1
        _stat_0 = _stat_0 + model_noise
        all_states.append(_stat_0)
        #if is_seq:  # just one time for filter
        #    break

    return np.array(all_states)


def determistic_forward_run(state0, F=8.0, t=1, error=[]):
    states1 = odeint(f, state0, t)
    return states1


def generate_truth(n, F):
    state0 = np.zeros(n) + F
    state0[19] += 0.01  # slightly perturb 20th variable
    t = np.arange(0.0, 20.0, 0.01)
    #k_true = forward_run(state0=state0, t=t)
    #return t, k_true
    state_truth = odeint(f, state0, t)
    #pd.DataFrame(state_truth)[14].plot()
    return state_truth


def forward_run():
    from scipy.integrate import odeint
    import numpy as np

    df_in = pd.read_csv(os.path.join("lorenz96_in.csv"))

    state0 = df_in.loc[df_in.name.str.startswith("x"), 'value']
    names = df_in.loc[df_in.name.str.startswith("x"), 'name'].apply(lambda x: x.split("_")[0])
    t_start = df_in.loc[df_in['name'] == 't_start', 'value'].values[0]
    t_end = df_in.loc[df_in['name'] == 't_end', 'value'].values[0]
    delt = df_in.loc[df_in['name'] == 'delt', 'value'].values[0]
    #nobs_t = df_in.loc[df_in['name'] == 'nobs_t', 'value'].values[0]
    #nobs_loc = df_in.loc[df_in['name'] == 'nobs_loc', 'value'].values[0]
    is_random = df_in.loc[df_in['name'] == 'is_random', 'value'].values[0]
    #is_seq = df_in.loc[df_in['name'] == 'is_seq', 'value'].values[0]
    t = np.arange(0.0, t_end - t_start, delt)  # arange exclusive gotcha!
    #dt_m = nobs_t / ((t_end - t_start) / delt)
    #ind_m = (np.linspace(int(dt_m / delt), int(t_end / delt) - 1, int(nobs_t))).astype(int)
    #ind_loc = np.linspace(0, len(state0) - 1, int(nobs_loc)).astype(int)

    # measurement noise
    sig_m = 0.1
    #cov_m = sig_m**2 * np.eye(m)

    # TODO: some people use physical start instead of 8, 8, .., 8.01, ... by running model first...


    ki = rand_evolve(state0=state0, t=t)#, is_random=is_random)
    print(ki)
    onames,ovals = [],[]
    for iname,name in enumerate(names):
        time_names = ["{0}_{1:06.3f}".format(name,otime) for otime in t]
        onames.extend(time_names)
        ovals.extend(list(ki[:,iname]))
    print(len(onames),len(ovals))
    df_out = pd.DataFrame({'obsname':onames,'simval':ovals},index=onames)
    df_out.to_csv(os.path.join('lorenz96_out.csv'), index=False)


    # names = []
    # obsvalue = []
    # for io, obs_t in enumerate(ind_m):
    #     vals = []
    #     for _io, obs_l in enumerate(ind_loc):
    #         vals.append(ki[obs_t][obs_l])
    #     names = names + ['x{0}_{1}'.format(x, obs_t) for x in ind_loc]
    #     obsvalue = obsvalue + vals
    #
    #
    # df_out = pd.DataFrame(columns=['obsname', 'simval'])
    # obs_name = names
    #
    # for indx, prefix in enumerate(['sim_x{}_'.format(x) for x in range(len(state0))]):
    #     dyn_name = [prefix + "{}".format(i) for i in range(len(t))]
    #     obs_name = obs_name + dyn_name
    #     current_value = ki[:, indx].flatten().tolist()
    #     obsvalue = obsvalue + current_value
    #
    # df_out['obsname'] = obs_name
    # df_out['simval'] = obsvalue
    #
    # df_out.to_csv(os.path.join(template_ws, 'lorenz96_out.csv'), index=False)



if __name__ == "__main__":
    forward_run()