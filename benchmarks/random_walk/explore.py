import numpy as np
import matplotlib.pyplot as plt

def gx(x):
    y = 7.0/12 * np.power(x, 3) - 7/2 * np.power(x, 2.0) + 8 * x
    return y

n_time_steps = 3000
nreal = 500
std = 1;
vals = [np.array([0]*nreal)]
for i in range(n_time_steps-1):
    if i == 0:
        newValue =  3 * np.random.randn() + vals[i] + std * np.random.randn(nreal)
    else:
        newValue = 3 * np.random.randn() + vals[i]   + std * np.random.randn(nreal)

    vals.append(newValue)
A = np.array(vals)
xx = 1