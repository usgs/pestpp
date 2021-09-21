import os, sys
import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv(r"cncnc.csv")
kh = [col for col in df.columns if 'kh' in col]
kk = df[kh].values
plt.imshow(kk.mean(axis = 0).reshape(50,50), cmap = 'jet')
plt.show()
xx = 1