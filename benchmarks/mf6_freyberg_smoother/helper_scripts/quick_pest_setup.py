
import os, sys
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyemu
sys.path.insert(0,r"D:\Workspace\Codes\curr_flopy\flopy")
import flopy

org_model_ws = os.path.join('freyberg_mf6')
os.listdir(org_model_ws)