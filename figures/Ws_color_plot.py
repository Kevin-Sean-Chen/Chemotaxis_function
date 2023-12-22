#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 19:50:14 2023

@author: kschen
"""
import numpy as np
from matplotlib import pyplot as plt
import scipy.io
import pandas as pd
import ssm
from scipy.optimize import minimize
from mpl_toolkits.axes_grid1 import make_axes_locatable
import json
import os

import matplotlib 
matplotlib.rc('xtick', labelsize=20) 
matplotlib.rc('ytick', labelsize=20) 

from scipy.io import loadmat

# %%
dir_ = '/projects/LEIFER/Kevin/Publications/Chen_learning_2023/data4plots/mutant_W_color.mat'
# Load the .mat file
mat_data = loadmat(dir_)

# Access the matrix (assuming it's stored with a specific variable name)
Ws = mat_data['Ws']

# %%
cmap = 'PiYG'
vmin,vmax = -0.3, 0.3

plt.figure()
plt.subplot(131)
plt.imshow(Ws[0,:,:], cmap=cmap, vmin=vmin, vmax=vmax)

plt.subplot(133)
plt.imshow(Ws[1,:,:], cmap=cmap, vmin=vmin, vmax=vmax)
plt.colorbar()

plt.subplot(132)
plt.imshow(Ws[2,:,:], cmap=cmap, vmin=vmin, vmax=vmax)

plt.savefig('Ws_color_plot.pdf', format='pdf', dpi=600)