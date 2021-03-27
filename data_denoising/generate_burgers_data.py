# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 13:59:28 2021
burgers_groundtruth.mat created using
burgers_equation.m

"""
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from PDE_FIND import *
import scipy.io as sio
import itertools

# load data
data=loadmat('data/burgers_groundtruth.mat')

eta=data['eta']
x=np.squeeze(data['x'])[1:-1]
t=np.squeeze(data['t'])[1:-1]
U=data['U'][1:-1,1:-1]
U_t=data['U_t']
U_x=data['U_x']
U_xx=data['U_xx']


# save as .npy
data = {}
data['eta']=eta
data['x'] = x
data['t'] = t
data['U'] = U
data['U_t'] = U_t
data['U_x'] = U_x
data['U_xx'] = U_xx
np.save('data/burgers_00',data)

# 
# noisy GLS solutions
# 

gamma = 1.0 # proportional error constant

# 1% error
noise_level = 0.01
U_noise = U + noise_level * np.abs(U)**gamma * np.random.normal(size=U.shape)
U_noise = (U_noise>0)*U_noise
data = {}
data['x'] = x
data['t'] = t
data['U'] = U_noise
data['gamma'] = gamma
np.save('data/burgers_01',data)

# 5% error
noise_level = 0.05
U_noise = U + noise_level * np.abs(U)**gamma * np.random.normal(size=U.shape)
U_noise = (U_noise>0)*U_noise
data = {}
data['x'] = x
data['t'] = t
data['U'] = U_noise
data['gamma'] = gamma
np.save('data/burgers_05',data)

# 10% error
noise_level = 0.10
U_noise = U + noise_level * np.abs(U)**gamma * np.random.normal(size=U.shape)
U_noise = (U_noise>0)*U_noise
data = {}
data['eta'] = eta
data['x'] = x
data['t'] = t
data['U'] = U_noise
data['gamma'] = gamma
np.save('data/burgers_10',data)

# 25% error
noise_level = 0.25
U_noise = U + noise_level * np.abs(U)**gamma * np.random.normal(size=U.shape)
U_noise = (U_noise>0)*U_noise
data = {}
data['x'] = x
data['t'] = t
data['U'] = U_noise
data['gamma'] = gamma
np.save('data/burgers_25',data)

# 50% error
noise_level = 0.50
U_noise = U + noise_level * np.abs(U)**gamma * np.random.normal(size=U.shape)
U_noise = (U_noise>0)*U_noise
data = {}
data['x'] = x
data['t'] = t
data['U'] = U_noise
data['gamma'] = gamma
np.save('data/burgers_50',data)



