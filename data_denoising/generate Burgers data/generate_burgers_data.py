# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 13:59:28 2021

@author: melis
"""
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from PDE_FIND import *
import scipy.io as sio
import itertools

# load data
data=sio.loadmat('./burgers.mat')
x= np.squeeze(np.real(data['x'][0]))[1:-1]
t =np.squeeze(np.real(data['t'][:,0]))[1:-1]
dt = t[1] - t[0]
dx = x[2]-x[1]

U = np.real(data['usol'][1:-1, 1:-1])

# first get derivatives U_t, U_x, U_xx

n, m = U.shape

Ut, R, rhs_des = build_linear_system(U, dt, dx, D=2, P=1, time_diff = 'FD', space_diff = 'FD')

U_t = Ut.reshape(len(x), len(t))
R = np.transpose(R)
U_x=R[2].reshape(len(x),len(t))
U_xx = R[4].reshape(len(x),len(t))

# need to transpose U, U_t, U_x, U_xx so in the right shape
# for the work being done
U = np.transpose(U)
U_t = np.transpose(U_t)
U_x = np.transpose(U_x)
U_xx = np.transpose(U_xx)

# save as .npy
data = {}
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



