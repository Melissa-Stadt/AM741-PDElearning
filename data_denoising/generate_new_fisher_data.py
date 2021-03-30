'''
Generates the new_fisher data that is created from the .m files
in "generate Fisher groundtruth" folder
'''

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.io import loadmat

prefix = input('What is prefix? (e.g. new_fisher)?')

# load solution data from .mat
filename = prefix + '_groundtruth.mat'
data=loadmat('data/'+filename)

D=data['D']
r=data['r']
K=data['K']
x=np.squeeze(data['x'])[1:-1]
t=np.squeeze(data['t'])[1:-1]
U=data['U'][1:-1,1:-1]
U_t=data['U_t']
U_x=data['U_x']
U_xx=data['U_xx']

# save as .npy
data = {}
data['D'] = D
data['r'] = r
data['K'] = K
data['x'] = x
data['t'] = t
data['U'] = U
data['U_t'] = U_t
data['U_x'] = U_x
data['U_xx'] = U_xx
savefile = 'data/'+prefix+'_00'
np.save(savefile,data)


# 
# noisy GLS solutions
# 

gamma = 1.0 # proportional error constant

# 1% error
noise_level = 0.01
U_noise = U + noise_level * np.abs(U)**gamma * np.random.normal(size=U.shape)
U_noise = (U_noise>0)*U_noise
data = {}
data['D'] = D
data['r'] = r
data['K'] = K
data['x'] = x
data['t'] = t
data['U'] = U_noise
data['gamma'] = gamma
savefile = 'data/'+prefix+'_01'
np.save(savefile,data)

# 5% error
noise_level = 0.05
U_noise = U + noise_level * np.abs(U)**gamma * np.random.normal(size=U.shape)
U_noise = (U_noise>0)*U_noise
data = {}
data['D'] = D
data['r'] = r
data['K'] = K
data['x'] = x
data['t'] = t
data['U'] = U_noise
data['gamma'] = gamma
savefile = 'data/'+prefix+'_05'
np.save(savefile,data)

# 10% error
noise_level = 0.10
U_noise = U + noise_level * np.abs(U)**gamma * np.random.normal(size=U.shape)
U_noise = (U_noise>0)*U_noise
data = {}
data['D'] = D
data['r'] = r
data['K'] = K
data['x'] = x
data['t'] = t
data['U'] = U_noise
data['gamma'] = gamma
savefile = 'data/'+prefix+'_10'
np.save(savefile,data)

# 25% error
noise_level = 0.25
U_noise = U + noise_level * np.abs(U)**gamma * np.random.normal(size=U.shape)
U_noise = (U_noise>0)*U_noise
data = {}
data['D'] = D
data['r'] = r
data['K'] = K
data['x'] = x
data['t'] = t
data['U'] = U_noise
data['gamma'] = gamma
savefile ='data/'+prefix+'_25'
np.save(savefile,data)

# 50% error
noise_level = 0.50
U_noise = U + noise_level * np.abs(U)**gamma * np.random.normal(size=U.shape)
U_noise = (U_noise>0)*U_noise
data = {}
data['D'] = D
data['r'] = r
data['K'] = K
data['x'] = x
data['t'] = t
data['U'] = U_noise
data['gamma'] = gamma
savefile = 'data/'+prefix+'_50'
np.save(savefile,data)

