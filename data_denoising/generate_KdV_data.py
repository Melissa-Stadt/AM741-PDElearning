import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.io import loadmat

# load solution data from .mat
data=loadmat('generate_KdV data/KdVgroundtruth.mat')

x=np.squeeze(data['x'])[1:-1]
t=np.squeeze(data['t'])[1:-1]
U=data['U'][1:-1,1:-1]
U_t=data['U_t']
U_x=data['U_x']
U_xxx=data['U_xxx']

# save as .npy
data = {}
data['x'] = x
data['t'] = t
data['U'] = U
data['U_t'] = U_t
data['U_x'] = U_x
data['U_xxx'] = U_xxx
np.save('data/KdV_00',data)

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
np.save('data/KdV_01',data)

# 5% error
noise_level = 0.05
U_noise = U + noise_level * np.abs(U)**gamma * np.random.normal(size=U.shape)
U_noise = (U_noise>0)*U_noise
data = {}
data['x'] = x
data['t'] = t
data['U'] = U_noise
data['gamma'] = gamma
np.save('data/KdV_05',data)

# 10% error
noise_level = 0.10
U_noise = U + noise_level * np.abs(U)**gamma * np.random.normal(size=U.shape)
U_noise = (U_noise>0)*U_noise
data = {}
data['x'] = x
data['t'] = t
data['U'] = U_noise
data['gamma'] = gamma
np.save('data/KdV_10',data)

# 25% error
noise_level = 0.25
U_noise = U + noise_level * np.abs(U)**gamma * np.random.normal(size=U.shape)
U_noise = (U_noise>0)*U_noise
data = {}
data['x'] = x
data['t'] = t
data['U'] = U_noise
data['gamma'] = gamma
np.save('data/KdV_25',data)

# 50% error
noise_level = 0.50
U_noise = U + noise_level * np.abs(U)**gamma * np.random.normal(size=U.shape)
U_noise = (U_noise>0)*U_noise
data = {}
data['x'] = x
data['t'] = t
data['U'] = U_noise
data['gamma'] = gamma
np.save('data/KdV_50',data)