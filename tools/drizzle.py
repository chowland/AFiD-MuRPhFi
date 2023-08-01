import h5py
from scipy.special import lambertw
import os
import afidtools as afid
import numpy as np

fld = os.environ['WORK']+'/RainyBenard/drizzle_IC'
inputs = afid.InputParams(fld)
nxm = inputs.nxm
print(nxm)

alpha = 3.0
beta = 2.0
gamma = beta/(1 - np.exp(-alpha))
print(alpha, beta, gamma)

def W(z):
    return np.real(lambertw(z))

# Define vertical coordinate grid
zc = np.linspace(0, 1, nxm+1)
z = 0.5*(zc[1:] + zc[:-1])

# Set moist static energy
P = gamma
Q = beta - 1 - gamma*(1 - np.exp(-alpha))
m = P + Q*z

# Define temperature profile
C = P + (Q - beta)*z
T = C - W(alpha*gamma*np.exp(alpha*C))/alpha

# Infer buoyancy profile
b = T + beta*z

# Compute humidity profile
q = (m - b)/gamma

import matplotlib.pyplot as plt
fig, ax = plt.subplots(layout='constrained')
ax.plot(b, z)
ax.plot(q, z)
fig.savefig('drizzle.png')

with h5py.File(fld+'/drizzle.h5','w') as f:
    f['b'] = b
    f['q'] = q