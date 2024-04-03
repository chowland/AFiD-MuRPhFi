import h5py
from scipy.special import lambertw
import afidtools as afid
import numpy as np

fld = '../examples/RainyBenard'
inputs = afid.InputParams(fld)
nxm = inputs.nxm
print(nxm)

mpar = afid.MoistParams(fld)
alpha = mpar.alpha
beta = mpar.beta
gamma = mpar.gamma

saturated_top = (mpar.qfixN==1)
if gamma < 0:
    if saturated_top:
        gamma = beta/(1 - np.exp(-alpha))
    else:
        gamma = beta
print(alpha, beta, gamma)

def W(z):
    return np.real(lambertw(z))

# Define vertical coordinate grid
zc = np.linspace(0, 1, nxm+1)
z = 0.5*(zc[1:] + zc[:-1])

# Set moist static energy
P = gamma
if saturated_top:
    Q = beta - 1 - gamma*(1 - np.exp(-alpha))
else:
    Q = beta - 1 - gamma
m = P + Q*z

# Define temperature profile
C = P + (Q - beta)*z
T = C - W(alpha*gamma*np.exp(alpha*C))/alpha

# Infer buoyancy profile
b = T + beta*z

# Compute humidity profile
q = (m - b)/gamma

if not saturated_top:
    # Calculate height above which humidity not at saturation
    A = beta - (1 + gamma)/(1 + gamma*alpha*np.exp(1-alpha))
    B = alpha*np.exp(1-alpha)*(1 + gamma)/(1 + gamma*alpha*np.exp(1-alpha))
    zc = 1 - (1 + gamma*alpha*np.exp(1-alpha))/alpha/(1+gamma)
    # Prescribe linear profiles above this level
    b[z > zc] = beta - 1 + A*(z[z > zc] - 1)
    q[z > zc] = B*(1 - z[z > zc])

import matplotlib.pyplot as plt
fig, ax = plt.subplots(layout='constrained')
ax.plot(b, z)
ax.plot(q, z)
ax.grid(True)
fig.savefig('drizzle.png')

with h5py.File(fld+'/drizzle.h5','w') as f:
    f['b'] = b
    f['q'] = q