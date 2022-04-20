import numpy as np
import csv
import matplotlib.pyplot as plt
from scipy.special import erfc

### INPUTS
# Prandtl & Schmidt numbers
Pr, Sc = 1.0, 10.0
# Density ratio
Rrho = 7e-3
###

# Diffusivity ratio
tau = Pr/Sc
# Far-field temperature (degree C)
DT = 5.0
# Stefan number
S = 2.5 # 80.0/DT
# Dimensionless liquidus slope
Lambda = 2.8e-3/Rrho
# Far-field salinity (g/kg)
DS = 3.87e-5*DT/7.86e-4/Rrho
print("Effective far-field salinity: ",DS," g/kg")

# Similarity variable values at interface (to be solved for)
alpha = np.arange(0.0,0.45,1e-6)
gamma = alpha/np.sqrt(tau)

lhs = alpha*np.exp(alpha**2)*erfc(-alpha)

qpi = np.sqrt(np.pi)

def f(gamma):
    return gamma*qpi*erfc(-gamma)/(
        np.exp(-gamma**2) + qpi*gamma*erfc(-gamma)
    )

rhs = (1 + Lambda)/S/qpi*(1 - Lambda/(1 + Lambda)*f(gamma))

# Find value of alpha which satisfies equation
a0 = np.interp(0, lhs - rhs, alpha)

# Calculate temperature and salinity coefficients
A = a0*np.exp(a0**2)*S*qpi
B = (1 + Lambda - A*erfc(-a0))/Lambda/erfc(-a0/np.sqrt(tau))

real_ratio = Rrho*(A*erfc(-a0)/B/erfc(-a0/np.sqrt(tau)))
print("Effective density ratio: ", real_ratio)
print("Interface melting temperature: ", DT*(1.0 - A*erfc(-a0)), "deg C")
print("Interface salinity: ", DS*(1.0 - B*erfc(-a0/np.sqrt(tau))), "g/kg")

x0 = 0.8
t0 = 1e-4
h0 = x0 + 2*a0*np.sqrt(t0)
x = np.linspace(0,1,1025)
T = 1.0 - A*erfc((x0 - x)/2.0/np.sqrt(t0))
S = 1.0 - B*erfc((x0 - x)/2.0/np.sqrt(tau*t0))
T[x>h0] = 1.0 - A*erfc(-a0)
# S[x>h0] = 0.0
plt.plot(x, T)
plt.plot(x, S)
plt.plot(h0, T[-1], 'o')
plt.plot(h0, 0.0, 'o')
plt.plot(h0, 1.0 - B*erfc(-a0/np.sqrt(tau)), 'o')
plt.show()

with open('pfparam.in', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow([A, B, a0])