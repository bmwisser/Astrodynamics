# Import libraries
import numpy as np
from numpy import linalg as LA
from Astro_Functions.v_to_E import v_to_E

#Constants
mu = 1;

#Givens
r0 = [1,0,0] #DU
v0 = [0,0,1.1] #DU/TU

#Calculations
r_mag = LA.norm(r0)
v_mag = LA.norm(v0)
t = 2

h0 = np.cross(r0,v0)
h0_mag = LA.norm(h0)
p = h0_mag**2/mu

energy = (v_mag**2)/2 - mu/(r_mag)
a = -1/(2*energy)
e = np.sqrt(1 - p/a)

ecosv = ((p/r_mag)-1)/e
nu_0 = np.arccos(round(ecosv,10))
E00 = v_to_E(nu_0,e)

E0 = np.pi
EE = []
M = np.sqrt(mu/a**3)*t
error = 1
#if e < 0:

while np.abs(error) > 1e-25:

    Mn = E0 - e*np.sin(E0)
    dM_dE = 1 - e*np.cos(E0)
    
    E1 = E0 + (M - Mn)/dM_dE
    
    error = M - Mn

    E0 = E1
    EE.append(E1)
    #print(E1)
    print(error)


E = EE[-1]
dE = E - E00

f = 1 - (a/r_mag)*(1 - np.cos(dE))
g = t - np.sqrt(a**3/mu) * (dE - np.sin(dE))

r_a = [i*f for i in r0]
r_b = [i*g for i in v0]

r = [i + j for i,j in zip(r_a,r_b)]
r_new = LA.norm(r)

f_dot = - (np.sqrt(mu*a)*np.sin(dE))/(r_new*r_mag)
g_dot = 1 - (a/r_new)*(1 - np.cos(dE))

v_a = [i*f_dot for i in r0]
v_b = [i*g_dot for i in v0]

v = [i + j for i,j in zip(v_a,v_b)]
v_new = LA.norm(v)

print(r)
print(v)













