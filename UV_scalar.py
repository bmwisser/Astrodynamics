# Import libraries
import numpy as np
from math import log
from sympy import Symbol, Eq, solve
from numpy.linalg import norm
from math import sqrt
from Astro_Functions.Stumpff import Stumpff

# Constants
mu = 3.986012e5
#mu = 1

# Givens
r_mag = 10000 #km
v_mag = 10 #km/s
nu0 = 30 * (np.pi/180) #radians
del_t = 1 * 3600 #seconds

# Calculations
energy = (v_mag ** 2) / 2 - mu / (r_mag)
a = - mu / (2 * energy) #print(a)

rp = r_mag*np.cos(nu0) #print(rp)
rq = r_mag*np.sin(nu0) #print(rq)
r0 =[rp, rq, 0]
#print('r0 = ' + str(r0))

e = Symbol('e')
f_e = Eq(r_mag*(1+e*np.cos(nu0)), a*(1-e**2))
ee = solve(f_e)

if type(ee) is list:
    e = [i for i in ee if i>0]
    e = float(e[0])
#print('e = ' + str(e))

p = a*(1-e**2)
coeff = sqrt(mu/p)

vp = - coeff * np.sin(nu0)
vq = coeff * (e + np.cos(nu0))
v0 = [vp, vq, 0]
#print('v0 = ' + str(v0) + '\n')

t = del_t

#print('t = ' + str(t) + ' sec')

# Approximation for x1
if a > 0:
    x1 = sqrt(mu) * (t) / a
else:
    if del_t>0:
        sign = 1
    else:
        sign = -1

    num = -2 * mu * (t)
    den = a * (np.dot(r0, v0) + sign * sqrt(-mu*a) * (1 - (r_mag/a)))

    x1 = sign * sqrt(-a)*log(num/den)

x = x1
#t0 = 0
dt = 1
i = 0
k = 0

data = []
xx = []
CC = []
SS = []

while (np.abs(dt) > 1e-12) and k < 1000:
    k += 1

    z = x ** 2 / a

    SC = Stumpff(z,20)   #; print(SC); print('\n')
    S = SC[0]
    C = SC[1]

    t1 = ((np.dot(r0, v0) / sqrt(mu)) * x ** 2 * C + (1 - (r_mag / a)) * x ** 3 * S + r_mag * x) / sqrt(mu)
    dt_dx = (x ** 2 * C + (np.dot(r0, v0) / sqrt(mu)) * x * (1 - z * S) + r_mag * (1 - z * C)) / sqrt(mu)

    x2 = x + (t - t1) / dt_dx

    dt = t - t1

    x = x2
    t0 = t1

    #print(dt)

    data.append(dt)
    xx.append(x)
    CC.append(C)
    SS.append(S)

#print('k = ' + str(k))

# pulling data from loop
x = xx[-1]
C = CC[-1]
S = SS[-1]
z = x ** 2 / a

# Calculations
f = 1 - (x**2/r_mag)*C
g = t - (x ** 3) * S / sqrt(mu)

r_a = [i * f for i in r0]
r_b = [i * g for i in v0]

r = [i + j for i, j in zip(r_a, r_b)]

r_new = norm(r)

f_dot = (sqrt(mu) / (r_new * r_mag)) * x * (z * S - 1)
g_dot = 1 - (x ** 2) * C / r_new

v_a = [i * f_dot for i in r0]
v_b = [i * g_dot for i in v0]

v = [i + j for i, j in zip(v_a, v_b)]
v_new = norm(v)

#import pdb; pdb.set_trace()

print('\n'+ 'New r = ' + str(r))

print('\n'+ 'New v = ' + str(v))

# print('\n' + str(f*g_dot - f_dot*g) + '  =  1?')




