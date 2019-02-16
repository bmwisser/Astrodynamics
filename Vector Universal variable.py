# Import libraries
import numpy as np
from math import log
from numpy.linalg import norm
from math import sqrt
from Astro_Functions.Stumpff import Stumpff

#Constants
mu = 3.986012e5
#mu = 1

#Givens

# r0 = [1,0,0] #DU
# v0 = [0,0,1.1] #DU/TU
# del_t = 2 #TU

r0 = [7000.0, -12124, 0]
v0 = [2.6679, 4.6210, 0]
del_t = 60*60

#Calculations
r_mag = norm(r0)
v_mag = norm(v0)
h0 = np.cross(r0,v0)
h_mag = norm(h0)
energy = (v_mag**2)/2 - mu/(r_mag)
a = -mu/(2*energy)
p = h_mag**2/mu
e = sqrt(1 - (p/a))

nu0 = np.arccos( ((p/r_mag) - 1) / e ); #print(nu0)

t = del_t

# Approximation for x1
if a > 0:
    x1 = sqrt(mu) * (t) / a
else:
    if del_t>0:
        sign = 1

    else:
        sign = -1

    num = -2 * mu * (t-t0)
    den = a * (np.dot(r0, v0) + sign * sqrt(-mu*a) * (1 - (r_mag/a)))

    x1 = sign * sqrt(-a)*log(num/den)

x = x1
dt = 1
i = 0
k = 0

data = []
xx = []
CC = []
SS = []

while (np.abs(dt) > 1e-12) and k<1000:
    k += 1

    z = x ** 2 / a
    # if z > 0:
    #     C1 = (1 - np.cos(sqrt(z))) / z
    #     S1 = (sqrt(z) - np.sin(sqrt(z))) / sqrt(z ** 3)
    # else:
    #     C1 = (1 - np.cosh(sqrt(-z))) / z
    #     S1 = (np.sinh(sqrt(-z)) - sqrt(-z)) / sqrt((-z) ** 3)

    #print([S1, C1])

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

print('k = ' + str(k))

# pulling data from loop
x = xx[-1]
C = CC[-1]
S = SS[-1]
z = x ** 2 / a

#import pdb; pdb.set_trace()

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

print('\n'+ 'New r = ' + str(r))

print('\n'+ 'New v = ' + str(v))

print('\n' + str(f*g_dot - f_dot*g) + '  =  1?')

