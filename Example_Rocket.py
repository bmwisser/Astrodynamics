# Import libraries
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from datetime import datetime
startTime = datetime.now()

###################################
#            Constants            #
###################################
R_e = 6378e3
# mu = 3.986e5*1e9 #m^3/s^2
g0 = 9.81

Cd0 = 0.17
A = 20  # m^2

Ve1 = 3000
Ve2 = 4500

ep1 = 0.10
ep2 = 0.09

m_star = 10875

m02 = 40000
m01 = 240000

delta = .999999999945

####################################################
#           Stage One Initial Conditions           #
####################################################
B1 = -1200
gamma0 = (np.pi/2)*delta
H0 = 0
X0 = 0
V0 = 0.001

m0 = m01
B = B1

t_b1 = -(m0 - m02)*(1 - ep1)/B1

a1 = (-B1*Ve1/m01 - g0)
a2 = (-B1*Ve1/((m01-m02)*ep1 + m02) - g0)


def f(y,t):
    import numpy as np

    rho = 1.2*np.exp(-y[1]/8000)
    m = m0+B*t

    # vm = y[2]/342
    #
    # Cd = 10*np.sqrt(Cd0)*(vm**9 - (np.sqrt(Cd0)/12))*np.exp(-4.2*vm**3) + Cd0

    Cd = 0.155

    T1 = -Ve1*B
    D1 = .5*Cd*A*rho*y[2]**2

    f0 = y[2]*np.cos(y[3])  # dX/dt
    f1 = y[2]*np.sin(y[3])  # dH/dt
    f2 = (T1/m) - (D1/m) - (g0 - (y[2]*np.cos(y[3]))**2/(R_e + y[1]))*np.sin(y[3])  # dV/dt
    f3 = - (g0/y[2] - (y[2]*np.cos(y[3]))**2/((R_e + y[1])*y[2]))*np.cos(y[3])  # d(gamma)/dt

    return [f0, f1, f2, f3]

t1 = np.linspace(0,t_b1,1000)
init_cond = [X0, H0, V0, gamma0]

soln = odeint(f,init_cond,t1)

X = soln[:, 0]
H = soln[:, 1]
V = soln[:, 2]
gam = soln[:, 3]

# if H[-1]/1e3 > 68 and gam[-1] <.9 and V[-1]/1e3 > 3.4:
#     print([(H[-1]) / 1e3, (V[-1]) / 1e3, gam[-1]])

print([X[-1], H[-1], V[-1], gam[-1]])

# print(datetime.now() - startTime)

abc = plt.figure(1)
plt.subplot(2,2,1)
plt.plot(t1, H)
# plt.xlabel('Time (sec)')
plt.ylabel('Height (km)')

plt.subplot(2,2,2)
plt.plot(t1, X)
# plt.xlabel('Time (sec)')
plt.ylabel('Downrange Distance (km)')

plt.subplot(2,2,3)
plt.plot(t1, V)
plt.xlabel('Time (sec)')
plt.ylabel('Velocity (km/s)')

plt.subplot(2,2,4)
plt.plot(t1, gam) #*180/np.pi)
plt.xlabel('Time (sec)')
plt.ylabel('Flight Path Angle (rad)')

plt.show(abc)

###################################################
#          Stage Two Initial Conditions           #
###################################################
#
# del2 = 0.08 #np.linspace(0.005,0.015,100)
#
# BBB = np.linspace(-80,-400,20)
# for j in range(len(BBB)):
#     B2 = BBB[j]
#
#     m0 = m02
#     B = B2
#
#     gamma0 = gam[-1] - del2
#     H0 = H[-1]
#     X0 = X[-1]
#     V0 = V[-1]
#
#     t_b2 = -(m02 - m_star)*(1 - ep2)/B2
#
#     t2 = np.linspace(0,t_b2,2000)
#     init_cond = [X0, H0, V0, gamma0]
#
#     soln2 = odeint(f,init_cond,t2)
#
#     X2 = soln2[:, 0]
#     H2 = soln2[:, 1]
#     V2 = soln2[:, 2]
#     gam2 = soln2[:, 3]
#
#     # import pdb; pdb.set_trace()
#
#     #Actual orbit apogee
#     #################################
#
#     h_f = ((H[-1]+R_e)*V2[-1]*np.cos(gam2[-1]))
#     p_f = h_f**2/mu
#     Em_f = V2[-1]**2/2 - mu/(H2[-1]+R_e)
#     a_f = -mu/(2*Em_f)
#     e_f = np.sqrt(1 - p_f/a_f)
#     r_a = a_f*(1 + e_f)/1e3
#
#     dr_a = r_a - (35786 + 6378)

    #print([(H2[-1])/1e3,(V2[-1])/1e3,gam2[-1], r_a, dr_a])

#for i in range(len(t2)):
#t2[i] = t2[i] + t1[-1]

#
# abcd = plt.figure(1)
# plt.subplot(2,2,1)
# plt.plot(t2, H2/1e3)
# plt.plot(t1,H/1e3)
# # plt.xlabel('Time (sec)')
# plt.ylabel('Height (km)')
#
# plt.subplot(2,2,2)
# plt.plot(t2, X2/1e3)
# plt.plot(t1,X/1e3)
# # plt.xlabel('Time (sec)')
# plt.ylabel('Downrange Distance (km)')
#
# plt.subplot(2,2,3)
# plt.plot(t2, V2/1e3)
# plt.plot(t1,V/1e3)
# plt.xlabel('Time (sec)')
# plt.ylabel('Velocity (km/s)')
#
# plt.subplot(2,2,4)
# plt.plot(t2, gam2) #*180/np.pi)
# plt.plot(t1,gam)
# plt.xlabel('Time (sec)')
# plt.ylabel('Flight Path Angle (rad)')
#
# plt.show(abcd)








