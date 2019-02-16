# Import libraries
import numpy as np
from numpy import cross,dot
from numpy.linalg import norm
from sympy import solveset, Symbol, Eq, S
from Astro_Functions.ECI import ECI
from Astro_Functions.Stumpff import Stumpff
import matplotlib.pyplot as plt
from Astro_Functions.Orbital_Elements import rv2oe

###################################
#            Inputs               #
###################################

# mu = 3.9860e5
# # H = 1
# # Lat = np.radians(40)
# # time = [0, 118.104, 237.577]
# #
# # RA = np.radians([43.5365, 54.4196, 64.3178])
# # DEC = np.radians([-8.78334, -12.0739, -15.1054])
# # LST = np.radians([44.5065, 45.000, 45.4992])

mu = 3.9860e5
H = 0
Lat = np.radians(29)
time = [0, 60.0, 120.0]

RA = np.radians([0.0, 65.9279, 79.850])
DEC = np.radians([51.5110, 27.9911, 14.6609])
LST = np.radians([0.0, 0.250684, 0.501369])

tau1 = time[0] - time[1]
tau3 = time[2] - time[1]
tau = time[2] - time[0]

###################################
#             R vectors           #
###################################

R_list = []
for i in range(0,len(LST)):
    R_list.append(ECI(H,Lat,LST[i]))

R1 = np.array(R_list[0]); R2 = np.array(R_list[1]); R3 = np.array(R_list[2])

###################################
#         Rho_hat vectors         #
###################################

rho_list = []
for i in range(0,len(LST)):
    Lx = np.cos(RA[i])*np.cos(DEC[i])
    Ly = np.sin(RA[i])*np.cos(DEC[i])
    Lz = np.sin(DEC[i])
    #print('\n' + str([Lx,Ly,Lz]) + '\n')

    rho_list.append([Lx,Ly,Lz])

rho1 = np.array(rho_list[0]); rho2 = np.array(rho_list[1]); rho3 = np.array(rho_list[2])

p1 = cross(rho2,rho3)
p2 = cross(rho1,rho3)
p3 = cross(rho1,rho2)

###################################
#             D Matrix            #
###################################

D0 = dot(rho1,p1)

D11 = dot(R1,p1)
D12 = dot(R1,p2)
D13 = dot(R1,p3)

D21 = dot(R2,p1)
D22 = dot(R2,p2)
D23 = dot(R2,p3)

D31 = dot(R3,p1)
D32 = dot(R3,p2)
D33 = dot(R3,p3)

###################################
#           r_2 and v_2           #
###################################

A = (-D12*(tau3/tau) + D22 + D32*(tau1/tau))/D0
B = (D12*(tau3**2 - tau**2)*(tau3/tau) + D32*(tau**2 - tau1**2)*(tau1/tau))/(6*D0)
E = dot(R2,rho2)
R2_mag2 = norm(R2)**2

a = -(A**2 + 2*A*E + R2_mag2)
b = -2*mu*B*(A+E)
c = -(mu*B)**2

x = Symbol('x')
func = Eq(x**8 + a*x**6 + b*x**3 + c, 0)
xx = solveset(func,x,domain=S.Reals);

R_2 = 0;
for i in xx:
    if i > 6378:
        R_2 = float(i)
        #print(R_2)

rho1_mag = (((6*(D31*tau1/tau3 + D21*tau/tau3)*R_2**3 + mu*D31*(tau**2 - tau1**2)*tau1/tau3))/(6*R_2**3 + mu*(tau**2 - tau3**2)) - D11)/D0 # print(rho1_mag)
rho2_mag = A + mu*B/R_2**3 #; print(rho2_mag)
rho3_mag = (((6*(D13*tau3/tau1 - D23*tau/tau1)*R_2**3 + mu*D13*(tau**2 - tau3**2)*tau3/tau1))/(6*R_2**3 + mu*(tau**2 - tau1**2)) - D33)/D0 #; print(rho3_mag)

c1 = tau3/tau * (1 + (mu/(6*R_2**3)*(tau**2 - tau3**2))) #C1
c3 = -tau1/tau * (1 + (mu/(6*R_2**3)*(tau**2 - tau1**2))) #C3

f1 = 1 - mu*tau1**2/(2*R_2**3)
f3 = 1 - mu*tau3**2/(2*R_2**3)

g1 = tau1 - mu*tau1**3/(6*R_2**3)
g3 = tau3 - mu*tau3**3/(6*R_2**3)

r_1 = R1 + rho1_mag*rho1 ; print('r_1 = ' + str(r_1))
r_2 = R2 + rho2_mag*rho2 ; print('r_2 = ' + str(r_2))
r_3 = R3 + rho3_mag*rho3 ; print('r_3 = ' + str(r_3))

v_2 = (-f3*r_1 + f1*r_3)/(f1*g3 - f3*g1) ; print('v_2 = ' + str(v_2))

###################################
#           Improvement           #
###################################
r_list = [norm(r_2)]
v_list = [norm(v_2)]
c1_list = [c1]
c3_list = [c3]
rho1_list = [rho1_mag]
rho2_list = [rho2_mag]
rho3_list = [rho3_mag]
d_rho = 1
k = 0
while d_rho > 1e-5:
    k+=1; print('\nN = ' + str(k))

    r2 = norm(r_2)
    v2 = norm(v_2)

    alpha = 2/r2 - v2**2/mu
    X = Symbol('X')
    SC = Stumpff(alpha*X**2,5); SS = SC[0]; CC = SC[1]

    ###################################
    #            X1 and X3            #
    ###################################

    func2 = Eq(np.sqrt(mu)*tau1,(dot(r_2,v_2)*X**2*CC)/np.sqrt(mu) + (1 - alpha*r2)*X**3*SS + r2*X)
    XX = solveset(func2,X,domain=S.Reals)
    for i in XX:
        X1 = i
        #print(X1)

    func3 = Eq(np.sqrt(mu)*tau3,(dot(r_2,v_2)*X**2*CC)/np.sqrt(mu) + (1 - alpha*r2)*X**3*SS + r2*X)
    XX = solveset(func3,X,domain=S.Reals);
    for i in XX:
        X3 = i
        #print(X3)

    SC1 = Stumpff(alpha*X1**2,5); S1 = SC1[0]; C1 = SC1[1]
    SC3 = Stumpff(alpha*X3**2,5); S3 = SC3[0]; C3 = SC3[1]

    f1n = 1 - (X1**2*C1)/r2 #; print(f1n)
    f3n = 1 - (X3**2*C3)/r2 #; print(f3n)
    g1n = tau1 - (X1**3*S1)/np.sqrt(mu) #; print(g1n)
    g3n = tau3 - (X3**3*S3)/np.sqrt(mu)# ; print(g3n)

    f1_avg = (f1+f1n)/2
    f3_avg = (f3+f3n)/2
    g1_avg = (g1+g1n)/2
    g3_avg = (g3+g3n)/2

    c1n = g3_avg/(f1_avg*g3_avg - f3_avg*g1_avg); c1_list.append(c1n)
    c3n = -g1_avg/(f1_avg*g3_avg - f3_avg*g1_avg); c3_list.append(c3n)

    rho1n = (-D11 + D21/c1n - c3n*D31/c1n)/D0; rho1_list.append(rho1n)
    rho2n = (-c1n*D12 + D22 - c3n*D32)/D0; rho2_list.append(rho2n) #; print(rho2n)
    rho3n = (-c1n*D13/c3n + D23/c3n - D33)/D0; rho3_list.append(rho3n)

    r_1n = R1 + rho1n*rho1 #; print(r_1n)
    r_2n = R2 + rho2n*rho2 ; print('r2 = ' + str(r_2n))
    r_3n = R3 + rho3n*rho3 #; print(r_3n)

    v_2n = (-f3_avg*r_1n + f1_avg*r_3n)/(f1_avg*g3_avg - f3_avg*g1_avg); print('v2 = ' + str(v_2n))

    d_rho = abs(rho2n - rho2_list[k-1])
    d_c1 = abs(c1n - c1_list[k-1])
    d_c3 = abs(c3n - c3_list[k-1])

    r_2 = np.array(r_2n).astype(np.float64); r_list.append(norm(r_2))
    v_2 = np.array(v_2n).astype(np.float64); v_list.append(norm(v_2))

    #print(X1)
    #print(X3)
    print('r2 = ' + str(norm(r_2)))
    print('v2 = ' + str(norm(v_2)))
    print('\u03c12 = ' + str(rho2n))
    print('\u0394\u03c1 = ' + str(d_rho))
    print('\u0394c1 = ' + str(d_c1))
    print('\u0394c3 = ' + str(d_c3))


print('----------------------------------------------------------------------')
print('\nr_2 = ' + str(r_2))
print('v_2 = ' + str(v_2))

print('\nr2 = ' + str(r2))
print('v2 = ' + str(v2))
print('\n----------------------------------------------------------------------\n\n')


###################################
#         Orbital Elements        #
###################################

rv2oe(r_2,v_2)

###################################
#         Plot Convergence        #
###################################

iter = []
for i in range(0,len(rho2_list)):
    iter.append(i)
#print(iter)

for i in iter:
    plt.plot(iter,rho2_list)

plt.xlim([0,k])
plt.ylim([570,572])
plt.xlabel('Iteration')
plt.ylabel('Magnitude of \u03c12')
#plt.legend()
plt.grid(True)
plt.title('Convergence of \u03c12')
plt.show()


