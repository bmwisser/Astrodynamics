# Import libraries
import numpy as np
from Astro_Functions.Gibbs_Method import Gibbs

###################################
#units = input('DU or km?   \n')

# if units == 'DU':
#     mu = 1
#     a_e = 6378.145 / 6378.145
#     b_e = 6356.785 / 6378.145
# elif units == 'km':
mu = 3.986012e5
a_e = 6378.145  # km
b_e = 6356.785  # km

#Givens
# ~~~~~~~CHECK UNITS ~~~~~~~~
time = [0, 2, 4] ; time = [60*i for i in time] #min
LST = np.radians([60.0, 60.5014, 61.0027])
Az = np.radians([165.931, 145.967, 2.40962])
El = np.radians([9.53549, 45.7711, 21.8825])
rho = [1214.89, 421.441, 732.079] #km

Lat = np.radians(-20)
H = 0.5 #km
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Reference Ellipsoid eccentricity
e = 0.08182

x = ((a_e)/(np.sqrt(1 - e**2*(np.sin(Lat))**2)) + H)*np.cos(Lat)
z = ((a_e*(1-e**2))/(np.sqrt(1 - (e*np.sin(Lat))**2)) + H)*np.sin(Lat)

# Right ascension and declination

dec = []
h = []
alpha = []
for i in range(0,len(LST)):

    dec.append(np.arcsin(np.cos(Lat) * np.cos(Az[i]) * np.cos(El[i]) + np.sin(Lat) * np.sin(El[i])))
    
    if 0 < Az[i] < np.pi:
        h.append(2*np.pi - np.arccos( (np.cos(Lat)*np.sin(El[i]) - np.sin(Lat)*np.cos(Az[i])*np.cos(El[i]))/np.cos(dec[i]) ))
    else:
        h.append(np.arccos( (np.cos(Lat)*np.sin(El[i]) - np.sin(Lat)*np.cos(Az[i])*np.cos(El[i]))/np.cos(dec[i]) ))
    
    alpha.append(LST[i] - h[i])


#Observation point position vectors
r_list = []; R_list = []
for i in range(0,len(LST)):
    Lx = np.cos(alpha[i])*np.cos(dec[i])
    Ly = np.sin(alpha[i])*np.cos(dec[i])
    Lz = np.sin(dec[i])

    # print('\n' + str([Lx,Ly,Lz]) + '\n')

    r_list.append([k * rho[i] for k in [Lx,Ly,Lz]])

    R_list.append([x * np.cos(LST[i]), x * np.sin(LST[i]), z])

# print(r_list)
rho_1 = r_list[0]; rho_2 = r_list[1]; rho_3 = r_list[2]

# print(R_list)
R_1 = R_list[0]; R_2 = R_list[1]; R_3 = R_list[2]

r_1 = [i + j for i,j in zip(R_1,rho_1)]
r_2 = [i + j for i,j in zip(R_2,rho_2)]
r_3 = [i + j for i,j in zip(R_3,rho_3)]

# import pdb; pdb.set_trace()

Gibbs(r_1,r_2,r_3)



