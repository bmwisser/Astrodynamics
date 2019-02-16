import numpy as np
import math

y = float(input('Input year (1901 - 2099): '))
m = float(input('Input month (1-12): '))
d = float(input('Input day (1-31): '))

h = float(input('Input hour (24hr): '))
min = float(input('Input min: '))
sec = float(input('Input sec: '))

J0 = 367*y - math.floor((7/4)*(y + math.floor((m+9)/12))) + math.floor(275*m/9) + d + 1721013.5

#print(J0)

UTC = (h + min/60 + sec/3600)/24

print('Julian Date = ',J0 + UTC)

if 2025 <= y < 2075:
    T0 = (J0 - 2469807.5) / 36525
elif 1975 <= y < 2025:
    T0 = (J0-2451545.0)/36525
elif 1925 <= y < 1975:
    T0 = (J0-2433282.5)/36525
else:
    print('Shit, out of range')

print(T0)


#if  y >= 2000:
#    T0 = (J0-2451545.0)/36525
#else: # 1950<=y<2000
#    T0 = (J0-2433282.5)/36525

Theta_G0 = math.fmod((100.4606184 + 36000.77004*T0 + 0.000387933*T0**2 - 2.583e-8*T0**3),360)

Theta_G = math.fmod((Theta_G0 + 360.98564724*UTC),360)

if Theta_G0 <0:
    Theta_G0 = Theta_G0 + 360
if Theta_G <0:
    Theta_G = Theta_G + 360

print('\n'+'\u03B8_G0 = ' + str(round(Theta_G0,4)) + '\u00b0')
print('The Greenwich time is: ' + str(round(Theta_G,4)) + '\u00b0 \n')

Long = float(input('Input Longitude (East positive): '))
Theta = math.fmod((Theta_G + Long),360)

if Theta <0:
    Theta = Theta + 360

print('\nThe local Sidereal time is: ' + str(round(Theta,4)) + '\u00b0')








