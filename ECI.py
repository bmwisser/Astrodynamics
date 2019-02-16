def ECI(H,Lat,LST):

    import numpy as np

    a_e = 6378.15  # km
    f = 1/298.26  #0.003353

    x = (a_e/(np.sqrt(1 - (2*f - f**2)*(np.sin(Lat))**2)) + H)*np.cos(Lat)
    z = (((a_e*(1-f)**2)/np.sqrt(1 - (2*f - f**2)*(np.sin(Lat))**2)) + H)*np.sin(Lat)

    R = [x * np.cos(LST), x * np.sin(LST), z]

    return R

