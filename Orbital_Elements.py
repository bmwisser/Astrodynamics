def rv2oe(r_2,v_2):
    import numpy as np
    from numpy.linalg import norm
    import math

    if norm(r_2) <= 500:
        mu = 1
    else:
        mu = 3.986012e5

    r2 = norm(r_2)
    v2 = norm(v_2)

    h0 = np.cross(r_2, v_2)
    h_mag = norm(h0)
    energy = (v2**2)/2 - mu/r2
    a = - mu / (2 * energy)
    p = h_mag**2 / mu
    e = np.sqrt(1 - (p / a))

    vva = [(v2**2 - mu/r2)/mu*i for i in r_2]
    rdvv = [np.dot(r_2,v_2)/mu*i for i in v_2]

    e_vec = [i-j for i,j in zip(vva,rdvv)] #; print(e_vec, e_vec[2])

    n = np.cross([0,0,1],h0)

    i = math.degrees(np.arccos(np.dot([0,0,1],h0)/h_mag))

    w = math.degrees(np.arccos(np.dot(n, e_vec) / (e * norm(n))))

    RAAN = math.degrees(np.arccos(np.dot([1, 0, 0], n) / norm(n)))

    nu2 = math.degrees(np.arccos(np.dot(e_vec, r_2) / (e * r2)))

    if n[1] < 0:
        RAAN = 360 - math.degrees(np.arccos(np.dot([1,0,0],n)/norm(n)))

    if e_vec[2] < 0:
        w = 360 - math.degrees(np.arccos(np.dot(n,e_vec)/(e*norm(n))))

    if np.dot(r_2,v_2) < 0:
        nu2 = 360 - nu2


    print('\na = '+str(round(a,5)) + '\n')
    print('e = '+str(round(e,5))+'\n')
    print('i = '+str(round(i,5))+u'\xb0 \n')
    print('\u03A9 = ' + str(round(RAAN, 5)) + u'\xb0 \n')
    print('\u03C9 = '+str(round(w,5))+u'\xb0 \n')
    print('\u03BD = '+str(round(nu2,5))+u'\xb0 \n')
    print('State vector = ' + str([a,e,i,RAAN,w,nu2]))

