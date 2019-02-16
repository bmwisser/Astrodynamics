def Gibbs(r_1,r_2,r_3):
    import numpy as np
    import math
    from numpy.linalg import norm

    # r_1 = [1.41422511, 0, 1.414202]
    # r_2 = [1.81065659, 1.06066883, 0.3106515]
    # r_3 = [1.35353995, 1.41422511, -0.6464495]

    if norm(r_1) <= 500 and norm(r_2) <= 500:
        mu = 1
    else:
        mu = 3.986012e5

    r1 = norm(r_1)
    r2 = norm(r_2)
    r3 = norm(r_3)

    #print(np.dot(r_1, np.cross(r_2, r_3)))
    #print(mu)

    C23 = np.array(np.cross(r_2, r_3))/(r2*r3) #; print(C23)
    n = np.array(r_1)/r1 #; print(n)

    # print(np.dot(C23, n))

    if np.dot(C23,n) < 1e-5:

        N = r1*np.array(np.cross(r_2,r_3)) + r2*np.array(np.cross(r_3,r_1)) + r3*np.array(np.cross(r_1,r_2))
        D = np.array(np.cross(r_1,r_2)) + np.array(np.cross(r_2,r_3)) + np.array(np.cross(r_3,r_1))
        S = np.array(r_1)*(r2-r3) + np.array(r_2)*(r3-r1) + np.array(r_3)*(r1-r2)

        N_mag = norm(N)
        D_mag = norm(D)

        v_2 = np.sqrt(mu/(N_mag*D_mag)) * ((np.array(np.cross(D,r_2))/r2) + S); print(v_2)
        v2 = norm(v_2); print(v2)

    else:
        print('r1, r2, and r3 are not coplanar vectors')
        exit()

    h0 = np.cross(r_2, v_2)
    h_mag = norm(h0)
    energy = (v2** 2) / 2 - mu / (r2)
    a = - mu / (2 * energy)
    p = h_mag** 2 / mu
    e = np.sqrt(1 - (p / a))

    vva = [(v2**2 - mu/r2)/mu*i for i in r_2]
    rdvv = [np.dot(r_2,v_2)/mu*i for i in v_2]

    e_vec = [i-j for i,j in zip(vva,rdvv)] ; print(e_vec, e_vec[2])

    n = np.cross([0,0,1],h0)

    i = math.degrees(np.arccos(np.dot([0,0,1],h0)/h_mag))

    w = math.degrees(np.arccos(np.dot(n, e_vec) / (e * norm(n))))#; print(w)

    RAAN = math.degrees(np.arccos(np.dot([1, 0, 0], n) / norm(n)))#; print(RAAN)

    nu2 = math.degrees(np.arccos(np.dot(e_vec, r_2) / (e * r2)))  # ; print(nu2)

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

    # import pdb; pdb.set_trace()

# r111 = [-9025.0458584, -1773.6338262, 1083.3613982]
# r222 = [-8638.6007898, -4987.4984914, 0]
# r333 = [-6909.9359212, -7915.9089353, -1237.6481231]

# r11 = [5887, -3520, -1204]
# r22 = [5572, -3457, -2376]
# r33 = [5088, -3289, -3480]

# r1 = [-294.32, 4265.1, 5986.7]
# r2 = [-1365.5, 3637.6, 6346.8]
# r3 = [-2940.3, 2473.7, 6555.8]
#
# Gibbs(r1,r2,r3)








