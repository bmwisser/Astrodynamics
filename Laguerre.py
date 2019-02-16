def Laguerre(r0,v0,dt):

    # Import libraries
    import numpy as np
    from numpy.linalg import norm
    from numpy import sqrt
    from Astro_Functions.Stumpff import Stumpff

    #Constants
    if norm(r0) > 500:
        mu = 3.986012e5
    else:
        mu = 1.0

    #Calculations
    r_mag = norm(r0)
    v_mag = norm(v0)
    h0 = np.cross(r0,v0)
    h_mag = norm(h0)
    energy = (v_mag**2)/2 - mu/r_mag
    a = -mu/(2*energy)
    alpha = -(2*energy)/mu
    p = h_mag**2/mu
    e = sqrt(1 - (p*alpha))

    nu0 = np.arccos(round( ((p/r_mag) - 1) / e, 12)) #; print(nu0)

    t = dt

    Period = 2*np.pi*sqrt(a**3 / mu)

    if dt > 0:
        sign = 1
    else:
        sign = -1

    if alpha < 0:
        dt = dt - sign * math.floor(abs(dt)/Period) * p

    #Loop Setup
    x0 = sqrt(mu)* dt*abs(alpha)
    x = x0
    error = 1
    k = 0
    tol = 1e-12

    data = []
    xx = []
    CC = []
    SS = []

    while (abs(error) > tol) and (k < 3000):
        k += 1

        z = x ** 2 * alpha

        SC = Stumpff(z, 20)
        S = float(SC[0])
        C = float(SC[1])

        F = (1 - (r_mag*alpha))*S*x**3 + np.dot(r0,v0)/sqrt(mu)*C*x**2 + r_mag*x - sqrt(mu)*t
        Fp = C*x**2 + np.dot(r0,v0)/sqrt(mu)*(1 - S*z)*x + r_mag*(1 - C*z)
        Fpp = (1 - (r_mag*alpha))*(1 - S*z)*x + np.dot(r0,v0)/sqrt(mu)*(1-C*z)

        delta = 2*sqrt(4*Fp**2 - 5*F*Fpp)

        if Fp > 0:
            sign = 1
        else:
            sign = -1

        dx = 5*F / (Fp + sign*delta)
        error = dx**2*alpha

        x = x - dx

        data.append(dx)
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

    # print('\n'+ 'New r = ' + str(r))

    # print('\n'+ 'New v = ' + str(v))

    # print('\n' + str(f*g_dot - f_dot*g) + '  =  1?')

    return r,v


