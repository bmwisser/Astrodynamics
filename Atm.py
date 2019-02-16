def alt():

    # Import libraries
    import numpy as np
    import matplotlib.pyplot as plt

    ###################################
    #            Constants            #
    ###################################
    g0 = 9.81

    ###################################
    #           Calculations          #
    ###################################
    T0 = 288.16
    P0 = 101325
    rho0 = 1.225
    R = 287

    TT = [288.16]
    PP = [101325]
    den = [1.225]
    den2 = [1.225]

    for i in range(1,11001):
        TT.append(T0 + -6.5e-3*i)
        PP.append(P0*(TT[-1]/T0) ** ((-g0/(R*-6.5e-3))))
        den.append(rho0*(TT[-1]/T0) ** -((g0/(R*-6.5e-3)) +1))
        den2.append(PP[i] / (R * TT[i]))

    T0 = TT[-1]; P0 = PP[-1]; rho0 = den[-1]
    for i in range(1,(25000-11000)+1):
        TT.append(TT[-1])
        PP.append(P0 * (np.exp(-(g0 / (R * TT[-1]))))**(i))
        den.append(rho0 * (np.exp(-(g0 / (R * TT[-1])))) ** (i))
        den2.append(PP[i+11000] / (R * TT[i+11000]))

    T0 = TT[-1]; P0 = PP[-1]; rho0 = den[-1]
    for i in range(1,(47000-25000)+1):
        T = TT[-1] + 3e-3
        TT.append(T)
        PP.append(P0 * (TT[-1] / T0) ** ((-g0 / (R * 3e-3))))
        den.append(rho0 * (TT[-1] / T0) ** -((g0 / (R * 3e-3)) + 1))
        den2.append(PP[i+25000] / (R * TT[i+25000]))

    T0 = TT[-1]; P0 = PP[-1]; rho0 = den[-1]
    for i in range(1,(53000-47000)+1):
        TT.append(TT[-1])
        PP.append(P0 * np.exp((-g0 / (R * TT[-1])))**(i))
        den.append(rho0 * (np.exp(-(g0 / (R * TT[-1])))) ** (i))
        den2.append(PP[i+47000] / (R * TT[i+47000]))

    T0 = TT[-1]; P0 = PP[-1]; rho0 = den[-1]
    for i in range(1,(79000-53000)+1):
        T = TT[-1] + -4.5e-3
        TT.append(T)
        PP.append(P0 * (TT[-1] / T0) ** ((-g0 / (R * -4.5e-3))))
        den.append(rho0 * (TT[-1] / T0) ** -((g0 / (R * -4.5e-3)) + 1))
        den2.append(PP[i+53000] / (R * TT[i+53000]))

    T0 = TT[-1]; P0 = PP[-1]; rho0 = den[-1]
    for i in range(1,(90000-79000)+1):
        TT.append(TT[-1])
        PP.append(P0 * np.exp((-g0 / (R * TT[-1])))**(i))
        den.append(rho0 * (np.exp(-(g0 / (R * TT[-1])))) ** (i))
        den2.append(PP[i+79000] / (R * TT[i+79000]))

    T0 = TT[-1]; P0 = PP[-1]; rho0 = den[-1]
    for i in range(1,(105000-90000)+1):
        T = TT[-1] + 4e-3
        TT.append(T)
        PP.append(P0 * (TT[-1] / T0) ** ((-g0 / (R * 4e-3))))
        den.append(rho0 * (TT[-1] / T0) ** ((-g0 / (R * 4e-3)) + 1))
        den2.append(PP[i+90000] / (R * TT[i+90000]))

    return TT,PP,den

# import numpy as np
# import matplotlib.pyplot as plt


# h = list(np.linspace(0,105000,105001))
#
# plt.plot(TT, h)
#
# #plt.figure()
# plt.xlim([160,320])
# plt.ylim([0,110000])
# plt.xlabel('Temperature (K)')
# plt.ylabel('Altitude (m)')
# # plt.legend()
# plt.grid(True)
# plt.title('Standard Atmosphere')
# plt.show()
#
#
# plt.plot(PP, h)
#
# # plt.figure()
# plt.xlim([-5000,101325])
# plt.ylim([0,110000])
# plt.xlabel('Pressure (Pa)')
# plt.ylabel('Altitude (m)')
# # plt.legend()
# plt.grid(True)
# plt.title('Standard Atmosphere')
# plt.show()
#
# DEN = np.zeros(len(h))
# for i in range(len(h)):
#     DEN[i] = 1.225*np.exp(-h[i]/8000)

# plt.plot(den, h)
#
# # plt.figure()
# plt.xlim([-.05,1.225])
# plt.ylim([0,110000])
# plt.xlabel('Density (kg/m^3)')
# plt.ylabel('Altitude (m)')
# # plt.legend()
# plt.grid(True)
# plt.title('Standard Atmosphere')
# plt.show()


