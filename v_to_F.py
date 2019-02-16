def v_to_F(v,e):
    import numpy as np
    
    F = np.arccosh((e+np.cos(v))/(1 + e*np.cos(v)))

    return F


