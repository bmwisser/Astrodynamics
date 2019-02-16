def v_to_E(v,e):
    import numpy as np

    a = np.sqrt((1-e)/(1+e))
    
    E = 2*np.arctan(a*np.tan(v/2))
    return E




    
