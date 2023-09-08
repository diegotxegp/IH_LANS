import numpy as np

def wmoore_calc(D50):
    wmoore = np.zeros_like(D50)
    
    for i in range(len(D50)):
        if D50[i] <= 0.1e-3:
            wmoore[i] = 1.1e6 * D50[i]**2
        elif 0.1e-3 < D50[i] <= 1e-3:
            wmoore[i] = 273 * D50[i]**1.1
        elif D50[i] >= 1e-3:
            wmoore[i] = 4.36 * D50[i]**0.5
    
    return wmoore