import numpy as np

def verifica_courant(dx, Q, dc, kcerc):
    tmin = (np.min(dx) ** 2) * (np.min(dc) / (4 * np.max(np.abs(Q) * np.max(kcerc))))
    return tmin