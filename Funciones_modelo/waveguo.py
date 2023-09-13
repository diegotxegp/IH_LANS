import numpy as np

def waveguo(T, h):
    # Guo 2002
    w = 2 * np.pi / T

    kh = w ** 2 * h / 9.81 * (1 - np.exp(- (w * np.sqrt(h / 9.81)) ** 2.5)) ** (-0.4)

    k = kh / h

    L = 2 * np.pi / k

    return L, k