import numpy as np

def calcula_jacobiano_longitudinal(kcerc, Dc, dQdx, dt, sigmaK):
    n = len(kcerc)
    Jacobito = np.zeros((3, 3, n))

    # Rellenamos la diagonal con 1's
    for i in range(n):
        Jacobito[:, :, i] = np.eye(3)

    # Calculamos las derivadas parciales
    for i in range(n):
        Jacobito[0, 2, i] = dt
        Jacobito[0, 1, i] = -1.0 / Dc[i] * sigmaK * kcerc[i] * dQdx[i] * dt

    return Jacobito