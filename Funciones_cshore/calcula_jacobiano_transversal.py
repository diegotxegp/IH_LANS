import numpy as np

def calcula_jacobiano_transversal(Yeq, Yst, posero, kacr, kero, PLCS_CS, dt, dy0):
    num_elements = len(PLCS_CS)
    Jacobito = np.tile(np.diag(np.ones(5)), (1, 1, num_elements))

    for i in range(num_elements):
        j = PLCS_CS[i]
        if posero[i] == 1:
            Jacobito[0, 0, i] = 1 - kero[j] * dt
            Jacobito[0, 1, i] = dt * (dy0[i] + Yeq[i] - Yst[j])
            Jacobito[0, 4, i] = kero[j] * dt
            Jacobito[2, :, i] = 0
            Jacobito[3, :, i] = 0
        else:
            Jacobito[2, 2, i] = 1 - kacr[j] * dt
            Jacobito[2, 3, i] = dt * (dy0[i] + Yeq[i] - Yst[j])
            Jacobito[2, 4, i] = kacr[j] * dt
            Jacobito[0, :, i] = 0
            Jacobito[1, :, i] = 0

    return Jacobito