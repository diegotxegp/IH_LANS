# Nombre del archivo: bahias.py
# Autor: Diego García Prieto (diegotxegp @ Github)
# Fecha de creación: septiembre de 2023
# 

import numpy as np

def calcula_jacobiano_transversal(Yeq, Yst, posero, kacr, kero, PLCS_CS, dt, dy0):

    Yeq = np.array(Yeq)

    num_elements = len(PLCS_CS)

    unitario_5x5 = np.eye(5)
    Jacobito = np.array([unitario_5x5] * num_elements)

    for i in range(num_elements):
        j = PLCS_CS[i]
        if posero[i] == 1:
            Jacobito[i,0,0] = 1 - kero[j] * dt
            Jacobito[i,0,1] = dt * (dy0[i] + Yeq[i] - Yst[j])
            Jacobito[i,0,4] = kero[j] * dt
            Jacobito[i,2,:] = 0
            Jacobito[i,3,:] = 0
        else:
            Jacobito[i,2,2] = 1 - kacr[j] * dt
            Jacobito[i,2,3] = dt * (dy0[i] + Yeq[i] - Yst[j])
            Jacobito[i,2,4] = kacr[j] * dt
            Jacobito[i,0,:] = 0
            Jacobito[i,1,:] = 0

    return Jacobito