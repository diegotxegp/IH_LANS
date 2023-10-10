import numpy as np

from Funciones_modelo.calcula_jacobiano_longitudinal import calcula_jacobiano_longitudinal
from Funciones_modelo.filtro_kalman_longitudinal import filtro_kalman_longitudinal

def kalman_longitudinal(Ylt, kcerc, vlt, dQdx, dt, DALCS, it, PLCS, Dc, Ber, sigmaK):

    kcerc = np.array(kcerc)
    vlt = np.array(vlt)
    Dc = np.array(Dc)
    Ber = np.array(Ber)

    estado_ant = np.array([Ylt[PLCS], kcerc[PLCS], vlt[PLCS]])
    Jacobito = calcula_jacobiano_longitudinal(kcerc[PLCS], Dc[PLCS] + Ber[PLCS], dQdx[PLCS], dt, sigmaK)
    estado_post, DALCS = filtro_kalman_longitudinal(estado_ant, Jacobito, DALCS, it)
    
    # Actualizamos el estado
    saltoYlt = estado_post[0, :] - estado_ant[0]
    Ylt[PLCS] = estado_post[0, :]
    kcerc[PLCS] = estado_post[1, :]
    vlt[PLCS] = estado_post[2, :]
    
    return Ylt, kcerc, vlt, saltoYlt, DALCS