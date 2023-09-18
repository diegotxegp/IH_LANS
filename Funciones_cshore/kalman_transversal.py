import numpy as np

from Funciones_cshore.calcula_jacobiano_transversal import calcula_jacobiano_transversal
from Funciones_cshore.filtro_kalman_transversal import filtro_kalman_transversal

def kalman_transversal(Yct, YCTi, Yeq, kacr, kero, dt, DACS, it, PLCS_CS, posero, dy0):
    estado_ant = np.vstack([Yct[PLCS_CS], kero[PLCS_CS] * posero, Yct[PLCS_CS], kacr[PLCS_CS] * (1 - posero), dy0])
    Jacobito = calcula_jacobiano_transversal(Yeq, YCTi, posero, kacr, kero, PLCS_CS, dt, dy0)
    estado_post, DACSi, saltoYct = filtro_kalman_transversal(estado_ant, Jacobito, DACS, it, posero, Yct)
    DACS = DACSi
    Yct[PLCS_CS] = estado_post[0, :] + estado_post[2, :]
    kero[PLCS_CS[posero]] = estado_post[1, posero]
    kacr[PLCS_CS[np.where(1 - posero)]] = estado_post[3, np.where(1 - posero)]
    dy0 = estado_post[4, :]
    return Yct, kacr, kero, saltoYct, DACS, dy0