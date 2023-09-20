import numpy as np

from Funciones_modelo.calcula_jacobiano_longitudinal import calcula_jacobiano_longitudinal
from Funciones_modelo.ensambla_lc import ensambla_lc
from Funciones_modelo.filtro_kalman_longitudinal_transversal import filtro_kalman_longitudinal_transversal

from Funciones_cshore.calcula_jacobiano_transversal import calcula_jacobiano_transversal

def kalman_longitudinal_transversal(Ylt, kcerc, vlt, dQdx, dt, DALCS, it, PLCS, PLCS_CS_LCS, Dc, Ber, sigmaK, Yct, YCTi, Yeq, kacr, kero, posero, dy0):
    # Vectores de estado longitudinales
    estado_ant_l = np.array([Ylt[PLCS], kcerc[PLCS], vlt[PLCS]])
    
    # Vectores de estado transversales
    estado_ant_c = np.array([Yct[PLCS], kero[PLCS] * posero[PLCS_CS_LCS], Yct[PLCS], kacr[PLCS] * (1 - posero[PLCS_CS_LCS]), dy0[PLCS_CS_LCS]])
    
    # Cálculo de jacobianos
    Jacobito_l = calcula_jacobiano_longitudinal(kcerc[PLCS], Dc[PLCS] + Ber[PLCS], dQdx[PLCS], dt, sigmaK)
    Jacobito_c = calcula_jacobiano_transversal(Yeq[PLCS_CS_LCS], YCTi, posero, kacr, kero, PLCS, dt, dy0[PLCS_CS_LCS])
    Jacobito_lc = ensambla_lc(Jacobito_l, Jacobito_c)
    
    # Calculamos el estado posterior
    estado_ant = np.concatenate((estado_ant_l, estado_ant_c))
    estado_post, DALCS, saltoYlt, saltoYct = filtro_kalman_longitudinal_transversal(estado_ant, Jacobito_lc, DALCS, it, posero, Ylt, Yct)
    
    # Actualizamos el estado
    Ylt[PLCS] = estado_post[0, :]
    kcerc[PLCS] = estado_post[1, :]
    vlt[PLCS] = estado_post[2, :]
    Yct[PLCS] = estado_post[3, :] + estado_post[5, :]
    kero[PLCS[posero[PLCS_CS_LCS]]] = estado_post[4, posero[PLCS_CS_LCS]]
    kacr[PLCS[~posero[PLCS_CS_LCS]]] = estado_post[6, ~posero[PLCS_CS_LCS]]
    dy0[PLCS_CS_LCS] = estado_post[7, :]
    
    return Ylt, kcerc, vlt, saltoYlt, DALCS, Yct, kacr, kero, saltoYct, dy0