import numpy as np

def calcula_cshore_md(YSTi, wbd, Hbd, SS0, AT0, kacr, kero, dy0, PCS, dt, Hberm):
    DYeq = -wbd[PCS] * (0.106 * Hbd[PCS] + SS0[PCS] + AT0[PCS]) / (Hberm[PCS] + 2 * Hbd[PCS])
    Yeq = dy0 + DYeq

    # Discretización hacia adelante en YST y hacia atrás en Yeq
    ystii = YSTi[PCS] + kacr[PCS] * dt * (Yeq - YSTi[PCS])

    # Corregimos erosiones
    posero = np.where(ystii < YSTi[PCS])
    posero_out = ystii < YSTi[PCS]
    ystii[posero] = YSTi[PCS][posero] + kero[PCS][posero] * dt * (Yeq[posero] - YSTi[PCS][posero])

    YSTii = np.zeros(YSTi.shape)
    YSTii[PCS] = ystii

    return YSTii, posero_out, DYeq