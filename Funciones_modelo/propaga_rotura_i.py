import numpy as np

from Funciones_modelo.snell_shoalrefrot import snell_shoalrefrot
from Funciones_modelo.szonewidth import szonewidth

def propaga_rotura_i(H0, D0, Tp0, AT0, SS0, SLR0, refNMM, h0, dinperf, nbati_calc, cotasZ, gamma, PERF):
    H = H0[dinperf]
    T = Tp0[dinperf]
    DIR = D0[dinperf]
    nivel = h0[dinperf] + AT0[dinperf] + SS0[dinperf] + SLR0[dinperf] + refNMM
    Hi, _, _, Di, hb, _, _ = snell_shoalrefrot(H, T, DIR, nivel, gamma, nbati_calc)
    wdean = szonewidth(hb - SS0 - AT0, cotasZ, PERF)
    Hi = Hi.flatten()
    Di = Di.flatten()
    D0i = DIR.flatten()
    w0 = wdean.flatten()

    return Hi, D0i, Di, w0