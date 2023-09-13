import numpy as np
from scipy.optimize import fsolve

from Funciones_modelo.waveguo import waveguo

def snell_shoalrefrot(H, T, DIR, h0, gamma, nbati):
    # Función a hacer cero: H0*kr*ks-hbgamma=0
    def prot(H0, T, alpha0, c0, cg0, h1, gamma):
        angle_arg = np.arcsin(np.sin(np.deg2rad(alpha0)) / (c0 / T * (2 * np.pi / T) ** 2 * h1 / 9.81 *
            (1 - np.exp(-(2 * np.pi / T * np.sqrt(h1 / 9.81)) ** 2.5)) ** (-0.4) / h1) ** (-1) * 2 * np.pi)
        kr_term = np.sqrt(np.cos(np.deg2rad(alpha0)) / (np.cos(angle_arg)))
        ks_term = np.sqrt(cg0 / (0.5 * (1 + (4 * np.pi * h1 / (
            (2 * np.pi / T) ** 2 * h1 / 9.81 * (1 - np.exp(-(2 * np.pi / T * np.sqrt(h1 / 9.81)) ** 2.5)) ** (-0.4) / h1) ** (-1) * 2 * np.pi)) / (
            np.sinh(4 * np.pi * h1 / (
            (2 * np.pi / T) ** 2 * h1 / 9.81 * (1 - np.exp(-(2 * np.pi / T * np.sqrt(h1 / 9.81)) ** 2.5)) ** (-0.4) / h1) ** (-1) * 2 * np.pi)))) * \
            ((2 * np.pi / T) ** 2 * h1 / 9.81 * (1 - np.exp(-(2 * np.pi / T * np.sqrt(h1 / 9.81)) ** 2.5)) ** (-0.4) / h1) ** (-1) * 2 * np.pi / T
        return H0 * kr_term * ks_term - gamma * h1

    # Calculamos longitud de onda offshore
    L0, k0 = waveguo(T, h0)
    c0 = L0 / T
    n0 = 0.5 * (1 + 2 * k0 * h0 / (np.sinh(2 * k0 * h0)))
    cg0 = c0 * n0
    
    # Calculamos dirección wrt normal saliente de la batimetría
    Dir_dif = (DIR - nbati + 180) % 360 - 180

    hsol = np.zeros_like(h0)
    
    options = {'xtol': 0.001}
    for it in range(len(H)):
        if abs(Dir_dif[it]) < 90 and H[it] > 0.15:  # Oleaje entrante grande (el pequeño da errores numéricos)
            try:
                hsol[it] = fsolve(lambda h1: prot(H[it], T[it], Dir_dif[it], c0[it], cg0[it], h1, gamma), H[it] / gamma, **options)
            except:
                try:
                    hsol[it] = fsolve(lambda h1: prot(H[it], T[it], Dir_dif[it], c0[it], cg0[it], h1, gamma), H[it] / gamma * 0.8, **options)
                except:
                    try:
                        hsol[it] = fsolve(lambda h1: prot(H[it], T[it], Dir_dif[it], c0[it], cg0[it], h1, gamma), H[it] / gamma * 0.5, **options)
                    except:
                        hsol[it] = H[it] / gamma
            if np.isnan(hsol[it]):
                # Por ola pequeña
                hsol[it] = h0[it]
        else:
            hsol[it] = h0[it]  # Oleaje saliente, no propagamos

    # Calculamos a partir de la profundidad de rotura el resto de parámetros
    # Hb, alphab, Ks y Kr
    L1, k1 = waveguo(T, hsol)
    c1 = L1 / T
    n1 = 0.5 * (1 + 2 * k1 * hsol / (np.sinh(2 * k1 * hsol)))
    cg1 = c1 * n1
    Ks = np.sqrt(cg0 / cg1)
    Dir_difb = np.arcsin(np.sin(np.deg2rad(Dir_dif)) / (c0 / T * c1))
    alphab = (np.rad2deg(Dir_difb) + nbati) % 360
    Kr = np.sqrt(np.cos(np.deg2rad(Dir_dif)) / np.cos(Dir_difb))
    
    # Calculamos Hb
    Hb = H * Kr * Ks
    
    # Los no propagados, les afectamos las propiedades offshore
    pnoprop = np.where(hsol == h0)[0]
    pprop = np.setdiff1d(np.arange(len(hsol)), pnoprop)
    Ks[pnoprop] = 1
    Kr[pnoprop] = 1
    alphab[pnoprop] = DIR[pnoprop]
    Hb[pnoprop] = H[pnoprop]
    hb = hsol

    return Hb, Ks, Kr, alphab, hb, pprop, pnoprop