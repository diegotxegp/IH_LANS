import numpy as np

from Generales.anglecalc import anglecalc

def snell_shoalref(Hs0, T0, Dir0, h0, hpoint, nbati):
    # Propagación Ley-Snell + Conservación de Flujo de Energía.
    
    # Obtención de la longitud de onda en el punto desde el que se propaga
    L0, k0 = waveguo(T0, h0)
    c0 = L0 / T0
    cg0 = (c0 / 2) * (1 + (2 * k0 * h0) / (np.sinh(2 * k0 * h0)))
    Dir_dif = circ_dist2((Dir0) * np.pi / 180, nbati * np.pi / 180) * 180 / np.pi
    
    # Obtención de la longitud de onda en el punto OBJETIVO
    L1, k1 = waveguo(T0, hpoint)
    c1 = L1 / T0
    n1 = 0.5 * (1 + (2 * k1 * hpoint) / (np.sinh(2 * k1 * hpoint)))
    cg1 = n1 * c1
    Ks = np.sqrt(cg0 / cg1)
    Dir_difb = np.arcsin(np.sin(Dir_dif) / c0 * c1)  # respecto a normal saliente
    
    # Deshacemos la resta y normalizamos de 0 a 360 grados
    alphab = Dir_difb + nbati
    alphab = anglecalc(np.sin(alphab), np.cos(alphab))
    
    Kr = np.sqrt(np.cos(Dir_dif) / np.cos(Dir_difb))
    
    # Calculamos Hb
    Hb = Hs0 * Kr * Ks
    
    # Los no propagados, les afectamos las propiedades offshore
    pnoprop = np.where(hpoint == h0)[0]
    pprop = np.setdiff1d(np.arange(len(hpoint)), pnoprop)
    
    Ks[pnoprop] = 1
    Kr[pnoprop] = 1
    alphab[pnoprop] = Dir0[pnoprop]
    Hb[pnoprop] = Hs0[pnoprop]
    
    return Hb, alphab