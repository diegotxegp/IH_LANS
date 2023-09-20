import numpy as np

from Funciones_modelo.snell_shoalrefrot import snell_shoalrefrot

def millerydean_dy0(tstab, tcent, Adean, Kacr0, Kero0, gammab, DYN, dinperf, t, Hberm, dt, PCALC, calcularotura, nbati, refNMM, gamma):
    # MODELO DE MILLER Y DEAN CENTRADO Y ESTABILIZADO para añadir directamente al longshore
    # Se considera un tiempo de estabilización tstab en el que se desprecia la parte transitoria en la que se alcanza la posición de equilibrio.
    # La posición de equilibrio se calcula promediando un tiempo de centrado, a partir de la fase transitoria tstab

    # Cargamos las dinámicas brutas
    Hs0 = np.zeros((len(DYN[dinperf[0]]['Hs']), len(dinperf)))
    SS0 = np.zeros((len(DYN[dinperf[0]]['SS']), len(dinperf)))
    AT0 = np.zeros((len(DYN[dinperf[0]]['AT']), len(dinperf)))
    TP0 = np.zeros((len(DYN[dinperf[0]]['Tp']), len(dinperf)))
    DIR0 = np.zeros((len(DYN[dinperf[0]]['Dir']), len(dinperf)))
    h0 = np.zeros(len(dinperf))

    for i in range(len(dinperf)):
        Hs0[:, i] = DYN[dinperf[i]]['Hs']
        SS0[:, i] = DYN[dinperf[i]]['SS']
        AT0[:, i] = DYN[dinperf[i]]['AT']
        TP0[:, i] = DYN[dinperf[i]]['Tp']
        DIR0[:, i] = DYN[dinperf[i]]['Dir']
        h0[i] = DYN[dinperf[i]]['h0']

    # Seleccionamos los perfiles
    HsP = Hs0[:, PCALC]
    SSP = SS0[:, PCALC]
    ATP = AT0[:, PCALC]
    TpP = TP0[:, PCALC]
    DIRP = DIR0[:, PCALC]
    nbatiP = [nbati[pcalc] for pcalc in PCALC]
    hP = [h0[pcalc] for pcalc in PCALC]

    # Extraer fecha común de cálculo
    pcomun = np.where(np.isin(np.array([DYN[0]['t']]), np.array([t])).tolist()[0])[0].tolist()
    HsC = HsP[pcomun, :]
    SSC = SSP[pcomun, :]
    ATC = ATP[pcomun, :]
    DIRC = DIRP[pcomun, :]
    TPC = TpP[pcomun, :]
    HsCP = np.zeros(HsC.shape)

    # Duplicamos la parte de estabilización
    if calcularotura == 0:
        Hs_s = np.vstack((HsC[0:tstab, :], HsC))
    else:
        Hs_s = np.copy(HsCP)
        for i in range(len(PCALC)):
            print('Perfil', i, 'de', len(PCALC))
            AT = ATC[:, i]
            SS = SSC[:, i]
            hi = hP[i]
            nivel = AT + SS + refNMM + hi
            Hs = HsC[:, i]
            Dir = DIRC[:, i]
            Tp = TPC[:, i]
            Hb, _, _, _, _, _, pnoprop = snell_shoalrefrot(Hs, Tp, Dir, nivel, gamma, nbatiP[i])
            Hb[pnoprop] = 0.1  # No demasiado pequeña por posibles problemas numéricos
            Hs_s[:, i] = Hb

    SS_s = np.vstack((SSC[0:tstab, :], SSC))
    AT_s = np.vstack((ATC[0:tstab, :], ATC))

    Adean_dd = [Adean[pcalc] for pcalc in PCALC]
    Hberm_dd = [Hberm[pcalc] for pcalc in PCALC]
    Kacr = [Kacr0[pcalc] for pcalc in PCALC]
    Kero = [Kero0[pcalc] for pcalc in PCALC]

    Wast_s = ((Hs_s / gammab) / Adean_dd) ** 1.5

    Dyeq = -Wast_s * ((0.106 * Hs_s + SS_s + AT_s) / (Hberm_dd + 2 * Hs_s))
    Yeq_nsl = Dyeq
    YstMDt = np.zeros(Hs_s.shape)

    # Calculamos la media para el centrado
    Aacr = [(ka / 2 * dt) for ka in Kacr]
    Aero = [(ke / 2 * dt) for ke in Kero]

    for n in range(tstab + tcent):
        A = Aacr
        YstMD1 = (YstMDt[n, :] + A * (Yeq_nsl[n + 1, :] + Yeq_nsl[n, :] - YstMDt[n, :])) / (np.ones(len(A)) + A)
        posneg = np.where(YstMD1 < YstMDt[n, :])[0]
        # Corrección
        for i in range(len(posneg)):
            A = Aero
            YstMD1[posneg[i]] = (YstMDt[n, posneg[i]] + A[posneg[i]] * (Yeq_nsl[n + 1, posneg[i]] + Yeq_nsl[n, posneg[i]] - YstMDt[n, posneg[i]])) / (1 + A[posneg[i]])
        YstMDt[n + 1, :] = YstMD1

    yst = YstMDt[0:tstab + tcent, :]
    dy0 = np.mean(YstMDt[tstab + 1:tstab + tcent + 1, :], axis=0)

    return yst, dy0, Dyeq