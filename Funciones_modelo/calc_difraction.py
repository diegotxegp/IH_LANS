import numpy as np

from Funciones_modelo.calc_zonasombra import calc_zonasombra

def calc_difraction(PERF, ACT, estdif, Hi, Di, wi, YLTi, it, gamma):
    xzr = [p['xon'] + p['nx'] * (wi + YLTi) for p in PERF]
    yzr = [p['yon'] + p['ny'] * (wi + YLTi) for p in PERF]
    xlt = [p['xon'] + p['nx'] * YLTi for p in PERF]
    ylt = [p['yon'] + p['ny'] * YLTi for p in PERF]
    Dall = np.empty((len(estdif), len(Hi)))
    Dall[:] = np.nan
    Kdall = np.ones((len(estdif), len(Hi)))

    for i in range(len(estdif)):
        psombra, pcesion = calc_zonasombra(PERF, ACT[estdif[i]], YLTi, wi, it)

        ACTi = ACT[estdif[i]]
        if ACTi['Tipo'] == 'dper' and psombra:
            if psombra['d1']:
                APx = xzr[psombra['d1']] - ACTi['X'][-1]
                APy = yzr[psombra['d1']] - ACTi['Y'][-1]
                AP = np.angle(APx + APy * 1j, deg=True)
                deltaP = np.degrees((270 - ACTi['Dir'][it]) - AP)
                kdsombra = 0.71 - 0.0093 * abs(deltaP) + 0.000025 * deltaP**2
                Kdall[i, psombra['d1']] = kdsombra

                if pcesion['d1']['p1']:
                    ASx = xzr[pcesion['d1']['p1']] - ACTi['X'][-1]
                    ASy = yzr[pcesion['d1']['p1']] - ACTi['Y'][-1]
                    AS = np.angle(ASx + ASy * 1j, deg=True)
                    deltaS = np.degrees((270 - ACTi['Dir'][it]) - AS)
                    posmax = np.argmax(np.abs(deltaS))
                    kdcesion = 0.71 + 0.29 * np.sin(90 * deltaS / deltaS[posmax])
                    Kdall[i, pcesion['d1']['p1']] = kdcesion

                nbati_cesion = [PERF[psombra['d1'][j]]['nbati'] for j in range(len(psombra['d1']))]
                nbation = 270 - [PERF[j]['nbati'] for j in range(len(PERF))]
                alpha_b = np.degrees((270 - Di) - nbation)
                nxest = ACTi['X'][0] - ACTi['X'][1]
                nyest = ACTi['Y'][0] - ACTi['Y'][1]
                modest = np.hypot(nxest, nyest)
                alphaest = np.angle(nxest + nyest * 1j, deg=True)
                alpha_s = np.degrees((270 - ACTi['Dir'][it]) - alphaest)
                condicion = 0.5 * (np.tan(np.radians(alpha_s)) + np.tan(np.radians(0.88 * alpha_b[ACTi['PERF']])))
                
                nxpbs = xzr[psombra['d1']] - ACTi['X'][1]
                nypbs = yzr[psombra['d1']] - ACTi['Y'][1]
                modpbs = np.hypot(nxpbs, nypbs)
                apbs = np.angle(nxpbs + nypbs * 1j, deg=True)
                anglepbs = np.degrees(np.angle(np.exp(1j * (apbs - alphaest)), deg=True))
                condpbs = np.tan(np.radians(anglepbs))
                correg = np.where(np.abs(condpbs) < np.abs(condicion))[0]
                kdsombra1 = Kdall[i, psombra['d1']]
                kdcesion1 = Kdall[i, pcesion['d1']['p1']]
                Dcesion = 270 - (nbation[psombra['d1']] + alpha_b[psombra['d1']] * kdsombra1**0.375)
                Dcesion[correg] = 270 - (nbation[psombra['d1'][correg]] + alpha_b[psombra['d1'][correg]] * kdsombra1[correg]**0.375
                * (np.tan(np.radians(anglepbs[correg])) / condicion))

                Dall[i, psombra['d1']] = Dcesion
                if pcesion['d1']['p1']:
                    Dall[i, pcesion['d1']['p1']] = 270 - (nbation[pcesion['d1']['p1']] + alpha_b[pcesion['d1']['p1']] * kdcesion1**0.375)
        
        elif ACTi['Tipo'] == 'dpar' and psombra:
            xzr = [p['xon'] + p['nx'] * (wi + p['yc']) for p in PERF]
            yzr = [p['yon'] + p['ny'] * (wi + p['yc']) for p in PERF]
            if psombra['d1'] or psombra['d2']:
                Kdpar0 = np.zeros((2, len(Hi)))
                for idif in range(2):
                    if psombra[f'd{idif + 1}']:
                        APx = xzr[psombra[f'd{idif + 1}']] - ACTi['X'][idif]
                        APy = yzr[psombra[f'd{idif + 1}']] - ACTi['Y'][idif]
                        AP = np.angle(APx + APy * 1j, deg=True)
                        deltaP = np.degrees((270 - ACTi['Dir'][it]) - AP)
                        kdsombra = 0.71 - 0.0093 * abs(deltaP) + 0.000025 * deltaP**2
                        Kdpar0[idif, psombra[f'd{idif + 1}']] = kdsombra

                Dall[i, :] = 270 - (ACTi['nbati'] + (Kdpar0[0, :] * Kdpar0[1, :])**0.375 * (np.tan(np.radians(90 * (Kdpar0[1, :] - Kdpar0[0, :])))))

        elif ACTi['Tipo'] == 'oper':
            Dall[i, :] = 270 - ACTi['nbati'] - (ACTi['Disp'] / 100) * Hi

    return Dall, Kdall