import numpy as np
from Funciones_modelo.calc_zonasombra import calc_zonasombra

def reparto_bypass(Q, EA, ACT, PERF, YLTi, wi, it):
    Qc = Q.copy()
    if 'dper' in EA:
        ndper = len(EA['dper'])
        for i in range(ndper):
            ACTi = ACT[EA['dper'][i]]
            psombra, pcesion = calc_zonasombra(PERF, ACTi, YLTi, wi, it, 0)
            Qbypass = Q[ACTi['PERF'][0]]
            if Qbypass != 0:
                cede = [psombra['d1'], pcesion['d1']['p1']] + np.sign(Qbypass)
                Qc[cede] = Q[cede] + Qbypass / len(cede)
                Qc[ACTi['PERF'][0]] = Qbypass / len(cede)
    elif 'dperpar' in EA:
        pass
    else:
        pass
    return Qc