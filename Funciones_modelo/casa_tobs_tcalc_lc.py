import math
import numpy as np

from Funciones_modelo.ensambla_lc import ensambla_lc

def casa_tobs_tcalc_lc(t, PERF, data_asim_l, data_asim_c, data_asim_lc, **kwargs):
    if data_asim_l:
        sigmaK = kwargs.get('sigmaK', 0.5)
        kcerc = kwargs.get('kcerc', [])
    if data_asim_c:
        kacr = kwargs.get('kacr', [])
        kero = kwargs.get('kero', [])

    DA = []
    for ip in range(len(PERF)):
        tobs = PERF[ip]['date_obs']
        posin_t = np.array([i for i, obs in enumerate(tobs) if t[0] < obs <= t[-1]])
        tasim = [tobs[i] for i in posin_t]
        nasim = []
        for i in range(len(tasim)):
            tasimi = tasim[i]
            possup = [j for j, time in enumerate(t) if time >= tasimi]
            nasim.append(min(possup))

        # escogemos sólo valores únicos
        pos_unique = []
        nasim_unique = []
        for pos, value in enumerate(sorted(set(nasim))):
            pos_unique.append(pos)
            nasim_unique.append(value)

        DA_entry = {
            'tasim': np.array([element[0] for element in tasim])[pos_unique],
            'nasim': np.array(nasim_unique),
            'itmax': max(nasim_unique),
            'pos': 1,
            'stop': 0,
            'pos_l': 1,
            'stop_l': 0,
            'pos_c': 1,
            'stop_c': 0,
            'pos_lc': 1,
            'stop_lc': 0
        }

        if data_asim_l:
            DA_entry['Ylt'] = [PERF[ip]['Y_obs_lt'][i] for i in posin_t[pos_unique]]
            DA_entry['kcerc0'] = kcerc[ip]
            DA_entry['sigmaK'] = sigmaK
            P0 = PERF[ip]['rP0'] ** 2
            Q = PERF[ip]['rQ'] ** 2
            P0[1, 1] = (1 / sigmaK * math.log((kcerc[ip] + PERF[ip]['rP0'][1, 1]) / kcerc[ip])) ** 2
            Q[1, 1] = (1 / sigmaK * math.log((kcerc[ip] + PERF[ip]['rQ'][1, 1]) / kcerc[ip])) ** 2
            DA_entry['P0'] = P0
            DA_entry['Q'] = Q
            DA_entry['R'] = PERF[ip]['R']

        if data_asim_c:
            DA_entry['Yct'] = [PERF[ip]['Y_obs_ct'][i] for i in posin_t[pos_unique]]
            DA_entry['kacr0'] = kacr[ip]
            DA_entry['kero0'] = kero[ip]
            DA_entry['Pacr0_c'] = PERF[ip]['rPacr0_c'] ** 2
            DA_entry['Pero0_c'] = PERF[ip]['rPero0_c'] ** 2
            DA_entry['Qacr_c'] = PERF[ip]['rQacr_c'] ** 2
            DA_entry['Qero_c'] = PERF[ip]['rQero_c'] ** 2
            DA_entry['R_c'] = PERF[ip]['R_c']

        if data_asim_lc:
            Ylt = [PERF[ip]['Y_obs_lt'][i] for i in posin_t[pos_unique]]
            Yct = [PERF[ip]['Y_obs_ct'][i] for i in posin_t[pos_unique]]
            DA_entry['YY'] = [lt + ct for lt, ct in zip(Ylt, Yct)]
            DA_entry['kcerc0'] = kcerc[ip]
            DA_entry['sigmaK'] = sigmaK
            DA_entry['kacr0'] = kacr[ip]
            DA_entry['kero0'] = kero[ip]
            P0 = PERF[ip]['rP0'] ** 2
            Q = PERF[ip]['rQ'] ** 2
            Pl = P0.copy()
            Ql = Q.copy()
            Pl[1, 1] = (1 / sigmaK * math.log((kcerc[ip] + PERF[ip]['rP0'][1, 1]) / kcerc[ip])) ** 2
            Ql[1, 1] = (1 / sigmaK * math.log((kcerc[ip] + PERF[ip]['rQ'][1, 1]) / kcerc[ip])) ** 2
            Pacr0_c = PERF[ip]['rPacr0_c'] ** 2
            Pero0_c = PERF[ip]['rPero0_c'] ** 2
            Qacr_c = PERF[ip]['rQacr_c'] ** 2
            Qero_c = PERF[ip]['rQero_c'] ** 2
            DA_entry['Pero0_lc'] = ensambla_lc(Pl, Pero0_c)
            DA_entry['Pacr0_lc'] = ensambla_lc(Pl, Pacr0_c)
            DA_entry['Qero_lc'] = ensambla_lc(Ql, Qero_c)
            DA_entry['Qacr_lc'] = ensambla_lc(Ql, Qacr_c)
            DA_entry['R_lc'] = PERF[ip]['R_c'] + PERF[ip]['R']

        DA.append(DA_entry)
        print("Nuevo DA_entry")
        print(ip)
        print((len(PERF)))

    DA = np.array(DA)

    return DA