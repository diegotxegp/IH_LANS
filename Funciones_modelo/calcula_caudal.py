import numpy as np

from Funciones_modelo.point_to_line_distance import point_to_line_distance

def calcula_caudal(Hbd, Dbd, wbd, dx, EA, YLTi, ACT, PERF, it, Tipo):
    xlt = [p['xon'] + p['nx'] * YLTi[idx] for idx, p in enumerate(PERF)]
    ylt = [p['yon'] + p['ny'] * YLTi[idx] for idx, p in enumerate(PERF)]
    xzr = [p['xon'] + p['nx'] * (YLTi[idx] + w) for idx, (p, w) in enumerate(zip(PERF, wbd))]
    yzr = [p['yon'] + p['ny'] * (YLTi[idx] + w) for idx, (p, w) in enumerate(zip(PERF, wbd))]

    Hcelda = 0.5 * (Hbd[1:] + Hbd[:-1])

    alpha_wave0 = (270 - 0.5 * (Dbd[1:] + Dbd[:-1])) * np.pi / 180
    alpha_shoreline = np.arctan2(np.diff(ylt), np.diff(xlt)) - np.pi / 2

    alpha = np.angle(np.exp(1j * (alpha_wave0 - alpha_shoreline)))

    Kest = np.ones_like(Hcelda)

    if EA:
        estper = 'dper' in EA
        estrig = 'drig' in EA

        if any([estper, estrig]):
            if estrig:
                indrig = EA['drig']

                for i in range(len(indrig)):
                    ACTi = ACT[indrig[i]]
                    posret = np.abs(np.array(YLTi[ACTi['PERF']]) - ACTi['YAFC']) > 0.5
                    promcelda = (posret[1:] + posret[:-1]) / 2
                    indest = ACTi['PERF'][:-1]
                    Kest[indest] = promcelda
                    alpha[ACTi['PERF'][0] - 1] = np.pi / 2 if abs(2 * alpha[ACTi['PERF'][0] - 1]) > np.pi else alpha[ACTi['PERF'][0] - 1]
                    alpha[ACTi['PERF'][-1]] = np.pi / 2 if abs(2 * alpha[ACTi['PERF'][-1]]) > np.pi else alpha[ACTi['PERF'][-1]]
            if estper:
                indper = EA['dper']

                for i in range(len(indper)):
                    ACTi = ACT[indper[i]]
                    Dirdique = ACTi['Dir'][it]
                    alpha_wavedique = (270 - Dirdique) * np.pi / 180
                    xest, yest = ACTi['X'], ACTi['Y']
                    nx0 = xest[0] - xest[1]
                    ny0 = yest[0] - yest[1]
                    mod = np.hypot(nx0, ny0)
                    angle_est = np.angle(complex(ny0 / mod, nx0 / mod))
                    alpharef = alpha_wavedique * np.pi / 180 - angle_est * np.pi / 180
                    indest = ACTi['PERF']
                    xest3 = xest[1] + 1000 * np.cos(angle_est - np.pi / 2)
                    yest3 = yest[1] + 1000 * np.sin(angle_est - np.pi / 2)
                    signoy = np.sign(yest[1] - ylt[indest[1]])
                    Ldique = point_to_line_distance([xlt[indest[1]], ylt[indest[1]]], [xest[1], yest[1]], [xest3, yest3]) * signoy
                    Lzr = wbd[indest[1]]
                    Kest[indest[0]] = 1 - 0.8 * Ldique / Lzr if alpharef < 0 else 1
                    alpha = np.angle(np.exp(1j * (alpha_wave0 - alpha_shoreline)))

    pneg = np.where(Kest < 0)
    ppos = np.where(Kest > 1)
    Kest[pneg] = 0
    Kest[ppos] = 1

    if Tipo == 'CERC':
        Q = Kest * Hcelda**2.5 * np.sin(2 * alpha)
    elif Tipo == 'KAMP':
        pass
    elif Tipo == 'CERC_GRAD':
        nslope = [p['near_slope'] for p in PERF]
        mb = (nslope[1:] + nslope[:-1]) / 2
        dHbdx = np.gradient(Hbd, dx)
        dHbdx[Kest != 1] = 0
        Q = Kest * Hcelda**2.5 * (np.sin(2 * alpha) - 1 / mb * np.cos(alpha) * dHbdx)
    elif Tipo == 'KAMP_GRAD':
        pass

    return Q, Kest