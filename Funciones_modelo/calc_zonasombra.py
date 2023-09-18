import numpy as np
import matplotlib.pyplot as plt

from Generales.anglecalc import anglecalc
from CircSat.circ_dist import circ_dist

def calc_zonasombra(PERF, ACTi, YLTi, wi, it, ploteo):
    # Calcula las zonas de sombra y de cesión de estructuras

    # ACTi = ACT[0]
    # YLTi = YLT[it, :]

    # Aplicamos criterio de 2.5 L de onda para la zona de cesión
    tipo = ACTi['Tipo']
    xlc = [perf['xon'] + perf['nx'] * YLTi for perf in PERF]
    ylc = [perf['yon'] + perf['ny'] * YLTi for perf in PERF]

    xzr = [perf['xon'] + perf['nx'] * (YLTi + wi) for perf in PERF]
    yzr = [perf['yon'] + perf['ny'] * (YLTi + wi) for perf in PERF]

    # Pintamos estructura
    if tipo == 'dper':
        # Miramos si tenemos que hacer cálculo al comprobar si existe
        # intersección entre la línea de rompientes y la estructura
        pafect = ACTi['PERF']
        xi, _ = plt.polyxpoly(ACTi['X'], ACTi['Y'], xzr, yzr)
        
        # Evitar problemas de inestabilidades.
        xi1, _ = plt.polyxpoly([PERF[pafect[0]]['xon'], PERF[pafect[0]]['xon']], ACTi['Y'], xzr, yzr)
        xi2, _ = plt.polyxpoly([PERF[pafect[1]]['xon'], PERF[pafect[1]]['xon']], ACTi['Y'], xzr, yzr)
        
        if not xi:  # and not xi1 and not xi2
            psombra = []
            pcesion = []
        else:
            xest = ACTi['X'][-1]
            yest = ACTi['Y'][-1]
            zest = ACTi['Z'][-1]
            Dest = ACTi['Dir'][it]
            Test = ACTi['Tp'][it]  # Aproximamos a indefinidas

            # L0est, _ = waveguo(Test, abs(zest))
            L0est = 1.56 * Test ** 2
            Dproy = 50 * np.hypot(ACTi['X'][1] - ACTi['X'][0], ACTi['Y'][1] - ACTi['Y'][0])

            # Proyectamos sobre la LC desde la estructura según la dirección del oleaje
            xptemp = xest + np.cos(np.radians(270 - Dest)) * Dproy
            yptemp = yest + np.sin(np.radians(270 - Dest)) * Dproy

            # Calculamos la intersección con LC
            xpint, ypint = plt.polyxpoly([xest, xptemp], [yest, yptemp], xlc, ylc)

            if not xpint:
                Dlc = np.min(np.hypot(np.array(xlc) - xest, np.array(ylc) - yest))
                xpint = xest + np.cos(np.radians(270 - Dest)) * Dlc
                ypint = yest + np.sin(np.radians(270 - Dest)) * Dlc
            else:
                xpint = xpint[0]
                ypint = ypint[0]

            # Calculamos el polígono que recoge los puntos afectados
            xpolydif = [xest, xpint, ACTi['X'][0], xest]
            ypolydif = [yest, ypint, ACTi['Y'][0], yest]

            # Extraemos los puntos interiores
            psombra_d1 = [i for i, (x, y) in enumerate(zip(xzr, yzr)) if plt.pnpoly(x, y, xpolydif, ypolydif)]

            # Calculamos zona de cesión
            # Ángulo estructura
            sinest = ACTi['Y'][0] - ACTi['Y'][1]
            cosest = ACTi['X'][0] - ACTi['X'][1]
            angest = anglecalc(sinest, cosest)

            # Ángulo estructura-oleaje incidente
            angdif = circ_dist((270 - Dest) * np.pi / 180, angest * np.pi / 180) * 180 / np.pi

            if angdif > 0:
                # Olas siguen dirección de la LC (izquierda-derecha)
                nxest = np.cos(np.radians(angest + 90))
                nyest = np.sin(np.radians(angest + 90))
            else:
                # Olas en dirección contraria a la LC
                nxest = np.cos(np.radians(angest - 90))
                nyest = np.sin(np.radians(angest - 90))

            pestcesx = xest + nxest * L0est * 2.5
            pestcesy = yest + nyest * L0est * 2.5

            # Proyectamos sobre la LC desde la estructura según la dirección del oleaje
            xptemp2 = pestcesx + np.cos(np.radians(270 - Dest)) * Dproy
            yptemp2 = pestcesy + np.sin(np.radians(270 - Dest)) * Dproy

            # Calculamos la intersección con LC
            xpint2, ypint2 = plt.polyxpoly([pestcesx, xptemp2], [pestcesy, yptemp2], xlc, ylc)

            if not xpint2:
                Dlc = np.min(np.hypot(np.array(xlc) - xest, np.array(ylc) - yest))
                xpint2 = pestcesx + np.cos(np.radians(270 - Dest)) * Dlc
                ypint2 = pestcesy + np.sin(np.radians(270 - Dest)) * Dlc
            else:
                xpint2 = xpint2[0]
                ypint2 = ypint2[0]

            # Calculamos el polígono que recoge los puntos afectados
            xpolyces = [xest, xpint2, xpint, xest]
            ypolyces = [yest, ypint2, ypint, yest]

            pcesion_d1_p1 = [i for i, (x, y) in enumerate(zip(xzr, yzr)) if plt.pnpoly(x, y, xpolyces, ypolyces)]

            # Plot de verificación
            if ploteo:
                # Pintamos
                plt.figure()
                for perf in PERF:
                    plt.plot([perf['xon'], perf['xof']], [perf['yon'], perf['yof']], 'k:')
                plt.plot(xlc, ylc, 'o')
                plt.plot(xzr, yzr, '*')
                plt.plot(xpint, ypint, '*')
                plt.plot(ACTi['X'], ACTi['Y'], 'k-', linewidth=5)
                plt.quiver(ACTi['X'][-1], ACTi['Y'][-1], np.cos(np.radians(270 - Dest)), np.sin(np.radians(270 - Dest)), autoscale=200)
                plt.plot(xpolydif, ypolydif)
                plt.plot([xzr[i] for i in psombra_d1], [yzr[i] for i in psombra_d1], 'rd', markersize=10)

                # Pintamos cesión
                plt.plot(pestcesx, pestcesy, '*')
                plt.quiver(pestcesx, pestcesy, np.cos(np.radians(270 - Dest)), np.sin(np.radians(270 - Dest)), autoscale=200)
                plt.plot(xpolyces, ypolyces)
                plt.plot([xzr[i] for i in pcesion_d1_p1], [yzr[i] for i in pcesion_d1_p1], 'gd', markersize=10)

                # Pintamos para figura
                # plt.fill(xpolyces, ypolyces, color=[54, 237, 17] / 255, facealpha=0.5, edgecolor=[1, 1, 1], linewidth=2)
                # plt.fill(xpolydif, ypolydif, color=[17, 107, 237] / 255, facealpha=0.5, edgecolor=[1, 1, 1], linewidth=2)


    elif tipo == 'dperpar':
        # Miramos si tenemos que hacer cálculo al comprobar si existe
        # intersección entre la línea de rompientes y la estructura
        xi, _ = plt.polyxpoly(ACTi['X'], ACTi['Y'], xzr, yzr)
        if not xi:
            psombra = []
            pcesion = []
        else:
            Dest = ACTi['Dir'][it]

            # Miramos si continuamos los cálculos en el punto intersección dique
            # perpendicular y paralelo
            if len(ACTi['X']) == 3:
                nxper0 = ACTi['X'][0] - ACTi['X'][1]
                nyper0 = ACTi['Y'][0] - ACTi['Y'][1]
                mper = np.hypot(nxper0, nyper0)
                nxper1 = nxper0 / mper
                nyper1 = nyper0 / mper

                nxpar0 = ACTi['X'][2] - ACTi['X'][1]
                nypar0 = ACTi['Y'][2] - ACTi['Y'][1]
                mpar = np.hypot(nxpar0, nypar0)
                nxpar1 = nxpar0 / mpar
                nypar1 = nypar0 / mpar

                angleper = anglecalc(nyper1, nxper1)
                anglepar = anglecalc(nypar1, nxpar1)
                deltaest = circ_dist(angleper * np.pi / 180, anglepar * np.pi / 180) * 180 / np.pi
                deltawaves = circ_dist(angleper * np.pi / 180, (270 - Dest[0]) * np.pi / 180) * 180 / np.pi

                if np.sign(deltaest) != np.sign(deltawaves):  # Afectan los dos focos de difracción
                    # Calculamos el primero igual que un dique perpendicular
                    xest1 = ACTi['X'][-2]
                    yest1 = ACTi['Y'][-2]
                    zest1 = ACTi['Z'][-2]
                    Dest1 = ACTi['Dir'][it, -2]
                    Test1 = ACTi['Tp'][it, -2]
                    L0est = 1.56 * Test ** 2

                    Dproy = 50 * np.hypot(ACTi['X'][1] - ACTi['X'][0], ACTi['Y'][1] - ACTi['Y'][0])

                    # Proyectamos sobre la LC desde la estructura según la dirección del oleaje
                    xptemp = xest1 + np.cos(np.radians(270 - Dest1)) * Dproy
                    yptemp = yest1 + np.sin(np.radians(270 - Dest1)) * Dproy

                    # Calculamos la intersección con LC
                    xpint, ypint = plt.polyxpoly([xest1, xptemp], [yest1, yptemp], xlc, ylc)

                    if not xpint:
                        Dlc = np.min(np.hypot(np.array(xlc) - xest1, np.array(ylc) - yest1))
                        xpint = xest1 + np.cos(np.radians(270 - Dest1)) * Dlc
                        ypint = yest1 + np.sin(np.radians(270 - Dest1)) * Dlc
                    else:
                        xpint = xpint[0]
                        ypint = ypint[0]

                    # Calculamos el polígono que recoge los puntos afectados
                    xpolydif = [xest1, xpint, ACTi['X'][0], xest1]
                    ypolydif = [yest1, ypint, ACTi['Y'][0], yest1]

                    # Extraemos los puntos interiores
                    psombra_d1 = [i for i, (x, y) in enumerate(zip(xzr, yzr)) if plt.pnpoly(x, y, xpolydif, ypolydif)]

                    # Calculamos zona de cesión
                    # Ángulo estructura
                    sinest = ACTi['Y'][0] - ACTi['Y'][1]
                    cosest = ACTi['X'][0] - ACTi['X'][1]
                    angest = anglecalc(sinest, cosest)

                    # Ángulo estructura-oleaje incidente
                    angdif = circ_dist((270 - Dest1) * np.pi / 180, angest * np.pi / 180) * 180 / np.pi

                    if angdif > 0:
                        # Olas siguen dirección de la LC (izquierda-derecha)
                        nxest = np.cos(np.radians(angest + 90))
                        nyest = np.sin(np.radians(angest + 90))
                    else:
                        # Olas en dirección contraria a la LC
                        nxest = np.cos(np.radians(angest - 90))
                        nyest = np.sin(np.radians(angest - 90))

                    pestcesx = xest1 + nxest * L0est * 2.5
                    pestcesy = yest1 + nyest * L0est * 2.5

                    # Proyectamos sobre la LC desde la estructura según la dirección del oleaje
                    xptemp2 = pestcesx + np.cos(np.radians(270 - Dest1)) * Dproy
                    yptemp2 = pestcesy + np.sin(np.radians(270 - Dest1)) * Dproy

                    # Calculamos la intersección con LC
                    xpint2, ypint2 = plt.polyxpoly([pestcesx, xptemp2], [pestcesy, yptemp2], xlc, ylc)

                    if not xpint2:
                        Dlc = np.min(np.hypot(np.array(xlc) - xest1, np.array(ylc) - yest1))
                        xpint2 = pestcesx + np.cos(np.radians(270 - Dest1)) * Dlc
                        ypint2 = pestcesy + np.sin(np.radians(270 - Dest1)) * Dlc
                    else:
                        xpint2 = xpint2[0]
                        ypint2 = ypint2[0]


                    # Calculamos el polígono que recoge los puntos afectados
                    xpolyces = [xest1, xpint2, xpint, xest1]
                    ypolyces = [yest1, ypint2, ypint, yest1]
                    pcesion_d1_p1 = [i for i, (x, y) in enumerate(zip(xzr, yzr)) if plt.pnpoly(x, y, xpolyces, ypolyces)]

                    # Verificamos que el rayo no corte a la estructura
                    nx1 = ACTi['X'][0] - ACTi['X'][2]
                    ny1 = ACTi['Y'][0] - ACTi['Y'][2]
                    mod = np.hypot(nx1, ny1)
                    nxm = nx1 / mod
                    nym = ny1 / mod
                    angestmo = anglecalc(nym, nxm)
                    angwaves = 270 - Dest[1]
                    delta = circ_dist(angestmo * np.pi / 180, angwaves * np.pi / 180) * 180 / np.pi

                    if delta > 0:  # El rayo corta a la estructura
                        xest2 = ACTi['X'][-1]
                        yest2 = ACTi['Y'][-1]
                        zest2 = ACTi['Z'][-1]
                        Test2 = ACTi['Tp'][it, -1]
                        L0est2 = 1.56 * Test ** 2

                        Dproy = 50 * np.hypot(ACTi['X'][1] - ACTi['X'][0], ACTi['Y'][1] - ACTi['Y'][0])

                        # Proyectamos sobre la LC desde la estructura según la dirección del oleaje
                        xptemp3 = xest2 + np.cos(np.radians(270 - Dest[1])) * Dproy
                        yptemp3 = yest2 + np.sin(np.radians(270 - Dest[1])) * Dproy

                        # Calculamos la intersección con LC
                        xpint3, ypint3 = plt.polyxpoly([xest2, xptemp3], [yest2, yptemp3], ACTi['X'][:2], ACTi['Y'][:2])

                        # Calculamos el polígono que recoge los puntos afectados
                        xpolydif2 = [xest2, xpint3, ACTi['X'][1], xest2]
                        ypolydif2 = [yest2, ypint3, ACTi['Y'][1], yest2]

                        psombra_d2 = [i for i, (x, y) in enumerate(zip(xzr, yzr)) if plt.pnpoly(x, y, xpolydif2, ypolydif2)]

                        if psombra_d2:
                            # Calculamos cesión
                            sinpar = ACTi['Y'][2] - ACTi['Y'][1]
                            cospar = ACTi['X'][2] - ACTi['X'][1]
                            angpar = anglecalc(sinpar, cospar)
                            deltaest = circ_dist(angest * np.pi / 180, angpar * np.pi / 180) * 180 / np.pi
                            angcesion = angest - np.sign(deltaest) * 90
                            Tces = ACTi['Tp'][-1]
                            Lces = 1.56 * Tces**2
                            xfces = ACTi['X'][2] + 2.5 * Lces * np.cos(np.radians(angcesion))
                            yfces = ACTi['Y'][2] + 2.5 * Lces * np.sin(np.radians(angcesion))

                            # Cortamos LC
                            xptemp4 = xfces + Dproy * np.cos(np.radians(270 - Dest[1]))
                            yptemp4 = yfces + Dproy * np.sin(np.radians(270 - Dest[1]))

                            # Calculamos la intersección con LC
                            xpint4, ypint4 = plt.polyxpoly([xfces, xptemp4], [yfces, yptemp4], xlc, ylc)

                            if not xpint4:
                                Dlc = np.min(np.hypot(np.array(xlc) - xest2, np.array(ylc) - yest2))
                                xpint4 = xfces + np.cos(np.radians(270 - Dest[1])) * Dlc
                                ypint4 = yfces + np.sin(np.radians(270 - Dest[1])) * Dlc
                            else:
                                xpint4 = xpint4[0]
                                ypint4 = ypint4[0]

                            xpolices2 = [xfces, xpint4, xpint3, xest2, xfces]
                            ypolices2 = [yfces, ypint4, ypint3, yest2, yfces]

                            pcesion_d2_p2 = [i for i, (x, y) in enumerate(zip(xzr, yzr)) if plt.pnpoly(x, y, xpolices2, ypolices2)]
                        else:
                            pcesion_d2_p2 = []

                    else:  # El rayo no corta a la estructura
                        # Calculamos el segundo punto
                        xest2 = ACTi['X'][-1]
                        yest2 = ACTi['Y'][-1]
                        zest2 = ACTi['Z'][-1]
                        Test2 = ACTi['Tp'][it, -1]
                        L0est2 = 1.56 * Test ** 2

                        Dproy = 50 * np.hypot(ACTi['X'][1] - ACTi['X'][0], ACTi['Y'][1] - ACTi['Y'][0])

                        # Proyectamos sobre la LC desde la estructura según la dirección del oleaje
                        xptemp3 = xest2 + np.cos(np.radians(270 - Dest[1])) * Dproy
                        yptemp3 = yest2 + np.sin(np.radians(270 - Dest[1])) * Dproy

                        # Calculamos la intersección con LC
                        xpint3, ypint3 = plt.polyxpoly([xest2, xptemp3], [yest2, yptemp3], xlc, ylc)

                        if not xpint3:
                            Dlc = np.min(np.hypot(np.array(xlc) - xest1, np.array(ylc) - yest1))
                            xpint3 = pestcesx + np.cos(np.radians(270 - Dest[1])) * Dlc
                            ypint3 = pestcesy + np.sin(np.radians(270 - Dest[1])) * Dlc
                        else:
                            xpint3 = xpint3[0]
                            ypint3 = ypint3[0]

                        # Calculamos el polígono que recoge los puntos afectados
                        xpolydif2 = [xest2, xpint3, ACTi['X'][:2], xest2]
                        ypolydif2 = [yest2, ypint3, ACTi['Y'][:2], yest2]

                        psombra_d2 = [i for i, (x, y) in enumerate(zip(xzr, yzr)) if plt.pnpoly(x, y, xpolydif2, ypolydif2)]

                        if psombra_d2:
                            # Calculamos cesión
                            sinpar = ACTi['Y'][2] - ACTi['Y'][1]
                            cospar = ACTi['X'][2] - ACTi['X'][1]
                            angpar = anglecalc(sinpar, cospar)
                            deltaest = circ_dist(angest * np.pi / 180, angpar * np.pi / 180) * 180 / np.pi
                            angcesion = angest - np.sign(deltaest) * 90
                            Tces = ACTi['Tp'][-1]
                            Lces = 1.56 * Tces**2
                            xfces = ACTi['X'][2] + 2.5 * Lces * np.cos(np.radians(angcesion))
                            yfces = ACTi['Y'][2] + 2.5 * Lces * np.sin(np.radians(angcesion))

                            # Cortamos LC
                            xptemp4 = xfces + Dproy * np.cos(np.radians(270 - Dest[1]))
                            yptemp4 = yfces + Dproy * np.sin(np.radians(270 - Dest[1]))

                            # Calculamos la intersección con LC
                            xpint4, ypint4 = plt.polyxpoly([xfces, xptemp4], [yfces, yptemp4], xlc, ylc)

                            if not xpint4:
                                Dlc = np.min(np.hypot(np.array(xlc) - xest2, np.array(ylc) - yest2))
                                xpint4 = xfces + np.cos(np.radians(270 - Dest[1])) * Dlc
                                ypint4 = yfces + np.sin(np.radians(270 - Dest[1])) * Dlc
                            else:
                                xpint4 = xpint4[0]
                                ypint4 = ypint4[0]

                            xpolices2 = [xfces, xpint4, xpint3, xest2, xfces]
                            ypolices2 = [yfces, ypint4, ypint3, yest2, yfces]

                            pcesion_d2_p2 = [i for i, (x, y) in enumerate(zip(xzr, yzr)) if plt.pnpoly(x, y, xpolices2, ypolices2)]
                        else:
                            pcesion_d2_p2 = []

                    if ploteo:
                        # Pintamos
                        plt.plot([PERF[i]['xon'], PERF[i]['xof'] for i in range(len(PERF))], [PERF[i]['yon'], PERF[i]['yof'] for i in range(len(PERF))], 'k:')
                        plt.plot(xlc, ylc, 'o')
                        plt.plot(xzr, yzr, '*')
                        plt.plot(xpint, ypint, '*')
                        plt.plot(ACTi['X'], ACTi['Y'], 'k-', linewidth=5)
                        plt.quiver(ACTi['X'][-2:], ACTi['Y'][-2:], np.cos(np.radians(270 - Dest)), np.sin(np.radians(270 - Dest)), autoscale=1)
                        plt.plot(xpolydif, ypolydif)
                        plt.plot(xzr[psombra_d1], yzr[psombra_d1], 'rd', markersize=10)

                        # Pintamos cesión
                        plt.plot(pestcesx, pestcesy, '*')
                        plt.quiver(pestcesx, pestcesy, np.cos(np.radians(270 - Dest[0])), np.sin(np.radians(270 - Dest[0])), autoscale=200)
                        plt.plot(xpolyces, ypolyces)
                        plt.plot(xzr[pcesion_d1_p1], yzr[pcesion_d1_p1], 'gd', markersize=10)

                        # Punto 2
                        plt.plot(xpint3, ypint3, 's')
                        plt.plot(xpolydif2, ypolydif2)
                        plt.plot(xzr[psombra_d2], yzr[psombra_d2], 'rd', markersize=10)

                        # Cesión
                        if psombra_d2:
                            plt.plot(xfces, yfces, '*')
                            plt.quiver(xfces, yfces, np.cos(np.radians(270 - Dest[1])), np.sin(np.radians(270 - Dest[1])), autoscale=200)
                            plt.plot(xpint4, ypint4, 'o')
                            plt.plot(xpolices2, ypolices2)
                            plt.plot(xzr[pcesion_d2_p2], yzr[pcesion_d2_p2], 'gd', markersize=10)


                else: # afecta el foco principal
                    xest = ACTi.X[-1]
                    yest = ACTi.Y[-1]
                    zest = ACTi.Z[-1]
                    Dest = ACTi.Dir[it, -1]
                    Test = ACTi.Tp[it, -1]
                    L0est = 1.56 * Test ** 2

                    # L0est = waveguo(Test, abs(zest))  # Descomenta si tienes una función waveguo definida

                    Dproy = 50 * np.hypot(ACTi.X[1] - ACTi.X[0], ACTi.Y[1] - ACTi.Y[0])

                    # proyectamos sobre la LC desde la estructura según la dirección del oleaje
                    xptemp = xest + np.cos(np.radians(270 - Dest)) * Dproy
                    yptemp = yest + np.sin(np.radians(270 - Dest)) * Dproy

                    # calculamos la intersección con LC
                    mask = ~np.isnan(xlc) & ~np.isnan(ylc)
                    xpint, ypint = polyxpoly([xest, xptemp], [yest, yptemp], xlc[mask], ylc[mask])

                    # calculamos el polígono que recoge los puntos afectados
                    xpolydif = [xest, xpint, ACTi.X[0:2], xest]
                    ypolydif = [yest, ypint, ACTi.Y[0:2], yest]

                    # extraemos los puntos interiores
                    psombra_d2 = np.where(inpolygon(xzr, yzr, xpolydif, ypolydif))

                    # calculamos zona de cesión
                    # angulo estructura
                    sinest = ACTi.Y[0] - ACTi.Y[1]
                    cosest = ACTi.X[0] - ACTi.X[1]
                    angest = np.degrees(np.arctan2(sinest, cosest))

                    if deltaest < 0:
                        # olas siguen dirección de la LC (izquierda-derecha)
                        nxest = np.cos(np.radians(angest + 90))
                        nyest = np.sin(np.radians(angest + 90))
                    else:
                        # olas en dirección contraria a la LC
                        nxest = np.cos(np.radians(angest - 90))
                        nyest = np.sin(np.radians(angest - 90))

                    # Calculamos la posición de cesión
                    pestcesx = xest + nxest * L0est * 2.5
                    pestcesy = yest + nyest * L0est * 2.5

                    # Proyectamos sobre la LC desde la estructura según la dirección del oleaje
                    xptemp2 = pestcesx + np.cos(np.radians(270 - Dest)) * Dproy
                    yptemp2 = pestcesy + np.sin(np.radians(270 - Dest)) * Dproy

                    # Calculamos la intersección con LC
                    mask = ~np.isnan(xlc) & ~np.isnan(ylc)
                    xpint2, ypint2 = polyxpoly([pestcesx, xptemp2], [pestcesy, yptemp2], xlc[mask], ylc[mask])

                    # Calculamos el polígono que recoge los puntos afectados
                    xpolyces = [xest, xpint2, xpint, xest]
                    ypolyces = [yest, ypint2, ypint, yest]

                    # Encontramos los puntos afectados en la zona de cesión
                    pcesion_d2_p2 = np.where(inpolygon(xzr, yzr, xpolyces, ypolyces))

                    if ploteo:
                        # Pintamos los resultados si es necesario
                        plt.figure()
                        for i in range(len(PERF)):
                            plt.plot([PERF[i].xon, PERF[i].xof], [PERF[i].yon, PERF[i].yof], 'k:')
                        plt.plot(xlc, ylc, 'o')
                        plt.plot(xzr, yzr, '*')
                        plt.plot(xpint, ypint, '*')
                        plt.plot(ACTi.X, ACTi.Y, 'k-', linewidth=5)
                        plt.quiver(ACTi.X[-1], ACTi.Y[-1], np.cos(np.radians(270 - Dest)), np.sin(np.radians(270 - Dest)), scale=200)
                        plt.plot(xpolydif, ypolydif)
                        plt.plot(xzr[psombra_d2], yzr[psombra_d2], 'rd', markersize=10)
                        # Pintamos cesión
                        plt.plot(pestcesx, pestcesy, '*')
                        plt.quiver(pestcesx, pestcesy, np.cos(np.radians(270 - Dest)), np.sin(np.radians(270 - Dest)), scale=200)
                        plt.plot(xpolyces, ypolyces)
                        plt.plot(xzr[pcesion_d2_p2], yzr[pcesion_d2_p2], 'gd', markersize=10)
                        plt.show()


    elif tipo == 'dpar':
        # Calculamos xi, la intersección entre la línea de rompientes y la estructura
        # Consideramos ZR como el valor inicial para evitar inestabilidades
        # xzr = [PERF.xon] + [PERF.nx] * ([PERF.yc] + wi)
        # yzr = [PERF.yon] + [PERF.ny] * ([PERF.yc] + wi)

        pref = ACTi.PERF[0]
        xi, _ = polyxpoly([PERF[pref].xon, ACTi.X, PERF[pref].xon], [PERF[pref].yon, ACTi.Y, PERF[pref].yon], xzr, yzr)

        if len(xi) == 0:
            psombra = {}
            pcesion = {}
        else:
            if ploteo:
                plt.figure()
                # for i in range(len(PERF)):
                #     plt.plot([PERF[i].xon, PERF[i].xof], [PERF[i].yon, PERF[i].yof], 'k:')
                # plt.plot(xlc, ylc, 'o')
                # plt.plot(xzr, yzr, '*')

            for ip in range(2):
                # Punto de difracción 1
                xest1 = ACTi.X[0]
                yest1 = ACTi.Y[0]
                # Punto de difracción 2
                xest2 = ACTi.X[1]
                yest2 = ACTi.Y[1]
                zest = ACTi.Z[ip]
                # Dirección del punto ip
                Dest = ACTi.Dir(it, ip)
                Tpest = ACTi.Tp(it, ip)
                L0est = 1.56 * Tpest ** 2
                # [L0est, _] = waveguo(Tpest, abs(zest))
                # Distancia de proyección
                Dproy = 50 * hypot(ACTi.X[1] - ACTi.X[0], ACTi.Y[1] - ACTi.Y[0])
                # Rayo desde el punto 1
                xptemp1 = xest1 + np.cos(np.radians(270 - Dest)) * Dproy
                yptemp1 = yest1 + np.sin(np.radians(270 - Dest)) * Dproy
                # Rayo desde el punto 2
                xptemp2 = xest2 + np.cos(np.radians(270 - Dest)) * Dproy
                yptemp2 = yest2 + np.sin(np.radians(270 - Dest)) * Dproy
                # Intersección con LC 1
                xpint1, ypint1 = polyxpoly([xest1, xptemp1], [yest1, yptemp1], xlc, ylc)
                # Intersección con LC 2
                xpint2, ypint2 = polyxpoly([xest2, xptemp2], [yest2, yptemp2], xlc, ylc)

                if len(xpint1) == 0:
                    Dlc = min(hypot(xlc - xest1, ylc - yest1))
                    xpint1 = xest1 + np.cos(np.radians(270 - Dest)) * Dlc
                    ypint1 = yest1 + np.sin(np.radians(270 - Dest)) * Dlc
                else:
                    xpint1 = xpint1[0]
                    ypint1 = ypint1[0]

                if len(xpint2) == 0:
                    Dlc = min(hypot(xlc - xest2, ylc - yest2))
                    xpint2 = xest2 + np.cos(np.radians(270 - Dest)) * Dlc
                    ypint2 = yest2 + np.sin(np.radians(270 - Dest)) * Dlc
                else:
                    xpint2 = xpint2[0]
                    ypint2 = ypint2[0]

                # Calculamos el polígono que recoge los puntos afectados
                xpolydif = [xest1, xpint1, xpint2, xest2, xest1]
                ypolydif = [yest1, ypint1, ypint2, yest2, yest1]

                # Extraemos los puntos interiores
                psombraT = np.where(inpolygon(xzr, yzr, xpolydif, ypolydif))

                # Si no hay intersección entre LC y estructura, hacemos continua la sombra
                intx, _ = polyxpoly(ACTi.X, ACTi.Y, xlc, ylc)
                if len(intx) == 0:
                    psombra['d' + str(ip)] = list(range(psombraT[0], psombraT[-1] + 1))
                else:
                    psombra['d' + str(ip)] = list(psombraT)

                # Calculamos zona de cesión
                # Ángulo de la estructura
                sinest = ACTi.Y[0] - ACTi.Y[1]
                cosest = ACTi.X[0] - ACTi.X[1]
                angest = anglecalc(sinest, cosest)

                # Ángulo entre la estructura y el oleaje incidente
                angdif = circ_dist((270 - Dest) * np.pi / 180, angest * np.pi / 180) * 180 / np.pi

                pestcesx1 = ACTi.X[0] + np.cos(np.radians(angest)) * L0est * 2.5
                pestcesy1 = ACTi.Y[0] + np.sin(np.radians(angest)) * L0est * 2.5
                pestcesx2 = ACTi.X[1] - np.cos(np.radians(angest)) * L0est * 2.5
                pestcesy2 = ACTi.Y[1] - np.sin(np.radians(angest)) * L0est * 2.5

                # Proyectamos sobre la LC desde la estructura según la dirección del oleaje
                xptemp3 = pestcesx1 + np.cos(np.radians(270 - Dest)) * Dproy
                yptemp3 = pestcesy1 + np.sin(np.radians(270 - Dest)) * Dproy
                xptemp4 = pestcesx2 + np.cos(np.radians(270 - Dest)) * Dproy
                yptemp4 = pestcesy2 + np.sin(np.radians(270 - Dest)) * Dproy

                # Calculamos la intersección con LC 1
                xpint3, ypint3 = polyxpoly([pestcesx1, xptemp3], [pestcesy1, yptemp3], xlc, ylc)
                # Calculamos la intersección con LC 2
                xpint4, ypint4 = polyxpoly([pestcesx2, xptemp4], [pestcesy2, yptemp4], xlc, ylc)

                if len(xpint3) == 0:
                    Dlc = min(hypot(xlc - pestcesx1, ylc - pestcesy1))
                    xpint3 = pestcesx1 + np.cos(np.radians(270 - Dest)) * Dlc
                    ypint3 = pestcesy1 + np.sin(np.radians(270 - Dest)) * Dlc
                else:
                    xpint3 = xpint3[0]
                    ypint3 = ypint3[0]

                if len(xpint4) == 0:
                    Dlc = min(hypot(xlc - pestcesx2, ylc - pestcesy2))
                    xpint4 = pestcesx2 + np.cos(np.radians(270 - Dest)) * Dlc
                    ypint4 = pestcesy2 + np.sin(np.radians(270 - Dest)) * Dlc
                else:
                    xpint4 = xpint4[0]
                    ypint4 = ypint4[0]

                # Calculamos el polígono que recoge los puntos afectados
                xpolyces1 = [xest1, xpint1, xpint3, xest1]
                ypolyces1 = [yest1, ypint1, ypint3, yest1]

                pcesion['d' + str(ip)] = {}
                pcesion['d' + str(ip)]['p1'] = np.where(inpolygon(xzr, yzr, xpolyces1, ypolyces1))

                xpolyces2 = [xest2, xpint2, xpint4, xest2]
                ypolyces2 = [yest2, ypint2, ypint4, yest2]

                pcesion['d' + str(ip)]['p2'] = np.where(inpolygon(xzr, yzr, xpolyces2, ypolyces2))

                if ploteo:
                    plt.plot([xpint1, xpint2], [ypint1, ypint2], '*')
                    plt.plot(ACTi.X, ACTi.Y, 'k-', linewidth=5)
                    plt.quiver(ACTi.X, ACTi.Y, np.cos(np.radians(270 - Dest)) * np.array([1, 1]), np.sin(np.radians(270 - Dest)) * np.array([1, 1]), scale=1)
                    plt.plot(xpolydif, ypolydif, 'd-')
                    plt.plot(xzr[psombra['d' + str(ip)]], yzr[psombra['d' + str(ip)]], 'rd', markersize=10)
                    # Cesion
                    plt.plot([pestcesx1, pestcesx2], [pestcesy1, pestcesy2], '*')
                    plt.quiver([pestcesx1, pestcesx2], [pestcesy1, pestcesy2], np.cos(np.radians(270 - Dest)) * np.array([1, 1]), np.sin(np.radians(270 - Dest)) * np.array([1, 1]), scale=1)
                    plt.plot(xpolyces1, ypolyces1)
                    plt.plot(xpolyces2, ypolyces2, 'd-')
                    plt.plot(xzr[pcesion['d' + str(ip)]['p1']], yzr[pcesion['d' + str(ip)]['p1']], 'gd', markersize=10)
                    plt.plot(xzr[pcesion['d' + str(ip)]['p2']], yzr[pcesion['d' + str(ip)]['p2']], 'gd', markersize=10)
        # plt.fill(xpolydif, ypolydif, color=[17, 107, 237] / 255, alpha=0.25, edgecolor=[1, 1, 1], linewidth=2)
        # plt.fill(xpolyces1, ypolyces1, color=[54, 237, 17] / 255, alpha=0.25, edgecolor=[1, 1, 1], linewidth=2)
        # plt.fill(xpolyces2, ypolyces2, color=[54, 237, 17] / 255, alpha=0.25, edgecolor=[1, 1, 1], linewidth=2)
        plt.show()