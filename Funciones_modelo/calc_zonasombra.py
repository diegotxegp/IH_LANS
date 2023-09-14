import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from shapely.geometry import Polygon, Point
from shapely.ops import cascaded_union

def calc_zonasombra(PERF, ACTi, YLTi, wi, it, ploteo):
    # Calcula las zonas de sombra y de cesión de estructuras
    # Definimos variables previamente

    # ACTi = ACT(1);
    # YLTi = YLT[it, :];

    # Aplicamos criterio de 2.5 L de onda para la zona de cesión
    tipo = ACTi.Tipo
    xlc = [p['xon'] + p['nx'] * YLTi for p in PERF]
    ylc = [p['yon'] + p['ny'] * YLTi for p in PERF]

    xzr = [p['xon'] + p['nx'] * (YLTi + wi) for p in PERF]
    yzr = [p['yon'] + p['ny'] * (YLTi + wi) for p in PERF]

    # Pintamos la estructura
    if tipo == 'dper':
        # Miramos si tenemos que hacer cálculo al comprobar si existe
        # intersección entre la línea de rompientes y la estructura
        pafect = ACTi.PERF
        xi, _ = polyxpoly(ACTi.X, ACTi.Y, xzr, yzr)
        # Evitar problemas de inestabilidades.
        xi1, _ = polyxpoly([PERF[pafect[0]]['xon'], PERF[pafect[0]]['xon']], ACTi.Y, xzr, yzr)
        xi2, _ = polyxpoly([PERF[pafect[1]]['xon'], PERF[pafect[1]]['xon']], ACTi.Y, xzr, yzr)
        if len(xi) == 0:  # and len(xi1) == 0 and len(xi2) == 0:
            psombra = []
            pcesion = []
        else:
            xest = ACTi.X[-1]
            yest = ACTi.Y[-1]
            zest = ACTi.Z[-1]
            Dest = ACTi.Dir[it]
            Test = ACTi.Tp[it]  # Aproximamos a indefinidas
            # L0est, _ = waveguo(Test, abs(zest))
            L0est = 1.56 * Test ** 2
            Dproy = 50 * np.hypot(ACTi.X[1] - ACTi.X[0], ACTi.Y[1] - ACTi.Y[0])
            # Proyectamos sobre la LC desde la estructura según la dirección del oleaje
            xptemp = xest + np.cos(np.deg2rad(270 - Dest)) * Dproy
            yptemp = yest + np.sin(np.deg2rad(270 - Dest)) * Dproy
            # Calculamos la intersección con LC
            xpint, ypint = polyxpoly([xest, xptemp], [yest, yptemp], xlc, ylc)
            if len(xpint) == 0:
                Dlc = min(np.hypot(np.array(xlc) - xest, np.array(ylc) - yest))
                xpint = xest + np.cos(np.deg2rad(270 - Dest)) * Dlc
                ypint = yest + np.sin(np.deg2rad(270 - Dest)) * Dlc
            else:
                xpint = xpint[0]
                ypint = ypint[0]
            # Calculamos el polígono que recoge los puntos afectados
            xpolydif = [xest, xpint, ACTi.X[0], xest]
            ypolydif = [yest, ypint, ACTi.Y[0], yest]
            # Extraemos los puntos interiores
            psombra_d1 = []
            for i in range(len(xzr)):
                point = Point(xzr[i], yzr[i])
                polygon = Polygon([(x, y) for x, y in zip(xpolydif, ypolydif)])
                if polygon.contains(point):
                    psombra_d1.append(i)
            psombra = {'d1': psombra_d1}
            # Calculamos zona de cesión
            # Ángulo estructura
            sinest = ACTi.Y[0] - ACTi.Y[1]
            cosest = ACTi.X[0] - ACTi.X[1]
            angest = anglecalc(sinest, cosest)
            # Ángulo estructura-oleaje incidente
            angdif = circ_dist((270 - Dest) * np.pi / 180, angest * np.pi / 180) * 180 / np.pi
            if angdif > 0:
                # Olas siguen dirección de la LC (izquierda-derecha)
                nxest = np.cos(np.deg2rad(angest + 90))
                nyest = np.sin(np.deg2rad(angest + 90))
            else:
                # Olas en dirección contraria a la LC
                nxest = np.cos(np.deg2rad(angest - 90))
                nyest = np.sin(np.deg2rad(angest - 90))
            pestcesx = xest + nxest * L0est * 2.5
            pestcesy = yest + nyest * L0est * 2.5
            # Proyectamos sobre la LC desde la estructura según la dirección del oleaje
            xptemp2 = pestcesx + np.cos(np.deg2rad(270 - Dest)) * Dproy
            yptemp2 = pestcesy + np.sin(np.deg2rad(270 - Dest)) * Dproy
            # Calculamos la intersección con LC
            xpint2, ypint2 = polyxpoly([pestcesx, xptemp2], [pestcesy, yptemp2], xlc, ylc)
            if len(xpint2) == 0:
                Dlc = min(np.hypot(np.array(xlc) - xest, np.array(ylc) - yest))
                xpint2 = pestcesx + np.cos(np.deg2rad(270 - Dest)) * Dlc
                ypint2 = pestcesy + np.sin(np.deg2rad(270 - Dest)) * Dlc
            else:
                xpint2 = xpint2[0]
                ypint2 = ypint2[0]
            # Calculamos el polígono que recoge los puntos afectados
            xpolyces = [xest, xpint2, xpint, xest]
            ypolyces = [yest, ypint2, ypint, yest]
            pcesion_d1_p1 = []
            for i in range(len(xzr)):
                point = Point(xzr[i], yzr[i])
                polygon = Polygon([(x, y) for x, y in zip(xpolyces, ypolyces)])
                if polygon.contains(point):
                    pcesion_d1_p1.append(i)
            pcesion = {'d1': {'p1': pcesion_d1_p1}}
            # Plot de verificación
            if ploteo:
                # Pintamos
                plt.figure()
                for i in range(len(PERF)):
                    plt.plot([PERF[i]['xon'], PERF[i]['xof']], [PERF[i]['yon'], PERF[i]['yof']], 'k:')
                plt.plot(xlc, ylc, 'o')
                plt.plot(xzr, yzr, '*')
                plt.plot(xpint, ypint, '*')
                plt.plot(ACTi.X, ACTi.Y, 'k-', linewidth=5)
                plt.quiver(ACTi.X[-1], ACTi.Y[-1], np.cos(np.deg2rad(270 - Dest)), np.sin(np.deg2rad(270 - Dest)), autoscale=200)
                plt.plot(xpolydif, ypolydif)
                plt.plot(np.array(xzr)[psombra_d1], np.array(yzr)[psombra_d1], 'rd', markersize=10)
                # Pintamos cesión
                plt.plot(pestcesx, pestcesy, '*')
                plt.quiver(pestcesx, pestcesy, np.cos(np.deg2rad(270 - Dest)), np.sin(np.deg2rad(270 - Dest)), autoscale=200)
                plt.plot(xpolyces, ypolyces)
                plt.plot(np.array(xzr)[pcesion_d1_p1], np.array(yzr)[pcesion_d1_p1], 'gd', markersize=10)
                # Pintamos para la figura
                # patch(xpolyces,ypolyces,[54 237 17]./255,'Facealpha',0.5,'EdgeColor',[1 1 1],'LineWidth',2)
                # patch(xpolydif,ypolydif,[17 107 237]./255,'Facealpha',0.5,'EdgeColor',[1 1 1],'LineWidth',2)