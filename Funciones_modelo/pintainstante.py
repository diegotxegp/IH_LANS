import numpy as np
import matplotlib.pyplot as plt

def pintainstante(PERF, YLTi, ACT, EA, Hi, D0, Di, wi, t, it, escalaprin):
    # Crea una figura y ajusta su posición
    # h = plt.figure()
    # plt.title(t[it])
    xlci = [perf["xon"] + perf["nx"] * YLTi for perf in PERF]
    ylci = [perf["yon"] + perf["ny"] * YLTi for perf in PERF]

    # Pinta la LC y el área del mar
    politierrax = [perf["xon"] for perf in PERF] + list(reversed(xlci)) + [PERF[0]["xon"]]
    politierray = [perf["yon"] for perf in PERF] + list(reversed(ylci)) + [PERF[0]["yon"]]
    tierra = plt.Polygon(np.array([politierrax, politierray]).T, facecolor=(247/255, 243/255, 141/255))
    plt.gca().add_patch(tierra)

    minx = float("inf")
    maxx = float("-inf")
    miny = float("inf")
    maxy = float("-inf")

    for perf in PERF:
        xon = perf["xon"]
        xof = perf["xof"]
        yon = perf["yon"]
        yof = perf["yof"]
        plt.plot([xon, xof], [yon, yof], 'k.:', linewidth=2)
        if xon < minx or xof < minx:
            minx = min(xon, xof)
        if xon > maxx or xof > maxx:
            maxx = max(xon, xof)
        if yon < miny or yof < miny:
            miny = min(yon, yof)
        if yon > maxy or yof > maxy:
            maxy = max(yon, yof)

    # Pinta el oleaje en la cabeza de los perfiles
    quiver_params = {
        'angles': 270 - np.array(D0),
        'scale': escalaprin,
        'scale_units': 'xy',
        'pivot': 'tail',
        'color': 'k',
        'linewidth': 2
    }
    for perf, q_params in zip(PERF, quiver_params):
        plt.quiver(perf.xof, perf.yof, *q_params)

    # Pinta ZR0
    xr = [perf["xon"] + perf["nx"] * (wi + YLTi) for perf in PERF]
    yr = [perf["yon"] + perf["ny"] * (wi + YLTi) for perf in PERF]
    xr2 = [perf["xon"] + perf["nx"] * (1.25 * wi + YLTi) for perf in PERF]
    yr2 = [perf["yon"] + perf["ny"] * (1.25 * wi + YLTi) for perf in PERF]
    plt.plot(xr, yr, 'b--', linewidth=2)
    plt.plot(xr, yr, 'bo', markerfacecolor='b')
    plt.plot(xr2, yr2, 'r--', linewidth=2)

    # Pinta dirección sin corregir
    quiver_params = {
        'angles': 270 - np.array(Di),
        'scale': escalaprin,
        'scale_units': 'xy',
        'pivot': 'tail',
        'color': 'b',
        'linewidth': 2
    }
    for (xi, yi), q_params in zip(zip(xr, yr), quiver_params):
        plt.quiver(xi, yi, *q_params)

    nxp = [np.gradient(xlci) / np.hypot(np.gradient(ylci), np.gradient(xlci))]
    nyp = [np.gradient(ylci) / np.hypot(np.gradient(ylci), np.gradient(xlci))]

    # Pinta LC0
    xlc0 = [perf["xon"] + perf["nx"] * perf["yc"] for perf in PERF]
    ylc0 = [perf["yon"] + perf["ny"] * perf["yc"] for perf in PERF]
    plt.plot(xlc0, ylc0, 'k--', linewidth=2)

    # Pinta estructuras activas
    if EA is not None:
        tiposact = EA.keys()
        for acttipo in tiposact:
            if acttipo == 'noest':
                continue

            nest = len(EA[acttipo])

            for ii in range(nest):
                index = EA[acttipo][ii]
                ACT0 = ACT[index]

                if acttipo == 'dper':
                    plt.plot(ACT0.X, ACT0.Y, linewidth=5, color=(126/255, 126/255, 129/255))
                    Dest = ACT0.Dir[it]
                    plt.quiver(ACT0.X[-1], ACT0.Y[-1], np.cos(np.deg2rad(270 - Dest)) * escalaprin, np.sin(np.deg2rad(270 - Dest)) * escalaprin, 0, 'b', linewidth=2)
                elif acttipo == 'dpar':
                    plt.plot(ACT0.X, ACT0.Y, linewidth=5, color=(126/255, 126/255, 129/255))
                    Dest = ACT0.Dir[it]
                    plt.quiver(ACT0.X, ACT0.Y, np.cos(np.deg2rad(270 - Dest)) * escalaprin, np.sin(np.deg2rad(270 - Dest)) * escalaprin, 0, 'b', linewidth=2)
                elif acttipo == 'dperpar':
                    plt.plot(ACT0.X, ACT0.Y, linewidth=5, color=(126/255, 126/255, 129/255))
                    Dest = ACT0.Dir[it]
                    if len(ACT0.X) == 4:
                        plt.quiver([ACT0.X[-2], ACT0.X[-1]], [ACT0.Y[-2], ACT0.Y[-1]], np.cos(np.deg2rad(270 - Dest)) * escalaprin, np.sin(np.deg2rad(270 - Dest)) * escalaprin, 0, 'b', linewidth=2)
                    elif len(ACT0.X) == 3:
                        nxest0 = ACT0.X[0] - ACT0.X[1]
                        nyest0 = ACT0.Y[0] - ACT0.Y[1]
                        mod = np.hypot(nxest0, nyest0)
                        nxest = nxest0 / mod
                        nyest = nyest0 / mod
                        angper = np.degrees(np.arctan2(nyest, nxest))
                        nxest1 = ACT0.X[-1] - ACT0.X[-2]
                        nyest1 = ACT0.Y[-1] - ACT0.Y[-2]
                        mod1 = np.hypot(nxest1, nyest1)
                        nxest2 = nxest1 / mod1
                        nyest2 = nyest1 / mod1
                        angpar = np.degrees(np.arctan2(nyest2, nxest2))
                        deltaest = np.rad2deg(np.arccos(np.cos(np.deg2rad(angper)) * np.cos(np.deg2rad(angpar)) + np.sin(np.deg2rad(angper)) * np.sin(np.deg2rad(angpar))))
                        Destcar = 270 - Dest[0]
                        deltaola = np.rad2deg(np.arccos(np.cos(np.deg2rad(angper)) * np.cos(np.deg2rad(Destcar)) + np.sin(np.deg2rad(angper)) * np.sin(np.deg2rad(Destcar))))
                        dref = np.rad2deg(np.arccos(np.cos(np.deg2rad(angper)) * np.cos(np.deg2rad(Destcar)) + np.sin(np.deg2rad(angper)) * np.sin(np.deg2rad(Destcar))))

                        if np.sign(deltaola) != np.sign(deltaest):
                            plt.quiver([ACT0.X[-2], ACT0.X[-1]], [ACT0.Y[-2], ACT0.Y[-1]], np.cos(np.deg2rad(270 - Dest)) * escalaprin, np.sin(np.deg2rad(270 - Dest)) * escalaprin, 0, 'b', linewidth=2)
                        else:
                            plt.quiver(ACT0.X[-1], ACT0.Y[-1], np.cos(np.deg2rad(270 - Dest[-1])) * escalaprin, np.sin(np.deg2rad(270 - Dest[-1])) * escalaprin, 0, 'b', linewidth=2)
                elif acttipo == 'drig':
                    plt.plot(ACT0.X, ACT0.Y, linewidth=5, color=(126/255, 126/255, 129/255))
                elif acttipo == 'nour':
                    xnour0 = ACT0.X
                    xnour = xnour0 + [xnour0[0]]
                    ynour0 = ACT0.Y
                    ynour = ynour0 + [ynour0[0]]
                    plt.plot(xnour, ynour, linewidth=2, color=(232/255, 177/255, 19/255))

    plt.axis('equal')
    plt.xlim(minx, maxx)
    plt.ylim(miny, maxy)
    plt.show()