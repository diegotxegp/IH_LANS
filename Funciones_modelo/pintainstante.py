import numpy as np
import matplotlib.pyplot as plt
import time
from datetime import datetime, timedelta

from Funciones_modelo.pinta_perfiles import pinta_perfiles

def pintainstante(PERF, YLTi, ACT, EA, Hi, D0, Di, wi, t, it, escalaprin):

    tref = datetime(1,1,1)
    delta = timedelta(t[it])

    instante = tref + delta

    # Pintamos instante de lc
    plt.title(str(instante.date()))
    
    xon = [perf['xon'] for perf in PERF]
    nx = [perf['nx'] for perf in PERF]
    yon = [perf['yon'] for perf in PERF]
    ny = [perf['ny'] for perf in PERF]

    xlci = xon + nx * YLTi
    ylci = yon + ny * YLTi
    
    # Pintamos LC y mar con área
    politierrax = xon + list(reversed(xlci)) + xon[0]
    politierray = yon + list(reversed(ylci)) + yon[0]
    
    #tierra = plt.Polygon([(1,2),(3,4)], facecolor=(247 / 255, 243 / 255, 141 / 255), edgecolor='none')
    #plt.gca().add_patch(tierra)
    
    minx = 1e9
    maxx = -1e9
    miny = 1e9
    maxy = -1e9

    pinta_perfiles(PERF)
    
    for i in range(len(PERF)):
        xon = PERF[i]['xon']
        xof = PERF[i]['xof']
        yon = PERF[i]['yon']
        yof = PERF[i]['yof']
        
        plt.plot([xon, xof], [yon, yof], 'k.:', linewidth=0.5)
        
        if xon < minx or xof < minx:
            minx = min(xon, xof)
        if xon > maxx or xof > maxx:
            maxx = max(xon, xof)
        if yon < miny or yof < miny:
            miny = min(yon, yof)
        if yon > maxy or yof > maxy:
            maxy = max(yon, yof)
    
    # Pintamos oleaje en cabeza de perfiles
    quiver_scale = escalaprin
    quiver_x = [perf['xof'] for perf in PERF]
    quiver_y = [perf['yof'] for perf in PERF]
    quiver_dir = 270 - np.array(D0)
    quiver_dx = quiver_scale * np.cos(np.deg2rad(quiver_dir))
    quiver_dy = quiver_scale * np.sin(np.deg2rad(quiver_dir))
    
    plt.quiver(quiver_x, quiver_y, quiver_dx, quiver_dy, angles='xy', scale_units='xy', scale=1, color='k', linewidth=2)
    
    # Pintamos ZR0
    xon = [perf['xon'] for perf in PERF]
    nx = [perf['nx'] for perf in PERF]
    yon = [perf['yon'] for perf in PERF]
    ny = [perf['ny'] for perf in PERF]

    xr = np.array(xon) + np.array(nx) * (np.array(wi) + YLTi)
    yr = np.array(yon) + np.array(ny) * (np.array(wi) + YLTi)
    xr2 = np.array(xon) + np.array(nx) * (1.25 * np.array(wi) + YLTi)
    yr2 = np.array(yon) + np.array(ny) * (1.25 * np.array(wi) + YLTi)
    
    plt.plot(xr, yr, 'b--', linewidth=2)
    plt.plot(xr, yr, 'bo', markerfacecolor='b')
    plt.plot(xr2, yr2, 'r--', linewidth=2)
    
    # Pintamos dirección sin corregir
    quiver_x = xr
    quiver_y = yr
    quiver_dir = 270 - np.array(Di)
    quiver_dx = quiver_scale * np.cos(np.deg2rad(quiver_dir))
    quiver_dy = quiver_scale * np.sin(np.deg2rad(quiver_dir))
    
    plt.quiver(quiver_x, quiver_y, quiver_dx, quiver_dy, angles='xy', scale_units='xy', scale=1, color='b', linewidth=2)
    
    # Calculamos los vectores normales
    nxp = np.gradient(xlci) / np.hypot(np.gradient(ylci), np.gradient(xlci))
    nyp = np.gradient(ylci) / np.hypot(np.gradient(ylci), np.gradient(xlci))
    
    # Pintamos LC0
    xlc0 = [perf['xon'] + perf['nx'] * perf['yc'] for perf in PERF]
    ylc0 = [perf['yon'] + perf['ny'] * perf['yc'] for perf in PERF]
    
    plt.plot(xlc0, ylc0, 'k--', linewidth=2)
    
    # Pintamos estructuras activas
    if EA:
        tiposact = list(EA.keys())
        for acttipo in tiposact:
            if acttipo == 'noest':
                continue
            
            acttipo_estructuras = EA[acttipo]
            for i in range(len(acttipo_estructuras)):
                index = acttipo_estructuras[i]
                ACT0 = ACT[index]
                
                if acttipo == 'dper':
                    plt.plot(ACT0['X'], ACT0['Y'], linewidth=5, color=(126/255, 126/255, 129/255))
                    
                    Dest = ACT0['Dir'][it]
                    quiver_x = ACT0['X'][-1]
                    quiver_y = ACT0['Y'][-1]
                    quiver_dir = 270 - Dest
                    quiver_dx = quiver_scale * np.cos(np.deg2rad(quiver_dir))
                    quiver_dy = quiver_scale * np.sin(np.deg2rad(quiver_dir))
                    
                    plt.quiver(quiver_x, quiver_y, quiver_dx, quiver_dy, angles='xy', scale_units='xy', scale=1, color='b', linewidth=2)
                elif acttipo == 'dpar':
                    plt.plot(ACT0['X'], ACT0['Y'], linewidth=5, color=(126/255, 126/255, 129/255))
                    
                    Dest = ACT0['Dir'][it]
                    quiver_x = ACT0['X']
                    quiver_y = ACT0['Y']
                    quiver_dir = 270 - Dest
                    quiver_dx = quiver_scale * np.cos(np.deg2rad(quiver_dir))
                    quiver_dy = quiver_scale * np.sin(np.deg2rad(quiver_dir))
                    
                    plt.quiver(quiver_x, quiver_y, quiver_dx, quiver_dy, angles='xy', scale_units='xy', scale=1, color='b', linewidth=2)
                elif acttipo == 'dperpar':
                    plt.plot(ACT0['X'], ACT0['Y'], linewidth=5, color=(126/255, 126/255, 129/255))
                    
                    Dest = ACT0['Dir'][it]
                    
                    if len(ACT0['X']) == 4:
                        quiver_x = [ACT0['X'][-2], ACT0['X'][-1]]
                        quiver_y = [ACT0['Y'][-2], ACT0['Y'][-1]]
                        quiver_dir = 270 - Dest
                        quiver_dx = quiver_scale * np.cos(np.deg2rad(quiver_dir))
                        quiver_dy = quiver_scale * np.sin(np.deg2rad(quiver_dir))
                        
                        plt.quiver(quiver_x, quiver_y, quiver_dx, quiver_dy, angles='xy', scale_units='xy', scale=1, color='b', linewidth=2)
                    elif len(ACT0['X']) == 3:
                        nxest0 = ACT0['X'][0] - ACT0['X'][1]
                        nyest0 = ACT0['Y'][0] - ACT0['Y'][1]
                        mod = np.hypot(nxest0, nyest0)
                        nxest = nxest0 / mod
                        nyest = nyest0 / mod
                        angper = np.degrees(np.angle(nxest + nyest * 1j))
                        
                        nxest1 = ACT0['X'][-1] - ACT0['X'][-2]
                        nyest1 = ACT0['Y'][-1] - ACT0['Y'][-2]
                        mod1 = np.hypot(nxest1, nyest1)
                        nxest2 = nxest1 / mod1
                        nyest2 = nyest1 / mod1
                        angpar = np.degrees(np.angle(nxest2 + nyest2 * 1j))
                        
                        deltaest = (angper - angpar) % 360
                        deltaola = (angper - Dest) % 360
                        
                        dref = (angper - Dest) % 360
                        
                        if np.sign(deltaola) != np.sign(deltaest):
                            quiver_x = [ACT0['X'][-2], ACT0['X'][-1]]
                            quiver_y = [ACT0['Y'][-2], ACT0['Y'][-1]]
                            quiver_dir = 270 - Dest
                            quiver_dx = quiver_scale * np.cos(np.deg2rad(quiver_dir))
                            quiver_dy = quiver_scale * np.sin(np.deg2rad(quiver_dir))
                            
                            plt.quiver(quiver_x, quiver_y, quiver_dx, quiver_dy, angles='xy', scale_units='xy', scale=1, color='b', linewidth=2)
                        else:
                            quiver_x = ACT0['X'][-1]
                            quiver_y = ACT0['Y'][-1]
                            quiver_dir = 270 - Dest[-1]
                            quiver_dx = quiver_scale * np.cos(np.deg2rad(quiver_dir))
                            quiver_dy = quiver_scale * np.sin(np.deg2rad(quiver_dir))
                            
                            plt.quiver(quiver_x, quiver_y, quiver_dx, quiver_dy, angles='xy', scale_units='xy', scale=1, color='b', linewidth=2)
                elif acttipo == 'drig':
                    # Resaltamos la LC en la zona rigidizada
                    plt.plot(ACT0['X'], ACT0['Y'], linewidth=5, color=(126/255, 126/255, 129/255))
                elif acttipo == 'nour':
                    xnour0 = ACT0['X']
                    xnour = np.append(xnour0, xnour0[0])
                    ynour0 = ACT0['Y']
                    ynour = np.append(ynour0, ynour0[0])
                    
                    plt.plot(xnour, ynour, linewidth=2, color=(232/255, 177/255, 19/255))
    
    plt.xlabel('Longitudinal, [m]')
    plt.ylabel('Transversal, [m]')
    plt.gca().set_aspect('equal', adjustable='box')
    
    # Ajusta los límites de acuerdo a los valores mínimos y máximos calculados
    plt.xlim(minx-50, maxx+50)
    plt.ylim(miny-50, maxy+500)
    
    mng = plt.get_current_fig_manager()
    mng.window.state('zoomed')
    plt.show()

    time.sleep(1)
    plt.close('all')
