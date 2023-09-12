import numpy as np

def perfiles_actuacion(ACT, PERF):
    # Identificamos los perfiles de afectación de las estructuras
    # Asumimos que los perfiles son normales a la LC, al igual que los diques
    xon = [perf['xon'] for perf in PERF]
    yon = [perf['yon'] for perf in PERF]
    xof = [perf['xof'] for perf in PERF]
    yof = [perf['yof'] for perf in PERF]
    xlc = xon + np.array([perf['nx'] for perf in PERF]) * np.array([perf['yc'] for perf in PERF])
    ylc = yon + np.array([perf['ny'] for perf in PERF]) * np.array([perf['yc'] for perf in PERF])
    
    for i in range(len(ACT)):
        ACT0 = ACT[i]
        if ACT0['Tipo'] == 'dper' or ACT0['Tipo'] == 'dperpar':
            # Como punto de referencia de la estructura, cogemos la intersección
            # del dique con la LC inicial
            intersection = np.poly1d([ACT0['X'][0], ACT0['Y'][0]]).roots(np.poly1d([xlc, ylc]))
            xest, yest = intersection[0], intersection[1]
            d1 = np.hypot(np.array(xon) - xest, np.array(yon) - yest)
            posmin1 = np.argmin(d1)
            
            # Buscamos en las celdas contiguas
            if np.point_inside_polygon(xest, yest, [xon[posmin1], xon[posmin1+1], xof[posmin1+1], xof[posmin1], xon[posmin1]],
                                        [yon[posmin1], yon[posmin1+1], yof[posmin1+1], yof[posmin1], yon[posmin1]]):
                pafc = [posmin1, posmin1+1]
            elif np.point_inside_polygon(xest, yest, [xon[posmin1-1], xon[posmin1], xof[posmin1], xof[posmin1-1], xon[posmin1-1]],
                                          [yon[posmin1-1], yon[posmin1], yof[posmin1], yof[posmin1-1], yon[posmin1-1]]):
                pafc = [posmin1-1, posmin1]
            else:
                # Recorremos todas las celdas y vemos dónde cae el punto inicial
                count = 0
                noin = 0
                while noin == 0:
                    count += 1
                    noin = np.point_inside_polygon(xest, yest, [xon[count], xon[count+1], xof[count+1], xof[count], xon[count]],
                                                   [yon[count], yon[count+1], yof[count+1], yof[count], yon[count]])
                    pafc = [count, count+1]
            ACT[i]['PERF'] = pafc
        elif ACT0['Tipo'] == 'drig' or ACT0['Tipo'] == 'dpar':
            # Cortamos la estructura con los perfiles
            pafc = []
            yafc = []
            xest = ACT0['X']
            yest = ACT0['Y']
            count = 0
            for ip in range(len(PERF)):
                xperf = [PERF[ip]['xon'], PERF[ip]['xof']]
                yperf = [PERF[ip]['yon'], PERF[ip]['yof']]
                # Cortamos los perfiles con las estructuras
                intersection = np.poly1d(xperf).roots(np.poly1d([xest, yest]))
                if len(intersection) > 0:
                    count += 1
                    pafc.append(ip)
                    # Distancia con el punto onshore del perfil
                    yafc.append(np.hypot(xperf[0] - intersection[0], yperf[0] - intersection[1]))
            ACT[i]['PERF'] = pafc
            ACT[i]['YAFC'] = yafc
    
    return ACT