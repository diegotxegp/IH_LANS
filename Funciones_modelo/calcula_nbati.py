import numpy as np
from scipy.spatial import distance
from shapely.geometry import LineString, Point

def calcula_nbati(xprop, yprop, PERF):
    # Calculamos las normales salientes de la línea
    dx = np.gradient(xprop)
    dy = np.gradient(yprop)
    
    # Angulo normal
    angnorm = np.angle(dx - 1j * dy, deg=True)
    
    # Calculamos perfiles
    for i in range(len(PERF)):
        xp = [PERF[i]['xon'], PERF[i]['xof']]
        yp = [PERF[i]['yon'], PERF[i]['yof']]
        
        # Intersectamos el contorno offshore
        line1 = LineString([(xp[0], yp[0]), (xp[1], yp[1])])
        line2 = LineString(list(zip(xprop, yprop)))
        intersection = line1.intersection(line2)
        
        if intersection.is_empty:
            continue
        
        if intersection.geom_type == 'Point':
            xint, yint = intersection.x, intersection.y
        else:
            # Si hay múltiples intersecciones, elegimos la más cercana
            closest_point = min(intersection, key=lambda p: distance.euclidean((xprop[0], yprop[0]), (p.x, p.y)))
            xint, yint = closest_point.x, closest_point.y
        
        # Calculamos la distancia entre el punto de intersección y los puntos del perfil
        dist = [distance.euclidean((xint, yint), (xprop[j], yprop[j])) for j in range(len(xprop))]
        
        # Ordenamos las distancias en orden creciente
        sorted_indices = np.argsort(dist)
        
        # Promediamos el ángulo entre puntos adyacentes ponderado por la distancia
        idx1, idx2 = sorted_indices[:2]
        nbati = 90 - (angnorm[idx1] * dist[idx1] + angnorm[idx2] * dist[idx2]) / (dist[idx1] + dist[idx2])
        
        PERF[i]['nbati'] = nbati
        
    return PERF