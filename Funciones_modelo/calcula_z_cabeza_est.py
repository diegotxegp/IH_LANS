import numpy as np

from Funciones_modelo.perfiles_difraccion import perfiles_difraccion
from Funciones_modelo.point_to_line_distance import point_to_line_distance
from Funciones_modelo.calcula_z import calcula_z

def calcula_z_cabeza_est(ACT, PERF, tipo):
    # Calculamos la z en la cabeza de la estructura y la normal a considerar en la propagación
    if not ACT:
        return ACT
    
    iper = [act['Tipo'] == 'dper' for act in ACT]
    ipar = [act['Tipo'] == 'dpar' for act in ACT]
    iparper = [act['Tipo'] == 'dperpar' for act in ACT]
    idif = [i for i, (p, a, ap) in enumerate(zip(iper, ipar, iparper)) if p or a or ap]
    
    if tipo == 'FIX':
        return ACT
    
    for i in idif:
        est = i
        ACT0 = ACT[est]
        tipo = ACT0['Tipo']
        
        if tipo == 'dper':
            pafc = perfiles_difraccion(ACT0, PERF)
            z = np.zeros_like(pafc, dtype=float)
            d1 = np.zeros_like(pafc, dtype=float)
            nbatiP = np.zeros_like(pafc, dtype=float)
            
            for ifila in range(pafc.shape[0]):
                for icolumn in range(pafc.shape[1]):
                    xonp = PERF[pafc[ifila, icolumn]]['xon']
                    yonp = PERF[pafc[ifila, icolumn]]['yon']
                    xofp = PERF[pafc[ifila, icolumn]]['xof']
                    yofp = PERF[pafc[ifila, icolumn]]['yof']
                    xest = ACT0['X'][-1]
                    yest = ACT0['Y'][-1]
                    d = point_to_line_distance([xest, yest], [xonp, yonp], [xofp, yofp])
                    d1[ifila, icolumn] = d
                    dondif = np.hypot(xonp - xest, yonp - yest)
                    dperfil = np.sqrt(dondif**2 - d**2)
                    z[ifila, icolumn] = calcula_z(dperfil, PERF[pafc[ifila, icolumn]], tipo)
                    nbatiP[ifila, icolumn] = PERF[pafc[ifila, icolumn]]['nbati']
                    
            ACT[est]['Z'] = np.sum(z * d1) / np.sum(d1)
            ACT[est]['nbati'] = np.sum(nbatiP * d1) / np.sum(d1)
        
        elif tipo == 'dperpar':
            pafc = perfiles_difraccion(ACT0, PERF)
            z = np.zeros_like(pafc, dtype=float)
            d1 = np.zeros_like(pafc, dtype=float)
            nbatiP = np.zeros_like(pafc, dtype=float)
            
            for ifila in range(pafc.shape[0]):
                for icolumn in range(pafc.shape[1]):
                    ind = icolumn - 2
                    xonp = PERF[pafc[ifila, icolumn]]['xon']
                    yonp = PERF[pafc[ifila, icolumn]]['yon']
                    xofp = PERF[pafc[ifila, icolumn]]['xof']
                    yofp = PERF[pafc[ifila, icolumn]]['yof']
                    xest = ACT0['X'][len(ACT0['X']) + ind]
                    yest = ACT0['Y'][len(ACT0['Y']) + ind]
                    d = point_to_line_distance([xest, yest], [xonp, yonp], [xofp, yofp])
                    d1[ifila, icolumn] = d
                    dondif = np.hypot(xonp - xest, yonp - yest)
                    dperfil = np.sqrt(dondif**2 - d**2)
                    z[ifila, icolumn] = calcula_z(dperfil, PERF[pafc[ifila, icolumn]], tipo)
                    nbatiP[ifila, icolumn] = PERF[pafc[ifila, icolumn]]['nbati']
                    
            ACT[est]['Z'] = np.sum(z * d1) / np.sum(d1)
            ACT[est]['nbati'] = np.sum(nbatiP * d1) / np.sum(d1)
        
        elif tipo == 'dpar':
            pafc = perfiles_difraccion(ACT[est], PERF)
            z = np.zeros_like(pafc, dtype=float)
            d1 = np.zeros_like(pafc, dtype=float)
            nbatiP = np.zeros_like(pafc, dtype=float)
            
            for ifila in range(pafc.shape[0]):
                for icolumn in range(pafc.shape[1]):
                    xonp = PERF[pafc[ifila, icolumn]]['xon']
                    yonp = PERF[pafc[ifila, icolumn]]['yon']
                    xofp = PERF[pafc[ifila, icolumn]]['xof']
                    yofp = PERF[pafc[ifila, icolumn]]['yof']
                    xest = ACT0['X'][ifila]
                    yest = ACT0['Y'][ifila]
                    d = point_to_line_distance([xest, yest], [xonp, yonp], [xofp, yofp])
                    d1[ifila, icolumn] = d
                    dondif = np.hypot(xonp - xest, yonp - yest)
                    dperfil = np.sqrt(dondif**2 - d**2)
                    z[ifila, icolumn] = calcula_z(dperfil, PERF[pafc[ifila, icolumn]], tipo)
                    nbatiP[ifila, icolumn] = PERF[pafc[ifila, icolumn]]['nbati']
                    
            ACT[est]['Z'] = np.sum(z * d1) / np.sum(d1)
            ACT[est]['nbati'] = np.sum(nbatiP * d1) / np.sum(d1)
    
    return ACT