import numpy as np

from Generales.wmoore_calc import wmoore_calc

def calcula_z(dperfil, PERF0, tipo):
    # Función para calcular el calado de un punto a d del origen del perfil

    if tipo == 'DEAN_D50':
        xdean = dperfil - PERF0['Yberma']
        Adean = 0.51 * wmoore_calc([PERF0['d50']]) ** 0.44
        z = PERF0['Berma'] - Adean * xdean ** (2/3)
    elif tipo == 'DEAN_A':
        xdean = dperfil - PERF0['Yberma']
        z = PERF0['Berma'] - PERF0['Adean'] * xdean ** (2/3)
    elif tipo == 'PROFILE':
        x1 = PERF0['xall']
        y1 = PERF0['yall']
        dall = [0] + np.cumsum(np.hypot(x1[1:] - x1[:-1], y1[1:] - y1[:-1])).tolist()
        zall = PERF0['zall']
        _, z = polyxpoly(dall, zall, [dperfil, dperfil], [min(zall), max(zall)])
    else:
        raise ValueError('Not a valid option to calculate the depth of the difraction point')

    return z