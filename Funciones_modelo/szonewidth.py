import numpy as np

from Generales.wmoore_calc import wmoore_calc

def szonewidth(hb, type, PERF0):
    # Anchura de la ZR desde la línea de costa
    if type == 'DEAN_D50':
        Adean = 0.51 * wmoore_calc([PERF0['d50']]) ** 0.44
        wb = ((hb) / Adean) ** (3/2)
    elif type == 'DEAN_A':
        wb = ((hb) / PERF0['Adean']) ** (3/2)
    elif type == 'PROFILE':
        xall = np.array(PERF0['xall'])
        yall = np.array(PERF0['yall'])
        zall = np.array(PERF0['zall'])
        dall = [0] + list(np.cumsum(np.sqrt((xall[1:] - xall[:-1]) ** 2 + (yall[1:] - yall[:-1]) ** 2)))
        
        # Cogemos la parte sumergida del perfil
        posem = np.where(zall > 0)[0]
        da = np.array(dall)
        da[posem] = []
        za = np.array(zall)
        za[posem] = []
        
        try:
            # Si funciona
            wb = np.interp(-hb, za, da) - PERF0['yc']
        except:
            wb = np.zeros_like(hb)
            for i in range(len(hb)):
                xi, _ = polyxpoly(dall, zall, [min(dall), max(dall)], [-hb[i], -hb[i]])
                if len(xi) > 1:
                    # Escogemos primero o último en función del lado de mar
                    if zall[0] > zall[-1]:  # Lado mar a la derecha
                        wb[i] = xi[-1] - PERF0['yc']
                    else:
                        wb[i] = xi[0] - PERF0['yc']
                else:
                    wb[i] = xi[0] - PERF0['yc']
    else:
        raise ValueError('Opción no válida para calcular la anchura de la zona de rompientes')

    return wb