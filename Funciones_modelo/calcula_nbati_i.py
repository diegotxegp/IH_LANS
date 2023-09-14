import numpy as np

from Generales.anglecalc import anglecalc

def calcula_nbati_i(xprop, yprop, *varargin):
    # Calculamos normales salientes de la línea
    # Gradient devuelve la diferencia central entre el siguiente y el anterior
    # En las esquinas, la diferencia es de un solo lado
    nx = np.gradient(xprop)
    ny = np.gradient(yprop)

    # Angulo normal
    angnorm = anglecalc(nx, -ny)
    nbati = 90 - angnorm

    # Corregimos saltos en tramos LCS y estructuras
    posBC = [i for i, x in enumerate(varargin) if x == 'PBC'] + 1
    posACT = [i for i, x in enumerate(varargin) if x == 'ACT'] + 1

    if posBC:
        PBC = varargin[posBC[0]]
        ntramos = PBC.shape[0]

        for i in range(ntramos):
            if PBC[i, 0] != 1 and PBC[i, 0] > 2:
                nbati[PBC[i, 0] - 1] = nbati[PBC[i, 0] - 2]
                nbati[PBC[i, 0]] = nbati[PBC[i, 0] + 1]

            if PBC[i, 1] < len(xprop) - 1:
                nbati[PBC[i, 1]] = nbati[PBC[i, 1] - 1]
                nbati[PBC[i, 1] + 1] = nbati[PBC[i, 1] + 2]

    if posACT:
        ACT = varargin[posACT[0]]
        # Completar en los diques si es necesario
        # Agregar código aquí si es necesario

    return nbati