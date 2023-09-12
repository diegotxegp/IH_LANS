import numpy as np

def detecta_lcs_bcpun(PERF):
    tipo = [p['Tipo'] for p in PERF]
    lcs = [i for i, t in enumerate(tipo) if 'lcs' in t]
    PNL = [i for i, t in enumerate(tipo) if 'lcs' not in t]
    dlcs = np.diff(lcs)
    ntramos = np.sum(dlcs != 1) + 1
    pcambios = np.where(dlcs != 1)[0]

    PBC = np.empty((ntramos, 2))
    
    if pcambios.size > 0:
        for i in range(ntramos):
            if i == 0:
                PBC[i, :] = [lcs[0], lcs[pcambios[i]]]
            elif i == ntramos - 1:
                PBC[i, :] = [lcs[pcambios[i - 1] + 1], lcs[-1]]
            else:
                PBC[i, :] = [lcs[pcambios[i - 1] + 1], lcs[pcambios[i]]]
    else:
        PBC = np.array([lcs[0], lcs[-1]])

    return PBC, PNL