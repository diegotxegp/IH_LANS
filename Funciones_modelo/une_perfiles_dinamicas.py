import numpy as np

def une_perfiles_dinamicas(PERF, DYN):
    # Une los perfiles a sus dinámicas forzadoras por cercanía con el punto
    # offshore del perfil
    dinperf = np.zeros(len(PERF), dtype=int)
    xdyn = [dyn['X'] for dyn in DYN]
    ydyn = [dyn['Y'] for dyn in DYN]
    
    for i in range(len(PERF)):
        xof = PERF[i]['xof']
        yof = PERF[i]['yof']
        xon = PERF[i]['xon']
        yon = PERF[i]['yon']
        xm = (xof + xon) / 2
        ym = (yof + yon) / 2
        d = np.hypot(xm - xdyn, ym - ydyn)
        posmin = np.argmin(d)
        dinperf[i] = posmin
    
    return dinperf