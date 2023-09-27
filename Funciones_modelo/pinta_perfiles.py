import matplotlib.pyplot as plt
import numpy as np

from Funciones_modelo.detecta_lcs_bcpun import detecta_lcs_bcpun

def pinta_perfiles(PERF):
    PBC, PNL = detecta_lcs_bcpun(PERF)
    
    for i in range(len(PERF)):
        if 'Tipo' in PERF[i]:
            if PERF[i]['Tipo'] == 'lcs':
                color = 'g'
            elif PERF[i]['Tipo'] == 'cs':
                color = 'y'
            elif PERF[i]['Tipo'] == 's':
                color = 'r'
            elif PERF[i]['Tipo'] == 'ac':
                color = 'k'
        else:
            color = 'k'
        
        if np.isin([i],PBC):
            linestyle = '--'
            linewidth = 2
        else:
            linestyle = ':'
            linewidth = 1
        
        plt.plot([PERF[i]['xon'], PERF[i]['xof']], [PERF[i]['yon'], PERF[i]['yof']], linestyle, color=color, linewidth=linewidth)
        xm = (PERF[i]['xon'] + PERF[i]['xof']) / 2
        ym = (PERF[i]['yon'] + PERF[i]['yof']) / 2
        # plt.text(xm, ym, str(i))  # Descomentar para etiquetar los perfiles
    plt.gca().set_aspect('equal', adjustable='box')
    