import numpy as np

from Funciones_modelo.une_perfiles_dinamicas import une_perfiles_dinamicas
from Funciones_modelo.snell_shoalrefrot import snell_shoalrefrot
from Funciones_modelo.szonewidth import szonewidth
from Generales.anglecalc import anglecalc

def propaga_rotura(PERF, DYN, gamma, t, refNMM, cotasZ, calcularotur):
    # Unimos perfiles y dinámicas
    dinperf = une_perfiles_dinamicas(PERF, DYN)
    
    # Buscamos fechas coincidentes entre el tiempo de análisis y las dinámicas
    poscalc = np.where(np.isin([d["t"] for d in DYN], t).tolist()[0])[0].tolist()

    sum = []
    for i in range(len(poscalc)-1):
        sum.append(poscalc[i+1]-poscalc[i])
    
    if np.sum(sum) != len(t) - 1:
        raise ValueError('No hay dinámicas en las fechas del análisis')
    
    DYNP = []
    
    if calcularotur == 1:
        for i in range(len(PERF)):
            print(f'Perfil {i+1} de {len(PERF)}')
            id = dinperf[i]  # Unión perfil dinámicas
            AT = DYN[id]["AT"][poscalc]
            SS = DYN[id]["SS"][poscalc]
            SLR = DYN[id]["SLR"][poscalc]
            nivel = AT + SS + SLR + refNMM + DYN[id]["h0"]
            Hs = np.array(DYN[id]["Hs"])[poscalc]
            Dir = DYN[id]["Dir"][poscalc]
            Tp = DYN[id]["Tp"][poscalc]
            nbati = PERF[i]["nbati"]
            Hb, _, _, alphab, hb, _, pnoprop = snell_shoalrefrot(Hs, Tp, Dir, nivel, gamma, nbati)
            # Cocinamos los no propagados
            Hb[pnoprop] = 0.1  # No demasiado pequeña por posibles problemas numéricos
            alphab[pnoprop] = anglecalc(np.sin(np.deg2rad(nbati)), np.cos(np.deg2rad(nbati)))  # Incidencia normal en naúticas de donde vienen
            DYNP.append({
                'Hb': Hb,
                'Dirb': alphab,
                'hb': hb,  # La profundidad de las no propagadas mantenemos la inicial para que no difracten
                'Dir0': Dir,
                'SS': SS,
                'AT': AT,
                'SLR': SLR,
                'wb': szonewidth(hb - SS - AT, cotasZ, PERF[i])  # Trabajamos en nmm local
            })
    else:
        for i in range(len(PERF)):
            print(f'Perfil {i+1} de {len(PERF)}')
            id = dinperf[i]  # Unión perfil dinámicas
            AT = DYN[id]["AT"][:len(poscalc)]
            SS = DYN[id]["SS"][:len(poscalc)]
            SLR = DYN[id]["SLR"][:len(poscalc)]
            nivel = AT + SS + SLR + refNMM + DYN[id]["h0"]
            Hs = DYN[id]["Hs"][:len(poscalc)]
            Dir = DYN[id]["Dir"][:len(poscalc)]
            Tp = DYN[id]["Tp"][:len(poscalc)]
            nbati = PERF[i]["nbati"]
            hb = [hs / gamma for hs in Hs]
            DYNP.append({
                'Hb': Hs,
                'Dirb': Dir,
                'hb': hb,
                'Dir0': Dir,
                'SS': SS,
                'AT': AT,
                'SLR': SLR,
                'wb': szonewidth(np.array(hb) - SS - AT, cotasZ, PERF[i])  # Trabajamos en nmm local
            })

    if 'RSLR' in DYN[0]:
        for i in range(len(PERF)):
            id = dinperf[i]
            DYNP[i]['RSLR'] = DYN[id].RSLR[poscalc]
    
    return DYNP, poscalc