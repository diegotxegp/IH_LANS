import numpy as np

def estructuras_clasifica(t, ACT):
    # Función para detectar las estructuras activas y clasificarlas
    index = np.zeros(len(ACT))
    for i in range(len(ACT)):
        t0e = ACT[i]['Tini']
        tfe = ACT[i]['Tfin']
        if t >= t0e and t < tfe:
            index[i] = 1
    ind = np.where(index)[0]

    # Clasificamos las estructuras activas
    EA = {}
    if 'dper' in [ACT[i]['Tipo'] for i in ind]:
        EA['dper'] = [i for i in ind if ACT[i]['Tipo'] == 'dper']
    if 'dpar' in [ACT[i]['Tipo'] for i in ind]:
        EA['dpar'] = [i for i in ind if ACT[i]['Tipo'] == 'dpar']
    if 'dperpar' in [ACT[i]['Tipo'] for i in ind]:
        EA['dperpar'] = [i for i in ind if ACT[i]['Tipo'] == 'dperpar']
    if 'drig' in [ACT[i]['Tipo'] for i in ind]:
        EA['drig'] = [i for i in ind if ACT[i]['Tipo'] == 'drig']
    if 'nour' in [ACT[i]['Tipo'] for i in ind]:
        EA['nour'] = [i for i in ind if ACT[i]['Tipo'] == 'nour']
    if not ind:
        EA['noest'] = []

    return EA