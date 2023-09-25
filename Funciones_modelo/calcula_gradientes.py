import numpy as np

def calcula_gradientes(Q, dx, PNL, PBC, bctype, bctypeval):
    Qall = np.empty(len(dx) + 1)
    Qall[1:-1] = Q

    # Aplicar condiciones de contorno
    ntramos = len(PBC)
    for i in range(ntramos):
        Tipoi = bctype[i]
        PBCi = PBC[i]
        vali = bctypeval[i]

        if Tipoi[0] == 'Dirichlet':  # En caudales frontera
            Qall[PBCi[0]] = vali[0]

        if Tipoi[1] == 'Dirichlet':
            Qall[PBCi[1] + 1] = vali[1]

    Qbc = Qall.copy()
    dQdx = (Qall[1:] - Qall[:-1]) / dx

    for i in range(ntramos):
        Tipoi = bctype[i]
        PBCi = PBC[i]
        vali = bctypeval[i]

        if Tipoi[0] == 'Neumann':  # En gradientes
            dQdx[PBCi[0]] = vali[0]
            Qbc[PBCi[0]] = -dQdx[PBCi[0]] * dx[PBCi[0]] + Qbc[PBCi[0] + 1]

        if Tipoi[1] == 'Neumann':
            dQdx[PBCi[1]] = vali[1]
            Qbc[PBCi[1] + 1] = dQdx[PBCi[1]] * dx[PBCi[1]] + Qbc[PBCi[1]]

    # Anular gradientes en perfiles no lcs
    dQdx[PNL] = 0

    # Corregir Qbc con los límites de PBC
    Qbc[:PBC[0, 0]] = 0
    Qbc[PBC[-1, 1] + 2:] = 0

    if ntramos > 1:
        for j in range(ntramos - 1):
            Qbc[PBC[j, 1] + 2:PBC[j + 1, 0] - 1] = 0

    return dQdx, Qbc