import numpy as np

def ensambla_lc(Jacobito_l, Jacobito_c):
    dim1 = Jacobito_l.shape[0] + Jacobito_c.shape[0]
    dim2 = Jacobito_l.shape[1] + Jacobito_c.shape[1]
    dim3 = Jacobito_l.shape[2]
    Jacobito_lc = np.zeros((dim1, dim2, dim3))

    for i in range(dim3):
        Jacobito_lc[0:Jacobito_l.shape[0], 0:Jacobito_l.shape[1], i] = Jacobito_l[:, :, i]
        Jacobito_lc[Jacobito_l.shape[0]:dim1, Jacobito_l.shape[1]:dim2, i] = Jacobito_c[:, :, i]

    return Jacobito_lc