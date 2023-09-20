import numpy as np

def filtro_kalman_longitudinal_transversal(estado_ant, Jacobito_lc, DALCS, it, posero, Ylt, Yct):
    n = estado_ant.shape[1]
    saltoYlt = np.zeros(n)
    saltoYct = np.zeros(n)
    Yctn = np.zeros(Yct.shape)
    estado_post = np.zeros_like(estado_ant)

    Hero = np.array([1, 0, 0, 1, 0, 0, 0, 0])
    Hacr = np.array([1, 0, 0, 0, 0, 1, 0, 0])

    for i in range(n):
        if posero[i] == 1:
            P0 = DALCS[i]["Pero0_lc"]
            Q = DALCS[i]["Qero_lc"]
            P = np.dot(np.dot(Jacobito_lc[:, :, i], P0), Jacobito_lc[:, :, i].T) + Q
            DALCS[i]["Pero0_lc"] = P
            DALCS[i]["Pacr0_lc"][8, :] = P[8, :]
            DALCS[i]["Pacr0_lc"][:, 8] = P[:, 8]
            DALCS[i]["Pacr0_lc"][0:3, :] = P[0:3, :]
            DALCS[i]["Pacr0_lc"][:, 0:3] = P[:, 0:3]
            H = Hero
            HM = np.array([1, 1, 1, 1, 1, 0, 1, 1])
        else:
            P0 = DALCS[i]["Pacr0_lc"]
            Q = DALCS[i]["Qacr_lc"]
            P = np.dot(np.dot(Jacobito_lc[:, :, i], P0), Jacobito_lc[:, :, i].T) + Q
            DALCS[i]["Pacr0_lc"] = P
            DALCS[i]["Pero0_lc"][8, :] = P[8, :]
            DALCS[i]["Pero0_lc"][:, 8] = P[:, 8]
            DALCS[i]["Pero0_lc"][0:3, :] = P[0:3, :]
            DALCS[i]["Pero0_lc"][:, 0:3] = P[:, 0:3]
            H = Hacr
            HM = np.array([1, 1, 1, 0, 1, 1, 1, 1])

        contador = DALCS[i]["pos_lc"]

        if it + 1 == DALCS[i]["nasim"][contador] and DALCS[i]["stop_lc"] == 0:
            if it + 1 >= DALCS[i]["itmax"]:
                DALCS[i]["stop_lc"] = 1
            else:
                DALCS[i]["pos_lc"] = contador + 1

            Yobs = DALCS[i]["YY"][contador]
            K = np.dot(np.dot(P, H.T), np.linalg.inv(np.dot(np.dot(H, P), H.T) + DALCS[i]["R_lc"]))
            estado_ant[1, i] = np.log(estado_ant[1, i] / DALCS[i]["kcerc0"]) / DALCS[i]["sigmaK"]
            modificacionKalman = np.dot(K, (Yobs - np.dot(H, estado_ant[:, i])))

            estado_post[:, i] = HM * (estado_ant[:, i] + modificacionKalman)
            estado_post[1, i] = DALCS[i]["kcerc0"] * np.exp(DALCS[i]["sigmaK"] * estado_post[1, i])

            Pnew = np.dot(np.dot(np.eye(len(H)) - np.dot(np.dot(K, H), P), (1.0 / len(H)) - K), P)
            saltoYlt[i] = estado_post[0, i] - Ylt[i]

            if posero[i] == 1:
                DALCS[i]["Pero0_lc"] = Pnew
                DALCS[i]["Pacr0_lc"][8, :] = Pnew[8, :]
                DALCS[i]["Pacr0_lc"][:, 8] = Pnew[:, 8]
                DALCS[i]["Pacr0_lc"][0:3, :] = Pnew[0:3, :]
                DALCS[i]["Pacr0_lc"][:, 0:3] = Pnew[:, 0:3]
                saltoYct[i] = estado_post[3, i] - Yct[i]
            else:
                DALCS[i]["Pacr0_lc"] = Pnew
                DALCS[i]["Pero0_lc"][8, :] = Pnew[8, :]
                DALCS[i]["Pero0_lc"][:, 8] = Pnew[:, 8]
                DALCS[i]["Pero0_lc"][0:3, :] = Pnew[0:3, :]
                DALCS[i]["Pero0_lc"][:, 0:3] = Pnew[:, 0:3]
                saltoYct[i] = estado_post[5, i] - Yct[i]
        else:
            estado_post[:, i] = HM * estado_ant[:, i]

    return estado_post, DALCS, saltoYlt, saltoYct