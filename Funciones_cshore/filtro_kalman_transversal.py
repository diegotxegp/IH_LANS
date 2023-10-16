import numpy as np

def filtro_kalman_transversal(estado_ant, Jacobito, DACS, it, pero, Yct):
    Hero = np.array([1, 0, 0, 0, 0])
    Hacr = np.array([0, 0, 1, 0, 0])
    nel = estado_ant.shape[0]
    saltoYct = np.zeros(Jacobito.shape[2])
    Yctn = np.zeros(Yct.shape)
    estado_post = np.zeros(estado_ant.shape)

    for i in range(len(DACS)):
        if pero[i] == 1:
            P = np.dot(np.dot(Jacobito[i][:][:], DACS[i]["Pero0_c"]), np.transpose(Jacobito[i][:][:])) + DACS[i]["Qero_c"]
            DACS[i]["Pero0_c"] = P
            DACS[i]["Pacr0_c"][4][:] = P[4][:]
            DACS[i]["Pacr0_c"][:][4] = P[:][4]
            Q = DACS[i]["Qero_c"]
            H = Hero
            HM = np.array([1, 1, 0, 1, 1])
        else:
            P = np.dot(np.dot(Jacobito[i][:][:], DACS[i]["Pacr0_c"]), np.transpose(Jacobito[i][:][:])) + DACS[i]["Qacr_c"]
            DACS[i]["Pacr0_c"] = P
            DACS[i]["Pero0_c"][4, :] = P[4, :]
            DACS[i]["Pero0_c"][:, 4] = P[:, 4]
            Q = DACS[i]["Qacr_c"]
            H = Hacr
            HM = np.array([0, 1, 1, 1, 1])

        contador = DACS[i]["pos_c"]

        if "nasim" in DACS[i]:
            if it + 1 == DACS[i]["nasim"][contador] and DACS[i]["stop_c"] == 0:
                if it + 1 >= DACS[i]["itmax"]:
                    DACS[i]["stop_c"] = 1
                else:
                    DACS[i]["pos_c"] = contador + 1

                Yobs = DACS[i]["Yct"][contador]
                K = np.dot(np.dot(P, H.T), np.linalg.inv(np.dot(np.dot(H, P), H.T) + DACS[i]["R_c"]))
                modificacionKalman = np.dot(K, Yobs - np.dot(H, estado_ant[:][i]))
                estado_post[:][i] = HM * (estado_ant[:][i] + modificacionKalman)

                if pero[i] == 1:
                    DACS[i]["Pero0_c"] = np.dot(np.eye(len(H)) - np.dot(K, H), P)
                    DACS[i]["Pacr0_c"][4][:] = DACS[i]["Pero0_c"][4][:]
                    DACS[i]["Pacr0_c"][:][4] = DACS[i]["Pero0_c"][:][4]
                    saltoYct[i] = estado_post[0, i] - Yct[i]
                else:
                    DACS[i]["Pacr0_c"] = np.dot(np.eye(len(H)) - np.dot(K, H), P)
                    DACS[i]["Pero0_c"][4][:] = DACS[i]["Pacr0_c"][4][:]
                    DACS[i]["Pero0_c"][:][4] = DACS[i]["Pacr0_c"][:][4]
                    saltoYct[i] = estado_post[2][i] - Yct[i]
            else:
                estado_post[:,i] = HM * estado_ant[:, i]
        else:
            estado_post[:,i] = HM * estado_ant[:,i]

    return estado_post, DACS, saltoYct