import numpy as np

def filtro_kalman_longitudinal(estado_ant, Jacobito, DA, it):
    H = np.array([1, 0, 0])
    n = estado_ant.shape[1]
    estado_post = np.zeros_like(estado_ant)

    for i in range(n):
        DA[i]["P0"] = np.dot(np.dot(Jacobito[:, :, i], DA[i]["P0"]), np.transpose(Jacobito[:, :, i])) + DA[i]["Q"]
        contador = DA[i]["pos_l"]

        if it + 1 == DA[i]["nasim"][contador] and DA[i]["stop_l"] == 0:
            if it + 1 >= DA[i]["itmax"]:
                DA[i]["stop_l"] = 1
            else:
                DA[i]["pos_l"] = contador + 1

            Yobs = DA[i]["Ylt"][contador]
            K = DA[i]["P0"]*H.T / (H*DA[i]["P0"]*H.T + DA[i]["R"])
            modificacionKalman = np.dot(K, (Yobs - np.dot(H, estado_ant[:, i])))
            estado_ant_m = estado_ant[:, i]

            estado_ant_m[1] = np.log(estado_ant[1, i] / DA[i]["kcerc0"]) / DA[i]["sigmaK"]
            estado_post[:, i] = estado_ant_m + modificacionKalman
            estado_post[1, i] = DA[i]["kcerc0"] * np.exp(DA[i]["sigmaK"] * estado_post[1, i])
            DA[i]["P0"] = (np.eye(len(H)) - K*H) * DA[i]["P0"]

        else:
            estado_post[:, i] = estado_ant[:, i]

    return estado_post, DA