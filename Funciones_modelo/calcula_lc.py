import numpy as np

def calcula_lc(YLTi, dt, Dc, kcerc, dQdx, Q, dx, ACT, EA):
    if EA is not None:
        camposEA = EA.keys()
        estpar = 'dpar' in camposEA
        estrig = 'drig' in camposEA
        estcamb = [estpar, estrig]

        if any(estcamb):  # Existencia de estructuras/actuaciones que cambian la línea de costa
            ylt0 = YLTi - dt / Dc * kcerc * dQdx

            if estpar:
                indpar = EA['dpar']
                for i in range(len(indpar)):
                    ACTi = ACT[indpar[i]]
                    ylt = ylt0.copy()  # Copiar ylt0 para trabajar con él

            if estrig:
                indrig = EA['drig']
                for i in range(len(indrig)):
                    ACTi = ACT[indrig[i]]
                    pafc = ACTi['PERF']

                    # Verificar si la línea de costa retrocede detrás de la estructura
                    posretant = np.where(YLTi[pafc] < ACTi['YAFC'])[0]

                    if len(posretant) > 0:
                        YLTi[pafc[posretant]] = ACTi['YAFC'][posretant]

                    posret = np.where(ylt0[pafc] < ACTi['YAFC'])[0]

                    if len(posret) > 0:  # Se viola la condición de la estructura
                        for iv in range(len(posret)):
                            posreti = pafc[posret[iv]]
                            afci = posret[iv]
                            nocump = 1
                            Qnext = Q[posreti + 1]
                            ysi = ACTi['YAFC'][afci]
                            yi = YLTi[posreti]
                            Qant = Q[posreti]
                            Dci = Dc[posreti]
                            dxi = dx[posreti]
                            kcerci = kcerc[posreti]

                            if Qnext > 0 and Qant < 0 and posreti > pafc[0] and posreti + 1 <= pafc[-1]:
                                q2q1ast = (ysi - yi) / dt * -Dci / kcerci * dxi
                                q2q1 = Qnext - Qant
                                Qnext = q2q1ast / q2q1 * Qnext
                                Qant = q2q1ast / q2q1 * Qant

                                Q[posreti + 1] = Qnext
                                Q[posreti] = Qant
                                dQdxi = (Q[1:] - Q[:-1]) / dx
                                yltact = YLTi - dt / Dc * kcerc * dQdxi

                                posretnext = posreti + 1
                                posretant = posreti - 1
                                afcnext = afci + 1
                                afcant = afci - 1

                                nocumpnext = np.where(yltact[posretnext] < ACTi['YAFC'][afcnext])[0]
                                nocumpant = np.where(yltact[posretant] < ACTi['YAFC'][afcant])[0]

                                while len(nocumpnext) > 0:
                                    posreti = posretnext
                                    afci = afcnext
                                    Qnext = Q[posreti + 1]
                                    ysi = ACTi['YAFC'][afci]
                                    yi = YLTi[posreti]
                                    Qant = Q[posreti]
                                    Dci = Dc[posreti]
                                    dxi = dx[posreti]
                                    kcerci = kcerc[posreti]

                                    if Qnext >= 0 and Qant >= 0:
                                        posretnext = posreti + 1
                                        afcnext = afci + 1
                                        Qnext_ast = (yi - ysi) * Dci * dxi / kcerci / dt + Qant
                                        Q[posreti + 1] = Qnext_ast
                                    elif Qnext <= 0 and Qant <= 0:
                                        posretnext = posreti - 1
                                        afcnext = afci - 1
                                        Qant_ast = (ysi - yi) * Dci * dxi / kcerci / dt + Qnext
                                        Q[posreti] = Qant_ast

                                    if posretnext in pafc:
                                        dQdxi = (Q[1:] - Q[:-1]) / dx
                                        yltact = YLTi - dt / Dc * kcerc * dQdxi
                                        nocumpnext = np.where(yltact[posretnext] < ACTi['YAFC'][afcnext])[0]
                                    else:
                                        nocumpnext = []

                                while len(nocumpant) > 0:
                                    posreti = posretant
                                    afci = afcant
                                    Qnext = Q[posreti + 1]
                                    ysi = ACTi['YAFC'][afci]
                                    yi = YLTi[posreti]
                                    Qant = Q[posreti]
                                    Dci = Dc[posreti]
                                    dxi = dx[posreti]
                                    kcerci = kcerc[posreti]

                                    if Qnext >= 0 and Qant >= 0:
                                        posretnext = posreti + 1
                                        afcnext = afci + 1
                                        Qnext_ast = (yi - ysi) * Dci * dxi / kcerci / dt + Qant
                                        Q[posreti + 1] = Qnext_ast
                                    elif Qnext <= 0 and Qant <= 0:
                                        posretnext = posreti - 1
                                        afcnext = afci - 1
                                        Qant_ast = (ysi - yi) * Dci * dxi / kcerci / dt + Qnext
                                        Q[posreti] = Qant_ast

                                    if posretnext in pafc:
                                        dQdxi = (Q[1:] - Q[:-1]) / dx
                                        yltact = YLTi - dt / Dc * kcerc * dQdxi
                                        nocumpant = np.where(yltact[posretnext] < ACTi['YAFC'][afcnext])[0]
                                        posretant = posretnext
                                        afcant = afcnext
                                    else:
                                        nocumpant = []

                                ylt = yltact
                            else:
                                c = 0
                                nocump = 1

                                while len(nocump) > 0:
                                    Qnext = Q[posreti + 1]
                                    ysi = ACTi['YAFC'][afci]
                                    yi = YLTi[posreti]
                                    Qant = Q[posreti]
                                    Dci = Dc[posreti]
                                    dxi = dx[posreti]
                                    kcerci = kcerc[posreti]

                                    if Qnext >= 0 and Qant >= 0:
                                        posretnext = posreti + 1
                                        afcnext = afci + 1
                                        Qnext_ast = (yi - ysi) * Dci * dxi / kcerci / dt + Qant
                                        Q[posreti + 1] = Qnext_ast
                                    elif Qnext <= 0 and Qant <= 0:
                                        posretnext = posreti - 1
                                        afcnext = afci - 1
                                        Qant_ast = (ysi - yi) * Dci * dxi / kcerci / dt + Qnext
                                        Q[posreti] = Qant_ast
                                    elif afci == 0 and Qnext > 0 and Qant < 0:
                                        posretnext = posreti + 1
                                        afcnext = afci + 1
                                        q2q1ast = (ysi - yi) / dt * -Dci / kcerci * dxi
                                        q2q1 = Qnext - Qant
                                        Qnext = q2q1ast / q2q1 * Qnext
                                        Qant = q2q1ast / q2q1 * Qant

                                        Q[posreti + 1] = Qnext
                                        Q[posreti] = Qant
                                    elif afci == len(ACTi['PERF']) - 1 and Qnext > 0 and Qant < 0:
                                        posretnext = posreti - 1
                                        afcnext = afci - 1
                                        q2q1ast = (ysi - yi) / dt * -Dci / kcerci * dxi
                                        q2q1 = Qnext - Qant
                                        Qnext = q2q1ast / q2q1 * Qnext
                                        Qant = q2q1ast / q2q1 * Qant

                                        Q[posreti + 1] = Qnext
                                        Q[posreti] = Qant
                                    elif Qnext < 0 and Qant > 0:
                                        posretnext = []  # No se propaga más

                                    if posretnext in pafc:
                                        dQdxi = (Q[1:] - Q[:-1]) / dx
                                        yltact = YLTi - dt / Dc * kcerc * dQdxi
                                        nocump = np.where(yltact[posretnext] < ACTi['YAFC'][afcnext])[0]

                                        posreti = posretnext
                                        afci = afcnext
                                    else:
                                        nocump = []

                                ylt = yltact
                    else:  # No se viola la condición de la estructura
                        ylt = ylt0
            else:
                ylt = YLTi - dt / Dc * kcerc * dQdx
        else:
            ylt = YLTi - dt / Dc * kcerc * dQdx
    else:  # No hay estructuras
        ylt = YLTi - dt / Dc * kcerc * dQdx

    return ylt