import numpy as np

def calcula_lc(YLTi, dt, Dc, kcerc, dQdx, Q, dx, ACT, EA):

    # no hay estructuras que la hacen cambiar
    if len(ACT) > 0:
        camposEA = list(EA.keys())
        estpar = 'dpar' in camposEA
        estrig = 'drig' in camposEA
        estcamb = [estpar, estrig]

        if any(estcamb):  # existen estructuras/actuaciones que cambian LC
            # avanzamos la línea de costa
            ylt0 = YLTi - dt / Dc * kcerc * dQdx
            # calculamos la posición de la lc
            # comprobamos si existe intersección con la estructura
            # Estructuras tipo drig y nour
            if estpar:
                indpar = EA['dpar']
                for i in range(len(indpar)):
                    ACTi = ACT[indpar[i]]
                    ylt = ylt0.copy()

            if estrig:
                indrig = EA['drig']
                for i in range(len(indrig)):
                    ACTi = ACT[indrig[i]]
                    pafc = ACTi['PERF']
                    # verificamos si la lc retrocede por detrás de la escollera
                    posretant = np.where(YLTi[pafc] < ACTi['YAFC'])[0]
                    if len(posretant) > 0:
                        YLTi[pafc[posretant]] = ACTi['YAFC'][posretant]

                    posret = np.where(ylt0[pafc] < ACTi['YAFC'])[0]
                    if len(posret) > 0:  # se viola la condición de la escollera
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
                                # actualizamos caudales y recalculamos
                                Q[posreti + 1] = Qnext
                                Q[posreti] = Qant
                                dQdxi = (Q[1:] - Q[:-1]) / dx
                                yltact = YLTi - dt / Dc * kcerc * dQdxi
                                # ahora posreti cumple, el siguiente?
                                posretnext = posreti + 1
                                posretant = posreti - 1
                                afcnext = afci + 1
                                afcant = afci - 1
                                nocumpnext = np.where(yltact[posretnext] < ACTi['YAFC'][afcnext])[0]
                                nocumpant = np.where(yltact[posretant] < ACTi['YAFC'][afcant])[0]
                                # corregimos nocumpnext
                                while len(nocumpnext) > 0:
                                    # tte a derechas
                                    posreti = posretnext
                                    afci = afcnext
                                    Qnext = Q[posreti + 1]
                                    ysi = ACTi['YAFC'][afci]
                                    yi = YLTi[posreti]
                                    Qant = Q[posreti]
                                    Dci = Dc[posreti]
                                    dxi = dx[posreti]
                                    kcerci = kcerc[posreti]
                                    if Qnext >= 0 and Qant >= 0:  # a derechas (debería ser este siempre)
                                        posretnext = posreti + 1
                                        afcnext = afci + 1
                                        Qnext_ast = (yi - ysi) * Dci * dxi / kcerci / dt + Qant
                                        Q[posreti + 1] = Qnext_ast
                                    elif Qnext <= 0 and Qant <= 0:  # a izquierda
                                        posretnext = posreti - 1
                                        afcnext = afci - 1
                                        Qant_ast = (ysi - yi) * Dci * dxi / kcerci / dt + Qnext
                                        Q[posreti] = Qant_ast
                                    if posretnext in pafc:
                                        # recalculamos LC, BC ya impuestas
                                        dQdxi = (Q[1:] - Q[:-1]) / dx
                                        yltact = YLTi - dt / Dc * kcerc * dQdxi
                                        # ahora posreti cumple, el siguiente?
                                        nocumpnext = np.where(yltact[posretnext] < ACTi['YAFC'][afcnext])[0]
                                    else:  # salimos del bucle
                                        nocumpnext = []

                                while len(nocumpant) > 0:
                                    # tte a izquierdas
                                    posreti = posretant
                                    afci = afcant
                                    Qnext = Q[posreti + 1]
                                    ysi = ACTi['YAFC'][afci]
                                    yi = YLTi[posreti]
                                    Qant = Q[posreti]
                                    Dci = Dc[posreti]
                                    dxi = dx[posreti]
                                    kcerci = kcerc[posreti]
                                    if Qnext >= 0 and Qant >= 0:  # a derechas
                                        posretnext = posreti + 1
                                        afcnext = afci + 1
                                        Qnext_ast = (yi - ysi) * Dci * dxi / kcerci / dt + Qant
                                        Q[posreti + 1] = Qnext_ast
                                    elif Qnext <= 0 and Qant <= 0:  # a izquierda (debería ser este siempre)
                                        posretnext = posreti - 1
                                        afcnext = afci - 1
                                        Qant_ast = (ysi - yi) * Dci * dxi / kcerci / dt + Qnext
                                        Q[posreti] = Qant_ast
                                    if posretnext in pafc:
                                        # recalculamos LC, BC ya impuestas
                                        dQdxi = (Q[1:] - Q[:-1]) / dx
                                        yltact = YLTi - dt / Dc * kcerc * dQdxi
                                        # ahora posreti cumple, el siguiente?
                                        nocumpant = np.where(yltact[posretnext] < ACTi['YAFC'][afcnext])[0]
                                        # actualizamos el que no cumple
                                        posretant = posretnext
                                        afcant = afcnext
                                    else:  # salimos del bucle
                                        nocumpant = []
                                        yltact = YLTi - dt / Dc * kcerc * (Q[1:] - Q[:-1]) / dx
                                ylt = yltact
                            else:  # no es celda vaciante de dos lados
                                c = 0
                                while len(nocump) > 0:
                                    Qnext = Q[posreti + 1]
                                    ysi = ACTi['YAFC'][afci]
                                    yi = YLTi[posreti]
                                    Qant = Q[posreti]
                                    Dci = Dc[posreti]
                                    dxi = dx[posreti]
                                    kcerci = kcerc[posreti]
                                    if Qnext >= 0 and Qant >= 0:  # a derechas
                                        posretnext = posreti + 1
                                        afcnext = afci + 1
                                        Qnext_ast = (yi - ysi) * Dci * dxi / kcerci / dt + Qant
                                        Q[posreti + 1] = Qnext_ast
                                    elif Qnext <= 0 and Qant <= 0:  # a izquierda
                                        posretnext = posreti - 1
                                        afcnext = afci - 1
                                        Qant_ast = (ysi - yi) * Dci * dxi / kcerci / dt + Qnext
                                        Q[posreti] = Qant_ast
                                    elif afci == 0 and Qnext > 0 and Qant < 0:  # vaciante primer perfil
                                        # corrección hacia delante
                                        posretnext = posreti + 1
                                        afcnext = afci + 1
                                        # tocamos los dos caudales
                                        q2q1ast = (ysi - yi) / dt * -Dci / kcerci * dxi
                                        q2q1 = Qnext - Qant
                                        Qnext = q2q1ast / q2q1 * Qnext
                                        Qant = q2q1ast / q2q1 * Qant
                                        # actualizamos caudales y recalculamos
                                        Q[posreti + 1] = Qnext
                                        Q[posreti] = Qant
                                    elif afci == len(ACTi['PERF']) - 1 and Qnext > 0 and Qant < 0:  # vaciante ultimo perfil
                                        # corrección hacia atrás
                                        posretnext = posreti - 1
                                        afcnext = afci - 1
                                        # tocamos los dos caudales
                                        q2q1ast = (ysi - yi) / dt * -Dci / kcerci * dxi
                                        q2q1 = Qnext - Qant
                                        Qnext = q2q1ast / q2q1 * Qnext
                                        Qant = q2q1ast / q2q1 * Qant
                                        # actualizamos caudales y recalculamos
                                        Q[posreti + 1] = Qnext
                                        Q[posreti] = Qant
                                    elif Qnext < 0 and Qant > 0:  # celda llenante
                                        posretnext = []  # no se propaga más
                                    if posretnext in pafc:
                                        # recalculamos LC, BC ya impuestas
                                        dQdxi = (Q[1:] - Q[:-1]) / dx
                                        yltact = YLTi - dt / Dc * kcerc * dQdxi
                                        # ahora posreti cumple, el siguiente?
                                        nocump = np.where(yltact[posretnext] < ACTi['YAFC'][afcnext])[0]
                                        # actualizamos el que no cumple
                                        posreti = posretnext
                                        afci = afcnext
                                    else:  # salimos del bucle
                                        nocump = []
                                        yltact = YLTi - dt / Dc * kcerc * (Q[1:] - Q[:-1]) / dx
                                ylt = yltact
                        ylt = ylt0
            else:  # no se viola la condición de escollera
                ylt = ylt0
        else:
            ylt = YLTi - dt / Dc * kcerc * dQdx
    else:  # no hay estructuras
        ylt = YLTi - dt / Dc * kcerc * dQdx

    return ylt, kcerc