def corrige_instante_inicial_escollera(YLTi, ACT, EA, t, it):
    if not ACT:
        return YLTi

    camposEA = EA.keys()
    estrig = 'drig' in camposEA

    if estrig:  # Hay escolleras
        indrig = EA['drig']
        
        for i in range(len(indrig)):
            ACTi = ACT[indrig[i]]

            # Es su instante inicial
            if t[it] == ACTi['Tini']:
                # Es su instante inicial, verificamos si alguna posición
                # retrocede por detrás de la estructura
                pafc = ACTi['PERF']  # Perfiles afectados
                yafc = ACTi['YAFC']  # Posiciones en los perfiles afectados
                posret = [idx for idx, val in enumerate(YLTi[pafc]) if val < yafc[idx]]  # Posiciones que retroceden

                # Actualizamos
                for idx in posret:
                    YLTi[pafc[idx]] = yafc[idx]

    return YLTi