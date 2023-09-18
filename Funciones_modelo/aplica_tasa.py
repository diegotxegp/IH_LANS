import numpy as np

def aplica_tasa(Yltii, ACT, EA, vlt, RSLR):
    # Suma tasa nourishments
    Yltnour = np.zeros(len(Yltii))

    if EA:
        camposEA = EA.keys()
        estnour = 'nour' in camposEA

        if estnour:
            for i in range(len(EA['nour'])):
                ACTi = ACT[EA['nour'][i]]
                tasanour = ACTi['tasadt']
                Yltnour = Yltnour + tasanour

    # Suma tasa permanente y de RSLR
    Ylt = Yltii + Yltnour + vlt + RSLR

    return Ylt