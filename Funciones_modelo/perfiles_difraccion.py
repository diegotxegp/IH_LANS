import numpy as np

def perfiles_difraccion(ACT0, PERF):
    xon = [p['xon'] for p in PERF]
    yon = [p['yon'] for p in PERF]
    xof = [p['xof'] for p in PERF]
    yof = [p['yof'] for p in PERF]
    pafc = []

    if ACT0['Tipo'] == 'dper':
        xest = ACT0['X'][-1]
        yest = ACT0['Y'][-1]
        d1 = [np.hypot(x - xest, y - yest) for x, y in zip(xon, yon)]
        d2 = [np.hypot(x - xest, y - yest) for x, y in zip(xof, yon)]
        posmin1 = np.argmin(d1)
        
        if inpolygon(xest, yest, [xon[posmin1], xon[posmin1+1], xof[posmin1+1], xof[posmin1], xon[posmin1]],
                     [yon[posmin1], yon[posmin1+1], yof[posmin1+1], yof[posmin1], yon[posmin1]]):
            pafc = [posmin1, posmin1+1]
        elif inpolygon(xest, yest, [xon[posmin1-1], xon[posmin1], xof[posmin1], xof[posmin1-1], xon[posmin1-1]],
                       [yon[posmin1-1], yon[posmin1], yof[posmin1], yof[posmin1-1], yon[posmin1-1]]):
            pafc = [posmin1-1, posmin1]
        else:
            count = 0
            noin = 0
            while noin == 0:
                count = count + 1
                noin = inpolygon(xest, yest, [xon[count], xon[count+1], xof[count+1], xof[count], xon[count]],
                                 [yon[count], yon[count+1], yof[count+1], yof[count], yon[count]])
                pafc = [count, count+1]
    elif ACT0['Tipo'] == 'dperpar':
        pact = np.zeros((2, 2))
        for i in range(2):
            ind = i - 2
            xest = ACT0['X'][len(ACT0['X']) + ind]
            yest = ACT0['Y'][len(ACT0['Y']) + ind]
            d1 = [np.hypot(x - xest, y - yest) for x, y in zip(xon, yon)]
            posmin1 = np.argmin(d1)
            
            if inpolygon(xest, yest, [xon[posmin1], xon[posmin1+1], xof[posmin1+1], xof[posmin1], xon[posmin1]],
                         [yon[posmin1], yon[posmin1+1], yof[posmin1+1], yof[posmin1], yon[posmin1]]):
                pact[:, i] = [posmin1, posmin1+1]
            elif inpolygon(xest, yest, [xon[posmin1-1], xon[posmin1], xof[posmin1], xof[posmin1-1], xon[posmin1-1]],
                           [yon[posmin1-1], yon[posmin1], yof[posmin1], yof[posmin1-1], yon[posmin1-1]]):
                pact[:, i] = [posmin1-1, posmin1]
            else:
                count = 0
                noin = 0
                while noin == 0:
                    count = count + 1
                    noin = inpolygon(xest, yest, [xon[count], xon[count+1], xof[count+1], xof[count], xon[count]],
                                     [yon[count], yon[count+1], yof[count+1], yof[count], yon[count]])
                    pact[:, i] = [count, count+1]
        pafc = pact.tolist()
    elif ACT0['Tipo'] == 'dpar':
        pact = np.zeros((2, 2))
        for i in range(len(ACT0['X'])):
            xest = ACT0['X'][i]
            yest = ACT0['Y'][i]
            d1 = [np.hypot(x - xest, y - yest) for x, y in zip(xon, yon)]
            d2 = [np.hypot(x - xest, y - yest) for x, y in zip(xof, yon)]
            posmin1 = np.argmin(d1)
            
            if inpolygon(xest, yest, [xon[posmin1], xon[posmin1+1], xof[posmin1+1], xof[posmin1], xon[posmin1]],
                         [yon[posmin1], yon[posmin1+1], yof[posmin1+1], yof[posmin1], yon[posmin1]]):
                pact[:, i] = [posmin1, posmin1+1]
            elif inpolygon(xest, yest, [xon[posmin1-1], xon[posmin1], xof[posmin1], xof[posmin1-1], xon[posmin1-1]],
                           [yon[posmin1-1], yon[posmin1], yof[posmin1], yof[posmin1-1], yon[posmin1-1]]):
                pact[:, i] = [posmin1-1, posmin1]
            else:
                count = 0
                noin = 0
                while noin == 0:
                    count = count + 1
                    noin = inpolygon(xest, yest, [xon[count], xon[count+1], xof[count+1], xof[count], xon[count]],
                                     [yon[count], yon[count+1], yof[count+1], yof[count], yon[count]])
                    pact[:, i] = [count, count+1]
        pafc = pact.tolist()
    else:
        pafc = 0
    
    return pafc