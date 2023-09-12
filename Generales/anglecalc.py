import numpy as np

def anglecalc(sinn, coss):
    angle = np.zeros(sinn.size)

    for i in range(len(sinn)):
        if coss[i] > 0 and sinn[i] > 0:
            angle[i] = np.degrees(np.arctan(sinn[i] / coss[i]))
        elif coss[i] < 0 and sinn[i] > 0:
            angle[i] = np.degrees(np.arctan(sinn[i] / coss[i])) + 180
        elif coss[i] < 0 and sinn[i] < 0:
            angle[i] = np.degrees(np.arctan(sinn[i] / coss[i])) + 180
        elif coss[i] > 0 and sinn[i] < 0:
            angle[i] = np.degrees(np.arctan(sinn[i] / coss[i])) + 360
        elif coss[i] == 0 and sinn[i] > 0:
            angle[i] = 90
        elif coss[i] == 0 and sinn[i] < 0:
            angle[i] = 270
        elif sinn[i] == 0 and coss[i] > 0:
            angle[i] = 0
        elif sinn[i] == 0 and coss[i] < 0:
            angle[i] = 180
        else:
            angle[i] = np.nan

    return angle.tolist()