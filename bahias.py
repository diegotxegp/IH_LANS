import os
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

from Generales.wmoore_calc import wmoore_calc
from Funciones_modelo.calcula_linea_xy import calcula_linea_xy

## SET-UP DEL MODELO

# Tiempo del análisis
t0 = datetime(1995, 1, 1)
tfin = datetime(2020, 1, 1)
dt = 1 # Intervalo de tiempo en horas

# Lista de tiempos entre ambas fecha cada hora
t = [i for i in range(0, int((tfin - t0).total_seconds() / 3600) + 1)]


## Perfiles
nperf = 100
d50 = [0.5e-3]
dtierra = 200
dc = 10
dmax = 10
berma = 1
Adean = 0.51 * wmoore_calc(d50) ** 0.44

# definición LC
# calculamos LC inicial
# aqui indica la geometría de los perfiles
Lbahia1 = 2000
Hbahia1 = 500
A1 = 4 * Hbahia1 / (Lbahia1 ** 2)

Lbahia2 = 500
Hbahia2 = 150
A2 = 4 * Hbahia2 / (Lbahia2 ** 2)

Lbahia3 = 1000
Hbahia3 = 250
A3 = 4 * Hbahia3 / (Lbahia3 ** 2)

# % hacemos bahias independientes y concatenamos
dint = np.linspace(0, Lbahia1 + 100 + Lbahia2 + 100 + Lbahia3, nperf)
PERF = []

for i in range(nperf):
    profile = {}
    profile["xon"] = 0 + 1 * dint[i]
    profile["yon"] = 0 + 0 * dint[i]

    if dint[i] <= Lbahia1:
        xc = profile["xon"] - Lbahia1 / 2
        dparab = A1 * abs(xc) ** 2
        profile["Tipo"] = 'lcs'
    elif Lbahia1 < dint[i] <= Lbahia1 + 100:
        dparab = Hbahia1
        profile["Tipo"] = 's'
    elif Lbahia1 + 100 < dint[i] <= Lbahia1 + 100 + Lbahia2:
        xc = dint[i] - Lbahia1 - 100 - Lbahia2 / 2
        dparab = A2 * abs(xc) ** 2 + Hbahia1 - Hbahia2
        profile["Tipo"] = 'cs'
    elif Lbahia1 + 100 + Lbahia2 < dint[i] <= Lbahia1 + 100 + Lbahia2 + 100:
        profile["Tipo"] = 'ac'
        dparab = Hbahia1
    elif Lbahia1 + 100 + Lbahia2 + 100 < dint[i] <= Lbahia1 + 100 + Lbahia2 + 100 + Lbahia3:
        xc = dint[i] - Lbahia1 - 100 - Lbahia2 - 100 - Lbahia3 / 2
        dparab = A3 * abs(xc) ** 2 + Hbahia1 - Hbahia3
        profile["Tipo"] = 'lcs'

    profile["xof"] = profile["xon"] - 0 * 1000
    profile["yof"] = profile["yon"] + 1 * 1000
    nx0 = profile["xof"] - profile["xon"]
    ny0 = profile["yof"] - profile["yon"]
    profile["nx"] = nx0 / np.hypot(nx0, ny0)
    profile["ny"] = ny0 / np.hypot(nx0, ny0)
    profile["yc"] = dparab
    profile["dc"] = dc
    profile["d50"] = d50
    profile["Adean"] = Adean
    profile["Berma"] = berma
    profile["Yberma"] = dtierra
    profile["nbati"] = 0  # normal saliente en nauticas

    PERF.append(profile)

# Cálculo de la línea
xlc, ylc = calcula_linea_xy(PERF,[PERF.yc])
for profile in PERF:
    xlc.append(profile["xon"] + profile["nx"] * profile["yc"])
    ylc.append(profile["yon"] + profile["ny"] * profile["yc"])

# Gráficos
plt.plot(xlc, ylc)
plt.hold(True)