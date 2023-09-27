import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import scipy.io as sio
import time

from Generales.wmoore_calc import wmoore_calc
from Funciones_modelo.calcula_linea_xy import calcula_linea_xy
from Funciones_modelo.calcula_nbati import calcula_nbati
from Funciones_modelo.pinta_perfiles import pinta_perfiles
from Funciones_modelo.IH_LANS import IH_LANS

## SET-UP DEL MODELO

# Tiempo del análisis
t0 = datetime(1995, 1, 1)
tfin = datetime(2020, 1, 1)
dt = 1 # Intervalo de tiempo en horas

# Lista de tiempos entre ambas fecha cada hora
t = [i for i in range((tfin - t0).days + 1)]


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
yc = [profile['yc'] for profile in PERF]
xlc, ylc = calcula_linea_xy(PERF,yc)

PERF=calcula_nbati(xlc,ylc,PERF)

# Gráficos
plt.figure()
plt.plot(xlc, ylc)

pinta_perfiles(PERF)

plt.show()

time.sleep(2)
plt.close('all')


# Definir estructuras
ACT = []

# Oleaje y nivel
# Creamos serie de dinámicas en cabeza de perfil (batimétrica de 10 m)
# para flexibilizar las dinámicas de la cabeza de los perfiles
# cogemos diez puntos de dinámicas distribuidos uniformemente entre los perfiles
path_d0 = r'D:\IHLANS'
path_d = os.path.join(path_d0, "data_estudio1.mat")
DD = sio.loadmat(path_d)
npundyn = 1
nperfmaestro = np.round(np.linspace(nperf-1, 0, npundyn, dtype=int))
nyears = 30

DYN = []

for i in range(len(nperfmaestro)):
    id = nperfmaestro[i]
    DYN.append({
        'X': PERF[id]['xof'],
        'Y': PERF[id]['yof'],
        'Hs': np.ones(nyears * 365) * 1,
        'Tp': np.ones(nyears * 365) * 7,
        'Hs': [fila[0] for fila in DD['Hs'][:nyears * 365]],  # incluimos oleaje real
        'Dir': np.ones(nyears * 365) * 0,
        'SS': np.ones(nyears * 365) * 0,
        'AT': np.ones(nyears * 365) * 0,
        'SLR': np.ones(nyears * 365) * 0,
        'h0': 10  # calado positivo
    })

t0dyn = pd.to_datetime('01-Jan-1995', format='%d-%b-%Y')
tfin = t[0] + len(DYN[0]['Hs'])
tdyn = [i for i in range(t[0], tfin, 1)]

DYN[0]['t'] = tdyn


# Creamos estructura de datos de entrada
# Datos base
INPUT = {}
INPUT['PERF'] = PERF
INPUT['DYN'] = DYN
INPUT['ACT'] = ACT
INPUT['t'] = t
# INPUT['DA'] = DA
INPUT['dt'] = dt  # horario

# Datos propagación
INPUT['calcularotura'] = 0
INPUT['inthidromorfo'] = 0
INPUT['alpha_int'] = 0
INPUT['gamma'] = 0.50

# Datos longshore
INPUT['kcerc'] = [20]
INPUT['bctype'] = [['Dirichlet', 'Dirichlet'], ['Dirichlet', 'Dirichlet']]
INPUT['bctypeval'] = np.array([[0, 0], [0, 0]])  # Por defecto Dirichlet!! Q=0 en extremos (playa encajada)
INPUT['fcourant'] = 0.1  # Entre 0 y 1; 1 cumple courant estricto

# Datos cross-shore
INPUT['kacr'] = [3.89E-3 * 3]
INPUT['kero'] = [1.74e-2 * 8]
INPUT['tstab'] = 365 * 2  # Tiempo de estabilización (inicialización modelo)
INPUT['tcent'] = 365 * 5  # Tiempo tras el tiempo de estabilización en el que se calcula dy0

# Datos ploteo
INPUT['plott'] = 1
INPUT['tplot'] = 365 / 2
INPUT['toutp'] = 1


# Ejecutamos el modelo
RES = IH_LANS(INPUT)

# Graficamos RES.YLT[:,99]
plt.plot(RES['YLT'][:, 99])

# Graficamos RES.YLT[:,99] + RES.YCT[:,99]
plt.plot(RES['YLT'][:, 99] + RES['YCT'][:, 99])

# Graficamos RES.YCT[:,[0, 49, 99]]
plt.plot(RES['YCT'][:, [0, 49, 99]])

plt.show()


kcerc = [40]
kacr = [3.89E-3]
kero = [1.74E-2]
plott = 0

# Copia los valores de INPUT en INPUT2 y actualiza los parámetros relevantes
INPUT2 = INPUT.copy()
INPUT2['kcerc'] = kcerc
INPUT2['kacr'] = kacr
INPUT2['kero'] = kero
INPUT2['plott'] = plott

# Llama a la función IH_LANS con los nuevos parámetros
RES2 = IH_LANS(INPUT2)

# Grafica los resultados
plt.plot(RES['YLT'][:, 99], label='RES.YLT[:,99]')
plt.plot(RES2['YLT'][:, 99], label='RES2.YLT[:,99]')
plt.plot(RES['YLT'][:, 99] + RES['YCT'][:, 99], label='RES.YLT[:,99] + RES.YCT[:,99]')
plt.plot(RES2['YLT'][:, 99] + RES2['YCT'][:, 99], label='RES2.YLT[:,99] + RES2.YCT[:,99]')
plt.plot(RES['YCT'][:, 99], label='RES.YCT[:,99]')
plt.plot(RES2['YCT'][:, 99], label='RES2.YCT[:,99]')

# Añade leyendas y muestra el gráfico
plt.legend()
plt.show()


# Preparamos datos de asimilación cross-shore
errY = 1  # Error en la posición inicial, [m]
errK = 50  # Error en la constante longshore
errvlt = 0.001 / 365 / 24  # Error en la tendencia
qY = 1e-5  # Ruido en la posición
qerrK = 1e-1  # Ruido en la constante longshore
qerrvlt = 1e-9  # Ruido en la tendencia a largo plazo
R = 1  # Error de la observación
# Preparamos asimilación cross-shore
errYct = 1
qerrYct = 0.2
errKacr = 1e-3
errKero = 1e-2
qerrKacr = 1e-4
qerrKero = 1e-2
errDy0 = 1
qerrDy0 = 0.2
R_c = 1
posguarda = list(range(2 * 365, 15 * 365 + 1, 15))  # Posiciones de guardado

for i in range(len(PERF)):
    # Longshore
    PERF[i]['rP0'] = np.diag([errY, errK, errvlt])
    PERF[i]['rQ'] = np.diag([qY, qerrK, qerrvlt])
    PERF[i]['R'] = R
    PERF[i]['date_obs'] = RES['t_output'][posguarda]
    PERF[i]['Y_obs_lt'] = RES['YLT'][posguarda, i]

    # Cross-shore
    PERF[i]['rPero0_c'] = np.diag([errYct, errKero, 0, 0, errDy0])
    PERF[i]['rPacr0_c'] = np.diag([0, 0, errYct, errKacr, errDy0])
    PERF[i]['rQero_c'] = np.diag([qerrYct, qerrKero, 0, 0, qerrDy0])
    PERF[i]['rQacr_c'] = np.diag([0, 0, qerrYct, qerrKacr, qerrDy0])
    PERF[i]['R_c'] = R_c
    PERF[i]['Y_obs_ct'] = RES['YCT'][posguarda, i]

# Crea una copia de INPUT2 y actualiza los parámetros necesarios
INPUT3 = INPUT2.copy()
INPUT3['plott'] = 0
INPUT3['PERF'] = PERF
INPUT3['data_asim_l'] = 1
INPUT3['data_asim_c'] = 1
INPUT3['OUTPUTLIST'] = ['Hbd', 'kcerc', 'vlt']  # Outputs por defecto: t_output e YLT

# Llama a la función IH_LANS con los nuevos parámetros
RES3 = IH_LANS(INPUT3)


tr = 30

# Gráfico Longshore
plt.figure(figsize=(10, 5))
plt.plot(RES['t_output'], RES['YLT'][:, tr], label='RES')
plt.plot(RES['t_output'], RES2['YLT'][:, tr], label='RES2')
plt.plot(RES['t_output'], RES3['YLT'][:, tr], label='RES3')
plt.plot(PERF[tr]['date_obs'], PERF[tr]['Y_obs_lt'], 'o', label='Observación')
plt.gca().xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))  # Formato de fechas en el eje x
plt.xlabel('Fecha')
plt.ylabel('YLT')
plt.legend()
plt.title('Gráfico Longshore')

# Gráfico Cross-shore
plt.figure(figsize=(10, 5))
plt.plot(RES['t_output'], RES['YCT'][:, tr], label='RES')
plt.plot(RES['t_output'], RES2['YCT'][:, tr], label='RES2')
plt.plot(RES['t_output'], RES3['YCT'][:, tr], label='RES3')
plt.plot(PERF[tr]['date_obs'], PERF[tr]['Y_obs_ct'], 'o', label='Observación')
plt.gca().xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))  # Formato de fechas en el eje x
plt.xlabel('Fecha')
plt.ylabel('YCT')
plt.legend()
plt.title('Gráfico Cross-shore')
plt.ylim([-50, 50])

plt.show()
