import os
import numpy as np
import datetime
import matplotlib.pyplot as plt

from Funciones_modelo.detecta_lcs_bcpun import detecta_lcs_bcpun
from Funciones_modelo.casa_tobs_tcalc_lc import casa_tobs_tcalc_lc
from Funciones_modelo.une_perfiles_dinamicas import une_perfiles_dinamicas
from Funciones_modelo.perfiles_actuacion import perfiles_actuacion
from Funciones_modelo.calcula_z_cabeza_est import calcula_z_cabeza_est
from Funciones_modelo.propaga_cabeza import propaga_cabeza
from Funciones_modelo.propaga_rotura import propaga_rotura
from Funciones_modelo.calcula_linea_xy import calcula_linea_xy
from Funciones_modelo.calcula_nbati_i import calcula_nbati_i
from Funciones_modelo.reduce_estructuras import reduce_estructuras
from Funciones_modelo.szonewidth import szonewidth
from Funciones_modelo.propaga_rotura_i import propaga_rotura_i
from Funciones_modelo.estructuras_clasifica import estructuras_clasifica
from Funciones_modelo.calc_difraction import calc_difraction
from Funciones_modelo.pintainstante import pintainstante
from Funciones_modelo.aplica_tasa import aplica_tasa

from Funciones_cshore.millerydean_dy0 import millerydean_dy0
from Funciones_cshore.calcula_cshore_md import calcula_cshore_md
from Funciones_cshore.kalman_transversal import kalman_transversal

from CircSat.circ_mean import circ_mean

def IH_LANS(INPUT):

    # Identificamos inputs/cargamos valores por defecto
    gamma = INPUT.get('gamma', 0.55)
    toutp = INPUT.get('toutp', 24 * 30)
    kcerc = INPUT.get('kcerc', [200] * len(INPUT['PERF']))
    fcourant = INPUT.get('fcourant', 0.5)
    kacr = INPUT.get('kacr', [3.89E-3] * len(INPUT['PERF']))
    kero = INPUT.get('kero', [1.74e-2] * len(INPUT['PERF']))
    tstab = INPUT.get('tstab', 365 * 2)
    tcent = INPUT.get('tcent', 365 * 5)
    vlt = INPUT.get('vlt', [0] * len(INPUT['PERF']))
    RSLR = INPUT.get('RSLR', [0] * len(INPUT['PERF']))
    tipotrans = INPUT.get('tipotrans', 'CERC')
    data_asim_l = INPUT.get('data_asim_l', 0)
    data_asim_c = INPUT.get('data_asim_c', 0)
    data_asim_lc = INPUT.get('data_asim_lc', 0)
    plott = INPUT.get('plott', 1)
    tplot = INPUT.get('tplot', 24)
    OUTPUTLIST = INPUT.get('OUTPUTLIST', [])
    inthidromorfo = INPUT.get('inthidromorfo', 0)
    refNMM = INPUT.get('refNMM', 0)
    cotasZ = INPUT.get('cotasZ', 'DEAN_A')
    calcularotura = INPUT.get('calcularotura', 1)


    # Detectamos perfiles
    PERF = INPUT['PERF']
    PLCS = [i for i, p in enumerate(PERF) if 'lcs' in p['Tipo']]  # lcs
    PLCS_CS = [i for i, p in enumerate(PERF) if 'cs' in p['Tipo']]  # lcs y cs
    PLCS_CS_LCS = [i for i in PLCS_CS if 'lcs' in PERF[i]['Tipo']]
    PCS = [i for i, p in enumerate(PERF) if 'cs' in p['Tipo'] and 'lcs' not in p['Tipo']]  # solo cs
    PS = [i for i, p in enumerate(PERF) if 's' in p['Tipo'] and 'cs' not in p['Tipo']]
    PAC = [i for i, p in enumerate(PERF) if 'ac' in p['Tipo']]

    # Iniciamos condiciones de contorno longshore
    if PLCS:
        PBC, PNL = detecta_lcs_bcpun(INPUT['PERF'])
        ntramos = len(PBC)
        if 'bctype' not in INPUT:
            bctype = [['Dirichlet', 'Dirichlet']] * ntramos
        else:
            if len(INPUT['bctype']) != ntramos:
                print('Las filas de bctype no coinciden con el número de tramos')
                print('Asignamos el resto de tramos igual que la primera fila')
                bctype = [INPUT['bctype'][0]] * ntramos
            else:
                bctype = INPUT['bctype']

        if 'bctypeval' not in INPUT:
            bctypeval = [[0, 0]] * ntramos
        else:
            if len(INPUT['bctypeval']) != ntramos:
                print('Las filas de bctypeval no coinciden con el número de tramos')
                print('Asignamos el resto de tramos igual que la primera fila')
                bctypeval = [INPUT['bctypeval'][0]] * ntramos
            else:
                bctypeval = INPUT['bctypeval']

        if len(kcerc) < len(PERF):
            print('Longitud de kcerc no coincide con el número de perfiles')
            print('Completamos con el primer valor')
            kcerc = [kcerc[0]] * len(PERF)

        if len(vlt) < len(PERF):
            print('Longitud de vlt no coincide con el número de perfiles')
            print('Completamos con el primer valor')
            vlt = [vlt[0]] * len(PERF)

    if len(RSLR) < len(PERF):
        print('Longitud de RSLR no coincide con el número de perfiles')
        print('Completamos con el primer valor')
        RSLR = [RSLR[0]] * len(PERF)

    if len(kacr) < len(PERF):
        print('Longitud de kcerc no coincide con el número de perfiles')
        print('Completamos con el primer valor')
        kacr = [kacr[0]] * len(PERF)

    if len(kero) < len(PERF):
        print('Longitud de kcerc no coincide con el número de perfiles')
        print('Completamos con el primer valor')
        kero = [kero[0]] * len(PERF)

    if data_asim_l and data_asim_lc:
        print('data_asim_c y data_asim_lc activados, se aplica data_asim_lc')

    if data_asim_c and data_asim_lc:
        print('data_asim_l y data_asim_lc activados, se aplica data_asim_lc')

    # Cargamos inputs
    PERF = INPUT['PERF']
    DYN = INPUT['DYN']
    ACT = INPUT['ACT']
    t = INPUT['t']
    Dc = [p['dc'] for p in PERF]
    Ber = [p['Berma'] for p in PERF]
    dt = INPUT['dt']

    # verificamos si hay datos para asimilar
    if data_asim_l or data_asim_c or data_asim_lc:
        if data_asim_l and not data_asim_c and hasattr(PERF[0], 'rP0') and hasattr(PERF[0], 'rQ') and hasattr(PERF[0], 'R') and hasattr(PERF[0], 'date_obs') and hasattr(PERF[0], 'Y_obs_lt'):
            print('Preparando asimilación longshore')
            sigmaK = INPUT.get('sigmaK', 0.5)
            DA = casa_tobs_tcalc_lc(t, PERF, data_asim_l, data_asim_c, data_asim_lc, sigmaK=sigmaK, kcerc=kcerc)
            # Borramos información de PERF transferida a DA
            for key in ['date_obs', 'Y_obs_lt', 'rP0', 'rQ', 'R']:
                PERF[0].pop(key, None)
            # DALCS = DA[PLCS]
        elif data_asim_c and not data_asim_l and hasattr(PERF[0], 'rPero0_c') and hasattr(PERF[0], 'rQero_c') and hasattr(PERF[0], 'R_c') and hasattr(PERF[0], 'date_obs') and hasattr(PERF[0], 'Y_obs_ct'):
            print('Preparando asimilación cross-shore')
            DA = casa_tobs_tcalc_lc(t, PERF, data_asim_l, data_asim_c, data_asim_lc, kacr=kacr, kero=kero)
            for key in ['date_obs', 'Y_obs_ct', 'rPero0_c', 'rPacr0_c', 'rQero_c', 'rQacr_c', 'R_c']:
                PERF[0].pop(key, None)
            # DACS = DA[PLCS_CS]
        elif (data_asim_lc or (data_asim_c and data_asim_l)) and hasattr(PERF[0], 'rPero0_c') and hasattr(PERF[0], 'rQero_c') and hasattr(PERF[0], 'rPacr0_c') and hasattr(PERF[0], 'rQacr_c') and hasattr(PERF[0], 'R_c') and hasattr(PERF[0], 'date_obs') and hasattr(PERF[0], 'Y_obs_ct') and \
                hasattr(PERF[0], 'rP0') and hasattr(PERF[0], 'rQ') and hasattr(PERF[0], 'R') and hasattr(PERF[0], 'Y_obs_lt'):
            print('Preparando asimilación longshore y cross-shore')
            sigmaK = INPUT.get('sigmaK', 0.5)
            DA = casa_tobs_tcalc_lc(t, PERF, data_asim_l, data_asim_c, data_asim_lc, sigmaK=sigmaK, kcerc=kcerc, kacr=kacr, kero=kero)
            for key in ['date_obs', 'Y_obs_lt', 'rP0', 'rQ', 'R', 'Y_obs_ct', 'rPero0_c', 'rQero_c', 'rPacr0_c', 'rQacr_c', 'R_c']:
                PERF[0].pop(key, None)
            # DALCS = DA[PLCS]
            # DACS = DA[PLCS_CS]
        else:
            raise ValueError('No se puede asimilar, inputs incompletos')
        

    # Calculamos el espaciado entre perfiles
    x = np.array([perf['xon'] for perf in PERF]) + np.array([perf['nx'] for perf in PERF]) * np.array([perf['yc'] for perf in PERF])
    y = np.array([perf['yon'] for perf in PERF]) + np.array([perf['ny'] for perf in PERF]) * np.array([perf['yc'] for perf in PERF])
    dx = np.empty(len(PERF))
    dx[1:-1] = 0.5 * (np.sqrt((x[2:] - x[1:-1])**2 + (y[2:] - y[1:-1])**2)
        + np.sqrt((x[1:-1] - x[:-2])**2 + (y[1:-1] - y[:-2])**2))
    dx[0] = np.sqrt((x[1] - x[0])**2 + (y[1] - y[0])**2)
    dx[-1] = np.sqrt((x[-1] - x[-2])**2 + (y[-1] - y[-2])**2)

    # Inicializamos resultados
    len_t = len(t)
    Ntr = len(PERF)
    len_interval = len(np.where(np.mod(np.arange(1, len_t), toutp) == 0)[0]) + 1
    RES = {
        'YLT': np.full((len_interval, Ntr), np.nan),
        't_output': np.full((len_interval, 1), np.nan)
    }
    RES['YLT'] = np.full((len_interval, Ntr), np.nan)
    YLT = [perf['yc'] for perf in PERF]
    RES['t_output'] = np.full((len_interval, 1), np.nan)
    t_output = t[0]

    if data_asim_l and not data_asim_c:
        if 'OUTPUTPERF_ASIM' in INPUT:
            psaveasim = INPUT['OUTPUTPERF_ASIM']
            Ntra = len(psaveasim)
        else:
            Ntra = Ntr
            psaveasim = list(range(1, Ntr + 1))

        output_list = ['Hbd', 'Dbd', 'wbd', 'Q', 'kest', 'dQdx', 'Qbc', 'saltoYlt', 'kcerc', 'vlt', 'rP']
        output_length = [Ntr, Ntr, Ntr, Ntr - 1, Ntr - 1, Ntr, Ntr + 1, Ntra, Ntra, Ntra, Ntra]

        # Inicializamos el salto
        saltoYlt = np.zeros(Ntra)

        if 'OUTPUTLIST' in INPUT:
            output_var = INPUT['OUTPUTLIST']
            # Extraemos las variables comunes a output_list
            posvar = [output_list.index(var) for var in output_var if var in output_list]
            # Inicializamos los resultados para las variables solicitadas
            for iou in posvar:
                if output_list[iou] == 'rP':
                    RES['rP'] = np.full((len_interval, 3, output_length[iou]), np.nan)
                else:
                    RES[output_list[iou]] = np.full((len_interval, output_length[iou]), np.nan)

    elif data_asim_c and not data_asim_l:
        if 'OUTPUTPERF_ASIM' in INPUT:
            psaveasim = INPUT['OUTPUTPERF_ASIM']
            Ntra = len(psaveasim)
        else:
            Ntra = Ntr
            psaveasim = list(range(1, Ntr + 1))

        output_list = ['Hbd', 'Dbd', 'wbd', 'Q', 'kest', 'dQdx', 'Qbc', 'saltoYct', 'kero', 'kacr', 'dy0']
        output_length = [Ntr, Ntr, Ntr, Ntr - 1, Ntr - 1, Ntr, Ntr + 1, Ntra, Ntra, Ntra, Ntra]

        # Inicializamos el salto
        saltoYct = np.zeros(Ntra)

        if 'OUTPUTLIST' in INPUT:
            output_var = INPUT['OUTPUTLIST']
            # Extraemos las variables comunes a output_list
            posvar = [output_list.index(var) for var in output_var if var in output_list]
            # Inicializamos los resultados para las variables solicitadas
            for iou in posvar:
                if output_list[iou] == 'rP':
                    RES['rP'] = np.full((len_interval, 3, output_length[iou]), np.nan)
                else:
                    RES[output_list[iou]] = np.full((len_interval, output_length[iou]), np.nan)

    elif (data_asim_l and data_asim_c) or data_asim_lc:
        if 'OUTPUTPERF_ASIM' in INPUT:
            psaveasim = INPUT['OUTPUTPERF_ASIM']
            Ntra = len(psaveasim)
        else:
            Ntra = Ntr
            psaveasim = list(range(1, Ntr + 1))

        output_list = ['Hbd', 'Dbd', 'wbd', 'Q', 'kest', 'dQdx', 'Qbc', 'saltoYlt', 'kcerc', 'vlt', 'rP', 'saltoYct', 'kero', 'kacr', 'dy0']
        output_length = [Ntr, Ntr, Ntr, Ntr - 1, Ntr - 1, Ntr, Ntr + 1, Ntra, Ntra, Ntra, Ntra, Ntra, Ntra, Ntra, Ntra]

        # Inicializamos el salto
        saltoYlt = np.zeros(Ntra)
        saltoYct = np.zeros(Ntra)

        if 'OUTPUTLIST' in INPUT:
            output_var = INPUT['OUTPUTLIST']
            # Extraemos las variables comunes a output_list
            posvar = [output_list.index(var) for var in output_var if var in output_list]
            # Inicializamos los resultados para las variables solicitadas
            for iou in posvar:
                RES[output_list[iou]] = np.full((len_interval, output_length[iou]), np.nan)

    else:
        output_list = ['Hbd', 'Dbd', 'wbd', 'Q', 'kest', 'dQdx', 'Qbc']
        output_length = [Ntr, Ntr, Ntr, Ntr - 1, Ntr - 1, Ntr, Ntr + 1]

        if 'OUTPUTLIST' in INPUT:
            output_var = INPUT['OUTPUTLIST']
            # Extraemos las variables comunes a output_list
            posvar = [output_list.index(var) for var in output_var if var in output_list]
            # Inicializamos los resultados para las variables solicitadas
            for iou in posvar:
                RES[output_list[iou]] = np.full((len_interval, output_length[iou]), np.nan)

    count_output = 1
    cdif = ['dper', 'dpar', 'dperpar']  # Estructuras que difractan
            

    # Propagamos a rotura en estructuras
    dinperf = une_perfiles_dinamicas(PERF, DYN)
    tdyncomun = np.where(np.isin(DYN[0]["t"], t))[0].tolist()
    ACT = perfiles_actuacion(ACT, PERF)
    ACT = calcula_z_cabeza_est(ACT, PERF, cotasZ)
    ACT = propaga_cabeza(ACT, DYN, tdyncomun)

   # Propagamos a rotura sin difracción si no hay interacción
    if inthidromorfo == 0:
        print('Calculando rotura sin interacción')
        DYNP, poscalc = propaga_rotura(PERF, DYN, gamma, t, refNMM, cotasZ, calcularotura)
    else:
        print('Calculando rotura con interacción')

    if inthidromorfo == 1 and calcularotura == 0:
        warning_msg = (
            "No se puede calcular la interacción hidro-morfo sin propagar a rotura\n"
            "Propagación a rotura activada"
        )
        print(warning_msg)
        calcula_rotura = 1


    ## Inicializamos Cshore calculado posicion equilibrio
    # Usamos dinámicas sin propagar a rotura
    # Buscamos perfiles a calcular
    adean = [perf["Adean"][0] for perf in PERF]
    berma = [perf["Berma"] for perf in PERF]
    nbati = [perf["nbati"] for perf in PERF]
    _, dy00, _ = millerydean_dy0(tstab, tcent, adean, kacr, kero, gamma, DYN, dinperf, t, berma, dt, PLCS_CS, calcularotura, nbati, refNMM, gamma)
    dy0 = -dy00
    # dy0 = np.zeros_like(dy0)  # Comentado porque no se usa en el código
    # yst0 = millerydean_na(-dy0, [PERF.Adean], kacr, kero, gamma, DYN, dinperf, t, [PERF.Berma], dt, PCS, calcularotura, [PERF.nbati], refNMM, gamma)


    # Cuerpo del modelo
    YLTi = [perf["yc"] for perf in PERF]
    YCTi = np.zeros(len(YLTi))
    DY0 = np.zeros(len(YLTi))

    if len(PLCS) > 0:
        if plott == 1:
            h = plt.figure()
            h.set_size_inches(18.12, 6.66)

        nbati0 = [perf["nbati"] for perf in PERF]

        for it in range(len(t) - 1):
            if inthidromorfo == 1:  # con interacción hidromorfo
                if 'alpha_int' not in INPUT:
                    alpha_int = 1
                else:
                    alpha_int = INPUT.alpha_int

                xlc, ylc = calcula_linea_xy(PERF, YLTi)
                nbatii = calcula_nbati_i(xlc, ylc, 'PBC', PBC)
                nbati_calc = circ_mean(
                    [nbatii, nbati0] * np.pi / 180, np.tile([alpha_int, 1 - alpha_int], (1, len(nbatii))), 1) * 180 / np.pi
                poscalc = np.argmin(np.abs(t[it] - DYN[0]["t"]))
                H0 = reduce_estructuras(DYN, 'Hs', poscalc)
                D0 = reduce_estructuras(DYN, 'Dir', poscalc)
                Tp0 = reduce_estructuras(DYN, 'Tp', poscalc)
                AT0 = reduce_estructuras(DYN, 'AT', poscalc)
                SS0 = reduce_estructuras(DYN, 'SS', poscalc)
                SLR0 = reduce_estructuras(DYN, 'SLR', poscalc)
                h0 = [DYN.h0]
                Hi, D0, Di, w0 = propaga_rotura_i(H0, D0, Tp0, AT0, SS0, SLR0, refNMM, h0, dinperf, nbati_calc, cotasZ, gamma, PERF)
            else:
                Hi = reduce_estructuras(DYNP, 'Hb', it)
                D0 = reduce_estructuras(DYNP, 'Dir0', it)
                Di = reduce_estructuras(DYNP, 'Dirb', it)
                w0 = reduce_estructuras(DYNP, 'wb', it)
                AT0 = reduce_estructuras(DYNP, 'AT', it)
                SS0 = reduce_estructuras(DYNP, 'SS', it)
                SLR0 = reduce_estructuras(DYNP, 'SLR', it)

            # Detectamos estructuras activas
            if len(ACT) > 0:
                EA = estructuras_clasifica(t[it], ACT)

                camposdif = [i for i, x in enumerate(EA["columns"]) if x in cdif]

                if camposdif:  # existen estructuras que difractan
                    estdif = []

                    for tt in camposdif:
                        estdif.extend(EA[cdif[tt]])

                    Hbd, Dbd = calc_difraction(PERF, ACT, estdif, Hi, Di, w0, YLTi, it, gamma)
                    wbd = szonewidth(Hbd / gamma, cotasZ, PERF)
                else:
                    wbd = w0
                    Hbd = Hi
                    Dbd = Di
            else:  # no hay estructuras en ningún momento
                EA = []
                wbd = w0
                Hbd = Hi
                Dbd = Di

            ATi = np.zeros(len(PERF))
            SSi = np.zeros(len(PERF))

            # Ojo, añadir resto dinámicas
            # Añadimos Bruun a YLT aquí

            escalaprin = 50

            if it % tplot == 0 and plott == 1:
                pintainstante(PERF, YLTi + YCTi, ACT, EA, Hbd, D0, Dbd, wbd, t, it, escalaprin)
                drawnow
                cla()

            # Calculamos cshore M&D
            Yct, posero, Yeq = calcula_cshore_md(YCTi, wbd, Hbd, SSi, ATi, kacr, kero, dy0, PLCS_CS, dt, [PERF.Berma])

            # Actualizamos
            YLTi = corrige_instante_inicial_escollera(YLTi, ACT, EA, t, it)

            # Calculamos caudales longitudinales
            Q, kest = calcula_caudal(Hbd, Dbd, wbd, dx, EA, YLTi, ACT, PERF, it, tipotrans)
            Qc = reparto_bypass(Q, EA, ACT, PERF, YLTi, wbd, it)

            # Calcula gradiente Q con BC en extremos
            dQdx, Qbc = calcula_gradientes(Qc, dx, PNL, PBC, bctype, bctypeval)

            # Verifica Courant
            tmin = verifica_courant(dx, Qbc, [PERF.dc] + [PERF.Berma], kcerc)
            if tmin < dt:
                npartes = min(100, int(1 / tmin / fcourant))
                warnings.warn(f'Se viola el Courant, Tmin = {tmin:.2f}\nSubdividimos en {npartes} el step')
                for ipartes in range(npartes):
                    Yltii = YLTi
                    Q, _ = calcula_caudal(Hbd, Dbd, wbd, dx, EA, Yltii, ACT, PERF, it, tipotrans)
                    Qc = reparto_bypass(Q, EA, ACT, PERF, Yltii, wbd, it)
                    dQdx, Qbc = calcula_gradientes(Qc, dx, PNL, PBC, bctype, bctypeval)
                    Yltii, _ = calcula_lc(Yltii, dt, Dc + Ber, kcerc / npartes, dQdx, Qbc, dx, ACT, EA)
                    YLTi = Yltii
            else:
                # No se viola el Courant
                Yltii, kcerc = calcula_lc(YLTi, dt, Dc + Ber, kcerc, dQdx, Qbc, dx, ACT, EA)
                YLTi = Yltii

            # Aplico tasas y tendencias
            Ylt = aplica_tasa(Yltii, ACT, EA, vlt, RSLR)

            # Asimilación de observaciones
            if data_asim_l and not data_asim_lc:
                Ylt, kcerc, vlt, saltoYlt, DALCSi = kalman_longitudinal(Ylt, kcerc, vlt, dQdx, dt, DA(PLCS), it, PLCS, Dc, Ber, sigmaK)
                DA(PLCS) = DALCSi

            if data_asim_c and not data_asim_lc:
                Yct, kacr, kero, saltoYct, DACSi, dy0 = kalman_transversal(Yct, YCTi, Yeq, kacr, kero, dt, DA(PLCS_CS), it, PLCS_CS, posero, dy0)
                DA(PLCS_CS) = DACSi
                DY0(PLCS_CS) = dy0

            if data_asim_lc:
                Ylt1, kcerc1, vlt1, saltoYlt1, DALCS1, Yct1, kacr1, kero1, saltoYct1, dy01 = kalman_longitudinal_transversal(
                    Ylt, kcerc, vlt, dQdx, dt, DA(PLCS), it, PLCS, PLCS_CS_LCS, Dc, Ber, sigmaK, Yct, YCTi, Yeq, kacr, kero, posero, dy0)
                Ylt = Ylt1
                kcerc = kcerc1
                vlt = vlt1
                saltoYlt = saltoYlt1
                DA(PLCS) = DALCS1
                kacr = kacr1
                kero = kero1
                saltoYct = saltoYct1
                dy0 = dy01


            # Guardamos resultados
            if it % toutp == 0:
                print(f"{it / (len(t) - 1) * 100:.2f}% completado")
                count_output += 1
                RES["YLT"][count_output, :] = Ylt
                RES["YCT"][count_output, :] = Yct
                RES["t_output"][count_output, :] = t[it + 1]

                if "posvar" in locals():
                    for ivar in range(len(posvar)):
                        if output_list[posvar[ivar]] == "rP":
                            rP = np.empty((3, len(psaveasim)))
                            for is0 in range(len(psaveasim)):
                                tsave = psaveasim[is0]
                                temp = DALCS[tsave].P0
                                erYlt = np.sqrt(temp[0, 0])
                                erK1 = np.sqrt(temp[1, 1])
                                K1 = np.log(kcerc[tsave] / DALCS[tsave].kcerc0) / DALCS[tsave].sigmaK
                                erK = -kcerc[tsave] + DALCS[tsave].kcerc0 * np.exp(DALCS[tsave].sigmaK * (K1 + erK1))
                                ervlt = np.sqrt(temp[2, 2])
                                rP[:, is0] = [erYlt, erK, ervlt]
                            RES["rP"][count_output, :, :] = rP

                        elif output_list[posvar[ivar]] == "saltoYlt":
                            RES[output_list[posvar[ivar]]][count_output, 0:output_length[posvar[ivar]]] = saltoYlt[psaveasim]
                        elif output_list[posvar[ivar]] == "kcerc":
                            RES[output_list[posvar[ivar]]][count_output, 0:output_length[posvar[ivar]]] = kcerc[psaveasim]
                        elif output_list[posvar[ivar]] == "vlt":
                            RES[output_list[posvar[ivar]]][count_output, 0:output_length[posvar[ivar]]] = vlt[psaveasim]
                        elif output_list[posvar[ivar]] == "saltoYct":
                            RES[output_list[posvar[ivar]]][count_output, 0:output_length[posvar[ivar]]] = saltoYct[psaveasim]
                        elif output_list[posvar[ivar]] == "kacr":
                            RES[output_list[posvar[ivar]]][count_output, 0:output_length[posvar[ivar]]] = kacr[psaveasim]
                        elif output_list[posvar[ivar]] == "kero":
                            RES[output_list[posvar[ivar]]][count_output, 0:output_length[posvar[ivar]]] = kero[psaveasim]
                        elif output_list[posvar[ivar]] == "dy0":
                            RES[output_list[posvar[ivar]]][count_output, 0:output_length[posvar[ivar]]] = DY0[psaveasim]
                        else:
                            RES[output_list[posvar[ivar]]][count_output, 0:output_length[posvar[ivar]]] = eval(output_list[posvar[ivar]])

            # Actualizamos YLTi para el siguiente time step
            YLTi = Ylt
            YCTi = Yct

    if not lcs:  # no lcs
        if plott == 1:
            h = plt.figure()
            h.set_position([50, 50, 1812, 666])

        nbati0 = [PERF.nbati]
        for it in range(len(t) - 1):
            if inthidromorfo == 1:  # con interacción hidromorfo
                if "alpha_int" not in INPUT:
                    alpha_int = 1
                else:
                    alpha_int = INPUT["alpha_int"]
                xlc, ylc = calcula_linea_xy(PERF, YLTi)
                nbatii = calcula_nbati_i(xlc, ylc, 'PBC', PBC)
                nbati_calc = circ_mean([nbatii, nbati0] * np.pi / 180, np.tile([alpha_int, 1 - alpha_int], (1, len(nbatii), 1)), 1) * 180 / np.pi
                poscalc = np.argmin(np.abs(t[it] - DYN[0]["t"]))
                H0 = reduce_estructuras(DYN, 'Hs', poscalc)
                D0 = reduce_estructuras(DYN, 'Dir', poscalc)
                Tp0 = reduce_estructuras(DYN, 'Tp', poscalc)
                AT0 = reduce_estructuras(DYN, 'AT', poscalc)
                SS0 = reduce_estructuras(DYN, 'SS', poscalc)
                SLR0 = reduce_estructuras(DYN, 'SLR', poscalc)
                RSLR = reduce_estructuras(DYN, 'RSLR', it)
                h0 = [DYN[0]["h0"]]
                Hi, D0, Di, w0 = propaga_rotura_i(H0, D0, Tp0, AT0, SS0, SLR0, refNMM, h0, dinperf, nbati_calc, cotasZ, gamma, PERF)

            else:
                Hi = reduce_estructuras(DYNP, 'Hb', it)
                D0 = reduce_estructuras(DYNP, 'Dir0', it)
                Di = reduce_estructuras(DYNP, 'Dirb', it)
                w0 = reduce_estructuras(DYNP, 'wb', it)
                AT0 = reduce_estructuras(DYNP, 'AT', it)
                SS0 = reduce_estructuras(DYNP, 'SS', it)
                SLR0 = reduce_estructuras(DYNP, 'SLR', it)
                RSLR = reduce_estructuras(DYNP, 'RSLR', it)

            # Detectamos estructuras activas
            if not isinstance(ACT, list) or len(ACT) == 0:
                EA = []
            else:
                EA = estructuras_clasifica(t[it], ACT)

            camposdif = [campo for campo in EA if campo in cdif]

            if len(camposdif) > 0:
                estdif = []
                for tt in range(len(camposdif)):
                    estdif.extend(EA[camposdif[tt]])
                Hbd, Dbd = calc_difraction(PERF, ACT, estdif, Hi, Di, w0, YLTi, it, gamma)
                wbd = szonewidth(Hbd / gamma, cotasZ, PERF)
            else:
                wbd = w0
                Hbd = Hi
                Dbd = Di

            ATi = np.zeros(len(PERF))
            SSi = np.zeros(len(PERF))
            # Ojo añadir resto dinámicas
            # Añadimos Bruun a YLT aquí
            escalaprin = 50
            if it % tplot == 0 and plott == 1:
                pintainstante(PERF, YLTi + YCTi, ACT, EA, Hbd, D0, Dbd, wbd, t, it, escalaprin)
                plt.draw()
                plt.clf()

            # Calculamos cshore M&D
            Yct, posero, Yeq = calcula_cshore_md(YCTi, wbd, Hbd, SSi, ATi, kacr, kero, dy0, PLCS_CS, dt, [PERF["Berma"]])

            # Aplicamos tasas y tendencias
            Ylt = aplica_tasa(YLTi, ACT, EA, vlt, RSLR)

            # Asimilación de observaciones
            if data_asim_c and not data_asim_lc:
                Yct, kacr, kero, saltoYct, DACSi, dy0 = kalman_transversal(Yct, YCTi, Yeq, kacr, kero, dt, DA[PLCS_CS], it, PLCS_CS, posero, dy0)
                DA[PLCS_CS] = DACSi
                DY0[PLCS_CS] = dy0

    # Cierre de todas las figuras (plots)
    if "path_save" in INPUT:
        if "nsave" in INPUT:
            nsave = INPUT["nsave"]
        else:
            nsave = ""
        if not os.path.exists(INPUT["path_save"]):
            os.mkdir(INPUT["path_save"])
        namesave = os.path.join(INPUT["path_save"], f"TEST_{nsave}_date_{datetime.now().strftime('%H-%M-%S')}.mat")
        # Borramos DYNP del guardado, pesa mucho
        if "DYNP" in INPUT:
            del INPUT["DYNP"]
        np.savez(namesave, RES=RES, INPUT=INPUT)