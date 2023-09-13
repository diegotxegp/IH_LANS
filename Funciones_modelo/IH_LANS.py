import numpy as np

from Funciones_modelo.detecta_lcs_bcpun import detecta_lcs_bcpun
from Funciones_modelo.casa_tobs_tcalc_lc import casa_tobs_tcalc_lc
from Funciones_modelo.une_perfiles_dinamicas import une_perfiles_dinamicas
from Funciones_modelo.perfiles_actuacion import perfiles_actuacion
from Funciones_modelo.calcula_z_cabeza_est import calcula_z_cabeza_est
from Funciones_modelo.propaga_cabeza import propaga_cabeza
from Funciones_modelo.propaga_rotura import propaga_rotura

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
        DYNP = propaga_rotura(PERF, DYN, gamma, t, refNMM, cotasZ, calcularotura)
    else:
        print('Calculando rotura con interacción')

    if inthidromorfo == 1 and calcularotura == 0:
        warning_msg = (
            "No se puede calcular la interacción hidro-morfo sin propagar a rotura\n"
            "Propagación a rotura activada"
        )
        print(warning_msg)
        calcula_rotura = 1


    print("Completed")
