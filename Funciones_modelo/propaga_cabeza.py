import numpy as np

def propaga_cabeza(ACT, DYN, t):
    tdyncomun = np.where([d['t'] == t for d in DYN])[0]
    
    if not ACT:
        return ACT
    
    iper = [act['Tipo'] == 'dper' for act in ACT]
    ipar = [act['Tipo'] == 'dpar' for act in ACT]
    iparper = [act['Tipo'] == 'dperpar' for act in ACT]
    idif = [i for i, (p, a, ap) in enumerate(zip(iper, ipar, iparper)) if p or a or ap]
    xdyn = [d['X'] for d in DYN]
    ydyn = [d['Y'] for d in DYN]
    
    for i in idif:
        ACT0 = ACT[i]
        if ACT0['Tipo'] == 'dper':
            xest = ACT0['X'][-1]
            yest = ACT0['Y'][-1]
            d = [np.hypot(x - xest, y - yest) for x, y in zip(xdyn, ydyn)]
            posmin = np.argmin(d)
            h1 = abs(ACT0['Z']) + DYN[posmin]['SS'][tdyncomun] + DYN[posmin]['AT'][tdyncomun] + DYN[posmin]['SLR'][tdyncomun]
            Hs0 = DYN[posmin]['Hs'][tdyncomun]
            T0 = DYN[posmin]['Tp'][tdyncomun]
            Dir0 = DYN[posmin]['Dir'][tdyncomun]
            nivel = DYN[posmin]['SS'][tdyncomun] + DYN[posmin]['AT'][tdyncomun] + DYN[posmin]['SLR'][tdyncomun] + DYN[posmin]['h0']
            nbati = ACT0['nbati']
            Hb, alphab = snell_shoalref(Hs0, T0, Dir0, nivel, h1, nbati)
            ACT[i]['Hb'] = Hb
            ACT[i]['Dir'] = alphab
            ACT[i]['Tp'] = T0
        elif ACT0['Tipo'] == 'dperpar':
            Hb1 = np.zeros((len(tdyncomun), 2))
            alphab1 = np.zeros((len(tdyncomun), 2))
            Tb = np.zeros((len(tdyncomun), 2))
            for ip in range(2):
                ind = ip - 2
                xest = ACT0['X'][len(ACT0['X']) + ind]
                yest = ACT0['Y'][len(ACT0['Y']) + ind]
                d = [np.hypot(x - xest, y - yest) for x, y in zip(xdyn, ydyn)]
                posmin = np.argmin(d)
                h1 = abs(ACT0['Z'][len(ACT0['Z']) + ind]) + DYN[posmin]['SS'][tdyncomun] + DYN[posmin]['AT'][tdyncomun] + DYN[posmin]['SLR'][tdyncomun]
                Hs0 = DYN[posmin]['Hs'][tdyncomun]
                T0 = DYN[posmin]['Tp'][tdyncomun]
                Dir0 = DYN[posmin]['Dir'][tdyncomun]
                nivel = DYN[posmin]['SS'][tdyncomun] + DYN[posmin]['AT'][tdyncomun] + DYN[posmin]['SLR'][tdyncomun] + DYN[posmin]['h0']
                nbati = ACT0['nbati'][ip]
                Hb, alphab = snell_shoalref(Hs0, T0, Dir0, nivel, h1, nbati)
                Hb1[:, ip] = Hb
                alphab1[:, ip] = alphab
                Tb[:, ip] = T0
            ACT[i]['Hb'] = Hb1
            ACT[i]['Dir'] = alphab1
            ACT[i]['Tp'] = Tb
        elif ACT0['Tipo'] == 'dpar':
            Hb1 = np.zeros((len(tdyncomun), 2))
            alphab1 = np.zeros((len(tdyncomun), 2))
            Tb = np.zeros((len(tdyncomun), 2))
            for ip in range(2):
                xest = ACT0['X'][ip]
                yest = ACT0['Y'][ip]
                d = [np.hypot(x - xest, y - yest) for x, y in zip(xdyn, ydyn)]
                posmin = np.argmin(d)
                h1 = abs(ACT0['Z'][ip]) + DYN[posmin]['SS'][tdyncomun] + DYN[posmin]['AT'][tdyncomun] + DYN[posmin]['SLR'][tdyncomun]
                Hs0 = DYN[posmin]['Hs'][tdyncomun]
                T0 = DYN[posmin]['Tp'][tdyncomun]
                Dir0 = DYN[posmin]['Dir'][tdyncomun]
                nivel = DYN[posmin]['SS'][tdyncomun] + DYN[posmin]['AT'][tdyncomun] + DYN[posmin]['SLR'][tdyncomun] + DYN[posmin]['h0']
                nbati = ACT0['nbati'][ip]
                Hb, alphab = snell_shoalref(Hs0, T0, Dir0, nivel, h1, nbati)
                Hb1[:, ip] = Hb
                alphab1[:, ip] = alphab
                Tb[:, ip] = T0
            ACT[i]['Hb'] = Hb1
            ACT[i]['Dir'] = alphab1
            ACT[i]['Tp'] = Tb
    return ACT