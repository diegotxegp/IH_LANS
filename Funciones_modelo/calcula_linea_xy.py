def calcula_linea_xy(PERF, YLTi):
    x = [perf['xon'] + perf['nx'] * YLTi for perf in PERF]
    y = [perf['yon'] + perf['ny'] * YLTi for perf in PERF]
    
    return x, y