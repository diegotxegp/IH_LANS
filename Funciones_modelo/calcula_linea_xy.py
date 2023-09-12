def calcula_linea_xy(PERF, YLTi):
    xon = [profile["xon"] for profile in PERF]
    nx = [profile["nx"] for profile in PERF]

    yon = [profile["yon"] for profile in PERF]
    ny = [profile["ny"] for profile in PERF]

    x = [a+b*c for a,b,c in zip(xon, nx, YLTi)]
    y = [a+b*c for a,b,c in zip(yon, ny, YLTi)]
    
    return x, y