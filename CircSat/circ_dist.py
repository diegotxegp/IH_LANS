import numpy as np

def circ_dist(x, y):
    """
    Compute pairwise differences x_i - y_i around the circle efficiently.

    Parameters:
    - x: array-like
        Sample of linear random variable.
    - y: array-like
        Sample of linear random variable or one single angle.

    Returns:
    - r: ndarray
        Matrix with differences.

    References:
    - Biostatistical Analysis, J. H. Zar, p. 651

    Circular Statistics Toolbox for Matlab
    By Philipp Berens, 2009
    berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html
    """
    if x.shape[1] > x.shape[0]:
        x = x.T

    if y.shape[1] > y.shape[0]:
        y = y.T

    if len(x) != len(y) and not len(y) == 1:
        raise ValueError('Input dimensions do not match.')

    r = np.angle(np.exp(1j * x) / np.exp(1j * y))

    return r