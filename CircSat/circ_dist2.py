import numpy as np

def circ_dist2(x, y=None):
    """
    Compute pairwise circular differences between angles.

    Args:
        x (array-like): A 1-D array of angles in radians.
        y (array-like, optional): A 1-D array of angles in radians. If not provided, `y` will be set to `x`.

    Returns:
        numpy.ndarray: A matrix with pairwise differences between angles.

    References:
        Biostatistical Analysis, J. H. Zar, p. 651

    PHB 3/19/2009

    Circular Statistics Toolbox for Matlab
    By Philipp Berens, 2009
    berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html
    """

    if y is None:
        y = x

    if x.ndim > 1:
        x = x.flatten()
    if y.ndim > 1:
        y = y.flatten()

    x = np.exp(1j * x)
    y = np.exp(1j * y)

    r = np.angle(np.outer(x, 1 / y))

    return r