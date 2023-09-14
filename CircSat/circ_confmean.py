import numpy as np
from scipy.stats import chi2

from CircSat.circ_r import circ_r

def circ_confmean(alpha, xi=0.05, w=None, d=0):
    """
    Compute the confidence limits on the mean for circular data.

    Parameters:
    alpha: array_like
        Sample of angles in radians.
    xi: float, optional
        (1-xi)-confidence limits are computed, default 0.05.
    w: array_like, optional
        Number of incidences in case of binned angle data.
    d: float, optional
        Spacing of bin centers for binned data. If supplied, a correction factor is used to correct for bias
        in estimation of r, in radians.

    Returns:
    t: float
        Mean +- d yields upper/lower (1-xi)% confidence limit.
    """
    # Check vector size
    if alpha.shape[1] > alpha.shape[0]:
        alpha = alpha.T

    # Set confidence limit size to default
    if xi is None:
        xi = 0.05

    if w is None:
        # If no specific weighting has been specified, assume no binning has taken place
        w = np.ones(alpha.shape)
    else:
        if w.shape[1] > w.shape[0]:
            w = w.T

        if len(alpha) != len(w):
            raise ValueError('Input dimensions do not match.')

    if d is None:
        # Per default, do not apply correction for binned data
        d = 0

    # Compute ingredients for conf. lim.
    r = circ_r(alpha, w, d)
    n = np.sum(w)
    R = n * r
    c2 = chi2.ppf(1 - xi, 1)

    # Check for resultant vector length and select appropriate formula
    if 0.9 > r > np.sqrt(c2 / (2 * n)):
        t = np.sqrt((2 * n * (2 * R ** 2 - n * c2)) / (4 * n - c2))  # Equ. 26.24
    elif r >= 0.9:
        t = np.sqrt(n ** 2 - (n ** 2 - R ** 2) * np.exp(c2 / n))  # Equ. 26.25
    else:
        t = np.nan
        print('Resultant vector does not allow specifying confidence limits on mean. '
              'Results may be wrong or inaccurate.')

    # Apply final transform
    t = np.arccos(t / R)
    
    return t