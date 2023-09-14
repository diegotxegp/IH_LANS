import numpy as np

def circ_r(alpha, w=None, d=0):
    """
    Compute mean resultant vector length for circular data.

    Parameters:
    alpha: array_like
        Sample of angles in radians.
    w: array_like, optional
        Number of incidences in case of binned angle data.
    d: float, optional
        Spacing of bin centers for binned data. If supplied, a correction factor is used to correct for bias
        in estimation of r, in radians.

    Returns:
    r: float
        Mean resultant length.
    """
    # Check vector size
    if alpha.shape[1] > alpha.shape[0]:
        alpha = alpha.T

    if w is None:
        # If no specific weighting has been specified, assume no binning has taken place
        w = np.ones(alpha.shape)
    else:
        if w.shape[1] > w.shape[0]:
            w = w.T

    if d is None:
        # Per default, do not apply correction for binned data
        d = 0

    # Compute weighted sum of cos and sin of angles
    r = np.dot(w, np.exp(1j * alpha))

    # Obtain length
    r = abs(r) / np.sum(w)

    # For data with known spacing, apply correction factor to correct for bias
    # in the estimation of r (see Zar, p. 601, equ. 26.16)
    if d != 0:
        c = d / (2 * np.sin(d / 2))
        r = c * r

    return r