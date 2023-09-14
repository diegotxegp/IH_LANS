import numpy as np

from CircSat.circ_confmean import circ_confmean

def circ_mean(alpha, w=None, dim=0):
    """
    Compute the mean direction for circular data.

    Parameters:
    alpha: array_like
        Sample of angles in radians.
    w: array_like, optional
        Weightings in case of binned angle data.
    dim: int, optional
        Compute along this dimension, default is 0.

    Returns:
    mu: float
        Mean direction.
    ul: float, optional
        Upper 95% confidence limit.
    ll: float, optional
        Lower 95% confidence limit.
    """
    if dim is None:
        dim = 0

    if w is None:
        # If no specific weighting has been specified, assume no binning has taken place
        w = np.ones(alpha.shape)
    else:
        if w.shape != alpha.shape:
            raise ValueError('Input dimensions do not match')

    # Compute weighted sum of cos and sin of angles
    r = np.sum(w * np.exp(1j * alpha), axis=dim)

    # Obtain mean
    mu = np.angle(r)

    # Confidence limits if desired
    ul = None
    ll = None
    if ul is not None:
        t = circ_confmean(alpha, 0.05, w, None, dim)
        ul = mu + t
        ll = mu - t

    return mu, ul, ll