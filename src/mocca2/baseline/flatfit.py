"""
Implementation of asymetric least squares with smoothness penalty
"""

from typing import Tuple
from numpy.typing import ArrayLike, NDArray

import warnings

import numpy as np
from scipy.signal import savgol_filter # type: ignore
from scipy import sparse # type: ignore
from scipy.sparse.linalg import spsolve # type: ignore


def flatfit(data: ArrayLike, smoothness: float, p: float) -> NDArray:
    """
    FlatFit: least squares weighted by inverse scale of 1st and 2nd derivatives with smoothness penalty

    Parameters
    ----------

    data: ArrayLike
        1D data

    smoothness: float
        size of smoothness penalty

    p: float
        relative size of Savitzky-Golay filter

    Returns
    -------
    NDArray
        values that minimize the asymmetric squared error with smoothness penalty

    Description
    -----------

    This routine finds vector `z` that minimized:

    `(y-z).T @ W @ (y-z) + smoothness * z.T @ D.T @ D @ z`

    where `W` is determined by slope and curvature at given point

    """
    y = np.array(data)
    assert len(y.shape) == 1, "Incorrect data input shape for AsLS"
    assert len(y) > 3, "At least 4 data points muts be provided for AsLS"

    L = len(y)
    
    # scale lambda to maintain invariance for sampling frequency
    lamb = smoothness * L**4

    # calculate w
    filter_window = max(4, int(L * p))
    slope = savgol_filter(y, filter_window, 3, deriv=1)
    curvature = np.gradient(slope)
    slope = slope**2
    slope /= np.sum(slope)
    curvature = curvature**2
    curvature /= np.sum(curvature)

    w = 1 / (slope + curvature + 1e-10)

    # calculate baseline

    D = sparse.diags([1, -2, 1], [0, -1, -2], shape=(L, L-2))
    H = lamb * D.dot(D.transpose())

    W = sparse.spdiags(w, 0, L, L)
    Z = W + H
    z = spsolve(Z, w*y)

    return z
