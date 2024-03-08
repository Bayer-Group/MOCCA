"""Wrapper for baseline estimation routines"""

from typing import Literal
from numpy.typing import NDArray

import numpy as np

from scipy.signal import savgol_filter # type: ignore

from mocca2.baseline.arpls import arpls
from mocca2.baseline.asls import asls
from mocca2.baseline.flatfit import flatfit
from mocca2.classes import Data2D


def estimate_baseline(
        data: NDArray | Data2D,
        method : Literal['asls', 'arpls', 'flatfit'] = 'arpls',
        smoothness: float = 1.,
        p: float | None = None,
        tol: float = 1e-7,
        max_iter: int | None = None,
        smooth_wl: int | None = None
    ) -> NDArray:
    """
    Estimates baseline using AsLS, arPLS or FlatFit algorithm

    Parameters
    ----------
    data: NDArray | Data2D
        Data with shape [N] or [sample, N]

    method: Literal['asls', 'arpls', 'flatfit']
        Possible baseline estimation methods are AsLS, arPLS and FlatFit. FlatFit and AsLS work well with smooth data, asPLS works better with noisy data
        
    smoothness: float
        size of smoothness penalty

    p: float | None
        Assymetry factor, different for AsLS and arPLS

    tol: float
        maximum relative change of `w` for convergence

    max_iter: int | None 
        maximum number of iterations. If not specified, guessed automatically
    
    smooth_wl: int | None
        if specified, applies Savitzky-Golay filter (order 2) accross wavelength axis with given window size

    Returns
    -------
    NDArray
        values that minimize the asymmetric squared error with smoothness penalty,
        same shape as `data`

    See details in the individual routines or at
    [StackOverflow](https://stackoverflow.com/a/50160920) and 
    [10.1039/C4AN01061B](https://doi.org/10.1039/C4AN01061B).

    """

    if method not in {'asls', 'arpls', 'flatfit'}:
        raise AttributeError(f"Invalid method for baseline correction '{method}'")
    
    if isinstance(data, Data2D):
        y = data.data
    else:
        y = data

    assert len(y.shape) in {1,2}, "The data for baseline correction has invalid shape"
    
    # default value of p depends on the method
    param = p if p is not None else {'asls': 0.001, 'arpls': 2., 'flatfit': 0.03}[method]

    if method == 'asls':
        routine = lambda d, b0: asls(d, smoothness, param, tol, max_iter, b0)
    elif method == 'arpls':
        routine = lambda d, b0: arpls(d, smoothness, param, tol, max_iter, b0)
    elif method == 'flatfit':
        routine = lambda d, b0: flatfit(d, smoothness, param)
    

    if len(y.shape) == 1:
        return routine(y, None)

    baselines_arr = []
    b = None
    for d in y:
        b = routine(d, b)
        baselines_arr.append(b)

    baselines = np.array(baselines_arr)

    if smooth_wl is not None:
        baselines = savgol_filter(baselines, smooth_wl, 2, axis=0)

    return baselines