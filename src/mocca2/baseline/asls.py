"""
Implementation of asymetric least squares with smoothness penalty
"""

from typing import Tuple
from numpy.typing import ArrayLike, NDArray

import warnings

import numpy as np
from scipy import sparse # type: ignore
from scipy.sparse.linalg import spsolve # type: ignore


def asls(data: ArrayLike, smoothness: float, p: float, tol: float = 1e-7, max_iter: int | None = None, baseline_guess: NDArray | None = None) -> NDArray:
    """
    AsLS: Asymmetric Least Squares with smoothness penalty

    Parameters
    ----------

    data: ArrayLike
        1D data

    smoothness: float
        size of smoothness penalty

    p: float
        asymetry factor, `w = p if  y_fit < data else (1-p)`

    tol: float
        maximum relative change of `w` for convergence

    max_iter: int | None
        maximum number of iterations

    baseline_guess: ArrayLike | None
        initial guess for baseline

    Returns
    -------
    NDArray
        values that minimize the asymmetric squared error with smoothness penalty

    Description
    -----------

    This routine finds vector `z` that minimized:

    `(y-z).T @ W @ (y-z) + smoothness * z.T @ D.T @ D @ z`

    where `w = p if  y < z else (1-p)` and `D` is finite differences for second derivative.

    See details at
    [StackOverflow](https://stackoverflow.com/a/50160920) and 
    [10.1039/C4AN01061B](https://doi.org/10.1039/C4AN01061B).

    """
    y = np.array(data)
    assert len(y.shape) == 1, "Incorrect data input shape for AsLS"
    assert len(y) > 3, "At least 4 data points muts be provided for AsLS"

    L = len(y)

    # default number of maximum iterations
    if max_iter is None:
        max_iter = L
    
    # scale lambda to maintain invariance for sampling frequency
    lamb = smoothness * L**4 * 1e-7

    D = sparse.diags([1, -2, 1], [0, -1, -2], shape=(L, L-2))
    H = lamb * D.dot(D.transpose())

    def get_baseline(w):
        W = sparse.spdiags(w, 0, L, L)
        Z = W + H
        z = spsolve(Z, w*y)
        return z

    def get_w(z):
        w = p * (y > z) + (1-p) * (y < z)
        return w
    
    def get_loss(z, w):
        d = z - y
        loss = d.dot(w*d)
        loss += H.dot(z).dot(z)
        return loss

    # initialize values
    if baseline_guess is None:
        z = np.ones_like(y) * np.mean(y)
    else:
        z = baseline_guess

    w = get_w(z)
    old_w = w.copy()

    loss = np.inf

    for iters in range(max_iter):
        # determine search direction
        search_dir = get_baseline(old_w) - z

        # calculate loss
        s = search_dir
        for _ in range(5):
            new_z = z + s
            new_w = get_w(new_z)
            new_loss = get_loss(new_z, new_w)
            if new_loss < loss:
                z = new_z
                w = new_w
                loss = new_loss
            else:
                break
        
            s /= 2.

        change = np.linalg.norm(w-old_w)/np.linalg.norm(w)
        if change < tol:
            z = get_baseline(w)
            return z
        old_w = w

    warnings.warn(f'Failed to converge: after {max_iter} iterations the change is {change}, tol is {tol}')

    return z
