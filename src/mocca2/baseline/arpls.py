"""
Implementation of asymmetrically reweighted penalized least squares with smoothness penalty

Taken from https://doi.org/10.1039/C4AN01061B
"""

from typing import Tuple
from numpy.typing import ArrayLike, NDArray

import warnings

import numpy as np
from scipy import sparse  # type: ignore
from scipy.sparse.linalg import spsolve  # type: ignore


def arpls(data: ArrayLike, smoothness: float, p: float = 2., tol: float = 1e-7, max_iter: int | None = None, baseline_guess: NDArray | None = None) -> NDArray:
    """
    arPLS: Asymmetrically Reweighted Penalized Least Squares

    Parameters
    ----------

    data: ArrayLike
        1D data

    smoothness: float
        size of smoothness penalty
    
    p: float
        lower values shift the baseline lower

    tol: float
        maximum relative change of `w` for convergence

    max_iter: int | None
        maximum number of iterations

    baseline_guess: ArrayLike | None
        initial guess for baseline

    Returns
    -------
    NDArray
        values that minimize the asymmetric squared error with smoothness penalty, and `w`

    Description
    -----------

    This routine finds vector `z` that minimized:

    `(y-z).T @ W @ (y-z) + smoothness * z.T @ D.T @ D @ z`

    where `w` is nonlinear weighting function and `D` is finite differences for second derivative.

    See details at
    [10.1039/C4AN01061B](https://doi.org/10.1039/C4AN01061B).

    """

    y = np.array(data)
    assert len(y.shape) == 1, "Incorrect data input shape for arPLS"
    assert len(y) > 3, "At least 4 data points muts be provided for arPLS"

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
        z = spsolve(W + H, w*y)
        return z

    def get_w(z):
        d = y - z
        d_neg = d[d < 0]
        if len(d_neg) < 2:
            return np.ones_like(d)
        m = np.mean(d_neg)
        s = np.std(d_neg)
        w = sigmoid(2*(d - (-m+p*s)) / s)
        return w
    
    def get_loss(z, w):
        d = z-y
        loss = d.dot(w*d)
        loss += H.dot(z).dot(z)
        return loss

    # initial guess for baseline
    if baseline_guess is not None:
        z = baseline_guess
    else:
        z = np.ones(L) * np.mean(y)
    old_w = get_w(z)
    w = old_w

    loss = np.inf
    old_search_dir = None
    step = None

    for _ in range(max_iter):
        # determine search direction
        search_dir = get_baseline(old_w) - z

        # adjust search direction based on momentum
        # the momentum is same as in conjugate gradient, but search_dir is used instead of gradient
        if old_search_dir is not None:
            beta = search_dir.dot(search_dir) / old_search_dir.dot(old_search_dir)
            step = search_dir + beta * step
        else:
            step = search_dir

        old_search_dir = search_dir

        # calculate loss
        s = step
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


def sigmoid(x: NDArray) -> NDArray:
    """Logistic function with prevented overflow"""
    return 1/(1+np.exp(np.clip(x, None, 100)))
