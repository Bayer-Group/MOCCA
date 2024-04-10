# taken from scipy, but added l2 norm for numerical stability

import numpy as np
from scipy.linalg import lstsq, LinAlgWarning
import warnings


def nnls(A, b, maxiter=None, tol=None, l2=0.0):
    """
    Non-negative least squares solver, adapted from SciPy's nnls
    """
    m, n = A.shape

    AtA = A.T @ A + np.eye(n) * l2
    Atb = b @ A  # Result is 1D - let NumPy figure it out

    if not maxiter:
        maxiter = 3 * n
    if tol is None:
        tol = 10 * max(m, n) * np.spacing(1.0)

    # Initialize vars
    x = np.zeros(n, dtype=np.float64)
    s = np.zeros(n, dtype=np.float64)
    # Inactive constraint switches
    P = np.zeros(n, dtype=bool)

    # Projected residual
    w = Atb.copy().astype(np.float64)  # x=0. Skip (-AtA @ x) term

    # Overall iteration counter
    # Outer loop is not counted, inner iter is counted across outer spins
    iter = 0

    while (not P.all()) and (w[~P] > tol).any():  # B
        # Get the "most" active coeff index and move to inactive set
        k = np.argmax(w * (~P))  # B.2
        P[k] = True  # B.3

        # Iteration solution
        s[:] = 0.0
        # B.4
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", message="Ill-conditioned matrix", category=LinAlgWarning
            )
            s[P] = lstsq(AtA[np.ix_(P, P)], Atb[P])[0]

        # Inner loop
        while (iter < maxiter) and (s[P].min() < 0):  # C.1
            iter += 1
            inds = P * (s < 0)
            alpha = (x[inds] / (x[inds] - s[inds])).min()  # C.2
            x *= 1 - alpha
            x += alpha * s
            P[x <= tol] = False
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore", message="Ill-conditioned matrix", category=LinAlgWarning
                )
                s[P] = lstsq(AtA[np.ix_(P, P)], Atb[P])[0]
            s[~P] = 0  # C.6

        x[:] = s[:]
        w[:] = Atb - AtA @ x

        if iter == maxiter:
            raise RuntimeError("NNLS: Iteration limit reached")

    return x, np.linalg.norm(A @ x - b)
