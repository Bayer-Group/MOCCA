from typing import Tuple
from numpy.typing import NDArray

import warnings

import numpy as np

from mocca2.classes import Data2D
from mocca2.deconvolution.nonnegative_lstsq import concentrations_from_spectra, spectra_from_concentrations
from mocca2.deconvolution.guess_spectra import guess_spectra

def alternating_lstsq(
        data: NDArray | Data2D,
        n_compounds: int | None = None,
        initial_spectra: NDArray | None = None,
        initial_concentrations: NDArray | None = None,
        max_iter: int = 100,
        rtol : float = 1e-7
    ) -> Tuple[NDArray, NDArray, float]:
    """
    Finds non-negative spectra and concentrations that minimize the mean square error.

    One of `n_compounds`, `initial_spectra`, `initial concentrations` must be provided.
    
    Iteratively calls `concentrations_from_spectra` and `spectra_from_concentrations` until convergence.

    Parameters
    ----------

    data: NDArray | Data2D
        The peak data that should be decomposed

    n_compounds: int | None
        Number of components that should be used for deconvolution (rank of data).
        If initial guess is not provided, uses `guess_spectra` to get initial spectra
    
    initial_spectra: NDArray | None
        Initial guess of spectra of pure compounds

    initial_concentrations: NDArray | None
        Initial guess of concentrations of pure compounds

    rtol : float
        Tolerance required for convergence

    Returns
    -------

    NDArray
        Optimized concentrations [component, time]

    NDArray
        Optimized spectra [component, wavelength], normalized such that mean = 1

    float
        Mean squared error

    """

    if isinstance(data, Data2D):
        y = data.data
    else:
        y = data

    if initial_spectra is not None:
        c = concentrations_from_spectra(y, initial_spectra)[0]
        s = initial_spectra
    elif initial_concentrations is not None:
        s = spectra_from_concentrations(y, initial_concentrations)[0]
        c = initial_concentrations
    elif n_compounds is not None:
        s = guess_spectra(y, n_compounds)
        c = concentrations_from_spectra(y, s)[0]
    else:
        raise Exception("One of `n_compounds`, `initial_spectra`, `initial concentrations` must be provided.")
    
    s = (s.T / np.mean(s, axis=1)).T

    old_s = s
    old_c = c

    converged = False

    for _ in range(max_iter):
        # Update spectra
        s = spectra_from_concentrations(y, c)[0]
        # Normalize spectra
        s = (s.T / np.mean(s, axis=1)).T
        # Update concentrations 
        c, mse = concentrations_from_spectra(y, s)

        # Check convergence
        change_s = np.linalg.norm(old_s - s) / np.linalg.norm(s)
        change_c = np.linalg.norm(old_c - c) / np.linalg.norm(c)
        if change_s < rtol and change_c < rtol:
            converged = True
            break
    
        old_s = s
        old_c = c

    if not converged:
        warnings.warn("alternating_lstsq() did not converge, maximum number of iterations reached")
    
    return c, s, mse