from typing import Tuple
from numpy.typing import NDArray

import numpy as np

# from scipy.linalg import solve, LinAlgWarning
from mocca2.deconvolution.nnls import nnls
import warnings


def concentrations_from_spectra(
    data: NDArray, spectra: NDArray
) -> Tuple[NDArray, float]:
    """
    Calculates positive concentrations such that the mean squared error is minimized. Returns concentrations and MSE

    Parameters
    ----------
    data: NDArray
        Absorbance data [wavelength, time]

    spectra: NDArray
        Spectra of the individual compounds [compound, wavelength]

    Returns
    -------
    NDArray
        Positive concentrations that minimize the MSE [compound, time]
    float
        MSE

    """

    assert (
        data.shape[0] == spectra.shape[1]
    ), "The shapes do not match in concentrations_from_spectra()"

    concs = []
    error = 0
    for t in range(data.shape[1]):
        try:
            conc, loss = nnls(spectra.T, data[:, t])
        except:
            # if failed to converge, solve with l2 norm
            conc, loss = nnls(spectra.T, data[:, t], l2=1e-3)

        concs.append(conc)
        error += loss**2

    concentrations = np.array(concs).T
    mse = error / data.shape[0] / data.shape[1]

    return concentrations, mse


def spectra_from_concentrations(
    data: NDArray, concentrations: NDArray
) -> Tuple[NDArray, float]:
    """
    Calculates positive spectra such that the mean squared error is minimized. Returns spectra and MSE

    Parameters
    ----------
    data: NDArray
        Absorbance data [wavelength, time]

    concentrations: NDArray
        Concentrations of the individual compounds [compound, time]

    Returns
    -------
    NDArray
        Positive spectra that minimize the MSE [compound, wavelength]
    float
        MSE
    """

    assert (
        data.shape[1] == concentrations.shape[1]
    ), "The shapes do not match in spectra_from_concentrations()"

    # Mathematically, this is identical problem, so reuse the other function
    return concentrations_from_spectra(data.T, concentrations)
