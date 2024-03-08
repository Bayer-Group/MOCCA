from typing import Tuple, Literal
from numpy.typing import NDArray

import numpy as np

from mocca2.deconvolution.nonnegative_lstsq import concentrations_from_spectra, spectra_from_concentrations
from mocca2.deconvolution.fit_peak_model import fit_peak_model
from mocca2.deconvolution.peak_models import PeakModel, BiGaussian, BiGaussianTailing, FraserSuzuki, Bemg


def deconvolve_adaptive(
        data: NDArray,
        model: PeakModel | Literal['BiGaussian', 'BiGaussianTailing', 'FraserSuzuki', 'Bemg'],
        max_mse: float,
        relaxe_concs: bool,
        min_comps: int,
        max_comps: int
) -> Tuple[NDArray, NDArray, float]:
    """
    Deconvolves data with increasingly more components until MSE limit is reached

    Parameters
    ----------
    data: NDArray
        2D data [wavelength, data]

    model: PeakModel | Literal['BiGaussian', 'BiGaussianTailing', 'FraserSuzuki']
        mathematical model used for fitting shapes of components of peaks

    max_mse: float
        Maximum allowed MSE for termination

    relaxe_concs: bool
        If False, the fitted peak model functions are returned

        Otherwise, the concentrations are refined with restricted least squares

    min_comps: int
        Minimum number of components that can be fitted

    max_comps: int
        Maximum number of components that can be fitted

    Returns
    -------
    NDArray
        concentrations [compound, time]

    NDArray
        spectra [compound, wavelength], normalized such that mean = 1

    float
        MSE

    """

    if isinstance(model, str):
        model = {
            'BiGaussian': BiGaussian(),
            'BiGaussianTailing': BiGaussianTailing(),
            'FraserSuzuki': FraserSuzuki(),
            'Bemg': Bemg(),
        }[model]

    for n_comps in range(min_comps, max_comps+1):
        # Deconvolve peak with some increasing number of components
        concs, spectra, mse = deconvolve_fixed(
            data, n_comps, model, relaxe_concs
        )
        # Check whether MSE is sufficiently small
        if mse < max_mse:
            break

    return concs, spectra, mse


def deconvolve_fixed(data: NDArray, n_comps: int, model: PeakModel, relaxe_concs: bool) -> Tuple[NDArray, NDArray, float]:
    """
    Deconvolves data with given number of components. Returns concentration, spectra and MSE

    Parameters
    ----------
    data: NDArray
        2D data [wavelength, data]

    n_comps: int
        how many components should be used for deconvolution

    model: PeakModel
        mathematical model used for fitting shapes of components of peaks

    relaxe_concs: bool
        If False, the fitted bigaussian functions are returned

        Otherwise, the concentrations are refined with restricted least squares

    Returns
    -------
    NDArray
        concentrations [compound, time]

    NDArray
        spectra [compound, wavelength], normalized such that mean = 1

    float
        MSE

    """

    # Fit peak model
    concs, mse, _ = fit_peak_model(data, model, n_compounds=n_comps)

    # Get spectra
    spectra, mse = spectra_from_concentrations(data, concs)

    if relaxe_concs:
        # Relaxe the constrain on peaks being bigaussian
        concs, mse = concentrations_from_spectra(data, spectra)
        spectra, mse = spectra_from_concentrations(data, concs)

    # Normalize the spectra and scale concentrations accordingly
    norm_factors = np.mean(spectra, axis=1) + 1e-7
    spectra = (spectra.T / norm_factors).T
    concs = (concs.T * norm_factors).T

    return concs, spectra, mse
