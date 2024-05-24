"""Routines for fitting peak models to peak data"""

from typing import Tuple
from numpy.typing import NDArray
import warnings

import numpy as np
from scipy.optimize import minimize  # type: ignore
from scipy.signal import find_peaks  # type: ignore
from scipy.signal import savgol_filter  # type: ignore

from mocca2.deconvolution.peak_models import PeakModel
from mocca2.deconvolution.nonnegative_lstsq import (
    concentrations_from_spectra,
    spectra_from_concentrations,
)
from mocca2.deconvolution.guess_spectra import guess_spectra


def fit_peak_model(
    data: NDArray,
    model: PeakModel,
    spectra: NDArray | None = None,
    n_compounds: int | None = None,
    adjust_spectra: bool = True,
    initial_params: None | NDArray = None,
    initial_concs: None | NDArray = None,
) -> Tuple[NDArray, float, NDArray]:
    """
    Fits peak shapes to data. Spectra can be fixed or fitted too.

    Parameters
    ----------
    data: NDArray
        2D chromatogram data, [wavelength, time]

    model: PeakModel
        mathematical model of the peak shape

    spectra: NDArray | None
        absorption spectra of individual components, [components, wavelength]. If not provided, guess_spectra() is used.

    n_compounds: int | None = None
        number of peaks to deconvolve, must be provided if spectra is None. Ignored if spectra are provided

    adjust_spectra: bool
        whether the spectra should be also fitted. When False, the spectra are kept constant

    initial_params: None | NDArray
        initial parameters for the models, flattened

    initial_concs: None | NDArray
        initial guess of the concentrations

    Returns
    -------
    NDArray:
        Concentrations [component, time]

    float
        Mean Squared Error

    NDArray:
        The optimized parameters, flattened

    """
    if spectra is not None:
        assert (
            data.shape[0] == spectra.shape[1]
        ), "The shapes of data and spectra don't match in fit_peak_model()"
    else:
        spectra = guess_spectra(data, n_compounds)

    # Some generally useful variables
    t = np.arange(data.shape[1], dtype=float)
    n_comps = spectra.shape[0]
    n_params = model.n_params()

    bounds = model.get_bounds(t[-1])

    # Get initial guess of the model parameters
    if initial_params is None:
        p0 = []

        if initial_concs is None:
            # adjust spectra guess to mean=1
            spectra = (spectra.T / np.mean(spectra, axis=1)).T
            # add peaks and remove their contribution one by one
            residual_data = data.copy()
            for peak_index in range(spectra.shape[0]):
                # get concentrations guess
                concentrations, _ = concentrations_from_spectra(
                    residual_data, spectra[peak_index:]
                )
                # smooth out the concentrations
                for idx in range(concentrations.shape[0]):
                    concentrations[idx] = savgol_filter(
                        concentrations[idx], max(data.shape[1] // 20, 5), 2
                    )

                # find highest peak
                component, _ = np.unravel_index(
                    np.argmax(concentrations), concentrations.shape
                )

                conc = concentrations[component]
                peaks, info = find_peaks(conc, width=2, height=0, rel_height=0.3)
                if len(peaks) > 0:
                    idx = np.argmax(info["peak_heights"])
                    height = info["peak_heights"][idx]
                    mu = peaks[idx]
                    s1 = peaks[idx] - info["left_ips"][idx]
                    s2 = info["right_ips"][idx] - peaks[idx]
                # otherwise use moments
                else:
                    height = np.max(conc)
                    mu = np.sum(conc * t) / (np.sum(conc) + 1e-7)
                    s1 = np.sqrt(
                        np.sum(conc * (t - mu) ** 2 * (t < mu))
                        / (np.sum(conc * (t < mu)) + 1e-7)
                    )
                    s2 = np.sqrt(
                        np.sum(conc * (t - mu) ** 2 * (t < mu))
                        / (np.sum(conc * (t < mu)) + 1e-7)
                    )
                # get the peak parameters
                params = model.init_guess(height, mu, s1 * 1.4, s2 * 1.4)
                p0.append(params)

                # substract data accounted for by the peak
                residual_data -= np.outer(
                    spectra[component + peak_index], model.val(t, *params)
                )

                # swap remove the spectrum guess
                spectra[[peak_index, component + peak_index]] = spectra[
                    [component + peak_index, peak_index]
                ]
        else:
            concentrations = initial_concs

            for conc in concentrations:
                # if peak is present use it to get initial parameters
                peaks, info = find_peaks(conc, width=2, height=0.0, rel_height=0.3)
                if len(peaks) > 0:
                    idx = np.argmax(info["peak_heights"])
                    height = info["peak_heights"][idx]
                    mu = peaks[idx]
                    s1 = peaks[idx] - info["left_ips"][idx]
                    s2 = info["right_ips"][idx] - peaks[idx]
                # otherwise use moments
                else:
                    height = np.max(conc)
                    mu = np.sum(conc * t) / (np.sum(conc) + 1e-7)
                    s1 = np.sqrt(
                        np.sum(conc * (t - mu) ** 2 * (t < mu))
                        / (np.sum(conc * (t < mu)) + 1e-7)
                    )
                    s2 = np.sqrt(
                        np.sum(conc * (t - mu) ** 2 * (t < mu))
                        / (np.sum(conc * (t < mu)) + 1e-7)
                    )

                # get initial guess
                params = model.init_guess(height, mu, s1 * 1.4, s2 * 1.4)

                # make sure the initial guess is within bounds
                for i, (low, high) in enumerate(bounds):
                    params[i] = np.clip(params[i], low, high)

                p0.append(params)

        initial_params = np.array(p0).flatten()

    def loss_and_grad(
        params: NDArray, data: NDArray, spectra: NDArray, t: NDArray
    ) -> Tuple[float, NDArray]:
        """Returns MSE and gradient. Would be good to test that this is correct."""
        # Reshape params to [component, param]
        if not adjust_spectra:
            params = params.reshape([n_comps, n_params])
        else:
            params = params.reshape([n_comps, n_params - 1])
            params = np.concatenate(
                [np.ones([n_comps, 1], dtype=float), params], axis=1
            )

        # Get concentrations [component, time]
        concs = np.array([model.val(t, *p) for p in params])

        # Get spectra
        if adjust_spectra:
            spectra = spectra_from_concentrations(data, concs)[0]

        # Get residues and mse
        residues = spectra.T @ concs - data
        mse = np.mean(residues**2)

        # Get gradients of reconstructed data [param, wl, time]
        if not adjust_spectra:
            sTc_grad = np.zeros([n_comps * n_params, data.shape[0], data.shape[1]])
            for comp_idx, param in enumerate(params):
                sTc_grad[comp_idx * n_params : (comp_idx + 1) * n_params] = np.einsum(
                    "w,pt->pwt", spectra[comp_idx], model.grad(t, *param)
                )
        else:
            sTc_grad = np.zeros(
                [n_comps * (n_params - 1), data.shape[0], data.shape[1]]
            )
            for comp_idx, param in enumerate(params):
                sTc_grad[
                    comp_idx * (n_params - 1) : (comp_idx + 1) * (n_params - 1)
                ] = np.einsum("w,pt->pwt", spectra[comp_idx], model.grad(t, *param)[1:])

        # Get gradients fo square residues [param, wavele, time]
        res_grad = 2 * sTc_grad * residues

        # Average over wavelength and time axis
        grad = np.mean(res_grad, axis=(1, 2))

        # Return MSE
        return mse, grad

    # Prepare for the optimization
    if not adjust_spectra:
        bounds = bounds * n_comps
    else:
        bounds = bounds[1:] * n_comps

    # remove height from params if spectra can be adjusted
    if adjust_spectra:
        initial_params = initial_params.reshape([n_comps, n_params])
        initial_params = initial_params[:, 1:]
        initial_params = initial_params.flatten()

    # These methods seem to work well: L-BFGS-B (with increased maxls), TNC, SLSQP
    result = minimize(
        lambda params: loss_and_grad(params, data, spectra, t),
        x0=initial_params,
        jac=True,
        bounds=bounds,
        method="L-BFGS-B",
        options={"maxls": 100},
        tol=1e-10,
    )

    # it is not a big issue if the optimization does not converged, if this function is called iteratively
    if not result.success:
        warnings.warn(
            "fit_peak_model(): Optimization failed to converge: " + result.message
        )

    params = result.x
    # Reshape params to [component, param]
    if not adjust_spectra:
        params = params.reshape([len(params) // n_params, n_params])
    else:
        n_peaks = len(params) // (n_params - 1)
        params = params.reshape([n_peaks, n_params - 1])
        params = np.concatenate([np.ones([n_peaks, 1], dtype=float), params], axis=1)
    # Get concentrations [component, time]
    concs = np.array([model.val(t, *p) for p in params])

    return concs, result.fun, params.flatten()
