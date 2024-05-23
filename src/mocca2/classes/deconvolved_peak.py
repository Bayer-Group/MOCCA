from __future__ import annotations
from typing import List, Dict, Any
from numpy.typing import NDArray

import numpy as np

from mocca2.classes import Peak, Component
from mocca2.serializing import dict_encoder


class DeconvolvedPeak(Peak):
    """Information about peak and its deconvolved components"""

    components: List[Component]
    """Deconvolved components of the peak"""

    residual_mse: float
    """Residual MSE after deconvolution"""

    r2: float
    """R2 after deconvolution"""

    resolved: bool
    """This specifies, whether the deconvolution sufficiently explains the peak"""

    def __init__(
        self,
        peak: Peak,
        concentrations: NDArray,
        spectra: NDArray,
        residual_mse: float,
        r2: float,
        resolved: bool,
    ):
        """
        Initializes DeconvolvedPeak based on the Peak and deconvolved data

        Parameters
        ----------
        peak: Peak
            The base Peak object that was deconvolved

        concentrations: NDArray = None
            Concentrations of the individual components [component, time]

        spectra: NDArray = None
            Spectra of the individual components [component, wavelength]

        residual_mse: float
            Residual MSE after deconvolution

        r2: float
            R2 after deconvolution

        resolved: bool
            This specifies, whether the deconvolution sufficiently explains the peak
        """

        self.__dict__ = peak.__dict__

        self.residual_mse = residual_mse
        self.r2 = r2
        self.components = []
        self.resolved = resolved

        assert (
            concentrations.shape[0] == spectra.shape[0]
        ), "The number of components in concentrations and spectra does not match"

        total_conc = np.sum(concentrations)

        for conc, spec in zip(concentrations, spectra):
            peak_frac = np.sum(conc) / total_conc
            self.components.append(Component(conc, spec, peak.left, peak_frac))

    def merge_same_components(self) -> None:
        """
        Merges all the components with identical ID (not None).

        Concentrations are added, spectra are averaged, weighted by concentration integral.
        """

        # Assign all components to a dict by compound_id
        components: Dict[int, List[Component]] = dict()
        for component in self.components:
            id = component.compound_id
            if id in components:
                components[id].append(component)
            elif id is not None:
                components[id] = [component]

        # Merge components with identical compound_id
        for compound_id, comps in components.items():
            if len(comps) == 1 or compound_id is None:
                continue

            integrals = np.array([c.integral for c in comps])
            spectra = np.array([c.spectrum for c in comps])

            conc = np.sum([c.concentration for c in comps], axis=0)
            spectrum = np.sum(spectra.T * integrals, axis=1) / np.sum(integrals)

            peak_frac = np.sum(c.peak_fraction for c in comps)

            merged = Component(conc, spectrum, self.left, peak_frac, compound_id)

            components[compound_id] = [merged]

        # Flatten the dict into list and save
        self.components = [
            component for comps in components.values() for component in comps
        ]

    def to_dict(self) -> Dict[str, Any]:
        """Converts the data to a dictionary for serialization"""
        data = super().to_dict().copy()
        data["components"] = [c.to_dict().copy() for c in self.components]
        data["residual_mse"] = self.residual_mse
        data["r2"] = self.r2
        data["resolved"] = self.resolved

        data["__classname__"] = "DeconvolvedPeak"
        return dict_encoder(data)

    @staticmethod
    def from_dict(data: Dict[str, Any]) -> DeconvolvedPeak:
        """Creates a DeconvolvedPeak object from a dictionary"""
        assert data["__classname__"] == "DeconvolvedPeak"

        peak = Peak.from_dict(data | {"__classname__": "Peak"})
        residual_mse = float(data["residual_mse"])
        r2 = float(data["r2"])
        resolved = bool(data["resolved"])
        components = [Component.from_dict(c) for c in data["components"]]
        concs = np.array([c.concentration for c in components])
        specs = np.array([c.spectrum for c in components])

        deconvolved_peak = DeconvolvedPeak(
            peak, concs, specs, residual_mse, r2, resolved
        )
        deconvolved_peak.components = components
        return deconvolved_peak
