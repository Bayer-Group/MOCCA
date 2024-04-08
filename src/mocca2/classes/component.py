from __future__ import annotations
from typing import Any, Dict
from numpy.typing import NDArray

import numpy as np


class Component:
    """Information about single deconvolved component of a peak"""

    concentration: NDArray
    """Concentration profile in the selected range"""

    spectrum: NDArray
    """Spectrum of the component. Normalized such that `mean = 1`"""

    compound_id: int | None
    """ID of compound, if assigned"""

    elution_time: int
    """Index of time point with maximum concentration"""

    integral: float
    """Integral (sum of individual time points) of the concentration of this component"""

    peak_fraction: float
    """Fraction of the peak area that this component represents"""

    def __init__(
        self,
        concentration: NDArray,
        spectrum: NDArray,
        time_offset: int = 0,
        peak_fraction: float = 1.0,
        compound_id: int | None = None,
    ):
        self.concentration = concentration
        self.spectrum = spectrum
        self.elution_time = int(np.argmax(concentration)) + time_offset
        self.integral = np.sum(concentration)
        self.compound_id = compound_id
        self.peak_fraction = peak_fraction

    def get_area(self, wl_idx: int) -> float:
        """Returns peak area at given wavelength (specified by index)"""

        return self.integral * self.spectrum[wl_idx]

    def to_dict(self) -> Dict[str, Any]:
        """Converts the data to a dictionary for serialization"""
        data = self.__dict__
        data["spectrum"] = data["spectrum"].tolist()
        data["concentration"] = data["concentration"].tolist()
        data["__classname__"] = "Component"
        return data

    @staticmethod
    def from_dict(data: Dict[str, Any]) -> Component:
        """Creates a Component object from a dictionary"""
        assert data["__classname__"] == "Component"

        component = Component(
            np.array(data["concentration"]),
            np.array(data["spectrum"]),
            int(data["elution_time"]),
            float(data["peak_fraction"]),
            data["compound_id"],
        )
        component.integral = data["integral"]
        component.elution_time = int(data["elution_time"])
        return component
