from __future__ import annotations
from typing import List, Tuple, Dict, Any
from numpy.typing import NDArray

from scipy.ndimage import gaussian_filter  # type: ignore
from scipy.signal import find_peaks  # type: ignore
import numpy as np

from mocca2.serializing import dict_encoder


class Compound:
    """Information about single chemical compound"""

    elution_time: int
    """Index of the elution time on the time scale"""

    spectrum: NDArray
    """Absorption spectrum of the compound, normalized to mean = 1"""

    name: str | None
    """Name of the compound"""

    concentration_factor: float | None
    """Conversion factor to get absolute concentration, such that `concentration = integral * concentration_factor`"""

    concentration_factor_vs_istd = float | None
    """Conversion factor to get absolute concentration relative to ISTD `rel_conc = istd_conc * concentration_factor_vs_istd`"""

    _absorption_maxima: List[Tuple[int, float]] | None
    """Cached indeces of absorption maxima. Access this by the absorption_maxima() function"""

    def __init__(
        self,
        elution_time: int,
        spectrum: NDArray,
        name: str | None = None,
        concentration_factor: float | None = None,
        concentration_factor_vs_istd: float | None = None,
    ):
        self.elution_time = elution_time
        self.spectrum = spectrum
        self.name = name
        self.concentration_factor = concentration_factor
        self.concentration_factor_vs_istd = concentration_factor_vs_istd
        self._absorption_maxima = None

    def absorption_maxima(self) -> List[Tuple[int, float]]:
        """
        Finds absorption maxima using 2nd derivatives. Returns indeces of absorption maxima and relative heights.
        """

        if self._absorption_maxima is None:
            self._absorption_maxima = []

            second_deriv = gaussian_filter(
                self.spectrum, len(self.spectrum) / 40, order=2
            )
            peaks, _ = find_peaks(-second_deriv, rel_height=0.3)
            peaks = sorted(peaks)
            max_height = np.max(self.spectrum)
            rel_heights = [self.spectrum[p] / max_height for p in peaks]
            for p, h in zip(peaks, rel_heights):
                if h > 0.03:
                    self._absorption_maxima.append((p, h))

        return self._absorption_maxima

    def to_dict(self) -> Dict[str, Any]:
        """Converts the data to a dictionary for serialization"""
        data = self.__dict__.copy()
        data["spectrum"] = data["spectrum"]
        return dict_encoder(data | {"__classname__": "Compound"})

    @staticmethod
    def from_dict(data: Dict[str, Any]) -> Compound:
        """Creates a Compound object from a dictionary"""
        assert data["__classname__"] == "Compound"

        return Compound(
            int(data["elution_time"]),
            np.array(data["spectrum"]),
            data["name"],
            (
                float(data["concentration_factor"])
                if data["concentration_factor"] is not None
                else None
            ),
            (
                float(data["concentration_factor_vs_istd"])
                if data["concentration_factor_vs_istd"] is not None
                else None
            ),
        )
