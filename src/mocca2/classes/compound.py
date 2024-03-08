from typing import List, Tuple
from numpy.typing import NDArray
from copy import deepcopy

from scipy.ndimage import gaussian_filter # type: ignore
from scipy.signal import find_peaks # type: ignore
import numpy as np 

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
            concentration_factor_vs_istd: float | None = None
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

            second_deriv = gaussian_filter(self.spectrum, len(self.spectrum)/40, order=2)
            peaks, _ = find_peaks(-second_deriv, rel_height=0.3)
            peaks = sorted(peaks)
            max_height = np.max(self.spectrum)
            rel_heights = [self.spectrum[p]/max_height for p in peaks]
            for p, h in zip(peaks, rel_heights):
                if h > 0.03:
                    self._absorption_maxima.append((p, h))
                
        return self._absorption_maxima

    def to_json(self):
        json_dict = deepcopy(self.__dict__)
        json_dict['spectrum'] = json_dict['spectrum'].tolist()
        return json_dict

    def from_json(json_dict):
        compound_init_vars = ['elution_time', 'spectrum', 'name', 'concentration_factor', 'concentration_factor_vs_istd']
        compound_init_dict = {key: json_dict[key] for key in compound_init_vars}
        compound_init_dict['spectrum'] = np.array(compound_init_dict['spectrum'])
        compound = Compound(**compound_init_dict)
        compound._absorption_maxima = json_dict['_absorption_maxima']
        return compound
