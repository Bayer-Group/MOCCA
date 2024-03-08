from numpy.typing import NDArray
from copy import deepcopy

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

    def __init__(self, concentration: NDArray, spectrum: NDArray, time_offset: int = 0, peak_fraction: float = 1., compound_id: int | None = None):
        self.concentration = concentration
        self.spectrum = spectrum
        self.elution_time = int(np.argmax(concentration)) + time_offset
        self.integral = np.sum(concentration)
        self.compound_id = compound_id
        self.peak_fraction = peak_fraction

    def get_area(self, wl_idx: int) -> float:
        """Returns peak area at given wavelength (specified by index)"""

        return self.integral * self.spectrum[wl_idx]

    def to_json(self):
        json_dict = deepcopy(self.__dict__)
        json_dict['concentration'] = json_dict['concentration'].tolist()
        json_dict['spectrum'] = json_dict['spectrum'].tolist()

        return json_dict

    def from_json(json_dict):
        component_init_vars = ['concentration', 'spectrum', 'time_offset', 'peak_fraction', 'compound_id']
        component_inst_vars = ['elution_time', 'integral']

        component_init_vars_dict = {}
        for var in component_init_vars:
            if var in json_dict.keys():
                component_init_vars_dict[var] = json_dict[var]

        component_inst_vars_dict = {}
        for var in component_inst_vars:
            if var in json_dict.keys():
                component_inst_vars_dict[var] = json_dict[var]

        component_init_vars_dict['concentration'] = np.array(component_init_vars_dict['concentration'])
        component_init_vars_dict['spectrum'] = np.array(component_init_vars_dict['spectrum'])

        component = Component(**component_init_vars_dict)
        component.elution_time = component_inst_vars_dict['elution_time']
        component.integral = component_inst_vars_dict['integral']

        return component
