from __future__ import annotations
from typing import Literal
from dataclasses import dataclass
from copy import deepcopy

import yaml # type: ignore

@dataclass(init=True)
class ProcessingSettings:
    """Collection of all settings required for automatic chromatogram processing in MoccaDataset"""

    baseline_model: Literal['asls', 'arpls', 'flatfit'] = 'flatfit'
    """Name of baseline estimator"""
    baseline_smoothness: float = 1.
    """Smoothness penalty for baseline"""
    min_rel_prominence: float = 0.01
    """Minimal relative peak height"""
    min_prominence: float = 1
    """Minimal peak height"""
    border_max_peak_cutoff: float = 0.1
    """Maximum relative peak height for peak cutoff"""
    split_threshold: float = 0.05
    """Maximum relative height of minima between peaks to split them"""
    explained_threshold: float = 0.995
    """Minimal R2 to consider peak resolved"""
    peak_model: Literal['BiGaussian', 'BiGaussianTailing', 'FraserSuzuki', 'Bemg'] = 'Bemg'
    """Model that describes the peak shape"""
    max_peak_comps: int = 4
    """Maximum number of deconvolved components in single peak"""
    max_peak_distance: float = 1.
    """Maximum peak distance deviation relative to peak width for one compound"""
    min_spectrum_correl: float = 0.99
    """Minimum correlation of spectra for one compound"""
    min_elution_time: float = 0.4
    """Peaks with maxima before this time will not be considered"""
    max_elution_time: float = 10.
    """Peaks with maxima after this time will not be considered"""
    min_wavelength: float = 220.
    """The data will be cropped such that lower wavelengths are not included"""
    max_wavelength: float = 400.
    """The data will be cropped such that higher wavelengths are not included"""
    min_rel_integral: float = 0.01
    """Minimum integral relative to the largest peak"""
    relaxe_concs: bool = False
    """If True, the concentrations will be relaxed to fit the calibration curve without any peak model"""

    def to_yaml(self) -> str:
        """Converts self to YAML string"""
        data = yaml.dump(self.__dict__)
        return data

    @staticmethod
    def from_yaml(data: str) -> ProcessingSettings:
        loaded = yaml.safe_load(data)
        settings = ProcessingSettings(**loaded)
        return settings

    def to_json(self):
        json_dict = deepcopy(self.__dict__)
        return json_dict

    def from_json(json_dict):
        settings = ProcessingSettings(**json_dict)
        return settings