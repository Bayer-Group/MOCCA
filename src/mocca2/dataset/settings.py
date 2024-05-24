from __future__ import annotations
from typing import Literal, Dict, Any
from dataclasses import dataclass
import yaml  # type: ignore

from mocca2.serializing import dict_encoder


@dataclass(init=True)
class ProcessingSettings:
    """Collection of all settings required for automatic chromatogram processing in MoccaDataset"""

    baseline_model: Literal["asls", "arpls", "flatfit"] = "flatfit"
    """Name of baseline estimator"""
    baseline_smoothness: float = 1.0
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
    peak_model: Literal["BiGaussian", "BiGaussianTailing", "FraserSuzuki", "Bemg"] = (
        "Bemg"
    )
    """Model that describes the peak shape"""
    max_peak_comps: int = 4
    """Maximum number of deconvolved components in single peak"""
    max_peak_distance: float = 1.0
    """Maximum peak distance deviation relative to peak width for one compound"""
    min_spectrum_correl: float = 0.99
    """Minimum correlation of spectra for one compound"""
    min_elution_time: float = 0.4
    """Peaks with maxima before this time will not be considered"""
    max_elution_time: float = 10.0
    """Peaks with maxima after this time will not be considered"""
    min_wavelength: float = 220.0
    """The data will be cropped such that lower wavelengths are not included"""
    max_wavelength: float = 400.0
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

    def to_dict(self) -> Dict[str, Any]:
        """Converts the data to a dictionary for serialization"""
        return dict_encoder(self.__dict__ | {"__classname__": "ProcessingSettings"})

    @staticmethod
    def from_dict(data: Dict[str, Any]) -> ProcessingSettings:
        """Creates a ProcessingSettings object from a dictionary"""
        assert data["__classname__"] == "ProcessingSettings"

        data.pop("__classname__")
        settings = ProcessingSettings(
            baseline_model=str(data["baseline_model"]),
            baseline_smoothness=float(data["baseline_smoothness"]),
            min_rel_prominence=float(data["min_rel_prominence"]),
            min_prominence=float(data["min_prominence"]),
            border_max_peak_cutoff=float(data["border_max_peak_cutoff"]),
            split_threshold=float(data["split_threshold"]),
            explained_threshold=float(data["explained_threshold"]),
            peak_model=str(data["peak_model"]),
            max_peak_comps=int(data["max_peak_comps"]),
            max_peak_distance=float(data["max_peak_distance"]),
            min_spectrum_correl=float(data["min_spectrum_correl"]),
            min_elution_time=float(data["min_elution_time"]),
            max_elution_time=float(data["max_elution_time"]),
            min_wavelength=float(data["min_wavelength"]),
            max_wavelength=float(data["max_wavelength"]),
            min_rel_integral=float(data["min_rel_integral"]),
            relaxe_concs=bool(data["relaxe_concs"]),
        )
        return settings
