"""The class Dataset provides very high-level interface for processing sets of chromatograms"""

from __future__ import annotations
from typing import Dict, Tuple, Any, List
from numpy.typing import NDArray

from multiprocessing import Pool
import pandas as pd  # type: ignore
import numpy as np

from mocca2.classes.chromatogram import Chromatogram
from mocca2.classes.data2d import Data2D
from mocca2.classes import Compound, Component
from mocca2.dataset.settings import ProcessingSettings
from mocca2.clustering.cluster_components import cluster_components
from mocca2.math import cosine_similarity
from mocca2.serializing import dict_encoder


class MoccaDataset:
    """Collection of chromatograms, compounds, and other information"""

    _raw_2d_data: Dict[int, Data2D]
    """2D chromatogram data without wavelength cropping"""

    chromatograms: Dict[int, Chromatogram]
    """Chromatograms in the dataset [id -> Chromatogram]"""

    compounds: Dict[int, Compound]
    """Compounds by ID [id -> Compound]"""

    compound_references: Dict[int, Tuple[str, float | None]]
    """If some chromatogram is reference for a compound, this we be stored here as [chromatogram id -> (compound name, concentration?)]"""

    istd_concentrations: Dict[int, float]
    """Concentrations of internal standard [chromatogram id -> concentration]"""

    _istd_chromatogram: int | None
    """ID of reference chromatogram for internal standard"""

    istd_compound: int | None
    """ID of compound that is internal standard"""

    settings: ProcessingSettings | None
    """Settings for automatic chromatogram processing"""

    def __init__(self):
        self.chromatograms = {}
        self._raw_2d_data = {}
        self.compounds = {}
        self.compound_references = {}
        self.istd_concentrations = {}
        self._istd_chromatogram = None
        self.istd_compound = None
        self.settings = ProcessingSettings()

    def add_chromatogram(
        self,
        chromatogram: Chromatogram,
        istd_concentration: float | None = None,
        reference_for_compound: str | None = None,
        compound_concentration: float | None = None,
        istd_reference: bool = False,
    ) -> int:
        """
        Adds chromatogram to the dataset, returns the assigned ID

        Parameters
        ----------

        chromatogram: Chromatogram
            the chromatogram that is added

        istd_concentration: float | None
            if internal standard is present, specify concentration

        reference_for_compound: str | None
            if this chromatogram contains compound reference, specify compound name

        compound_concentration: float | None
            if this chromatogram contains compound reference with known concentration, specify concentration

        istd_reference: bool = False
            specify whether this chromatogram is reference for internal standard

        """

        # Find new id for the chromatogram
        new_id = 0
        while new_id in self.chromatograms:
            new_id += 1

        # check that the sampling is same as in existing chromatograms
        if len(self.chromatograms) > 0:
            if not next(iter(self._raw_2d_data.values())).check_same_sampling(
                chromatogram
            ):
                raise Exception(
                    "Cannot add this chromatogram to the campaign, because the time or wavelength points are different"
                )

        # Save all data
        self.chromatograms[new_id] = chromatogram
        raw_data = Data2D(chromatogram.time, chromatogram.wavelength, chromatogram.data)
        self._raw_2d_data[new_id] = raw_data

        if istd_concentration is not None:
            self.istd_concentrations[new_id] = istd_concentration

        if reference_for_compound is not None:
            self.compound_references[new_id] = (
                reference_for_compound,
                compound_concentration,
            )

        if istd_reference:
            self._istd_chromatogram = new_id

        return new_id

    def time(self) -> NDArray | None:
        """Returns time axis of the chromatograms, if there are any"""
        if len(self.chromatograms) == 0:
            return None

        time = next(iter(self.chromatograms.values())).time
        return time

    def wavelength(self) -> NDArray | None:
        """Returns wavelength axis of the chromatograms, if there are any"""
        if len(self.chromatograms) == 0:
            return None

        wavelength = next(iter(self.chromatograms.values())).wavelength
        return wavelength

    def wavelength_raw(self) -> NDArray | None:
        """Returns wavelength axis of the raw data (without cropping), if there are any"""
        if len(self._raw_2d_data) == 0:
            return None

        wavelength = next(iter(self._raw_2d_data.values())).wavelength
        return wavelength

    def time_step(self) -> float | None:
        """Returns the sampling step of the time axis in the chromatograms, if there are any"""
        if len(self.chromatograms) == 0:
            return None

        time = next(iter(self.chromatograms.values())).time_step()
        return time

    def wavelength_step(self) -> float | None:
        """Returns the sampling step of the wavelength axis in the chromatograms, if there are any"""
        if len(self.chromatograms) == 0:
            return None

        wavelength = next(iter(self.chromatograms.values())).wavelength_step()
        return wavelength

    def closest_time(self, time: float) -> Tuple[int, float] | None:
        """Returns index and value of time point that is closest to specified `time`, if there are any chromatograms"""
        if len(self.chromatograms) == 0:
            return None

        data = next(iter(self.chromatograms.values())).closest_time(time)
        return data

    def closest_wavelength(self, wavelength: float) -> Tuple[int, float] | None:
        """Returns index and value of wavelength point that is closest to specified `wavelength`, if there are any chromatograms"""
        if len(self.chromatograms) == 0:
            return None

        data = next(iter(self.chromatograms.values())).closest_wavelength(wavelength)
        return data

    def _name_compounds(self):
        """Gives default name to all compounds and assigns concentration conversion factors to compounds"""

        def name_main_compound_in_chromatogram(
            chromatogram_id: int, name: str, conc: float | None
        ) -> int | None:
            chromatogram = self.chromatograms[chromatogram_id]
            components = chromatogram.all_components(sort_by=lambda c: -c.integral)

            for component in components:
                # Skip unresolved and named compounds
                if component.compound_id is None:
                    continue
                if self.compounds[component.compound_id].name is not None:
                    continue

                # Name the compound
                self.compounds[component.compound_id].name = name

                # Find concentration factor
                if conc is not None:
                    compound = self.compounds[component.compound_id]
                    integral = component.integral
                    # Absolute concentration factor
                    compound.concentration_factor = conc / integral
                    # Relative concentration factor
                    if (
                        self.istd_compound is not None
                        and chromatogram_id in self.istd_concentrations
                    ):
                        istd_integral = sum(
                            [
                                c.integral
                                for c in chromatogram.all_components()
                                if c.compound_id == self.istd_compound
                            ]
                        )
                        if istd_integral > 0:
                            istd_conc = self.istd_concentrations[chromatogram_id]
                            integral = component.integral

                            compound.concentration_factor_vs_istd = (
                                conc / integral * istd_integral / istd_conc
                            )

                            print(
                                f"Compound {compound.name} has conc factor vs ISTD {compound.concentration_factor_vs_istd:0.3f}"
                            )

                return component.compound_id
            return None

        # def name_impurities_in_chromatogram(chromatogram: Chromatogram, name: str):
        #     components = chromatogram.all_components(sort_by=lambda c: c.integral)
        #     count = 1
        #     for component in components:
        #         if component.compound_id is None:
        #             continue
        #         if self.compounds[component.compound_id].name is not None:
        #             continue
        #         self.compounds[component.compound_id].name = name + f" impurity {count}"
        #         count += 1

        # name ISTD
        if self._istd_chromatogram is not None:
            name, conc = self.compound_references[self._istd_chromatogram]
            istd_id = name_main_compound_in_chromatogram(
                self._istd_chromatogram, name, conc
            )
            self.istd_compound = istd_id

        # Set names and conversion factors for known compounds
        for idx, (name, conc) in self.compound_references.items():
            # make sure to skip ISTD
            if idx == self._istd_chromatogram:
                continue
            name_main_compound_in_chromatogram(idx, name, conc)

        # Set names for impurities in reference chromatograms
        # for idx, (name, _) in self.compound_references.items():
        #     name_impurities_in_chromatogram(self.chromatograms[idx], name)

        # Set default names for other compounds
        for compound in self.compounds.values():
            if compound.name is None:
                compound.name = f"@ {self.time()[compound.elution_time]:0.3f}"

    def process_all(
        self, settings: ProcessingSettings, verbose: bool = True, cores: int = 1
    ):
        """Processes all chromatograms: finds and deconvolves peaks, creates averaged compounds, and refines peaks"""
        self.settings = settings

        def check_nans(self):
            if any([c is None for c in self.chromatograms.values()]):
                raise Exception("Some chromatograms are missing")

        check_nans(self)

        # Reset some values
        self.compounds = {}
        for chromatogram in self.chromatograms.values():
            chromatogram.peaks = []

        if verbose:
            print("Cropping wavelengths")
        # Cropping wavelengths
        for idx, raw in self._raw_2d_data.items():
            cropped = raw.extract_wavelength(
                settings.min_wavelength, settings.max_wavelength
            )
            self.chromatograms[idx].data = cropped.data
            self.chromatograms[idx].time = cropped.time
            self.chromatograms[idx].wavelength = cropped.wavelength

        # Baseline correction
        if verbose:
            print("Correcting baseline")
        kwargs = {
            "method": settings.baseline_model,
            "smoothness": settings.baseline_smoothness,
            "smooth_wl": max(len(self.wavelength()) // 20, 4),
        }
        if cores == 1:
            for chromatogram in self.chromatograms.values():
                chromatogram.correct_baseline(**kwargs)
        else:
            with Pool(cores) as p:
                keys = list(self.chromatograms.keys())
                results = p.starmap(
                    _double_star_args,
                    [
                        (
                            Chromatogram.correct_baseline,
                            {"self": self.chromatograms[k], **kwargs},
                        )
                        for k in keys
                    ],
                )
                for id, chromatogram in zip(keys, results):
                    self.chromatograms[id] = chromatogram
        check_nans(self)

        # Peak picking
        if verbose:
            print("Picking peaks")
        kwargs = {
            "min_rel_height": settings.min_rel_prominence,
            "min_height": settings.min_prominence,
            "width_at": settings.border_max_peak_cutoff,
            "split_threshold": settings.split_threshold,
            "min_elution_time": settings.min_elution_time,
            "max_elution_time": settings.max_elution_time,
        }
        if cores == 1:
            for chromatogram in self.chromatograms.values():
                chromatogram.find_peaks(**kwargs)
        else:
            with Pool(cores) as p:
                keys = list(self.chromatograms.keys())
                results = p.starmap(
                    _double_star_args,
                    [
                        (
                            Chromatogram.find_peaks,
                            {"self": self.chromatograms[k], **kwargs},
                        )
                        for k in keys
                    ],
                )
                for id, chromatogram in zip(keys, results):
                    self.chromatograms[id] = chromatogram
        check_nans(self)

        # Initial deconvolution
        if verbose:
            print("Deconvolution")
        kwargs = {
            "model": settings.peak_model,
            "min_r2": settings.explained_threshold,
            "relaxe_concs": settings.relaxe_concs,
            "max_comps": settings.max_peak_comps,
        }
        if cores == 1:
            for idx, chromatogram in enumerate(self.chromatograms.values()):
                if verbose:
                    print(f"Chromatogram {idx+1}/{len(self.chromatograms)}")
                chromatogram.deconvolve_peaks(**kwargs)
        else:
            with Pool(cores) as p:
                keys = list(self.chromatograms.keys())
                results = p.starmap(
                    _double_star_args,
                    [
                        (
                            Chromatogram.deconvolve_peaks,
                            {"self": self.chromatograms[k], **kwargs},
                        )
                        for k in keys
                    ],
                )
                for id, chromatogram in zip(keys, results):
                    self.chromatograms[id] = chromatogram
        check_nans(self)

        # Cluster individual peaks to build averaged compounds
        if verbose:
            print("Clustering compounds")
        components = [
            component
            for chromatogram in self.chromatograms.values()
            for component in chromatogram.all_components()
        ]

        dt = self.time_step()
        assert dt is not None

        def are_same_compound(comp1: Component, comp2: Component) -> bool:
            # estimate peak width
            pw1 = (
                np.sum(
                    np.clip(
                        comp1.concentration - np.max(comp1.concentration) / 2, 0, np.inf
                    )
                    > 0
                )
                / 2
            )
            pw2 = (
                np.sum(
                    np.clip(
                        comp2.concentration - np.max(comp2.concentration) / 2, 0, np.inf
                    )
                    > 0
                )
                / 2
            )
            max_peak_dist = pw1 + pw2

            if (
                abs(comp1.elution_time - comp2.elution_time)
                > max_peak_dist * settings.max_peak_distance
            ):
                return False
            if (
                cosine_similarity(comp1.spectrum, comp2.spectrum)
                < settings.min_spectrum_correl
            ):
                return False

            return True

        def importance(comp: Component) -> float:
            return comp.integral * comp.peak_fraction**4

        self.compounds = cluster_components(
            components, are_same=are_same_compound, weights=importance
        )

        # Refine peaks
        if verbose:
            print("Refining peaks")
        kwargs = {
            "compounds": self.compounds,
            "model": settings.peak_model,
            "relaxe_concs": settings.relaxe_concs,
            "min_rel_integral": settings.min_rel_integral,
        }
        if cores == 1:
            for chromatogram in self.chromatograms.values():
                chromatogram.refine_peaks(**kwargs)
        else:
            with Pool(cores) as p:
                keys = list(self.chromatograms.keys())
                results = p.starmap(
                    _double_star_args,
                    [
                        (
                            Chromatogram.refine_peaks,
                            {"self": self.chromatograms[k], **kwargs},
                        )
                        for k in keys
                    ],
                )
                for id, chromatogram in zip(keys, results):
                    self.chromatograms[id] = chromatogram
        check_nans(self)

        # Remove compounds that are not present
        if verbose:
            print("Naming compounds")
        present = set()
        for chromatogram in self.chromatograms.values():
            for component in chromatogram.all_components():
                present.add(component.compound_id)
        self.compounds = {
            id: compound for id, compound in self.compounds.items() if id in present
        }

        # Name all compounds
        self._name_compounds()

        if verbose:
            print("Processing finished!")

    def get_area_percent(self, wl_idx: int) -> Tuple[pd.DataFrame, List[int]]:
        """
        Calculates area % of deconvolved peaks at given wavelength

        Parameters
        ----------

        wl_idx: int
            index of wavelength which will be used for calculating area %

        Returns
        -------

        pd.DataFrame
            The columns of the dataframe are: 'Chromatogram ID', 'Chromatogram', and names of compounds

        List[int]
            IDs of compounds in the same order as in DataFrame

        """

        columns = [
            list(self.chromatograms.keys()),
            [chrom.name for chrom in self.chromatograms.values()],
        ]

        column_names = ["Chromatogram ID", "Chromatogram"]

        area_percents = [
            chrom.get_area_percent(wl_idx) for chrom in self.chromatograms.values()
        ]

        compound_ids = sorted(
            self.compounds.keys(), key=lambda c_id: self.compounds[c_id].elution_time
        )

        for c_id in compound_ids:
            name = self.compounds[c_id].name
            column_names.append(name)  # type: ignore
            columns.append(
                [None if c_id not in ap else ap[c_id] for ap in area_percents]
            )

        df = pd.DataFrame(zip(*columns), columns=column_names)  # type: ignore

        return df, compound_ids

    def get_concentrations(self) -> Tuple[pd.DataFrame, List[int]]:
        """
        Calculates integrals or concentrations of deconvolved peaks.

        If compound has `concentration_factor` specified, the integrals are multiplied by this factor

        Returns
        -------

        pd.DataFrame
            The columns of the dataframe are: 'Chromatogram ID', 'Chromatogram', and names of compounds

        List[int]
            IDs of compounds in the same order as in DataFrame

        """

        columns = [
            list(self.chromatograms.keys()),
            [chrom.name for chrom in self.chromatograms.values()],
        ]

        column_names = ["Chromatogram ID", "Chromatogram"]

        integrals = [chrom.get_integrals() for chrom in self.chromatograms.values()]

        compound_ids = sorted(
            self.compounds.keys(), key=lambda c_id: self.compounds[c_id].elution_time
        )

        for c_id in compound_ids:
            name = self.compounds[c_id].name
            conc_factor = self.compounds[c_id].concentration_factor
            if conc_factor is None:
                conc_factor = 1.0
            column_names.append(name)
            columns.append(
                [
                    None if c_id not in ints else ints[c_id] * conc_factor
                    for ints in integrals
                ]
            )

        df = pd.DataFrame(zip(*columns), columns=column_names)

        return df, compound_ids

    def get_relative_concentrations(self) -> Tuple[pd.DataFrame, List[int]]:
        """
        Calculates integrals or concentrations of deconvolved peaks relative to internal standard.

        If compound has `concentration_factor` specified, the integrals are multiplied by this factor

        Returns
        -------

        pd.DataFrame
            The columns of the dataframe are: 'Chromatogram ID', 'Chromatogram', and names of compounds

        List[int]
            IDs of compounds in the same order as in DataFrame

        """

        assert (
            self.istd_compound is not None
        ), "Cannot calculate relative concentrations, the internal standard is not specified"

        columns = [
            list(self.chromatograms.keys()),
            [chrom.name for chrom in self.chromatograms.values()],
        ]

        column_names = ["Chromatogram ID", "Chromatogram"]

        # Get relative integrals (compound_integral / ISTD_integral)
        integrals = [
            chrom.get_relative_integrals(self.istd_compound)
            for chrom in self.chromatograms.values()
        ]

        # Scale the integrals by given ISTD concentration
        for idx, ch_id in enumerate(columns[0]):
            if ch_id in self.istd_concentrations:
                istd_conc = self.istd_concentrations[ch_id]
            else:
                istd_conc = 1.0
            for key in integrals[idx]:
                integrals[idx][key] *= istd_conc

        compound_ids = sorted(
            self.compounds.keys(), key=lambda c_id: self.compounds[c_id].elution_time
        )

        # put all results into dataframe
        for c_id in compound_ids:
            name = self.compounds[c_id].name
            conc_factor = self.compounds[c_id].concentration_factor_vs_istd
            if conc_factor is None:
                conc_factor = 1.0
            column_names.append(name)
            columns.append(
                [
                    None if c_id not in ints else ints[c_id] * conc_factor
                    for ints in integrals
                ]
            )

        df = pd.DataFrame(zip(*columns), columns=column_names)

        return df, compound_ids

    def to_dict(self) -> Dict[str, Any]:
        """Converts the data to a dictionary for serialization"""
        dic = {
            "chromatograms": {k: v.to_dict() for k, v in self.chromatograms.items()},
            "_raw_2d_data": {k: v.to_dict() for k, v in self._raw_2d_data.items()},
            "compounds": {k: v.to_dict() for k, v in self.compounds.items()},
            "compound_references": self.compound_references,
            "istd_concentrations": self.istd_concentrations,
            "istd_chromatogram": self._istd_chromatogram,
            "istd_compound": self.istd_compound,
            "settings": self.settings.to_dict(),
        }
        return dict_encoder(dic | {"__classname__": "MoccaDataset"})

    @staticmethod
    def from_dict(data: Dict[str, Any]) -> MoccaDataset:
        """Creates a MoccaDataset object from a dictionary"""
        assert data["__classname__"] == "MoccaDataset"

        dataset = MoccaDataset()
        dataset.chromatograms = {
            int(k): Chromatogram.from_dict(v) for k, v in data["chromatograms"].items()
        }
        dataset._raw_2d_data = {
            int(k): Data2D.from_dict(v) for k, v in data["_raw_2d_data"].items()
        }
        dataset.compounds = {
            int(k): Compound.from_dict(v) for k, v in data["compounds"].items()
        }
        dataset.compound_references = {
            int(k): (v[0], float(v[1]) if v[1] is not None else None)
            for k, v in data["compound_references"].items()
        }
        dataset.istd_concentrations = {
            int(k): float(v) for k, v in data["istd_concentrations"].items()
        }
        dataset._istd_chromatogram = (
            int(data["istd_chromatogram"])
            if data["istd_chromatogram"] is not None
            else None
        )
        dataset.istd_compound = (
            int(data["istd_compound"]) if data["istd_compound"] is not None else None
        )
        dataset.settings = (
            ProcessingSettings.from_dict(data["settings"])
            if data["settings"] is not None
            else None
        )
        return dataset


def _double_star_args(func, kwargs):
    return func(**kwargs)
