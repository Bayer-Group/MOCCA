"""High-level interface for processing single chromatogram"""

from __future__ import annotations
from typing import List, Literal, Dict, Callable, Any

import matplotlib.axes
import numpy as np
import matplotlib

from mocca2.classes import Data2D, Peak, DeconvolvedPeak, Component, Compound
from mocca2 import parsers
from mocca2.baseline import estimate_baseline
from mocca2.peaks import find_peaks
from mocca2.deconvolution.deconvolve import deconvolve_adaptive
from mocca2.deconvolution.fit_peak_model import fit_peak_model
from mocca2.deconvolution.nonnegative_lstsq import concentrations_from_spectra
from mocca2.deconvolution.peak_models import (
    PeakModel,
    Bemg,
    FraserSuzuki,
    BiGaussian,
    BiGaussianTailing,
)
from mocca2.serializing import dict_encoder


class Chromatogram(Data2D):
    """Information about single chromatogram, based on Data2D"""

    peaks: List[Peak | DeconvolvedPeak]
    """Peaks in the chromatogram"""

    sample_path: str | None
    """Filename of the raw chromatogram file"""

    blank_path: str | None
    """Filename of the raw chromatogram file with blank"""

    name: str | None
    """Name of this chromatogram"""

    def __init__(
        self,
        sample: Data2D | str,
        blank: Data2D | str | None = None,
        name: str | None = None,
        interpolate_blank=False,
    ):
        """
        Creates chromatogram from the given sample. Substracts blank if provided.

        Parameters
        ----------

        sample: Data2D | str
            Data of the sample, or path to the raw chromatogram file

        blank: Data2D | str | None
            Data of the blank, or path to the raw chromatogram file

        name: str | None
            Name of this chromatogram

        interpolate_blank: bool
            If True and the time sampling of sample and blank are different, the blank is interpolated to match the sampling of the sample

        """

        # Initialize sample data
        if isinstance(sample, Data2D):
            self.__dict__ = sample.__dict__
        else:
            self.__dict__ = parsers.load_data2d(sample).__dict__

        # Substract blank
        if blank is not None:
            if isinstance(blank, Data2D):
                blank_data = blank
            else:
                blank_data = parsers.load_data2d(blank)

            assert self.check_same_sampling(
                blank_data, time=False
            ), "The wavelength sampling of the sample and blank are different"

            if not self.check_same_sampling(blank_data, wavelength=False):
                assert (
                    interpolate_blank
                ), "The time sampling of the sample and blank are different, and `interpolate_blank` is False"
                blank_data = blank_data.interpolate_time(self.time)

            self.data -= blank_data.data

        # Other fields
        self.peaks = []
        self.name = name
        self.sample_path = sample if isinstance(sample, str) else None
        self.blank_path = blank if isinstance(blank, str) else None

    def correct_baseline(
        self,
        method: Literal["asls", "arpls", "flatfit"] = "flatfit",
        smoothness: float = 1.0,
        p: float | None = None,
        tol: float = 1e-7,
        max_iter: int | None = None,
        smooth_wl: int | None = None,
    ) -> Chromatogram:
        """
        Corrects the baseline using AsLS, arPLS or FlatFit algorithm

        Parameters
        ----------
        data: NDArray | Data2D
            Data with shape [N] or [sample, N]

        method: Literal['asls', 'arpls', 'flatfit']
            Possible baseline estimation methods are AsLS, arPLS and FlatFit. FlatFit and AsLS work well with smooth data, asPLS works better with noisy data

        smoothness: float
            size of smoothness penalty

        p: float | None
            Assymetry factor, different for AsLS and arPLS

        tol: float
            maximum relative change of `w` for convergence

        max_iter: int | None
            maximum number of iterations. If not specified, guessed automatically

        smooth_wl: int | None
            if specified, applies Savitzky-Golay filter (order 2) accross wavelength axis with given window size

        Returns
        -------
        Chromatogram
            Returns self
        """

        self.data -= estimate_baseline(
            self, method, smoothness, p, tol, max_iter, smooth_wl
        ).data

        return self

    def find_peaks(
        self,
        contraction: Literal["mean", "max", "weighted_mean"] = "mean",
        min_rel_height: float = 0.01,
        min_height: float = 10.0,
        width_at: float = 0.1,
        expand_borders: bool = True,
        merge_overlapping: bool = True,
        split_threshold: float | None = 0.05,
        min_elution_time: float | None = None,
        max_elution_time: float | None = None,
    ) -> Chromatogram:
        """
        Finds all peaks in contracted data. Assumes that baseline is flat and centered around 0.

        Parameters
        ----------
        contraction: Literal['mean', 'max', 'weighted_mean']
            Contraction method to project 2D data into 1D

        min_rel_height: float
            minimum relative prominence of the peaks (relative to highest peak)

        min_height: float
            minimum prominence of the peaks

        width_at: float
            the peak width will be measured at this fraction of peak height

        expand_borders: bool
            if True, tries to find peak borders. Otherwise borders from scipy are returned

        merge_overlapping: bool
            if True, also calls the merge_overlapping_peaks before returning the peaks

        split_threshold: float | None
            maximum height of a minimum separating two peaks for splitting, relative to smaller of the peaks

        min_elution_time: int | None
            if specified, peaks with maximum before `min_elution_time` are omitted

        max_elution_time: int | None
            if specified, peaks with maximum after `max_elution_time` are omitted

        Returns
        -------
        Chromatogram
            Returns self

        Description
        -----------

        1. The peaks are picked using scipy.signal.find_peaks and filtered based on `min_rel_height`
        2. If `min_elution_time` or `max_elution_time` are specified, the peaks are filtered
        3. If `expand_borders`, the borders of the peaks are expanded down to baseline (up to estimated background noise)
        4. If `merge_overlapping`, any overlapping peaks are merged. See `merge_overlapping_peaks`
        5. If `split_threshold` is provided, merged peaks with sufficient minimum separating them are split. See `split_peaks`

        """

        contracted = self.contract(contraction)

        self.peaks = find_peaks(
            contracted,
            min_rel_height,
            min_height,
            width_at,
            expand_borders,
            merge_overlapping,
            split_threshold,
            (
                self.closest_time(min_elution_time)[0]
                if min_elution_time is not None
                else None
            ),
            (
                self.closest_time(max_elution_time)[0]
                if max_elution_time is not None
                else None
            ),
        )

        return self

    def deconvolve_peaks(
        self,
        model: (
            PeakModel
            | Literal["BiGaussian", "BiGaussianTailing", "FraserSuzuki", "Bemg"]
        ),
        min_r2: float,
        relaxe_concs: bool,
        max_comps: int,
    ) -> Chromatogram:
        """
        Deconvolves peaks with increasingly more components until MSE limit is reached. See `deconvolve_adaptive()` for details.

        Parameters
        ----------
        model: PeakModel | Literal['BiGaussian', 'BiGaussianTailing', 'FraserSuzuki']
            mathematical model used for fitting shapes of components of peaks

        min_r2: float
            Minimum required R2 for deconvolution

        relaxe_concs: bool
            If False, the fitted peak model functions are returned. Otherwise, the concentrations are refined with restricted least squares

        max_comps: int
            Maximum number of components that can be fitted

        Returns
        -------
        Chromatogram
            Returns self
        """

        if len(self.peaks) == 0:
            return self

        base_ms = np.mean([np.mean(peak.data(self.data) ** 2) for peak in self.peaks])

        for idx, peak in enumerate(self.peaks):
            # prepare input data
            peak_data = peak.data(self)
            peak_ms = np.mean(peak_data**2)

            # if peak MS (mean square) is smaller than average, use average to avoid small peaks having too many components
            ms = max(peak_ms, base_ms)
            max_mse = (1 - min_r2) * ms
            min_comps = min(max(1, len(peak.all_maxima)), max_comps)

            # deconvolve
            concs, spectra, mse = deconvolve_adaptive(
                peak_data, model, max_mse, relaxe_concs, min_comps, max_comps
            )

            # save deconvolved peak
            r2 = 1 - mse / peak_ms
            self.peaks[idx] = DeconvolvedPeak(
                peak=peak,
                concentrations=concs,
                spectra=spectra,
                residual_mse=mse,
                r2=r2,
                resolved=max_mse > mse,
            )

        return self

    def all_components(
        self, sort_by: Callable[[Component], Any] | None = None
    ) -> List[Component]:
        """Returns all peak components from this chromatogram"""

        components = []
        for peak in self.peaks:
            if isinstance(peak, DeconvolvedPeak):
                components.extend(peak.components)

        if sort_by is not None:
            components = sorted(components, key=sort_by)

        return components

    def get_area_percent(self, wl_idx: int) -> Dict[int, float]:
        """
        Returns area % of individual components. Only deconvolved peaks are considered.

        Parameters
        ----------

        wl_idx: int
            Index of wavelength which will be used for calculating peak area

        Returns
        -------
        Dict[int, float]
            Area % of individual compounds [compound_id -> area %]

        """

        integrals: Dict[int, float] = {}

        # Get areas of individual compounds
        for component in self.all_components():
            id = component.compound_id
            if id is None:
                continue
            if id in integrals:
                integrals[id] += component.get_area(wl_idx)
            else:
                integrals[id] = component.get_area(wl_idx)

        if len(integrals) == 0:
            return {}

        # Normalize the area
        sum = np.sum([integ for integ in integrals.values()]) / 100.0

        integrals = {id: integ / sum for id, integ in integrals.items()}

        return integrals

    def get_integrals(self) -> Dict[int, float]:
        """
        Returns integrals of individual components [compound_id -> area %]. Only deconvolved peaks are considered.
        """

        integrals: Dict[int, float] = {}

        # Get integrals of individual compounds
        for component in self.all_components():
            id = component.compound_id
            if id is None:
                continue
            if id in integrals:
                integrals[id] += component.integral
            else:
                integrals[id] = component.integral

        return integrals

    def get_relative_integrals(self, relative_to: int) -> Dict[int, float]:
        """
        Returns integrals of individual components relative to the specified compound.
        If reference compound is not present, returns empty dictionary.
        Only deconvolved peaks are considered.

        Parameters
        ----------

        relative_to: int
            All integrals will be divided by integral of compound with this ID

        Returns
        -------
        Dict[int, float]
            Relative integrals of individual compounds [compound_id -> area %]

        """
        # Get integrals of all compounds
        integrals = self.get_integrals()

        # Check that reference compound is present
        if relative_to not in integrals or integrals[relative_to] <= 0:
            return {}

        # Divide integrals by the integral of the reference compound
        rel_integral = integrals[relative_to]
        integrals = {id: integ / rel_integral for id, integ in integrals.items()}

        return integrals

    def refine_peaks(
        self,
        compounds: Dict[int, Compound],
        model: (
            PeakModel
            | Literal["BiGaussian", "BiGaussianTailing", "FraserSuzuki", "Bemg"]
        ),
        relaxe_concs: bool,
        min_rel_integral: float,
    ) -> Chromatogram:
        """
        Refines the concentration profiles using averaged spectra of the compounds.

        Removes components with insuffucient integrals.
        """

        if isinstance(model, str):
            model = {
                "BiGaussian": BiGaussian,
                "BiGaussianTailing": BiGaussianTailing,
                "FraserSuzuki": FraserSuzuki,
                "Bemg": Bemg,
            }[model]()

        # fine tune concentration profiles
        for peak in self.peaks:
            spectra = np.array(
                [
                    compounds[component.compound_id].spectrum
                    for component in peak.components
                ]
            )
            data = peak.data(self.data)

            if relaxe_concs:
                concs, _ = concentrations_from_spectra(data, spectra)
            else:
                concs = np.array(
                    [component.concentration for component in peak.components]
                )
                concs, _, _ = fit_peak_model(
                    peak.data(self.data),
                    model,
                    spectra,
                    adjust_spectra=False,
                    initial_concs=concs,
                )

            total_conc = np.sum(concs)
            for conc, comp, spectrum in zip(concs, peak.components, spectra):
                comp.concentration = conc
                comp.integral = np.sum(conc)
                comp.spectrum = spectrum
                comp.peak_fraction = np.sum(conc) / total_conc

        # remove components with insufficient integrals
        min_integral = (
            max([0] + [component.integral for component in self.all_components()])
            * min_rel_integral
        )

        # remove small components
        for peak in self.peaks:
            peak.components = [
                component
                for component in peak.components
                if component.integral > min_integral
            ]

        # remove peaks with no components
        self.peaks = [peak for peak in self.peaks if len(peak.components) > 0]

        return self

    # serialization and deserialization
    def to_dict(self) -> Dict[str, Any]:
        """Converts the data to a dictionary for serialization"""
        data = super().to_dict()
        data["peaks"] = [peak.to_dict().copy() for peak in self.peaks]
        data["sample_path"] = self.sample_path
        data["blank_path"] = self.blank_path
        data["name"] = self.name

        return dict_encoder(data | {"__classname__": "Chromatogram"})

    @staticmethod
    def from_dict(data: Dict[str, Any]) -> Chromatogram:
        """Creates a Chromatogram object from a dictionary"""
        assert data["__classname__"] == "Chromatogram"
        data2d = Data2D.from_dict(data | {"__classname__": "Data2D"})
        peaks = [
            (
                Peak.from_dict(peak)
                if peak["__classname__"] == "Peak"
                else DeconvolvedPeak.from_dict(peak)
            )
            for peak in data["peaks"]
        ]
        sample_path = data["sample_path"]
        blank_path = data["blank_path"]
        name = data["name"]

        chrom = Chromatogram(data2d)
        chrom.peaks = peaks
        chrom.sample_path = sample_path
        chrom.blank_path = blank_path
        chrom.name = name
        return chrom

    # plotting
    def plot(
        self,
        ax: matplotlib.axes.Axes = None,
        color: str = "k",
        label: str | None = None,
        plot_peaks: bool = True,
        zero_line: bool = False,
    ) -> matplotlib.axes.Axes:
        ax = super().plot(ax=ax, color=color, label=label, zero_line=zero_line)

        if not plot_peaks:
            return ax

        # add peaks
        colors = [
            "#1f77b4",
            "#ff7f0e",
            "#2ca02c",
            "#d62728",
            "#9467bd",
            "#8c564b",
            "#e377c2",
            "#7f7f7f",
            "#bcbd22",
            "#17becf",
        ]

        idx = 0
        contracted = self.contract()

        for peak in self.peaks:
            if isinstance(peak, DeconvolvedPeak):
                for component in peak.components:
                    time = peak.time(self.time)
                    intensity = component.concentration
                    ax.plot(
                        time,
                        intensity,
                        color=colors[idx % len(colors)],
                        linestyle="-",
                    )
                    ax.fill_between(
                        time,
                        intensity,
                        color=colors[idx % len(colors)],
                        alpha=0.3,
                    )
                    idx += 1
            else:
                mid = (self.time[peak.left] + self.time[peak.right]) / 2
                stretch = (-self.time[peak.left] + self.time[peak.right]) / 2
                ax.errorbar(
                    [mid],
                    [0],
                    xerr=[stretch],
                    color=colors[idx % len(colors)],
                    capsize=3,
                )
                for max in peak.all_maxima:
                    ax.vlines(
                        self.time[max],
                        0,
                        contracted[max],
                        colors=colors[idx % len(colors)],
                        lw=1.0,
                    )
                idx += 1

        return ax
