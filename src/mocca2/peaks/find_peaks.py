from typing import List
from numpy.typing import NDArray

import numpy as np
from scipy.signal import find_peaks as scipy_find_peaks  # type: ignore
from scipy.ndimage import gaussian_filter  # type: ignore

from mocca2.classes import Peak
from mocca2.peaks import merge_overlapping_peaks
from mocca2.peaks import split_peaks


def find_peaks(
        data: NDArray,
        min_rel_height: float = 0.01,
        min_height: float = -np.inf,
        width_at: float = 0.1,
        expand_borders: bool = True,
        merge_overlapping: bool = True,
        split_threshold: float | None = 0.05,
        min_elution_time: int | None = None,
        max_elution_time: int | None = None
) -> List[Peak]:
    """
    Finds all peaks in given 1D data. Assumes that baseline is flat and centered around 0.

    Parameters
    ----------
    data: ArrayLike
        1D chromatogram data

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
    List[Peak]

    Description
    -----------

    1. The peaks are picked using scipy.signal.find_peaks and filtered based on `min_rel_height`
    2. If `min_elution_time` or `max_elution_time` are specified, the peaks are filtered
    3. If `expand_borders`, the borders of the peaks are expanded down to baseline (up to estimated background noise)
    4. If `merge_overlapping`, any overlapping peaks are merged. See `merge_overlapping_peaks`
    5. If `split_threshold` is provided, merged peaks with sufficient minimum separating them are split. See `split_peaks`

    """

    # Find peaks
    peaks = _initial_peak_picking(data, min_rel_height, min_height, width_at)

    if expand_borders:
        # estimate background noise before filtering peaks
        background_std = _estimate_background_noise(data, peaks)

    # Filter peaks by elution time
    if min_elution_time is None:
        min_elution_time = 0
    if max_elution_time is None:
        max_elution_time = len(data)

    peaks = [
        peak for peak in peaks
        if peak.maximum >= min_elution_time and peak.maximum <= max_elution_time
    ]

    # Expand peak borders
    if expand_borders:
        _expand_peaks(data, peaks, background_std)

    # Merge overlapping
    if merge_overlapping:
        peaks = merge_overlapping_peaks(data, peaks)

    # Try splitting merged
    if split_threshold is not None:
        peaks = split_peaks(data, peaks, split_threshold)

    return peaks


def _initial_peak_picking(data: NDArray, min_rel_height: float = 0.01, min_height: float = -np.inf, width_at: float = 0.1) -> List[Peak]:
    """Runs scipy.signal.find_peaks and filters peaks by relative prominence"""

    # Find peaks
    maxima, info = scipy_find_peaks(
        data, height=-np.inf, prominence=min_height, width=-np.inf, rel_height=1.-width_at)

    if len(maxima) == 0:
        return []

    # Filter by relative height
    max_prominence = np.max(info['prominences'])
    keep = info['prominences']/max_prominence > min_rel_height

    maxima = maxima[keep]
    for key in info:
        info[key] = info[key][keep]

    # Create instances of Peak
    peaks = []
    for idx in range(len(maxima)):
        peak = Peak(
            left=int(info['left_ips'][idx]),
            right=int(info['right_ips'][idx]),
            maximum=maxima[idx],
            height=info['peak_heights'][idx],
            prominence=info['prominences'][idx]
        )
        peaks.append(peak)

    return peaks


def _estimate_background_noise(data: NDArray, peaks: List[Peak], mult: float = 4) -> float:
    """Iteratively calculates std outside of peaks and crops data at multiple of std"""

    y = np.array(data)

    keep = np.ones(y.shape) == 1.
    for peak in peaks:
        keep[peak.left:peak.right] = False

    y = y[keep]

    # avoid empty
    if len(y) == 0:
        y = np.array(data)
    min_std = np.std(y) * 0.01

    old_std = 0
    while True:
        std = np.std(y)
        if std < min_std:
            return min_std
        if np.abs(old_std - std)/std < 1e-5:
            return std
        y = np.clip(y, -mult*std, mult*std)
        old_std = std


def _expand_peaks(data: NDArray, peaks: List[Peak], max_cutoff: float) -> None:
    """Expands peak borders until derivative changes sign and signal decreases under max_cutoff"""

    smooth = gaussian_filter(data, 1.)
    deriv = gaussian_filter(data, 1., order=1)

    for p in peaks:
        while p.left > 0 and (smooth[p.left] > max_cutoff or deriv[p.left] > 0):
            p.left -= 1

        while p.right < smooth.shape[0] - 1 and (smooth[p.right] > max_cutoff or deriv[p.right] < 0):
            p.right += 1
