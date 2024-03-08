from typing import List
from numpy.typing import NDArray

from mocca2.classes import Peak


def merge_overlapping_peaks(data: NDArray, peaks: List[Peak]) -> List[Peak]:
    """
    Merges any overlapping peaks

    Parameters
    ----------

    data: NDArray
        1D chromatogram data

    peaks: List[Peak]
        peaks that are checked for overlap and merged

    Returns
    -------

    List[Peak]
        List of new peaks. Non-overlapping peaks are copied as-is, and overlapping peaks are replaced by merged peaks
    
    """

    if len(peaks) < 2:
        return peaks

    sorted_peaks = sorted(peaks, key=lambda p:p.left)

    merged_peaks = [peaks[0]]
    for peak in sorted_peaks[1:]:
        if _overlap(merged_peaks[-1], peak) >= 0:
            merged_peaks[-1] = _merge_two_peaks(data, merged_peaks[-1], peak)
        else:
            merged_peaks.append(peak)

    return merged_peaks


def _overlap(p1: Peak, p2: Peak) -> int:
    """Returns overlap (or negative distance) of two peaks"""
    left = max(p1.left, p2.left)
    right = min(p1.right, p2.right)
    return right-left

def _merge_two_peaks(data: NDArray, p1: Peak, p2: Peak) -> Peak:
    """Merges two peaks"""
    larger = p1 if data[p1.maximum] > data[p2.maximum] else p2

    return Peak(
        left=min(p1.left, p2.left),
        right=max(p1.right, p2.right),
        maximum=larger.maximum,
        height=larger.height,
        prominence=larger.prominence,
        all_maxima=p1.all_maxima+p2.all_maxima
    )
