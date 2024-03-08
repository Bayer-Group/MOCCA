from typing import List
from numpy.typing import NDArray

import numpy as np

from mocca2.classes import Peak

def split_peaks(data: NDArray, peaks: List[Peak], max_height: float = 0.05) -> List[Peak]:
    """
    Splits any peaks that are have sufficient minimum separating them

    Parameters
    ----------
    data: NDArray
        1D chromatogram data
    
    peaks:
        List of peaks, only merged peaks can be split
    
    max_height:
        maximum height of a minimum separating two peaks, relative to smaller of the peaks
    
    Returns
    -------
    List[Peak]
        split peaks
    """

    split_peaks = []
    for peak in peaks:
        split_peaks += _split_peak(data, peak, max_height)

    return split_peaks
    
def _split_peak(data: NDArray, peak: Peak, max_height: float) -> List[Peak]:
    """Splits single peak, see description of split_peaks()"""

    if len(peak.all_maxima) < 2:
        return [peak]
    
    # get indeces of peak maxima
    maxima = np.unique(np.array(sorted(peak.all_maxima), dtype=int))

    # get indeces of minima between peaks
    minima = [np.argmin(data[peak.left:maxima[0]]) + peak.left]
    for idx in range(len(maxima) - 1):
        minimum = np.argmin(data[maxima[idx]:maxima[idx+1]]) + maxima[idx]
        minima.append(minimum)
    minima += [np.argmin(data[maxima[-1]:peak.right]) + maxima[-1]]

    # determine which minima are sufficient to separate peaks
    separate = [True]*(len(maxima)+1)
    abs_min = np.min(data[peak.left:peak.right])

    for idx in range(len(maxima)):
        height = data[maxima[idx]] - abs_min
        if data[minima[idx]] - abs_min > height * max_height:
            separate[idx] = False
        if data[minima[idx+1]] - abs_min > height * max_height:
            separate[idx+1] = False

    # The peak will always start and end
    separate[0] = True
    separate[-1] = True

    # create the separated peaks
    peaks = []
    sep_start = 0
    for idx in range(len(maxima)):
        if separate[idx+1]:
            maximum = np.argmax(data[maxima[sep_start:idx+1]]) + sep_start
            maximum = maxima[maximum]
            new_peak = Peak(
                left=minima[sep_start], # type: ignore
                right=minima[idx+1]-1, # type: ignore
                height=data[maximum],
                prominence=data[maximum] - abs_min,
                maximum=maximum, # type: ignore
                all_maxima=list(maxima[sep_start:idx+1]),
            )
            sep_start = idx + 1

            peaks.append(new_peak)

    return peaks