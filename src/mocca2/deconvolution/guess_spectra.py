"""Routines from here are used to get initial guess of spectra of individual compounds from single peak"""

from numpy.typing import NDArray

import numpy as np
from sklearn.cluster import SpectralClustering # type: ignore

from mocca2.math import cosine_similarity
from mocca2.peaks import find_peaks

def guess_spectra(data: NDArray, n_compounds: int) -> NDArray:
    """
    Clusters data in the peak to try to isolate dominant spectra.

    Parameters
    ----------
    data: NDArray
        2D peak data [wavelength, time]
    
    n_compounds: int
        how many spectra should be generated

    Returns
    -------
    NDArray
        The generated spectra [compound, wavelength]. Spectra are not normalized

    """

    # Find peaks
    avg = np.mean(data, axis=0)
    peaks = find_peaks(avg, expand_borders=False, merge_overlapping=False, split_threshold=None)
    peaks = sorted(peaks, key=lambda p: -p.prominence)

    # If there are enough peaks, initial guesses are just spectra at peak maxima
    if len(peaks) >= n_compounds:
        maxima = [p.maximum for p in peaks[:n_compounds]]
        return data[:,maxima].T
    
    # Trim data near base-line at 1/3 height of smallest peak
    if len(peaks) > 0:
        idx = len(peaks)-1
        trim_height = peaks[idx].height * 0.3
    else:
        trim_height = np.max(data) * 0.3
    y = data[:, avg > trim_height]

    # Get similarity matrix
    similarity_matrix = _get_similarity_matrix(y)

    # Get spectra of the averaged clusters
    spectra = _get_clustered_spectra(y, similarity_matrix, n_compounds)

    return spectra


def _get_similarity_matrix(data: NDArray) -> NDArray:
    """Creates similarity matrix with samples along last axis"""

    N = data.shape[1]

    similarity = np.zeros([N,N])
    for i in range(N):
        s = cosine_similarity(data[:,i], data[:, i+1:].T)
        similarity[i, i+1:] = s
    
    similarity += similarity.T

    return similarity + 1. + 1e-7

def _get_clustered_spectra(data: NDArray, similarity_matrix: NDArray, n_compounds: int) -> NDArray:
    """Clusters time points based on cosine similarity of spectra and returns the dominant spectra"""
    clustering = SpectralClustering(n_clusters=n_compounds,affinity='precomputed')
    clusters = clustering.fit_predict(similarity_matrix)

    spectra = []
    for cluster in set(clusters):
        sp = np.mean(data[:,clusters==cluster], axis=1)
        sp /= np.linalg.norm(sp)
        spectra.append(sp)

    return np.array(spectra)