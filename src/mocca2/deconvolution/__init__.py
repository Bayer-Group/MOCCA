"""
The approach for peak deconvolution is inspired by Shimadzu 10.1016/j.chroma.2016.09.037.

The entire procedure is as follows:

1. Get initial guess of spectra
    - if enough maxima are present, spectra at maxima are used
    - otherwise, the timepoints are clustered based on cosine similarity and averaged spectra of clusters are used

2. Refine spectra by fitting peak shapes to PeakModel

3. Unrestricted refinement of concentrations based on spectra

"""

from mocca2.deconvolution.deconvolve import deconvolve_adaptive, deconvolve_fixed