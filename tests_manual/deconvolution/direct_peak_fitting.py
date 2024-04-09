from typing import Tuple
from numpy.typing import NDArray
import warnings

import numpy as np
from scipy.optimize import minimize # type: ignore
from scipy.signal import find_peaks # type: ignore

from mocca2.deconvolution.peak_models import PeakModel
from mocca2.deconvolution.nonnegative_lstsq import concentrations_from_spectra, spectra_from_concentrations
from mocca2.deconvolution.fit_peak_model import fit_peak_model

from mocca2 import Chromatogram
from time import time
from matplotlib import pyplot as plt

from mocca2.deconvolution.peak_models import FraserSuzuki, Bemg

sample = Chromatogram('tests/test_data/tripeak.arw', 'tests/test_data/tripeak_blank.arw')

print("Correcting baseline")
# correct baseline
sample.correct_baseline('flatfit', smoothness=1.)

peak = sample.extract_time(1.44, 1.54)


print("Deconvolving peak")
start = time()
concs, mse, _ = fit_peak_model(peak.data, Bemg(), n_compounds=3)
print("MSE:",mse)
spectra = spectra_from_concentrations(peak.data, concs)[0]
concs = (concs.T * np.mean(spectra, axis=1)).T

print("Done deconvolving, elapsed", time() - start)

plt.figure()
plt.plot(peak.time, peak.contract())
for c in concs:
    plt.plot(peak.time, c)

plt.show()