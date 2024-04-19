from typing import Tuple
from numpy.typing import NDArray

from matplotlib import pyplot as plt # type: ignore
import numpy as np

from mocca2 import load_data2d, estimate_baseline
from mocca2.deconvolution.deconvolve import deconvolve_adaptive
from mocca2.deconvolution.peak_models import BiGaussian

print("Loading data")
# load data
blank = load_data2d('tests/test_data/09072021_gradient_97.txt')
sample = load_data2d('tests/test_data/09072021_sample_5.txt')
assert blank.check_same_sampling(sample)

print("Correcting baseline")
# correct baseline
sample -= blank
sample -= estimate_baseline(sample, method='arpls')

# Isolate region with two overlapping peaks
peak_data = sample.data[:,980:1060]
peak_time = sample.time[980:1060]

print("Deconvolving")
# Deconvolve
max_mse = np.mean(peak_data**2) * 0.01

concs, spectra, mse = deconvolve_adaptive(peak_data, BiGaussian(), max_mse, relaxe_concs=False, min_comps=1, max_comps=3)

print(concs.shape, spectra.shape, max_mse, mse)

print("Done")

for s in spectra:
    plt.plot(sample.wavelength, s)

plt.figure()
plt.plot(peak_time, np.mean(peak_data, axis=0))
for c in concs:
    plt.plot(peak_time, c)

plt.show()