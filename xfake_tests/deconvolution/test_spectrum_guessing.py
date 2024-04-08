from matplotlib import pyplot as plt # type: ignore
import numpy as np
# from sklearn.decomposition import PCA # type: ignore
from sklearn.cluster import SpectralClustering # type: ignore

from mocca2 import load_data2d, estimate_baseline
from mocca2.deconvolution.nonnegative_lstsq import concentrations_from_spectra
from mocca2.deconvolution.guess_spectra import guess_spectra

# load data
blank = load_data2d('tests/test_data/09072021_gradient_97.txt')
sample = load_data2d('tests/test_data/09072021_sample_5.txt')
assert blank.check_same_sampling(sample)

# correct baseline
sample -= blank
sample -= estimate_baseline(sample, method='arpls')

# Isolate region with two overlapping peaks
sample.data = sample.data[:,980:1060]
sample.time = sample.time[980:1060]

spectra = guess_spectra(sample.data, 2)

spectra = (spectra.T / np.mean(spectra, axis=1)).T

concentrations, mse = concentrations_from_spectra(sample.data, spectra)

plt.plot(sample.time, sample.contract())
for c in concentrations:
    plt.plot(sample.time, c)
plt.figure()

# plt.imshow(-sample.data, cmap='gist_rainbow', aspect='auto')
# plt.figure()

for s in spectra:
    plt.plot(sample.wavelength, s)
plt.show()
