from typing import Tuple
from numpy.typing import NDArray

from matplotlib import pyplot as plt # type: ignore
import numpy as np
from scipy.optimize import curve_fit, minimize # type: ignore

from mocca2 import load_data2d, estimate_baseline
from mocca2.deconvolution.nonnegative_lstsq import concentrations_from_spectra, spectra_from_concentrations
from mocca2.deconvolution.guess_spectra import guess_spectra
from mocca2.deconvolution.fit_peak_model_spectra import fit_peak_model_and_spectra
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

plt.imshow(peak_data, aspect='auto')
plt.show()

print("Guessing spectra")
spectra_guess = guess_spectra(peak_data, 2)
concs_guess, mse_guess = concentrations_from_spectra(peak_data, spectra_guess)

print("Fitting peak model")
model = BiGaussian()
concs, spectra, mse_bigaussian  = fit_peak_model_and_spectra(peak_data, spectra_guess, model)

print("MSE:", mse_guess, mse_bigaussian)

for s in spectra:
    plt.plot(sample.wavelength, s)
plt.show()

colors = ['b', 'r']
for cg, c, color in zip(concs_guess, concs, colors):
    plt.plot(peak_time, cg, color+'--')
    plt.plot(peak_time, c, color+'-')
plt.show()