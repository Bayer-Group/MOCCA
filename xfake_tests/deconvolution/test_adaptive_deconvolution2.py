from typing import Tuple
from numpy.typing import NDArray

from matplotlib import pyplot as plt # type: ignore
import numpy as np

from mocca2 import load_data2d, estimate_baseline, find_peaks, deconvolve_adaptive
from mocca2.peaks import split_peaks
from mocca2.deconvolution.peak_models import BiGaussian, BiGaussianTailing, FraserSuzuki

print("Loading data")
# load data
blank = load_data2d('tests/test_data/09072021_gradient_97.txt')
sample = load_data2d('tests/test_data/09072021_sample_7.txt')

assert blank.check_same_sampling(sample)

print("Correcting baseline")
# correct baseline
sample -= blank
start_idx, _ = sample.closest_time(2)
sample.data = sample.data[:,start_idx:]
sample.time = sample.time[start_idx:]
sample -= estimate_baseline(sample, method='arpls')

contracted = sample.contract()

print("Finding peaks")
# Find peaks
peaks = find_peaks(contracted, min_rel_height=0.02, merge_overlapping=True)
peaks = split_peaks(contracted, peaks, 0.05)

print(f"Found {len(peaks)} peaks")

print("Deconvolving peaks")
# Deconvolve all peaks

concentrations = []
times = []
for peak in peaks:
    peak_data = peak.data(sample)
    max_mse = np.mean(peak_data**2) * 0.03
    # max_mse = max(total_ms*0.05*0.01, peak_ms*0.01)
    min_comps = min(3, len(peak.all_maxima))

    concs, spectra, mse = deconvolve_adaptive(peak_data, FraserSuzuki(), max_mse, relaxe_concs=False, min_comps=min_comps, max_comps=4)
    if mse > max_mse:
        print("Failed to deconvolve peak at", sample.time[peak.maximum])
    
    concentrations += list(concs)
    times += [peak.time(sample)] * len(concs)

print("Done deconvolving")
print(f"Found {len(peaks)} peaks and {len(concentrations)} components")

plt.figure()
plt.plot(sample.time, contracted)
for t,c in zip(times, concentrations):
    plt.plot(t, c)

plt.show()