from matplotlib import pyplot as plt  # type: ignore

import numpy as np

from mocca2 import load_data2d, find_peaks, estimate_baseline
from mocca2.peaks import find_peaks, merge_overlapping_peaks

# Load data
ls_blank = load_data2d('tests/test_data/09072021_gradient_97.txt')
ls_sample = load_data2d('tests/test_data/09072021_istd_96.txt')
assert ls_blank.check_same_sampling(ls_sample)
# ls_blank = load_data2d('tests/test_data/raw_data3725.arw')
# ls_sample = load_data2d('tests/test_data/raw_data3734.arw')
# assert ls_blank.check_same_sampling(ls_sample)

ls_corrected = ls_sample - ls_blank

# Contract data to 1D
averaged = np.mean(ls_corrected.data, axis=0)

averaged -= estimate_baseline(averaged, method='arpls')

# Find peaks
h = 0.1
peaks = find_peaks(averaged, min_rel_height=0.01, width_at=h)
merged = merge_overlapping_peaks(averaged, peaks)

# Plot results
t = ls_corrected.time

plt.plot(t, averaged, 'k-')
plt.plot(t[[0,-1]], [0,0], 'k--', lw=0.5, alpha=0.7)
plt.vlines(t[[p.maximum for p in peaks]], [p.height - p.prominence for p in peaks], [p.height for p in peaks], colors='green', lw=1)
plt.hlines([p.height - (1-h)*p.prominence for p in peaks], t[[p.left for p in peaks]], t[[p.right for p in peaks]], colors='red', lw=1)
plt.hlines(range(len(merged)), t[[p.left for p in merged]], t[[p.right for p in merged]], colors='orange', lw=1)

plt.show()