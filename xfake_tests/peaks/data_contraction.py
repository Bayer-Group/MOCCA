from matplotlib import pyplot as plt  # type: ignore

import numpy as np

from mocca2 import load_data2d, estimate_baseline

# Load data
ls_blank = load_data2d('tests/test_data/09072021_gradient_97.txt')
# ls_sample = load_data2d('tests/test_data/09072021_istd_96.txt')
ls_sample = load_data2d('tests/test_data/09072021_sample_4.txt')
assert ls_blank.check_same_sampling(ls_sample)
# ls_blank = load_data2d('tests/test_data/raw_data3725.arw')
# ls_sample = load_data2d('tests/test_data/raw_data3734.arw')
# assert ls_blank.check_same_sampling(ls_sample)

# Correct baseline
ls_corrected = ls_sample - ls_blank

ls_corrected -= estimate_baseline(ls_corrected, method='arpls')

# Contract data to 1D
d1 = ls_corrected.contract(method='mean')
d3 = ls_corrected.contract(method='max')
d4 = ls_corrected.contract(method='weighted_mean', damping=0.2)


# Plot results
t = ls_corrected.time

plt.plot(t, d1/np.max(d1), label='mean')
# plt.plot(t, d3/np.max(d3), label='max')
plt.plot(t, d4/np.max(d4), label='weighted_mean')

plt.legend()
plt.show()