from matplotlib import pyplot as plt  # type: ignore

from mocca2 import load_data2d

# LabSolutions
ls_blank = load_data2d('tests/test_data/09072021_gradient_97.txt')
ls_sample = load_data2d('tests/test_data/09072021_istd_96.txt')
assert ls_blank.check_same_sampling(ls_sample)

# Empower
em_blank = load_data2d('tests/test_data/raw_data3725.arw')
em_sample = load_data2d('tests/test_data/raw_data3734.arw')
assert em_blank.check_same_sampling(em_sample)

# ChemStation
cs_blank = load_data2d('tests/test_data/2022-01-26_19-43-27_gradient.D/DAD1.CSV')
cs_sample = load_data2d('tests/test_data/2022-01-26_20-00-51_ba_0.5.D')
assert cs_blank.check_same_sampling(cs_sample)

# Plot some of the data
ls_corrected = ls_sample - ls_blank
wl_idx, _ = ls_corrected.closest_wavelength(230.)

t = ls_corrected.time
plt.plot(t, ls_corrected.data[wl_idx])
plt.plot(t, ls_sample.data[wl_idx])
plt.plot(t, ls_blank.data[wl_idx])
plt.show()
