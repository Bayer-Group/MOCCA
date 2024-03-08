# See detailed description for this example in docs

from mocca2 import Chromatogram, estimate_baseline
from matplotlib import pyplot as plt

# Load sample chromatogram and substract blank
chromatogram = Chromatogram(
    'tests/test_data/tripeak.arw',
    blank='tests/test_data/tripeak_blank.arw'
)

# Refine the baseline
chromatogram.correct_baseline()

# Plotting the chromatogram with corrected baseline
plt.figure()
plt.plot(chromatogram.time, chromatogram.contract())
plt.xlabel('Time [min]')
plt.ylabel('Average absorbance [mAU]')
plt.xlim(chromatogram.time[0], chromatogram.time[-1])

# But we can explore different baseline correction algorithms
# Let's use only averaged absorbance

# Load and average data
chromatogram = Chromatogram('tests/test_data/tripeak.arw')
averaged = chromatogram.contract()
averaged_blank = Chromatogram('tests/test_data/tripeak_blank.arw').contract()

# Estimate baseline using different methods
baseline_asls = estimate_baseline(averaged, method='asls')
baseline_arpls = estimate_baseline(averaged, method='arpls')
baseline_flatfit = estimate_baseline(averaged, method='flatfit')

# Plot the result
plt.figure()
plt.plot(chromatogram.time, averaged, label='Raw data', lw=1.)
plt.plot(chromatogram.time, averaged_blank, label='Blank', lw=1.)
plt.plot(chromatogram.time, baseline_asls, label='ASLS baseline', lw=1.)
plt.plot(chromatogram.time, baseline_arpls, label='arPLS baseline', lw=1.)
plt.plot(chromatogram.time, baseline_flatfit, label='FlatFit baseline', lw=1.)
plt.xlabel('Time [min]')
plt.ylabel('Average absorbance [mAU]')
plt.xlim(chromatogram.time[0], chromatogram.time[-1])
plt.ylim(-30, 100)
plt.legend()

plt.show()