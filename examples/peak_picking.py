from mocca2 import Chromatogram, find_peaks
from matplotlib import pyplot as plt

# Load data
chromatogram = Chromatogram('tests/test_data/tripeak.arw', 'tests/test_data/tripeak_blank.arw')
chromatogram.correct_baseline()

# Plot the chromatogram
plt.figure()
plt.plot(chromatogram.time, chromatogram.contract())
plt.xlabel('Time [min]')
plt.ylabel('Average absorbance [mAU]')
plt.xlim(chromatogram.time[0], chromatogram.time[-1])

# Cropping unrelevant parts of the chromatogram
chromatogram.extract_time(2, 10, inplace=True)
chromatogram.extract_wavelength(220, 400, inplace=True)

# Plot the chromatogram
plt.figure()
plt.plot(chromatogram.time, chromatogram.contract())
plt.xlabel('Time [min]')
plt.ylabel('Average absorbance [mAU]')
plt.xlim(chromatogram.time[0], chromatogram.time[-1])

# Pick peaks
contracted = chromatogram.contract()

peaks = find_peaks(
    contracted,
    min_rel_height=0.03,
    expand_borders=True,
    merge_overlapping=True,
)

# Plot chromatogram with the peaks
plt.figure()
plt.plot(chromatogram.time, contracted)
plt.xlabel('Time [min]')
plt.ylabel('Average absorbance [mAU]')
plt.xlim(chromatogram.time[0], chromatogram.time[-1])

t = chromatogram.time
for peak in peaks:
    plt.hlines(-10, t[peak.left], t[peak.right], colors='red')
    for max in peak.all_maxima:
        plt.vlines(t[max], -10, contracted[max], colors='green', lw=1.)

# Pick peaks
contracted = chromatogram.contract()

peaks = find_peaks(
    contracted,
    min_rel_height=0.01,
    expand_borders=True,
    merge_overlapping=True,
)

# Plot chromatogram with the peaks
plt.figure()
plt.plot(chromatogram.time, contracted)
plt.xlabel('Time [min]')
plt.ylabel('Average absorbance [mAU]')
plt.xlim(chromatogram.time[0], chromatogram.time[-1])

t = chromatogram.time
for peak in peaks:
    plt.hlines(-10, t[peak.left], t[peak.right], colors='red')
    for max in peak.all_maxima:
        plt.vlines(t[max], -10, contracted[max], colors='green', lw=1.)

plt.show()