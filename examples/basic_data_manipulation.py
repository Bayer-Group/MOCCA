# See detailed description for this example in docs

from mocca2 import Chromatogram
import numpy as np
from matplotlib import pyplot as plt

# Load raw data from file
chromatogram = Chromatogram('tests/test_data/tripeak.arw')

# Crop wavelength between 220 and 400 nm
chromatogram.extract_wavelength(220, 400, inplace=True)

# Plot chromatogram (averaged over wavelenghts)
plt.figure()
plt.plot(chromatogram.time, chromatogram.contract())
plt.xlabel('Time [min]')
plt.ylabel('Average absorbance [mAU]')
plt.xlim(chromatogram.time[0], chromatogram.time[-1])

# Plot 2D chromatogram as heatmap
plt.figure()
plt.imshow(
    chromatogram.data,
    cmap='gist_ncar',
    extent=[
        chromatogram.time[0],
        chromatogram.time[-1],
        chromatogram.wavelength[-1],
        chromatogram.wavelength[0]
    ],
    aspect='auto'
)
plt.xlabel('Time [min]')
plt.ylabel('Wavelength [nm]')
plt.colorbar()

# Extract raegion with peak
peak_region = chromatogram.extract_time(2.2, 2.35)

# Average spectrum
# The chromatogram data are just numpy arrays, axes are [wavelength, time]
avg_spectrum = np.mean(peak_region.data, axis=1)
plt.figure()
plt.plot(peak_region.wavelength, avg_spectrum)
plt.xlabel('Wavelength [nm]')
plt.ylabel('Average absorbance [mAU]')
plt.show()