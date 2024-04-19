# See detailed description for this example in docs

from mocca2 import Chromatogram, example_data
import numpy as np
from matplotlib import pyplot as plt

PATH_TO_SAMPLE_CHROMATOGRAM = "src/mocca2/example_data/data/examples/chrom1.arw"
PATH_TO_BLANK_CHROMATOGRAM = "src/mocca2/example_data/data/examples/blank1.arw"

# Load raw data from file
# Including the blank is optional, but makes the processing easier
chromatogram = Chromatogram(
    sample=PATH_TO_SAMPLE_CHROMATOGRAM,
    blank=PATH_TO_BLANK_CHROMATOGRAM,
)

# Alternatively, you can use example data
chromatogram = example_data.example_1()

# Crop wavelength between 220 and 400 nm
chromatogram.extract_wavelength(220, 400, inplace=True)

# Plot chromatogram (averaged over wavelenghts)
fig, ax = plt.subplots(figsize=(8, 5))
chromatogram.plot(ax)
# plt.savefig("docs/_static/ex_basic_chromatogram.svg")
plt.show()


# Plot 2D chromatogram as heatmap
fig, ax = plt.subplots(figsize=(8, 5))
chromatogram.plot_2d(ax)
# plt.savefig("docs/_static/ex_basic_chromatogram_2d.svg")
plt.show()

# Extract region with peak
peak_region = chromatogram.extract_time(2.2, 2.35)

# Average spectrum
# The chromatogram data are just numpy arrays, axes are [wavelength, time]
avg_spectrum = np.mean(peak_region.data, axis=1)
plt.subplots(figsize=(8, 5))
plt.plot(peak_region.wavelength, avg_spectrum)
plt.xlabel("Wavelength [nm]")
plt.ylabel("Average absorbance [mAU]")
# plt.savefig("docs/_static/ex_basic_spectrum.svg")
plt.show()
