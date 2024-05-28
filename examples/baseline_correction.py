# See detailed description for this example in docs

from mocca2 import example_data, estimate_baseline
from matplotlib import pyplot as plt

# Load example chromatogram
chromatogram = example_data.example_1(substract_blank=True)

# Refine the baseline
chromatogram.correct_baseline()

# Load the chromatogram without baseline correction and without blank subtraction
chromatogram_no_baseline = example_data.example_1(substract_blank=True)
chromatogram_no_blank = example_data.example_1(substract_blank=False)

# Plot the chromatogram with corrected baseline
fig, ax = plt.subplots(figsize=(8, 5))
chromatogram_no_blank.plot(ax, color="green", label="No blank subtraction")
chromatogram_no_baseline.plot(ax, color="red", label="No baseline correction")
chromatogram.plot(ax, label="Corrected")

plt.legend()
# plt.savefig("docs/_static/ex_baseline_corrected.svg")
plt.show()

# We can explore different baseline correction algorithms

# Load example chromatogram
chromatogram = example_data.example_1(substract_blank=False)

# To make things faster, lets average absorbance over all wavelengths
mean_absorbance = chromatogram.contract()

# Estimate baseline using different methods
baseline_asls = estimate_baseline(mean_absorbance, method="asls")
baseline_arpls = estimate_baseline(mean_absorbance, method="arpls")
baseline_flatfit = estimate_baseline(mean_absorbance, method="flatfit")

# Plot the result
fig, ax = plt.subplots(figsize=(8, 5))

chromatogram.plot(ax, label="Original")
ax.plot(chromatogram.time, baseline_arpls, label="AsLS")
ax.plot(chromatogram.time, baseline_asls, label="arPLS")
ax.plot(chromatogram.time, baseline_flatfit, label="FlatFit")

ax.set_ylim(-30, 75)
plt.legend()
# plt.savefig("docs/_static/ex_baseline_comparison.svg")
plt.show()
