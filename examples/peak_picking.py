from mocca2 import example_data
from matplotlib import pyplot as plt
import pandas as pd

# Load the chromatogram and correct baseline
chromatogram = example_data.example_2()
chromatogram.correct_baseline()

# Crop out the relevant part of the chromatogram
chromatogram.extract_time(1.5, None, inplace=True)

# Find the peaks
chromatogram.find_peaks()

# Plot chromatogram with the peaks
chromatogram.plot()
plt.savefig("docs/_static/ex_peak_picking_initial.svg")
plt.show()

# It is possible to adjust the threshold to pick even the smaller peaks
chromatogram.find_peaks(min_height=1)

# Let's see the chromatogram again
chromatogram.plot()
plt.savefig("docs/_static/ex_peak_picking_adjusted.svg")
plt.show()

# As you can see, some the algorithm detected that some peaks overlap
# Let's deconvolve them
chromatogram.deconvolve_peaks(
    model="FraserSuzuki",
    min_r2=0.99,
    relaxe_concs=False,
    max_comps=3,
)

# And plot the chromatogram again
chromatogram.plot()
plt.savefig("docs/_static/ex_peak_picking_deconvolved.svg")
plt.show()

# Finally, let's get all the components and calculate their area %
# (the integral is averaged over all wavelengths)
data = pd.DataFrame(
    [
        {
            "Elution Time [min]": chromatogram.time[component.elution_time],
            "Integral": component.integral,
        }
        for component in chromatogram.all_components()
    ]
)
data["Area %"] = data["Integral"] / data["Integral"].sum() * 100
data.sort_values("Elution Time [min]", inplace=True)
print(data.round(2).to_string(index=False))
