# See detailed description for this example in docs

from mocca2 import example_data, MoccaDataset, ProcessingSettings
from matplotlib import pyplot as plt
import numpy as np

# Load example data
chromatograms = example_data.cyanation()

# Create the MOCCA2 dataset
dataset = MoccaDataset()

# Specify chromatogram with with internal standard
tetralin_concentration = 0.06094
dataset.add_chromatogram(
    chromatograms["istd"],
    reference_for_compound="tetralin",
    istd_reference=True,
    compound_concentration=tetralin_concentration,
    istd_concentration=tetralin_concentration,
)

# Add standards for starting material and product
dataset.add_chromatogram(
    chromatograms["educt_1"],
    reference_for_compound="starting_material",
    compound_concentration=0.0603,
    istd_concentration=tetralin_concentration,
)
dataset.add_chromatogram(
    chromatograms["product_1"],
    reference_for_compound="product",
    compound_concentration=0.05955,
    istd_concentration=tetralin_concentration,
)

# Add the chromatograms for the reactions
for chromatogram in chromatograms["reactions"]:
    dataset.add_chromatogram(chromatogram, istd_concentration=tetralin_concentration)

# Specify the processing settings
# Default values are usually fine, but check the results and adjust if necessary
settings = ProcessingSettings(
    baseline_model="arpls",
    min_elution_time=2.5,
    max_elution_time=5,
    min_wavelength=230,
    # Some of the chromatograms contain very intense peaks, most likely decomposed reagents
    # To detect the smaller peaks of interest, disable filtering by relative height
    min_rel_prominence=0.0,
    min_prominence=1,
    # Increase the required peak purity
    explained_threshold=0.998,
)

# Process the dataset
dataset.process_all(settings, verbose=True, cores=15)

# Get concentrations relative to the internal standard
results = dataset.get_relative_concentrations()[0][
    ["Chromatogram", "starting_material", "product"]
]

# If a compound is not detected, the concentration is set to nan
# Convert nan to 0
results = results.fillna(0)

# Calculate conversion and yield
initial_concentration = 0.06
results["Conversion [%]"] = (
    100 * (initial_concentration - results["starting_material"]) / initial_concentration
)
results["Yield [%]"] = 100 * results["product"] / initial_concentration

# Print the results
print(
    results[["Chromatogram", "Conversion [%]", "Yield [%]"]]
    .round(0)
    .to_string(index=False)
)

# Plot the yields using a heatmap
# List of reagents in rows and columns
rows = "10a 10b 10c 10d 10e 10f 10g".split()
columns = "DBU/XPhos DBU/tBu-XPhos DBU/CM-Phos TMG/XPhos TMG/tBu-XPhos TMG/CM-Phos DMAP/XPhos DMAP/tBu-XPhos DMAP/CM-Phos DIPEA/XPhos DIPEA/tBu-XPhos DIPEA/CM-Phos".split()

# Extract the yields of the reaction and reshape
yields = results["Yield [%]"][
    results["Chromatogram"].apply(lambda s: s.startswith("reaction_"))
].values

yields = np.reshape(yields, [len(rows), len(columns)])

# Plot the heatmap
plt.imshow(yields, vmin=0, vmax=100, cmap="viridis")
plt.xticks(
    np.arange(len(columns)),
    labels=columns,
    rotation=45,
    ha="right",
    rotation_mode="anchor",
)
plt.yticks(np.arange(len(rows)), labels=rows)

# Add annotations
for i in range(len(rows)):
    for j in range(len(columns)):
        text = plt.text(
            j, i, f"{yields[i, j]:0.0f}", ha="center", va="center", color="w"
        )

# Show the plot
plt.title("Yields of the cyanation reaction")
plt.tight_layout()
# plt.savefig("docs/_static/ex_batch_processing.svg")
plt.show()
