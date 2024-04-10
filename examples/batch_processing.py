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
    min_elution_time=2.5,
    max_elution_time=4,
    min_wavelength=230,
    min_rel_prominence=0.05,
    explained_threshold=0.998,
)

# Process the dataset
dataset.process_all(settings, verbose=True, cores=14)

# Get concentrations relative to the internal standard
print(dataset.get_concentrations()[0].to_string(index=False))

results = dataset.get_relative_concentrations()[0][
    ["Chromatogram", "tetralin", "starting_material", "product"]
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
print(results[["Chromatogram", "Conversion [%]", "Yield [%]"]].to_string(index=False))

# for chrom in [dataset.chromatograms[i] for i in range(len(dataset.chromatograms))]:
#     chrom.plot()
#     plt.title(chrom.name)
#     plt.show()
# dataset.chromatograms[0].plot()
# plt.show(block=False)
# dataset.chromatograms[2].plot()
# plt.show(block=False)
# dataset.chromatograms[4].plot()
# plt.show(block=True)

rows = "10a 10b 10c 10d 10e 10f 10g".split()
columns = "DBU/XPhos DBU/tBu-XPhos DBU/CM-Phos TMG/XPhos TMG/tBu-XPhos TMG/CM-Phos DMAP/XPhos DMAP/tBu-XPhos DMAP/CM-Phos DIPEA/XPhos DIPEA/tBu-XPhos DIPEA/CM-Phos".split()

yields = results["Yield [%]"][
    results["Chromatogram"].apply(lambda s: s.startswith("reaction_"))
].values

yields = np.pad(yields, (0, len(rows) * len(columns) - len(yields)), constant_values=-1)

yields = np.reshape(yields, [len(rows), len(columns)])


plt.imshow(yields, vmin=0, vmax=100, cmap="viridis")
plt.xticks(
    np.arange(len(columns)),
    labels=columns,
    rotation=45,
    ha="right",
    rotation_mode="anchor",
)
plt.yticks(np.arange(len(rows)), labels=rows)

# Loop over data dimensions and create text annotations.
for i in range(len(rows)):
    for j in range(len(columns)):
        text = plt.text(
            j, i, f"{yields[i, j]:0.0f}", ha="center", va="center", color="w"
        )

plt.title("Yields of the cyanation reaction")
plt.colorbar(label="Yield [%]")
plt.tight_layout()
plt.show()
