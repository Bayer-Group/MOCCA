[![PyPI](https://img.shields.io/pypi/v/mocca2.svg)](https://pypi.org/project/mocca2/)
![pytest](https://github.com/bayer-group/MOCCA/actions/workflows/ci.yaml/badge.svg)
[![Docs Pages](https://github.com/bayer-group/MOCCA/actions/workflows/deploy_pages.yaml/badge.svg)](https://bayer-group.github.io/MOCCA/)
[![Example Data](https://github.com/bayer-group/MOCCA/actions/workflows/package_data.yaml/badge.svg)](https://github.com/bayer-group/MOCCA/tree/example-data)

# Welcome to MOCCA2

MOCCA2 is a Python package for automatic processing of HPLC chromatograms.

To automate your workflow and get accurate results, MOCCA2 features:
 - support for raw data files from Agilent, Shimadzu and Waters
 - automatic baseline correction
 - adaptive peak picking
 - automatic purity checking and peak deconvolution
 - compound tracking across chromatograms
 - fully automatic processing of any number of chromatograms


## Documentation

Examples and detailed documentation are documented at [https://bayer-group.github.io/MOCCA](https://bayer-group.github.io/MOCCA).

## Getting Started

The latest version of MOCCA2 can be installed simply using pip:

```
pip install mocca2
```

Example data can be then downloaded using the following command:

```
python -m mocca2 --download-data
```

Now you are ready to process your first chromatogram!

```
from mocca2 import example_data
from matplotlib import pyplot as plt

# Load example data
chromatogram = example_data.example_1()

# Correct the baseline
chromatogram.correct_baseline()

# Crop the chromatogram to the region of interest, 1.4 to 1.8 minutes
chromatogram.extract_time(1.4, 1.8, inplace=True)

# Exclude low wavelengths that tend to be noisy - ignore everything below 220 nm
chromatogram.extract_wavelength(220, None, inplace=True)

# Find peaks in the chromatogram
chromatogram.find_peaks(min_height=2)

# Deconvolve the peaks
print("Deconvolving peaks, this migth take a minute...")

chromatogram.deconvolve_peaks(
    model="FraserSuzuki", min_r2=0.999, relaxe_concs=False, max_comps=5
)

print("Deconvolved!")

# Plot the chromatogram
chromatogram.plot()
plt.show()
```

## Publications and MOCCA

This package is based on [MOCCA](https://github.com/HaasCP/mocca) package by [HaasCP](https://github.com/HaasCP). This work has been published by [Christian Haas et al. in 2023](https://doi.org/10.1021/acscentsci.2c01042).

Inspired by MOCCA, MOCCA2 features more Pythonic interface as well as adaptive and more accurate algorithms.

Publication featuring MOCCA2 is coming soon!

## Repository Details

This repository automates numerous workflows:

### Automatic testing
On push to `main`, all tests in the `tests` directory are automatically run. Currently, MOCCA2 is tested on Ubuntu with Python 3.10, 3.11 and 3.12.

### Docs
On push to `main`, the Sphinx docs are automatically compiled and published to [GitHub pages](https://bayer-group.github.io/mocca).

### Example data
The repository contains various example datasets:
 - Knoevenagel condensation ([Christian Haas et al., 2023](https://doi.org/10.1021/acscentsci.2c01042))
 - Cyanation screening ([Christian Haas et al., 2023](https://doi.org/10.1021/acscentsci.2c01042))
 - Diterpene esters from coffee extracts ([Erny et al., 2021](https://doi.org/10.5281/zenodo.5412345))
 - and various standalone chromatograms

Since these datasets don't fit into the PyPI package size limit, they are automatically compressed and published onto `example-data` branch on push to `main`. 

The data can be automatically downloaded using ``python -m mocca2 --download-data``.

### Publishing to PyPI and GitHub
On push to `main`, the MOCCA2 package is automatically published to [PyPI](https://pypi.org/project/mocca2/) and [GitHub Releases](https://github.com/bayer-group/MOCCA/releases).

## Contributing

The process for contributing is outlined in [CONTRIBUTING.md](https://github.com/bayer-group/MOCCA/blob/main/CONTRIBUTING.md).