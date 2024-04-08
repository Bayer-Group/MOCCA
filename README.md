# Welcome to MOCCA2

This package is based on [MOCCA](https://github.com/HaasCP/mocca) package by [HaasCP](https://github.com/HaasCP).

MOCCA2 provides powerful functions for 2D chromatogram processing. This facilitates data exploration and building automated processing pipelines.

Current features include:
 - data import using `Chromatogram()`
 - baseline correction using `estimate_baseline()` or `Chromatogram.correct_baseline()`
 - peak picking using `find_peaks()` or `Chromatogram.find_peaks()`
 - deconvolution using `deconvolve_fixed()`, `deconvolve_adaptive()` or `Chromatogram.deconvolve_peaks()`

There is also a rough automated processing pipeline. See the class `MoccaDataset`.

The docs contain concise examples that demonstrate basic features of MOCCA2.

### Documentation

The project is documented using Sphinx. As long as the repository is private, the docs cannot be published on GitHub pages, but the pages are built on the `gh-pages` branch, or you can build them yourself by running `make html` in the `docs` directory.

You can find examples in the `examples` folder.

### Tests

Currently, there is not any automated testing. The `tests` directory contains only code snippets for manual testing of some of the features. Some of the scripts in `tests` might not even work.