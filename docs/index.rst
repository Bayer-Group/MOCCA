Welcome to MOCCA2!
==================

MOCCA2 is a Python package for automatic processing of HPLC chromatograms.

To automate your workflow and get accurate results, MOCCA2 features:
 - support for raw data files from Agilent, Shimadzu and Waters
 - automatic baseline correction
 - adaptive peak picking
 - automatic purity checking and peak deconvolution
 - compound tracking across chromatograms
 - fully automatic processing of any number of chromatograms

Getting Started
---------------

`Get Started <installation.html>`_ by installing MOCCA2 and processing your first chromatogram!

Publications and MOCCA
----------------------

This package is based on `MOCCA <https://github.com/HaasCP/mocca>`_ package by `HaasCP <https://github.com/HaasCP>`_. This work has been published by `Christian Haas et al. in 2023 <https://doi.org/10.1021/acscentsci.2c01042>`_.

Inspired by MOCCA, MOCCA2 features more Pythonic interface as well as adaptive and more accurate algorithms.

Publication featuring MOCCA2 is coming soon!

Contents
========

.. toctree::
   :maxdepth: 1
   :hidden:

   self
   installation

.. toctree::
   :maxdepth: 1
   :includehidden:
   :caption: Examples

   ex_basic_data
   ex_baseline
   ex_peak_picking

.. toctree::
   :maxdepth: 1
   :includehidden:
   :caption: Reference

   ref_classes
   ref_dataset
   ref_parsers
   ref_baseline
   ref_peak
   ref_deconvolution
   ref_clustering

Contributing
============

Please see the `Contributing Guide <https://github.com/oboril/mocca/blob/main/CONTRIBUTING.md>`_.