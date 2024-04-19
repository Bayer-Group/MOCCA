Module **classes**
==================

This module defines light-weight classes for common data.

Data2D
~~~~~~
The Data2D class stores absorbance data at given times and wavelengths.

.. autoclass:: mocca2.classes.Data2D
    :members:

Chromatogram
~~~~~~~~~~~~
The Chromatogram class extends the Data2D class, adding more metadata.

.. autoclass:: mocca2.Chromatogram
    :inherited-members:
    :members:

Peak
~~~~
The Peak class stores information about single peak or multiple overlapping peaks.

.. autoclass:: mocca2.classes.Peak
    :members:

DeconvolvedPeak
~~~~~~~~~~~~~~~
The DeconvolvedPeak class extends the Peak class, containing information about all deconvolved components of the peak.

.. autoclass:: mocca2.classes.DeconvolvedPeak
    :inherited-members:
    :members:

Component
~~~~~~~~~
The Component class stores information about a single peak component, usually a pure compound.

.. autoclass:: mocca2.classes.Component
    :members:

Compound
~~~~~~~~
The Compound class stores information about a chemical compound.

.. autoclass:: mocca2.classes.Compound
    :members: