Basic Data Manipulation
=======================

This examples demonstrates how to load and manipulate chromatogram data using the :class:`Data2D<mocca2.classes.data2d.Data2D>` and :class:`Chromatogram<mocca2.Chromatogram>` classes.

First, let's import some essential packages:

.. code-block:: Python

    from mocca2 import Chromatogram
    import numpy as np
    from matplotlib import pyplot as plt

The data can be loaded simply using the :class:`Chromatogram<mocca2.Chromatogram>` constructor.

.. code-block:: Python

    # Load raw data from file
    chromatogram = Chromatogram('tests/test_data/tripeak.arw')

The :class:`Data2D<mocca2.classes.data2d.Data2D>` class provides basic functions for cropping, averaging, etc. Let's crop data to wavelengths between 220 - 400 nm.


.. code-block:: Python

    # Crop wavelength between 220 and 400 nm
    chromatogram.extract_wavelength(220, 400, inplace=True)

The :class:`Data2D<mocca2.classes.data2d.Data2D>` stores all data in numpy ndarrays, so plotting is straigh-forward.

.. code-block:: Python

    # Plot chromatogram (averaged over wavelenghts)
    plt.figure()
    plt.plot(chromatogram.time, chromatogram.contract())
    plt.xlabel('Time [min]')
    plt.ylabel('Average absorbance [mAU]')
    plt.xlim(chromatogram.time[0], chromatogram.time[-1])

.. image:: _static/ex_basic_chromatogram.svg

.. code-block:: Python

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

.. image:: _static/ex_basic_chromatogram_2d.svg

Let's take a closer look on the peak around 2.3 minutes. First, we extract the relevant region.

.. code-block:: Python

    # Extract raegion with peak
    peak_region = chromatogram.extract_time(2.2, 2.35)

Now, we could for example average the absorption spectrum over the peak.

.. code-block:: Python

    # Average spectrum
    # The chromatogram data are just numpy arrays, axes are [wavelength, time]
    avg_spectrum = np.mean(peak_region.data, axis=1)
    plt.figure()
    plt.plot(peak_region.wavelength, avg_spectrum)
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Average absorbance [mAU]')

.. image:: _static/ex_basic_spectrum.svg

For futher details see the reference for :class:`Data2D<mocca2.classes.data2d.Data2D>` and :class:`Chromatogram<mocca2.Chromatogram>`.