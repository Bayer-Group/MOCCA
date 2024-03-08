Peak Picking
============

This example covers the peak picking functionality provided by MOCCA2. The function :func:`pick_peaks() <mocca2.pick_peaks>` can be called directly, or using the :meth:`Chromatogram.find_peaks() <mocca2.Chromatogram.find_peaks>` wrapper.

First, let's import libraries and load the data.

.. code-block:: Python 

    from mocca2 import Chromatogram, find_peaks
    from matplotlib import pyplot as plt

    # Load data
    chromatogram = Chromatogram('tests/test_data/tripeak.arw', 'tests/test_data/tripeak_blank.arw')
    chromatogram.correct_baseline()

    # Plot the chromatogram
    plt.figure()
    plt.plot(chromatogram.time, chromatogram.contract())
    plt.xlabel('Time [min]')
    plt.ylabel('Average absorbance [mAU]')
    plt.xlim(chromatogram.time[0], chromatogram.time[-1])

.. image:: !!! TO DO !!!

The region until 2 minutes contains mostly signal from solvent, and it can be excluded from the analysis.

.. code-block:: Python

    # Cropping unrelevant parts of the chromatogram
    chromatogram.extract_time(2, 10, inplace=True)
    chromatogram.extract_wavelength(220, 400, inplace=True)

.. image:: !!! TODO !!!

The current peak picking method supports only 1D data, so the chromatogram must be first contracted. The peak picking is then straight-forward.

.. code-block:: Python 

    # Pick peaks
    contracted = chromatogram.contract()

    peaks = find_peaks(
        contracted,
        min_rel_height=0.03,
        expand_borders=True,
        merge_overlapping=True,
    )

The :class:`Peak<mocca2.classes.peak.Peak>` class contains basic information about the peaks, such as left and right borders and location of maxima.

By default, the peak picking method tries to find suitable peak borders and merges overlapping peaks, so that they can be deconvolved later.

Let's see the picked peaks.

.. code-block:: Python

    # Plot chromatogram with the peaks
    plt.figure()
    plt.plot(chromatogram.time, contracted)
    plt.xlabel('Time [min]')
    plt.ylabel('Average absorbance [mAU]')
    plt.xlim(chromatogram.time[0], chromatogram.time[-1])

    t = chromatogram.time
    for peak in peaks:
        plt.hlines(-10, t[peak.left], t[peak.right], colors='red')
        for max in peak.all_maxima:
            plt.vlines(t[max], -10, contracted[max], colors='green', lw=1.)

.. image:: !!! TO DO !!!

Maybe we want to pick the smaller peaks too, this can be done by adjusting the threshold.

.. code-block:: Python

    peaks = find_peaks(
        contracted,
        min_rel_height=0.01,
        expand_borders=True,
        merge_overlapping=True,
    )

.. image:: !!! TO DO !!!

It is often rather subjective which peaks to pick, which peaks overlap, and where are the peak borders. For this reason, it might be neccessary to tweak the parameters to your needs.