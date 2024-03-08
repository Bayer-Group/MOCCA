Baseline Correction
===================

Baseline correction is often overlooked, but it is crucial for accurate peak modelling and integration.

First let's import MOCCA2.

.. code-block:: Python

    from mocca2 import Chromatogram, estimate_baseline
    from matplotlib import pyplot as plt

The best way to correct baseline is using blank, and then by refining the baseline. This is a two-liner:

.. code-block:: Python

    # Load sample chromatogram and substract blank
    chromatogram = Chromatogram(
        'tests/test_data/tripeak.arw',
        blank='tests/test_data/tripeak_blank.arw'
    )

    # Refine the baseline
    chromatogram.correct_baseline()

    # Plotting the chromatogram with corrected baseline
    plt.figure()
    plt.plot(chromatogram.time, chromatogram.contract())
    plt.xlabel('Time [min]')
    plt.ylabel('Average absorbance [mAU]')
    plt.xlim(chromatogram.time[0], chromatogram.time[-1])

.. image:: _static/ex_baseline_corrected.svg

We can also take a look onto different algorightms for baseline estimation.
Let's pretend we don't have the blank run, so that we can compare it to the estimated baseline.

.. code-block:: Python

    # Load and average data
    chromatogram = Chromatogram('tests/test_data/tripeak.arw')
    averaged = chromatogram.contract()
    averaged_blank = Chromatogram('tests/test_data/tripeak_blank.arw').contract()

    # Estimate baseline using different methods
    baseline_asls = estimate_baseline(averaged, method='asls')
    baseline_arpls = estimate_baseline(averaged, method='arpls')
    baseline_flatfit = estimate_baseline(averaged, method='flatfit')

    # Plot the result
    plt.figure()
    plt.plot(chromatogram.time, averaged, label='Raw data', lw=1.)
    plt.plot(chromatogram.time, averaged_blank, label='Blank', lw=1.)
    plt.plot(chromatogram.time, baseline_asls, label='ASLS baseline', lw=1.)
    plt.plot(chromatogram.time, baseline_arpls, label='arPLS baseline', lw=1.)
    plt.plot(chromatogram.time, baseline_flatfit, label='FlatFit baseline', lw=1.)
    plt.xlabel('Time [min]')
    plt.ylabel('Average absorbance [mAU]')
    plt.xlim(chromatogram.time[0], chromatogram.time[-1])
    plt.ylim(-30, 100)
    plt.legend()

.. image:: _static/ex_baseline_comparison.svg

A few remarks:
 * The blank has slightly lower intensity than the sample baseline - for this reason it is better to refine baseline even if blank is substracted
 * AsLS significanly underestimates baseline in the 0.5 - 1.5 minute region
 * It is hard to interpret the 0 - 0.5 min region where are positive and negative peaks caused by the solvents. These obviously cannot be corrected by any general baseline correction algorithm
 * The region under the peak around 1.4 - 1.8 minutes shows how different methods approach baseline under peaks differently. The FlatFit seems to best capture the bend in the baseline.

For not-very-noisy data, I would recommend using FlatFit. It is very fast and stable algorithm which can also handle negative peaks.

The description of individual methods is described in the :ref:`baseline <ref_baseline>` reference.