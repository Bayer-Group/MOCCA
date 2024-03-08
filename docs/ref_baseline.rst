.. _ref_baseline:

Module **baseline**
===================

The baseline can be estimated using :func:`estimate_baseline()<mocca2.estimate_baseline>`.

Currently available algorithms are AsLS (asymetric least squares with smoothness penalty), arPLS (asymmetrically reweighted penalized least squares) and FlatFit. The AsLS and arPLS algorithms were adapted from `10.1039/C4AN01061B <https://doi.org/10.1039/C4AN01061B>`_.

As described in the paper, AsLS is biased to estimating baseline lower than it actually is. On the other hand, arPLS can estimate the baseline too high which can cut off some smallest peaks.

In general, I would recommend using FlatFit. If your data is very noisy, arPLS should be a good choice.

estimate_baseline()
~~~~~~~~~~~~~~~~~~~

.. autofunction:: mocca2.estimate_baseline

asls()
~~~~~~

.. autofunction:: mocca2.baseline.asls.asls

arpls()
~~~~~~~

.. autofunction:: mocca2.baseline.arpls.arpls


flatfit()
~~~~~~~~~

.. autofunction:: mocca2.baseline.flatfit.flatfit