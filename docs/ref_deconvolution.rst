Module **deconvolution**
========================

This module contains numerous routines for fitting peak shapes and deconvolution of 2D spectra.

**Automatic deconvolution**
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The deconvolution process is inspired by article by Shimadzu `10.1016/j.chroma.2016.09.037 <https://doi.org/10.1016/j.chroma.2016.09.037>`_.

The entire procedure is as follows:

1. Get initial guess of spectra
    - if enough maxima are present, spectra at maxima are used
    - otherwise, the timepoints are clustered based on cosine similarity and averaged spectra of clusters are used

2. Refine spectra by fitting peak shapes to some :class:`PeakModel<mocca2.deconvolution.peak_models.PeakModel>`

3. Optionally, unrestricted refinement of concentrations based on spectra

deconvolve_adaptive()
~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mocca2.deconvolve_adaptive


deconvolve_fixed()
~~~~~~~~~~~~~~~~~~

.. autofunction:: mocca2.deconvolve_fixed

**Peak Models**
~~~~~~~~~~~~~~~

Mathematical models that describe peak shapes. A nice summary of peak shapes is in `10.1016/S0922-3487(98)80024-5 <https://doi.org/10.1016/S0922-3487(98)80024-5>`_.

You can define your own peak shape model, but it must inherit from :class:`PeakModel<mocca2.deconvolution.peak_models.PeakModel>` and implement all relevant methods.

PeakModel
~~~~~~~~~

Abstract class that defines all methods that peak model must have.

.. autoclass:: mocca2.deconvolution.peak_models.PeakModel
    :members:

FraserSuzuki
~~~~~~~~~~~~

Fraser-Suzuki peak model works very well on high-quality experimental data. It is described by the following equation.

.. math::

    f(t, h, \mu, \sigma, a) = h \cdot \exp \left( - \frac{1}{a^2} \ln^2 \left( 1 + \frac{a(t-\mu)}{\sigma} \right) \right)


.. autoclass:: mocca2.deconvolution.peak_models.FraserSuzuki

BiGaussian
~~~~~~~~~~

BiGaussian peak model is simple yet efficient for describing peaks. It is described by the following equation.

.. math::

    f(t, h, \mu, \sigma_1, \sigma_2) =
    \begin{cases} 
        h \cdot \exp \left( -\frac{(t-\mu)^2}{2{\sigma_1}^2}\right) &\text{if } t < \mu \\
        h \cdot \exp \left( -\frac{(t-\mu)^2}{2{\sigma_2}^2}\right) &\text{if } t \ge \mu
    \end{cases}


.. autoclass:: mocca2.deconvolution.peak_models.BiGaussian


BiGaussianTailing
~~~~~~~~~~~~~~~~~

BiGaussian with exponential tailing. This is suitable if the peaks have slowly decaying tails. It is described by the following equation.

.. math::

    f(t, h, \mu, \sigma_1, \sigma_2, t_h, t_w) =
    \begin{cases} 
        h \cdot \exp \left( -\frac{(t-\mu)^2}{2{\sigma_1}^2}\right) &\text{if } t < \mu \\
        h \cdot (1-t_h)\cdot \exp \left( -\frac{(t-\mu)^2}{2{\sigma_2}^2}\right) + h \cdot t_h \cdot \exp \left(-\frac{t-\mu}{t_w}\right) &\text{if } t \ge \mu
    \end{cases}


.. autoclass:: mocca2.deconvolution.peak_models.BiGaussianTailing

BEMG
~~~~

Bi-Exponentially Modified Gaussian function as described in `10.1016/j.chroma.2016.09.037 <https://doi.org/10.1016/j.chroma.2016.09.037>`_. Very versatile model that can describe tailing and fronting.


.. autoclass:: mocca2.deconvolution.peak_models.Bemg

**Fitting Peak Models**
~~~~~~~~~~~~~~~~~~~~~~~

The peak shape, and optionally also the component spectra, can be fitted using :func:`fit_peak_model()<mocca2.deconvolution.fit_peak_model.fit_peak_model>` to best explain the 2D chromatogram data.

Initial guess for spectra of individual components can be obtained from :func:`guess_spectra()<mocca2.deconvolution.guess_spectra.guess_spectra>`.


fit_peak_model_and_spectra()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mocca2.deconvolution.fit_peak_model.fit_peak_model

guess_spectra()
~~~~~~~~~~~~~~~

.. autofunction:: mocca2.deconvolution.guess_spectra.guess_spectra

**Constrained Linear Problems**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This submodule provides methods for estimating the **non-negative** concentrations or spectra that best fit the 2D chromatogram data, assuming the other information (concentrations or spectra) is known.

concentrations_from_spectra()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mocca2.deconvolution.nonnegative_lstsq.concentrations_from_spectra

spectra_from_concentrations()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mocca2.deconvolution.nonnegative_lstsq.spectra_from_concentrations