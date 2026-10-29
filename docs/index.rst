lfast-etc: Exposure Time Calculator for LFAST
=============================================

.. image:: _static/logo_schmidt_observatory_system.png
   :alt: Schmidt Observatory System
   :width: 320px

``lfast_etc`` is the exposure time calculator for the **LFAST** telescope array at Kitt Peak:
nominally 20 unit telescopes of 0.726 m, each feeding a fiber, with the fibers merged into one
spectrograph input. Two spectrographs are modelled, the **WIYN Bench** (visible, 15 grating
modes) and **CSHELL** (1.08 to 1.75 um, R up to 43,000). The number of telescopes, the fiber,
and every instrument parameter are plain inputs, so the calculator follows the array design as
it evolves.

.. code-block:: python

    from lfast_etc import LFASTETC, Source

    etc = LFASTETC("bench_r1800", n_telescopes=20)
    star = Source.blackbody(5800, mag=15, band="V")
    res = etc.exptime(star, snr=50)       # total exposure time per resolution element
    res.plot(); res.plot_noise_budget(); res.summary()

.. plot::
    :include-source: False

    etc = LFASTETC("bench_r1800", n_telescopes=20)
    star = Source.blackbody(5800, mag=17, band="V").add_line(656.3, 0.3, 0.05)
    res = etc.snr(star, exptime=600)
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.4))
    res.plot_spectrum(ax=axes[0], seed=1); axes[0].set_title("simulated 600 s spectrum, V=17", fontsize=9)
    res.plot_noise_budget(ax=axes[1]); axes[1].set_title("what limits the SNR", fontsize=9)
    fig.tight_layout()

.. note::

   New here? Read :doc:`installation`, then work through :doc:`quickstart` and
   :doc:`concepts`. The :doc:`notebooks` are executed Jupyter notebooks that cover
   each science case end to end, with their figures and printed numbers saved in place.

What it models
--------------

- A telescope array: the number of unit telescopes, the fiber feeding each one and the
  combiner, giving an equivalent collecting diameter that follows the array design.
- Two spectrographs with 18 presets, the 15 WIYN Bench grating modes and three CSHELL
  modes, each a plain TOML file you can copy and edit.
- Sources: blackbodies, power laws, user spectra and added spectral lines, normalised to a
  magnitude in any band.
- The sky: the bundled SkyCalc table for Kitt Peak, radiance per resel and per pixel through
  the line-spread function, and a hook for your own table.
- Signal-to-noise per resolution element and per detector pixel, the exposure time that
  reaches a target SNR, limiting magnitudes, saturation, and maps over magnitude and time.
- Plots: simulated spectra, noise budgets, sensitivity maps and exposure-versus-magnitude
  curves, in one house style.

Contents
--------

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   installation
   quickstart
   concepts

.. toctree::
   :maxdepth: 2
   :caption: User guide

   instruments
   telescope
   sources
   sky
   sensitivity
   plotting
   web
   physics

.. toctree::
   :maxdepth: 2
   :caption: Examples

   notebooks

.. toctree::
   :maxdepth: 2
   :caption: Reference

   api
   changelog
   contributing
   review_report

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
