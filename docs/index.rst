lfast-etc
=========

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

.. grid:: 1 2 2 3
   :gutter: 2

   .. grid-item-card:: Install
      :link: installation
      :link-type: doc

      pip install, optional extras, checking the installation.

   .. grid-item-card:: Quick start
      :link: quickstart
      :link-type: doc

      Ten minutes from import to a simulated spectrum.

   .. grid-item-card:: Concepts
      :link: concepts
      :link-type: doc

      Resolution elements, SNR per resel versus per pixel, how exposure time and limiting magnitude are solved.

   .. grid-item-card:: Instruments and presets
      :link: instruments
      :link-type: doc

      The 15 Bench modes and three CSHELL modes, every preset field, how to add your own.

   .. grid-item-card:: Telescope array
      :link: telescope
      :link-type: doc

      Number of telescopes, the two fibers, the combiner model, equivalent diameter.

   .. grid-item-card:: Sources
      :link: sources
      :link-type: doc

      Blackbodies, power laws, user spectra, spectral lines, magnitudes in other bands.

   .. grid-item-card:: Sky
      :link: sky
      :link-type: doc

      The bundled SkyCalc table, what "grey time" means for your numbers, bringing your own table.

   .. grid-item-card:: Sensitivity functions
      :link: sensitivity
      :link-type: doc

      Maps over magnitude and time, limiting magnitudes, exposure versus magnitude, saturation.

   .. grid-item-card:: Plot gallery
      :link: plotting
      :link-type: doc

      Every plot helper, with the figure it makes and the one-line call behind it.

   .. grid-item-card:: Physics
      :link: physics
      :link-type: doc

      The noise model, throughput conventions and the list of known limitations.

   .. grid-item-card:: Notebooks
      :link: notebooks
      :link-type: doc

      Ten executed notebooks: a walkthrough and nine science cases versus the array size.

   .. grid-item-card:: API reference
      :link: api
      :link-type: doc

      ``LFASTETC``, ``ETCResult``, ``Source``, ``Telescope``, ``Instrument``, ``Sky``, ``plotting``.

.. toctree::
   :hidden:
   :caption: Getting started

   installation
   quickstart
   concepts

.. toctree::
   :hidden:
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
   :hidden:
   :caption: Examples

   notebooks

.. toctree::
   :hidden:
   :caption: Reference

   api
   changelog
   contributing
   review_report
