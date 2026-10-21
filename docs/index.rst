lfast-etc
=========

``lfast_etc`` is an exposure time calculator for the **LFAST** telescope array at Kitt Peak:
nominally 20 x 0.726 m unit telescopes feeding one spectrograph through an 18 um fiber combiner.
Two instruments are modelled, the **WIYN Bench** spectrograph (visible) and **CSHELL**
(1.08-1.75 um, R up to 43,000). The number of telescopes, the fiber and every instrument
parameter are plain inputs, so the calculator follows the array design as it evolves.

.. code-block:: python

    from lfast_etc import LFASTETC, Source

    etc = LFASTETC("bench_r1800", n_telescopes=20)
    star = Source.blackbody(5800, mag=15, band="V")
    res = etc.exptime(star, snr=50)       # per resolution element
    res.plot(); res.summary()

.. grid:: 1 2 2 2
   :gutter: 2

   .. grid-item-card:: Install & quick start
      :link: usage
      :link-type: doc

      pip install, the three objects you need, the web ETC.

   .. grid-item-card:: Physics
      :link: physics
      :link-type: doc

      Per-resel SNR model, throughput conventions, telescope combination.

   .. grid-item-card:: Notebooks
      :link: notebooks
      :link-type: doc

      Sensitivity maps, metallicity surveys, TDAMM follow-up, exoplanet atmospheres.

   .. grid-item-card:: API
      :link: api
      :link-type: doc

      ``LFASTETC``, ``Source``, ``Telescope``, ``Instrument``, ``Sky``.

.. toctree::
   :hidden:
   :maxdepth: 2

   usage
   physics
   notebooks
   api
   review_report
