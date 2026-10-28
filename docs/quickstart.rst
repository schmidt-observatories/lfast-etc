Quick start
===========

This page is the ten-minute version of :doc:`notebook 01 <notebooks/01_quickstart>`. Every
figure is produced by the code above it when the documentation is built.

1. Build a calculator
---------------------

A calculator is an instrument preset, a telescope array and a sky model. Printing it tells you
what you built.

.. code-block:: python

    from lfast_etc import LFASTETC, Source, Telescope, list_instruments

    list_instruments()                 # 'bench_r1800', 14 more bench_* modes, three cshell_* modes
    etc = LFASTETC("bench_r1800", n_telescopes=20)
    print(etc)
    # LFASTETC(bench_r1800: R=1800, 412-686 nm, 2.2 pix/resel, S_opt=0.3, slit=1.0, RN=3.7 e-,
    #   dark=0.0003 e-/s/pix | 20 x 0.726 m telescopes (total area 8.28 m^2, 3.2 m equivalent),
    #   18um fiber (1.48" on sky), combine=fiber, 1 input(s) | SkyCalc table skytable_secz1.3.fits x 1, 918 bins)

The grid has one point per resolution element, :math:`\delta\lambda = \lambda/R`. The calculator
is immutable: to change the instrument, the array or the sky, build a new one.

2. Define a target
------------------

A :class:`~lfast_etc.Source` is a spectral shape plus one magnitude in one band. The magnitude
sets the absolute scale.

.. plot::

    src = Source.blackbody(5800, mag=15, band="V").add_line(656.3, depth=0.3, sigma_nm=0.05)
    print(src.describe(), "| synthetic J =", round(src.magnitude("J"), 2))
    src.plot(wavelength_nm=np.linspace(380, 1000, 2000))

The shaded box is the band the magnitude refers to. Lines are averaged over each resolution
element, so a 0.05 nm line on a 0.36 nm resel shows up as a shallower, resel-wide dip with the
same equivalent width.

3. Exposure time for a target SNR
---------------------------------

.. plot::

    etc = LFASTETC("bench_r1800", n_telescopes=20)
    src = Source.blackbody(5800, mag=15, band="V").add_line(656.3, depth=0.3, sigma_nm=0.05)
    res = etc.exptime(src, snr=50)
    res.plot()
    print(res.summary())
    print(res.at(550.0))

``exptime`` returns an :class:`~lfast_etc.ETCResult` with one value per bin. ``summary()``
gives percentiles over the band, ``at()`` one wavelength, ``to_table()`` an astropy table for
export. The time rises steeply to the blue because both the throughput and this star's SED fall
there.

4. SNR for a fixed exposure, and what limits it
-----------------------------------------------

.. plot::

    etc = LFASTETC("bench_r1800", n_telescopes=20)
    src = Source.blackbody(5800, mag=17, band="V")
    res = etc.snr(src, exptime=600)
    fig, axes = plt.subplots(3, 1, figsize=(9, 9), sharex=True)
    res.plot(ax=axes[0]); axes[0].set_xlabel("")
    res.plot_rates(ax=axes[1]); axes[1].set_xlabel("")
    res.plot_noise_budget(ax=axes[2]); axes[2].set_title("")
    fig.tight_layout()

The noise budget shows the fraction of the variance from star, sky, dark current and read noise
at every wavelength. For a V = 17 star in 600 s the star still dominates; the sky takes over at
its emission lines and, three magnitudes fainter, everywhere.

5. What the spectrum would look like
------------------------------------

.. plot::

    etc = LFASTETC("bench_r1800", n_telescopes=20)
    src = Source.blackbody(5800, mag=17, band="V").add_line(656.3, depth=0.3, sigma_nm=0.05)
    etc.snr(src, exptime=600).plot_spectrum(seed=1)

One realisation of the extracted spectrum with the per-resel SNR applied. This is the plot to
show when the question is whether a feature will be visible.

6. Change the instrument, the array, the sky
--------------------------------------------

.. code-block:: python

    from lfast_etc import Instrument, Sky

    # CSHELL, one echelle setting around the He 10830 line
    he = LFASTETC("cshell_r43000", n_telescopes=20, wavelength_min_nm=1081, wavelength_max_nm=1086)

    # override any preset field
    quiet = Instrument.from_preset("cshell_r18000", read_noise=5.0)
    LFASTETC(quiet, wavelength_min_nm=1550, wavelength_max_nm=1554)

    # a different array: 40 telescopes with the 50 um fibers (one spectrum per telescope)
    LFASTETC("bench_r1800", telescope=Telescope(n_telescopes=40, fiber="50um"))

    # a crude dark-time sky
    LFASTETC("bench_r1800", sky=Sky(scale=0.4))

7. Sensitivity maps and limiting magnitudes
-------------------------------------------

.. plot::

    etc = LFASTETC("bench_r1800", n_telescopes=20)
    gstar = Source.blackbody(5800, mag=0, band="V")          # the magnitude is supplied by the grid
    mags = np.linspace(8, 22, 57); exps = np.logspace(0, 4, 49)
    snr = etc.snr_grid(gstar, mags, exps, wavelength=(500, 650))
    plot_snr_map(mags, exps, snr, band="V", title="Bench R~1800, 20 telescopes, median SNR/resel 500-650 nm")
    print("V limit, SNR 10 in 15 min:", round(etc.limiting_magnitude(gstar, 900, 10, wavelength=(500, 650)), 2))

Because the stellar rate scales exactly as :math:`10^{-0.4 m}`, these are analytic in magnitude
and fast. :doc:`sensitivity` describes every function; :doc:`plotting` every plot.

Where next
----------

* :doc:`concepts` explains the per-resel convention and the SNR equation.
* :doc:`notebooks` has the full walkthrough and nine science cases.
* :doc:`physics` lists what the model does not include.
