Plot gallery
============

Every helper takes an optional ``ax`` and returns it, so the plots compose with
``plt.subplots``. The result, calculator and source objects expose them as methods; the
functions live in :mod:`lfast_etc.plotting`. The figures below are generated when the
documentation is built.

Result plots
------------

``res.plot()`` draws the solved quantity (exposure time or SNR) versus wavelength; ``what=``
selects another column.

.. plot::

    etc = LFASTETC("bench_r1800")
    src = Source.blackbody(5800, 15).add_line(656.3, 0.3, 0.05)
    res = etc.exptime(src, snr=50)
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.4))
    res.plot(ax=axes[0]); res.plot(ax=axes[1], what="n_star")
    fig.tight_layout()

``res.plot_rates()``: star, sky and dark electrons per pixel per second.
``res.plot_noise_budget()``: fraction of the variance from each term, for the exposure stored
in the result. The second answers "what limits me here"; the first tells you by how much.

.. plot::

    res = LFASTETC("bench_r1800").snr(Source.blackbody(3500, 19, "V"), 900)
    fig, axes = plt.subplots(2, 1, figsize=(9, 6.5), sharex=True)
    res.plot_rates(ax=axes[0]); axes[0].set_xlabel("")
    res.plot_noise_budget(ax=axes[1]); axes[1].set_title("")
    fig.tight_layout()

``res.plot_spectrum(seed=...)``: one noisy realisation of the extracted spectrum
(``res.simulate`` returns the array). Use it to judge whether a feature will be visible.

.. plot::

    etc = LFASTETC("bench_316_63p4_o9")          # the R~23,000 H-alpha mode
    src = Source.blackbody(3300, 11, "V").add_line(656.3, -2.0, 0.03)
    res = etc.snr(src, 60)
    ax = res.plot_spectrum(seed=2); ax.set_xlim(654, 659)

Calculator plots
----------------

``etc.plot_throughput()``: :math:`K_{opt}`, :math:`\eta_{fib}`, spectrograph x QE and the end-to-end
product. ``etc.plot_sky()``: sky surface brightness per resel and the sky rate per pixel.

.. plot::

    etc = LFASTETC("cshell_r18000")
    fig, axes = plt.subplots(2, 1, figsize=(9, 6.5), sharex=True)
    etc.plot_throughput(ax=axes[0]); axes[0].set_xlabel("")
    etc.plot_sky(ax=axes[1])
    fig.tight_layout()

Source plot
-----------

``source.plot()``: the SED with the normalising band shaded and the lines marked;
``photon=True`` for photon flux.

.. plot::

    fig, ax = plt.subplots(figsize=(9, 3.6))
    for teff in (3500, 5800, 12000):
        Source.blackbody(teff, 15, "V").plot(wavelength_nm=np.linspace(380, 1800, 2000), ax=ax)
    ax.set_title("three blackbodies at V = 15", fontsize=9); ax.legend_.remove()

Detector profile and seeing
---------------------------

``ins.plot_profile()``: the fiber image of one resel on the pixel grid, with the peak-pixel
fraction used by the saturation flag (:ref:`aperture`). ``tel.plot_fiber_coupling()``: the
fraction of a point source inside the fiber versus seeing for the Moffat and Gaussian discs,
with the reference and current seeing marked (:ref:`seeing`).

.. plot::

    from lfast_etc import Instrument, Telescope
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), gridspec_kw=dict(width_ratios=[1, 1.4]))
    Instrument.from_preset("cshell_r43000").plot_profile(ax=axes[0])
    Telescope(seeing_arcsec=1.5).plot_fiber_coupling(ax=axes[1])
    fig.tight_layout()

Throughput cascade
------------------

``etc.plot_throughput_cascade()``: cumulative throughput after each stage of the Choi et al.
(2026) budget and the instrument (:ref:`budget`); ``path="sky"`` drops the fiber injection.

.. plot::

    from lfast_etc import LFASTETC, Telescope
    etc = LFASTETC("bench_r1800", telescope=Telescope(seeing_arcsec=1.5))
    etc.plot_throughput_cascade()

Sensitivity plots
-----------------

``plot_snr_map(mags, exptimes, snr)``: colour map with labelled SNR contours, from
:meth:`~lfast_etc.LFASTETC.snr_grid`. ``plot_limiting_magnitude(etc, source, exptimes, snr)``:
limiting magnitude versus wavelength. ``plot_vs_n_telescopes(make_etc, source, ns, exptimes,
snr)``: limiting magnitude versus array size with an equivalent-diameter axis.

.. plot::

    cshell = LFASTETC("cshell_r18000", wavelength_min_nm=1550, wavelength_max_nm=1560)
    kstar = Source.blackbody(4500, 0, "H")
    mags = np.linspace(2, 16, 57); exps = np.logspace(1, 4, 40)
    plot_snr_map(mags, exps, cshell.snr_grid(kstar, mags, exps), band="H", levels=(3, 10, 30, 100, 300),
                 title="CSHELL R=18k, 20 telescopes, 1.55-1.56 um")

.. plot::

    etc = LFASTETC("cshell_r18000")
    ax = plot_limiting_magnitude(etc, Source.blackbody(4500, 0, "H"), exptimes=[60, 600, 3600], snr=10)
    ax.set_title("CSHELL R=18k, SNR=10 per resel (gaps: K_opt = 0 in the 1.36-1.42 um water band)", fontsize=9)

Composing
---------

All helpers accept ``ax``, so a one-page summary of an observation is a few lines:

.. plot::

    etc = LFASTETC("bench_r1800")
    src = Source.blackbody(5800, 17, "V").add_line(656.3, 0.3, 0.05)
    res = etc.snr(src, 600)
    fig, axes = plt.subplots(2, 2, figsize=(12, 6.5))
    src.plot(ax=axes[0, 0], wavelength_nm=np.linspace(380, 1000, 1000))
    res.plot(ax=axes[0, 1])
    res.plot_noise_budget(ax=axes[1, 0]); axes[1, 0].set_title("")
    res.plot_spectrum(ax=axes[1, 1]); axes[1, 1].set_title("")
    fig.tight_layout()

Styling
-------

Importing ``lfast_etc.plotting`` (which every ``.plot*`` method does) applies the package style,
``lfast_etc/data/lfast.mplstyle``: STIX serif fonts, inward ticks with minor ticks, a light grid
and a six-colour palette. The helpers themselves set labels and titles only, so
``plt.style.use('default')`` after the import restores matplotlib's look, and any other style
sheet can be layered on top.
