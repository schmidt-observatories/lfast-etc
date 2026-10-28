Sensitivity functions
=====================

The two basic calls return a full :class:`~lfast_etc.ETCResult`:

.. code-block:: python

    res = etc.exptime(source, snr=50, n_exp=1)      # total exposure time [s] per resel to reach the SNR
    res = etc.snr(source, exptime=600, n_exp=3)     # SNR per resel for 600 s split into 3 reads

    res.wavelength, res.delta_lambda                # the grid [nm]
    res.snr, res.exptime                            # one of them solved, the other the input
    res.n_star, res.n_sky_pix, res.dark_pix, res.read_noise, res.n_pix, res.n_exp
    res.peak_e_per_pixel, res.saturated, res.below_min_exptime
    res.at(550.0)                                   # dict at one wavelength (interpolated across NaN bins)
    res.summary()                                   # percentiles, saturated fraction, number of bins
    res.to_table()                                  # astropy Table of every per-bin column
    res.simulate(seed=0)                            # one noisy realisation [e- per resel]

The functions below exploit :math:`N_* \propto 10^{-0.4 m}` and are analytic in magnitude. They
use only the *shape* of the source; its magnitude is ignored. All take ``n_exp`` and the
``wavelength`` selector described in :doc:`concepts` (``None`` median, float nearest bin,
``(lo, hi)`` window median, ``"all"`` per bin).

snr_grid
--------

.. code-block:: python

    mags = np.linspace(8, 22, 57); exps = np.logspace(0, 4, 49)
    snr = etc.snr_grid(source, mags, exps, wavelength=(500, 650))      # shape (len(mags), len(exps))
    plot_snr_map(mags, exps, snr, band="V")

The bins are selected before the (magnitude, time, bin) cube is formed, so a 50 x 50 grid over a
50 nm window costs a few megabytes; ``wavelength="all"`` on the full CSHELL range does not.

limiting_magnitude
------------------

.. code-block:: python

    etc.limiting_magnitude(source, exptime=900, snr=10, wavelength=(500, 650))   # one number
    etc.limiting_magnitude(source, exptime=900, snr=10, wavelength="all")        # per bin
    plot_limiting_magnitude(etc, source, exptimes=[60, 600, 3600], snr=10)

The magnitude in ``source.band`` at which the SNR per resel equals ``snr`` after ``exptime``
seconds. Solved per bin from the quadratic in :math:`10^{-0.4 m}`, then reduced.

exptime_vs_mag
--------------

.. code-block:: python

    t = etc.exptime_vs_mag(source, mags, snr=30, wavelength=(500, 560))   # seconds, shape (len(mags),)
    targets_per_night = 8 * 3600 / (t + 120)                              # with 2 min overhead

Survey speed in one line (notebooks 03 and 08).

reads_needed
------------

.. code-block:: python

    n = etc.reads_needed(source, exptime=300, margin=0.8)
    res = etc.snr(source, 300, n_exp=n)

The smallest number of reads that keeps the brightest pixel below ``margin * full_well``
(``0 < margin <= 1``) in a total exposure; 1 when the instrument's full well is unknown. The
brightest pixel is approximated as a uniform share of the resel on the most heavily loaded
detector, so the estimate is optimistic (see :doc:`physics`). Reads shorter than
``Instrument.min_exptime`` are not prevented: check ``ETCResult.below_min_exptime``. Notebook 05
uses it to split bright He 10830 hosts into reads.

Array scaling
-------------

The number of telescopes is a parameter of :class:`~lfast_etc.Telescope`, so a scaling study is
a loop over calculators:

.. plot::

    gstar = Source.blackbody(5800, 0, "V")
    plot_vs_n_telescopes(lambda n: LFASTETC("bench_r1800", n_telescopes=n), gstar,
                         ns=[1, 2, 5, 10, 20, 50, 100, 200], exptimes=[60, 600, 3600], snr=10,
                         wavelength=(500, 650))

``plot_vs_n_telescopes`` takes a function that builds the calculator for a given :math:`N`, so
you decide fiber, combiner and bandpass; the top axis is the equivalent single-telescope
diameter. Notebooks 06 to 10 use this pattern for 20, 200 and 1720 telescopes.

Rebinning to a lower resolution
-------------------------------

Co-adding :math:`k = R/R_{bin}` resels raises the SNR by exactly :math:`\sqrt{k}` in this noise
model. To ask "what exposure gives SNR 5 per R = 300 bin", ask for SNR
:math:`5/\sqrt{1800/300} = 2.04` per resel:

.. code-block:: python

    k = etc.instrument.resolving_power / 300
    t = etc.exptime_vs_mag(source, mags, snr=5 / np.sqrt(k), wavelength=(450, 650))
