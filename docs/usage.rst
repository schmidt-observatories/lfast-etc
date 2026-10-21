Install and quick start
=======================

Installation
------------

.. code-block:: bash

    git clone git@github.com:schmidt-observatories/lfast-etc.git
    cd lfast-etc
    pip install -e .            # core: numpy, scipy, astropy, matplotlib
    pip install -e ".[web]"     # + Flask/Bokeh browser ETC
    pip install -e ".[docs]"    # + Sphinx to build these pages

The throughput curves (Choi et al. 2026), the SkyCalc sky table and the instrument presets ship
inside the package; nothing needs to be downloaded.

Three objects
-------------

**Instrument** -- a spectrograph configuration. Bundled presets:

.. code-block:: python

    from lfast_etc import list_instruments, Instrument
    list_instruments()      # ['bench_r1800', 'cshell_r18000', 'cshell_r21500', 'cshell_r43000']
    ins = Instrument.from_preset("cshell_r43000", slit_throughput=0.6)   # override any field

Presets are TOML files in ``src/lfast_etc/data/config``; copy one to add a new configuration.
The CSHELL slit-loss values and the detector full wells are placeholders flagged in the files.

**Telescope** -- the array. ``n_telescopes`` is a free parameter.

.. code-block:: python

    from lfast_etc import Telescope
    Telescope(n_telescopes=20, fiber="18um")                     # combiner: all fibers -> one input
    Telescope(n_telescopes=20, fiber="50um")                     # no combiner: one spectrum per telescope
    Telescope(n_telescopes=40, fiber="18um", combine="independent")

**Source** -- an SED shape and one magnitude.

.. code-block:: python

    from lfast_etc import Source
    Source.blackbody(5800, mag=15, band="V")                     # Vega magnitudes, U..K box bands
    Source.flat(mag=20, band="R", system="AB")
    Source.powerlaw(alpha=-0.7, mag=18, band="V", system="AB")   # F_nu ~ nu^alpha
    Source.from_arrays(wavelength_nm, flux, flux_unit="flam", mag=16, band="I")
    Source.blackbody(5800, 15).add_line(656.3, depth=0.3, sigma_nm=0.05)

Calculations
------------

.. code-block:: python

    from lfast_etc import LFASTETC
    etc = LFASTETC("bench_r1800", n_telescopes=20)               # or telescope=Telescope(...)
    res = etc.exptime(star, snr=50)             # total exposure time per resel
    res = etc.snr(star, exptime=600, n_exp=3)   # SNR for 600 s split into 3 reads
    res.wavelength, res.snr, res.exptime, res.n_star, res.n_sky_pix, res.saturated
    res.at(550.0); res.summary(); res.plot(); res.plot_rates()

    # sensitivity maps (analytic in magnitude, so fast)
    etc.snr_grid(star, mags, exptimes, wavelength=(500, 650))    # (len(mags), len(exptimes))
    etc.limiting_magnitude(star, exptime=900, snr=10, wavelength=(500, 650))
    etc.exptime_vs_mag(star, mags, snr=30)

``wavelength`` selects how the per-bin arrays are reduced: ``None`` = median over the grid, a
float = nearest bin, ``(lo, hi)`` = median over a window, ``"all"`` = the per-bin array (for
``limiting_magnitude``).

CSHELL's free spectral range is ~2.5e-3 lambda (one setting is ~40 A wide at 1.55 um), so pass
``wavelength_min_nm`` / ``wavelength_max_nm`` to ``LFASTETC`` for the setting you care about.

Browser ETC
-----------

.. code-block:: bash

    pip install -e ".[web]"
    lfast-etc-web --port 5000        # then open http://127.0.0.1:5000

The page exposes the same inputs as the Python API (instrument preset with overrides, number of
telescopes, fiber, source, SNR or exposure time) and shows interactive Bokeh plots of the exposure
time or SNR per resel, the electron rates and the throughput curves. ``GET /api/calc?...`` returns
the same calculation as JSON.
