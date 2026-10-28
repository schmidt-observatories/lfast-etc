Browser ETC and JSON API
========================

.. code-block:: bash

    pip install -e ".[web]"
    lfast-etc-web --port 5000        # then open http://127.0.0.1:5000
    python -m lfast_etc.web.app --host 0.0.0.0 --port 8000

The page exposes the Python API: instrument preset with per-field overrides (blank means the
preset value, shown as the placeholder), number of telescopes, fiber and combiner, source type,
magnitude, band and system, and the calculation mode (exposure time for a target SNR, or SNR for
an exposure time split into ``n_exp`` reads). The result is a summary table (percentiles over the
band, saturated fraction, reads below the instrument minimum, zero-throughput bins, the SNR 10
limiting magnitude in one hour at the band centre) and three interactive Bokeh plots: the solved
quantity versus wavelength with saturated bins marked, the electron rates per pixel, and the
telescope throughput curves. Long grids are decimated for display and the page says so.

Every input is validated and bounded (resolving power up to 1e6, up to 1e5 telescopes, finite
numbers only), so bad input returns a message, never a 500.

JSON API
--------

``GET /api/calc?...`` takes the same field names as the form and returns the calculation:

.. code-block:: bash

    curl "http://127.0.0.1:5000/api/calc?instrument=bench_r1800&n_telescopes=20&teff=5800&mag=15&band=V&mode=exptime&snr=50"

.. code-block:: json

    {
      "mode": "exptime",
      "decimation": 1,
      "describe": "bench_r1800: R=1800, 412-686 nm, ... | 20 x 0.726 m telescopes ... | SkyCalc table ...",
      "source": "blackbody 5800 K, V=15.00 (Vega)",
      "wavelength": [412.0, 412.23, "..."],
      "exptime": [1012.3, 1008.9, "..."],
      "snr": [50.0, 50.0, "..."],
      "saturated": [false, false, "..."],
      "below_min_exptime": [false, false, "..."],
      "summary": {"exptime": {"10": 95.1, "50": 140.2, "90": 320.7}, "snr": {"10": 50.0, "50": 50.0, "90": 50.0},
                  "saturated_fraction": 0.0, "below_min_exptime_fraction": 0.0, "nan_fraction": 0.0,
                  "n_bins": 918, "limiting_mag_snr10_1hr": 19.9, "limiting_mag_band": "V", "limiting_mag_system": "Vega"}
    }

Zero-throughput bins are ``null`` (NaN is not JSON). Arrays are decimated to at most 2000
points; ``decimation`` says by how much. Invalid input returns HTTP 400 with ``{"error": "..."}``.

Fields: ``instrument``, ``resolving_power``, ``pixels_per_resel``, ``spectrograph_efficiency``,
``slit_throughput``, ``qe``, ``read_noise``, ``dark_current``, ``n_telescopes``, ``fiber``
(``18um``/``50um``), ``combine`` (``auto``/``fiber``/``independent``), ``source_type``
(``blackbody``/``flat``/``powerlaw``), ``teff``, ``alpha``, ``mag``, ``band``, ``system``,
``mode`` (``exptime``/``snr``), ``snr``, ``exptime``, ``n_exp``, ``wmin``, ``wmax``.

Deploying
---------

``lfast_etc.web.app:create_app`` is a Flask application factory, so any WSGI server works:

.. code-block:: bash

    pip install gunicorn
    gunicorn "lfast_etc.web.app:create_app()" --bind 0.0.0.0:8000 --workers 2

Each request builds its own calculator (a few milliseconds), so the app is stateless.
