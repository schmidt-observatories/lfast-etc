Telescope array
===============

A :class:`~lfast_etc.Telescope` is :math:`N` identical unit telescopes feeding one spectrograph.

.. code-block:: python

    from lfast_etc import Telescope
    Telescope(n_telescopes=20, fiber="18um")                       # nominal: combiner, one input
    Telescope(n_telescopes=20, fiber="50um")                       # one spectrum per telescope
    Telescope(n_telescopes=200, fiber="18um")                      # 10 arrays, 10 combiner inputs
    Telescope(n_telescopes=200, fiber="18um", telescopes_per_input=None)   # idealised single combiner
    Telescope(n_telescopes=20, diameter_m=0.8, fiber_diameter_arcsec=1.6)  # what-ifs

    t = Telescope(190)
    t.total_area_m2, t.equivalent_diameter_m, t.n_inputs          # 78.6 m^2, 10.0 m, 10

Fibers
------

==========  ================  =======================  ==========================================
fiber       on-sky diameter   default ``combine``      notes
==========  ================  =======================  ==========================================
``18um``    1.48 arcsec       ``"fiber"``              goes through the N-to-1 combiner
``50um``    4.10 arcsec       ``"independent"``        no combiner; about 2x the point-source
                                                       coupling and 7.7x the sky aperture
==========  ================  =======================  ==========================================

The two throughput curves per fiber (:math:`K_{opt}` and :math:`\eta_{fib}`, Choi et al. 2026)
ship with the package and are sampled at 1 nm from 400 to 1800 nm. Asking for the 50 um fiber
with ``combine="fiber"`` is allowed as a what-if and raises a warning, since that configuration
is not planned.

The combiner model
------------------

* ``combine="fiber"``: fibers are merged in groups of ``telescopes_per_input`` (20) into one
  spectrograph input each. Inside a group, star and sky electrons scale with the group size
  while dark current and read noise are paid once. Arrays larger than one group are
  :math:`n_{in} = \lceil N/20 \rceil` inputs co-added incoherently, so beyond 20 telescopes the
  SNR grows as :math:`\sqrt{n_{in}}` on top of the single-array value.
  ``telescopes_per_input=None`` models one combiner regardless of :math:`N`.
* ``combine="independent"``: every telescope yields its own spectrum. The effective dark current
  is :math:`N D` and the effective read noise :math:`\sqrt{N}\,RN`, so
  :math:`\mathrm{SNR} = \sqrt{N}\,\mathrm{SNR}_1` exactly. ``telescopes_per_input`` has no
  meaning in this mode.

The saturation estimate always refers to one physical detector: in independent mode the light
per pixel does not grow with :math:`N`; in fiber mode it grows with the group size but not with
the number of groups.

.. plot::

    gstar = Source.blackbody(5800, mag=0, band="V")
    ns = [1, 2, 4, 8, 12, 16, 20, 30, 40, 60, 80, 100, 150, 200]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for fiber, combine, tpi, ls in (("18um", "fiber", 20, "-"), ("18um", "fiber", None, ":"),
                                    ("18um", "independent", 20, "-"), ("50um", "independent", 20, "-")):
        lim = [LFASTETC("bench_r1800", telescope=Telescope(n, fiber, combine=combine, telescopes_per_input=tpi))
               .limiting_magnitude(gstar, 900, 10, wavelength=(500, 650)) for n in ns]
        ax.plot(ns, lim, marker="o", ms=3, ls=ls, label=f"{fiber}, {combine}" + (" (ideal combiner)" if tpi is None else ""))
    ax.set_xscale("log"); ax.invert_yaxis(); ax.axvline(20, color="0.4", ls="--", lw=0.8)
    ax.set_xlabel("Number of telescopes"); ax.set_ylabel("Limiting V (SNR=10/resel, 900 s)")
    ax.legend(fontsize=8); ax.set_title("Bench R~1800, median over 500-650 nm", fontsize=9)
    equivalent_diameter_axis(ax)

Geometry
--------

``diameter_m`` is the clear aperture of one unit telescope (0.726 m). The 15.5 cm central
obscuration is already in the throughput budget, so do not subtract it again. ``area_m2`` is one
telescope, ``total_area_m2`` the array, ``equivalent_diameter_m`` the single telescope with the
same collecting area (3.2 m for 20 units, 10.3 m for 200, 30.1 m for 1720).
