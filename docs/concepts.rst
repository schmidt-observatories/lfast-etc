Concepts
========

Resolution elements
-------------------

Everything in ``lfast_etc`` is computed per **resolution element** (resel). The wavelength grid
has one point per resel, :math:`\lambda_{j+1} = \lambda_j (1 + 1/R)`, and each bin is one resel
wide, :math:`\delta\lambda_j = \lambda_j/R`. A resel is the image of the fiber face,
``pixels_per_resel`` wide along the dispersion direction (2.2 for the Bench, 6 / 5 / 2.5 for the
CSHELL modes) and ``fiber_image_px`` high (circular by default), :math:`n_{pix}` pixels in all.

Three quantities carry units you should keep straight:

* :math:`N_*` [e\ :sup:`-`/s per resel]: stellar electrons in the whole resel, summed over the
  array. Not divided by :math:`n_{pix}`.
* :math:`N_{sky}` [e\ :sup:`-`/s per pixel]: sky electrons per pixel.
* :math:`D` [e\ :sup:`-`/s per pixel] and :math:`RN` [e\ :sup:`-` per pixel per read].

Per resel versus per pixel
--------------------------

Survey papers usually quote SNR per pixel. In the photon-limited regime

.. math:: \mathrm{SNR}_{pixel} \approx \mathrm{SNR}_{resel} / \sqrt{n_{pix}},

so SNR 30 per resel is about SNR 20 per pixel on the Bench. Convert before comparing with
LAMOST, SEGUE or SDSS numbers. Rebinning :math:`k` adjacent resels of equal SNR to a lower
resolution :math:`R/k` raises the SNR by exactly :math:`\sqrt{k}` in this noise model, because
every variance term scales with the number of pixels (notebooks 04 and 07 use this for
classification at R = 300).

The SNR equation
----------------

For a total exposure :math:`t` split into :math:`n_{exp}` reads (Howell 2006):

.. math:: \mathrm{SNR} = \frac{N_* t}{\sqrt{N_* t + n_{pix}\,(N_{sky} + D)\,t + n_{pix}\,n_{exp}\,RN^2}}.

``exptime`` inverts this quadratic for :math:`t`; ``snr`` evaluates it. ``n_exp`` means *the
total time is split into this many reads*, not "this many exposures of this length": the read
noise is paid ``n_exp`` times and the per-read time is ``exptime / n_exp``. Reads shorter than
the instrument's ``min_exptime`` are flagged in ``ETCResult.below_min_exptime``.

Where the source throughput is zero (telluric cores in the bundled :math:`K_{opt}` curve:
any sample below ``KOPT_EPS = 1e-6``, e.g. 687, 761, 1129 nm and the 1336 to 1424 nm water
band, plus the interpolation intervals next to those samples) the bin is invalid:
``ETCResult.valid`` is False and ``snr``/``exptime`` are NaN rather than an enormous number.
Medians and percentiles ignore those bins; ``at()`` never interpolates across them; the web
API renders them as ``null``.

Method names and result attributes
----------------------------------

``etc.exptime(source, snr=50)`` returns an :class:`~lfast_etc.ETCResult` whose ``.exptime``
array is the answer and whose ``.snr`` array is the constant target. ``etc.snr(source,
exptime=600)`` returns a result whose ``.snr`` is the answer and ``.exptime`` the constant input.
``result.solved_for`` says which one was solved, and ``result.plot()`` plots that one by default.

Limiting magnitude and sensitivity maps
---------------------------------------

Because :math:`N_* \propto 10^{-0.4 m}`, the SNR at fixed :math:`(t, n_{exp})` is a closed-form
function of magnitude. :meth:`~lfast_etc.LFASTETC.limiting_magnitude` solves the quadratic for
the magnitude that reaches a target SNR, :meth:`~lfast_etc.LFASTETC.snr_grid` evaluates the SNR on
a (magnitude, exposure) grid and :meth:`~lfast_etc.LFASTETC.exptime_vs_mag` the exposure time
along a magnitude axis. All three ignore ``source.mag`` (they use the shape only) and all three
take the same ``wavelength`` selector:

================  ==================================================  ==========================
``wavelength``    meaning                                             shape of the result
================  ==================================================  ==========================
``None``          median over the whole grid                          one number per grid point
float             the bin nearest that wavelength (must be on grid)   one number per grid point
``(lo, hi)``      median over the window                              one number per grid point
``"all"``         no reduction                                        an extra wavelength axis
================  ==================================================  ==========================

Telescope combination
---------------------

Two ways the unit telescopes are combined, chosen by the fiber (see :doc:`telescope`):

* ``combine="fiber"`` (18 um fibers through the 20-to-1 combiner): :math:`N_*` and
  :math:`N_{sky}` scale with :math:`N`, detector noise is paid once per combiner input. The gain
  is better than :math:`\sqrt{N}` whenever detector noise matters.
* ``combine="independent"`` (50 um fibers, one spectrum per telescope, co-added afterwards):
  :math:`\mathrm{SNR} = \sqrt{N}\,\mathrm{SNR}_1` exactly.

Immutability
------------

A :class:`~lfast_etc.LFASTETC` precomputes the throughputs and the sky rate when it is built.
Do not assign to ``etc.instrument``, ``etc.telescope`` or ``etc.sky`` afterwards; build a new
calculator. Building one costs a few milliseconds, so loops over array sizes or presets are
cheap.
