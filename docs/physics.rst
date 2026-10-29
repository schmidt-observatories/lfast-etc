Physics
=======

Everything is computed per **resolution element** on a grid with one point per resel,
:math:`\lambda_{j+1} = \lambda_j (1 + 1/R)`, :math:`\delta\lambda_j = \lambda_j / R`.
This is the per-bin logic of the reference LFAST ETC notebook (2026-10-12), which the test-suite
reproduces to 1e-6.

Source
------

A :class:`~lfast_etc.Source` is a relative SED shape :math:`s(\lambda)` and one magnitude. The
band flux of the magnitude in a box band of centre :math:`\lambda_0`, width :math:`\Delta\lambda`
and Vega zero point :math:`f_{\nu,0}` (Bessell et al. 1998) is

.. math:: F_{band} = f_{\nu,0}\,10^{-0.4 m}\,\frac{c\,\Delta\lambda}{\lambda_0^2},

and :math:`F_\lambda(\lambda) = s(\lambda)\,F_{band}/\int_{band} s\,d\lambda` in
W m\ :sup:`-2` nm\ :sup:`-1`. The frequency width :math:`c\Delta\lambda/\lambda_0^2` is the
reference-notebook approximation; the exact :math:`c\Delta\lambda/(\lambda_{lo}\lambda_{hi})`
differs by ~1% (0.01 mag). Lines are multiplicative Gaussians added *after* normalisation, so
``mag`` is the continuum magnitude.

Stellar electron rate
---------------------

With :math:`N` telescopes of area :math:`A` each, point-source telescope throughput
:math:`K_{opt}(\lambda)`, spectrograph efficiency :math:`S_{opt}` (times the slit throughput for
slit modes) and detector QE,

.. math:: N_*(\lambda_j) = N\,A\,K_{opt}\,S_{opt}\,QE\;\frac{F_\lambda\,\delta\lambda_j}{hc/\lambda_j}
   \quad [\mathrm{e^-\,s^{-1}\ per\ resel}].

:math:`N_*` is **not** divided by the number of pixels: the source is one process spread across
:math:`n_{pix}` pixels by the optics.

Throughput conventions
----------------------

The bundled curves follow Choi et al. (2026, SPIE 14147). :math:`K_{opt}` is the end-to-end
*point-source* throughput of the telescope + fiber at the reference airmass; it **already
contains** the fiber injection efficiency :math:`\eta_{fib}` (the 50 um and 18 um
:math:`K_{opt}` curves differ by exactly :math:`\eta_{fib,50}/\eta_{fib,18}`) and the site
extinction. Hence

* star (point source): :math:`K_{opt}`,
* sky (fills the fiber uniformly, no injection loss): :math:`K_{opt}/\eta_{fib}`.

``apply_extinction`` is therefore ``False`` by default; set it to ``True`` only with a
:math:`K_{opt}` that excludes the atmosphere.

Sky
---

The SkyCalc table (ESO sky model at 2640 m, Paranal-like aerosols since SkyCalc has no Kitt
Peak site; airmass 1.3, PWV 2.5 mm, **moon up** 45 degrees away at 45 degrees altitude, R=20,000
sampling) gives the radiance :math:`B_{sky}` in ph s\ :sup:`-1` m\ :sup:`-2` um\ :sup:`-1`
arcsec\ :sup:`-2`. Every SkyCalc sample inside a resel is averaged (nearest sample if the resel is
narrower than the table sampling). With the fiber solid angle :math:`\Omega`,

.. math:: N_{sky}(\lambda_j) = \frac{N\,A\,\Omega\,(K_{opt}/\eta_{fib})\,S_{opt}\,QE\,
   B_{sky}\,\delta\lambda_j}{n_{pix}} \quad [\mathrm{e^-\,s^{-1}\,pixel^{-1}}].

Pass ``Sky(path, scale)`` to use another SkyCalc file or a crude moon scaling (see :doc:`sky`).

SNR and exposure time
---------------------

For a total exposure :math:`t` split into :math:`n_{exp}` reads (Howell 2006),

.. math:: \mathrm{SNR} = \frac{N_* t}{\sqrt{N_* t + n_{pix}\,(N_{sky} + D)\,t + n_{pix}\,n_{exp}\,RN^2}}.

``exptime`` inverts this quadratic for :math:`t` at a target SNR. Bins where :math:`K_{opt} < 10^{-6}`
(telluric cores, ``lfast_etc.telescope.KOPT_EPS``) are invalid and return NaN. The sky-subtraction noise term :math:`(1 + n_{pix}/n_{sky})` of
Howell is omitted, as in the reference notebook (a sky fiber is assumed well sampled).
:math:`n_{pix}` is Howell's full extraction aperture, ``pixels_per_resel`` (the dispersion width
of one resel: 2.2 pixels on the Bench, 6 / 5 / 2.5 for the CSHELL modes at 2.7 km/s per pixel)
times ``fiber_image_px`` (the spatial height of the fiber image; see :ref:`aperture`). Because :math:`N_* \propto
10^{-0.4 m}` the limiting magnitude at fixed :math:`(t, \mathrm{SNR})` is also a quadratic and is
solved analytically; ``snr_grid`` evaluates the SNR on a (magnitude, exposure) grid the same way.

Telescope combination
---------------------

* ``combine="fiber"`` (18 um fibers through the 20-to-1 combiner): :math:`N_*` and :math:`N_{sky}`
  scale with :math:`N`, dark current and read noise are paid once per combiner input,
  :math:`n_{in} = \lceil N / 20 \rceil` (``telescopes_per_input``). Beyond one array the SNR
  therefore grows as :math:`\sqrt{n_{in}}` on top of the single-array value; set
  ``telescopes_per_input=None`` for an idealised single combiner.
* ``combine="independent"`` (50 um fibers, one spectrum per telescope, co-added afterwards):
  :math:`D \to N D` and :math:`RN^2 \to N\,RN^2`, so :math:`\mathrm{SNR} = \sqrt{N}\,\mathrm{SNR}_1`
  exactly.

.. _aperture:

Extraction aperture and detector profile
----------------------------------------

A fiber-fed spectrograph images the fiber face onto the detector, so one resel is a disc
``pixels_per_resel`` wide along the dispersion direction and ``fiber_image_px`` high across it
(#13). ``fiber_image_px=None`` (the default) assumes a circular image, ``fiber_image_px =
pixels_per_resel``; the CSHELL slit modes set it to the native 6 px fiber image explicitly because
the slit narrows the dispersion direction only. ``fiber_image_px=1.0`` reproduces the reference
notebook, which counted the dispersion direction alone. Howell's :math:`n_{pix}` is the product,
so the sky, dark and read-noise terms all scale with the fiber-image height; for the Bench at
3.7 e- this is the dominant faint-end sensitivity (a factor 2.2 in :math:`n_{pix}`). Replace the
default with the measured 1x1-binned fiber-image height when the LFAST team supplies it.

The *shape* of that disc sets the brightest pixel (#19). ``Instrument.profile`` chooses it:

* ``"tophat"`` (default): an elliptical top hat of those two diameters, the image of a uniformly
  illuminated fiber face, convolved with a 0.35 px Gaussian camera blur;
* ``"uniform"``: the light spread evenly over :math:`n_{pix}` pixels (peak fraction
  :math:`1/n_{pix}`, the pre-#19 convention);
* ``"gaussian"``: a Gaussian with those FWHMs (lower peak than the top hat because of its wings).

The image is rendered on a 20x sub-pixel grid centred on a pixel centre (the worst case) and binned
to pixels; ``Instrument.peak_fraction`` is the largest pixel's share and ``plot_profile()`` draws
it. Seeing does not enter: the fiber scrambles the input, so the detector profile is the same in
any seeing. For the Bench's 2.2 x 2.2 px image the top hat puts 25% of the resel in the peak
pixel, versus 21% uniform and 17% Gaussian.

Saturation
----------

The brightest pixel is :math:`(f\,(N_* + n_{pix} N_{sky})\,p + D)\,t/n_{exp}` with :math:`p` the
``peak_fraction`` above and :math:`f` the fraction of the telescopes on the most heavily loaded
input (a full group of 20 for ``N = 21``, one telescope in independent mode), flagged against
``full_well``; ``reads_needed`` returns the smallest :math:`n_{exp}` that keeps it under a margin of
the well. The sub-pixel centring is the worst case, so the flag is slightly conservative. An
``Instrument`` without a ``full_well`` (``inf``) is never flagged and reports
``saturated_fraction = None``.

.. _seeing:

Seeing
------

The shipped :math:`K_{opt}` and :math:`\eta_{fib}` curves are for one seeing (Choi et al. 2026;
the value is not stated in the data, and ``Telescope.seeing_ref_arcsec = 1.0`` is a placeholder
until it is confirmed). ``Telescope(seeing_arcsec=s)`` rescales the point-source throughput by the
geometric fraction of a seeing disc of FWHM :math:`s` that falls inside the 1.48 arcsec fiber,
relative to the same fraction at the reference seeing (#14):

.. math:: K_{opt} \to K_{opt}\,\frac{E(s(\lambda))}{E(s_{ref}(\lambda))}, \qquad
   E_{\rm Moffat}(s) = 1 - \left[1 + (r_f/\alpha)^2\right]^{1-\beta}, \quad
   \alpha = \frac{s}{2\sqrt{2^{1/\beta} - 1}},

with :math:`r_f` the fiber radius, :math:`\beta = 3.5` by default (``moffat_beta``; 4.77 for pure
Kolmogorov turbulence, 2.5 for profiles with strong wings) and a Gaussian alternative
(``seeing_psf="gaussian"``). The seeing is scaled with wavelength as :math:`\lambda^{-1/5}`
from its 500 nm value. The sky, which fills the fiber whatever the seeing, is unchanged. The
model is geometric (no diffraction, no focal-ratio degradation), so it is a *relative* correction
and only as good as the reference seeing. ``Telescope.plot_fiber_coupling()`` draws
:math:`E(s)` for both shapes.

Known limitations
-----------------

* QE is a flat placeholder per instrument (pass a callable ``qe(wavelength_nm)`` for a curve).
* CSHELL slit losses and full wells are placeholders.
* One sky table at airmass 1.3, sampled at R=20,000; no airmass or moon-phase model beyond a
  scale factor. For the CSHELL presets (R > 20,000) each resel takes the nearest sample, so
  airglow lines are pre-smoothed to R=20,000; regenerate a finer SkyCalc table for line-level work.
* :math:`K_{opt}` is sampled at 1 nm; the telluric treatment is coarse and belongs to the
  airmass assumed by Choi et al., not necessarily 1.3.
* The fiber-image height (``fiber_image_px``) defaults to a circular image; the measured 1x1 value
  is still needed from the LFAST team, and the reference seeing of the throughput curves is a
  placeholder (see :ref:`aperture` and :ref:`seeing`).
* Reads shorter than ``Instrument.min_exptime`` (8 s Fowler minimum for CSHELL) are flagged in
  ``ETCResult.below_min_exptime``; the read noise is not adjusted.
* No telluric absorption in the *source* path unless ``apply_extinction=True`` (and then the
  extinction already in :math:`K_{opt}` is double counted).
* Blackbody SEDs only, unless you supply a spectrum.
* No spectral line-spread function: spectral lines are averaged over the resel (equivalent width
  preserved) but not convolved with the instrument profile; the detector ``profile`` serves the
  peak-pixel estimate only.
* :math:`K_{opt}` is sampled every 1 nm. Samples below ``KOPT_EPS = 1e-6`` (687, 761, 1129 nm,
  the 1336-1424 nm water band, 1797 nm) are treated as exactly zero, and so is every wavelength
  whose interpolation interval touches such a sample, so a telluric core is NaN regardless of
  where the bin centres fall (#12). Small but finite samples (the 687 nm B-band wing, 950 nm)
  still give large, finite exposure times.
* The default sky is grey time (moon up); see :doc:`sky` for a dark-time proxy.
* ``grid="uniform"`` samples the band every ``step_nm`` but every bin stays one resel wide, so
  it over- or undersamples; it exists for plotting on a regular grid, not for summing counts.
