Sky
===

The sky background comes from an ESO SkyCalc FITS table wrapped by :class:`~lfast_etc.Sky`.

The bundled table is grey time
------------------------------

``skytable_secz1.3.fits`` is the ESO sky model at 2640 m (SkyCalc has no Kitt Peak site;
Paranal-like aerosols), airmass 1.3, PWV 2.5 mm, **moon up** at 45 degrees altitude and 45
degrees from the target, sampled at R = 20,000. Scattered moonlight is about 60% of the visible
sky: around 550 nm the table gives V ~ 20.3 mag per square arcsecond, roughly one magnitude
brighter than dark Kitt Peak (V ~ 21.7 to 22.0 at zenith, Neugent & Massey 2010).

Consequences:

* Bright targets do not care. The Bench becomes sky limited near V ~ 19.5 in 15 minutes.
* Faint-end limiting magnitudes are 0.3 to 0.5 mag pessimistic for dark time.
* ``Sky(scale=0.4)`` is a crude dark-time proxy. It scales the airglow lines as well as the
  continuum, which is wrong in detail.

.. plot::

    etc = LFASTETC("bench_r1800")
    etc.plot_sky()

The left axis is the resel-averaged sky surface brightness in AB mag per square arcsecond
(brighter is up); the grey curve is the sky electron rate per pixel in the 18 um fiber. The
spikes are airglow lines.

How the sky enters the noise
----------------------------

Every SkyCalc sample inside a resel is averaged (nearest sample if the resel is narrower than the
table sampling, which is the case for the CSHELL presets at R > 20,000, so airglow lines are
pre-smoothed to R = 20,000 there). With the fiber solid angle :math:`\Omega` and the
extended-source throughput :math:`K_{opt}/\eta_{fib}` (the sky fills the fiber and pays no
injection loss),

.. math:: N_{sky} = \frac{N\,A\,\Omega\,(K_{opt}/\eta_{fib})\,S_{opt}\,QE\,B_{sky}\,\delta\lambda}{n_{pix}}
   \quad [\mathrm{e^-\,s^{-1}\,pixel^{-1}}].

The sky-subtraction noise term :math:`(1 + n_{pix}/n_{sky})` is omitted (a sky fiber is assumed
well measured).

Bringing your own table
-----------------------

For real planning generate a table with the conditions you want at
https://www.eso.org/observing/etc/skycalc/ or with ``skycalc_cli``, keep the default FITS output,
and pass it in:

.. code-block:: python

    from lfast_etc import Sky, LFASTETC
    dark = Sky("skytable_dark_secz1.2.fits")
    etc = LFASTETC("bench_r1800", sky=dark)

The ETC reads the columns ``lam`` [nm, vacuum], ``flux`` [ph/s/m\ :sup:`2`/um/arcsec\ :sup:`2`]
and ``trans``. Use a wavelength grid at least as fine as your resel (R >= 50,000 for CSHELL) if
airglow lines matter to you. ``Sky.transmission(wavelength_nm)`` exposes the atmospheric
transmission; it is applied to the source only with ``LFASTETC(apply_extinction=True)``, which
double counts the extinction already inside :math:`K_{opt}` (see :doc:`physics`).
