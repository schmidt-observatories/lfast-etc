Sources
=======

A :class:`~lfast_etc.Source` is a relative spectral shape :math:`s(\lambda)` and one magnitude
in one band; the magnitude fixes the absolute scale. ``source.flux(wavelength_nm)`` returns
:math:`F_\lambda` at Earth in W m\ :sup:`-2` nm\ :sup:`-1`, ``photon_flux`` the same in photons
s\ :sup:`-1` m\ :sup:`-2` nm\ :sup:`-1`.

Constructors
------------

.. code-block:: python

    from lfast_etc import Source
    Source.blackbody(5800, mag=15, band="V")                       # Vega magnitudes by default
    Source.flat(mag=20, band="R", system="AB")                     # constant F_nu (kind="flam" for constant F_lambda)
    Source.powerlaw(alpha=-0.7, mag=18, band="R", system="AB")     # F_nu ~ nu^alpha
    Source.from_arrays(wavelength_nm, flux, flux_unit="flam", mag=16, band="I")
    Source.from_arrays(wavelength_nm, flux_W_m2_nm, flux_unit="W/m2/nm", band="V")   # absolute: mag is derived

Note the defaults: ``blackbody`` is Vega, ``flat`` and ``powerlaw`` are AB. Pass ``system=``
explicitly when it matters.

Bands
-----

Magnitudes are defined through box bands (centre, width, Vega zero point; Bessell et al. 1998,
V as in the reference notebook). "AB" uses 3631 Jy in every band.

=====  ===========  ==========  ===============
band   centre [nm]  width [nm]  Vega zero [Jy]
=====  ===========  ==========  ===============
U      366          66          1790
B      438          94          4063
V      550          90          3640
R      641          138         3064
I      798          149         2416
J      1220         213         1589
H      1630         307         1021
K      2190         390         640
=====  ===========  ==========  ===============

The band flux of magnitude :math:`m` is :math:`f_{\nu,0}\,10^{-0.4 m}\,c\,\Delta\lambda/\lambda_0^2`;
the exact :math:`c\Delta\lambda/(\lambda_{lo}\lambda_{hi})` differs by about 1% (0.01 mag).
``source.magnitude("J")`` integrates the same SED over another box band, which is how you get a
J magnitude for a star you only know in V (notebook 05 does this for He 10830 hosts known in J
but observed in H).

Spectral lines
--------------

.. code-block:: python

    src = Source.blackbody(5800, 15).add_line(656.3, depth=0.3, sigma_nm=0.05)    # 30% absorption
    src.add_line(486.1, depth=-2.0, sigma_nm=0.05)                                 # emission, 3x the continuum at the core

Lines are multiplicative Gaussians applied *after* the normalisation, so ``mag`` stays the
continuum magnitude. ``depth`` is the fraction of the continuum removed at the line centre
(at most 1); negative values are emission. On the ETC grid each line is averaged over the bin
:math:`[\lambda - \delta\lambda/2, \lambda + \delta\lambda/2]`, which preserves the equivalent
width: a 0.05 nm line on a 0.36 nm resel is a 10% dip, not the 30% you would see at infinite
resolution. There is no line-spread function on top of that averaging.

.. plot::

    src = Source.blackbody(5800, mag=15, band="V").add_line(656.3, 0.3, 0.05).add_line(486.1, 0.2, 0.05)
    etc = LFASTETC("bench_r1800")
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.4))
    src.plot(wavelength_nm=np.linspace(400, 700, 3000), ax=axes[0])
    w = etc.wavelength
    axes[1].plot(w, src.line_factor(w), lw=0.8, label="point sampled")
    axes[1].plot(w, src.line_factor(w, etc.delta_lambda), lw=1.2, label="averaged over the resel (what the ETC uses)")
    axes[1].set_xlim(650, 662); axes[1].set_xlabel("Wavelength (nm)"); axes[1].set_ylabel("line factor"); axes[1].legend(fontsize=8)
    fig.tight_layout()

User spectra
------------

``from_arrays`` accepts any sampled spectrum. ``flux_unit="flam"`` means the array is
proportional to :math:`F_\lambda` (any unit), ``"fnu"`` proportional to :math:`F_\nu`, and
``"W/m2/nm"`` an absolute :math:`F_\lambda` from which the magnitude is derived. Outside the
array's wavelength range the spectrum is zero; the normalising band must lie inside it or a
``ValueError`` is raised. Notebook 09 builds a flare spectrum this way (a 3300 K star plus a
10,000 K component, then an H-alpha emission line).

Blackbodies are a rough stand-in for real stars, especially M dwarfs at R ~ 1800 where TiO bands
dominate the visible SED; see the open issue on spectral templates.

Other magnitudes, copies
------------------------

.. code-block:: python

    src.magnitude("H")          # synthetic magnitude in another box band (same system)
    src.magnitude("H", "AB")
    src.with_mag(18.0)          # same shape, band, system and lines at another magnitude
    src.describe()              # 'blackbody 5800 K, V=15.00 (Vega)'

The sensitivity functions (:doc:`sensitivity`) call ``with_mag(0)`` internally and ignore the
magnitude you built the source with; only the shape matters there.
