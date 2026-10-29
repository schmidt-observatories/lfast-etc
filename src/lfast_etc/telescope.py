"""LFAST telescope array: collecting area, fiber, throughput curves.

Throughput model follows Choi et al. 2026, SPIE 14147, 1414736
(https://doi.org/10.1117/12.3103179; copy in input/literature/). Two curves
ship with the package for each fiber option:

* ``K_opt(lambda)`` -- end-to-end telescope throughput for a *point source*.
  This **already includes** the fiber injection efficiency ``eta_fib`` (the
  50 um and 18 um ``K_opt`` curves differ by exactly ``eta_fib_50/eta_fib_18``)
  and the site atmospheric extinction at the reference airmass.
* ``eta_fib(lambda)`` -- point-source fiber injection efficiency.

So a point source sees ``K_opt`` and an extended source (sky) that fills the
fiber uniformly sees ``K_opt / eta_fib``.
"""
from dataclasses import dataclass, field
import warnings
import numpy as np

KOPT_EPS = 1e-6   # K_opt below this is a telluric core: zero throughput, bin invalid (#12)

from .constants import THROUGHPUT_DIR, positive

# Choi et al. 2026 (SPIE 14147, 1414736) throughput budget for one unit telescope, star path.
# These are the wavelength-independent stages the paper quotes; the wavelength-dependent
# remainder (mirror reflectivity, PFC MgF2 coatings and bulk, fiber attenuation) is recovered
# from the shipped K_opt as a residual in ``Telescope.throughput_cascade``.
CHOI_STAGES = (
    ("obscuration: PFC and spider vanes (5.1%)", 0.949),                     # Sec. 2.2
    ("fiber puck: laser-ablation damage allowance", 0.95),                    # Sec. 2.5
    ("fiber feed: Fresnel loss, two circular segments at 92%", 0.92 ** 2),    # Sec. 2.5
    ("fiber combiner: microlens, dichroic, lens, output fiber (78.5%)", 0.785),  # Sec. 2.6
    ("alignment, polarization, vignetting, contamination (4 x 1%)", 0.99 ** 4),  # Sec. 2.7
)
CHOI_SCALAR = float(np.prod([v for _, v in CHOI_STAGES]))
SEEING_PSFS = ("choi", "moffat", "gaussian")
_CHOI_TABLE = {}


def _choi_table():
    """Lazy-load the Choi-calibrated seeing-factor table (scripts/build_seeing_table.py)."""
    if not _CHOI_TABLE:
        from scipy.interpolate import RegularGridInterpolator
        d = np.load(THROUGHPUT_DIR / "seeing_factor_choi.npz")
        for k in ("18", "50"):
            _CHOI_TABLE[k] = RegularGridInterpolator((d["seeing"], d["wavelength"]), d[f"factor_{k}"],
                                                     bounds_error=False, fill_value=None)
        _CHOI_TABLE["seeing"] = (float(d["seeing"].min()), float(d["seeing"].max()))
        _CHOI_TABLE["wavelength"] = (float(d["wavelength"].min()), float(d["wavelength"].max()))
    return _CHOI_TABLE

# fiber core -> (on-sky diameter [arcsec], throughput-file suffix, default combine mode)
FIBERS = {
    "18um": dict(diameter_arcsec=1.48, suffix="18", combine="fiber"),
    "50um": dict(diameter_arcsec=4.1, suffix="50", combine="independent"),
}


@dataclass
class Telescope:
    """An array of identical LFAST unit telescopes feeding one spectrograph.

    Parameters
    ----------
    n_telescopes : int
        Number of unit telescopes. Nominal Kitt Peak LFAST-20X is 20.
    fiber : {"18um", "50um"}
        Fiber core feeding the spectrograph. The 18 um fibers go through the
        N-to-1 fiber combiner (``combine="fiber"``); the 50 um fiber does not,
        so each telescope delivers its own spectrum (``combine="independent"``).
    combine : {"fiber", "independent", None}
        How the unit telescopes are combined. ``"fiber"``: fibers are merged
        in groups of ``telescopes_per_input`` into one spectrograph input each
        (one input for the nominal 20-telescope array), so signal *and* sky per
        pixel scale with the group size while dark current and read noise are
        paid once per input. ``"independent"``: every telescope produces its
        own spectrum which are co-added after the fact; SNR scales as
        ``sqrt(n_telescopes)``. ``None`` picks the default for the fiber.
    telescopes_per_input : int or None
        Size of one fiber combiner (20 for LFAST). Arrays larger than this
        are modelled as ``ceil(n_telescopes / telescopes_per_input)``
        independent inputs co-added incoherently. ``None`` means one input
        regardless of ``n_telescopes`` (an idealised single combiner).
    diameter_m : float
        Unit telescope clear-aperture diameter. The 15.5 cm central obscuration
        is already in the throughput budget; do not double count it.
    fiber_diameter_arcsec : float or None
        Override the on-sky fiber diameter.
    seeing_arcsec : float or None
        Seeing FWHM at 500 nm. ``None`` leaves ``K_opt`` as shipped. Otherwise
        the point-source throughput is scaled by the fraction of a seeing disc
        inside the fiber relative to the same fraction at ``seeing_ref_arcsec``
        (#14); the sky path is unchanged. The seeing is scaled with
        wavelength as ``(lambda / 500 nm) ** -0.2``.
    seeing_ref_arcsec : float
        Seeing assumed by the shipped ``K_opt``/``eta_fib`` curves (Choi et al.
        2026). 1.0 arcsec is the LFAST design seeing (Choi et al. 2024) and is kept until
        the curves' reference is confirmed;
        ``implied_seeing_ref()`` (about 1.5" at 550 nm) is the value at which the
        geometric coupling alone would explain the shipped ``eta_fib``.
    seeing_psf : {"choi", "moffat", "gaussian"}
        How the seeing factor is computed. ``"choi"`` (default) uses a table
        calibrated to the shipped ``eta_fib`` curves: a 1" Moffat (beta 4.765)
        seeing disc convolved with the non-seeing blur (mirror figure, wind
        shake, vibration) that Choi et al. (2026) include, fitted per wavelength
        to both fiber sizes (``scripts/build_seeing_table.py``). ``"moffat"``
        and ``"gaussian"`` are the bare geometric discs, which overstate the
        seeing dependence because they ignore that blur.
    moffat_beta : float
        Moffat index (4.77 for pure Kolmogorov turbulence, 2.5 to 3.5 for
        measured profiles with wings).
    """
    n_telescopes: int = 20
    fiber: str = "18um"
    combine: str | None = None
    diameter_m: float = 0.726
    fiber_diameter_arcsec: float | None = None
    telescopes_per_input: int | None = 20
    seeing_arcsec: float | None = None
    seeing_ref_arcsec: float = 1.0
    seeing_psf: str = "choi"
    moffat_beta: float = 3.5
    _kopt: tuple = field(init=False, default=None, repr=False, compare=False)
    _eta: tuple = field(init=False, default=None, repr=False, compare=False)

    def __post_init__(self):
        if self.fiber not in FIBERS:
            raise ValueError(f"fiber must be one of {list(FIBERS)}, got {self.fiber!r}")
        if int(self.n_telescopes) != self.n_telescopes or self.n_telescopes < 1:
            raise ValueError(f"n_telescopes must be a positive integer, got {self.n_telescopes!r}")
        self.n_telescopes = int(self.n_telescopes)
        self.diameter_m = positive("diameter_m", self.diameter_m)
        if self.telescopes_per_input is not None:
            if int(self.telescopes_per_input) != self.telescopes_per_input or self.telescopes_per_input < 1:
                raise ValueError("telescopes_per_input must be a positive integer or None")
            self.telescopes_per_input = int(self.telescopes_per_input)
        if self.fiber_diameter_arcsec is not None:
            self.fiber_diameter_arcsec = positive("fiber_diameter_arcsec", self.fiber_diameter_arcsec)
        spec = FIBERS[self.fiber]
        if self.combine is None:
            self.combine = spec["combine"]
        if self.combine not in ("fiber", "independent"):
            raise ValueError("combine must be 'fiber' or 'independent'")
        if self.combine == "fiber" and spec["combine"] != "fiber":
            warnings.warn(f"the {self.fiber} fiber does not go through the fiber combiner; "
                          "combine='fiber' is a what-if, not a planned LFAST configuration", stacklevel=2)
        if self.fiber_diameter_arcsec is None:
            self.fiber_diameter_arcsec = spec["diameter_arcsec"]
        if self.seeing_arcsec is not None:
            self.seeing_arcsec = positive("seeing_arcsec", self.seeing_arcsec)
        self.seeing_ref_arcsec = positive("seeing_ref_arcsec", self.seeing_ref_arcsec)
        self.moffat_beta = positive("moffat_beta", self.moffat_beta)
        if self.moffat_beta <= 1:
            raise ValueError("moffat_beta must exceed 1 for a normalisable profile")
        if self.seeing_psf not in SEEING_PSFS:
            raise ValueError(f"seeing_psf must be one of {SEEING_PSFS}, got {self.seeing_psf!r}")
        s = spec["suffix"]
        self._kopt = (np.load(THROUGHPUT_DIR / "Kopt_wavelength.npy"),
                      np.load(THROUGHPUT_DIR / f"Kopt_value_{s}.npy"))
        self._eta = (np.load(THROUGHPUT_DIR / "Eta_fib_wavelength.npy"),
                     np.load(THROUGHPUT_DIR / f"Eta_fib_value_{s}.npy"))

    # ---- geometry -----------------------------------------------------------
    @property
    def area_m2(self):
        """Collecting area of one unit telescope [m^2]."""
        return np.pi * (self.diameter_m / 2.0) ** 2

    @property
    def n_inputs(self):
        """Number of independent spectrograph inputs (detector read-outs)."""
        if self.combine == "independent":
            return self.n_telescopes
        if self.telescopes_per_input is None:
            return 1
        return -(-self.n_telescopes // self.telescopes_per_input)

    @property
    def telescopes_on_brightest_input(self):
        """Telescopes feeding the most heavily loaded input: inputs fill in packed
        groups of ``telescopes_per_input``, so N=21 is one full group plus one."""
        if self.combine == "independent":
            return 1
        if self.telescopes_per_input is None:
            return self.n_telescopes
        return min(self.n_telescopes, self.telescopes_per_input)

    @property
    def equivalent_diameter_m(self):
        """Diameter of a single telescope with the same collecting area [m]."""
        return self.diameter_m * np.sqrt(self.n_telescopes)

    @property
    def total_area_m2(self):
        """Collecting area of the whole array [m^2]."""
        return self.n_telescopes * self.area_m2

    @property
    def fiber_solid_angle_arcsec2(self):
        """On-sky solid angle of the fiber [arcsec^2]."""
        return np.pi * (self.fiber_diameter_arcsec / 2.0) ** 2

    @property
    def wavelength_range_nm(self):
        w = self._kopt[0]
        return float(w.min()), float(w.max())

    # ---- throughput ---------------------------------------------------------
    def kopt(self, wavelength_nm):
        """Point-source telescope throughput (includes fiber injection).

        Samples below ``KOPT_EPS`` are telluric cores and count as exactly
        zero, and so does any wavelength whose interpolation interval touches
        such a sample (#12): the rule does not depend on where the 1 nm
        sampling happens to land relative to the bin centres.
        """
        w, k = self._kopt
        clear = np.interp(wavelength_nm, w, (k >= KOPT_EPS).astype(float)) >= 1.0
        return np.where(clear, np.interp(wavelength_nm, w, k), 0.0)

    def eta_fib(self, wavelength_nm):
        """Point-source fiber injection efficiency."""
        return np.interp(wavelength_nm, *self._eta)

    # ---- seeing (#14) -------------------------------------------------------
    def fiber_coupling(self, seeing_arcsec, wavelength_nm=500.0, psf=None):
        """Fraction of a point source's light inside the fiber for a seeing disc of
        FWHM ``seeing_arcsec`` at 500 nm, scaled to ``wavelength_nm`` as
        ``(lambda / 500 nm) ** -0.2``. Geometry only: no telescope diffraction, no
        focal-ratio degradation, so this is the *relative* seeing dependence.
        """
        seeing = np.asarray(seeing_arcsec, float)
        if not np.all(np.isfinite(seeing)) or np.any(seeing <= 0):
            raise ValueError(f"seeing_arcsec must be > 0, got {seeing_arcsec!r}")
        fwhm = seeing * (np.asarray(wavelength_nm, float) / 500.0) ** -0.2
        r = self.fiber_diameter_arcsec / 2.0
        psf = psf or self.seeing_psf
        if psf == "choi":
            psf = "moffat"            # the geometric disc; the calibrated model lives in seeing_factor
        if psf == "gaussian":
            sigma = fwhm / 2.3548
            return 1.0 - np.exp(-0.5 * (r / sigma) ** 2)
        if psf != "moffat":
            raise ValueError("psf must be 'moffat' or 'gaussian'")
        b = self.moffat_beta
        alpha = fwhm / (2.0 * np.sqrt(2.0 ** (1.0 / b) - 1.0))
        return 1.0 - (1.0 + (r / alpha) ** 2) ** (1.0 - b)

    def implied_seeing_ref(self, wavelength_nm=500.0, psf=None):
        """Seeing [arcsec at 500 nm] at which the geometric ``fiber_coupling`` equals the
        shipped ``eta_fib`` at ``wavelength_nm``. If ``eta_fib`` were pure seeing loss this
        would be the reference seeing of the throughput curves; it is an upper limit,
        since focal-ratio degradation, alignment and diffraction (1.22 lambda/D = 0.57"
        at 1650 nm for 0.726 m, growing with wavelength while seeing shrinks) also
        reduce ``eta_fib``. About 1.5" at 550 nm for the 18 um fiber with the Moffat disc.
        """
        from scipy.optimize import brentq
        target = float(self.eta_fib(wavelength_nm))
        f = lambda s: float(self.fiber_coupling(s, wavelength_nm, psf=psf)) - target
        return float(brentq(f, 0.05, 20.0))

    def seeing_factor(self, wavelength_nm):
        """Point-source throughput relative to the reference seeing; 1 when ``seeing_arcsec`` is None.

        ``"choi"``: interpolated from the calibrated table (seeing clamped to its
        0.4-3.0" range, wavelength to 400-1800 nm). Otherwise
        ``fiber_coupling(seeing) / fiber_coupling(seeing_ref)``.
        """
        if self.seeing_arcsec is None:
            return 1.0
        if self.seeing_psf != "choi":
            return self.fiber_coupling(self.seeing_arcsec, wavelength_nm) / self.fiber_coupling(self.seeing_ref_arcsec, wavelength_nm)
        t = _choi_table(); key = FIBERS[self.fiber]["suffix"]
        w = np.clip(np.atleast_1d(np.asarray(wavelength_nm, float)), *t["wavelength"])
        s, s0 = (float(np.clip(x, *t["seeing"])) for x in (self.seeing_arcsec, self.seeing_ref_arcsec))
        pts = lambda x: np.stack([np.full_like(w, x), w], axis=-1)
        out = t[key](pts(s)) / t[key](pts(s0))
        return float(out[0]) if np.ndim(wavelength_nm) == 0 else out

    def throughput_cascade(self, wavelength_nm, path="star", sky=None):
        """The Choi et al. (2026) throughput budget as an ordered cascade of stages.

        Returns ``{"wavelength", "path", "stages": [{"name", "transmission", "source"}, ...],
        "cumulative": [array, ...]}``; the last cumulative entry equals
        ``throughput_point_source`` (``path="star"``) or ``throughput_extended``
        (``path="sky"``, no fiber injection, no seeing). Stages the paper quotes as
        wavelength-independent are ``CHOI_STAGES``; the wavelength-dependent remainder
        (mirror reflectivity, PFC coatings and bulk, fiber attenuation) is the residual
        ``K_opt / eta_fib / T_atm / CHOI_SCALAR``. ``K_opt`` carries the paper's own
        LBLRTM atmosphere (airmass 1.3, 5 mm PWV); the SkyCalc ``sky`` transmission
        (airmass 1.3, 2.5 mm) is only a proxy for showing it as a stage, so the
        residual absorbs the difference between the two atmospheres.
        """
        from .sky import Sky
        scalar = np.ndim(wavelength_nm) == 0
        w = np.atleast_1d(np.asarray(wavelength_nm, float))
        kopt, eta = self.kopt(w), self.eta_fib(w)
        H = np.where(eta > 0, kopt / np.where(eta > 0, eta, 1.0), 0.0)             # hardware incl. atmosphere
        T = np.maximum((sky or Sky()).transmission(w), 1e-3)
        residual = np.minimum(1.0, H / (T * CHOI_SCALAR))
        T_stage = np.where(residual > 0, H / np.where(residual > 0, residual, 1.0) / CHOI_SCALAR, T)
        stages = [dict(name="atmosphere (in K_opt: LBLRTM, airmass 1.3, 5 mm PWV; SkyCalc shown)", transmission=T_stage, source="Choi 2026 Sec. 2.1 / SkyCalc proxy"),
                  dict(name=CHOI_STAGES[0][0], transmission=np.full_like(w, CHOI_STAGES[0][1]), source="Choi 2026 Sec. 2.2"),
                  dict(name="mirror reflectivity + PFC coatings and bulk + fiber attenuation (derived residual)", transmission=residual, source="K_opt / eta_fib / T_atm / scalar stages")]
        if path == "star":
            stages.append(dict(name="fiber injection eta_fib at 1\" seeing (mirror figure, seeing, wind, vibration)", transmission=eta, source="Choi 2026 Sec. 2.4 (shipped curve)"))
            if self.seeing_arcsec is not None:
                stages.append(dict(name=f"seeing factor: {self.seeing_arcsec:g}\" vs {self.seeing_ref_arcsec:g}\" ({self.seeing_psf})", transmission=self.seeing_factor(w) * np.ones_like(w), source="lfast_etc seeing model"))
        elif path != "sky":
            raise ValueError("path must be 'star' or 'sky'")
        for name, v in CHOI_STAGES[1:]:
            stages.append(dict(name=name, transmission=np.full_like(w, v), source="Choi 2026"))
        cum, run = [], np.ones_like(w)
        for s in stages:
            run = run * s["transmission"]; cum.append(run)
        if scalar:
            for s in stages:
                s["transmission"] = float(s["transmission"][0])
            cum = [float(c[0]) for c in cum]
        return dict(wavelength=float(w[0]) if scalar else w, path=path, stages=stages, cumulative=cum)

    def plot_fiber_coupling(self, ax=None, seeing=None):
        """Fraction of light inside the fiber versus seeing: the geometric Moffat and Gaussian
        discs and the Choi-calibrated curve (``eta_fib`` at 500 nm times the table factor),
        with the current and reference seeing marked."""
        from .plotting import plot_fiber_coupling
        return plot_fiber_coupling(self, ax=ax, seeing=seeing)

    def throughput_point_source(self, wavelength_nm):
        """Telescope throughput for a point source: ``K_opt`` times the seeing factor."""
        return self.kopt(wavelength_nm) * self.seeing_factor(wavelength_nm)

    def throughput_extended(self, wavelength_nm):
        """Telescope throughput for a source filling the fiber: ``K_opt / eta_fib``."""
        return self.kopt(wavelength_nm) / self.eta_fib(wavelength_nm)

    def describe(self):
        return (f"{self.n_telescopes} x {self.diameter_m:.3f} m telescopes "
                f"(total area {self.total_area_m2:.2f} m^2, {self.equivalent_diameter_m:.1f} m equivalent), "
                f"{self.fiber} fiber ({self.fiber_diameter_arcsec:.2f}\" on sky), "
                f"combine={self.combine}, {self.n_inputs} input(s)"
                + (f", seeing {self.seeing_arcsec:g}\" ({self.seeing_psf}, ref {self.seeing_ref_arcsec:g}\")"
                   if self.seeing_arcsec is not None else ""))
