"""LFAST telescope array: collecting area, fiber, throughput curves.

Throughput model follows Choi et al. 2026, SPIE 14147
(https://doi.org/10.1117/12.3103179). Two curves ship with the package for
each fiber option:

* ``K_opt(lambda)`` -- end-to-end telescope throughput for a *point source*.
  This **already includes** the fiber injection efficiency ``eta_fib`` (the
  50 um and 18 um ``K_opt`` curves differ by exactly ``eta_fib_50/eta_fib_18``)
  and the site atmospheric extinction at the reference airmass.
* ``eta_fib(lambda)`` -- point-source fiber injection efficiency.

So a point source sees ``K_opt`` and an extended source (sky) that fills the
fiber uniformly sees ``K_opt / eta_fib``.
"""
from dataclasses import dataclass, field
import numpy as np

from .constants import THROUGHPUT_DIR, positive

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
        How the unit telescopes are combined. ``"fiber"``: all fibers are
        merged into one spectrograph input, so signal *and* sky per pixel scale
        with ``n_telescopes`` while dark current and read noise do not.
        ``"independent"``: every telescope produces its own spectrum which are
        co-added after the fact; SNR scales as ``sqrt(n_telescopes)``.
        ``None`` picks the default for the fiber.
    diameter_m : float
        Unit telescope clear-aperture diameter. The 15.5 cm central obscuration
        is already in the throughput budget; do not double count it.
    fiber_diameter_arcsec : float or None
        Override the on-sky fiber diameter.
    """
    n_telescopes: int = 20
    fiber: str = "18um"
    combine: str | None = None
    diameter_m: float = 0.726
    fiber_diameter_arcsec: float | None = None
    _kopt: tuple = field(init=False, default=None, repr=False, compare=False)
    _eta: tuple = field(init=False, default=None, repr=False, compare=False)

    def __post_init__(self):
        if self.fiber not in FIBERS:
            raise ValueError(f"fiber must be one of {list(FIBERS)}, got {self.fiber!r}")
        if int(self.n_telescopes) != self.n_telescopes or self.n_telescopes < 1:
            raise ValueError(f"n_telescopes must be a positive integer, got {self.n_telescopes!r}")
        self.n_telescopes = int(self.n_telescopes)
        self.diameter_m = positive("diameter_m", self.diameter_m)
        if self.fiber_diameter_arcsec is not None:
            self.fiber_diameter_arcsec = positive("fiber_diameter_arcsec", self.fiber_diameter_arcsec)
        spec = FIBERS[self.fiber]
        if self.combine is None:
            self.combine = spec["combine"]
        if self.combine not in ("fiber", "independent"):
            raise ValueError("combine must be 'fiber' or 'independent'")
        if self.fiber_diameter_arcsec is None:
            self.fiber_diameter_arcsec = spec["diameter_arcsec"]
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
        """Point-source telescope throughput (includes fiber injection)."""
        return np.interp(wavelength_nm, *self._kopt)

    def eta_fib(self, wavelength_nm):
        """Point-source fiber injection efficiency."""
        return np.interp(wavelength_nm, *self._eta)

    def throughput_point_source(self, wavelength_nm):
        """Telescope throughput for a point source: ``K_opt``."""
        return self.kopt(wavelength_nm)

    def throughput_extended(self, wavelength_nm):
        """Telescope throughput for a source filling the fiber: ``K_opt / eta_fib``."""
        return self.kopt(wavelength_nm) / self.eta_fib(wavelength_nm)

    def describe(self):
        return (f"{self.n_telescopes} x {self.diameter_m:.3f} m telescopes "
                f"(total area {self.total_area_m2:.2f} m^2), {self.fiber} fiber "
                f"({self.fiber_diameter_arcsec:.2f}\" on sky), combine={self.combine}")
