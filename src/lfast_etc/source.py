"""Target spectral energy distributions, normalised to a broadband magnitude.

A :class:`Source` returns ``F_lambda`` at Earth in W m^-2 nm^-1. The shape
comes from a blackbody, a power law, a flat spectrum or a user array; the
absolute scale is set by one magnitude in one band.
"""
from dataclasses import dataclass, field
import numpy as np

from .constants import H_PLANCK, C_LIGHT, K_BOLTZMANN

# Box bandpasses: centre [nm], width [nm], Vega zero point [Jy].
# V follows the reference LFAST notebook (Bessell et al. 1998 rounded);
# the rest are Bessell et al. 1998 / Bessell 2005. "AB" uses 3631 Jy in every band.
# The band width in frequency is approximated as c*dlam/lam0^2 (as in the
# reference notebook); the exact c*dlam/(lam_lo*lam_hi) differs by ~1%.
BANDS = {
    "U": (366.0, 66.0, 1790.0),
    "B": (438.0, 94.0, 4063.0),
    "V": (550.0, 90.0, 3640.0),
    "R": (641.0, 138.0, 3064.0),
    "I": (798.0, 149.0, 2416.0),
    "J": (1220.0, 213.0, 1589.0),
    "H": (1630.0, 307.0, 1021.0),
    "K": (2190.0, 390.0, 640.0),
}
AB_ZP_JY = 3631.0


def band_flux_W_m2(mag, band, system="Vega"):
    """Band-integrated flux [W/m^2] for an apparent magnitude in a box band."""
    center, width, zp_vega = BANDS[band]
    zp = AB_ZP_JY if system.upper() == "AB" else zp_vega
    f_nu = zp * 1e-26 * 10 ** (-0.4 * mag)                  # W/m^2/Hz
    delta_nu = C_LIGHT * (width * 1e-9) / (center * 1e-9) ** 2
    return f_nu * delta_nu


def planck_W_m2_nm_sr(wavelength_nm, teff):
    """Planck specific intensity B_lambda [W m^-2 nm^-1 sr^-1]."""
    wl = np.asarray(wavelength_nm, dtype=float) * 1e-9
    x = H_PLANCK * C_LIGHT / (wl * K_BOLTZMANN * teff)
    with np.errstate(over="ignore"):
        b = 2 * H_PLANCK * C_LIGHT ** 2 / (wl ** 5 * np.expm1(x))
    return b * 1e-9


@dataclass
class Source:
    """A point source with a known SED shape and a magnitude.

    Build with the class methods :meth:`blackbody`, :meth:`flat`,
    :meth:`powerlaw` or :meth:`from_arrays`, then optionally add
    Gaussian spectral lines with :meth:`add_line`.
    """
    shape: object                     # callable wavelength_nm -> relative F_lambda
    mag: float                        # continuum magnitude (lines added later are not renormalised)
    band: str = "V"
    system: str = "Vega"
    label: str = "source"
    lines: list = field(default_factory=list)
    _scale: float = field(default=None, repr=False)

    def __post_init__(self):
        if self.band not in BANDS:
            raise ValueError(f"band must be one of {list(BANDS)}")
        self.mag = float(self.mag)
        c, w, _ = BANDS[self.band]
        grid = np.linspace(c - w / 2, c + w / 2, 2001)
        integral = np.trapezoid(self.shape(grid), grid)
        if not np.isfinite(integral) or integral <= 0:
            raise ValueError("source produces no flux in the normalising band")
        self._scale = band_flux_W_m2(self.mag, self.band, self.system) / integral

    # ---- constructors -------------------------------------------------------
    @classmethod
    def blackbody(cls, teff, mag, band="V", system="Vega"):
        """Blackbody of temperature ``teff`` [K]."""
        return cls(lambda wl: planck_W_m2_nm_sr(wl, teff), mag, band, system,
                   label=f"blackbody {teff:.0f} K")

    @classmethod
    def flat(cls, mag, band="V", system="AB", kind="fnu"):
        """Flat spectrum: constant ``F_nu`` (``kind="fnu"``) or ``F_lambda``."""
        if kind == "fnu":
            shape = lambda wl: 1.0 / np.asarray(wl, float) ** 2
        else:
            shape = lambda wl: np.ones_like(np.asarray(wl, float))
        return cls(shape, mag, band, system, label=f"flat {kind}")

    @classmethod
    def powerlaw(cls, alpha, mag, band="V", system="AB"):
        """Power law ``F_nu ~ nu^alpha`` (so ``F_lambda ~ lambda^(-alpha-2)``)."""
        return cls(lambda wl: np.asarray(wl, float) ** (-alpha - 2), mag, band, system,
                   label=f"powerlaw alpha={alpha}")

    @classmethod
    def from_arrays(cls, wavelength_nm, flux, flux_unit="flam", mag=None, band="V", system="Vega"):
        """User spectrum. ``flux_unit``: ``"flam"`` (any F_lambda unit),
        ``"fnu"`` (any F_nu unit) or ``"W/m2/nm"`` (absolute; ``mag`` ignored
        and the band-integrated flux defines the magnitude)."""
        wl = np.asarray(wavelength_nm, float)
        fl = np.asarray(flux, float)
        if wl.ndim != 1 or wl.shape != fl.shape:
            raise ValueError("wavelength_nm and flux must be 1-D arrays of equal length")
        order = np.argsort(wl)
        wl, fl = wl[order], fl[order]
        if flux_unit == "fnu":
            fl = fl / wl ** 2
        shape = lambda w: np.interp(w, wl, fl, left=0.0, right=0.0)
        if flux_unit == "W/m2/nm":
            c, w, _ = BANDS[band]
            grid = np.linspace(c - w / 2, c + w / 2, 2001)
            f_band = np.trapezoid(shape(grid), grid)
            zp = AB_ZP_JY if system.upper() == "AB" else BANDS[band][2]
            delta_nu = C_LIGHT * (w * 1e-9) / (c * 1e-9) ** 2
            mag = -2.5 * np.log10(f_band / delta_nu / (zp * 1e-26))
        elif mag is None:
            raise ValueError("mag is required unless flux_unit='W/m2/nm'")
        return cls(shape, mag, band, system, label="user spectrum")

    # ---- lines --------------------------------------------------------------
    def add_line(self, center_nm, depth, sigma_nm):
        """Add a Gaussian line. ``depth`` > 0 is absorption (fraction of
        continuum removed at the core, at most 1), < 0 is emission. The
        magnitude stays that of the continuum."""
        if depth > 1:
            raise ValueError("absorption depth cannot exceed 1")
        self.lines.append((float(center_nm), float(depth), float(sigma_nm)))
        return self

    # ---- evaluation ---------------------------------------------------------
    def flux(self, wavelength_nm):
        """``F_lambda`` at Earth [W m^-2 nm^-1], lines included."""
        wl = np.asarray(wavelength_nm, float)
        f = self._scale * self.shape(wl)
        for c, d, s in self.lines:
            f = f * (1 - d * np.exp(-0.5 * ((wl - c) / s) ** 2))
        return f

    def photon_flux(self, wavelength_nm):
        """Photon flux density [photons s^-1 m^-2 nm^-1]."""
        wl = np.asarray(wavelength_nm, float)
        return self.flux(wl) / (H_PLANCK * C_LIGHT / (wl * 1e-9))

    def with_mag(self, mag):
        """Copy of this source at another magnitude (same band/system/lines)."""
        return Source(self.shape, mag, self.band, self.system, self.label, list(self.lines))

    def magnitude(self, band, system=None):
        """Synthetic magnitude of this source in another box band."""
        system = system or self.system
        c, w, zp_vega = BANDS[band]
        grid = np.linspace(c - w / 2, c + w / 2, 2001)
        f_band = np.trapezoid(self.flux(grid), grid)
        zp = AB_ZP_JY if system.upper() == "AB" else zp_vega
        delta_nu = C_LIGHT * (w * 1e-9) / (c * 1e-9) ** 2
        return -2.5 * np.log10(f_band / delta_nu / (zp * 1e-26))

    def describe(self):
        return f"{self.label}, {self.band}={self.mag:.2f} ({self.system})"
