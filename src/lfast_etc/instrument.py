"""Spectrograph + detector descriptions, loaded from TOML presets."""
import dataclasses
from dataclasses import dataclass
from functools import lru_cache
import tomllib

import numpy as np

from .constants import CONFIG_DIR, positive

PROFILES = ("uniform", "tophat", "gaussian")
BLUR_SIGMA_PX = 0.35      # camera blur applied to the top-hat fiber image [pixels], ~0.8 px FWHM
_SUB = 20                 # sub-pixels per pixel when rendering the profile


@lru_cache(maxsize=256)
def pixel_image(width_px, height_px, profile="tophat"):
    """The fiber image on the detector as fractions of the resel's light per pixel.

    ``width_px`` is the extent along the dispersion direction (one resel) and
    ``height_px`` the spatial extent; the image is centred on a pixel centre,
    which maximises the peak. ``"uniform"`` spreads the light evenly over
    ``width_px * height_px`` pixels (the pre-#19 convention); ``"tophat"`` is
    an elliptical disc, the image of the fiber face, convolved with a Gaussian
    camera blur of ``BLUR_SIGMA_PX``; ``"gaussian"`` has those FWHMs.
    """
    if profile not in PROFILES:
        raise ValueError(f"profile must be one of {PROFILES}, got {profile!r}")
    w, h = positive("width_px", width_px), positive("height_px", height_px)
    if profile == "uniform":
        out = np.full((1, 1), 1.0 / (w * h)); out.setflags(write=False)
        return out
    margin = 3 * max(w, h) / 2.3548 if profile == "gaussian" else 3 * max(BLUR_SIGMA_PX, 0.5)   # Gaussian wings need room
    half = int(np.ceil(max(w, h) / 2 + margin)) + 1
    x = (np.arange(-half * _SUB, half * _SUB) + 0.5) / _SUB - 0.5    # sub-pixel centres; pixel k spans [k-0.5, k+0.5)
    xx, yy = np.meshgrid(x, x)
    if profile == "tophat":
        from scipy.ndimage import gaussian_filter
        img = ((2 * xx / w) ** 2 + (2 * yy / h) ** 2 <= 1).astype(float)
        img = gaussian_filter(img, BLUR_SIGMA_PX * _SUB, mode="constant")
    else:
        sx, sy = w / 2.3548, h / 2.3548
        img = np.exp(-0.5 * ((xx / sx) ** 2 + (yy / sy) ** 2))
    img /= img.sum()
    n = 2 * half
    out = img.reshape(n, _SUB, n, _SUB).sum(axis=(1, 3))
    out.setflags(write=False)                       # cached: callers must not mutate it
    return out


def peak_pixel_fraction(width_px, height_px, profile="tophat"):
    """Largest fraction of a resel's light landing in one pixel (see :func:`pixel_image`)."""
    return float(pixel_image(float(width_px), float(height_px), profile).max())


def list_instruments():
    """Names of the bundled instrument presets."""
    return sorted(p.stem for p in CONFIG_DIR.glob("*.toml"))


@dataclass
class Instrument:
    """A spectrograph configuration.

    All parameters are held constant with wavelength; ``qe`` may be a scalar
    or a callable ``qe(wavelength_nm)`` for a wavelength-dependent curve.

    Parameters
    ----------
    name : str
    wavelength_min_nm, wavelength_max_nm : float
        Bandpass covered by this configuration.
    resolving_power : float
        ``R = lambda / delta_lambda``.
    pixels_per_resel : float
        Pixels spanned by one resolution element along the dispersion
        direction (the width of the fiber image, or of the slit, in pixels).
    fiber_image_px : float or None
        Spatial height of the fiber image on the detector in pixels at 1x1
        binning. ``None`` (default) assumes a circular fiber image, i.e.
        ``pixels_per_resel``; slit modes that narrow the dispersion direction
        only should set the native fiber-image height explicitly. Howell's
        extraction aperture is ``n_pix = pixels_per_resel * fiber_image_px``
        (#13); ``1.0`` reproduces the reference notebook's dispersion-only count.
    profile : {"tophat", "uniform", "gaussian"}
        Shape of the fiber image on the detector, used for the peak-pixel
        (saturation) estimate only (#19). ``"tophat"`` (default) is the image
        of the fiber face with a small camera blur; ``"uniform"`` is the
        pre-#19 convention (peak = ``1 / n_pix``); ``"gaussian"`` has the same
        FWHMs. Seeing does not enter: the fiber scrambles the input.
    spectrograph_efficiency : float
        Spectrograph optical throughput (``S_opt``), excluding slit losses and QE.
    slit_throughput : float
        Extra multiplicative throughput for slit modes (1 for the bare fiber).
        Applied to the sky as well as the source: the slit clips the uniformly
        illuminated fiber face, so both lose the same fraction.
    qe : float or callable
        Detector quantum efficiency.
    dark_current : float
        e-/s/pixel.
    read_noise : float
        e-/pixel per read.
    full_well : float
        e-/pixel saturation level (used for the saturation flag only). ``inf``
        or ``None`` means unknown: nothing is flagged and
        ``ETCResult.summary()["saturated_fraction"]`` is ``None``.
    min_exptime : float
        Shortest exposure [s] for which ``read_noise`` applies (e.g. the
        Fowler-sampling minimum). Results flag reads shorter than this in
        ``ETCResult.below_min_exptime``.
    free_spectral_range_frac : float or None
        ``FSR / lambda`` for echelle instruments; a single setting covers only
        ``lambda * free_spectral_range_frac``. Informational.
    """
    name: str
    wavelength_min_nm: float
    wavelength_max_nm: float
    resolving_power: float
    pixels_per_resel: float
    spectrograph_efficiency: float
    qe: object = 0.8
    slit_throughput: float = 1.0
    dark_current: float = 0.0
    read_noise: float = 0.0
    full_well: float | None = float("inf")
    min_exptime: float = 0.0
    free_spectral_range_frac: float | None = None
    description: str = ""
    fiber_image_px: float | None = None
    profile: str = "tophat"

    def __post_init__(self):
        for f in ("wavelength_min_nm", "wavelength_max_nm", "resolving_power", "pixels_per_resel",
                  "spectrograph_efficiency", "slit_throughput"):
            setattr(self, f, positive(f, getattr(self, f)))
        if self.full_well is None or self.full_well == float("inf"):
            self.full_well = float("inf")                 # unknown: saturation not assessed
        else:
            self.full_well = positive("full_well", self.full_well)
        for f in ("dark_current", "read_noise", "min_exptime"):
            setattr(self, f, positive(f, getattr(self, f), allow_zero=True))
        if not callable(self.qe):
            self.qe = positive("qe", self.qe)
        for f in ("spectrograph_efficiency", "slit_throughput") + (() if callable(self.qe) else ("qe",)):
            if getattr(self, f) > 1:
                raise ValueError(f"{f} is an efficiency and must be <= 1, got {getattr(self, f)}")
        if self.wavelength_max_nm <= self.wavelength_min_nm:
            raise ValueError("wavelength_max_nm must exceed wavelength_min_nm")
        if self.fiber_image_px is not None:
            self.fiber_image_px = positive("fiber_image_px", self.fiber_image_px)
        if self.profile not in PROFILES:
            raise ValueError(f"profile must be one of {PROFILES}, got {self.profile!r}")

    # ---- extraction aperture and detector profile (#13, #19) ----------------
    @property
    def spatial_px(self):
        """Fiber-image height in pixels (``fiber_image_px``, or circular: ``pixels_per_resel``)."""
        return self.pixels_per_resel if self.fiber_image_px is None else self.fiber_image_px

    @property
    def n_pix(self):
        """Howell's extraction aperture: ``pixels_per_resel * spatial_px`` pixels per resel."""
        return self.pixels_per_resel * self.spatial_px

    @property
    def peak_fraction(self):
        """Fraction of a resel's light in its brightest pixel for this ``profile``."""
        return peak_pixel_fraction(self.pixels_per_resel, self.spatial_px, self.profile)

    def pixel_image(self):
        """Fiber image as per-pixel fractions of the resel's light (see :func:`pixel_image`)."""
        return pixel_image(float(self.pixels_per_resel), float(self.spatial_px), self.profile)

    def plot_profile(self, ax=None):
        """Draw the pixelated fiber image with the peak fraction and ``n_pix``."""
        from .plotting import plot_profile
        return plot_profile(self, ax=ax)

    @classmethod
    def from_preset(cls, preset, **overrides):
        """Load a bundled preset (see :func:`list_instruments`), optionally overriding any field (``name`` included)."""
        path = CONFIG_DIR / f"{preset}.toml"
        if preset not in list_instruments():
            raise ValueError(f"Unknown instrument {preset!r}; available: {list_instruments()}")
        with open(path, "rb") as f:
            cfg = tomllib.load(f)["instrument"]
        cfg.update(overrides)
        return cls(**cfg)

    def replace(self, **changes):
        """Copy with some fields changed (``dataclasses.replace``)."""
        return dataclasses.replace(self, **changes)

    def qe_at(self, wavelength_nm):
        return self.qe(wavelength_nm) if callable(self.qe) else self.qe

    @property
    def efficiency(self):
        """Spectrograph efficiency including slit losses (``S_opt * slit_throughput``)."""
        return self.spectrograph_efficiency * self.slit_throughput

    def setting_window(self, center_nm):
        """``(lo, hi)`` [nm] covered by one setting centred on ``center_nm``.

        For an echelle with ``free_spectral_range_frac`` set (CSHELL) this is
        one free spectral range, ``center_nm * (1 -/+ frac / 2)``; otherwise
        the full instrument range. Pass the result as
        ``LFASTETC(..., wavelength_min_nm=lo, wavelength_max_nm=hi)``.
        """
        c = positive("center_nm", center_nm)
        if not self.wavelength_min_nm <= c <= self.wavelength_max_nm:
            raise ValueError(f"center_nm={c} is outside {self.wavelength_min_nm}-{self.wavelength_max_nm} nm")
        if self.free_spectral_range_frac is None:
            return self.wavelength_min_nm, self.wavelength_max_nm
        half = 0.5 * c * self.free_spectral_range_frac
        return max(c - half, self.wavelength_min_nm), min(c + half, self.wavelength_max_nm)

    def describe(self):
        return (f"{self.name}: R={self.resolving_power:.0f}, "
                f"{self.wavelength_min_nm:.0f}-{self.wavelength_max_nm:.0f} nm, "
                f"{self.pixels_per_resel}x{self.spatial_px:g} pix/resel ({self.profile}), S_opt={self.spectrograph_efficiency}, "
                f"slit={self.slit_throughput}, RN={self.read_noise} e-, dark={self.dark_current} e-/s/pix")
