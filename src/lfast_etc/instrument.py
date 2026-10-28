"""Spectrograph + detector descriptions, loaded from TOML presets."""
import dataclasses
from dataclasses import dataclass
import tomllib

from .constants import CONFIG_DIR, positive


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
        Number of detector pixels in the extraction aperture of one
        resolution element (Howell's ``n_pix``). The bundled presets follow
        the reference notebook and count the dispersion direction only
        (2.2 / 2.5); multiply by the fiber-image height in pixels once that
        is known for 1x1 binning.
    spectrograph_efficiency : float
        Spectrograph optical throughput (``S_opt``), excluding slit losses and QE.
    slit_throughput : float
        Extra multiplicative throughput for slit modes (1 for the bare fiber).
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
                f"{self.pixels_per_resel} pix/resel, S_opt={self.spectrograph_efficiency}, "
                f"slit={self.slit_throughput}, RN={self.read_noise} e-, dark={self.dark_current} e-/s/pix")
