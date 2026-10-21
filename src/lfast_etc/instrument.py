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
        e-/pixel saturation level (used for the saturation flag only).
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
    full_well: float = float("inf")
    min_exptime: float = 0.0
    free_spectral_range_frac: float | None = None
    description: str = ""

    def __post_init__(self):
        for f in ("wavelength_min_nm", "wavelength_max_nm", "resolving_power", "pixels_per_resel",
                  "spectrograph_efficiency", "slit_throughput", "full_well"):
            setattr(self, f, positive(f, getattr(self, f)))
        for f in ("dark_current", "read_noise", "min_exptime"):
            setattr(self, f, positive(f, getattr(self, f), allow_zero=True))
        if not callable(self.qe):
            self.qe = positive("qe", self.qe)
        if self.wavelength_max_nm <= self.wavelength_min_nm:
            raise ValueError("wavelength_max_nm must exceed wavelength_min_nm")

    @classmethod
    def from_preset(cls, name, **overrides):
        """Load a bundled preset (see :func:`list_instruments`), optionally overriding fields."""
        path = CONFIG_DIR / f"{name}.toml"
        if name not in list_instruments():
            raise ValueError(f"Unknown instrument {name!r}; available: {list_instruments()}")
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

    def describe(self):
        return (f"{self.name}: R={self.resolving_power:.0f}, "
                f"{self.wavelength_min_nm:.0f}-{self.wavelength_max_nm:.0f} nm, "
                f"{self.pixels_per_resel} pix/resel, S_opt={self.spectrograph_efficiency}, "
                f"slit={self.slit_throughput}, RN={self.read_noise} e-, dark={self.dark_current} e-/s/pix")
