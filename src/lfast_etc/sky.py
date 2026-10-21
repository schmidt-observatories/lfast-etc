"""Sky background from an ESO SkyCalc table."""
from functools import lru_cache
from pathlib import Path
import numpy as np
from astropy.io import fits

from .constants import SKY_DIR

DEFAULT_SKY = SKY_DIR / "skytable_secz1.3.fits"


@lru_cache(maxsize=8)
def _load(path):
    with fits.open(path) as hdul:
        d = hdul[1].data
        lam = np.asarray(d["lam"], float)
        flux = np.asarray(d["flux"], float)
        trans = np.asarray(d["trans"], float)
    o = np.argsort(lam)
    return lam[o], flux[o], trans[o]


class Sky:
    """Sky radiance and transmission from a SkyCalc FITS table.

    The bundled table is **grey time**: the ESO SkyCalc model at 2640 m (no
    Kitt Peak site exists in SkyCalc), August, airmass 1.3, PWV 2.5 mm, moon
    up at 45 deg altitude and 45 deg from the target (about 60% of the
    visible sky), R=20,000 sampling. ``scale=0.4`` is a crude dark-time proxy. Resels
    narrower than the table sampling (R > 20,000, i.e. the CSHELL presets)
    fall back to the nearest sample, so airglow lines are pre-smoothed to
    R=20,000. Generate your own table at
    https://www.eso.org/observing/etc/skycalc/ for other conditions or a
    finer grid.

    Parameters
    ----------
    path : str or Path
        SkyCalc FITS file with columns ``lam`` [nm], ``flux``
        [ph/s/m^2/um/arcsec^2] and ``trans``.
    scale : float
        Multiplicative factor on the sky radiance (e.g. crude moon scaling).
    """

    def __init__(self, path=DEFAULT_SKY, scale=1.0):
        self.path = str(path)
        self.scale = float(scale)
        self.wavelength_nm, self._flux, self._trans = _load(self.path)

    def radiance_per_bin(self, wavelength_nm, delta_lambda_nm):
        """Mean sky radiance in each bin [ph/s/m^2/nm/arcsec^2].

        Every SkyCalc sample inside ``[lambda - dl/2, lambda + dl/2]`` is
        averaged; a bin narrower than the table sampling falls back to the
        nearest sample.
        """
        wl = np.asarray(wavelength_nm, float)
        half = np.asarray(delta_lambda_nm, float) / 2.0
        lam, flux = self.wavelength_nm, self._flux
        lo = np.searchsorted(lam, wl - half, side="left")
        hi = np.searchsorted(lam, wl + half, side="right")
        csum = np.concatenate([[0.0], np.cumsum(flux)])
        n = hi - lo
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = (csum[hi] - csum[lo]) / n
        nearest = np.clip(np.searchsorted(lam, wl), 1, len(lam) - 1)
        use_prev = np.abs(wl - lam[nearest - 1]) <= np.abs(lam[nearest] - wl)
        nearest = np.where(use_prev, nearest - 1, nearest)
        out = np.where(n > 0, mean, flux[nearest])
        return out / 1000.0 * self.scale   # per micron -> per nm

    def transmission(self, wavelength_nm):
        """Atmospheric transmission interpolated onto ``wavelength_nm``."""
        return np.interp(wavelength_nm, self.wavelength_nm, self._trans)

    def describe(self):
        return f"SkyCalc table {Path(self.path).name} x {self.scale:g}"
