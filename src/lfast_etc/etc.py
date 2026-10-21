"""The LFAST exposure time calculator.

Per-bin logic (identical to the reference LFAST ETC notebook, 2026-10-12):

* The wavelength grid has one point per resolution element,
  ``delta_lambda = lambda / R``.
* ``N_star`` is the total stellar electron rate [e-/s] inside one resolution
  element, summed over the whole array.
* ``N_sky`` is the sky electron rate per *pixel* [e-/s/pix].
* With ``n_pix`` pixels per resel, dark current ``D`` [e-/s/pix], read noise
  ``RN`` [e-/pix] and ``n_exp`` reads over a total time ``t`` (Howell 2006)::

      SNR = N_star t / sqrt( N_star t + n_pix (N_sky + D) t + n_pix n_exp RN^2 )

Telescope combination (:class:`~lfast_etc.telescope.Telescope`):

* ``combine="fiber"``: ``N`` fibers merge into one spectrograph input, so
  ``N_star`` and ``N_sky`` scale with ``N`` and the detector terms do not.
* ``combine="independent"``: each telescope yields its own spectrum; the
  detector noise terms also scale with ``N`` and ``SNR ~ sqrt(N)``.
"""
from dataclasses import dataclass, field
import numpy as np

from .instrument import Instrument
from .telescope import Telescope
from .sky import Sky
from .source import Source


@dataclass
class ETCResult:
    """Per-resolution-element output of an ETC calculation.

    Arrays are all on ``wavelength`` [nm]. ``exptime`` is the *total*
    exposure time [s] (sum of ``n_exp`` reads). Bins with zero throughput
    (telluric cores in ``K_opt``) are NaN in ``snr``/``exptime``.
    """
    wavelength: np.ndarray
    delta_lambda: np.ndarray
    snr: np.ndarray
    exptime: np.ndarray
    n_star: np.ndarray          # e-/s per resel, whole array
    n_sky_pix: np.ndarray       # e-/s per pixel (summed over detectors in independent mode)
    dark_pix: float             # e-/s per pixel (effective: N*D in independent mode)
    read_noise: float           # e-/pixel per read (effective: sqrt(N)*RN in independent mode)
    n_pix: float
    n_exp: int
    peak_e_per_pixel: np.ndarray   # brightest physical pixel per read (approx.)
    saturated: np.ndarray
    below_min_exptime: np.ndarray  # per-read time shorter than Instrument.min_exptime
    source: str = ""
    setup: str = ""
    solved_for: str = ""

    def at(self, wavelength_nm):
        """Interpolated ``snr``, ``exptime`` and rates at one wavelength [nm] (NaN outside the grid)."""
        w = self.wavelength
        return {k: float(np.interp(wavelength_nm, w, getattr(self, k), left=np.nan, right=np.nan))
                for k in ("snr", "exptime", "n_star", "n_sky_pix", "peak_e_per_pixel")}

    def summary(self, percentiles=(10, 50, 90)):
        """Percentiles of ``snr`` and ``exptime`` over the band, plus a saturation flag."""
        return {
            "snr": dict(zip(percentiles, np.nanpercentile(self.snr, percentiles))),
            "exptime": dict(zip(percentiles, np.nanpercentile(self.exptime, percentiles))),
            "saturated_fraction": float(np.mean(self.saturated)),
            "below_min_exptime_fraction": float(np.mean(self.below_min_exptime)),
            "n_bins": int(len(self.wavelength)),
        }

    def plot(self, ax=None, what=None, **kw):
        """Plot ``exptime`` (if solved for) or ``snr`` versus wavelength."""
        import matplotlib.pyplot as plt
        if ax is None:
            _, ax = plt.subplots(figsize=(9, 4))
        what = what or self.solved_for or "snr"
        y = getattr(self, what)
        kw.setdefault("lw", 0.8)
        ax.plot(self.wavelength, y, **kw)
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Exposure time (s)" if what == "exptime" else "SNR per resel")
        ax.set_yscale("log")
        ax.set_title(f"{self.source} | {self.setup.split(':')[0]}, {self.setup.split('|')[1].split('(')[0].strip()}", fontsize=9)
        return ax

    def plot_rates(self, ax=None):
        """Plot stellar, sky and dark electron rates per pixel versus wavelength
        (effective rates, summed over detectors in independent mode)."""
        import matplotlib.pyplot as plt
        if ax is None:
            _, ax = plt.subplots(figsize=(9, 4))
        ax.plot(self.wavelength, self.n_star / self.n_pix, lw=0.8, label="star (per pixel)")
        ax.plot(self.wavelength, self.n_sky_pix, lw=0.8, label="sky (per pixel)")
        ax.axhline(self.dark_pix, color="0.4", ls="--", label="dark")
        ax.set_yscale("log")
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("e$^-$ s$^{-1}$ pixel$^{-1}$")
        ax.legend()
        return ax


class LFASTETC:
    """Exposure time calculator for one instrument on the LFAST array.

    Parameters
    ----------
    instrument : str or Instrument
        Preset name (see :func:`~lfast_etc.list_instruments`) or an
        :class:`~lfast_etc.Instrument`.
    telescope : Telescope, optional
        Defaults to ``Telescope(n_telescopes=n_telescopes, fiber=fiber)``.
    sky : Sky, optional
        Defaults to the bundled SkyCalc table.
    n_telescopes, fiber :
        Shortcuts for building the default telescope.
    wavelength_min_nm, wavelength_max_nm : float, optional
        Restrict the grid (e.g. one CSHELL echelle setting). Defaults to the
        instrument bandpass clipped to the throughput/sky data range. For an
        echelle with a small free spectral range (CSHELL) the per-bin results
        are what *one setting* would give at that wavelength; covering the
        full range takes many settings, which are not summed.
    grid : {"resel", "uniform"}
        One bin per resolution element (default) or a uniform ``step_nm`` grid.
    apply_extinction : bool
        Multiply the source by the SkyCalc atmospheric transmission. Off by
        default because ``K_opt`` already includes extinction at the
        reference airmass.
    """

    def __init__(self, instrument="bench_r1800", telescope=None, sky=None,
                 n_telescopes=20, fiber="18um", wavelength_min_nm=None,
                 wavelength_max_nm=None, grid="resel", step_nm=1.0,
                 apply_extinction=False):
        self.instrument = (Instrument.from_preset(instrument)
                           if isinstance(instrument, str) else instrument)
        self.telescope = telescope or Telescope(n_telescopes=n_telescopes, fiber=fiber)
        self.sky = sky or Sky()
        self.apply_extinction = apply_extinction
        lo, hi = self.telescope.wavelength_range_nm
        slo, shi = float(self.sky.wavelength_nm.min()), float(self.sky.wavelength_nm.max())
        wmin = max(wavelength_min_nm or self.instrument.wavelength_min_nm, lo, slo)
        wmax = min(wavelength_max_nm or self.instrument.wavelength_max_nm, hi, shi)
        if wmax <= wmin:
            raise ValueError(f"empty wavelength range {wmin}-{wmax} nm")
        self.wavelength = self._build_grid(wmin, wmax, grid, step_nm)
        self.delta_lambda = self.wavelength / self.instrument.resolving_power
        # Everything below is fixed for the life of the object: build a new
        # LFASTETC rather than mutating instrument/telescope/sky attributes.
        self._throughput = self._compute_throughputs()
        self._sky_rate_pix = (self.sky.radiance_per_bin(self.wavelength, self.delta_lambda)
                              * self._throughput["sky"] / self.n_pix)

    # ---- grid ---------------------------------------------------------------
    def _build_grid(self, wmin, wmax, grid, step_nm):
        if grid == "uniform":
            # bins are still one resel wide (delta_lambda = lambda/R); step_nm only sets the sampling
            return np.arange(wmin, wmax + step_nm / 2, step_nm)
        if grid != "resel":
            raise ValueError("grid must be 'resel' or 'uniform'")
        # lambda_{j+1} = lambda_j (1 + 1/R)  -> geometric grid
        R = self.instrument.resolving_power
        n = int(np.floor(np.log(wmax / wmin) / np.log1p(1.0 / R))) + 1
        return wmin * (1.0 + 1.0 / R) ** np.arange(n)

    # ---- fixed (source-independent) pieces ----------------------------------
    @property
    def n_pix(self):
        return self.instrument.pixels_per_resel

    @property
    def _n_detector(self):
        """How many independent detector read-outs contribute noise."""
        return self.telescope.n_telescopes if self.telescope.combine == "independent" else 1

    def _compute_throughputs(self):
        wl = self.wavelength
        tel, ins = self.telescope, self.instrument
        qe = ins.qe_at(wl) * np.ones_like(wl)
        trans = self.sky.transmission(wl) if self.apply_extinction else 1.0
        common = tel.area_m2 * tel.n_telescopes * ins.efficiency * qe * self.delta_lambda
        return {
            "star": common * tel.throughput_point_source(wl) * trans,
            "sky": common * tel.throughput_extended(wl) * tel.fiber_solid_angle_arcsec2,
            "qe": qe,
        }

    def throughputs(self):
        """Per-bin conversion factors for the whole array, QE included.
        ``star``: photons/s/m^2/nm -> e-/s per resel [m^2 nm e-/ph].
        ``sky``: ph/s/m^2/nm/arcsec^2 -> e-/s per resel [m^2 nm arcsec^2 e-/ph]
        (the fiber solid angle is folded in)."""
        return self._throughput

    def sky_rate_per_pixel(self):
        """Sky electron rate per pixel [e-/s/pix], whole array."""
        return self._sky_rate_pix

    def star_rate(self, source):
        """Stellar electron rate per resel [e-/s], whole array."""
        return source.photon_flux(self.wavelength, self.delta_lambda) * self.throughputs()["star"]

    def rates(self, source):
        """All noise-model ingredients for ``source`` on this setup."""
        ins = self.instrument
        nd = self._n_detector
        return dict(
            n_star=self.star_rate(source),
            n_sky_pix=self.sky_rate_per_pixel(),
            dark_pix=ins.dark_current * nd,
            read_noise=ins.read_noise * np.sqrt(nd),
            n_pix=self.n_pix,
        )

    # ---- SNR <-> exposure time ---------------------------------------------
    @staticmethod
    def _snr_formula(n_star, n_sky_pix, dark, rn, n_pix, t, n_exp):
        """Howell SNR; NaN where the source throughput is zero."""
        var = n_star * t + n_pix * (n_sky_pix + dark) * t + n_pix * n_exp * rn ** 2
        with np.errstate(divide="ignore", invalid="ignore"):
            snr = np.where(var > 0, n_star * t / np.sqrt(var), 0.0)
        return np.where(n_star > 0, snr, np.nan)

    @staticmethod
    def _exptime_formula(n_star, n_sky_pix, dark, rn, n_pix, snr, n_exp):
        """Positive root of the SNR quadratic in t; NaN where the source throughput is zero."""
        s = np.where(n_star > 0, n_star, 1.0)
        b = n_pix * (n_sky_pix + dark)
        a = s ** 2
        bq = -(snr ** 2) * (s + b)
        c = -(snr ** 2) * n_pix * n_exp * rn ** 2
        t = (-bq + np.sqrt(bq ** 2 - 4 * a * c)) / (2 * a)
        return np.where(n_star > 0, t, np.nan)

    def _result(self, source, r, t, snr, n_exp, solved_for):
        t = np.broadcast_to(np.asarray(t, float), self.wavelength.shape).copy()
        snr = np.broadcast_to(np.asarray(snr, float), self.wavelength.shape).copy()
        # ponytail: peak pixel = uniform share of the resel on ONE physical
        # detector; a real LSF peaks higher. In independent mode the rates in
        # ``r`` are summed over the N detectors, so divide them back out.
        nd = self._n_detector
        peak = ((r["n_star"] / r["n_pix"] + r["n_sky_pix"]) / nd
                + self.instrument.dark_current) * t / n_exp
        with np.errstate(invalid="ignore"):
            saturated = peak > self.instrument.full_well
            below_min = (t / n_exp) < self.instrument.min_exptime
        return ETCResult(
            wavelength=self.wavelength, delta_lambda=self.delta_lambda, snr=snr, exptime=t,
            n_star=r["n_star"], n_sky_pix=r["n_sky_pix"], dark_pix=r["dark_pix"],
            read_noise=r["read_noise"], n_pix=r["n_pix"], n_exp=n_exp,
            peak_e_per_pixel=peak, saturated=saturated, below_min_exptime=below_min,
            source=source.describe(), setup=self.describe(), solved_for=solved_for)

    def snr(self, source, exptime, n_exp=1):
        """SNR per resel for a total exposure ``exptime`` [s] split into ``n_exp`` reads."""
        r = self.rates(source)
        s = self._snr_formula(r["n_star"], r["n_sky_pix"], r["dark_pix"], r["read_noise"],
                              r["n_pix"], np.asarray(exptime, float), n_exp)
        return self._result(source, r, exptime, s, n_exp, "snr")

    def exptime(self, source, snr, n_exp=1):
        """Total exposure time [s] to reach ``snr`` per resel in every bin."""
        r = self.rates(source)
        t = self._exptime_formula(r["n_star"], r["n_sky_pix"], r["dark_pix"], r["read_noise"],
                                  r["n_pix"], float(snr), n_exp)
        return self._result(source, r, t, snr, n_exp, "exptime")

    # ---- sensitivity maps ---------------------------------------------------
    def _select(self, arr, wavelength):
        """Reduce a per-bin array (last axis) to one number per leading index.

        ``None``: median over the band; float: nearest bin; ``(lo, hi)``:
        median over the window; ``"all"``: the full per-bin array.
        """
        w = self.wavelength
        if isinstance(wavelength, str):
            if wavelength == "all":
                return arr
            raise ValueError("wavelength must be None, a float, (lo, hi) or 'all'")
        if wavelength is None:
            return np.nanmedian(arr, axis=-1)
        if np.ndim(wavelength) == 0:
            return arr[..., np.argmin(np.abs(w - wavelength))]
        lo, hi = wavelength
        m = (w >= lo) & (w <= hi)
        if not m.any():
            raise ValueError(f"no bins in {lo}-{hi} nm")
        return np.nanmedian(arr[..., m], axis=-1)

    def snr_grid(self, source, mags, exptimes, wavelength=None, n_exp=1):
        """SNR per resel on a (magnitude x exposure time) grid.

        ``wavelength``: ``None`` -> median over the band; a float -> nearest
        bin; ``(lo, hi)`` -> median over that window; ``"all"`` -> per bin.
        Returns an array of shape ``(len(mags), len(exptimes))`` (plus a
        wavelength axis for ``"all"``).
        """
        mags = np.atleast_1d(np.asarray(mags, float))
        ts = np.atleast_1d(np.asarray(exptimes, float))
        r = self.rates(source.with_mag(0.0))
        scale = 10 ** (-0.4 * (mags - 0.0))
        s = r["n_star"][None, None, :] * scale[:, None, None]
        t = ts[None, :, None]
        snr = self._snr_formula(s, r["n_sky_pix"], r["dark_pix"], r["read_noise"],
                                r["n_pix"], t, n_exp)
        return self._select(snr, wavelength)

    def limiting_magnitude(self, source, exptime, snr, wavelength=None, n_exp=1):
        """Magnitude (in ``source.band``) reaching ``snr`` per resel in ``exptime`` [s].

        Solved analytically per bin, then reduced with ``wavelength`` as in
        :meth:`snr_grid` (``"all"`` returns the per-bin array).
        """
        r = self.rates(source.with_mag(0.0))
        t, q = float(exptime), float(snr)
        c = r["n_pix"] * (r["n_sky_pix"] + r["dark_pix"]) * t + r["n_pix"] * n_exp * r["read_noise"] ** 2
        # s^2 t^2 - q^2 s t - q^2 c = 0  ->  s = (q^2 t + sqrt(q^4 t^2 + 4 t^2 q^2 c)) / (2 t^2)
        s = (q ** 2 * t + np.sqrt(q ** 4 * t ** 2 + 4 * t ** 2 * q ** 2 * c)) / (2 * t ** 2)
        with np.errstate(divide="ignore", invalid="ignore"):
            m = np.where(r["n_star"] > 0, -2.5 * np.log10(s / np.where(r["n_star"] > 0, r["n_star"], 1.0)), np.nan)
        return self._select(m, wavelength)

    def exptime_vs_mag(self, source, mags, snr, wavelength=None, n_exp=1):
        """Total exposure time [s] to reach ``snr`` versus magnitude, reduced with ``wavelength`` as in :meth:`snr_grid`."""
        mags = np.atleast_1d(np.asarray(mags, float))
        r = self.rates(source.with_mag(0.0))
        s = r["n_star"][None, :] * 10 ** (-0.4 * mags)[:, None]
        t = self._exptime_formula(s, r["n_sky_pix"], r["dark_pix"], r["read_noise"],
                                  r["n_pix"], float(snr), n_exp)
        return self._select(t, wavelength)

    def reads_needed(self, source, exptime, margin=0.8):
        """Smallest ``n_exp`` keeping the brightest pixel below ``margin * full_well``
        in a total exposure ``exptime`` [s] (1 if nothing saturates)."""
        peak = self.snr(source, exptime, n_exp=1).peak_e_per_pixel
        return int(max(1, np.ceil(np.nanmax(peak) / (margin * self.instrument.full_well))))

    # ---- misc ---------------------------------------------------------------
    def describe(self):
        return f"{self.instrument.describe()} | {self.telescope.describe()} | {self.sky.describe()}"

    def __repr__(self):
        return f"LFASTETC({self.describe()}, {len(self.wavelength)} bins)"
