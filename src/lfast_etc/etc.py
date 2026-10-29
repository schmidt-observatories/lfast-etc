"""The LFAST exposure time calculator.

Per-bin logic (identical to the reference LFAST ETC notebook, 2026-10-12):

* The wavelength grid has one point per resolution element,
  ``delta_lambda = lambda / R``.
* ``N_star`` is the total stellar electron rate [e-/s] inside one resolution
  element, summed over the whole array.
* ``N_sky`` is the sky electron rate per *pixel* [e-/s/pix].
* With ``n_pix`` pixels per resel (dispersion width x fiber-image height), dark current ``D`` [e-/s/pix], read noise
  ``RN`` [e-/pix] and ``n_exp`` reads over a total time ``t`` (Howell 2006)::

      SNR = N_star t / sqrt( N_star t + n_pix (N_sky + D) t + n_pix n_exp RN^2 )

Telescope combination (:class:`~lfast_etc.telescope.Telescope`):

* ``combine="fiber"``: fibers merge in groups of ``telescopes_per_input``
  (20) into one spectrograph input each, so ``N_star`` and ``N_sky`` scale
  with ``N`` while the detector terms scale with the number of inputs
  ``ceil(N / 20)``; beyond one array ``SNR ~ sqrt(n_inputs)``.
* ``combine="independent"``: each telescope yields its own spectrum; the
  detector noise terms scale with ``N`` and ``SNR ~ sqrt(N)``.
"""
from dataclasses import dataclass, field
import warnings
import numpy as np

from .instrument import Instrument
from .telescope import Telescope
from .sky import Sky
from .source import Source
from .constants import positive

MAX_BINS = 5_000_000   # guard against absurd resolving powers (memory)


@dataclass(eq=False)
class ETCResult:
    """Per-resolution-element output of an ETC calculation.

    Arrays are all on ``wavelength`` [nm]. ``exptime`` is the *total*
    exposure time [s] (sum of ``n_exp`` reads). Bins with zero throughput
    (telluric cores in ``K_opt``, see ``Telescope.kopt``) are invalid:
    ``valid`` is False and ``snr``/``exptime`` are NaN there, whichever was
    solved for.
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
    full_well: float = float("inf")   # inf: unknown, saturation not assessed
    source: str = ""
    setup: str = ""
    title: str = ""
    solved_for: str = ""

    @property
    def valid(self):
        """Bins with non-zero source throughput (False in telluric cores)."""
        return self.n_star > 0

    def at(self, wavelength_nm):
        """``snr``, ``exptime`` and rates at one wavelength [nm].

        An exact bin centre returns that bin's values (NaN ``snr``/``exptime``
        if the bin is invalid). Anywhere else is a linear interpolation between
        the two bracketing bins, NaN if either is invalid or the point is
        outside the grid (#37).
        """
        w, x = self.wavelength, float(wavelength_nm)
        keys = ("snr", "exptime", "n_star", "n_sky_pix", "peak_e_per_pixel")
        i = int(np.argmin(np.abs(w - x)))
        if np.isclose(w[i], x, rtol=1e-12, atol=0.0):
            return {k: float(getattr(self, k)[i]) for k in keys}
        j = int(np.searchsorted(w, x))
        if j == 0 or j == len(w) or not (self.valid[j - 1] and self.valid[j]):
            return {k: np.nan for k in keys}
        f = (x - w[j - 1]) / (w[j] - w[j - 1])
        return {k: float((1 - f) * getattr(self, k)[j - 1] + f * getattr(self, k)[j]) for k in keys}

    def summary(self, percentiles=(10, 50, 90)):
        """Percentiles of ``snr`` and ``exptime`` over the band, plus a saturation flag
        (``saturated_fraction`` is ``None`` when the instrument's full well is unknown)."""
        def pct(a):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)       # all-NaN band -> NaN, quietly
                return {int(p): float(v) for p, v in zip(percentiles, np.nanpercentile(a, percentiles))}
        return {
            "snr": pct(self.snr),
            "exptime": pct(self.exptime),
            "saturated_fraction": float(np.mean(self.saturated)) if np.isfinite(self.full_well) else None,
            "below_min_exptime_fraction": float(np.mean(self.below_min_exptime)),
            "n_bins": int(len(self.wavelength)),
        }

    def plot(self, ax=None, what=None, **kw):
        """Plot ``exptime`` (if solved for) or ``snr`` versus wavelength."""
        import matplotlib.pyplot as plt
        if ax is None:
            _, ax = plt.subplots(figsize=(9, 4))
        labels = {"exptime": "Total exposure time (s)", "snr": "SNR per resel",
                  "n_star": "Star (e$^-$/s per resel)", "n_sky_pix": "Sky (e$^-$/s per pixel)",
                  "peak_e_per_pixel": "Peak e$^-$ per pixel per read"}
        what = what or self.solved_for or "snr"
        y = getattr(self, what)
        kw.setdefault("lw", 0.8)
        ax.plot(self.wavelength, y, **kw)
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel(labels.get(what, what))
        ax.set_yscale("log")
        ax.set_title(self.title, fontsize=9)
        return ax

    def plot_rates(self, ax=None):
        """Plot stellar, sky and dark electron rates per pixel versus wavelength
        (effective rates, summed over detectors in independent mode)."""
        import matplotlib.pyplot as plt
        if ax is None:
            _, ax = plt.subplots(figsize=(9, 4))
        ax.plot(self.wavelength, self.n_star / self.n_pix, lw=0.8, label="star (mean per pixel in aperture)")
        ax.plot(self.wavelength, self.n_sky_pix, lw=0.8, label="sky (per pixel)")
        ax.axhline(self.dark_pix, color="0.4", ls="--", label="dark")
        ax.set_yscale("log")
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("e$^-$ s$^{-1}$ pixel$^{-1}$")
        ax.legend()
        return ax

    def plot_noise_budget(self, ax=None):
        """Stacked fractions of the noise variance (star, sky, dark, read) per resel."""
        from .plotting import plot_noise_budget
        return plot_noise_budget(self, ax=ax)

    def simulate(self, seed=0):
        """One noisy realisation of the extracted spectrum [e- per resel]:
        ``N_* t (1 + g/SNR)`` with ``g`` standard normal; NaN where ``snr`` is NaN."""
        rng = np.random.default_rng(seed)
        signal = self.n_star * self.exptime
        with np.errstate(invalid="ignore", divide="ignore"):
            return signal * (1 + rng.standard_normal(signal.shape) / self.snr)

    def plot_spectrum(self, ax=None, seed=0, show_noiseless=True):
        """Plot one simulated realisation of the spectrum (see :meth:`simulate`)."""
        from .plotting import plot_spectrum
        return plot_spectrum(self, ax=ax, seed=seed, show_noiseless=show_noiseless)

    def to_table(self):
        """All per-bin arrays as an :class:`astropy.table.Table` (for CSV/FITS export)."""
        from astropy.table import Table
        cols = ("wavelength", "delta_lambda", "snr", "exptime", "n_star", "n_sky_pix",
                "peak_e_per_pixel", "saturated", "below_min_exptime")
        t = Table({k: getattr(self, k) for k in cols})
        t.meta.update(source=self.source, setup=self.setup, n_exp=self.n_exp, n_pix=self.n_pix,
                      dark_pix=self.dark_pix, read_noise=self.read_noise, solved_for=self.solved_for)
        return t


class LFASTETC:
    """Exposure time calculator for one instrument on the LFAST array.

    The object is immutable after construction (throughputs and the sky rate
    are precomputed): build a new ``LFASTETC`` instead of changing
    ``instrument``/``telescope``/``sky`` attributes.

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
        One bin per resolution element (default), or bins *sampled* every
        ``step_nm``. Bins are always one resel wide (``delta_lambda = lambda/R``),
        so a uniform grid oversamples (step < resel) or undersamples the band.
    apply_extinction : bool
        Multiply the source by the SkyCalc atmospheric transmission. Off by
        default because ``K_opt`` already includes extinction at the
        reference airmass.
    """

    def __init__(self, instrument="bench_r1800", telescope=None, sky=None,
                 n_telescopes=None, fiber=None, wavelength_min_nm=None,
                 wavelength_max_nm=None, grid="resel", step_nm=1.0,
                 apply_extinction=False):
        self.instrument = (Instrument.from_preset(instrument)
                           if isinstance(instrument, str) else instrument)
        if telescope is not None and (n_telescopes is not None or fiber is not None):
            raise ValueError("give either telescope= or n_telescopes=/fiber=, not both")
        if telescope is None:
            telescope = Telescope(n_telescopes=20 if n_telescopes is None else n_telescopes,
                                  fiber="18um" if fiber is None else fiber)
        self.telescope = telescope
        self.sky = sky or Sky()
        self.apply_extinction = apply_extinction
        lo, hi = self.telescope.wavelength_range_nm
        slo, shi = float(self.sky.wavelength_nm.min()), float(self.sky.wavelength_nm.max())
        wmin = self.instrument.wavelength_min_nm if wavelength_min_nm is None else positive("wavelength_min_nm", wavelength_min_nm)
        wmax = self.instrument.wavelength_max_nm if wavelength_max_nm is None else positive("wavelength_max_nm", wavelength_max_nm)
        ins = self.instrument
        wmin = max(wmin, lo, slo, ins.wavelength_min_nm)      # never evaluate a preset outside its bandpass
        wmax = min(wmax, hi, shi, ins.wavelength_max_nm)
        if wmax <= wmin:
            raise ValueError(f"empty wavelength range {wmin}-{wmax} nm (instrument covers "
                             f"{ins.wavelength_min_nm}-{ins.wavelength_max_nm} nm, data {max(lo, slo):.0f}-{min(hi, shi):.0f} nm)")
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
            step_nm = positive("step_nm", step_nm)
            n = int((wmax - wmin) / step_nm) + 1
        elif grid == "resel":
            # lambda_{j+1} = lambda_j (1 + 1/R)  -> geometric grid
            n = int(np.floor(np.log(wmax / wmin) / np.log1p(1.0 / self.instrument.resolving_power))) + 1
        else:
            raise ValueError("grid must be 'resel' or 'uniform'")
        if n > MAX_BINS:
            raise ValueError(f"grid would have {n} bins (> {MAX_BINS}); lower the resolving power or narrow the range")
        if grid == "uniform":
            return wmin + step_nm * np.arange(n)
        return wmin * (1.0 + 1.0 / self.instrument.resolving_power) ** np.arange(n)

    # ---- fixed (source-independent) pieces ----------------------------------
    @property
    def n_pix(self):
        """Extraction aperture per resel: dispersion width x fiber-image height (#13)."""
        return self.instrument.n_pix

    @property
    def _n_detector(self):
        """How many independent detector read-outs contribute noise."""
        return self.telescope.n_inputs

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
    # ``n_star > 0`` is the one validity rule: ``Telescope.kopt`` returns an
    # exact 0 below ``KOPT_EPS``, so the formulas never see a near-zero rate.
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
        with np.errstate(over="ignore"):             # overflow -> inf, flagged in _result
            t = (-bq + np.sqrt(bq ** 2 - 4 * a * c)) / (2 * a)
        return np.where(n_star > 0, t, np.nan)

    @staticmethod
    def _check_n_exp(n_exp):
        if int(n_exp) != n_exp or n_exp < 1:
            raise ValueError(f"n_exp must be a positive integer, got {n_exp!r}")
        return int(n_exp)

    def _result(self, source, r, t, snr, n_exp, solved_for):
        valid = r["n_star"] > 0
        t = np.where(valid, np.broadcast_to(np.asarray(t, float), self.wavelength.shape), np.nan)
        snr = np.where(valid, np.broadcast_to(np.asarray(snr, float), self.wavelength.shape), np.nan)
        # Peak pixel: the resel's star + sky light times the profile's peak
        # fraction (#19; 1/n_pix for profile="uniform"), image centred on a pixel
        # (the maximum for one resel; a continuum source sums neighbouring resels
        # along the dispersion, which is ~5% lower for the top hat and ~10% higher
        # for a Gaussian), on the most heavily
        # loaded detector. The rates in ``r`` are whole-array, so scale to the
        # brightest input's share of the telescopes (a full 20-group for N=21: #36).
        tel = self.telescope
        share = tel.telescopes_on_brightest_input / tel.n_telescopes
        peak = ((r["n_star"] + r["n_sky_pix"] * r["n_pix"]) * self.instrument.peak_fraction * share
                + self.instrument.dark_current) * t / n_exp
        with np.errstate(invalid="ignore"):
            saturated = peak > self.instrument.full_well
            below_min = (t / n_exp) < self.instrument.min_exptime
        if np.any((r["n_star"] > 0) & ~(np.isfinite(t) & np.isfinite(snr))):
            raise ValueError("numerical overflow: snr or exptime is out of float range")
        return ETCResult(
            wavelength=self.wavelength, delta_lambda=self.delta_lambda, snr=snr, exptime=t,
            n_star=r["n_star"], n_sky_pix=r["n_sky_pix"], dark_pix=r["dark_pix"],
            read_noise=r["read_noise"], n_pix=r["n_pix"], n_exp=n_exp,
            peak_e_per_pixel=peak, saturated=saturated, below_min_exptime=below_min,
            full_well=self.instrument.full_well, source=source.describe(), setup=self.describe(), solved_for=solved_for,
            title=f"{source.describe()} | {self.instrument.name}, {self.telescope.n_telescopes} telescopes ({self.telescope.combine})")

    def snr(self, source, exptime, n_exp=1):
        """SNR per resel for a total exposure ``exptime`` [s] split into ``n_exp`` reads."""
        n_exp = self._check_n_exp(n_exp)
        t = np.asarray(exptime, float)
        if np.any(t < 0) or not np.all(np.isfinite(t)):
            raise ValueError("exptime must be finite and >= 0")
        r = self.rates(source)
        s = self._snr_formula(r["n_star"], r["n_sky_pix"], r["dark_pix"], r["read_noise"],
                              r["n_pix"], np.asarray(exptime, float), n_exp)
        return self._result(source, r, exptime, s, n_exp, "snr")

    def exptime(self, source, snr, n_exp=1):
        """Total exposure time [s] to reach ``snr`` per resel in every bin."""
        n_exp = self._check_n_exp(n_exp)
        snr = positive("snr", snr, allow_zero=True)
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
        return self._reduce(arr[..., self._bin_mask(wavelength)], wavelength)

    def snr_grid(self, source, mags, exptimes, wavelength=None, n_exp=1):
        """SNR per resel on a (magnitude x exposure time) grid.

        ``wavelength``: ``None`` -> median over the band; a float -> nearest
        bin; ``(lo, hi)`` -> median over that window; ``"all"`` -> per bin.
        Returns an array of shape ``(len(mags), len(exptimes))`` (plus a
        wavelength axis for ``"all"``).
        """
        n_exp = self._check_n_exp(n_exp)
        mags = np.atleast_1d(np.asarray(mags, float))
        ts = np.atleast_1d(np.asarray(exptimes, float))
        r = self.rates(source.with_mag(0.0))
        # ponytail: select the wavelength bins first so the (mag, t, bin) cube stays small
        sel = self._bin_mask(wavelength)
        s = r["n_star"][sel][None, None, :] * 10 ** (-0.4 * mags)[:, None, None]
        t = ts[None, :, None]
        snr = self._snr_formula(s, r["n_sky_pix"][sel], r["dark_pix"], r["read_noise"],
                                r["n_pix"], t, n_exp)
        return self._reduce(snr, wavelength)

    def _bin_mask(self, wavelength):
        """Bins that :meth:`_select` would use for this ``wavelength`` selector."""
        w = self.wavelength
        if wavelength is None or (isinstance(wavelength, str) and wavelength == "all"):
            return slice(None)
        if isinstance(wavelength, str):
            raise ValueError("wavelength must be None, a float, (lo, hi) or 'all'")
        if np.ndim(wavelength) == 0:
            if not (w.min() <= wavelength <= w.max()):
                raise ValueError(f"wavelength {wavelength} nm is outside the grid {w.min():.1f}-{w.max():.1f} nm")
            i = int(np.argmin(np.abs(w - wavelength)))
            return slice(i, i + 1)
        if len(wavelength) != 2:
            raise ValueError("wavelength window must be (lo, hi)")
        m = (w >= wavelength[0]) & (w <= wavelength[1])
        if not m.any():
            raise ValueError(f"no bins in {wavelength[0]}-{wavelength[1]} nm")
        return m

    @staticmethod
    def _reduce(arr, wavelength):
        """Companion to :meth:`_bin_mask`: collapse the already-selected last axis."""
        if isinstance(wavelength, str) and wavelength == "all":
            return arr
        if wavelength is not None and np.ndim(wavelength) == 0:
            return arr[..., 0]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)           # all-NaN window -> NaN, quietly
            return np.nanmedian(arr, axis=-1)

    def limiting_magnitude(self, source, exptime, snr, wavelength=None, n_exp=1):
        """Magnitude (in ``source.band``) reaching ``snr`` per resel in ``exptime`` [s].

        Solved analytically per bin, then reduced with ``wavelength`` as in
        :meth:`snr_grid` (``"all"`` returns the per-bin array).
        """
        n_exp = self._check_n_exp(n_exp)
        r = self.rates(source.with_mag(0.0))
        t, q = positive("exptime", exptime), positive("snr", snr)
        c = r["n_pix"] * (r["n_sky_pix"] + r["dark_pix"]) * t + r["n_pix"] * n_exp * r["read_noise"] ** 2
        # s^2 t^2 - q^2 s t - q^2 c = 0  ->  s = (q^2 t + sqrt(q^4 t^2 + 4 t^2 q^2 c)) / (2 t^2)
        s = (q ** 2 * t + np.sqrt(q ** 4 * t ** 2 + 4 * t ** 2 * q ** 2 * c)) / (2 * t ** 2)
        with np.errstate(divide="ignore", invalid="ignore"):
            m = np.where(r["n_star"] > 0, -2.5 * np.log10(s / np.where(r["n_star"] > 0, r["n_star"], 1.0)), np.nan)
        return self._select(m, wavelength)

    def exptime_vs_mag(self, source, mags, snr, wavelength=None, n_exp=1):
        """Total exposure time [s] to reach ``snr`` versus magnitude, reduced with ``wavelength`` as in :meth:`snr_grid`."""
        n_exp = self._check_n_exp(n_exp)
        mags = np.atleast_1d(np.asarray(mags, float))
        r = self.rates(source.with_mag(0.0))
        sel = self._bin_mask(wavelength)
        s = r["n_star"][sel][None, :] * 10 ** (-0.4 * mags)[:, None]
        t = self._exptime_formula(s, r["n_sky_pix"][sel], r["dark_pix"], r["read_noise"],
                                  r["n_pix"], positive("snr", snr, allow_zero=True), n_exp)
        return self._reduce(t, wavelength)

    def reads_needed(self, source, exptime, margin=0.8):
        """Smallest ``n_exp`` keeping the brightest pixel below ``margin * full_well``
        (``0 < margin <= 1``) in a total exposure ``exptime`` [s]; 1 if nothing
        saturates or the full well is unknown."""
        margin = positive("margin", margin)
        if margin > 1:
            raise ValueError(f"margin is a fraction of the full well and must be <= 1, got {margin}")
        peak = self.snr(source, exptime, n_exp=1).peak_e_per_pixel
        if not np.isfinite(peak).any():                   # every bin invalid: nothing to saturate
            return 1
        return int(max(1, np.ceil(np.nanmax(peak) / (margin * self.instrument.full_well))))

    # ---- plots --------------------------------------------------------------
    def plot_throughput(self, ax=None):
        """Telescope, fiber, instrument and end-to-end throughput versus wavelength."""
        from .plotting import plot_throughput
        return plot_throughput(self, ax=ax)

    def plot_sky(self, ax=None):
        """Sky surface brightness per resel and the sky electron rate per pixel."""
        from .plotting import plot_sky
        return plot_sky(self, ax=ax)

    # ---- misc ---------------------------------------------------------------
    def describe(self):
        return f"{self.instrument.describe()} | {self.telescope.describe()} | {self.sky.describe()}"

    def __repr__(self):
        return f"LFASTETC({self.describe()}, {len(self.wavelength)} bins)"
