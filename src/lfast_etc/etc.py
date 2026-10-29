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

from .instrument import Instrument, BLUR_SIGMA_PX, peak_pixel_fraction
from .telescope import Telescope, KOPT_EPS
from .sky import Sky
from .source import Source
from .constants import positive

MAX_BINS = 5_000_000   # guard against absurd resolving powers (memory)
LSFS = ("gaussian", "tophat")


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
    pixel: "PixelResult | None" = None   # per-detector-pixel view (``lsf`` runs only)

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


@dataclass(eq=False)
class PixelResult:
    """Per-detector-pixel companion of :class:`ETCResult` (``result.pixel``), along the dispersion.

    One entry per pixel of the mode (``pixels_per_resel`` per resolution element, on the same
    geometric grid as the resels, from the same first wavelength). Rates are per *extracted*
    pixel with the ``spatial_px`` pixels of the aperture summed in: ``n_star`` [e-/s] is the
    source in that pixel, ``n_sky_pix`` the sky per physical pixel (as in ``ETCResult``), so one
    extracted pixel's variance is ``n_star t + spatial_px (n_sky_pix + dark_pix) t + spatial_px
    n_exp RN^2``. ``exptime`` is the per-resel solve carried to the pixels of each resel;
    ``snr`` is per extracted pixel, which for a continuum source equals ``ETCResult.snr /
    sqrt(pixels_per_resel)``. Pixels whose resolution element has zero throughput are invalid
    (``n_star`` 0, ``snr``/``exptime`` NaN), as for the resels.
    """
    wavelength: np.ndarray
    delta_lambda: np.ndarray
    snr: np.ndarray
    exptime: np.ndarray
    n_star: np.ndarray          # e-/s per extracted pixel, whole array
    n_sky_pix: np.ndarray       # e-/s per physical pixel
    dark_pix: float
    read_noise: float
    spatial_px: float
    pixels_per_resel: float
    n_exp: int
    peak_e_per_pixel: np.ndarray
    saturated: np.ndarray
    full_well: float = float("inf")

    @property
    def valid(self):
        return self.n_star > 0

    def variance(self):
        """Noise variance of each extracted pixel [e-^2] for ``exptime``."""
        t = self.exptime
        return self.n_star * t + self.spatial_px * ((self.n_sky_pix + self.dark_pix) * t + self.n_exp * self.read_noise ** 2)

    def simulate(self, seed=0):
        """One noisy realisation [e- per extracted pixel]: the signal plus a Gaussian draw of
        the pixel's own variance (source, sky, dark and read noise); NaN where invalid."""
        rng = np.random.default_rng(seed)
        signal = self.n_star * self.exptime
        with np.errstate(invalid="ignore"):
            return np.where(self.valid, signal + rng.standard_normal(signal.shape) * np.sqrt(self.variance()), np.nan)

    def to_table(self):
        """All per-pixel arrays as an :class:`astropy.table.Table`."""
        from astropy.table import Table
        cols = ("wavelength", "delta_lambda", "snr", "exptime", "n_star", "n_sky_pix", "peak_e_per_pixel", "saturated")
        t = Table({k: getattr(self, k) for k in cols})
        t.meta.update(n_exp=self.n_exp, spatial_px=self.spatial_px, pixels_per_resel=self.pixels_per_resel,
                      dark_pix=self.dark_pix, read_noise=self.read_noise)
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
        reference airmass. Ignored with an ``lsf``, which always applies the
        resolved transmission once (see below).
    lsf : {"gaussian", "tophat", None}
        Line-spread function. With one, the source, the sky and the atmospheric
        transmission are evaluated on a fine working grid (``oversample``
        samples per detector pixel), the source and sky electron-rate densities
        are convolved with the LSF (``"gaussian"``: FWHM one resel, lambda/R;
        ``"tophat"``: the fiber image, a top hat ``pixels_per_resel`` wide blurred
        by ``BLUR_SIGMA_PX``), and integrated over each detector pixel
        (``result.pixel``) and over each resolution element (the per-resel API,
        which is thus the sum of its pixels). Narrow lines and airglow then have
        the instrument's width instead of being averaged into one bin. ``K_opt``
        carries a 1 nm-sampled atmosphere: on the fine grid it is divided by the
        SkyCalc transmission box-averaged over 1 nm (floored at 0.2) and
        multiplied by the resolved transmission, so telluric lines are resolved
        too (within the sky table's own sampling). ``None`` is the reference
        notebook's boxcar per resel with no pixel view; a ``"uniform"`` grid
        is always that.
    oversample : int
        Fine-grid samples per detector pixel for the LSF path (default 8).
    """

    def __init__(self, instrument="bench_r1800", telescope=None, sky=None,
                 n_telescopes=None, fiber=None, wavelength_min_nm=None,
                 wavelength_max_nm=None, grid="resel", step_nm=1.0,
                 apply_extinction=False, lsf="gaussian", oversample=8):
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
        if lsf is not None and lsf not in LSFS:
            raise ValueError(f"lsf must be one of {LSFS} or None, got {lsf!r}")
        self.lsf = lsf if grid == "resel" else None       # pixels are laid out from the resel grid: a uniform grid is boxcar-only
        self.oversample = None
        if self.lsf is not None:
            self.oversample = int(positive("oversample", oversample))
            self._lsf_setup()

    # ---- line-spread function and pixel grid --------------------------------
    def _lsf_setup(self):
        """Pixel grid, fine working grid, LSF kernel and the source-independent pieces."""
        ins, tel, sky = self.instrument, self.telescope, self.sky
        R, ppr, os_ = ins.resolving_power, ins.pixels_per_resel, self.oversample
        w0, w1 = float(self.wavelength[0]), float(self.wavelength[-1])
        rp = 1.0 + 1.0 / (R * ppr)                               # one detector pixel, geometric like the resels
        n_px = int(np.floor(np.log(w1 / w0) / np.log(rp))) + 1
        self.pixel_wavelength = w0 * rp ** np.arange(n_px)
        self.pixel_delta_lambda = self.pixel_wavelength * (np.sqrt(rp) - 1.0 / np.sqrt(rp))
        self._px_edges = (self.pixel_wavelength / np.sqrt(rp), self.pixel_wavelength * np.sqrt(rp))
        self._resel_edges = (self.wavelength - self.delta_lambda / 2.0, self.wavelength + self.delta_lambda / 2.0)
        # fine grid: os_ samples per pixel, padded so the kernel never sees an edge
        sigma = ppr * os_ / 2.3548                                # Gaussian LSF, FWHM = one resel [fine samples]
        blur = BLUR_SIGMA_PX * os_
        pad = int(np.ceil(4 * sigma + ppr * os_ / 2 + 4 * blur)) + 1
        n_fine = n_px * os_ + 2 * pad + 1
        if n_fine > MAX_BINS:
            raise ValueError(f"fine grid would have {n_fine} samples (> {MAX_BINS}); lower oversample or narrow the range")
        rf = rp ** (1.0 / os_)
        lf = w0 * rf ** np.arange(-pad, n_px * os_ + pad + 1)
        self._fine_w, self._fine_dl = lf, lf * (np.sqrt(rf) - 1.0 / np.sqrt(rf))
        x = np.arange(-pad, pad + 1, dtype=float)
        if self.lsf == "gaussian":
            k = np.exp(-0.5 * (x / sigma) ** 2)
        else:                                                     # fiber image: top hat ppr pixels wide, camera blur
            box = np.clip(ppr * os_ / 2.0 - np.abs(x) + 0.5, 0.0, 1.0)   # fractional edge samples
            g = np.exp(-0.5 * (x / blur) ** 2)
            k = np.convolve(box, g, mode="same")
        self._kernel = k / k.sum()
        # throughput on the fine grid: K_opt's 1 nm atmosphere out, the resolved transmission in (#32)
        t_1nm = np.clip(sky.transmission_smoothed(lf, 1.0), 0.2, 1.0)
        t_res = sky.transmission(lf)
        common = tel.area_m2 * tel.n_telescopes * ins.efficiency * ins.qe_at(lf) * np.ones_like(lf)
        k_eff = tel.throughput_point_source(lf) / t_1nm * t_res                              # resolved point-source throughput
        self._fine_clear = (tel.kopt(lf) > 0) & (k_eff >= KOPT_EPS)                          # kopt's rule, on the resolved curve
        self._fine_star_factor = common * k_eff                                               # e-/s/nm per ph/s/m^2/nm
        sky_density = common * tel.throughput_extended(lf) / t_1nm * t_res * tel.fiber_solid_angle_arcsec2 \
            * sky.radiance_per_bin(lf, self._fine_dl)                                       # e-/s/nm, whole array
        sky_c = self._convolve(sky_density)
        # a resel or pixel is valid only if every fine sample in it has throughput: the same
        # neighbour-zeroing as Telescope.kopt's rule (#12), so no bin lives on a 1e-6 wing
        clear = self._fine_clear.astype(float)
        self._resel_clear = self._integrate(clear, *self._resel_edges) / self.delta_lambda > 0.999
        w, wp = self.wavelength, self.pixel_wavelength
        j = np.clip(np.searchsorted(w, wp), 1, len(w) - 1)              # the resel each pixel belongs to
        self._px_resel = np.where(wp - w[j - 1] <= w[j] - wp, j - 1, j)
        self._px_clear = (self._integrate(clear, *self._px_edges) / self.pixel_delta_lambda > 0.999) & self._resel_clear[self._px_resel]
        self._sky_resel_pix = self._integrate(sky_c, *self._resel_edges) / self.n_pix
        self._sky_px_pix = self._integrate(sky_c, *self._px_edges) / ins.spatial_px

    def _convolve(self, density):
        return np.convolve(density, self._kernel, mode="same")

    def _integrate(self, density, lo, hi):
        """Integral of a fine-grid density over ``[lo, hi]`` (trapezoid, via the cumulative integral)."""
        lf = self._fine_w
        cum = np.concatenate([[0.0], np.cumsum(0.5 * (density[1:] + density[:-1]) * np.diff(lf))])
        return np.interp(hi, lf, cum) - np.interp(lo, lf, cum)

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
        """Sky electron rate per pixel [e-/s/pix], whole array (LSF-convolved with an ``lsf``)."""
        return self._sky_rate_pix if self.lsf is None else self._sky_resel_pix

    def star_rate(self, source):
        """Stellar electron rate per resel [e-/s], whole array."""
        if self.lsf is not None:
            return self.rates(source)["n_star"]
        return source.photon_flux(self.wavelength, self.delta_lambda) * self.throughputs()["star"]

    def rates(self, source):
        """All noise-model ingredients for ``source`` on this setup. With an ``lsf`` the
        ``"pixel"`` entry carries the per-extracted-pixel ``n_star`` and per-physical-pixel
        ``n_sky_pix`` on ``pixel_wavelength``."""
        ins = self.instrument
        nd = self._n_detector
        r = dict(dark_pix=ins.dark_current * nd, read_noise=ins.read_noise * np.sqrt(nd), n_pix=self.n_pix)
        if self.lsf is None:
            return dict(n_star=source.photon_flux(self.wavelength, self.delta_lambda) * self.throughputs()["star"],
                        n_sky_pix=self._sky_rate_pix, **r)
        star = self._convolve(source.photon_flux(self._fine_w, self._fine_dl) * self._fine_star_factor)
        return dict(n_star=np.where(self._resel_clear, self._integrate(star, *self._resel_edges), 0.0),
                    n_sky_pix=self._sky_resel_pix,
                    pixel=dict(n_star=np.where(self._px_clear, self._integrate(star, *self._px_edges), 0.0),
                               n_sky_pix=self._sky_px_pix), **r)

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
        pixel = None
        if "pixel" in r:
            ins, rp, wp = self.instrument, r["pixel"], self.pixel_wavelength
            sp, j = ins.spatial_px, self._px_resel
            valid_px = rp["n_star"] > 0                                      # pixels follow their resel's validity
            t_px = np.where(valid_px, t[j], np.nan)
            snr_px = self._snr_formula(rp["n_star"], rp["n_sky_pix"], r["dark_pix"], r["read_noise"], sp,
                                       np.nan_to_num(t_px), n_exp)
            snr_px = np.where(np.isfinite(t_px), snr_px, np.nan)
            f_sp = peak_pixel_fraction(1.0, sp, ins.profile)                  # brightest pixel of one extracted column
            peak_px = ((rp["n_star"] + rp["n_sky_pix"] * sp) * f_sp * share + ins.dark_current) * t_px / n_exp
            with np.errstate(invalid="ignore"):
                sat_px = peak_px > ins.full_well
            pixel = PixelResult(wavelength=wp, delta_lambda=self.pixel_delta_lambda, snr=snr_px, exptime=t_px,
                                n_star=rp["n_star"], n_sky_pix=rp["n_sky_pix"], dark_pix=r["dark_pix"],
                                read_noise=r["read_noise"], spatial_px=sp, pixels_per_resel=ins.pixels_per_resel,
                                n_exp=n_exp, peak_e_per_pixel=peak_px, saturated=sat_px, full_well=ins.full_well)
        return ETCResult(
            wavelength=self.wavelength, delta_lambda=self.delta_lambda, snr=snr, exptime=t,
            n_star=r["n_star"], n_sky_pix=r["n_sky_pix"], dark_pix=r["dark_pix"],
            read_noise=r["read_noise"], n_pix=r["n_pix"], n_exp=n_exp,
            peak_e_per_pixel=peak, saturated=saturated, below_min_exptime=below_min,
            full_well=self.instrument.full_well, source=source.describe(), setup=self.describe(), solved_for=solved_for,
            title=f"{source.describe()} | {self.instrument.name}, {self.telescope.n_telescopes} telescopes ({self.telescope.combine})",
            pixel=pixel)

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

    # ---- throughput cascade ---------------------------------------------------
    def throughput_cascade(self, path="star"):
        """``Telescope.throughput_cascade`` on this grid, extended with the spectrograph optics,
        slit and detector QE; the last cumulative entry is the end-to-end throughput."""
        ins = self.instrument
        c = self.telescope.throughput_cascade(self.wavelength, path=path, sky=self.sky)
        w = self.wavelength
        extra = [dict(name=f"spectrograph optics S_opt ({ins.name})", transmission=np.full_like(w, ins.spectrograph_efficiency), source="instrument preset")]
        if ins.slit_throughput != 1.0:
            extra.append(dict(name="slit throughput", transmission=np.full_like(w, ins.slit_throughput), source="instrument preset"))
        extra.append(dict(name="detector QE", transmission=ins.qe_at(w) * np.ones_like(w), source="instrument preset"))
        run = c["cumulative"][-1]
        for s in extra:
            run = run * s["transmission"]; c["stages"].append(s); c["cumulative"].append(run)
        return c

    def plot_throughput_cascade(self, ax=None, path="star"):
        """Cumulative throughput after each stage of the Choi et al. (2026) budget plus the instrument."""
        from .plotting import plot_throughput_cascade
        return plot_throughput_cascade(self, ax=ax, path=path)

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
        lsf = "boxcar per resel (no LSF)" if self.lsf is None else f"{self.lsf} LSF, {self.oversample}x oversampled pixels"
        return f"{self.instrument.describe()} | {self.telescope.describe()} | {self.sky.describe()} | {lsf}"

    def __repr__(self):
        return f"LFASTETC({self.describe()}, {len(self.wavelength)} bins)"
