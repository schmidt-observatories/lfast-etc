"""Plot helpers: sensitivity maps, noise budgets, throughput, sky, simulated spectra.

Every function takes an optional ``ax`` and returns it, so the plots compose
with ``plt.subplots``. The ``ETCResult``, ``LFASTETC`` and ``Source`` classes
expose the relevant ones as methods (``res.plot_noise_budget()`` etc.).
"""
import numpy as np
import matplotlib.pyplot as plt

from .constants import H_PLANCK, C_LIGHT

UNIT_DIAMETER_M = 0.726


def _ax(ax, figsize):
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    return ax


def plot_noise_budget(res, ax=None):
    """Stacked fractions of the noise variance (star, sky, dark, read) per resel.

    Uses the exposure time stored in ``res`` (per bin), so it answers "what
    limits me at this wavelength for this observation".
    """
    ax = _ax(ax, (9, 3.6))
    t = res.exptime
    terms = {
        "star": res.n_star * t,
        "sky": res.n_pix * res.n_sky_pix * t,
        "dark": res.n_pix * res.dark_pix * t * np.ones_like(t),
        "read": res.n_pix * res.n_exp * res.read_noise ** 2 * np.ones_like(t),
    }
    total = sum(terms.values())
    with np.errstate(invalid="ignore", divide="ignore"):
        fracs = [v / total for v in terms.values()]
    ax.stackplot(res.wavelength, *fracs, labels=list(terms), alpha=0.85, lw=0)
    ax.set_ylim(0, 1)
    ax.set_xlim(res.wavelength.min(), res.wavelength.max())
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Fraction of noise variance")
    ax.legend(loc="upper right", ncol=4, fontsize=8)
    ax.set_title(res.title, fontsize=9)
    return ax


def plot_spectrum(res, ax=None, seed=0, show_noiseless=True):
    """One noisy realisation of the extracted spectrum (electrons per resel).

    The realisation is ``N_* t (1 + g / SNR)`` with ``g`` standard normal, i.e.
    the per-resel SNR of ``res`` applied to the noiseless signal.
    """
    ax = _ax(ax, (9, 3.6))
    signal = res.n_star * res.exptime
    noisy = res.simulate(seed=seed)
    ax.plot(res.wavelength, noisy, lw=0.6, color="0.2", label=f"simulated (seed {seed})")
    if show_noiseless:
        ax.plot(res.wavelength, signal, lw=1.2, label="noiseless")
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("e$^-$ per resel")
    ax.legend(fontsize=8)
    ax.set_title(res.title, fontsize=9)
    return ax


def plot_throughput(etc, ax=None):
    """Telescope ``K_opt``, fiber injection, instrument (``S_opt x QE``) and the end-to-end product."""
    ax = _ax(ax, (9, 3.6))
    w = etc.wavelength
    tel, ins = etc.telescope, etc.instrument
    kopt = tel.kopt(w)
    instr = ins.efficiency * ins.qe_at(w) * np.ones_like(w)
    ax.plot(w, kopt, label="telescope + fiber, point source ($K_{opt}$)")
    ax.plot(w, tel.eta_fib(w), ls="--", label=r"fiber injection ($\eta_{fib}$, inside $K_{opt}$)")
    ax.plot(w, instr, ls=":", label=f"spectrograph x QE ({ins.name})")
    ax.plot(w, kopt * instr, lw=2.2, color="k", label="end to end")
    ax.set_ylim(0, 1.4 * max(np.nanmax(kopt), np.nanmax(tel.eta_fib(w)), np.nanmax(instr)))
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Throughput")
    ax.legend(fontsize=8, loc="upper left", ncol=2)
    return ax


def sky_surface_brightness_ab(etc):
    """Resel-averaged sky surface brightness on the ETC grid [AB mag / arcsec^2]."""
    rad = etc.sky.radiance_per_bin(etc.wavelength, etc.delta_lambda)      # ph/s/m^2/nm/arcsec^2
    lam_m = etc.wavelength * 1e-9
    f_lambda = rad * H_PLANCK * C_LIGHT / lam_m * 1e9                       # W/m^2/m/arcsec^2
    f_nu_jy = f_lambda * lam_m ** 2 / C_LIGHT / 1e-26                       # Jy/arcsec^2
    with np.errstate(divide="ignore"):
        return -2.5 * np.log10(f_nu_jy / 3631.0)


def plot_sky(etc, ax=None):
    """Sky surface brightness per resel (AB mag/arcsec^2) and the fiber's sky electron rate per pixel."""
    ax = _ax(ax, (9, 3.6))
    ax.plot(etc.wavelength, sky_surface_brightness_ab(etc), lw=0.8)
    ax.invert_yaxis()
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Sky (AB mag arcsec$^{-2}$ per resel)")
    ax2 = ax.twinx()
    ax2.plot(etc.wavelength, etc.sky_rate_per_pixel(), lw=0.8, color="0.5", alpha=0.7)
    ax2.set_yscale("log")
    ax2.set_ylabel("Sky e$^-$ s$^{-1}$ pixel$^{-1}$ (grey)", color="0.4")
    ax.set_title(etc.sky.describe() + f" | {etc.telescope.fiber} fiber", fontsize=9)
    return ax


def plot_source(source, wavelength_nm=None, ax=None, photon=False):
    """The source SED (``F_lambda`` or photon flux) with the normalising band shaded."""
    from .source import BANDS
    ax = _ax(ax, (9, 3.6))
    if wavelength_nm is None:
        wavelength_nm = np.linspace(350, 1800, 4000)
    w = np.asarray(wavelength_nm, float)
    y = source.photon_flux(w) if photon else source.flux(w)
    ax.plot(w, y, lw=1)
    c, width, _ = BANDS[source.band]
    ax.axvspan(c - width / 2, c + width / 2, color="0.85", zorder=0,
               label=f"{source.band} band ({source.band}={source.mag:.2f} {source.system})")
    for center, depth, _ in source.lines:
        ax.axvline(center, color="0.5", ls=":", lw=0.8)
    ax.set_yscale("log")
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("photons s$^{-1}$ m$^{-2}$ nm$^{-1}$" if photon else "$F_\\lambda$ (W m$^{-2}$ nm$^{-1}$)")
    ax.legend(fontsize=8)
    ax.set_title(source.describe(), fontsize=9)
    return ax


def equivalent_diameter_axis(ax, diameter_m=UNIT_DIAMETER_M, ticks=(1, 2, 3, 5, 10, 20, 30)):
    """Top axis in single-telescope equivalent diameter for an x axis in number of telescopes."""
    from matplotlib.ticker import FixedLocator, FixedFormatter, NullFormatter
    sec = ax.secondary_xaxis("top", functions=(lambda n: diameter_m * np.sqrt(n),
                                               lambda d: (d / diameter_m) ** 2))
    lo, hi = ax.get_xlim()
    ticks = [t for t in ticks if lo <= (t / diameter_m) ** 2 <= hi]
    sec.xaxis.set_major_locator(FixedLocator(ticks))
    sec.xaxis.set_major_formatter(FixedFormatter([f"{t:g}" for t in ticks]))
    sec.xaxis.set_minor_formatter(NullFormatter())
    sec.set_xlabel("Equivalent diameter (m)")
    return sec


def plot_vs_n_telescopes(make_etc, source, ns, exptimes, snr, wavelength=None, ax=None,
                         nominal=20, label_fmt="{t:g} s"):
    """Limiting magnitude versus number of telescopes for several exposure times.

    ``make_etc(n)`` must return the :class:`~lfast_etc.LFASTETC` for ``n``
    telescopes (so the caller decides fiber, combiner and bandpass).
    """
    ax = _ax(ax, (7, 4.5))
    ns = np.asarray(ns, int)
    for t in np.atleast_1d(exptimes):
        lim = [make_etc(int(n)).limiting_magnitude(source, t, snr, wavelength=wavelength) for n in ns]
        ax.plot(ns, lim, marker="o", ms=3, label=label_fmt.format(t=t))
    ax.set_xscale("log")
    ax.invert_yaxis()
    if nominal is not None:
        ax.axvline(nominal, color="0.4", ls="--", lw=0.8)
    ax.set_xlabel("Number of telescopes")
    ax.set_ylabel(f"Limiting {source.band} mag (SNR={snr:g}/resel)")
    ax.legend(title="exposure", fontsize=8)
    equivalent_diameter_axis(ax)
    return ax


def plot_snr_map(mags, exptimes, snr, ax=None, levels=(5, 10, 30, 100), band="V",
                 cmap="viridis", title=None):
    """Colour map of SNR over (exposure time, magnitude) with labelled SNR contours.

    ``snr`` has shape ``(len(mags), len(exptimes))`` as returned by
    :meth:`~lfast_etc.LFASTETC.snr_grid`.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 5))
    snr = np.asarray(snr)
    pcm = ax.pcolormesh(exptimes, mags, np.log10(np.maximum(snr, 1e-3)), cmap=cmap, shading="auto")
    cs = ax.contour(exptimes, mags, snr, levels=levels, colors="w", linewidths=1)
    ax.clabel(cs, fmt="SNR=%g", fontsize=8)
    ax.set_xscale("log")
    ax.invert_yaxis()
    ax.set_xlabel("Total exposure time (s)")
    ax.set_ylabel(f"{band} magnitude")
    plt.colorbar(pcm, ax=ax, label="log$_{10}$ SNR per resel")
    if title:
        ax.set_title(title, fontsize=9)
    return ax


def plot_limiting_magnitude(etc, source, exptimes, snr, ax=None, label_fmt="{t:g} s"):
    """Limiting magnitude versus wavelength for several exposure times."""
    if ax is None:
        _, ax = plt.subplots(figsize=(9, 4))
    for t in np.atleast_1d(exptimes):
        m = etc.limiting_magnitude(source, t, snr, wavelength="all")
        ax.plot(etc.wavelength, m, lw=0.8, label=label_fmt.format(t=t))
    ax.invert_yaxis()
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel(f"Limiting {source.band} mag (SNR={snr:g}/resel)")
    ax.legend(title="exposure")
    return ax


def plot_profile(instrument, ax=None):
    """Pixelated fiber image of one resel on the detector (#19) with the peak fraction."""
    ax = _ax(ax, (4.2, 4))
    img = instrument.pixel_image()
    n = img.shape[0]
    ext = (-n / 2 - 0.5, n / 2 - 0.5, -n / 2 - 0.5, n / 2 - 0.5)      # pixel k spans [k-0.5, k+0.5)
    im = ax.imshow(img, origin="lower", extent=ext, cmap="viridis", interpolation="nearest")
    if n > 1:
        for k in np.arange(-n / 2 - 0.5, n / 2 + 0.5):
            ax.axhline(k, color="w", lw=0.3, alpha=0.5); ax.axvline(k, color="w", lw=0.3, alpha=0.5)
    ax.set_xlabel("dispersion (pixels)"); ax.set_ylabel("spatial (pixels)")
    ax.set_title(f"{instrument.name}: {instrument.profile}, {instrument.pixels_per_resel:g} x "
                 f"{instrument.spatial_px:g} px, n_pix={instrument.n_pix:.3g}\npeak pixel holds "
                 f"{100 * instrument.peak_fraction:.1f}% of the resel", fontsize=9)
    plt.colorbar(im, ax=ax, fraction=0.046, label="fraction of resel light per pixel")
    return ax


def plot_fiber_coupling(telescope, ax=None, seeing=None):
    """Geometric fraction of a point source inside the fiber versus seeing FWHM at
    500 nm (#14), for both PSF shapes, with the reference and current seeing marked."""
    ax = _ax(ax, (6.5, 4))
    seeing = np.linspace(0.4, 3.0, 131) if seeing is None else np.asarray(seeing, float)
    for psf, ls in (("moffat", "-"), ("gaussian", "--")):
        label = f"{psf} disc (beta={telescope.moffat_beta:g})" if psf == "moffat" else "gaussian disc"
        ax.plot(seeing, telescope.fiber_coupling(seeing, psf=psf), ls, label=label)
    try:                                       # Choi-calibrated: eta_fib at 500 nm times the table factor
        import dataclasses
        eta0 = float(telescope.eta_fib(500.0))
        f = np.array([dataclasses.replace(telescope, seeing_arcsec=float(x), seeing_psf="choi", seeing_ref_arcsec=1.0).seeing_factor(500.0) for x in seeing])
        ax.plot(seeing, eta0 * f, "-.", label="Choi 2026 calibrated (eta_fib x factor)")
    except FileNotFoundError:
        pass
    ax.axvline(telescope.seeing_ref_arcsec, color="0.5", lw=0.8, ls=":", label=f"reference {telescope.seeing_ref_arcsec:g}\"")
    if telescope.seeing_arcsec is not None:
        ax.axvline(telescope.seeing_arcsec, color="C3", lw=0.8, label=f"current {telescope.seeing_arcsec:g}\" (x{telescope.seeing_factor(500.0):.2f})")
    ax.set_xlabel("seeing FWHM at 500 nm (arcsec)")
    ax.set_ylabel(f"fraction inside the {telescope.fiber_diameter_arcsec:.2f}\" fiber")
    ax.set_ylim(0, 1); ax.set_xlim(seeing.min(), seeing.max())
    ax.legend(fontsize=8)
    return ax


def plot_throughput_cascade(etc, ax=None, path="star"):
    """Cumulative throughput after each stage of the Choi et al. (2026) budget, then the
    spectrograph and detector, like their Figure 11 (``path="sky"`` drops the fiber injection)."""
    ax = _ax(ax, (10, 6.2))
    c = etc.throughput_cascade(path=path)
    w = c["wavelength"]
    n = len(c["stages"])
    colors = plt.cm.viridis(np.linspace(0.05, 0.9, n))
    for i, (s, cum) in enumerate(zip(c["stages"], c["cumulative"])):
        short = s["name"].split(":")[0].split(" (")[0]
        last = i == n - 1
        ax.plot(w, cum, color="k" if last else colors[i], lw=2.4 if last else 1.1,
                label=("end to end" if last else f"after {short}") + f" ({100 * np.nanmedian(cum):.0f}%)")
    ax.set_ylim(0, 1); ax.set_xlim(w.min(), w.max())
    ax.set_xlabel("Wavelength (nm)"); ax.set_ylabel("Cumulative throughput")
    ax.set_title(f"Throughput cascade, {path} path: {etc.telescope.fiber} fiber, {etc.instrument.name} (Choi et al. 2026 budget)", fontsize=9)
    ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2)
    return ax
