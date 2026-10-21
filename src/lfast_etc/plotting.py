"""Plot helpers for sensitivity maps."""
import numpy as np
import matplotlib.pyplot as plt


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
