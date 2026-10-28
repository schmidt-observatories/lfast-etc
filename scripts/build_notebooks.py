"""Build and execute the example notebooks 01-05 in notebooks/.

Usage: .venv/bin/python scripts/build_notebooks.py [name ...]

Each notebook is a list of ("md", text) / ("code", text) cells. Edit here, rerun, commit the
executed .ipynb. The science-survey notebooks 06-10 live in build_science_notebooks.py.
"""
import sys
from pathlib import Path
import nbformat as nbf
from nbclient import NotebookClient

OUT = Path(__file__).resolve().parents[1] / "notebooks"
PREAMBLE = """import numpy as np
import matplotlib.pyplot as plt
from lfast_etc import LFASTETC, Source, Telescope, Instrument, Sky, list_instruments
from lfast_etc.plotting import plot_snr_map, plot_limiting_magnitude, plot_vs_n_telescopes, equivalent_diameter_axis
try:
    plt.style.use('gks')      # the author's matplotlib style; skipped if you do not have it
except OSError:
    pass"""

NOTEBOOKS = {}


def nb(name, cells):
    NOTEBOOKS[name] = cells


def md(s):
    return ("md", s.strip())


def code(s):
    return ("code", s.strip())


# --------------------------------------------------------------------------- 01
nb("01_quickstart", [
md("""
# LFAST ETC quick start

LFAST is an array of small telescopes: nominally 20 units of 0.726 m at Kitt Peak, each feeding
a fiber, with the fibers merged into one spectrograph input. `lfast_etc` answers the two
questions every observer asks of such a system. How long do I need to expose to reach a given
signal-to-noise ratio (SNR) on this target, and what SNR do I get in a given time?

Two spectrographs are modelled: the WIYN Bench (visible; the preset used below is the R ~ 1800
sample configuration) and CSHELL (1.08 to 1.75 um, R = 18,000 to 43,000). The throughput
model follows Choi et al. (2026, SPIE 14147). The test suite reproduces the reference LFAST ETC
notebook of 2026-10-12 to one part in a million, so the numbers here are those of the reference
calculation.

After this notebook you will be able to build an ETC for any preset and array size, define a
target, read exposure times and SNRs off the result object, see which noise term limits you, and
plot what the extracted spectrum would look like.

## One idea to hold on to: everything is per resolution element

The spectrum is cut into bins one resolution element (resel) wide, $\\delta\\lambda = \\lambda/R$.
In each bin the star delivers $N_*$ electrons per second (summed over the whole array), the sky
delivers $N_{sky}$ electrons per second per pixel, and the resel is spread over $n_{pix}$ detector
pixels (2.2 for the Bench, 2.5 for CSHELL). With dark current $D$, read noise $RN$ and a total
exposure $t$ split into $n_{exp}$ reads, the SNR per resel is (Howell 2006)

$$\\mathrm{SNR} = \\frac{N_* t}{\\sqrt{N_* t + n_{pix}(N_{sky} + D)\\,t + n_{pix}\\, n_{exp}\\, RN^2}}.$$

Survey papers often quote SNR per pixel. In the photon-limited regime that is lower by about
$\\sqrt{n_{pix}}$, so SNR 30 per resel is roughly SNR 20 per pixel on the Bench. Keep that
conversion in mind when you compare with the literature.
"""),
code(PREAMBLE),
md("""
## Build an ETC

An `LFASTETC` is an instrument preset plus a telescope array plus a sky model. The presets are
TOML files inside the package; `list_instruments()` shows what ships. Printing the object tells
you what you built: resolving power, bandpass, pixels per resel, detector numbers, the array,
and the sky table.
"""),
code("""print(list_instruments())
etc = LFASTETC("bench_r1800", n_telescopes=20)
etc"""),
md("""
The grid has one point per resel. Here that is about 900 bins between 412 and 686 nm. The preset
stops at 686 nm because the telescope throughput curve has a zero at the 687 nm O$_2$ band.

## Where the photons go: throughput

$K_{opt}$ is the end-to-end telescope throughput for a point source, and it already contains the
fiber injection efficiency $\\eta_{fib}$ and the atmospheric extinction at the reference airmass
(that is how Choi et al. define it). The sky fills the fiber uniformly and pays no injection loss,
so it sees $K_{opt}/\\eta_{fib}$. The spectrograph efficiency and the detector QE are flat
placeholders in the presets. The black curve is the product of everything, and it is low: a few
percent of the photons that hit the primary mirrors end up as detected electrons.
"""),
code("""etc.plot_throughput();"""),
md("""
## Define a source

A `Source` is a spectral shape plus one magnitude in one band. The magnitude sets the absolute
scale; the shape comes from a blackbody, a power law, a flat spectrum or your own arrays. You can
add Gaussian lines. They are averaged over each resel so the equivalent width is preserved: a
0.05 nm wide H$\\alpha$ line with a 30% core, placed on a 0.36 nm resel, appears as a resel-wide
dip of about 10%. The magnitude you gave stays the continuum magnitude.
"""),
code("""src = Source.blackbody(5800, mag=15, band="V").add_line(656.3, depth=0.3, sigma_nm=0.05)
print(src.describe(), "| synthetic J =", round(src.magnitude("J"), 2))
src.plot(wavelength_nm=np.linspace(380, 1000, 2000));"""),
md("""
The shaded box is the normalising band. The synthetic J magnitude comes from integrating the same
SED over the J box band, which is how you get a magnitude in another band for a target you only
know in one.

## Exposure time to reach a target SNR

`exptime` returns an `ETCResult` with one value per bin. The plot shows the total exposure time
needed for SNR 50 per resel at each wavelength. Time rises steeply to the blue because both the
throughput and this star's SED fall there. The H$\\alpha$ line shows up as a spike: fewer photons
in the line, more time needed.
"""),
code("""res = etc.exptime(src, snr=50)
res.plot()
res.summary()"""),
md("""
`summary()` gives the 10th, 50th and 90th percentiles over the band. The median is the number to
quote if you need one exposure time for the whole spectrum. Pick one wavelength with `at`, or
export every column to an astropy table (`to_table().write("spectrum.csv")` for a file).
"""),
code("""print(res.at(550.0))
res.to_table()[:5]"""),
md("""
## SNR for a fixed exposure time

The inverse question. The dashed line marks SNR 50 so you can see where the 600 s exposure
reaches the target set above.
"""),
code("""res = etc.snr(src, exptime=600)
ax = res.plot()
ax.axhline(50, color='0.4', ls='--')
res.summary()['snr']"""),
md("""
## Which noise term limits you?

Two views of the same budget. The rates plot shows star, sky and dark electrons per pixel per
second; the noise-budget plot shows the fraction of the variance each term contributes for this
exposure. For a V = 15 star in 600 s the star dominates almost everywhere and the sky takes over
only at the bright airglow lines. Go three magnitudes fainter and the picture changes (try it).

The bundled sky table is grey time: the moon is up, 45 degrees from the target, and scattered
moonlight is about 60% of the visible sky. Faint-end results are therefore pessimistic by about
0.3 to 0.5 mag for dark time. Notebook 02 shows a crude dark-time scaling.
"""),
code("""fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
res.plot_rates(ax=axes[0])
res.plot_noise_budget(ax=axes[1])
axes[0].set_xlabel('');"""),
md("""
## What the spectrum would look like

`simulate` draws one realisation of the extracted spectrum with the per-resel SNR of the result,
and `plot_spectrum` overlays it on the noiseless signal. The H$\\alpha$ dip is clearly detected at
this SNR. This is the plot to show when someone asks whether a feature will be visible.
"""),
code("""res.plot_spectrum(seed=1);"""),
md("""
## The sky

Sky surface brightness per resel on the ETC grid, in AB magnitudes per square arcsecond (left
axis, brighter is up), with the sky electron rate per pixel in the 18 um fiber in grey. Around
550 nm the grey-time sky is about 20.3 mag per square arcsecond; dark Kitt Peak is about a
magnitude fainter. The spikes are airglow lines.
"""),
code("""etc.plot_sky();"""),
md("""
## Switching instruments and overriding parameters

CSHELL is an echelle with a small free spectral range, about 2.5e-3 of the wavelength, so one
setting covers about 40 A at 1.55 um. Restrict the grid to the setting you care about. Any
instrument field can be overridden through `Instrument.from_preset(name, **overrides)` or
`instrument.replace(...)`; the example lowers the read noise to see how much it matters for a
bright K star in 300 s.
"""),
code("""cshell = LFASTETC("cshell_r18000", wavelength_min_nm=1550, wavelength_max_nm=1554)
kstar = Source.blackbody(4500, mag=9, band="H")
r = cshell.snr(kstar, exptime=300)
print("CSHELL R=18k, H=9 K star, 300 s: median SNR/resel =", round(float(np.nanmedian(r.snr)), 1),
      "| peak e-/pix =", int(np.nanmax(r.peak_e_per_pixel)), "| saturated:", bool(r.saturated.any()))

lowRN = Instrument.from_preset("cshell_r18000", read_noise=5.0)
r2 = LFASTETC(lowRN, wavelength_min_nm=1550, wavelength_max_nm=1554).snr(kstar, exptime=300)
print("same with RN = 5 e-:", round(float(np.nanmedian(r2.snr)), 1))"""),
md("""
## Scaling the array

The number of telescopes is an ordinary parameter. With the 18 um fibers every telescope feeds the
same spectrograph input (`combine="fiber"`): signal and sky per pixel grow with $N$ while the
detector noise is paid once. The 50 um fiber does not go through the combiner, so each telescope
produces its own spectrum that is co-added afterwards (`combine="independent"`), and the SNR
grows as $\\sqrt{N}$. For a faint star the fiber mode wins once read noise matters.
"""),
code("""faint = Source.blackbody(5800, mag=19, band="V")
for n in (1, 5, 20, 50):
    for fiber in ("18um", "50um"):
        e = LFASTETC("bench_r1800", telescope=Telescope(n_telescopes=n, fiber=fiber))
        t = e.exptime(faint, snr=20)
        print(f"N={n:3d} {fiber}: median t(SNR=20) = {np.nanmedian(t.exptime)/60:7.1f} min  [{e.telescope.combine}]")"""),
md("""
## Your own spectrum

Any `(wavelength_nm, flux)` pair works. `flux_unit="flam"` means the array is proportional to
$F_\\lambda$ in whatever unit; the magnitude sets the scale. Use `"fnu"` for an $F_\\nu$ array and
`"W/m2/nm"` for an absolute spectrum, in which case the magnitude is derived rather than given.
"""),
code("""wl = np.linspace(380, 1800, 3000)
flux = wl**-2 * (1 + 0.5*np.exp(-0.5*((wl-850)/20)**2))   # F_nu-flat continuum with a bump
user = Source.from_arrays(wl, flux, flux_unit="flam", mag=16, band="I")
r = LFASTETC("bench_r1800").snr(user, 900)
r.plot();"""),
md("""
## Check yourself

The round trip below must close: the SNR computed for the exposure times that `exptime` returned
has to be the target SNR in every bin. The second check is the $\\sqrt{N}$ law for independent
spectra. If either fails, something in your installation is off.
"""),
code("""t = etc.exptime(src, snr=50)
back = etc.snr(src, exptime=t.exptime)
assert np.allclose(back.snr[np.isfinite(back.snr)], 50, rtol=1e-8)

one = LFASTETC("bench_r1800", telescope=Telescope(1, combine="independent")).snr(src, 100).snr
twenty = LFASTETC("bench_r1800", telescope=Telescope(20, combine="independent")).snr(src, 100).snr
assert np.allclose(twenty, np.sqrt(20) * one, rtol=1e-10, equal_nan=True)
print("round trip and sqrt(N) scaling OK")"""),
md("""
## Things to try

* Change the star to `Source.blackbody(3500, mag=18, band="V")` and rerun the noise budget. At
  which wavelengths does the sky take over?
* Set `n_exp=10` in the SNR call. How much SNR does the extra read noise cost at V = 15, and at
  V = 19?
* Replace `"bench_r1800"` with `"bench_316_63p4_o9"` (the R ~ 23,000 H$\\alpha$ mode). What
  happens to the exposure time for SNR 50 and why?
"""),
])

# --------------------------------------------------------------------------- 02
nb("02_sensitivity_maps_and_array_scaling", [
md("""
# Sensitivity maps and array scaling

The stellar electron rate scales exactly as $10^{-0.4 m}$ with magnitude, so the ETC can evaluate
the SNR on a (magnitude, exposure time) grid without redoing the spectrum each time, and it can
solve the limiting magnitude for a given SNR and exposure analytically, bin by bin. This notebook
builds those maps for both instruments and then asks how they move with the number of
telescopes, which is still an open design choice for LFAST.
"""),
code(PREAMBLE),
md("""
## Bench: SNR per resel versus V magnitude and exposure time

The map is the median SNR over 500 to 650 nm. Read it along a horizontal line to see how the SNR
grows with time at fixed magnitude, or along a contour to see how much fainter a longer exposure
reaches. The magnitude of the source object does not matter here; the grid supplies it.
"""),
code("""bench = LFASTETC("bench_r1800", n_telescopes=20)
gstar = Source.blackbody(5800, mag=0, band="V")
mags = np.linspace(8, 22, 57)
exps = np.logspace(0, 4, 49)                     # 1 s to about 3 h
snr = bench.snr_grid(gstar, mags, exps, wavelength=(500, 650))
plot_snr_map(mags, exps, snr, band="V", title="Bench R~1800, 20 telescopes, median SNR/resel 500-650 nm");"""),
md("""
Where the contours are straight lines of slope $-2.5$ mag per decade of time the measurement is
photon limited by the star. They bend at the faint end where sky and read noise enter and the gain
per unit time slows down.

## CSHELL: SNR per resel versus H magnitude

Same map for CSHELL in its native R = 18,000 mode, for one echelle setting at 1.55 to 1.56 um.
The grid starts at 8 s because that is the shortest Fowler-sampled read for which the quoted
read noise applies.
"""),
code("""cshell = LFASTETC("cshell_r18000", n_telescopes=20, wavelength_min_nm=1550, wavelength_max_nm=1560)
kstar = Source.blackbody(4500, mag=0, band="H")
mags_h = np.linspace(2, 16, 57)
exps_h = exps[exps >= cshell.instrument.min_exptime]
snr_h = cshell.snr_grid(kstar, mags_h, exps_h, wavelength=(1550, 1560))
plot_snr_map(mags_h, exps_h, snr_h, band="H", levels=(3, 10, 30, 100, 300),
             title="CSHELL R=18k, 20 telescopes, 1.55-1.56 um");"""),
md("""
## Limiting magnitude versus wavelength

`plot_limiting_magnitude` solves for the magnitude that reaches the target SNR in each bin. For the
Bench the curve follows the throughput: deepest in the red, about a magnitude shallower at 420 nm.
"""),
code("""fig, ax = plt.subplots(figsize=(9, 4))
plot_limiting_magnitude(bench, gstar, exptimes=[60, 600, 3600], snr=10, ax=ax)
ax.set_title("Bench, SNR=10 per resel, 20 telescopes", fontsize=9);"""),
code("""fig, ax = plt.subplots(figsize=(9, 4))
full = LFASTETC("cshell_r18000", n_telescopes=20)   # the whole CSHELL range, many echelle settings
plot_limiting_magnitude(full, kstar, exptimes=[60, 600, 3600], snr=10, ax=ax)
ax.set_title("CSHELL R=18k, SNR=10 per resel, 20 telescopes", fontsize=9);"""),
md("""
Each wavelength in the CSHELL plot is what one echelle setting, about 2.5e-3 of the wavelength
wide, would give; the curve is not a single exposure. The gaps between 1.36 and 1.42 um are bins
where $K_{opt}$ is zero in the water band, and the ETC returns NaN there rather than a number.
The remaining structure is the throughput curve (tellurics and fiber coupling) on top of the sky.

## Scaling the number of telescopes

The 18 um fibers go through a 20-to-1 combiner. Within one array both $N_*$ and $N_{sky}$ grow
with $N$ while the detector terms do not, so the gain is better than $\\sqrt{N}$ whenever read
noise or dark current matter. Beyond 20 telescopes the ETC assumes one combiner input per 20
telescopes, co-added incoherently (`telescopes_per_input=20`); the dotted curve is an idealised
single combiner for comparison. The 50 um fiber bypasses the combiner, so each telescope yields a
spectrum and the SNR grows as $\\sqrt{N}$, but it couples about twice as much point-source light
and admits a 7.7 times larger patch of sky.
"""),
code("""ns = np.array([1, 2, 4, 8, 12, 16, 20, 30, 40, 60, 80, 100, 150, 200])
configs = (("18um", "fiber", 20), ("18um", "fiber", None), ("18um", "independent", 20), ("50um", "independent", 20))
fig, ax = plt.subplots(figsize=(7, 4.5))
for fiber, combine, tpi in configs:
    lim = [LFASTETC("bench_r1800", telescope=Telescope(n_telescopes=int(n), fiber=fiber, combine=combine, telescopes_per_input=tpi))
           .limiting_magnitude(gstar, exptime=900, snr=10, wavelength=(500, 650)) for n in ns]
    ax.plot(ns, lim, marker='o', ms=3, ls=':' if tpi is None else '-',
            label=f"{fiber}, {combine}" + (" (ideal single combiner)" if tpi is None else ""))
ax.set_xscale('log'); ax.invert_yaxis()
ax.set_xlabel("Number of telescopes"); ax.set_ylabel("Limiting V (SNR=10/resel, 900 s)")
ax.axvline(20, color='0.4', ls='--'); ax.text(21, ax.get_ylim()[1]+0.3, "nominal 20", fontsize=8)
ax.legend(fontsize=8); ax.set_title("Bench R~1800, median over 500-650 nm", fontsize=9)
equivalent_diameter_axis(ax);"""),
md("""
The top axis converts the number of telescopes into the diameter of a single telescope with the
same collecting area. The same comparison for CSHELL, where the 22 e$^-$ read noise and 0.5
e$^-$/s dark current are large, shows the fiber combiner mattering more:
"""),
code("""fig, ax = plt.subplots(figsize=(7, 4.5))
for fiber, combine in (("18um", "fiber"), ("18um", "independent"), ("50um", "independent")):
    lim = [LFASTETC("cshell_r18000", telescope=Telescope(n_telescopes=int(n), fiber=fiber, combine=combine),
                    wavelength_min_nm=1550, wavelength_max_nm=1560).limiting_magnitude(kstar, 900, 10)
           for n in ns]
    ax.plot(ns, lim, marker='o', ms=3, label=f"{fiber}, {combine}")
ax.set_xscale('log'); ax.invert_yaxis()
ax.set_xlabel("Number of telescopes"); ax.set_ylabel("Limiting H (SNR=10/resel, 900 s)")
ax.axvline(20, color='0.4', ls='--'); ax.legend(fontsize=8); ax.set_title("CSHELL R=18k, 1.55-1.56 um", fontsize=9)
equivalent_diameter_axis(ax);"""),
md("""
For the common case of one configuration and several exposure times there is a one-call helper,
`plot_vs_n_telescopes`. You pass a function that builds the ETC for a given $N$, so you decide
the fiber, the combiner and the bandpass.
"""),
code("""plot_vs_n_telescopes(lambda n: LFASTETC("bench_r1800", n_telescopes=n), gstar,
                     ns=[1, 2, 5, 10, 20, 50, 100, 200], exptimes=[60, 600, 3600], snr=10, wavelength=(500, 650));"""),
md("""
## Dark time versus the bundled grey-time sky

The bundled SkyCalc table has the moon up (45 degrees away, 45 degrees altitude). Scattered
moonlight is about 57% of the 540 to 560 nm sky, giving V ~ 20.5 mag per square arcsecond, about
one magnitude brighter than dark Kitt Peak (V ~ 21.7 to 22.0 at zenith, Neugent & Massey 2010).
`Sky(scale=0.4)` is a crude dark-time estimate. It scales the airglow lines too, which is wrong
in detail; generate a proper SkyCalc table for real planning. The Bench becomes sky limited at
V > ~19.5, so this matters for the faint end and not at all for bright stars.
"""),
code("""for label, sky in (("grey (bundled)", Sky()), ("dark (x0.4)", Sky(scale=0.4))):
    e = LFASTETC("bench_r1800", n_telescopes=20, sky=sky)
    print(f"{label:15s}: V_lim(SNR=10, 900 s) = {e.limiting_magnitude(gstar, 900, 10, wavelength=(500, 650)):.2f}")"""),
md("""
## Exposure time versus magnitude at fixed SNR

The most compact way to read off survey speed. The dashed line is one hour.
"""),
code("""mags = np.linspace(8, 21, 53)
fig, ax = plt.subplots(figsize=(7, 4.5))
for snr_t in (10, 30, 100):
    t = bench.exptime_vs_mag(gstar, mags, snr=snr_t, wavelength=(500, 650))
    ax.plot(mags, t/60, label=f"SNR={snr_t}")
ax.set_yscale('log'); ax.set_xlabel("V magnitude"); ax.set_ylabel("Exposure time (min)")
ax.axhline(60, color='0.4', ls='--'); ax.legend(); ax.set_title("Bench R~1800, 20 telescopes, 5800 K star", fontsize=9);"""),
md("""
## Things to try

* Redo the Bench sensitivity map for a 3500 K star normalised in I. Which contour moves most?
* Add `n_exp=10` to `snr_grid`. Where on the map does the extra read noise show?
* Change `telescopes_per_input` to 40 in the scaling plot. Why does the curve change only beyond
  N = 20?
"""),
])

# --------------------------------------------------------------------------- 03
nb("03_stellar_metallicity", [
md("""
# Science case: stellar metallicities with the Bench

At R ~ 1800 the Bench sample configuration is LAMOST-like (LAMOST LRS is R ~ 1800 over 370 to
900 nm), so LAMOST tells us what to expect. Survey pipelines at that resolution deliver [Fe/H] to
about 0.1 dex for FGK stars once the SNR per pixel exceeds 20 to 30 (LAMOST LASP: Wu et al.
2014, Luo et al. 2015; SEGUE SSPP: about 0.11 dex for high-SNR F and G stars, Allende Prieto et
al. 2008), with a floor near 0.1 dex. About 0.05 dex is reached only by data-driven methods
trained on high-resolution labels (Ho et al. 2017; Xiang et al. 2019), and those need SNR > 50.
The ETC works per resel (2.2 pixels here), so SNR 20 per pixel is about 30 per resel. The 412 to
686 nm window contains Ca I 422.7, the G band, the Fe and Mg b lines at 516 to 519 nm and Na D at
589 nm, which is where the metallicity sensitivity comes from.

Assumptions, all of them changeable in the first code cell:

* SNR per resel targets of 30 (about 0.1 dex, LAMOST-like) and 75 (what the 0.05 dex data-driven
  methods need).
* Metric window 500 to 560 nm (the Mg b region).
* Stars as blackbodies normalised in V: 5800 K (G dwarf), 4500 K (K dwarf), 3500 K (M dwarf).
  The M dwarf is included only for the photon budget. M-dwarf metallicities at R ~ 1800 do not
  come from Mg b or Na D, and LAMOST M-dwarf [Fe/H] values are poor, so read that row as a
  throughput number and nothing more.
"""),
code(PREAMBLE),
code("""bench = LFASTETC("bench_r1800", n_telescopes=20)
stars = {"G dwarf 5800 K": 5800, "K dwarf 4500 K": 4500, "M dwarf 3500 K": 3500}
mags = np.linspace(9, 20, 45)
window = (500, 560)"""),
md("""
## Exposure time to reach the SNR thresholds

Left panel SNR 30, right panel SNR 75. The horizontal lines mark 15 minutes and one hour. The
three curves lie almost on top of each other. That is expected: the metric window sits inside the
V band, where the magnitude pins the flux, so spectral type only enters through the slope of the
SED across the 60 nm window. Spectral type would matter if the window were far from the band in
which the magnitude is known (try I-band normalisation).
"""),
code("""fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
for ax, snr_t in zip(axes, (30, 75)):
    for label, teff in stars.items():
        t = bench.exptime_vs_mag(Source.blackbody(teff, 0, "V"), mags, snr=snr_t, wavelength=window)
        ax.plot(mags, t/60, label=label)
    ax.set_yscale('log'); ax.set_xlabel("V magnitude"); ax.set_title(f"SNR = {snr_t} per resel, 500-560 nm", fontsize=9)
    for h, txt in ((15, "15 min"), (60, "1 h")):
        ax.axhline(h, color='0.4', ls='--'); ax.text(mags[0], h*1.1, txt, fontsize=8, color='0.4')
axes[0].set_ylabel("Exposure time (min)"); axes[0].legend();"""),
md("## Survey depth: limiting magnitude per exposure length"),
code("""for exp in (300, 900, 3600):
    row = []
    for label, teff in stars.items():
        m30 = bench.limiting_magnitude(Source.blackbody(teff, 0, "V"), exp, 30, wavelength=window)
        m75 = bench.limiting_magnitude(Source.blackbody(teff, 0, "V"), exp, 75, wavelength=window)
        row.append(f"{label.split()[0]}: V<{m30:.1f} (SNR30) / V<{m75:.1f} (SNR75)")
    print(f"{exp:5d} s  " + " | ".join(row))"""),
md("""
## Where the information comes from: SNR across the band

The metallicity indicators are spread over the window, and the SNR per resel is far from flat
across it because the throughput and the stellar SED both change strongly between 412 and 686 nm.
The dotted lines mark the main features. For the cool stars the blue end, where Ca I and the G
band sit, is where the SNR is lowest.
"""),
code("""fig, ax = plt.subplots(figsize=(9, 4))
for label, teff in stars.items():
    r = bench.snr(Source.blackbody(teff, 16, "V"), exptime=900)
    ax.plot(r.wavelength, r.snr, lw=0.8, label=label)
ymax = ax.get_ylim()[1]
for w, name in ((422.7, "Ca I"), (430.8, "G band"), (486.1, "H-beta"), (518.4, "Mg b"), (589.3, "Na D"), (656.3, "H-alpha")):
    ax.axvline(w, color='0.6', ls=':', lw=0.8); ax.text(w, ymax*0.9, name, fontsize=7, rotation=90, va='top')
ax.set_xlabel("Wavelength (nm)"); ax.set_ylabel("SNR per resel"); ax.set_title("V = 16, 900 s, 20 telescopes", fontsize=9)
ax.legend(fontsize=8);"""),
md("""
## Survey throughput

Stars per 8 hour night as a function of magnitude for the SNR 30 target, with 2 minutes of
overhead per target. All fibers point at one star, so there is no multiplexing: this is a
one-target-at-a-time instrument.
"""),
code("""overhead = 120.0
t = bench.exptime_vs_mag(Source.blackbody(5800, 0, "V"), mags, snr=30, wavelength=window)
n_per_night = 8*3600 / (t + overhead)
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(mags, n_per_night)
ax.set_yscale('log'); ax.set_xlabel("V magnitude"); ax.set_ylabel("Stars per 8 h night (SNR=30)")
ax.set_title("G dwarfs, Bench R~1800, 20 telescopes, 2 min overhead", fontsize=9);"""),
md("## Array size and survey depth at SNR 30 in 15 minutes"),
code("""for n in (5, 10, 20, 40, 80):
    e = LFASTETC("bench_r1800", n_telescopes=n)
    print(f"N={n:3d}: V_lim(SNR=30, 15 min) = {e.limiting_magnitude(Source.blackbody(5800,0,'V'), 900, 30, wavelength=window):.2f}")"""),
md("""
## Things to try

* Move the metric window to the G band (425 to 435 nm) and compare the depth. This is the window
  that matters for carbon-enhanced metal-poor stars.
* Switch to the `bench_1200_28p7_o2` preset (R ~ 3700 around 400 nm). Does the higher resolution
  pay for itself in exposure time at fixed SNR per resel?
"""),
])

# --------------------------------------------------------------------------- 04
nb("04_tdamm_followup", [
md("""
# Science case: transient classification for time-domain and multi-messenger follow-up

Spectroscopic follow-up of transients is mostly classification: a low-SNR spectrum, rebinned to
R ~ 100 to 300, that tells a SN Ia from a SN II, a kilonova from an afterglow, or a tidal
disruption event from an AGN flare. For LFAST the question is how faint a transient can be
classified in a given exposure, and how soon after the alert.

Assumptions, changeable in the first code cell:

* Classification needs SNR 5 per rebinned element at $R_{bin}$ = 300 across 450 to 650 nm, the
  regime low-resolution classifiers are built for (SEDM: Blagorodnova et al. 2018; 4MOST TiDES:
  Swann et al. 2019). Co-adding $k = R/R_{bin}$ adjacent resels of equal SNR raises the SNR by
  exactly $\\sqrt{k}$ in this noise model, because every variance term scales with the number of
  pixels.
* Transient SEDs: a 12,000 K blackbody (young SN or shock cooling), a 6000 K blackbody (SN Ia
  near peak), a 3500 K blackbody (late kilonova or late SN II), and a power law
  $F_\\nu \\propto \\nu^{-0.7}$ for a GRB afterglow. Magnitudes are R (Cousins; within 0.1 to 0.2
  mag of SDSS r for these colours), so the SED colour matters inside the 450 to 650 nm window.
* Kilonova: AT2017gfo-like, r = 17.3 at 0.5 d at 40 Mpc, fading 0.95 mag per day in r over the
  first 5 days (Drout et al. 2017, Table S2), cooling from ~10,000 K at 0.5 d to ~5000 K at
  1.5 d and approaching ~2500 K after ~5 d (Drout et al. 2017; Shappee et al. 2017).
* Sky: the bundled grey-time table (moon up). Dark time is 0.3 to 0.5 mag deeper at the faint
  end (notebook 02).
"""),
code(PREAMBLE),
code("""# ---- knobs ----
N_TEL = 20
WINDOW = (450, 650)          # nm, classification window
R_BIN = 300                  # rebinned resolution for classification
SNR_CLASS = 5                # SNR per R_BIN element needed to classify
BAND = "R"                   # magnitudes in Cousins R (~ SDSS r)

bench = LFASTETC("bench_r1800", n_telescopes=N_TEL, wavelength_min_nm=WINDOW[0], wavelength_max_nm=WINDOW[1])
k = bench.instrument.resolving_power / R_BIN   # resels per rebinned element
snr_resel_needed = SNR_CLASS / np.sqrt(k)      # SNR per native resel that gives SNR_CLASS at R_BIN
print(f"{k:.1f} resels per R={R_BIN} bin -> need SNR {snr_resel_needed:.2f} per resel")

seds = {"blue 12,000 K": Source.blackbody(12000, 0, BAND),
        "SN Ia-like 6000 K": Source.blackbody(6000, 0, BAND),
        "red 3500 K": Source.blackbody(3500, 0, BAND),
        "afterglow Fnu~nu^-0.7": Source.powerlaw(-0.7, 0, BAND, system="AB")}"""),
md("""
## Exposure time to classify versus magnitude

At the same R magnitude the blue SED reaches deepest and the red one shallowest, by about 0.8
mag. The magnitude is fixed at 641 nm, at the red edge of the classification window, so a blue
SED puts more light into 450 to 650 nm than a red one, and that outweighs the lower throughput
in the blue. The afterglow power law behaves like the blue blackbody.
"""),
code("""mags = np.linspace(14, 22, 41)
fig, ax = plt.subplots(figsize=(7.5, 4.5))
for label, src in seds.items():
    t = bench.exptime_vs_mag(src, mags, snr=snr_resel_needed, wavelength=WINDOW)
    ax.plot(mags, t/60, label=label)
ax.set_yscale('log'); ax.set_xlabel(f"{BAND} magnitude"); ax.set_ylabel("Exposure time (min)")
ax.axhline(30, color='0.4', ls='--'); ax.text(14.1, 33, "30 min", fontsize=8, color='0.4')
ax.axhline(5, color='0.55', ls=':'); ax.text(14.1, 5.5, "5 min", fontsize=8, color='0.55')
ax.legend(fontsize=8); ax.set_title(f"Classification: SNR={SNR_CLASS} per R={R_BIN} bin, {WINDOW[0]}-{WINDOW[1]} nm, {N_TEL} telescopes", fontsize=9);"""),
md("## Classification depth versus exposure time"),
code("""for exp in (300, 900, 1800, 3600):
    lims = [f"{label}: {BAND}<{bench.limiting_magnitude(src, exp, snr_resel_needed, wavelength=WINDOW):.1f}"
            for label, src in seds.items()]
    print(f"{exp:5d} s   " + " | ".join(lims))"""),
md("""
## Kilonova: how many days after merger can we still classify?

AT2017gfo at 40 Mpc: r = 17.3 at 0.5 d, fading 0.95 mag per day in r, cooling 10,000 K to 5000 K
to 2500 K at 0.5, 1.5 and 5.5 d (linear ramps). We scale to other distances with the distance
modulus and ask for a classification spectrum (SNR 5 at R = 300) in a 30 minute exposure. The
blackbody is normalised in R, so the cooling changes how much of the light falls inside the
classification window.
"""),
code("""days = np.linspace(0.5, 7, 27)
def kn_mag(d, D_Mpc):
    return 17.3 + 0.95*(d-0.5) + 5*np.log10(D_Mpc/40)
def kn_teff(d):
    return np.interp(d, [0.5, 1.5, 5.5], [10000, 5000, 2500])

fig, ax = plt.subplots(figsize=(7.5, 4.5))
for D in (40, 100, 200):
    t_needed = np.array([bench.exptime_vs_mag(Source.blackbody(kn_teff(d), 0, BAND), [kn_mag(d, D)],
                                              snr=snr_resel_needed, wavelength=WINDOW)[0] for d in days])
    ax.plot(days, t_needed/60, label=f"{D} Mpc")
ax.axhline(30, color='0.4', ls='--'); ax.text(0.6, 33, "30 min", fontsize=8, color='0.4')
ax.set_yscale('log'); ax.set_xlabel("Days after merger"); ax.set_ylabel("Exposure for classification (min)")
ax.legend(title="distance"); ax.set_title("AT2017gfo-like kilonova, Bench, 20 telescopes", fontsize=9);"""),
md("""
## Rapid response: a fading afterglow

GRB afterglows fade as $F \\propto t^{-\\alpha}$ with $\\alpha \\sim 1$. For an afterglow at
R = 17 (AB) one hour after the burst, the plot shows how the time needed to classify grows with
the delay before LFAST is on target. A fiber-fed array with no instrument changes can in
principle respond within minutes, which is the regime where the curve is flat.
"""),
code("""delay_h = np.logspace(0, np.log10(48), 25)
mag_t = 17 + 2.5*1.0*np.log10(delay_h/1.0)
t_cls = bench.exptime_vs_mag(seds["afterglow Fnu~nu^-0.7"], mag_t, snr=snr_resel_needed, wavelength=WINDOW)
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(delay_h, t_cls/60)
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel("Hours after burst"); ax.set_ylabel("Exposure for classification (min)")
ax.set_title("Afterglow R=17 at 1 h, alpha=1, Bench, 20 telescopes", fontsize=9);"""),
md("## Array size and classification depth in 30 minutes"),
code("""for n in (5, 10, 20, 40, 80):
    e = LFASTETC("bench_r1800", n_telescopes=n, wavelength_min_nm=WINDOW[0], wavelength_max_nm=WINDOW[1])
    print(f"N={n:3d}: classification depth in 30 min, 6000 K: {BAND}<{e.limiting_magnitude(seds['SN Ia-like 6000 K'], 1800, snr_resel_needed, wavelength=WINDOW):.2f}")"""),
md("""
## Things to try

* Raise `SNR_CLASS` to 10, the level at which a SN Ia subtype can be told apart. How much depth
  does that cost?
* Use `Sky(scale=0.4)` for dark time and redo the kilonova plot. Does the 200 Mpc curve cross the
  30 minute line?
"""),
])

# --------------------------------------------------------------------------- 05
nb("05_exoplanet_atmospheres_cshell", [
md("""
# Science case: exoplanet atmospheres at high resolution with CSHELL

CSHELL at R = 43,000 (100 um slit) or R = 18,000 (bare fiber) covers 1.08 to 1.75 um, but its
free spectral range is only about 2.5e-3 of the wavelength: one setting is about 27 A wide at
1.08 um and about 40 A at 1.55 um. That rules out cross-correlation over hundreds of molecular
lines in one shot and points to single strong features. The best one is the metastable He I
triplet at 1083.3 nm, the standard probe of escaping atmospheres (Spake et al. 2018; Allart et
al. 2018), where transit depths of 0.5 to 5% over a ~1 A core are common for hot Jupiters and
warm Neptunes.

Assumptions, changeable in the first code cell:

* He 10830 setting 1081 to 1086 nm (vacuum), at R = 43,000 and R = 18,000.
* Detection metric: the excess absorption averaged over `n_core` resels centred on the planet
  line (3 resels = 0.75 A at R = 43k, 2 resels = 1.2 A at R = 18k), measured in transit minus
  out of transit with equal in and out exposure. Precision per transit:
  $\\sigma = \\sqrt{2}/(\\mathrm{SNR}_{resel}\\sqrt{n_{core}\\,N_{exp,in}})$.
* The planet line is a Gaussian of 0.7 A FWHM (typical for He 10830). The literature depth is the
  peak depth, so the measured depth is the mean of that profile over the `n_core` window (0.79 at
  R = 43k, 0.59 at R = 18k). Instrumental broadening is not added on top of that.
* Individual exposures of 300 s, split into reads only if the brightest pixel would saturate
  (`full_well` from the preset, 40 ke$^-$ inferred from the IRTF page, with an 80% margin).
* Transit duration 2 h (24 in-transit exposures) and the same out of transit. The real durations
  are 1.8 h for HD 189733 b and 2.7 h for WASP-107 b.
* Not included, each a 20 to 40% effect: the planet signal sits in the stellar He I line core,
  which is 30 to 50% deep in K dwarfs, while the continuum SNR is used here; the OH airglow
  doublet at 10832 and 10834 A and the water tellurics near the line; the array-wide line spread
  function.
* The slit throughput at R = 43,000 is a placeholder (0.5) until measured.
"""),
code(PREAMBLE),
code("""# ---- knobs ----
N_TEL = 20
EXP = 300.0                  # s per exposure
T_TRANSIT = 2*3600           # s in transit (same out of transit)
LINE_FWHM_A = 0.7            # He 10830 planet line FWHM
HE_NM = 1083.3               # vacuum
MARGIN = 0.8                 # fraction of full well allowed
N_CORE = {"cshell_r43000": 3, "cshell_r18000": 2}   # resels averaged over the line core
from scipy.special import erf

N_IN = int(T_TRANSIT/EXP)    # in-transit exposures per transit

def he_setup(name, n_tel=N_TEL):
    return LFASTETC(name, n_telescopes=n_tel, wavelength_min_nm=1081, wavelength_max_nm=1086)

def core_fraction(etc, n_core):
    \"\"\"Mean of a Gaussian line (peak 1, FWHM LINE_FWHM_A) over n_core resels centred on it.\"\"\"
    sigma = LINE_FWHM_A/2.3548/10          # nm
    half = n_core*HE_NM/etc.instrument.resolving_power/2
    return sigma*np.sqrt(2*np.pi)/(2*half) * erf(half/(sigma*np.sqrt(2)))

def precision_per_transit(etc, src, exp=EXP, n_in=N_IN):
    \"\"\"(sigma on the measured depth per transit, continuum SNR/resel/exposure, reads per exposure).\"\"\"
    n_exp = etc.reads_needed(src, exp, margin=MARGIN)
    r = etc.snr(src, exp, n_exp=n_exp)
    snr = r.at(HE_NM)["snr"]
    n_core = N_CORE[etc.instrument.name]
    return np.sqrt(2) / (snr*np.sqrt(n_core*n_in)), snr, n_exp

def significance(etc, src, depth):
    \"\"\"Detection significance per transit of a line of peak depth `depth`.\"\"\"
    sig = precision_per_transit(etc, src)[0]
    return depth*core_fraction(etc, N_CORE[etc.instrument.name])/sig

targets = {  # (Teff [K], 2MASS J, peak He 10830 depth): well-known detections
    "HD 189733 b": (5000, 6.07, 0.0104),   # Salz et al. 2018
    "WASP-69 b":   (4715, 8.03, 0.0359),   # Nortmann et al. 2018
    "HAT-P-11 b":  (4780, 7.61, 0.0108),   # Allart et al. 2018
    "GJ 3470 b":   (3600, 8.79, 0.015),    # Palle et al. 2020 (Ninan et al. 2020: EW 0.012 A)
    "WASP-107 b":  (4430, 9.38, 0.0554),   # Allart et al. 2019 (Kirk et al. 2020: 7.3%)
}"""),
md("""
## Per-transit precision and detection significance

For each known detection the table gives the continuum SNR per resel in one 300 s exposure, how
many reads that exposure has to be split into to avoid saturation, the resulting precision on the
measured depth per transit, and the significance of the published line. Only HD 189733 at
R = 18k needs a second read; every other case fits in one.
"""),
code("""for mode in ("cshell_r43000", "cshell_r18000"):
    etc = he_setup(mode)
    frac = core_fraction(etc, N_CORE[mode])
    print(f"--- {etc.instrument.description} | core fraction {frac:.2f} over {N_CORE[mode]} resels ---")
    for name, (teff, jmag, depth) in targets.items():
        src = Source.blackbody(teff, jmag, "J")
        sig, snr, n_exp = precision_per_transit(etc, src)
        print(f"{name:12s} J={jmag:5.2f}: SNR/resel/300s={snr:6.1f} (reads={n_exp})  "
              f"sigma={100*sig:.3f}% per transit  ->  {significance(etc, src, depth):5.1f} sigma on a {100*depth:.1f}% line")"""),
md("""
## Detection significance versus host J magnitude

Solid lines are R = 43,000, dashed R = 18,000. The lower resolution comes out ahead at every
host magnitude. Its resels collect 2.4 times more light and it pays no slit loss, and that
outweighs the smaller fraction of the line depth it retains (0.59 against 0.79). Much of the gap
is the placeholder slit throughput of 0.5 at R = 43k; the last section shows what happens when
you set it to 1.
"""),
code("""jmags = np.linspace(4, 13, 37)
fig, ax = plt.subplots(figsize=(7.5, 4.5))
for mode, ls in (("cshell_r43000", '-'), ("cshell_r18000", '--')):
    etc = he_setup(mode)
    for depth, color in ((0.005, 'C0'), (0.02, 'C1')):
        s = np.array([significance(etc, Source.blackbody(4800, j, "J"), depth) for j in jmags])
        ax.plot(jmags, s, ls=ls, color=color, label=f"R={etc.instrument.resolving_power:.0f}, depth {100*depth:g}%")
ax.axhline(5, color='0.4', ls='--'); ax.text(4.1, 5.5, "5 sigma", fontsize=8, color='0.4')
ax.set_yscale('log'); ax.set_xlabel("Host J magnitude"); ax.set_ylabel("He 10830 significance per transit")
ax.legend(fontsize=8, ncol=2); ax.set_title(f"CSHELL, {N_TEL} telescopes, 2 h transit, 300 s exposures, 4800 K host (line smearing included)", fontsize=9);"""),
md("## Number of transits for a 5 sigma detection of a 1% line"),
code("""etc = he_setup("cshell_r43000")
n_transits = np.ceil((5/np.array([significance(etc, Source.blackbody(4800, j, "J"), 0.01) for j in jmags]))**2)
fig, ax = plt.subplots(figsize=(7, 4))
ax.step(jmags, n_transits, where='mid')
ax.set_yscale('log'); ax.set_xlabel("Host J magnitude"); ax.set_ylabel("Transits for 5 sigma (1% He line)")
ax.set_title("CSHELL R=43k, 20 telescopes", fontsize=9);"""),
md("""
## Why read noise and the array size matter here

CSHELL's 22 e$^-$ read noise and 0.5 e$^-$/s dark current are large. With all fibers combined the
detector terms are paid once per exposure; with independent spectra they are paid $N$ times.
"""),
code("""ns = [5, 10, 20, 40, 80]
src = Source.blackbody(4800, 9.0, "J")
for combine in ("fiber", "independent"):
    row = []
    for n in ns:
        etc = LFASTETC("cshell_r43000", telescope=Telescope(n, "18um", combine=combine), wavelength_min_nm=1081, wavelength_max_nm=1086)
        row.append(f"N={n}: {100*precision_per_transit(etc, src)[0]:.2f}%")
    print(f"{combine:12s} " + "  ".join(row))"""),
md("""
## Cross-correlation within one echelle setting, as an order of magnitude

A cross-correlation over $N_{lines}$ planetary lines of individual contrast $\\epsilon$ (relative
to the star) in one 40 A window at 1.55 um gives, per transit,
$S/N_{CCF} \\approx \\epsilon\\, \\mathrm{SNR}_{resel} \\sqrt{N_{lines} N_{exp}}$. The values
$\\epsilon$ = 5e-4 and $N_{lines}$ = 10 to 30 usable water or CO lines per 40 A window are
guesses, not measurements. With so few lines per setting the He line is the better bet.
"""),
code("""etc = LFASTETC("cshell_r43000", n_telescopes=20, wavelength_min_nm=1550, wavelength_max_nm=1554)
for name, (teff, jmag, _) in list(targets.items())[:3]:
    hmag = Source.blackbody(teff, jmag, "J").magnitude("H")
    r = etc.snr(Source.blackbody(teff, hmag, "H"), EXP)
    snr = float(np.nanmedian(r.snr))
    for n_lines in (10, 30):
        print(f"{name:12s} H={hmag:4.1f}: SNR/resel={snr:5.0f}; CCF S/N per transit (eps=5e-4, {n_lines} lines) = {5e-4*snr*np.sqrt(n_lines*N_IN):.1f}")"""),
md("""
## Things to try

* Set `slit_throughput=1.0` on the R = 43k preset (`Instrument.from_preset("cshell_r43000",
  slit_throughput=1.0)`) to see the best case once the slit loss is measured.
* Change `EXP` to 60 s. The bright hosts no longer need to be split into reads, but each exposure
  now pays read noise against less signal. Which host magnitude is the break-even point?
"""),
])


def build(name, cells):
    nbk = nbf.v4.new_notebook()
    nbk.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
    for kind, src in cells:
        nbk.cells.append(nbf.v4.new_markdown_cell(src) if kind == "md" else nbf.v4.new_code_cell(src))
    client = NotebookClient(nbk, timeout=600, kernel_name="python3", resources={"metadata": {"path": str(OUT)}})
    client.execute()
    path = OUT / f"{name}.ipynb"
    nbf.write(nbk, path)
    print("wrote", path)


if __name__ == "__main__":
    names = sys.argv[1:] or list(NOTEBOOKS)
    for n in names:
        build(n, NOTEBOOKS[n])
