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
import pandas as pd
import matplotlib.pyplot as plt
from lfast_etc import LFASTETC, Source, Telescope, Instrument, Sky, list_instruments
from lfast_etc.plotting import plot_snr_map, plot_limiting_magnitude, plot_vs_n_telescopes, equivalent_diameter_axis"""

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

The LFAST prototype is an array of small telescopes: nominally 20 units of 0.726 m at Kitt Peak,
each feeding an 18 micron fiber, with the fibers merged into one spectrograph input by a fiber
combiner. `lfast_etc` answers the two
questions every observer asks of such a system. How long do I need to expose to reach a given
signal-to-noise ratio (SNR) on this target, and what SNR do I get in a given time?

Two spectrographs are modelled: the WIYN Bench (visible; the preset used below is the R ~ 1800
sample configuration) and CSHELL (near-infrared, 1.08 to 1.75 um, R = 18,000 to 43,000). The throughput
model calculations are based on the work of Choi et al. (2026, SPIE 14147).

A description of the Bench spectrograph, including modes and wavelength coverage, can be found here: https://www.wiyn.org/Instruments/wiynbench.html.

A description of the CSHELL spectrograph, including modes and wavelength coverage, can be found here: https://irtfweb.ifa.hawaii.edu/~cshell/info.html

After this notebook you will be able to build an ETC for any preset and array size, define a
target, read exposure times and SNRs off the result object, see which noise term limits you, and
plot what the extracted spectrum would look like.

## One idea to hold on to: everything is per resolution element

The spectrum is cut into bins one resolution element (resel) wide, $\\delta\\lambda = \\lambda/R$.
In each bin the source delivers $N_*$ electrons per second (summed over the whole array), the sky
delivers $N_{sky}$ electrons per second per pixel, and the resel is spread over $n_{pix}$ detector
pixels: the pixels per resel along the dispersion direction times the height of the fiber image on
the detector (2.2 x 2.2 = 4.8 for the Bench; 6 x 6 = 36 for CSHELL at R = 18,000, whose fiber
image is 6 px high). With dark current $D$, read noise $RN$ and a total
exposure $t$ split into $n_{exp}$ reads, the SNR per resel is (Howell 2006)

$$\\mathrm{SNR} = \\frac{N_* t}{\\sqrt{N_* t + n_{pix}(N_{sky} + D)\\,t + n_{pix}\\, n_{exp}\\, RN^2}}.$$

Survey papers often quote SNR per pixel. In the photon-limited regime that is lower by about
$\\sqrt{2.2}$, the pixels per resel along the dispersion direction, so SNR 30 per resel is roughly
SNR 20 per pixel on the Bench. Keep that
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
md("""## Available instrument modes

Each instrument mode is a preset: a set of spectrograph and detector parameters loaded by name. The
table lists the bundled presets. *One setting* is the wavelength range a single exposure covers:
the whole range for the Bench, but only one free spectral range (~2.5e-3 $\\lambda$) for CSHELL, so
a CSHELL observation covers ~3 nm at 1.08 um and ~4 nm at 1.55 um (see below). `n_pix` is the
extraction aperture per resel (dispersion pixels times the fiber-image height). The CSHELL slit
throughputs at R = 21,500 and 43,000 and the flat QEs are placeholders."""),
code("""rows = {}
for name in list_instruments():
    ins = Instrument.from_preset(name)
    lo, hi = ins.setting_window(0.5*(ins.wavelength_min_nm + ins.wavelength_max_nm))
    rows[name] = {"description": ins.description, "R": int(ins.resolving_power),
                  "range (nm)": f"{ins.wavelength_min_nm:.0f}-{ins.wavelength_max_nm:.0f}",
                  "one setting at centre (nm)": round(hi - lo, 1),
                  "efficiency": ins.efficiency, "QE": ins.qe, "pix/resel": ins.pixels_per_resel, "n_pix": ins.n_pix,
                  "read noise (e-)": ins.read_noise, "dark (e-/s/pix)": ins.dark_current,
                  "min exptime (s)": ins.min_exptime}
modes = pd.DataFrame(rows).T.rename_axis("preset")
modes"""),
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
code("""snr = 50
res = etc.exptime(src, snr=snr)
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
ax.axhline(snr, color='0.4', ls='--')
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
md("""## Switching instruments: one CSHELL setting

CSHELL is a single-order echelle: its free spectral range is only ~2.5e-3 $\\lambda$, so one setting
covers a narrow window (~4 nm at 1.55 um). `Instrument.setting_window(center_nm)` returns the
window of one setting, which is then passed to the ETC to restrict the wavelength grid."""),
code("""cs = Instrument.from_preset("cshell_r18000")
lo, hi = cs.setting_window(1555.0)
print(f"One CSHELL setting centred on 1555 nm covers {lo:.1f}-{hi:.1f} nm ({10*(hi - lo):.0f} A)")

cshell = LFASTETC(cs, wavelength_min_nm=lo, wavelength_max_nm=hi)
kstar = Source.blackbody(4500, mag=9, band="H")
r = cshell.snr(kstar, exptime=300)
print(f"H=9 K star, 300 s: median SNR/resel = {np.nanmedian(r.snr):.1f} | peak = {np.nanmax(r.peak_e_per_pixel):.0f} e-/pix | saturated: {bool(r.saturated.any())}")"""),
md("""## Overriding instrument parameters

Every preset field can be changed, which is how to test design choices: `Instrument.from_preset(name, **overrides)`
builds a modified preset and `instrument.replace(...)` copies an existing one with some fields changed.
For example, what would a lower-noise detector buy CSHELL? Here its read noise drops from 22 to 5 e-:"""),
code("""low_rn = cs.replace(read_noise=5.0)
r2 = LFASTETC(low_rn, wavelength_min_nm=lo, wavelength_max_nm=hi).snr(kstar, exptime=300)
for label, rr, ins in (("preset", r, cs), ("override", r2, low_rn)):
    print(f"{label:8s}: read noise {ins.read_noise:4.0f} e-  ->  median SNR/resel = {np.nanmedian(rr.snr):.1f}")"""),
md("""At H = 9 the star gives ~6600 e$^-$ per resel in 300 s, but the read noise, 22 e$^-$ on each of
the 36 pixels of the extraction aperture, adds a variance of ~17,000 e$^{-2}$, more than the star's
own photon noise. Lowering it to 5 e$^-$ therefore raises the SNR by ~50%; for fainter targets,
where read noise is an even larger share of the noise, the gain is larger still."""),
md("""
## Scaling the array

The number of telescopes is an ordinary parameter. With the 18 um fibers the telescopes of one
array feed the same spectrograph input through the combiner (`combine="fiber"`): signal and sky
per pixel grow with $N$ while the detector noise is paid once per input. Beyond 20 telescopes the
ETC assumes one combiner input per 20-telescope array (`telescopes_per_input=20`, the default), so
100 and 200 telescopes are 5 and 10 inputs, read out separately and co-added. The 50 um fiber does
not go through the combiner, so each telescope produces its own spectrum that is co-added
afterwards (`combine="independent"`), and the SNR grows as $\\sqrt{N}$. For a faint star the fiber
mode wins once read noise matters.
"""),
code("""faint_source = Source.blackbody(5800, mag=19, band="V")
snr_faint = 10
rows = {}
for n in (20, 100, 200):
    row = {}
    for fiber in ("18um", "50um"):
        e = LFASTETC("bench_r1800", telescope=Telescope(n_telescopes=n, fiber=fiber))
        row[f"{fiber} inputs"] = e.telescope.n_inputs
        row[f"{fiber} time (min)"] = round(float(np.nanmedian(e.exptime(faint_source, snr=snr_faint).exptime))/60, 1)
    rows[n] = row
pd.DataFrame(rows).T.rename_axis("N telescopes").astype({"18um inputs": int, "50um inputs": int})"""),
md("""With the 18 um fiber, 5x more telescopes cut the time by ~4x and 10x by ~6x. Co-adding the
inputs raises the SNR by $\\sqrt{n_{inputs}}$ at fixed time; that is worth a factor $n_{inputs}$ in
time only when the sky dominates the noise, and less here because read noise still contributes for
this faint star. The 50 um fiber is slower at every N: it collects more sky per telescope and pays
read noise once per telescope."""),
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
* Replace `"bench_r1800"` with `"bench_316_63p4_o9"` (the R ~ 31,000 H$\\alpha$ mode). What
  happens to the exposure time for SNR 50 and why?
"""),
])

# --------------------------------------------------------------------------- 02
nb("02_sensitivity_maps_and_array_scaling", [
md("""
# Sensitivity maps and array scaling

The source electron rate scales exactly as $10^{-0.4 m}$ with magnitude, so the ETC can evaluate
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
per unit time slows down. Following the SNR = 10 contour, a V ~ 18.9 star reaches it in 15 min; the
gain drops from 2.5 mag per decade at the bright end to 1.5 mag between 1000 and 10,000 s.

## CSHELL: SNR per resel versus H magnitude

Same map for CSHELL in its native R = 18,000 mode, for one echelle setting centred on 1.555 um (about 3.9 nm wide).
The grid starts at 8 s because that is the shortest Fowler-sampled read for which the quoted
read noise applies.
"""),
code("""cs = Instrument.from_preset("cshell_r18000")
lo, hi = cs.setting_window(1555.0)               # one echelle setting (free spectral range ~2.5e-3 lambda)
cshell = LFASTETC(cs, n_telescopes=20, wavelength_min_nm=lo, wavelength_max_nm=hi)
kstar = Source.blackbody(4500, mag=0, band="H")
mags_h = np.linspace(2, 16, 57)
exps_h = exps[exps >= cshell.instrument.min_exptime]
snr_h = cshell.snr_grid(kstar, mags_h, exps_h, wavelength=(lo, hi))
plot_snr_map(mags_h, exps_h, snr_h, band="H", levels=(3, 10, 30, 100, 300),
             title=f"CSHELL R=18k, 20 telescopes, one setting {lo:.1f}-{hi:.1f} nm");"""),
md("""
## Limiting magnitude versus wavelength

`plot_limiting_magnitude` solves for the magnitude that reaches the target SNR in each bin. For the
Bench the curve follows the throughput: deepest in the red, about a magnitude shallower at 420 nm.
The spikes are night-sky emission lines ([O I] 557.7 nm, Na D near 589 nm, [O I] 630.0/636.4 nm)
where the sky is locally much brighter.
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
and admits a 7.7 times larger patch of sky. The nominal 20-telescope prototype and the
100- and 200-telescope upgrade options are marked.
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
for n, txt, ha in ((20, "nominal 20", 'left'), (100, "100", 'left'), (200, "200", 'right')):
    ax.axvline(n, color='0.4', ls='--'); ax.text(n*1.05 if ha == 'left' else n/1.05, ax.get_ylim()[1]+0.3, txt, fontsize=8, ha=ha)
ax.legend(fontsize=8); ax.set_title("Bench R~1800, median over 500-650 nm", fontsize=9)
equivalent_diameter_axis(ax);"""),
md("""
The top axis converts the number of telescopes into the diameter of a single telescope with the
same collecting area. The same comparison for CSHELL, where the 22 e$^-$ read noise and the 0.5
e$^-$/s dark current, both over a 36-pixel extraction aperture, are large, shows the fiber combiner
mattering more:
"""),
code("""fig, ax = plt.subplots(figsize=(7, 4.5))
for fiber, combine in (("18um", "fiber"), ("18um", "independent"), ("50um", "independent")):
    lim = [LFASTETC("cshell_r18000", telescope=Telescope(n_telescopes=int(n), fiber=fiber, combine=combine),
                    wavelength_min_nm=lo, wavelength_max_nm=hi).limiting_magnitude(kstar, 900, 10)
           for n in ns]
    ax.plot(ns, lim, marker='o', ms=3, label=f"{fiber}, {combine}")
ax.set_xscale('log'); ax.invert_yaxis()
ax.set_xlabel("Number of telescopes"); ax.set_ylabel("Limiting H (SNR=10/resel, 900 s)")
for n in (20, 100, 200):
    ax.axvline(n, color='0.4', ls='--')
ax.legend(fontsize=8); ax.set_title(f"CSHELL R=18k, one setting {lo:.1f}-{hi:.1f} nm", fontsize=9)
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
in detail; generate a proper SkyCalc table for real planning. In a 15-minute exposure the sky
matches a V ~ 19 star in 500-650 nm, so this matters for the faint end and not at all for bright stars.
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
* Per-target overhead (slew, acquisition, read-out) of 10 min, a realistic allowance for a new pointing of the array.
"""),
code(PREAMBLE),
code("""bench = LFASTETC("bench_r1800", n_telescopes=20)
stars = {"G dwarf 5800 K": 5800, "K dwarf 4500 K": 4500, "M dwarf 3500 K": 3500}
mags = np.linspace(9, 20, 45)
window = (500, 560)
OVERHEAD = 600.0    # s per target (slew, acquisition, read-out)"""),
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
md("""## Survey depth: limiting magnitude per exposure length

The previous section asked how long a star of a given magnitude takes; this one turns the question
around and asks how faint we can go for a fixed time budget per star. The *limiting magnitude* is the
faintest V at which the median SNR per resel over 500-560 nm reaches the target in that exposure.
It is solved analytically from the noise equation, one bin at a time, so no grid of magnitudes is
needed.

The table gives the limiting magnitude for both SNR targets (30 for ~0.1 dex, 75 for ~0.05 dex) and
all three stars, in 5, 15 and 60 minutes of on-target time. Read it as: *with this exposure, every
star brighter than this V gets a metallicity at this precision.*"""),
code("""exps = (300, 900, 3600)
depth = pd.DataFrame({(label.split()[0], f"SNR {snr_t}"): [bench.limiting_magnitude(Source.blackbody(teff, 0, "V"), t, snr_t, wavelength=window) for t in exps]
                      for label, teff in stars.items() for snr_t in (30, 75)},
                     index=pd.Index([t // 60 for t in exps], name="exposure (min)"))
depth.round(1)"""),
md("""With 20 telescopes a G dwarf reaches the ~0.1 dex target at V ~ 16.9 in 15 min and V ~ 18.2 in
1 h; the ~0.05 dex target stops ~1.8 mag brighter. As before, the three stars differ by at most
~0.2 mag, and only the SED matters, so the 4500 K column applies equally to K giants.

In distance (ignoring extinction), V = 16.9 and 18.2 correspond to ~2.7 and ~4.7 kpc for a G dwarf
($M_V \\approx +4.8$), and V = 16.9 and 18.1 to ~19 and ~33 kpc for a K giant
($M_V \\approx +0.5$): dwarfs probe the thick disc, giants reach into the halo.

How the depth grows with time tells us what limits it. Going from 5 to 15 min (x3 in time) gains ~1.1 mag
and from 15 to 60 min (x4) ~1.25 mag. That lies between the source-limited scaling
($2.5\\log_{10} t$: 1.2 and 1.5 mag) and the sky-limited one ($1.25\\log_{10} t$: 0.6 and 0.75 mag),
closer to the former: at these magnitudes the star's own photon noise still dominates."""),
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
md("""## Survey throughput

LFAST observes one star at a time, since all telescopes point at the same target. The number of
stars per night therefore follows directly from what each star costs: its exposure time plus a fixed
overhead for slewing, acquisition and read-out,

$$N_{\\rm night} = \\frac{8\\ \\mathrm{h}}{t_{\\rm exp}(V) + t_{\\rm overhead}}.$$

The plot shows this for G dwarfs at the ~0.1 dex target (SNR = 30 per resel), with 20 telescopes and
a 10-minute overhead. Two regimes frame the curve:

* **Bright stars** ($t_{\\rm exp} \\ll t_{\\rm overhead}$): the overhead sets the yield, and the curve is flat.
* **Faint stars** ($t_{\\rm exp} \\gg t_{\\rm overhead}$): the exposure sets the yield, and it falls steeply with magnitude."""),
code("""t = bench.exptime_vs_mag(Source.blackbody(5800, 0, "V"), mags, snr=30, wavelength=window)
n_per_night = 8*3600 / (t + OVERHEAD)
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(mags, n_per_night)
ax.set_yscale('log'); ax.set_xlabel("V magnitude"); ax.set_ylabel("Stars per 8 h night (SNR=30)")
ax.set_title(f"G dwarfs, Bench R~1800, 20 telescopes, {OVERHEAD/60:.0f} min overhead", fontsize=9);"""),
md("""Bright stars are overhead-limited: at V ~ 13 the exposure is well under a minute, so the night yields
at most 48 stars (8 h / 10 min) whatever the array size, and a bigger array buys nothing. The
exposure equals the overhead at V ~ 16.5; fainter than that the exposure dominates, and the yield
falls to ~18 per night at V = 17 and ~8 at V = 18. A larger array helps only in this faint regime,
and a faster acquisition helps most for bright targets. With one target at a time LFAST cannot compete with multiplexed surveys (LAMOST observes
~4000 stars per pointing) on numbers; its niche is depth on individually selected targets (see the
Galactic-archaeology notebook)."""),
md("""## Array size and survey depth at SNR 30 in 15 minutes

The 20-telescope prototype is the nominal array; 100 and 200 telescopes are two possible upgrades.
In every case the fibers from all telescopes are combined into **one** spectrograph input
(`Telescope(..., telescopes_per_input=None)`). The terms in the noise budget therefore do not all
scale with the number of telescopes N in the same way:

| term | scales with N as | why |
|---|---|---|
| source signal | $N$ | every telescope collects light from the star |
| sky | $N$ | every telescope has its own fiber, and each fiber sees its own patch of sky: 20 fibers bring 20x the sky of one |
| read noise, dark current | 1 | there is one spectrograph, read out once |

Going from 20 to N telescopes multiplies the signal and the sky by N/20 and leaves the detector terms
unchanged. How much depth that buys depends on which term dominates the noise:

* source-dominated: $\\Delta m = 2.5\\log_{10}(N/20)$, i.e. 1.75 mag for 100 telescopes and 2.5 mag for 200;
* sky-dominated: half that, $1.25\\log_{10}(N/20)$, because the sky noise grows along with the signal;
* detector-dominated: the full $2.5\\log_{10}(N/20)$ again, because the read noise does not grow at all."""),
code("""rows = {}
for n in (20, 100, 200):
    e = LFASTETC("bench_r1800", telescope=Telescope(n_telescopes=n, telescopes_per_input=None))
    rows[n] = {"D_eq (m)": e.telescope.equivalent_diameter_m,
               "V_lim (SNR=30, 15 min)": e.limiting_magnitude(Source.blackbody(5800, 0, 'V'), 900, 30, wavelength=window)}
trade = pd.DataFrame(rows).T.rename_axis("N telescopes")
trade["gain vs N=20 (mag)"] = trade["V_lim (SNR=30, 15 min)"] - trade.loc[20, "V_lim (SNR=30, 15 min)"]
trade.round(2)"""),
md("""The gain (+1.4 mag for 100 telescopes, +2.0 for 200) sits between the two limits. At the
200-telescope limit (V ~ 18.9) the star and the sky each contribute ~1800 e$^-$ per resel in 15 min,
while the detector adds only ~70 e$^{-2}$ of variance, so for the Bench the detector noise hardly matters
either way: feeding each 20-telescope array into its own input instead would cost only ~0.1 mag
here. It matters much more for CSHELL, whose read noise is 6x higher.

Caveat: merging 100-200 telescopes into one spectrograph input needs a larger output fiber or a
second combiner stage. A larger fiber image covers more detector pixels per resel, which raises the
detector noise; the ETC keeps the pixels per resel fixed, so these numbers are an upper limit on the
gain."""),
md("""
## Things to try

* Move the metric window to the G band (425 to 435 nm) and compare the depth. This is the window
  that matters for carbon-enhanced metal-poor stars.
* Switch to the `bench_1200_28p7_o2` preset (R ~ 7600 around 400 nm). Does the higher resolution
  pay for itself in exposure time at fixed SNR per resel?
"""),
])

# --------------------------------------------------------------------------- 04
nb("04_tdamm_followup", [
md("""# Science case: TDAMM follow-up (transient classification)

Time-domain and multi-messenger (TDAMM) follow-up with a spectrograph is about
**classification** and **follow-up**: we will focus on the former in this notebook.
Classification typically uses a low-SNR spectrum, rebinned to R ~ 100-300, that can distinguish a SN Ia from a
SN II or a TDE from an AGN flare. The relevant question for LFAST
is: *what is the faint limit to classify a transient in a given exposure?* All times below are
on-target integration; slew, acquisition and other response-time overheads are not modelled.

**Assumptions used in this notebook**:

* Classification needs SNR ~ 5 per **rebinned** element at R_bin = 300 across 450-650 nm, the
  regime low-resolution classifiers are designed for (SEDM: Blagorodnova et al. 2018; 4MOST/TiDES:
  Swann et al. 2019). Co-adding $k = R/R_{bin}$ adjacent resels of equal SNR raises the SNR by
  $\\sqrt{k}$ in this noise model (every variance term scales with the number of pixels).
* 60 minutes is taken as the practical upper limit on the time spent on one classification spectrum.
* Instrument mode: one Bench setting, the R ~ 1800 preset (`bench_r1800`, 412-686 nm), with
  classification over 450-650 nm. This is only one of the Bench's grating and order settings
  (https://www.wiyn.org/Instruments/wiynbench.html); others cover redder or wider ranges at different
  resolution and ship as `bench_*` presets (`list_instruments()`). Several statements below about the 450-650 nm window are properties of this mode,
  and exploring other modes is left as an exercise (see the end of the notebook).
* Transient SEDs (blackbodies): a hot 25,000 K source (young SN II / TDE), a 10,000 K source
  (SN Ia near peak, where $B-V \\approx 0$; Phillips et al. 1999) and a red 3500 K source (kilonova at ~3 d, luminous red
  nova). Their magnitudes are **R** (Cousins, Vega; within ~0.1-0.2 mag of SDSS r for these colours),
  so the SED colour matters in the 450-650 nm window.
* GRB afterglow: a power law $F_\\nu \\propto \\nu^{-0.7}$ (mean optical slope 0.66; Kann et al. 2010), fading as $t^{-\\alpha}$ with $\\alpha \\approx 1$.
  Afterglow magnitudes are **AB**, as is usual in the GRB literature. For afterglows the metric is a
  redshift, not a classification (see that section).
* Kilonova: AT2017gfo-like, r = 17.3 at 0.5 d at 40 Mpc, fading 0.95 mag per day in r over the
  first 5 days (Drout et al. 2017, Table S2), cooling from ~10,000 K at 0.5 d to ~5000 K at
  1.5 d and approaching ~2500 K after ~5 d (Drout et al. 2017; Shappee et al. 2017). We only use the
  first 5 days, where this description holds.
* The sky background assumption is grey time (moon up); dark-time depths are ~0.3-0.5 mag deeper at the
  faint end (see the sensitivity-map notebook). Grey time is a fair average for target-of-opportunity
  follow-up, where the moon phase cannot be chosen."""),
code(PREAMBLE),
code("""# ---- knobs ----
N_TEL = 20
WINDOW = (450, 650)          # nm, classification window
R_BIN = 300                  # rebinned resolution for classification
SNR_CLASS = 5                # SNR per R_BIN element needed to classify
EXP_MAX = 3600               # s, practical upper limit for one classification spectrum
BAND = "R"                   # magnitudes in Cousins R (~ SDSS r)

bench = LFASTETC("bench_r1800", n_telescopes=N_TEL, wavelength_min_nm=WINDOW[0], wavelength_max_nm=WINDOW[1])
k = bench.instrument.resolving_power / R_BIN   # resels per rebinned element
snr_resel_needed = SNR_CLASS / np.sqrt(k)      # SNR per native resel that gives SNR_CLASS at R_BIN
print(f"{k:.1f} resels per R={R_BIN} bin -> need SNR {snr_resel_needed:.2f} per resel")

seds = {"hot 25,000 K": Source.blackbody(25000, 0, BAND),
        "SN Ia 10,000 K": Source.blackbody(10000, 0, BAND),
        "red 3500 K": Source.blackbody(3500, 0, BAND),
        "afterglow (AB)": Source.powerlaw(-0.7, 0, BAND, system="AB")}"""),
md("""## Exposure time to classify versus magnitude

The curves give the integration time needed to reach the classification SNR for each SED. The
reference magnitudes are the Argus Array detection depths (g, 5 sigma, dark time; science
requirements from the Argus specifications, argus.unc.edu/specifications; Law et al. 2022), since
Argus is (one of) the discovery engine(s) LFAST would follow up:

* **g = 20.0**, the single-epoch (1-minute) depth;
* **g = 21.5**, the 15-minute co-add depth.

These are detection limits, not classification limits, and they are in g (AB) while the curves are in
R (Vega). For SN-like colours near peak $g \\approx r$, and R (Vega) is ~0.2 mag brighter than r (AB);
both offsets are ignored here."""),
code("""mags = np.linspace(14, 22, 41)
fig, ax = plt.subplots(figsize=(7.5, 4.5))
for label, src in seds.items():
    t = bench.exptime_vs_mag(src, mags, snr=snr_resel_needed, wavelength=WINDOW)
    ax.plot(mags, t/60, label=label, ls='--' if 'afterglow' in label else '-')
ax.set_yscale('log'); ax.set_xlabel(f"{BAND} magnitude"); ax.set_ylabel("Exposure time (min)")
for tmin, ls, c in ((60, '-', '0.4'), (30, '--', '0.4'), (5, ':', '0.55')):
    ax.axhline(tmin, color=c, ls=ls, lw=0.8); ax.text(14.1, tmin*1.1, f"{tmin} min", fontsize=8, color=c)
ARGUS = {"Argus 1 min": 20.0, "Argus 15 min": 21.5}   # g, 5 sigma, dark time (science requirements)
for txt, mag in ARGUS.items():
    ax.axvline(mag, color='0.5', ls='-.', lw=0.8); ax.text(mag-0.08, 0.5, txt, fontsize=8, color='0.4', rotation=90, ha='right', va='top')
ax.legend(fontsize=8, loc='upper left', bbox_to_anchor=(0.0, 0.52))
ax.set_title(f"Classification: SNR={SNR_CLASS} per R={R_BIN} bin, {WINDOW[0]}-{WINDOW[1]} nm, {N_TEL} telescopes", fontsize=9)

t_ref = bench.exptime_vs_mag(seds["SN Ia 10,000 K"], list(ARGUS.values()), snr=snr_resel_needed, wavelength=WINDOW)
print("SN Ia: " + ";  ".join(f"at {txt} depth ({mag}): {t/60:.1f} min" for (txt, mag), t in zip(ARGUS.items(), t_ref)))"""),
md("""With 20 telescopes a SN Ia-like SED at the Argus 1 minute depth needs about 7 minutes, and one at the
15-minute co-add depth needs about an hour, right at the classification budget. The nominal array
could therefore classify what Argus finds in single epochs and the bright end of its 15-minute co-adds.
Because transient counts rise towards fainter magnitudes, a large share of the Argus detections
will sit close to these limits."""),
md("""## Classification depth versus exposure time

The table gives the faintest R magnitude that can be classified in 5, 15, 30 and 60 minutes (i.e. assuming SNR=5)."""),
code("""exps = (300, 900, 1800, 3600)
depth = pd.DataFrame({label: [bench.limiting_magnitude(src, t, snr_resel_needed, wavelength=WINDOW) for t in exps]
                      for label, src in seds.items()},
                     index=pd.Index([t // 60 for t in exps], name="exposure (min)"))
depth.round(1)"""),
md("""Note:
* **Exposure time.** Doubling the exposure from 30 to 60 min gains only ~0.4 mag. That is the
  background-limited scaling ($\\Delta m = 1.25\\log_{10} 2 \\approx 0.38$): at these magnitudes the sky
  and detector noise dominate, so SNR grows as $\\sqrt{t}$ rather than as $t$. The step from 5 to
  15 min is somewhat steeper because read noise still matters in the shortest exposures. Going
  beyond 60 min buys little: the next half-magnitude would cost another ~1.5 h."""),
md("""## How many days after merger can we still classify kilonovae?

Assumptions (AT2017gfo-like, at 40 Mpc):
* r = 17.3 at 0.5 d,
* fading 0.95 mag/day in r,
* cooling 10,000 K -> 5000 K -> 2500 K at 0.5, 1.5 and 5.5 d (linear ramps);

We stop at 5 d, the end of the range this light curve describes. 
We scale to other distances with the distance modulus and ask for a classification
spectrum (SNR=5 at R=300). 

40 Mpc is GW170817 itself; for reference the LIGO binary-neutron-star
range is 150-170 Mpc in O4 and 225-295 Mpc planned for O5 (the classification horizon versus
array size is in the LSST-era time-domain notebook). The blackbody is normalised in R."""),
code("""days = np.linspace(0.5, 5, 91)
def kn_mag(d, D_Mpc):
    return 17.3 + 0.95*(d-0.5) + 5*np.log10(D_Mpc/40)
def kn_teff(d):
    return np.interp(d, [0.5, 1.5, 5.5], [10000, 5000, 2500])

fig, ax = plt.subplots(figsize=(7.5, 4.5))
for D in (40, 100, 200):
    t_needed = np.array([bench.exptime_vs_mag(Source.blackbody(kn_teff(d), 0, BAND), [kn_mag(d, D)],
                                              snr=snr_resel_needed, wavelength=WINDOW)[0] for d in days]) / 60
    ax.plot(days, t_needed, label=f"{D} Mpc")
    last = {tmin: days[t_needed <= tmin].max() if (t_needed <= tmin).any() else np.nan for tmin in (30, 60)}
    print(f"{D:3d} Mpc: classifiable until {last[30]:.1f} d (30 min) / {last[60]:.1f} d (60 min)")
ax.axhline(60, color='0.4', lw=0.8); ax.text(0.6, 66, "60 min", fontsize=8, color='0.4')
ax.axhline(30, color='0.4', ls='--', lw=0.8); ax.text(0.6, 33, "30 min", fontsize=8, color='0.4')
ax.set_yscale('log'); ax.set_ylim(0.2, 1e3); ax.set_xlabel("Days after merger"); ax.set_ylabel("Exposure for classification (min)")
ax.legend(title="distance"); ax.set_title(f"AT2017gfo-like kilonova, Bench, {N_TEL} telescopes", fontsize=9);"""),
md("""A GW170817-like event at 40 Mpc stays classifiable for almost four days; at 100 Mpc the
window is about two nights, and at 200 Mpc (O5-like distances) only the first night, in this mode.
The kilonova fades by ~1 mag/day and also cools, moving its flux redward of the 450-650 nm window.
A redder Bench setting would follow it further; how much longer that keeps a kilonova classifiable
is left as an exercise (see the end of the notebook)."""),
md("""## GRB afterglows: time to a redshift

For an afterglow the gamma-ray or X-ray trigger typically already broadly indicates what the transient is; what the
spectrum adds is the **redshift**. It comes from absorption lines that the host-galaxy gas and any
intervening absorbers imprint on the afterglow continuum: Mg II 2796/2803 and Fe II at z <~ 2 and,
at higher redshift, the C IV and Si II lines, Ly-alpha / damped Ly-alpha absorption and the Ly-alpha
forest (Fynbo et al. 2009). Between them these put strong lines in the 450-650 nm window for
z ~ 0.6-4; other Bench settings shift or widen that redshift range. The strong lines have rest-frame equivalent widths of ~1-4 A (Mg II doublet ~3.8 A, C IV ~2.4 A, Fe II 2600 ~1.8 A, Si II 1527 ~1.0 A; de Ugarte Postigo et al. 2012), and the observed
equivalent width is larger by (1+z). At R = 1800 a resel is ~3 A at 550 nm, so a continuum SNR of
~10 per **native** resel (no rebinning) measures an equivalent width to ~0.3 A and detects Fe II 2600
at z = 1 (3.6 A observed) at ~10 sigma. That is the metric used here.

Afterglows fade as $F \\propto t^{-\\alpha}$ with $\\alpha \\approx 1$. The curves are for afterglows at
R = 17, 19 and 21 (AB) one hour after the burst. For afterglows detected within the first hour the
R magnitude at 1 h is distributed around R ~ 18.4 with a width of ~1.6 mag (Wang et al. 2013, Vega;
AB is ~0.2 mag fainter), so these curves span roughly +-1.5 sigma of that population. Optically dark
bursts are fainter still and are not reachable. The source is held at its
brightness at the start of the exposure, which is optimistic once the exposure becomes comparable
to the delay."""),
code("""SNR_Z = 10                    # SNR per native resel for an absorption-line redshift
delay_h = np.logspace(np.log10(0.25), np.log10(48), 40)
fig, ax = plt.subplots(figsize=(7.5, 4.5))
for m1h in (17, 19, 21):
    mag_t = m1h + 2.5*1.0*np.log10(delay_h/1.0)
    t_z = bench.exptime_vs_mag(seds["afterglow (AB)"], mag_t, snr=SNR_Z, wavelength=WINDOW)
    ax.plot(delay_h, t_z/60, label=f"R = {m1h} at 1 h")
ax.axhline(60, color='0.4', lw=0.8); ax.text(0.27, 66, "60 min", fontsize=8, color='0.4')
ax.axhline(30, color='0.4', ls='--', lw=0.8); ax.text(0.27, 33, "30 min", fontsize=8, color='0.4')
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_ylim(0.3, 1e3)
ax.set_xlabel("Hours after burst (start of exposure)"); ax.set_ylabel("Exposure for a redshift (min)")
ax.legend(fontsize=8)
ax.set_title(f"Afterglow, alpha=1, SNR={SNR_Z} per resel at R=1800, Bench, {N_TEL} telescopes", fontsize=9);"""),
md("""A bright afterglow (R = 17 at 1 h) gives a redshift in a couple of minutes at one hour and is
still within an hour's integration up to ~11 h after the burst. An R = 19 afterglow needs to be
observed within ~2 h of the burst, and an R = 21 afterglow is out of reach of the nominal array except in
the first ~20 minutes."""),
md("""## Array size trade for TDAMM depth (60 min, classification SNR)

The 20-telescope prototype is the nominal array; 100 and 200 telescopes are two possible upgrades.
In this model the telescopes are fiber-combined in groups of 20 (one combiner per array, as in
`Telescope(telescopes_per_input=20)`); each combiner output is a separate spectrograph input,
read out separately and co-added. Source and sky photons therefore grow with N, and so do dark
current and read noise (with the number of inputs, N/20). At the faint end every noise term grows in
proportion to the signal, so SNR $\\propto \\sqrt{N}$ and the depth improves by
$\\Delta m = 1.25\\log_{10}(N/20)$: about 0.9 mag for 100 telescopes and 1.25 mag for 200. That is half
the $2.5\\log_{10}(N/20)$ that the collecting area alone would suggest. A single combiner feeding all
N telescopes into one input would do better (see the array-scaling notebook). The table gives the
60-minute classification depth for each SED."""),
code("""rows = {}
for n in (20, 100, 200):
    e = LFASTETC("bench_r1800", n_telescopes=n, wavelength_min_nm=WINDOW[0], wavelength_max_nm=WINDOW[1])
    rows[n] = {"D_eq (m)": e.telescope.equivalent_diameter_m, "inputs": e.telescope.n_inputs,
               **{label: e.limiting_magnitude(src, EXP_MAX, snr_resel_needed, wavelength=WINDOW) for label, src in seds.items()}}
trade = pd.DataFrame(rows).T.rename_axis("N telescopes")
trade["inputs"] = trade["inputs"].astype(int)
trade["gain vs N=20 (mag)"] = trade["SN Ia 10,000 K"] - trade.loc[20, "SN Ia 10,000 K"]
trade.round(2)"""),
md("""Going from 20 to 200 telescopes pushes the 60-minute SN Ia classification depth from
R ~ 21.5, at the Argus 15-minute co-add depth, to R ~ 22.8, about a magnitude beyond it **assuming the Bench spectrograph**. For the kilonova, 1.25 mag is a factor ~1.8 in distance, so the
200-telescope array reaches at 200 Mpc roughly what the nominal array reaches at ~110 Mpc."""),
md("""## Things to try

**Other Bench modes.** Everything above uses one Bench setting. The Bench offers many grating and
order combinations (https://www.wiyn.org/Instruments/wiynbench.html), from coarse gratings covering
~500 nm in one setting to high-resolution orders covering ~25 nm, and they ship as `bench_*` presets
(`list_instruments()`). Redder settings suit cool sources (later phases in kilonovae, red transients)
and higher-redshift afterglow lines; bluer ones suit hot, young transients. Suggestions:

* how many more days a kilonova stays classifiable in a redder mode;
* how the colour ordering of the SEDs changes;
* which afterglow redshifts each setting covers.

The `bench_*` presets scale the spectrograph efficiency from the WIYN table, but the detector QE is
a flat 0.75 placeholder, so results at the red end are indicative.

**Other knobs.**

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
  (`full_well` from the preset, the 55 ke$^-$ 1%-linear well of Brown et al. 2026, with an 80% margin).
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

CSHELL's 22 e$^-$ read noise and 0.5 e$^-$/s dark current, over a 15-pixel extraction aperture at
R = 43k, are large. With all fibers combined the
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
