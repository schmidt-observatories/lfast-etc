"""Build and execute the example notebooks in notebooks/.

Usage: .venv/bin/python scripts/build_notebooks.py [name ...]
"""
import sys
from pathlib import Path
import nbformat as nbf
from nbclient import NotebookClient

OUT = Path(__file__).resolve().parents[1] / "notebooks"
PREAMBLE = """import numpy as np
import matplotlib.pyplot as plt
from lfast_etc import LFASTETC, Source, Telescope, Instrument, list_instruments
from lfast_etc.plotting import plot_snr_map, plot_limiting_magnitude
try:
    plt.style.use('gks')
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
# LFAST ETC quickstart

`lfast_etc` is an exposure time calculator for the LFAST telescope array at Kitt Peak
(nominally 20 x 0.726 m unit telescopes feeding one spectrograph through an 18 um fiber combiner).
Two instruments are modelled: the **WIYN Bench** spectrograph (visible, R ~ 1800 sample
configuration) and **CSHELL** (1.08-1.75 um, R = 18,000 / 21,500 / 43,000).

Everything is computed **per resolution element**: for each wavelength bin of width
$\\delta\\lambda = \\lambda/R$ the star gives $N_*$ e$^-$/s inside the bin (summed over the array),
the sky gives $N_{sky}$ e$^-$/s/pixel, and with $n_{pix}$ pixels per resel, dark current $D$,
read noise $RN$ and $n_{exp}$ reads of total length $t$:

$$\\mathrm{SNR} = \\frac{N_* t}{\\sqrt{N_* t + n_{pix}(N_{sky} + D)\\,t + n_{pix}\\, n_{exp}\\, RN^2}}$$

The throughput model follows Choi et al. (2026, SPIE 14147). This notebook reproduces the
reference LFAST ETC notebook (2026-10-12) and shows the main API.
"""),
code(PREAMBLE),
md("## Build an ETC\n\nPick an instrument preset and the number of telescopes."),
code("""print(list_instruments())
etc = LFASTETC("bench_r1800", n_telescopes=20)
etc"""),
md("## Define a source\n\nSources are an SED shape plus one magnitude. Spectral lines can be injected as Gaussians."),
code("""src = Source.blackbody(5800, mag=15, band="V").add_line(656.3, depth=0.3, sigma_nm=0.05)
print(src.describe(), "| synthetic J =", round(src.magnitude("J"), 2))"""),
md("## Exposure time to reach a target SNR per resel"),
code("""res = etc.exptime(src, snr=50)
res.plot()
res.summary()"""),
md("The result object carries every per-bin quantity. Query one wavelength:"),
code("""res.at(550.0)"""),
md("## SNR for a fixed exposure time"),
code("""res = etc.snr(src, exptime=600)
ax = res.plot()
ax.axhline(50, color='0.4', ls='--')
res.summary()['snr']"""),
md("## What limits us? Electron rates per pixel"),
code("""res.plot_rates();"""),
md("""## Switching instruments and overriding parameters

CSHELL's free spectral range is tiny (~2.5e-3 lambda, i.e. ~40 A at 1.55 um), so a single setting
covers one narrow window; restrict the grid accordingly. Any instrument field can be overridden
through `Instrument.from_preset(name, **overrides)` or `instrument.replace(...)`."""),
code("""cshell = LFASTETC("cshell_r18000", wavelength_min_nm=1550, wavelength_max_nm=1554)
kstar = Source.blackbody(4500, mag=9, band="H")
r = cshell.snr(kstar, exptime=300)
print("CSHELL R=18k, H=9 K star, 300 s: median SNR/resel =", round(float(np.median(r.snr)), 1),
      "| peak e-/pix =", int(r.peak_e_per_pixel.max()), "| saturated:", r.saturated.any())

lowRN = Instrument.from_preset("cshell_r18000", read_noise=5.0)
r2 = LFASTETC(lowRN, wavelength_min_nm=1550, wavelength_max_nm=1554).snr(kstar, exptime=300)
print("same with RN = 5 e-:", round(float(np.median(r2.snr)), 1))"""),
md("""## Scaling the array

The number of telescopes is a plain parameter. With the 18 um fibers all telescopes feed one
spectrograph input (`combine="fiber"`): signal and sky per pixel scale with $N$ while the detector
noise does not. The 50 um fiber does not go through the combiner, so each telescope delivers its
own spectrum (`combine="independent"`, SNR $\\propto\\sqrt{N}$)."""),
code("""faint = Source.blackbody(5800, mag=19, band="V")
for n in (1, 5, 20, 50):
    for fiber in ("18um", "50um"):
        e = LFASTETC("bench_r1800", telescope=Telescope(n_telescopes=n, fiber=fiber))
        t = e.exptime(faint, snr=20)
        print(f"N={n:3d} {fiber}: median t(SNR=20) = {np.median(t.exptime)/60:7.1f} min  [{e.telescope.combine}]")"""),
md("## Your own spectrum\n\nAny `(wavelength_nm, flux)` arrays work; the magnitude sets the absolute scale."),
code("""wl = np.linspace(380, 1800, 3000)
flux = wl**-2 * (1 + 0.5*np.exp(-0.5*((wl-850)/20)**2))   # F_nu-flat continuum with a bump
user = Source.from_arrays(wl, flux, flux_unit="flam", mag=16, band="I")
r = LFASTETC("bench_r1800").snr(user, 900)
r.plot();"""),
])

# --------------------------------------------------------------------------- 02
nb("02_sensitivity_maps_and_array_scaling", [
md("""
# Brightness sensitivity maps and array scaling

Because the stellar electron rate scales exactly as $10^{-0.4 m}$, the ETC can evaluate SNR on a
(magnitude x exposure time) grid and solve the limiting magnitude **analytically** per bin.
This notebook builds the sensitivity maps for both instruments and shows how they move with the
number of telescopes, which is not yet fixed.
"""),
code(PREAMBLE),
md("## Bench: SNR per resel as a function of V magnitude and exposure time"),
code("""bench = LFASTETC("bench_r1800", n_telescopes=20)
gstar = Source.blackbody(5800, mag=0, band="V")   # the magnitude here is irrelevant for the grid
mags = np.linspace(8, 22, 57)
exps = np.logspace(0, 4, 49)                     # 1 s .. ~3 h
snr = bench.snr_grid(gstar, mags, exps, wavelength=(500, 650))
plot_snr_map(mags, exps, snr, band="V", title="Bench R~1800, 20 telescopes, median SNR/resel 500-650 nm");"""),
md("## CSHELL: SNR per resel vs H magnitude (R = 18,000 native mode, 1.55-1.56 um)"),
code("""cshell = LFASTETC("cshell_r18000", n_telescopes=20, wavelength_min_nm=1550, wavelength_max_nm=1560)
kstar = Source.blackbody(4500, mag=0, band="H")
mags_h = np.linspace(2, 16, 57)
snr_h = cshell.snr_grid(kstar, mags_h, exps, wavelength=(1550, 1560))
plot_snr_map(mags_h, exps, snr_h, band="H", levels=(3, 10, 30, 100, 300),
             title="CSHELL R=18k, 20 telescopes, 1.55-1.56 um");"""),
md("## Limiting magnitude versus wavelength"),
code("""fig, ax = plt.subplots(figsize=(9, 4))
plot_limiting_magnitude(bench, gstar, exptimes=[60, 600, 3600], snr=10, ax=ax)
ax.set_title("Bench, SNR=10 per resel, 20 telescopes", fontsize=9);"""),
code("""fig, ax = plt.subplots(figsize=(9, 4))
full = LFASTETC("cshell_r18000", n_telescopes=20)   # whole CSHELL range (many echelle settings)
plot_limiting_magnitude(full, kstar, exptimes=[60, 600, 3600], snr=10, ax=ax)
ax.set_title("CSHELL R=18k, SNR=10 per resel, 20 telescopes (sky lines dominate the structure)", fontsize=9);"""),
md("""## Scaling the number of telescopes

The 18 um fibers go through an $N$-to-1 combiner: $N_*$ and $N_{sky}$ both scale with $N$ while the
detector terms do not, so the gain is *better* than $\\sqrt{N}$ when read noise or dark current
matter. The 50 um fiber bypasses the combiner (one spectrum per telescope, SNR $\\propto\\sqrt N$)
but has ~2x the point-source coupling and a ~7.7x larger sky aperture."""),
code("""ns = np.array([1, 2, 4, 8, 12, 16, 20, 30, 40, 60, 80, 100, 150, 200])
rows = {}
for fiber in ("18um", "50um"):
    for combine in ("fiber", "independent"):
        lim = []
        for n in ns:
            e = LFASTETC("bench_r1800", telescope=Telescope(n_telescopes=n, fiber=fiber, combine=combine))
            lim.append(e.limiting_magnitude(gstar, exptime=900, snr=10, wavelength=(500, 650)))
        rows[(fiber, combine)] = np.array(lim)

fig, ax = plt.subplots(figsize=(7, 4.5))
for (fiber, combine), lim in rows.items():
    ax.plot(ns, lim, marker='o', ms=3, label=f"{fiber}, {combine}")
ax.set_xscale('log'); ax.set_xlabel("Number of telescopes"); ax.set_ylabel("Limiting V (SNR=10/resel, 900 s)")
ax.axvline(20, color='0.4', ls='--'); ax.text(21, ax.get_ylim()[0]+0.2, "nominal 20", fontsize=8)
ax.legend(fontsize=8); ax.set_title("Bench R~1800, median over 500-650 nm", fontsize=9);"""),
md("The same for CSHELL, where read noise (22 e$^-$) and dark current (0.5 e$^-$/s) are large, so the fiber combiner matters more:"),
code("""fig, ax = plt.subplots(figsize=(7, 4.5))
for fiber, combine in (("18um", "fiber"), ("18um", "independent"), ("50um", "independent")):
    lim = [LFASTETC("cshell_r18000", telescope=Telescope(n_telescopes=n, fiber=fiber, combine=combine),
                    wavelength_min_nm=1550, wavelength_max_nm=1560).limiting_magnitude(kstar, 900, 10)
           for n in ns]
    ax.plot(ns, lim, marker='o', ms=3, label=f"{fiber}, {combine}")
ax.set_xscale('log'); ax.set_xlabel("Number of telescopes"); ax.set_ylabel("Limiting H (SNR=10/resel, 900 s)")
ax.axvline(20, color='0.4', ls='--'); ax.legend(fontsize=8); ax.set_title("CSHELL R=18k, 1.55-1.56 um", fontsize=9);"""),
md("## Exposure time versus magnitude at fixed SNR\n\nA compact way to read off survey speed."),
code("""mags = np.linspace(8, 21, 53)
fig, ax = plt.subplots(figsize=(7, 4.5))
for snr_t in (10, 30, 100):
    t = bench.exptime_vs_mag(gstar, mags, snr=snr_t, wavelength=(500, 650))
    ax.plot(mags, t/60, label=f"SNR={snr_t}")
ax.set_yscale('log'); ax.set_xlabel("V magnitude"); ax.set_ylabel("Exposure time (min)")
ax.axhline(60, color='0.4', ls='--'); ax.legend(); ax.set_title("Bench R~1800, 20 telescopes, 5800 K star", fontsize=9);"""),
])

# --------------------------------------------------------------------------- 03
nb("03_stellar_metallicity", [
md("""
# Science case: stellar metallicity characterisation with the Bench

At R ~ 1800 the Bench sample configuration is LAMOST-like (LAMOST LRS is R ~ 1800 over
370-900 nm). Survey pipelines at that resolution deliver [Fe/H] to ~0.1-0.15 dex once the
SNR per pixel exceeds ~20-30, and ~0.05 dex above SNR ~ 50 (e.g. the LAMOST LASP and
SDSS/SEGUE SSPP performance). The 412-688 nm window carries the Ca I 422.7, Fe/Mg b (516-519) and
Na D (589) features that drive the metallicity sensitivity.

**Assumptions used here** (knobs, change freely):

* SNR per resel targets: 20 (rough [Fe/H], ~0.15 dex) and 50 (~0.05 dex).
* Metric wavelength window: 500-560 nm (Mg b region).
* Stars: 5800 K (G dwarf), 4500 K (K dwarf), 3500 K (M dwarf) blackbodies normalised in V.
"""),
code(PREAMBLE),
code("""bench = LFASTETC("bench_r1800", n_telescopes=20)
stars = {"G dwarf 5800 K": 5800, "K dwarf 4500 K": 4500, "M dwarf 3500 K": 3500}
mags = np.linspace(9, 20, 45)
window = (500, 560)"""),
md("## Exposure time to reach the SNR thresholds"),
code("""fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
for ax, snr_t in zip(axes, (20, 50)):
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
        m20 = bench.limiting_magnitude(Source.blackbody(teff, 0, "V"), exp, 20, wavelength=window)
        m50 = bench.limiting_magnitude(Source.blackbody(teff, 0, "V"), exp, 50, wavelength=window)
        row.append(f"{label.split()[0]}: V<{m20:.1f} (SNR20) / V<{m50:.1f} (SNR50)")
    print(f"{exp:5d} s  " + " | ".join(row))"""),
md("""## Where does the information come from? SNR across the band

Metallicity indicators are spread over the whole window; the SNR per resel is far from flat
because the throughput and the stellar SED both change strongly across 412-688 nm."""),
code("""fig, ax = plt.subplots(figsize=(9, 4))
for label, teff in stars.items():
    r = bench.snr(Source.blackbody(teff, 16, "V"), exptime=900)
    ax.plot(r.wavelength, r.snr, lw=0.8, label=label)
for w, name in ((422.7, "Ca I"), (430.8, "G band"), (486.1, "Hb"), (518.4, "Mg b"), (589.3, "Na D"), (656.3, "Ha")):
    ax.axvline(w, color='0.6', ls=':', lw=0.8); ax.text(w, ax.get_ylim()[1]*0.9 if ax.get_ylim()[1]>0 else 1, name, fontsize=7, rotation=90, va='top')
ax.set_xlabel("Wavelength (nm)"); ax.set_ylabel("SNR per resel"); ax.set_title("V = 16, 900 s, 20 telescopes", fontsize=9)
ax.legend(fontsize=8);"""),
md("""## Survey throughput

Number of stars per 8 h night as a function of magnitude for the SNR = 30 target
(per-target overhead of 2 min assumed; 1 target at a time since all fibers point at one star)."""),
code("""overhead = 120.0
t = bench.exptime_vs_mag(Source.blackbody(5800, 0, "V"), mags, snr=30, wavelength=window)
n_per_night = 8*3600 / (t + overhead)
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(mags, n_per_night)
ax.set_yscale('log'); ax.set_xlabel("V magnitude"); ax.set_ylabel("Stars per 8 h night (SNR=30)")
ax.set_title("G dwarfs, Bench R~1800, 20 telescopes, 2 min overhead", fontsize=9);"""),
md("""## Effect of the array size on survey depth at fixed SNR = 30 in 15 min"""),
code("""ns = [5, 10, 20, 40, 80]
for n in ns:
    e = LFASTETC("bench_r1800", n_telescopes=n)
    print(f"N={n:3d}: V_lim(SNR=30, 15 min) = {e.limiting_magnitude(Source.blackbody(5800,0,'V'), 900, 30, wavelength=window):.2f}")"""),
])

# --------------------------------------------------------------------------- 04
nb("04_tdamm_followup", [
md("""
# Science case: TDAMM follow-up (transient classification)

Time-domain and multi-messenger (TDAMM) follow-up with a spectrograph is mostly about
**classification**: a low-SNR spectrum, rebinned to R ~ 100-300, that distinguishes a SN Ia from a
SN II, a kilonova from an afterglow, or a TDE from an AGN flare. The relevant question for LFAST
is: *how faint a transient can we classify in a given exposure, and how quickly?*

**Assumptions used here** (knobs):

* Classification needs SNR ~ 5 per **rebinned** element at R_bin = 300 across 450-650 nm. In the
  photon-limited regime co-adding $k = R/R_{bin}$ resels raises SNR by $\\sqrt{k}$; we use that
  approximation (it is exact when all noise terms scale with the number of pixels, which they do
  here except for the sqrt in the Howell formula, so it is slightly optimistic in the read-noise limit).
* Transient SEDs: a hot blue blackbody (12,000 K, young SN / shock cooling), a 6000 K blackbody
  (SN Ia near peak) and a cool red 3500 K blackbody (late kilonova / SN II late). A power-law
  $F_\\nu \\propto \\nu^{-0.7}$ represents a GRB afterglow.
* Kilonova: AT2017gfo-like, peak $M_{r}$ ~ -16 at 40 Mpc gives r ~ 17 at peak, fading ~1 mag/day
  in the optical.
"""),
code(PREAMBLE),
code("""bench = LFASTETC("bench_r1800", n_telescopes=20, wavelength_min_nm=450, wavelength_max_nm=650)
R_BIN = 300
k = bench.instrument.resolving_power / R_BIN   # resels per rebinned element
def snr_binned(res):
    return res.snr * np.sqrt(k)

seds = {"blue 12,000 K": Source.blackbody(12000, 0, "V"),
        "SN Ia-like 6000 K": Source.blackbody(6000, 0, "V"),
        "red 3500 K": Source.blackbody(3500, 0, "V"),
        "afterglow Fnu~nu^-0.7": Source.powerlaw(-0.7, 0, "V", system="AB")}"""),
md("## Exposure time to classify versus magnitude"),
code("""mags = np.linspace(14, 22, 41)
snr_resel_needed = 5 / np.sqrt(k)       # SNR per native resel that gives SNR=5 at R=300
fig, ax = plt.subplots(figsize=(7.5, 4.5))
for label, src in seds.items():
    t = bench.exptime_vs_mag(src, mags, snr=snr_resel_needed, wavelength=(450, 650))
    ax.plot(mags, t/60, label=label)
ax.set_yscale('log'); ax.set_xlabel("V magnitude"); ax.set_ylabel("Exposure time (min)")
ax.axhline(30, color='0.4', ls='--'); ax.text(14.1, 33, "30 min", fontsize=8, color='0.4')
ax.axhline(5, color='0.55', ls=':'); ax.text(14.1, 5.5, "5 min", fontsize=8, color='0.55')
ax.legend(fontsize=8); ax.set_title(f"Classification: SNR=5 per R={R_BIN} bin, 450-650 nm, 20 telescopes", fontsize=9);"""),
md("## Classification depth versus exposure time"),
code("""for exp in (300, 900, 1800, 3600):
    lims = [f"{label}: V<{bench.limiting_magnitude(src, exp, snr_resel_needed, wavelength=(450, 650)):.1f}"
            for label, src in seds.items()]
    print(f"{exp:5d} s   " + " | ".join(lims))"""),
md("""## Kilonova: how many days after merger can we still classify?

AT2017gfo at 40 Mpc: r ~ 17.0 at 0.5 d, fading ~1 mag/day over the first few days while cooling
from ~10,000 K to ~3000 K (simple linear ramps used here). We scale to other distances with the
distance modulus and ask for a classification spectrum (SNR=5 at R=300) in a 30-minute exposure."""),
code("""days = np.linspace(0.5, 7, 27)
def kn_mag(d, D_Mpc):
    return 17.0 + 1.0*(d-0.5) + 5*np.log10(D_Mpc/40)
def kn_teff(d):
    return np.interp(d, [0.5, 2, 7], [10000, 5000, 3000])

fig, ax = plt.subplots(figsize=(7.5, 4.5))
for D in (40, 100, 200):
    t_needed = np.array([bench.exptime_vs_mag(Source.blackbody(kn_teff(d), 0, "V"), [kn_mag(d, D)],
                                              snr=snr_resel_needed, wavelength=(450, 650))[0] for d in days])
    ax.plot(days, t_needed/60, label=f"{D} Mpc")
ax.axhline(30, color='0.4', ls='--'); ax.text(0.6, 33, "30 min", fontsize=8, color='0.4')
ax.set_yscale('log'); ax.set_xlabel("Days after merger"); ax.set_ylabel("Exposure for classification (min)")
ax.legend(title="distance"); ax.set_title("AT2017gfo-like kilonova, Bench, 20 telescopes", fontsize=9);"""),
md("""## Rapid response: a fading afterglow

GRB afterglows fade as $F \\propto t^{-\\alpha}$ with $\\alpha \\sim 1$. For an afterglow at
R = 17 one hour after the burst, how does the time-to-classify grow with the delay before LFAST is
on target?"""),
code("""delay_h = np.logspace(0, np.log10(48), 25)
mag_t = 17 + 2.5*1.0*np.log10(delay_h/1.0)
t_cls = bench.exptime_vs_mag(seds["afterglow Fnu~nu^-0.7"], mag_t, snr=snr_resel_needed, wavelength=(450, 650))
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(delay_h, t_cls/60)
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel("Hours after burst"); ax.set_ylabel("Exposure for classification (min)")
ax.set_title("Afterglow R=17 at 1 h, alpha=1, Bench, 20 telescopes", fontsize=9);"""),
md("## Array size trade for TDAMM depth (30 min, classification SNR)"),
code("""for n in (5, 10, 20, 40, 80):
    e = LFASTETC("bench_r1800", n_telescopes=n, wavelength_min_nm=450, wavelength_max_nm=650)
    print(f"N={n:3d}: classification depth in 30 min, 6000 K: V<{e.limiting_magnitude(seds['SN Ia-like 6000 K'], 1800, snr_resel_needed, wavelength=(450,650)):.2f}")"""),
])

# --------------------------------------------------------------------------- 05
nb("05_exoplanet_atmospheres_cshell", [
md("""
# Science case: exoplanet atmospheres at high resolution with CSHELL

CSHELL at R = 43,000 (100 um slit) or R = 18,000 (bare fiber) covers 1.08-1.75 um, but its free
spectral range is only ~2.5e-3 lambda: **one setting is ~27 A wide at 1.08 um and ~40 A at 1.55 um**.
That rules out broad-band cross-correlation over hundreds of lines in one shot, and points to
**single strong features**. The best one is the metastable **He I 1083.3 nm triplet**, the
workhorse of escaping-atmosphere studies (Spake et al. 2018; Allart et al. 2018), where transit
depths of 0.5-5% over a ~1 A core are common for hot Jupiters and warm Neptunes.

**Assumptions used here** (knobs):

* He 10830 setting: 1081-1086 nm, R = 43,000 (3 resels across the ~0.7 A FWHM core) and R = 18,000.
* Detection metric: the excess absorption averaged over `n_core` resels in the planet line, measured
  in-transit minus out-of-transit, with equal in/out exposure. Precision per transit:
  $\\sigma = \\sqrt{2}/(\\mathrm{SNR}_{resel}\\sqrt{n_{core}\\,N_{exp,in}})$.
* Individual exposures of 300 s; the exposure is split into reads only if the brightest pixel would
  saturate (`full_well` from the preset, placeholder 80 ke$^-$).
* Transit duration 2 h (so 24 in-transit exposures) and the same out of transit.
* Slit throughput for R = 43,000 is a placeholder (0.5) until measured.
"""),
code(PREAMBLE),
code("""EXP = 300.0
N_IN = int(2*3600/EXP)      # in-transit exposures per transit
def he_setup(name, n_tel=20):
    return LFASTETC(name, n_telescopes=n_tel, wavelength_min_nm=1081, wavelength_max_nm=1086)

def precision_per_transit(etc, src, n_core, exp=EXP, n_in=N_IN):
    r = etc.snr(src, exp)
    if r.saturated.any():                         # split the exposure to stay below full well
        n_exp = int(np.ceil(r.peak_e_per_pixel.max() / (0.8*etc.instrument.full_well)))
        r = etc.snr(src, exp, n_exp=n_exp)
    snr = r.at(1083.3)["snr"]
    return np.sqrt(2) / (snr*np.sqrt(n_core*n_in)), snr, r.n_exp

targets = {  # (Teff, J mag, He depth): representative well-known He 10830 detections
    "HD 189733 b": (5000, 5.6, 0.010),
    "WASP-69 b":   (4700, 8.0, 0.035),
    "HAT-P-11 b":  (4800, 7.6, 0.012),
    "GJ 3470 b":   (3650, 8.8, 0.015),
    "WASP-107 b":  (4400, 9.4, 0.055),
}"""),
md("## Per-transit precision and detection significance"),
code("""for mode, n_core in (("cshell_r43000", 3), ("cshell_r18000", 2)):
    etc = he_setup(mode)
    print(f"--- {etc.instrument.description} ---")
    for name, (teff, jmag, depth) in targets.items():
        sig, snr, n_exp = precision_per_transit(etc, Source.blackbody(teff, jmag, "J"), n_core)
        print(f"{name:12s} J={jmag:4.1f}: SNR/resel/300s={snr:6.1f} (reads={n_exp})  "
              f"sigma={100*sig:.3f}% per transit  ->  {depth/sig:5.1f} sigma on a {100*depth:.1f}% line")"""),
md("## Detection significance versus host J magnitude"),
code("""jmags = np.linspace(4, 13, 37)
fig, ax = plt.subplots(figsize=(7.5, 4.5))
for mode, n_core, ls in (("cshell_r43000", 3, '-'), ("cshell_r18000", 2, '--')):
    etc = he_setup(mode)
    for depth, color in ((0.005, 'C0'), (0.02, 'C1')):
        sig = np.array([precision_per_transit(etc, Source.blackbody(4800, j, "J"), n_core)[0] for j in jmags])
        ax.plot(jmags, depth/sig, ls=ls, color=color, label=f"{etc.instrument.resolving_power:.0f}, depth {100*depth:g}%")
ax.axhline(5, color='0.4', ls='--'); ax.text(4.1, 5.5, "5 sigma", fontsize=8, color='0.4')
ax.set_yscale('log'); ax.set_xlabel("Host J magnitude"); ax.set_ylabel("He 10830 significance per transit")
ax.legend(fontsize=8, ncol=2); ax.set_title("CSHELL, 20 telescopes, 2 h transit, 300 s exposures, 4800 K host", fontsize=9);"""),
md("## Number of transits for a 5-sigma detection of a 1% line"),
code("""etc = he_setup("cshell_r43000")
sig = np.array([precision_per_transit(etc, Source.blackbody(4800, j, "J"), 3)[0] for j in jmags])
n_transits = np.ceil((5*sig/0.01)**2)
fig, ax = plt.subplots(figsize=(7, 4))
ax.step(jmags, n_transits, where='mid')
ax.set_yscale('log'); ax.set_xlabel("Host J magnitude"); ax.set_ylabel("Transits for 5 sigma (1% He line)")
ax.set_title("CSHELL R=43k, 20 telescopes", fontsize=9);"""),
md("""## Why read noise and the array size matter here

CSHELL's 22 e$^-$ read noise and 0.5 e$^-$/s dark current are large. With all fibers combined the
detector terms are paid once; with independent spectra they are paid $N$ times."""),
code("""ns = [5, 10, 20, 40, 80]
src = Source.blackbody(4800, 9.0, "J")
for combine in ("fiber", "independent"):
    row = []
    for n in ns:
        etc = LFASTETC("cshell_r43000", telescope=Telescope(n, "18um", combine=combine), wavelength_min_nm=1081, wavelength_max_nm=1086)
        row.append(f"N={n}: {100*precision_per_transit(etc, src, 3)[0]:.2f}%")
    print(f"{combine:12s} " + "  ".join(row))"""),
md("""## Cross-correlation within one echelle setting (order-of-magnitude)

For completeness: a CCF over $N_{lines}$ planetary lines of individual contrast $\\epsilon$
(relative to the star) in one 40 A window at 1.55 um gives per transit
$S/N_{CCF} \\approx \\epsilon\\, \\mathrm{SNR}_{resel} \\sqrt{N_{lines} N_{exp}}$.
With only ~10-30 usable water/CO lines in 40 A, the He line is the better bet, as the numbers show."""),
code("""etc = LFASTETC("cshell_r43000", n_telescopes=20, wavelength_min_nm=1550, wavelength_max_nm=1554)
for name, (teff, jmag, _) in list(targets.items())[:3]:
    hmag = Source.blackbody(teff, jmag, "J").magnitude("H")
    r = etc.snr(Source.blackbody(teff, hmag, "H"), EXP)
    snr = float(np.median(r.snr))
    for n_lines in (10, 30):
        print(f"{name:12s} H={hmag:4.1f}: SNR/resel={snr:5.0f}; CCF S/N per transit (eps=5e-4, {n_lines} lines) = {5e-4*snr*np.sqrt(n_lines*N_IN):.1f}")"""),
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
