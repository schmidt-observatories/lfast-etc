"""Build the array-scaling science-survey notebooks (06-10).

Usage: .venv/bin/python scripts/build_science_notebooks.py [name ...]
"""
import sys
from build_notebooks import build, md, code, PREAMBLE

SCALES = '''# ---- knobs: array sizes compared throughout (replicated in arrays of 20) ----
SCALES = {20: "nominal 20x (3.2 m eq.)", 200: "10 arrays (10.3 m eq.)", 1720: "86 arrays (30.1 m eq.)"}
NS = np.array([20, 40, 60, 100, 200, 400, 600, 1000, 1720, 2500])
FIG = "../memo/figures/"      # memo figures are saved here as PDF
plt.rcParams.update({"font.size": 11, "axes.titlesize": 10, "legend.fontsize": 9})
def tel(n, telescopes_per_input=20):
    return Telescope(n_telescopes=int(n), fiber="18um", combine="fiber", telescopes_per_input=telescopes_per_input)
def deq(n):
    return 0.726*np.sqrt(n)
from matplotlib.ticker import FixedLocator, FixedFormatter, NullFormatter
def deq_axis(ax):
    """Top axis in equivalent diameter with plain tick labels."""
    sec = ax.secondary_xaxis('top', functions=(deq, lambda d: (d/0.726)**2))
    ticks = [3, 5, 10, 20, 30]
    sec.xaxis.set_major_locator(FixedLocator(ticks)); sec.xaxis.set_major_formatter(FixedFormatter([str(t) for t in ticks]))
    sec.xaxis.set_minor_formatter(NullFormatter()); sec.set_xlabel("Equivalent diameter (m)")
    return sec'''

NB = {}

# --------------------------------------------------------------------------- 06
NB["06_array_scaling_overview"] = [
md("""
# Array scaling overview: 20x, 10 m equivalent, 30 m equivalent

LFAST is replicated in arrays of 20 unit telescopes. The nominal Kitt Peak array (20 x 0.726 m,
8.3 m$^2$) has the collecting area of a **3.2 m** telescope; 10 arrays (200 telescopes) match a
**10.3 m** aperture and 86 arrays (1720 telescopes) a **30.1 m** aperture. This notebook sets up
the scaling model used by the science-case notebooks 07-10 and produces the headline sensitivity
curves.

**How arrays combine in the ETC.** Each 20-to-1 fiber combiner feeds one spectrograph input.
Within an array, star *and* sky electrons scale with the number of telescopes while detector
noise is paid once, so the SNR gain is better than $\\sqrt{N}$ when read noise or dark current
matter. Beyond 20, the ETC assumes one input per array, co-added incoherently, so the SNR grows as
$\\sqrt{n_{arrays}}$ on top of the single-array value (`Telescope(..., telescopes_per_input=20)`).
An idealised single combiner (`telescopes_per_input=None`) is shown dotted where it matters.

**Caveat.** This is a photon-budget scaling only. It assumes the spectrographs scale with the
arrays (one input per array on the same or replicated spectrographs), that the per-array
throughput stays at the Choi et al. (2026) values, and that sky subtraction, flat-fielding and
fiber-combiner stability do not limit faint-end performance. Multiplexing is not modelled: all
telescopes point at **one target** at a time.
"""),
code(PREAMBLE), code(SCALES),
md("## Equivalent apertures"),
code("""print(f"{'N tel':>6} {'arrays':>7} {'area m2':>8} {'D_eq m':>7} {'inputs':>7}")
for n in NS:
    t = tel(int(n))
    print(f"{n:6d} {-(-n//20):7d} {t.total_area_m2:8.1f} {t.equivalent_diameter_m:7.1f} {t.n_inputs:7d}")"""),
md("## Limiting magnitude versus array size, both instruments"),
code("""gstar = Source.blackbody(5800, 0, "V")
kstar = Source.blackbody(4500, 0, "H")
cases = {
    "Bench R~1800, V, 500-650 nm": ("bench_r1800", gstar, dict(), (500, 650)),
    "CSHELL R=18k, H, 1.55-1.56 um": ("cshell_r18000", kstar, dict(wavelength_min_nm=1550, wavelength_max_nm=1560), None),
    "CSHELL R=43k, J, He 10830": ("cshell_r43000", Source.blackbody(4800, 0, "J"), dict(wavelength_min_nm=1081, wavelength_max_nm=1086), 1083.3),
}
fig, axes = plt.subplots(1, 3, figsize=(12, 3.9))
summary = {}
for ax, (label, (preset, src, kw, win)) in zip(axes, cases.items()):
    for exp, color in ((60, 'C2'), (900, 'C0'), (3600, 'C1')):
        lim = [LFASTETC(preset, telescope=tel(int(n)), **kw).limiting_magnitude(src, exp, 10, wavelength=win) for n in NS]
        ideal = [LFASTETC(preset, telescope=tel(int(n), telescopes_per_input=None), **kw).limiting_magnitude(src, exp, 10, wavelength=win) for n in NS]
        ax.plot(NS, lim, marker='o', ms=3, color=color, label=f"{exp} s")
        ax.plot(NS, ideal, ls=':', color=color, lw=1)
        summary[(label, exp)] = dict(zip(NS, lim))
    for n in (20, 200, 1720):
        ax.axvline(n, color='0.6', ls='--', lw=0.8)
    ax.set_xscale('log'); ax.set_xlabel("Number of telescopes"); ax.set_ylabel(f"Limiting {src.band} (SNR=10/resel)")
    ax.set_title(label + "\\n(dotted: ideal single combiner)", fontsize=9); ax.legend(fontsize=8)
    deq_axis(ax)
plt.tight_layout(); fig.savefig(FIG + "limiting_mag_vs_N.pdf")"""),
code("""print("Limiting magnitude, SNR=10 per resel:")
for (label, exp), d in summary.items():
    print(f"{label:32s} {exp:5d} s: " + "  ".join(f"N={n}: {d[n]:5.2f}" for n in (20, 200, 1720)))"""),
md("""## How the gain is split: source-, sky- and read-noise-limited regimes

At fixed SNR the required exposure time scales as $1/N$ in every regime once the array count
grows (incoherent co-add of inputs), and *faster* than $1/N$ inside a single array when detector
noise matters. The magnitude gain at fixed time is 2.5 log$_{10} N$ when source-limited and
1.25 log$_{10} N$ when background-limited, so the faint-end gain from 20 to 1720 telescopes is
between 2.4 and 4.8 mag."""),
code("""fig, ax = plt.subplots(figsize=(7, 4.3))
for mag, label in ((12, "V=12 (source-limited)"), (18, "V=18"), (21, "V=21 (sky-limited)")):
    t = [LFASTETC("bench_r1800", telescope=tel(int(n))).exptime_vs_mag(gstar, [mag], snr=10, wavelength=(500, 650))[0] for n in NS]
    ax.plot(NS, np.array(t)/t[0], marker='o', ms=3, label=label)
ax.plot(NS, 20/NS, color='0.4', ls='--', label="1/N")
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel("Number of telescopes"); ax.set_ylabel("Exposure time relative to 20x")
ax.legend(fontsize=8); ax.set_title("Bench, SNR=10 per resel, 5800 K", fontsize=9)
fig.savefig(FIG + "exptime_scaling.pdf")"""),
md("""## Bright end: saturation

With one combiner input per 20-telescope array, the light per detector pixel is set by a single
array and does **not** grow with the number of arrays, so the bright limit is the same at 20, 200
and 1720 telescopes (an idealised single combiner would instead saturate N/20 times sooner). The
table gives the brightest star for which a 1 s (Bench) or 8 s (CSHELL, Fowler minimum) read stays
below 80% of full well (placeholder full wells)."""),
code("""for preset, band, teff, kw, tmin in (("bench_r1800", "V", 5800, {}, 1.0), ("cshell_r18000", "H", 4500, dict(wavelength_min_nm=1550, wavelength_max_nm=1560), 8.0)):
    row = []
    for n in (20, 200, 1720):
        e = LFASTETC(preset, telescope=tel(n), **kw)
        mags = np.arange(-2, 14, 0.1)
        peak = np.array([np.nanmax(e.snr(Source.blackbody(teff, m, band), tmin).peak_e_per_pixel) for m in mags])
        ok = mags[peak < 0.8*e.instrument.full_well]
        row.append(f"N={n}: {band}>{ok.min():.1f}")
    print(f"{preset:14s} ({tmin:.0f} s read): " + "   ".join(row))"""),
]

# --------------------------------------------------------------------------- 07
NB["07_time_domain_lsst_era"] = [
md("""
# Time-domain science in the LSST era: classification depth, kilonova horizon, SN Ia cosmology, TDEs and AGN

A dedicated spectroscopic array is most valuable where the bottleneck is **telescope time on
single targets**: transient classification, rapid multi-messenger response and high-cadence
monitoring. This notebook asks what magnitude and distance limits the nominal 20x array reaches
and how they move at 10 m and 30 m equivalent.

**Assumptions used here** (knobs):

* Classification: SNR 5 per R = 300 bin over 450-650 nm (SEDM-like, Blagorodnova et al. 2018),
  i.e. SNR 2.04 per native resel; magnitudes in R (Cousins ~ SDSS r).
* Kilonova: AT2017gfo-like, r = 17.3 at 1 d at 40 Mpc, fading 0.6 mag/day, cooling 10,000 K ->
  5000 K -> 3000 K at 0.5 / 1.5 / 3.5 d (Drout et al. 2017; Villar et al. 2017).
* SN Ia: $M_B$ = -19.3 at peak, a 6000 K blackbody, no K-correction; luminosity distance from a
  flat $\\Lambda$CDM with $H_0$ = 70, $\\Omega_m$ = 0.3.
* TDE: 20,000 K blackbody; monitoring needs SNR 20 per R = 300 bin.
* AGN reverberation: 5100 A continuum precision per epoch from the SNR per resel averaged over a
  50 A window (photon noise only; fiber spectrophotometric stability is *not* modelled).
* Sky: the bundled grey-time table (moon up); dark time is ~0.3-0.5 mag deeper at the faint end.
"""),
code(PREAMBLE), code(SCALES),
code("""WINDOW = (450, 650); R_BIN = 300; SNR_CLASS = 5; BAND = "R"
def bench(n):
    return LFASTETC("bench_r1800", telescope=tel(n), wavelength_min_nm=WINDOW[0], wavelength_max_nm=WINDOW[1])
k = 1800/R_BIN
snr_resel_class = SNR_CLASS/np.sqrt(k)
sn_like = Source.blackbody(6000, 0, BAND)"""),
md("## Classification depth versus array size"),
code("""fig, ax = plt.subplots(figsize=(7.5, 4.5))
depth = {}
for exp, color in ((300, 'C2'), (1800, 'C0'), (3600, 'C1')):
    d = np.array([bench(int(n)).limiting_magnitude(sn_like, exp, snr_resel_class, wavelength=WINDOW) for n in NS])
    depth[exp] = dict(zip(NS, d))
    ax.plot(NS, d, marker='o', ms=3, color=color, label=f"{exp/60:.0f} min")
for mag, txt in ((20.5, "ZTF single visit"), (24.5, "LSST single visit r")):
    ax.axhline(mag, color='0.5', ls=':'); ax.text(NS[0], mag-0.1, txt, fontsize=8, color='0.4')
for n in (20, 200, 1720): ax.axvline(n, color='0.6', ls='--', lw=0.8)
ax.set_xscale('log'); ax.set_xlabel("Number of telescopes"); ax.set_ylabel("Classification depth (R mag)")
ax.set_title("SNR=5 per R=300 bin, 450-650 nm, 6000 K SED, grey sky", fontsize=9); ax.legend()
deq_axis(ax)
fig.savefig(FIG + "classification_depth_vs_N.pdf")
for exp, d in depth.items():
    print(f"{exp:5d} s: " + "  ".join(f"N={n}: R<{d[n]:.1f}" for n in (20, 200, 1720)))"""),
md("""## Kilonova horizon

Distance at which an AT2017gfo-like kilonova can still be classified (30 min exposure) at
1, 2 and 3 days after merger, versus the array size. The LIGO-Virgo-KAGRA O5 binary-neutron-star
range is a few hundred Mpc (see memo for the cited values)."""),
code("""def kn_src(day):
    teff = np.interp(day, [0.5, 1.5, 3.5], [10000, 5000, 3000])
    return Source.blackbody(teff, 0, BAND)
def kn_mag_40(day):
    return 17.3 + 0.6*(day-1.0)
fig, ax = plt.subplots(figsize=(7.5, 4.5))
horizon = {}
for day, color in ((1, 'C0'), (2, 'C1'), (3, 'C3')):
    lim = np.array([bench(int(n)).limiting_magnitude(kn_src(day), 1800, snr_resel_class, wavelength=WINDOW) for n in NS])
    D = 40*10**((lim - kn_mag_40(day))/5)
    horizon[day] = dict(zip(NS, D))
    ax.plot(NS, D, marker='o', ms=3, color=color, label=f"{day} d after merger")
for n in (20, 200, 1720): ax.axvline(n, color='0.6', ls='--', lw=0.8)
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel("Number of telescopes"); ax.set_ylabel("Classification horizon (Mpc)")
ax.set_title("AT2017gfo-like kilonova, 30 min, SNR=5 per R=300 bin", fontsize=9); ax.legend()
deq_axis(ax)
fig.savefig(FIG + "kilonova_horizon_vs_N.pdf")
for day, d in horizon.items():
    print(f"day {day}: " + "  ".join(f"N={n}: {d[n]:5.0f} Mpc" for n in (20, 200, 1720)))"""),
md("## SN Ia classification redshift reach (peak and one magnitude past peak)"),
code("""from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
zs = np.linspace(0.01, 1.2, 600)
m_peak = -19.3 + cosmo.distmod(zs).value       # B ~ R at peak for a 6000 K SED, K-correction ignored
fig, ax = plt.subplots(figsize=(7.5, 4.5))
zreach = {}
for exp, color in ((1800, 'C0'), (3600, 'C1')):
    zr = []
    for n in NS:
        lim = bench(int(n)).limiting_magnitude(sn_like, exp, snr_resel_class, wavelength=WINDOW)
        zr.append((np.interp(lim, m_peak, zs), np.interp(lim - 1, m_peak, zs)))
    zr = np.array(zr); zreach[exp] = dict(zip(NS, zr[:, 0]))
    ax.plot(NS, zr[:, 0], marker='o', ms=3, color=color, label=f"{exp/60:.0f} min, at peak")
    ax.plot(NS, zr[:, 1], ls='--', color=color, label=f"{exp/60:.0f} min, peak + 1 mag")
for n in (20, 200, 1720): ax.axvline(n, color='0.6', ls='--', lw=0.8)
ax.set_xscale('log'); ax.set_xlabel("Number of telescopes"); ax.set_ylabel("Max redshift for classification")
ax.set_title("SN Ia, M_B=-19.3, SNR=5 per R=300 bin", fontsize=9); ax.legend(fontsize=8)
fig.savefig(FIG + "snia_zreach_vs_N.pdf")
for exp, d in zreach.items():
    print(f"{exp:5d} s at peak: " + "  ".join(f"N={n}: z<{d[n]:.2f}" for n in (20, 200, 1720)))"""),
md("""## TDE monitoring

Tidal disruption events stay bright (blue, ~20,000 K) for months; weekly spectra with SNR 20 per
R = 300 bin follow the line evolution. Time per epoch for TDEs at g ~ 18, 19 and 20."""),
code("""tde = Source.blackbody(20000, 0, BAND)
snr_resel_mon = 20/np.sqrt(k)
for n in (20, 200, 1720):
    t = bench(n).exptime_vs_mag(tde, [17, 18, 19, 20], snr=snr_resel_mon, wavelength=WINDOW)
    print(f"N={n:5d}: " + "  ".join(f"R={m}: {ti/60:6.1f} min" for m, ti in zip((17, 18, 19, 20), t)))"""),
md("""## AGN reverberation mapping: continuum precision per epoch

Nightly spectrophotometry at the ~1-2% level over months is what reverberation campaigns need.
Photon-noise precision on the 5100 A continuum (50 A window) per 10-minute epoch for Seyferts
at V = 14-17; the fiber-aperture and spectrophotometric stability floor is not included and is
the real question for this programme."""),
code("""agn = Source.powerlaw(-0.5, 0, "V", system="AB")
for n in (20, 200, 1720):
    e = LFASTETC("bench_r1800", telescope=tel(n))
    row = []
    for m in (14, 15, 16, 17):
        r = e.snr(agn.with_mag(m), 600)
        win = (r.wavelength >= 505) & (r.wavelength <= 515)
        snr50 = np.sqrt(np.nansum(r.snr[win]**2))       # co-add the resels in the window
        row.append(f"V={m}: {100/snr50:.2f}%")
    print(f"N={n:5d} (10 min): " + "  ".join(row))"""),
]

# --------------------------------------------------------------------------- 08
NB["08_galactic_archaeology"] = [
md("""
# Galactic archaeology: single-target spectroscopy of halo and metal-poor stars

LFAST is **not** a multiplexed survey instrument (all telescopes point at one star), so its niche
in Galactic archaeology is depth on individual targets selected by Gaia, LAMOST, LSST or
Pristine-like photometry: extremely metal-poor candidates, stream members, hypervelocity stars,
and faint halo giants. At R ~ 1800 the metallicity precision is LAMOST-like (~0.1 dex at SNR
20-30 per pixel; Wu et al. 2014, Luo et al. 2015); the question is how faint and how fast.

**Assumptions** (knobs): SNR 30 per resel (~0.1 dex) in the 500-560 nm window; 2 min overhead per
target; G dwarf $M_V$ = +4.8, K giant $M_V$ = +0.5, red-giant-branch tip $M_V$ = -2.5 for the
distance reach; a 5800 K blackbody for dwarfs and 4500 K for giants; grey-time sky.
"""),
code(PREAMBLE), code(SCALES),
code("""WINDOW = (500, 560); SNR_FEH = 30; OVERHEAD = 120.0
dwarf = Source.blackbody(5800, 0, "V"); giant = Source.blackbody(4500, 0, "V")
def bench(n):
    return LFASTETC("bench_r1800", telescope=tel(n))"""),
md("## Depth: limiting V for [Fe/H] in 15 min and 1 h"),
code("""fig, ax = plt.subplots(figsize=(7.5, 4.5))
lims = {}
for exp, color in ((900, 'C0'), (3600, 'C1')):
    d = np.array([bench(int(n)).limiting_magnitude(giant, exp, SNR_FEH, wavelength=WINDOW) for n in NS])
    lims[exp] = dict(zip(NS, d))
    ax.plot(NS, d, marker='o', ms=3, color=color, label=f"{exp/60:.0f} min")
for n in (20, 200, 1720): ax.axvline(n, color='0.6', ls='--', lw=0.8)
ax.set_xscale('log'); ax.set_xlabel("Number of telescopes"); ax.set_ylabel("Limiting V (SNR=30/resel, 500-560 nm)")
ax.set_title("Bench R~1800, 4500 K giant", fontsize=9); ax.legend()
deq_axis(ax)
fig.savefig(FIG + "metallicity_depth_vs_N.pdf")
for exp, d in lims.items():
    print(f"{exp:5d} s: " + "  ".join(f"N={n}: V<{d[n]:.1f}" for n in (20, 200, 1720)))"""),
md("## Distance reach by tracer (1 h, SNR 30)"),
code("""tracers = {"G dwarf (M_V=+4.8)": (dwarf, 4.8), "K giant (M_V=+0.5)": (giant, 0.5), "RGB tip (M_V=-2.5)": (giant, -2.5)}
print(f"{'tracer':22s}" + "".join(f"{'N=%d' % n:>14s}" for n in (20, 200, 1720)))
for label, (src, MV) in tracers.items():
    row = []
    for n in (20, 200, 1720):
        vlim = bench(n).limiting_magnitude(src, 3600, SNR_FEH, wavelength=WINDOW)
        d_kpc = 10**((vlim - MV + 5)/5)/1000
        row.append(f"{d_kpc:9.1f} kpc")
    print(f"{label:22s}" + "".join(f"{r:>14s}" for r in row))"""),
md("## Survey speed: targets per 8 h night versus magnitude"),
code("""mags = np.linspace(12, 21, 46)
fig, ax = plt.subplots(figsize=(7.5, 4.5))
for n in (20, 200, 1720):
    t = bench(n).exptime_vs_mag(giant, mags, snr=SNR_FEH, wavelength=WINDOW)
    ax.plot(mags, 8*3600/(t + OVERHEAD), label=SCALES[n])
ax.set_yscale('log'); ax.set_xlabel("V magnitude"); ax.set_ylabel("Targets per 8 h night (SNR=30/resel)")
ax.set_title("One target at a time, 2 min overhead, 4500 K giant", fontsize=9); ax.legend(fontsize=8)
fig.savefig(FIG + "survey_speed_vs_N.pdf")
for n in (20, 200, 1720):
    t = bench(n).exptime_vs_mag(giant, [15, 17, 19], snr=SNR_FEH, wavelength=WINDOW)
    print(f"N={n:5d}: " + "  ".join(f"V={m}: {8*3600/(ti+OVERHEAD):6.0f}/night ({ti/60:5.1f} min each)" for m, ti in zip((15, 17, 19), t)))"""),
md("""## Context: Gaia XP metal-poor candidates

Photometric metallicity catalogues from Gaia XP spectra contain ~10$^5$-10$^6$ candidates with
[Fe/H] < -2 at G ~ 13-17.5. At V ~ 16-17 the nominal array confirms a candidate in minutes, so the
20x array alone could verify thousands per season; the larger arrays push confirmation to the
faint, distant halo (V ~ 19-20) where multiplexed 4 m surveys run out of photons per fiber."""),
code("""for n in (20, 200, 1720):
    t = bench(n).exptime_vs_mag(giant, [16, 17.5, 19, 20], snr=SNR_FEH, wavelength=WINDOW)
    print(f"N={n:5d}: " + "  ".join(f"V={m}: {ti/60:6.1f} min" for m, ti in zip((16, 17.5, 19, 20), t)))"""),
]

# --------------------------------------------------------------------------- 09
NB["09_stellar_activity_flares"] = [
md("""
# Time-resolved stellar activity: M-dwarf flares and chromospheric monitoring

A fiber-combined array reads out once per exposure regardless of the number of telescopes, so
**second-to-minute cadence spectroscopy of bright-ish stars** is a natural niche: H$\\alpha$,
H$\\beta$ and He I D3 during M-dwarf flares, chromospheric variability of exoplanet hosts, and
accretion variability in young stars (visible lines here; Paschen-$\\beta$ with CSHELL in the
NIR notebook).

**Assumptions** (knobs): SNR 20 per resel in the continuum next to H$\\alpha$ (650-660 nm) as
the threshold for a usable line profile per exposure; cadences 10 s and 60 s; M dwarfs as
blackbodies of 3300 K normalised in V, with $M_V$ = +10 (M3), +13 (M5) and +16 (M8) for the
distance reach; grey-time sky.
"""),
code(PREAMBLE), code(SCALES),
code("""WINDOW = (650, 660); SNR_LINE = 20
mdwarf = Source.blackbody(3300, 0, "V")
def bench(n):
    return LFASTETC("bench_r1800", telescope=tel(n))"""),
md("## Faintest M dwarf with SNR 20 per resel at Halpha, versus cadence and array size"),
code("""fig, ax = plt.subplots(figsize=(7.5, 4.5))
lims = {}
for exp, color in ((10, 'C3'), (60, 'C0'), (600, 'C1')):
    d = np.array([bench(int(n)).limiting_magnitude(mdwarf, exp, SNR_LINE, wavelength=WINDOW) for n in NS])
    lims[exp] = dict(zip(NS, d))
    ax.plot(NS, d, marker='o', ms=3, color=color, label=f"{exp} s cadence")
for n in (20, 200, 1720): ax.axvline(n, color='0.6', ls='--', lw=0.8)
ax.set_xscale('log'); ax.set_xlabel("Number of telescopes"); ax.set_ylabel("Limiting V (SNR=20/resel at Halpha)")
ax.set_title("Bench R~1800, 3300 K", fontsize=9); ax.legend()
deq_axis(ax)
fig.savefig(FIG + "flare_cadence_depth_vs_N.pdf")
for exp, d in lims.items():
    print(f"{exp:4d} s: " + "  ".join(f"N={n}: V<{d[n]:.1f}" for n in (20, 200, 1720)))"""),
md("## Distance reach at 10 s cadence by spectral type"),
code("""types = {"M3 (M_V=+10)": 10, "M5 (M_V=+13)": 13, "M8 (M_V=+16)": 16}
print(f"{'type':14s}" + "".join(f"{'N=%d' % n:>12s}" for n in (20, 200, 1720)))
for label, MV in types.items():
    row = []
    for n in (20, 200, 1720):
        vlim = lims[10][n]
        row.append(f"{10**((vlim - MV + 5)/5):8.1f} pc")
    print(f"{label:14s}" + "".join(f"{r:>12s}" for r in row))"""),
md("""## Flare spectra: what a 1-magnitude flare looks like at 10 s cadence

During a large flare the blue continuum brightens by magnitudes while H$\\alpha$ goes into
strong emission. A flare is modelled as the quiescent 3300 K star plus a 10,000 K blackbody
component that doubles the V flux (a ~0.75 mag flare), with H$\\alpha$ in emission (depth -3 on
a 0.1 nm Gaussian). SNR per resel across the band in 10 s for a V = 12 M dwarf."""),
code("""wl = np.linspace(380, 1000, 4000)
quiet = Source.blackbody(3300, 12, "V")
hot = Source.blackbody(10000, 12, "V")
flare_flux = quiet.flux(wl) + hot.flux(wl)
flare = Source.from_arrays(wl, flare_flux, flux_unit="W/m2/nm", band="V").add_line(656.3, -3.0, 0.1)
fig, ax = plt.subplots(figsize=(9, 4))
for n, ls in ((20, '-'), (200, '--')):
    e = bench(n)
    ax.plot(e.wavelength, e.snr(quiet, 10).snr, ls=ls, color='C0', lw=0.8, label=f"quiescent, N={n}")
    ax.plot(e.wavelength, e.snr(flare, 10).snr, ls=ls, color='C3', lw=0.8, label=f"flare, N={n}")
ax.set_yscale('log'); ax.set_xlabel("Wavelength (nm)"); ax.set_ylabel("SNR per resel in 10 s")
ax.set_title(f"V=12 M dwarf, flare V={flare.mag:.2f}", fontsize=9); ax.legend(fontsize=8, ncol=2)
fig.savefig(FIG + "flare_snr_10s.pdf")"""),
md("## Exoplanet-host chromospheric monitoring: Ca II H&K is outside the Bench window, Halpha is not"),
code("""hosts = {"V=8 K dwarf": (4800, 8), "V=11 M dwarf": (3300, 11), "V=14 M dwarf": (3300, 14)}
for label, (teff, mag) in hosts.items():
    row = []
    for n in (20, 200, 1720):
        r = bench(n).snr(Source.blackbody(teff, mag, "V"), 60)
        row.append(f"N={n}: SNR={r.at(656.3)['snr']:6.0f}" + (" (sat!)" if r.saturated.any() else ""))
    print(f"{label:14s} 60 s at Halpha: " + "  ".join(row))"""),
]

# --------------------------------------------------------------------------- 10
NB["10_cshell_nir_science"] = [
md("""
# CSHELL science: He 10830 atmospheres, Paschen-beta accretion, and NIR radial velocities

CSHELL is a **bright-object, single-order** instrument: 1.08-1.75 um at R = 18,000-43,000, but one
setting is only ~27-44 A wide (FSR ~ 2.5e-3 lambda). The science that fits is a single strong line
per setting. Three cases are quantified here versus the array size:

1. **He I 10830 transmission spectroscopy** of transiting planets (escaping atmospheres).
2. **Paschen-beta (1282 nm) accretion monitoring** of T Tauri stars (nightly to hourly cadence).
3. **km/s radial velocities** at 1.55 um for Gaia astrometric / eclipsing binaries.

**Assumptions** (knobs): as in notebook 05 for He 10830 (0.7 A FWHM line, 300 s exposures, 2 h
transit, equal out-of-transit time, resel-averaged line core, placeholder slit throughput 0.5 at
R = 43k). Pa-beta: R = 18k, SNR 30 per resel needed for a 10% line-flux measurement of a line with
EW ~ 5 A. RVs: SNR 30 per resel at R = 43k in a 40 A window gives ~1 km/s for a late-type star
(Bouchy et al. 2001 scaling with a conservative quality factor; a knob, see below).
"""),
code(PREAMBLE), code(SCALES),
code("""from scipy.special import erf
EXP = 300.0; N_IN = int(2*3600/EXP); HE_NM = 1083.3; LINE_FWHM_A = 0.7
N_CORE = {"cshell_r43000": 3, "cshell_r18000": 2}
def he(n, preset="cshell_r43000"):
    return LFASTETC(preset, telescope=tel(n), wavelength_min_nm=1081, wavelength_max_nm=1086)
def core_fraction(etc, n_core):
    sigma = LINE_FWHM_A/2.3548/10; half = n_core*HE_NM/etc.instrument.resolving_power/2
    return sigma*np.sqrt(2*np.pi)/(2*half)*erf(half/(sigma*np.sqrt(2)))
def he_sigma(etc, src):
    n_exp = etc.reads_needed(src, EXP)
    snr = etc.snr(src, EXP, n_exp=n_exp).at(HE_NM)["snr"]
    n_core = N_CORE[etc.instrument.name]
    return np.sqrt(2)/(snr*np.sqrt(n_core*N_IN))/core_fraction(etc, n_core)   # sigma on the peak depth"""),
md("## He 10830: faintest host for a 5-sigma single-transit detection"),
code("""jmags = np.linspace(4, 15, 45)
fig, ax = plt.subplots(figsize=(7.5, 4.5))
reach = {}
for depth, color in ((0.005, 'C0'), (0.01, 'C1'), (0.03, 'C3')):
    lim = []
    for n in NS:
        etc = he(int(n))
        sig = np.array([he_sigma(etc, Source.blackbody(4800, j, "J")) for j in jmags])
        ok = jmags[depth/sig >= 5]
        lim.append(ok.max() if len(ok) else np.nan)
    reach[depth] = dict(zip(NS, lim))
    ax.plot(NS, lim, marker='o', ms=3, color=color, label=f"{100*depth:g}% line")
for n in (20, 200, 1720): ax.axvline(n, color='0.6', ls='--', lw=0.8)
ax.set_xscale('log'); ax.set_xlabel("Number of telescopes"); ax.set_ylabel("Faintest host J (5 sigma, 1 transit)")
ax.set_title("CSHELL R=43k, He 10830, 4800 K host, 2 h transit", fontsize=9); ax.legend()
deq_axis(ax)
fig.savefig(FIG + "he10830_reach_vs_N.pdf")
for depth, d in reach.items():
    print(f"{100*depth:4.1f}% line: " + "  ".join(f"N={n}: J<{d[n]:.1f}" for n in (20, 200, 1720)))"""),
md("## Paschen-beta accretion monitoring of T Tauri stars (Taurus, 140 pc: J ~ 8-12)"),
code("""def pab(n):
    return LFASTETC("cshell_r18000", telescope=tel(n), wavelength_min_nm=1280, wavelength_max_nm=1284)
for n in (20, 200, 1720):
    e = pab(n)
    t = e.exptime_vs_mag(Source.blackbody(4000, 0, "J"), [8, 10, 12, 14], snr=30, wavelength=1282.2)
    print(f"N={n:5d}: " + "  ".join(f"J={m}: {ti:7.0f} s" for m, ti in zip((8, 10, 12, 14), t)))"""),
md("""## Radial velocities at 1.55 um, R = 43,000

Photon-limited RV precision scales as $\\sigma_{RV} = c / (Q \\sqrt{N_{e}})$ (Bouchy et al. 2001)
with $N_e$ the total electrons in the spectrum and $Q$ the spectral quality factor. For a K/M dwarf
at R = 43k in a 40 A H-band window we take $Q$ = 300 (knob; Bouchy et al. give Q ~ 10$^4$ at
R = 100k across the full optical range, and a 40 A window with few lines is far poorer). The
result is only a scaling: km/s orbits for binaries are easy, m/s planet work is not."""),
code("""Q = 300.0
def rv(n):
    return LFASTETC("cshell_r43000", telescope=tel(n), wavelength_min_nm=1550, wavelength_max_nm=1590)
for n in (20, 200, 1720):
    e = rv(n); row = []
    for h in (6, 9, 12):
        r = e.snr(Source.blackbody(4000, h, "H"), 900, n_exp=e.reads_needed(Source.blackbody(4000, h, "H"), 900))
        ne = np.nansum(r.n_star*r.exptime)        # electrons across the window
        row.append(f"H={h}: {2.998e5/(Q*np.sqrt(ne)):7.2f} km/s")
    print(f"N={n:5d} (15 min): " + "  ".join(row))"""),
md("## Bright end: reads needed for a J = 3 star in 300 s (placeholder 80 ke- full well)"),
code("""for n in (20, 200, 1720):
    print(f"N={n:5d}: reads = {he(n, 'cshell_r18000').reads_needed(Source.blackbody(4800, 3, 'J'), 300)}")"""),
]


if __name__ == "__main__":
    names = sys.argv[1:] or list(NB)
    for n in names:
        build(n, NB[n])
