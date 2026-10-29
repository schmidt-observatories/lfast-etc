"""Build the array-scaling science-survey notebooks (06-10).

Usage: .venv/bin/python scripts/build_science_notebooks.py [name ...]

Figures are also written as PDF into ../memo/figures/ when that directory exists (the
author's memo build); otherwise the save is skipped, so the notebooks run from a plain clone.
"""
import sys
from build_notebooks import build, md, code, PREAMBLE

SCALES = '''# ---- knobs: array sizes compared throughout (replicated in arrays of 20) ----
from pathlib import Path
SCALES = {20: "nominal 20x (3.2 m eq.)", 200: "10 arrays (10.3 m eq.)", 1720: "86 arrays (30.1 m eq.)"}
NS = np.array([20, 40, 60, 100, 200, 400, 600, 1000, 1720, 2500])
FIG = Path("../memo/figures")          # memo figures land here if the directory exists
def save(fig, name):
    if FIG.is_dir():
        fig.savefig(FIG / name)
plt.rcParams.update({"font.size": 11, "axes.titlesize": 10, "legend.fontsize": 9})
def tel(n, telescopes_per_input=20):
    return Telescope(n_telescopes=int(n), fiber="18um", combine="fiber", telescopes_per_input=telescopes_per_input)
def mark_scales(ax):
    for n in (20, 200, 1720):
        ax.axvline(n, color='0.6', ls='--', lw=0.8)
    ax.set_xscale('log'); ax.set_xlabel("Number of telescopes")
    equivalent_diameter_axis(ax, ticks=(3, 5, 10, 20, 30))'''

NB = {}

# --------------------------------------------------------------------------- 06
NB["06_array_scaling_overview"] = [
md("""
# Array scaling overview: 20 telescopes, 10 m equivalent, 30 m equivalent

LFAST is replicated in arrays of 20 unit telescopes. The nominal Kitt Peak array (20 x 0.726 m,
8.3 m$^2$) has the collecting area of a 3.2 m telescope; 10 arrays (200 telescopes) match a
10.3 m aperture and 86 arrays (1720 telescopes) a 30.1 m aperture. This notebook sets up the
scaling model used by notebooks 07 to 10 and produces the headline sensitivity curves.

How arrays combine in the ETC: each 20-to-1 fiber combiner feeds one spectrograph input. Within an
array, star and sky electrons both scale with the number of telescopes while detector noise is
incurred once, so the SNR gain is better than $\\sqrt{N}$ when read noise or dark current matter.
Beyond 20, the ETC assumes one input per array, co-added incoherently, so the SNR grows as
$\\sqrt{n_{arrays}}$ on top of the single-array value (`Telescope(..., telescopes_per_input=20)`).
An idealised single combiner (`telescopes_per_input=None`) is shown dotted where it matters.

This is a photon-budget scaling and nothing more. It assumes the spectrographs scale with the
arrays (one input per array on the same or replicated spectrographs), that the per-array
throughput stays at the Choi et al. (2026) values, and that sky subtraction, flat fielding and
fiber-combiner stability do not limit faint-end performance. Multiplexing is not modelled: all
telescopes point at one target at a time.
"""),
code(PREAMBLE), code(SCALES),
md("## Equivalent apertures"),
code("""apertures = pd.DataFrame(
    {int(n): {"Arrays": -(-int(n)//20), "Area (m²)": tel(int(n)).total_area_m2,
              "D_eq (m)": tel(int(n)).equivalent_diameter_m, "Spectrographs": tel(int(n)).n_inputs}
     for n in NS}).T.rename_axis("N telescopes")
apertures.astype({"Arrays": int, "Spectrographs": int}).style.format(
    {"Area (m²)": "{:.1f}", "D_eq (m)": "{:.1f}"}).apply(
    lambda r: ["font-weight: bold" if r.name in SCALES else "" for _ in r], axis=1)"""),
md("Bold rows are the three scales compared throughout: the nominal 20x array, 10 m equivalent and 30 m equivalent. In this model every array has its own spectrograph, so the number of spectrographs equals the number of arrays."),
md("""
## Limiting magnitude versus array size, both instruments

Three cases: the Bench in V over 500 to 650 nm, CSHELL at R = 18k in H for one setting at 1.555 um,
and CSHELL at R = 43k in J at the He 10830 line. Solid lines use one combiner per 20 telescopes;
dotted lines are the idealised single combiner. The gap between them is the price of paying the
detector noise once per array instead of once in total, and it is largest for CSHELL, whose
detector is noisier.
"""),
code("""gstar = Source.blackbody(5800, 0, "V")
kstar = Source.blackbody(4500, 0, "H")
cs_lo, cs_hi = Instrument.from_preset("cshell_r18000").setting_window(1555.0)   # one CSHELL echelle setting
cases = {
    "Bench R~1800, V, 500-650 nm": ("bench_r1800", gstar, dict(), (500, 650)),
    "CSHELL R=18k, H, one setting at 1.555 um": ("cshell_r18000", kstar, dict(wavelength_min_nm=cs_lo, wavelength_max_nm=cs_hi), None),
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
    ax.set_ylabel(f"Limiting {src.band} (SNR=10/resel)")
    ax.set_title(label + "\\n(dotted: ideal single combiner)", fontsize=9); ax.legend(fontsize=8)
    mark_scales(ax)
plt.tight_layout(); save(fig, "limiting_mag_vs_N.pdf")"""),
code("""lim_table = pd.DataFrame({(label, exp): {f"N={n}": d[n] for n in SCALES} for (label, exp), d in summary.items()}).T
lim_table.index.names = ["case (SNR=10 per resel)", "exposure (s)"]
lim_table["gain 20->200"] = lim_table["N=200"] - lim_table["N=20"]
lim_table["gain 20->1720"] = lim_table["N=1720"] - lim_table["N=20"]
lim_table.style.format("{:.2f}")"""),
md("""* **Bench.** The gain is largest for short exposures: a V ~ 16.3 star in 60 s gives ~150 e$^-$ per
  resel against ~14 e$^-$ of sky, so the star dominates and the gain is 1.6 mag (20 -> 200). In
  1 h at V ~ 19.9 the sky (~820 e$^-$) outweighs the star (~350 e$^-$) and the gain drops to 1.4 mag.
* **CSHELL.** Detector noise dominates at every exposure (22 e$^-$ read noise and 0.5 e$^-$/s dark
  current over the 36-pixel extraction aperture: in 15 min at the faint limit the star gives
  ~1900 e$^-$ per resel against ~16,000 e$^-$ of dark and a read-noise variance of ~17,000 e$^{-2}$),
  so the gain stays at ~1.3 mag (20 -> 200) and ~2.4 mag (20 -> 1720), close to the
  background-limited value.
* **Single combiner (dotted).** Feeding all telescopes into one spectrograph would help the Bench
  only a little (~0.1 mag) but would help CSHELL a lot: ~1.2 mag deeper at 200 telescopes and
  ~2.3 mag at 1720 in 15 min, because its detector noise would then be paid only once."""),
md("""
## How the gain is split: source-, sky- and read-noise-limited regimes

At fixed SNR the required exposure time scales as $1/N$ in every regime once the array count
grows (incoherent co-add of inputs), and faster than $1/N$ inside a single array when detector
noise matters. The magnitude gain at fixed time is $2.5 \\log_{10} N$ when source limited and
$1.25 \\log_{10} N$ when background limited, so the faint-end gain from 20 to 1720 telescopes is
between 2.4 and 4.8 mag.
"""),
code("""fig, ax = plt.subplots(figsize=(7, 4.3))
for mag, label in ((12, "V=12 (source-limited)"), (18, "V=18"), (21, "V=21 (sky-limited)")):
    t = [LFASTETC("bench_r1800", telescope=tel(int(n))).exptime_vs_mag(gstar, [mag], snr=10, wavelength=(500, 650))[0] for n in NS]
    ax.plot(NS, np.array(t)/t[0], marker='o', ms=3, label=label)
ax.plot(NS, 20/NS, color='0.4', ls='--', label="1/N")
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel("Number of telescopes"); ax.set_ylabel("Exposure time relative to 20x")
ax.legend(fontsize=8); ax.set_title("Bench, SNR=10 per resel, 5800 K", fontsize=9)
save(fig, "exptime_scaling.pdf")"""),
md("""
## Bright end: saturation

With one combiner input per 20-telescope array, the light per detector pixel is set by a single
array and does not grow with the number of arrays, so the bright limit is the same at 20, 200
and 1720 telescopes. An idealised single combiner would instead saturate $N/20$ times sooner.
The table gives the brightest star for which a 1 s (Bench) or 8 s (CSHELL, Fowler minimum) read
stays below 80% of full well (STA1 87.5 ke$^-$ for the Bench, 55 ke$^-$ for CSHELL).
"""),
code("""bright = {}
for preset, band, teff, kw, tmin in (("bench_r1800", "V", 5800, {}, 1.0), ("cshell_r18000", "H", 4500, dict(wavelength_min_nm=cs_lo, wavelength_max_nm=cs_hi), 8.0)):
    row = {}
    for n in SCALES:
        e = LFASTETC(preset, telescope=tel(n), **kw)
        mags = np.arange(-2, 14, 0.1)
        peak = np.array([np.nanmax(e.snr(Source.blackbody(teff, m, band), tmin).peak_e_per_pixel) for m in mags])
        row[f"N={n}"] = f"{band} > {mags[peak < 0.8*e.instrument.full_well].min():.1f}"
    bright[f"{preset} ({tmin:.0f} s read)"] = row
pd.DataFrame(bright).T.rename_axis("brightest unsaturated star")"""),
md("""The reads here are already the shortest each detector allows (1 s for the Bench, the 8 s Fowler
minimum for CSHELL), so a brighter star has to be attenuated, e.g. by using fewer telescopes per
array or a neutral-density filter."""),
]

# --------------------------------------------------------------------------- 07
NB["07_time_domain_lsst_era"] = [
md("""
# Time-domain science in the LSST era

A dedicated spectroscopic array is most valuable where the bottleneck is telescope time on single
targets: transient classification, rapid multi-messenger response and high-cadence monitoring.
This notebook explores what magnitude and distance limits the nominal 20-telescope array reaches, and
how they move at 10 m and 30 m equivalent apertures.

Assumptions:

* Classification: SNR 5 per R = 300 bin over 450 to 650 nm (SEDM-like, Blagorodnova et al.
  2018), which is SNR 2.04 per native resel; magnitudes in R (Cousins, close to SDSS r).
* Kilonova: AT2017gfo-like, r = 17.3 at 0.5 d at 40 Mpc, fading 0.95 mag per day in r, cooling
  10,000 K to 5000 K to 2500 K at 0.5, 1.5 and 5.5 d (Drout et al. 2017; Shappee et al. 2017).
* SN Ia: $M_B$ = -19.3 at peak, a 10,000 K blackbody ($B-V \\approx 0$ at peak, so $B \\approx R$), no K-correction; luminosity distance from a
  flat $\\Lambda$CDM with $H_0$ = 70, $\\Omega_m$ = 0.3.
* TDE: 20,000 K blackbody (TDE optical temperatures are typically a few x 10$^4$ K); monitoring needs SNR 20 per R = 300 bin.
* AGN reverberation: 5100 A continuum precision per epoch from the SNR per resel co-added over a
  50 A window, photon noise only. Fiber spectrophotometric stability is not modelled.
* Sky: the bundled grey-time table (moon up). Dark time is 0.3 to 0.5 mag deeper at the faint end.
"""),
code(PREAMBLE), code(SCALES),
code("""WINDOW = (450, 650); R_BIN = 300; SNR_CLASS = 5; BAND = "R"
def bench(n):
    return LFASTETC("bench_r1800", telescope=tel(n), wavelength_min_nm=WINDOW[0], wavelength_max_nm=WINDOW[1])
k = 1800/R_BIN
snr_resel_class = SNR_CLASS/np.sqrt(k)
sn_like = Source.blackbody(10000, 0, BAND)   # SN Ia near peak, B-V ~ 0"""),
md("""## Classification depth versus array size

How faint of a transient can the array classify in a given exposure, and how does that change with the
array size? The figure below shows the classification depth (SNR = 5 per R = 300 bin over 450-650 nm) for a
SN Ia-like 10,000 K SED in 5, 30 and 60 minutes, against the number of telescopes (equivalent
diameter on the top axis). The horizontal lines mark the discovery depths of the surveys LFAST would
follow up:

* the Argus Array single-epoch (1 min, g = 20.0) and 15-minute co-add (g = 21.5) depths
  (science requirements, 5 sigma, dark time; argus.unc.edu/specifications; Law et al. 2022);
* the LSST single-visit depth (r ~ 24.5; Ivezic et al. 2019).

These are detection limits in g or r (AB). For SN Ia-like colours near peak
$g \\approx r$, and R (Vega) is ~0.2 mag brighter than r (AB); both offsets are ignored here.

A classifier that cannot reach a survey's discovery depth leaves its faintest transients
unclassified, so the distance between the curves and these lines is the number to watch."""),
code("""REFS = {"Argus 1 min (g)": 20.0, "Argus 15 min (g)": 21.5, "LSST single visit (r)": 24.5}
fig, ax = plt.subplots(figsize=(7.5, 4.5))
depth = {}
for exp, color in ((300, 'C2'), (1800, 'C0'), (3600, 'C1')):
    d = np.array([bench(int(n)).limiting_magnitude(sn_like, exp, snr_resel_class, wavelength=WINDOW) for n in NS])
    depth[exp] = dict(zip(NS, d))
    ax.plot(NS, d, marker='o', ms=3, color=color, label=f"{exp/60:.0f} min")
for txt, mag in REFS.items():
    ax.axhline(mag, color='0.5', ls=':'); ax.text(NS[-1], mag+0.05, txt, fontsize=8, color='0.4', ha='right')
ax.set_ylabel("Classification depth (R mag)")
ax.set_title("SNR=5 per R=300 bin, 450-650 nm, 10,000 K SED, grey sky", fontsize=9); ax.legend(loc='upper left')
mark_scales(ax)
save(fig, "classification_depth_vs_N.pdf")
pd.DataFrame({f"N={n}": {exp // 60: depth[exp][n] for exp in depth} for n in SCALES}).rename_axis("exposure (min)").style.format("{:.1f}").set_caption("Classification depth (R mag, SNR=5 per R=300 bin)")"""),
md("""The nominal array classifies to R ~ 21.1 in 30 min and R ~ 21.5 in 60 min, i.e. down to the Argus
15-minute co-add depth only with 1 hr exposure time. Ten arrays (200 telescopes) nearly reach that depth
in 5 min (R ~ 21.1) and go 1.3 mag beyond it in 60 min (R ~ 22.8); 86 arrays reach R ~ 24.0 in
60 min, within ~0.5 mag of the LSST single-visit depth."""),
md("""## Kilonova horizon

We can also investigate the distance at which an AT2017gfo-like kilonova can still be classified (30 min exposure) at
0.5, 1, 2 and 3 days after merger, as a function of the array size. The shaded bands are the LIGO binary-neutron-star ranges:
150-170 Mpc in O4 and 225-295 Mpc planned for O5 (LVK observing plan G2002127-v37, Oct 2026).
The right-hand axis gives the same distances as a distance modulus; the legend gives, for each
epoch, the kilonova's apparent R magnitude at the 20-telescope horizon (its classification limit)."""),
code("""def kn_src(day):
    teff = np.interp(day, [0.5, 1.5, 5.5], [10000, 5000, 2500])
    return Source.blackbody(teff, 0, BAND)
def kn_mag_40(day):
    return 17.3 + 0.95*(day-0.5)
fig, ax = plt.subplots(figsize=(7.5, 4.5))
horizon = {}
for day, color in ((0.5, 'C2'), (1, 'C0'), (2, 'C1'), (3, 'C3')):
    lim = np.array([bench(int(n)).limiting_magnitude(kn_src(day), 1800, snr_resel_class, wavelength=WINDOW) for n in NS])
    D = 40*10**((lim - kn_mag_40(day))/5)
    horizon[day] = dict(zip(NS, D))
    ax.plot(NS, D, marker='o', ms=3, color=color, label=f"{day} d after merger (R < {lim[0]:.1f} at N=20)")
ax.axhspan(225, 295, color='0.85', zorder=0); ax.text(NS[-1], 255, "LIGO O5 BNS range", fontsize=8, color='0.4', ha='right', va='center')
ax.axhspan(150, 170, color='0.92', zorder=0); ax.text(NS[-1], 159, "O4 BNS range", fontsize=8, color='0.4', ha='right', va='center')
ax.set_yscale('log'); ax.set_ylabel("Classification horizon (Mpc)")
from matplotlib.ticker import FixedLocator, FixedFormatter, NullFormatter, NullLocator
ticks_mpc = [60, 100, 150, 200, 300, 500, 800]
ax.yaxis.set_major_locator(FixedLocator(ticks_mpc)); ax.yaxis.set_major_formatter(FixedFormatter([str(t) for t in ticks_mpc]))
ax.yaxis.set_minor_formatter(NullFormatter())
mu_ax = ax.secondary_yaxis('right', functions=(lambda d: 5*np.log10(np.clip(d, 1e-3, None)) + 25,
                                               lambda mu: 10**((mu - 25)/5)))
mu_ax.yaxis.set_major_locator(FixedLocator(range(34, 40))); mu_ax.yaxis.set_major_formatter(FixedFormatter([str(m) for m in range(34, 40)]))
mu_ax.yaxis.set_minor_locator(NullLocator()); mu_ax.set_ylabel("Distance modulus (mag)")
ax.set_title("AT2017gfo-like kilonova, 30 min, SNR=5 per R=300 bin", fontsize=9); ax.legend(fontsize=8, loc='upper left')
mark_scales(ax)
save(fig, "kilonova_horizon_vs_N.pdf")
pd.DataFrame({f"N={n}": {f"{day:g}": horizon[day][n] for day in horizon} for n in SCALES}).rename_axis("days after merger").style.format("{:.0f}").set_caption("Kilonova classification horizon (Mpc, 30 min exposure)")"""),
md("""The nominal array classifies a kilonova out to ~230 Mpc at 0.5 d and ~180 Mpc at 1 d, which covers
the O4 range on the first night. Ten arrays push the 1-day horizon to ~320 Mpc, beyond the planned O5
(single facility) range. The kilonova fades and cools fast: at every array size the 3-day horizon is about a third of
the 1-day one."""),
md("""
## SN Ia classification redshift reach

Maximum redshift at which a SN Ia can be classified at peak (solid) and one magnitude past peak
(dashed). B is taken equal to R ($B-V \\approx 0$ at peak for the 10,000 K SED) and the K-correction is
ignored, so treat the
redshifts beyond about 0.3 as indicative.
"""),
code("""from astropy.cosmology import FlatLambdaCDM
from matplotlib.lines import Line2D
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
zs = np.linspace(0.01, 1.2, 600)
m_peak = -19.3 + cosmo.distmod(zs).value
fig, ax = plt.subplots(figsize=(7.5, 4.5))
zreach = {}
for exp, color in ((1800, 'C0'), (3600, 'C1')):
    zr = []
    for n in NS:
        lim = bench(int(n)).limiting_magnitude(sn_like, exp, snr_resel_class, wavelength=WINDOW)
        zr.append((np.interp(lim, m_peak, zs), np.interp(lim - 1, m_peak, zs)))
    zr = np.array(zr)
    zreach[(exp // 60, "at peak")] = dict(zip(NS, zr[:, 0]))
    zreach[(exp // 60, "peak + 1 mag")] = dict(zip(NS, zr[:, 1]))
    ax.plot(NS, zr[:, 0], marker='o', ms=3, color=color)
    ax.plot(NS, zr[:, 1], ls='--', color=color)
ax.set_ylabel("Max redshift for classification")
ax.set_title("SN Ia, M_B=-19.3, SNR=5 per R=300 bin", fontsize=9)
colour_legend = ax.legend([Line2D([], [], color='C0'), Line2D([], [], color='C1')], ["30 min", "60 min"],
                          title="exposure (colour)", fontsize=8, title_fontsize=8, loc='upper left')
ax.add_artist(colour_legend)
ax.legend([Line2D([], [], color='0.3', marker='o', ms=3), Line2D([], [], color='0.3', ls='--')], ["at peak", "1 mag past peak"],
          title="phase (line style)", fontsize=8, title_fontsize=8, loc='lower right')
mark_scales(ax)
save(fig, "snia_zreach_vs_N.pdf")
z_table = pd.DataFrame({f"N={n}": {key: zreach[key][n] for key in zreach} for n in SCALES}).round(2)
z_table.index.names = ["exposure (min)", "phase"]
z_table.style.format("{:.2f}").set_caption("SN Ia: maximum redshift for classification")"""),
md("""In 30 minutes the nominal array classifies SNe Ia at peak to z ~ 0.24 (z ~ 0.16 one magnitude past
peak). Ten arrays reach z ~ 0.40 and 86 arrays z ~ 0.62; doubling the exposure to 60 minutes adds
only ~0.05-0.1 in redshift."""),
md("""## TDE monitoring

Tidal disruption events stay bright and blue (~20,000 K) for months, and their broad H and He II
lines evolve over weeks, so a useful programme takes a spectrum every week or so. We ask for SNR 20
per R = 300 bin, enough to follow the shape and strength of the broad lines, and compute the time
per epoch for TDEs at R = 18-21."""),
code("""tde = Source.blackbody(20000, 0, BAND)
snr_resel_mon = 20/np.sqrt(k)
tde_mags = (18, 19, 20, 21)
pd.DataFrame({f"N={n}": dict(zip([f"R={m}" for m in tde_mags],
                                 bench(n).exptime_vs_mag(tde, tde_mags, snr=snr_resel_mon, wavelength=WINDOW)/60))
              for n in SCALES}).T.rename_axis("array size").style.format("{:.1f}").set_caption("TDE monitoring: exposure time per epoch (min, SNR=20 per R=300 bin)")"""),
md("""For the nominal array a TDE at R = 18-19 costs 5-14 minutes per epoch, so weekly spectra of several
TDEs are cheap. At R = 20 an epoch takes about 55 minutes, and R = 21 (~4.5 h) is out of reach. Ten
arrays bring R = 20 down to ~8 minutes and R = 21 to ~30 minutes; 86 arrays do R = 21 in ~6 minutes,
which opens up the much larger sample of fainter TDEs that LSST will find."""),
md("""
## AGN reverberation mapping: continuum precision per epoch

Reverberation campaigns need nightly spectrophotometry at the 1 to 2% level over months. The
table gives the photon-noise precision on the 5100 A continuum (50 A window) per 10 minute epoch
for Seyferts at V = 14 to 17.
"""),
code("""agn = Source.powerlaw(-0.5, 0, "V", system="AB")
WIN_A = 50; CENTER_NM = 510.0                               # the "50 A window at 5100 A"
lo, hi = CENTER_NM - WIN_A/20, CENTER_NM + WIN_A/20          # 10 A per nm, half-width each side
agn_mags = (14, 15, 16, 17)
rows = {}
for n in SCALES:
    e = LFASTETC("bench_r1800", telescope=tel(n))
    row = {}
    for m in agn_mags:
        r = e.snr(agn.with_mag(m), 600)
        win = (r.wavelength >= lo) & (r.wavelength <= hi)
        # information sum over the resels: the SNR of one common amplitude fitted against the
        # known continuum shape, (sum_i SNR_i^2)^(1/2); not the SNR of an unweighted flux sum
        snr50 = np.sqrt(np.nansum(r.snr[win]**2))
        row[f"V={m}"] = 100/snr50
    rows[f"N={n}"] = row
pd.DataFrame(rows).T.rename_axis("array size").style.format("{:.2f}").set_caption("AGN: 5100 A continuum precision per 10-min epoch (%, photon noise only)")"""),
md("""Photon noise is never the limit here: even the nominal array reaches 0.6% per 10 minutes at
V = 16, several times better than campaigns need. What limits reverberation mapping is the
night-to-night flux calibration. Fiber-coupling and seeing losses change from night to night, and
campaigns normally remove them by scaling every spectrum to the narrow [O III] 5007 A line, whose flux
is constant over the campaign (van Groningen & Wanders 1992). How well that works for a
fiber-combined array is the real question for this programme, but is beyond the scope of this notebook."""),
]

# --------------------------------------------------------------------------- 08
NB["08_galactic_archaeology"] = [
md("""
# Galactic archaeology: single-target spectroscopy of halo and metal-poor stars

LFAST is not a multiplexed survey instrument; all telescopes point at one star. Its niche in
Galactic archaeology is therefore depth on individual targets selected by Gaia, LAMOST, LSST or
Pristine-like photometry: extremely metal-poor candidates, stream members, hypervelocity stars
and faint halo giants. At R ~ 1800 the metallicity precision is LAMOST-like, about 0.1 dex at
SNR 20 to 30 per pixel (Wu et al. 2014; Luo et al. 2015). The question is how faint and how fast.

Assumptions: SNR 30 per resel (about 0.1 dex) in the 500 to 560 nm window; 2 min overhead per
target; G dwarf $M_V$ = +4.8, K giant $M_V$ = +0.5 and red-giant-branch tip $M_V$ = -2.5 for the
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
ax.set_ylabel("Limiting V (SNR=30/resel, 500-560 nm)")
ax.set_title("Bench R~1800, 4500 K giant", fontsize=9); ax.legend()
mark_scales(ax)
save(fig, "metallicity_depth_vs_N.pdf")
for exp, d in lims.items():
    print(f"{exp:5d} s: " + "  ".join(f"N={n}: V<{d[n]:.1f}" for n in (20, 200, 1720)))"""),
md("""
## Distance reach by tracer (1 h, SNR 30)

The limiting magnitude turned into a distance through the absolute magnitude of each tracer. A
giant at the tip of the red giant branch is 7 mag brighter than a G dwarf, so the same spectrum
reaches 25 times farther.
"""),
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
save(fig, "survey_speed_vs_N.pdf")
for n in (20, 200, 1720):
    t = bench(n).exptime_vs_mag(giant, [15, 17, 19], snr=SNR_FEH, wavelength=WINDOW)
    print(f"N={n:5d}: " + "  ".join(f"V={m}: {8*3600/(ti+OVERHEAD):6.0f}/night ({ti/60:5.1f} min each)" for m, ti in zip((15, 17, 19), t)))"""),
md("""
## Context: Gaia XP metal-poor candidates

Photometric metallicity catalogues from Gaia XP spectra (published to G < 17.65) contain about
2 x 10$^5$ candidates with [Fe/H] < -2 (Xylakis-Dornbusch et al. 2024) and more than 10$^4$ very
metal-poor giants (Andrae et al. 2023). At V ~ 16 to 17 the nominal array confirms a candidate in
minutes, so the 20-telescope array alone could verify thousands per season. The larger arrays
push confirmation to the faint, distant halo (V ~ 19 to 20), where multiplexed 4 m surveys run
out of photons per fiber.
"""),
code("""for n in (20, 200, 1720):
    t = bench(n).exptime_vs_mag(giant, [16, 17.5, 19, 20], snr=SNR_FEH, wavelength=WINDOW)
    print(f"N={n:5d}: " + "  ".join(f"V={m}: {ti/60:6.1f} min" for m, ti in zip((16, 17.5, 19, 20), t)))"""),
]

# --------------------------------------------------------------------------- 09
NB["09_stellar_activity_flares"] = [
md("""
# Time-resolved stellar activity: M-dwarf flares and chromospheric monitoring

A fiber-combined array reads out once per exposure regardless of the number of telescopes, so
second-to-minute cadence spectroscopy of bright-ish stars is a natural niche: H$\\alpha$, H$\\beta$
and He I D3 during M-dwarf flares, chromospheric variability of exoplanet hosts, and accretion
variability in young stars (visible lines here; Paschen-$\\beta$ with CSHELL in notebook 10).

Assumptions: SNR 20 per resel in the continuum next to H$\\alpha$ (650 to 660 nm) as the
threshold for a usable line profile per exposure; cadences of 10 s and 60 s; M dwarfs as 3300 K
blackbodies normalised in V, with $M_V$ = +10 (M3), +13 (M5) and +16 (M8) for the distance
reach; grey-time sky.
"""),
code(PREAMBLE), code(SCALES),
code("""WINDOW = (650, 660); SNR_LINE = 20
mdwarf = Source.blackbody(3300, 0, "V")
def bench(n):
    return LFASTETC("bench_r1800", telescope=tel(n))"""),
md("## Faintest M dwarf with SNR 20 per resel at H-alpha, versus cadence and array size"),
code("""fig, ax = plt.subplots(figsize=(7.5, 4.5))
lims = {}
for exp, color in ((10, 'C3'), (60, 'C0'), (600, 'C1')):
    d = np.array([bench(int(n)).limiting_magnitude(mdwarf, exp, SNR_LINE, wavelength=WINDOW) for n in NS])
    lims[exp] = dict(zip(NS, d))
    ax.plot(NS, d, marker='o', ms=3, color=color, label=f"{exp} s cadence")
ax.set_ylabel("Limiting V (SNR=20/resel at H-alpha)")
ax.set_title("Bench R~1800, 3300 K", fontsize=9); ax.legend()
mark_scales(ax)
save(fig, "flare_cadence_depth_vs_N.pdf")
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
md("""
## Flare spectra: a one-magnitude flare at 10 s cadence

During a large flare the blue continuum brightens by magnitudes while H$\\alpha$ goes into strong
emission. The flare is modelled as the quiescent 3300 K star plus a 10,000 K blackbody component
that doubles the V flux (a 0.75 mag flare), with H$\\alpha$ in emission (depth -3 on a 0.1 nm
Gaussian). The plot shows the SNR per resel across the band in 10 s for a V = 12 M dwarf, quiet and
in flare. The blue end gains the most, which is where the flare physics (the Balmer jump and the
hot continuum) lives.
"""),
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
save(fig, "flare_snr_10s.pdf")"""),
md("""
## Exoplanet-host chromospheric monitoring

Ca II H and K are outside the Bench window; H$\\alpha$ is inside it. SNR per resel at H$\\alpha$ in
a 60 s exposure for three typical hosts. The bright K dwarf saturates the STA1 full well
with 20 telescopes already, so it would need shorter reads.
"""),
code("""hosts = {"V=8 K dwarf": (4800, 8), "V=11 M dwarf": (3300, 11), "V=14 M dwarf": (3300, 14)}
for label, (teff, mag) in hosts.items():
    row = []
    for n in (20, 200, 1720):
        r = bench(n).snr(Source.blackbody(teff, mag, "V"), 60)
        row.append(f"N={n}: SNR={r.at(656.3)['snr']:6.0f}" + (" (sat!)" if r.saturated.any() else ""))
    print(f"{label:14s} 60 s at H-alpha: " + "  ".join(row))"""),
]

# --------------------------------------------------------------------------- 10
NB["10_cshell_nir_science"] = [
md("""
# CSHELL science: He 10830 atmospheres, Paschen-beta accretion, and NIR radial velocities

CSHELL is a bright-object, single-order instrument: 1.08 to 1.75 um at R = 18,000 to 43,000, but
one setting is only 27 to 44 A wide (FSR about 2.5e-3 of the wavelength). The science that fits is
a single strong line per setting. Three cases are quantified here versus the array size:

1. He I 10830 transmission spectroscopy of transiting planets (escaping atmospheres).
2. Paschen-beta (1282 nm) accretion monitoring of T Tauri stars (nightly to hourly cadence).
3. km/s radial velocities at 1.55 um for Gaia astrometric and eclipsing binaries.

Assumptions: as in notebook 05 for He 10830 (0.7 A FWHM line, 300 s exposures, 2 h transit, equal
out-of-transit time, resel-averaged line core, placeholder slit throughput 0.5 at R = 43k).
Pa-beta: R = 18k, SNR 30 per resel for a 10% line-flux measurement of a line with EW ~ 5 A. RVs:
a photon-only Bouchy et al. (2001) bound over one 40 A setting at R = 43k with an assumed
quality factor, set below.
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
md("## He 10830: faintest host for a 5 sigma single-transit detection"),
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
ax.set_ylabel("Faintest host J (5 sigma, 1 transit)")
ax.set_title("CSHELL R=43k, He 10830, 4800 K host, 2 h transit", fontsize=9); ax.legend()
mark_scales(ax)
save(fig, "he10830_reach_vs_N.pdf")
for depth, d in reach.items():
    print(f"{100*depth:4.1f}% line: " + "  ".join(f"N={n}: J<{d[n]:.1f}" for n in (20, 200, 1720)))"""),
md("## Paschen-beta accretion monitoring of T Tauri stars (Taurus, 140 pc: J ~ 8 to 12)"),
code("""def pab(n):
    return LFASTETC("cshell_r18000", telescope=tel(n), wavelength_min_nm=1280, wavelength_max_nm=1284)
for n in (20, 200, 1720):
    e = pab(n)
    t = e.exptime_vs_mag(Source.blackbody(4000, 0, "J"), [8, 10, 12, 14], snr=30, wavelength=1282.2)
    print(f"N={n:5d}: " + "  ".join(f"J={m}: {ti:7.0f} s" for m, ti in zip((8, 10, 12, 14), t)))"""),
md("""
## Radial velocities at 1.55 um, R = 43,000

Photon-limited RV precision scales as $\\sigma_{RV} = c / (Q \\sqrt{N_{e}})$ (Bouchy et al. 2001)
with $N_e$ the total electrons in the spectrum and $Q$ the spectral quality factor. For a K or M
dwarf at R = 43k in a 40 A window in the H band (one CSHELL setting) we take $Q$ = 300. Bouchy
et al. give Q ~ 10$^4$ at R = 100k across the full optical range, and a 40 A window with few
lines is far poorer, so $Q$ is the knob to revisit.

This is a photon-only bound with an assumed $Q$: $N_e$ counts source electrons only (no sky,
dark or read noise enters), and a blackbody carries no line information, so the number is only
as good as $Q$. The result is a scaling: km/s orbits for binaries are easy, m/s planet work is not.
"""),
code("""Q = 300.0
RV_WINDOW_A = 40; RV_START_NM = 1550.0                      # one CSHELL setting, 40 A wide
def rv(n):
    return LFASTETC("cshell_r43000", telescope=tel(n), wavelength_min_nm=RV_START_NM,
                    wavelength_max_nm=RV_START_NM + RV_WINDOW_A/10)
for n in (20, 200, 1720):
    e = rv(n); row = []
    for h in (6, 9, 12):
        src = Source.blackbody(4000, h, "H")
        r = e.snr(src, 900, n_exp=e.reads_needed(src, 900))
        ne = np.nansum(r.n_star*r.exptime)        # source electrons across the 40 A window
        row.append(f"H={h}: {2.998e5/(Q*np.sqrt(ne)):7.2f} km/s")
    print(f"N={n:5d} (15 min, photon-only, Q assumed): " + "  ".join(row))"""),
md("## Bright end: reads needed for a J = 3 star in 300 s (55 ke- full well)"),
code("""for n in (20, 200, 1720):
    print(f"N={n:5d}: reads = {he(n, 'cshell_r18000').reads_needed(Source.blackbody(4800, 3, 'J'), 300)}")"""),
]


if __name__ == "__main__":
    names = sys.argv[1:] or list(NB)
    for n in names:
        build(n, NB[n])
