![Logo](docs/logo_schmidt_observatory_system.png)

# lfast-etc

Exposure Time Calculator for the **LFAST** telescope array at Kitt Peak (nominally 20 x 0.726 m
unit telescopes feeding one spectrograph through an 18 um fiber combiner), with two instruments:
the **WIYN Bench** spectrograph (visible; the 15 grating/order modes of the
[WIYN configuration table](https://www.wiyn.org/Instruments/wiynbench.html) that reach 400 nm, scaled from the
reference 600 @ 10.1 deg mode by `scripts/build_bench_presets.py`) and **CSHELL** (1.08-1.75 um, R up to 43,000).

The number of telescopes, the fiber, and every instrument parameter are plain inputs, so the
calculator follows the array design as it evolves. The physics is the per-resolution-element model
of the reference LFAST ETC notebook (2026-10-12) with the Choi et al. (2026, SPIE 14147)
throughput model; the test-suite reproduces that notebook to 1e-6.

## Install

Requires Python >= 3.11 and numpy >= 2.0.

```sh
git clone git@github.com:schmidt-observatories/lfast-etc.git
cd lfast-etc
pip install -e .            # core
pip install -e ".[web]"     # + browser ETC
pip install -e ".[docs]"    # + Sphinx docs
```

## Quick start

```python
from lfast_etc import LFASTETC, Source, Telescope

etc = LFASTETC("bench_r1800", n_telescopes=20)          # presets: list_instruments()
star = Source.blackbody(5800, mag=15, band="V")
res = etc.exptime(star, snr=50)                         # per resolution element
res.plot(); res.summary(); res.at(550.0); res.to_table()

res = etc.snr(star, exptime=600, n_exp=3)               # SNR for 600 s in 3 reads
res.plot_noise_budget(); res.plot_spectrum()            # what limits you; what the spectrum looks like
etc.plot_throughput(); etc.plot_sky(); star.plot()
etc.limiting_magnitude(star, exptime=900, snr=10, wavelength=(500, 650))
etc.snr_grid(star, mags, exptimes, wavelength=(500, 650))   # sensitivity map

# scale the array
LFASTETC("cshell_r43000", telescope=Telescope(n_telescopes=40, fiber="18um"),
         wavelength_min_nm=1081, wavelength_max_nm=1086)
```

Everything is per resolution element; SNR per pixel is lower by about `sqrt(pixels_per_resel)`.
The bundled sky is grey time (moon up), so faint-end numbers are 0.3-0.5 mag pessimistic for
dark time. The docs (`docs/`, build with `make html`) have a quick start, the concepts, every
preset, the physics and its known limitations, and a plot gallery.

## Browser ETC

```sh
pip install -e ".[web]"
lfast-etc-web --port 5000     # http://127.0.0.1:5000  (JSON at /api/calc?...)
```

## Notebooks

`notebooks/` holds executed examples, written as a walkthrough for someone new to the package.
Start with 01 (the whole API, with checks and exercises at the end). 01-05 are rebuilt with
`python scripts/build_notebooks.py`: quick start, sensitivity maps and array scaling, stellar
metallicity, transient classification, exoplanet atmospheres (He 10830) with CSHELL. 06-10
(`python scripts/build_science_notebooks.py`) compare science cases versus the array size
(20 telescopes, 10 m equivalent = 200, 30 m equivalent = 1720): array-scaling overview,
LSST-era time domain, Galactic archaeology, stellar activity and flares, CSHELL NIR science.
The notebooks are generated from those scripts; edit the scripts, not the `.ipynb` files.

## Docs

```sh
pip install -e ".[docs]"
cd docs && make html         # needs pandoc on PATH for the notebooks
```

## Tests

```sh
pytest
```

CI (GitHub Actions) runs the tests on Python 3.11-3.13 and builds the docs on every push and
pull request. Bugs and planned features are tracked in the
[issues](https://github.com/schmidt-observatories/lfast-etc/issues), labelled by area and priority.

## Layout

```
src/lfast_etc/        package (etc.py, source.py, telescope.py, instrument.py, sky.py, plotting.py)
src/lfast_etc/data/   throughput curves, SkyCalc table, instrument presets (TOML)
src/lfast_etc/web/    Flask + Bokeh browser ETC
notebooks/            executed example notebooks
docs/                 Sphinx site
tests/                pytest, including the regression test against the reference notebook
```

## License

BSD 3-Clause. See [LICENSE](LICENSE).
