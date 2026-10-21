# lfast-etc

Exposure Time Calculator for the **LFAST** telescope array at Kitt Peak (nominally 20 x 0.726 m
unit telescopes feeding one spectrograph through an 18 um fiber combiner), with two instruments:
the **WIYN Bench** spectrograph (visible) and **CSHELL** (1.08-1.75 um, R up to 43,000).

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

etc = LFASTETC("bench_r1800", n_telescopes=20)          # presets: bench_r1800, cshell_r18000/21500/43000
star = Source.blackbody(5800, mag=15, band="V")
res = etc.exptime(star, snr=50)                         # per resolution element
res.plot(); res.summary(); res.at(550.0)

res = etc.snr(star, exptime=600, n_exp=3)               # SNR for 600 s in 3 reads
etc.limiting_magnitude(star, exptime=900, snr=10, wavelength=(500, 650))
etc.snr_grid(star, mags, exptimes, wavelength=(500, 650))   # sensitivity map

# scale the array
LFASTETC("cshell_r43000", telescope=Telescope(n_telescopes=40, fiber="18um"),
         wavelength_min_nm=1081, wavelength_max_nm=1086)
```

## Browser ETC

```sh
pip install -e ".[web]"
lfast-etc-web --port 5000     # http://127.0.0.1:5000  (JSON at /api/calc?...)
```

## Notebooks

`notebooks/` holds executed examples: quick start, brightness sensitivity maps and array
scaling, stellar metallicity surveys, TDAMM follow-up, and exoplanet atmospheres (He 10830)
with CSHELL. Rebuild them with `python scripts/build_notebooks.py`.

## Docs

```sh
pip install -e ".[docs]"
cd docs && make html         # needs pandoc on PATH for the notebooks
```

## Tests

```sh
pytest
```

## Layout

```
src/lfast_etc/        package (etc.py, source.py, telescope.py, instrument.py, sky.py, plotting.py)
src/lfast_etc/data/   throughput curves, SkyCalc table, instrument presets (TOML)
src/lfast_etc/web/    Flask + Bokeh browser ETC
notebooks/            executed example notebooks
docs/                 Sphinx site
input/                original hackathon materials (reference notebook; the data files are the
                      same bytes as src/lfast_etc/data, kept here for provenance)
tests/                pytest, including the regression test against the reference notebook
```
