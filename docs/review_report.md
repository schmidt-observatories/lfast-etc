# Independent review report

Three independent reviewer agents examined the package before release, each with a fresh
context and read-only access: (1) physics and numerics against the reference LFAST ETC notebook
(2026-10-12), (2) code quality, API, packaging and robustness, (3) science content of the
notebooks and documentation. Every item below was verified by the reviewer with a reproduction;
the status column records what was done afterwards. All 26 tests pass after the fixes.

## Bugs found and fixed

| # | Area | Finding | Fix |
|---|------|---------|-----|
| 1 | physics | Peak-pixel / saturation estimate was `N_telescopes` times too high in `combine="independent"` mode (the 50 um default), so every 50 um run with 20 telescopes was flagged saturated 20x too early. | Peak pixel is now computed per physical detector; regression test added. |
| 2 | physics | Zero-throughput bins (`K_opt = 0` at the 687 nm O2 B band and the 1355-1420 nm water band) returned exposure times of 1e20-1e42 s and limiting magnitudes of -750 instead of NaN, corrupting window medians, percentiles and the CSHELL limiting-magnitude plot. | NaN in `snr`/`exptime`/`limiting_magnitude`; `nanmedian` everywhere; Bench preset ends at 686 nm (the 688 nm edge sat on the 687 nm zero). |
| 3 | physics | Narrow spectral lines were point-sampled: a 0.05 nm H-alpha line on a 0.36 nm resel appeared as a 2.6% dip instead of 30% and changed with grid alignment (inherited from the reference notebook). | Lines are averaged over each resel (equivalent width preserved, `Source.line_factor`); test added. The regression test now uses a line-free source. |
| 4 | code | `requires-python = ">=3.10"` but `tomllib` needs 3.11; `numpy>=1.24` but `np.trapezoid` needs 2.0. The package could not import on either floor. | Floors raised to Python >= 3.11, numpy >= 2.0; README states them. |
| 5 | code | `/api/calc` emitted bare `NaN` tokens (invalid JSON) once telluric bins became NaN. | NaN -> `null` at the API boundary; summary guards; test added. |
| 6 | code | An unbounded `resolving_power` override in the web form (e.g. 1e9) allocated ~5e8 bins and hung the worker. | `LFASTETC` refuses grids above 5e6 bins; the web form caps every override. |
| 7 | science | NB05 used HD 189733's H magnitude (5.6) as J (2MASS J = 6.07); HAT-P-11 b depth 1.2% vs 1.08% (Allart et al. 2018). | Target table corrected with references. |
| 8 | science | NB05 compared R=18k and R=43k with the full literature line depth in both, ignoring that a 0.7 A line averaged over 1.2 A (R=18k, 2 resels) retains only 59% of its depth versus 79% at R=43k; R=18k was over-credited by ~1.3x. | `core_fraction()` applies the resel-averaged profile per mode. |
| 9 | science | NB04 fed AT2017gfo r-band magnitudes to a V-normalised blackbody; for the 3000 K late SED V-R ~ 1 mag, so days 3-7 were ~1 mag too bright. Fading rate (1 mag/day) and cooling ramp were also off. | NB04 normalises in R (Cousins ~ SDSS r), uses 0.6 mag/day and the Drout et al. 2017 temperature ramp. |
| 10 | science | NB04 claimed the sqrt(k) rebinning gain is "slightly optimistic in the read-noise limit"; in this noise model every variance term scales with the number of pixels so sqrt(k) is exact. | Text corrected. |

## Issues fixed

* `min_exptime` was loaded but unused: results now carry `below_min_exptime` (reads shorter than the Fowler minimum, where the quoted read noise does not apply); shown in the web summary.
* `wavelength="all"` only worked in `limiting_magnitude`; now accepted by all sensitivity functions.
* `ETCResult.at()` silently clamped outside the grid; now NaN.
* `summary()` leaked numpy scalars; now plain floats with int keys.
* `ETCResult` and `Source` dataclass equality crashed or compared lambdas by identity; `eq=False`.
* `Telescope` cache fields were exposed in `__init__`; now `init=False`.
* Silent garbage on `n_exp=0`/`2.5`, negative times/SNRs, fractional telescopes, negative diameters, zero resolving power, unknown magnitude systems and flux units; all now raise `ValueError` at the public boundary. `wavelength_min_nm=0` is rejected rather than treated as unset. Passing both `telescope=` and `n_telescopes=` raises.
* Throughput and sky rate caches went stale if attributes were mutated; everything is precomputed in `__init__` and the class is documented as immutable.
* `snr_grid` broadcast a (mag, time, bin) cube before selecting bins (1.3 GB for 50x50 on full-range CSHELL); bins are selected first.
* `ETCResult.plot` title overran the figure and parsed a human string; a short `title` is stored on the result, and y labels follow the plotted field.
* Web: fractional integers were truncated silently; preset placeholders are computed once; limiting magnitude of n/a is handled.
* `Instrument.replace` reimplemented `dataclasses.replace`; `from_preset` validates against the preset list (no path traversal).
* The `docs` extra listed the `pandoc` PyPI wrapper (not the binary); removed. Deprecated license metadata updated.
* `tests/` imports the reference oracle via `pythonpath = ["tests"]` so `--import-mode=importlib` works.
* NB02: the CSHELL sensitivity map starts at the 8 s read minimum; the non-existent "50 um through the combiner" configuration was dropped; the CSHELL limiting-magnitude caption now attributes the structure to `K_opt` tellurics, and states that each wavelength is a separate echelle setting.
* NB03: LAMOST/SEGUE anchors corrected (per-pixel SNR vs per-resel, ~0.1 dex floor for LASP/SSPP, ~0.05 dex only for data-driven methods); M-dwarf row caveated.
* NB04: classification SNR anchored to SEDM / 4MOST-TiDES; SED colour now matters (R normalisation); knobs gathered in the first code cell.
* NB05: missing caveats added (continuum SNR vs stellar He line core, OH doublet and tellurics, real transit durations); CCF numbers labelled as guesses; knobs gathered; `LFASTETC.reads_needed` replaces the hand-rolled saturation split.
* Docs: per-resel vs per-pixel paragraph, SkyCalc instructions and required columns, preset field list, grey-time sky statement, additional known limitations (no LSF, telluric zeros, box-band approximation, sky-subtraction term omitted, `n_pix` counts the dispersion direction only).

## Open items for the LFAST team (not fixable from the data in hand)

1. **`n_pix` is the dispersion direction only** (2.2 Bench, 2.5 CSHELL, from the reference notebook). Howell's `n_pix` is the full extraction aperture; multiply by the fiber-image height in pixels at 1x1 binning once known. For the Bench at 3.7 e- read noise this is the dominant uncertainty for faint targets.
2. **CSHELL slit throughputs** for the 200 um (R=21,500) and 100 um (R=43,000) slits are placeholders (0.75, 0.5). The reference notebook has "XX%"/"YY%".
3. **Detector full wells** (100 ke- Bench, 80 ke- CSHELL) and the flat QE values are placeholders.
4. **Sky conditions.** The bundled SkyCalc table is grey time (moon up, ~60% of the visible sky), ESO model at 2640 m (SkyCalc has no Kitt Peak site), sampled at R=20,000, so airglow lines are pre-smoothed for CSHELL. A dark-time table and an R >= 50,000 table should be generated for the CSHELL presets.
5. **Is the 50 um fiber available on all 20 telescopes?** The package assumes yes with independent spectra; the reference notebook suggests testing it with `N_TELESCOPES = 1`.
6. **Reference-notebook text vs code.** Its markdown says the star path uses `K_opt * eta_fib`; the code (and this package) use `K_opt` alone, which the throughput files prove correct (`K_opt_50/K_opt_18 == eta_50/eta_18`, i.e. `K_opt` already contains fiber injection). The notebook text should be corrected so nobody "fixes" the package later.
7. **Telluric treatment in `K_opt`** belongs to whatever airmass Choi et al. assumed, not necessarily the sky table's 1.3; `apply_extinction=True` double counts.

## Suggested future science notebooks (from the science reviewer)

Time-resolved flare spectroscopy of nearby M dwarfs with the Bench (1-10 s cadence, no per-telescope
readout penalty); Paschen-beta accretion monitoring of T Tauri stars with CSHELL (1282 nm fits in
one setting); radial velocities of Gaia astrometric binaries at R=43k (needs an RV-precision
scaling); spectrophotometric AGN reverberation monitoring (needs a fiber-stability budget).
