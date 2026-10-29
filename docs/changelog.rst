Changelog
==

Unreleased
----------

Code review, notebooks and documentation pass (see the GitHub issues #2 to #30 for the full list).

Fixed

* ``K_opt`` samples below ``KOPT_EPS = 1e-6`` are telluric cores and count as zero, together
  with the interpolation intervals next to them, so the 687 nm O2 B band on the
  ``bench_316_7p0_o1`` preset is NaN instead of a 5x finite exposure-time spike (#12).
* ``ETCResult.valid`` marks the usable bins. ``exptime()`` no longer broadcasts the requested
  SNR into invalid bins; ``at()`` returns the sample at an exact bin centre (also on a one-bin
  grid) and never interpolates across or into an invalid bin (#37).
* An explicit ``wavelength_min_nm``/``wavelength_max_nm`` is now clipped to the preset bandpass
  as well as to the data range, with a clear error when nothing is left (#2).
* A float ``wavelength`` selector outside the grid raises instead of snapping to the edge bin (#3).
* Web: the zero-throughput fraction follows the solved quantity (#4); summary percentiles are
  ``null`` rather than bare ``NaN`` for all-NaN ranges (#5).
* ``n_telescopes=0`` raises instead of falling back to 20 (#6).
* ``ETCResult.at`` interpolates across NaN bins instead of into them (#7).
* Finite-value checks on magnitudes, exposure times, line parameters, ``reads_needed`` margin;
  efficiencies bounded by 1; clearer web messages for ``inf`` (#8).
* ``Instrument.from_preset(preset, **overrides)`` can override ``name``; ``from_arrays`` with a
  bad band raises ``ValueError`` (#9).
* Science notebooks 06-10 no longer require the gitignored ``memo/figures`` directory (#10).
* ``Telescope`` warns for the 50 um fiber through the combiner (#11).
* The peak-pixel estimate follows the most heavily loaded combiner input (a full 20-group for
  ``N = 21``) instead of averaging over inputs (#36); ``Telescope.telescopes_on_brightest_input``.
* ``reads_needed`` rejects ``margin > 1`` (#39).
* ``Instrument`` accepts the default ``full_well = inf`` (or ``None``) as "unknown";
  ``summary()["saturated_fraction"]`` is then ``None`` (#46).
* Overflowed ``snr``/``exptime`` raise ``ValueError`` instead of returning ``inf``; the web API
  bounds ``snr`` (1e6) and ``exptime`` (1e9 s) and serialises every non-finite value as ``null``
  (#47).
* Notebook 07 integrates the stated 50 A AGN window (507.5 to 512.5 nm) and names the estimator
  (#35); notebook 10 integrates one 40 A CSHELL setting (1550 to 1554 nm) and labels the RV
  number as a photon-only, assumed-Q bound (#34).

Added

* ``.readthedocs.yaml``: the Sphinx site builds on Read the Docs with the ``docs`` and ``web``
  extras and pandoc (#26).
* Plot helpers: ``ETCResult.plot_noise_budget``, ``plot_spectrum`` and ``simulate``,
  ``to_table``; ``LFASTETC.plot_throughput``, ``plot_sky``; ``Source.plot``;
  ``plotting.plot_vs_n_telescopes``, ``equivalent_diameter_axis``, ``sky_surface_brightness_ab``.
* The ten notebooks rewritten as a walkthrough for new users, with checks and exercises.
* This documentation: concepts, instruments (preset table generated from the TOML files),
  telescope, sources, sky, sensitivity, plot gallery, web API, contributing.
* GitHub Actions: tests on Python 3.11-3.13 and a documentation build.

Changed

* CSHELL detector from Brown et al. (2026, SPIE 14149, Table 1; the IRTF measurements of the
  Hughes 256x256 SBRC InSb): dark current 0.5 e-/s/pix (the IRTF page's 5 e-/s was the K-band
  background through the IRTF), read noise 22 e- at 6 Fowler samples, 55,000 e- 1%-linear well at the
  default bias, QE 0.8. CSHELL faint limits recover about 1 mag relative to the 5 e-/s assumption.
* The extraction aperture is the full fiber image (#13): ``n_pix = pixels_per_resel *
  fiber_image_px``, with a circular image by default and the native 6 px height set explicitly
  for the CSHELL slit modes. Sky, dark and read-noise terms grow accordingly (x2.2 on the Bench);
  ``fiber_image_px=1.0`` restores the reference notebook's dispersion-only count.
* The peak pixel uses a rendered fiber-image profile (#19): ``Instrument.profile`` is
  ``"tophat"`` (default), ``"uniform"`` or ``"gaussian"``; ``peak_fraction``, ``pixel_image()``
  and ``plot_profile()`` expose it. Saturation flags and ``reads_needed`` follow.
* ``Telescope(seeing_arcsec=...)`` scales the point-source throughput with a Moffat (or Gaussian)
  seeing disc integrated over the fiber, relative to a placeholder 1.0 arcsec reference (#14);
  ``fiber_coupling()`` and ``plot_fiber_coupling()`` expose it. The sky path is unchanged.

* Full wells are sourced instead of placeholders: 87,500 e- for every Bench mode (WIYN STA1
  CCD page, 1% linearity to 63,000 e-) and, for CSHELL, 55,000 e- (Brown et al. 2026; superseding the 40,000 e- inferred from the IRTF page). Saturation flags, ``reads_needed`` and the web markers follow.
* Documentation theme is ``pydata_sphinx_theme`` (as for the WCC ETC).
* The duplicated wavelength-selector and box-band code paths are merged.

0.1.0 (2026-10-21)
------------------

* First release: per-resel ETC reproducing the reference LFAST ETC notebook (2026-10-12) to 1e-6;
  WIYN Bench reference mode plus 14 modes from the WIYN configuration table; three CSHELL modes;
  18 um / 50 um fibers with the combiner model and ``telescopes_per_input``; sensitivity maps,
  limiting magnitudes, exposure versus magnitude, ``reads_needed``; Flask/Bokeh browser ETC with
  a JSON API; notebooks 01-10; independent review (``review_report``).
