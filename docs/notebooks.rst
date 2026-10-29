Notebooks
=========

Ten executed notebooks live in the repository's ``notebooks/`` directory and are rendered here.
Start with 01, which walks through the whole API and ends with checks and exercises; 02 adds the
sensitivity maps and the array scaling; 03 to 05 are single-instrument science cases; 06 to 10
compare the nominal 20-telescope array with 10 m and 30 m equivalents.

Each science notebook states its assumptions in the first text cell and gathers the adjustable
numbers in the first code cell, so changing a threshold or a target is one edit. The notebooks
are generated from ``scripts/build_notebooks.py`` and ``scripts/build_science_notebooks.py``;
see :doc:`contributing` for how to rebuild them.

Walkthrough
-----------

- :doc:`notebooks/01_quickstart`: Build a calculator, define a source, exposure time and SNR, noise budget, simulated spectrum, sky, overrides, array scaling, your own spectrum, checks and exercises.
- :doc:`notebooks/02_sensitivity_maps_and_array_scaling`: SNR over (magnitude, exposure), limiting magnitude versus wavelength and versus the number of telescopes, dark versus grey sky, exposure versus magnitude.

.. toctree::
   :maxdepth: 1
   :hidden:

   notebooks/01_quickstart
   notebooks/02_sensitivity_maps_and_array_scaling

Science cases, nominal array
----------------------------

- :doc:`notebooks/03_stellar_metallicity`: LAMOST-like R~1800 metallicity depth and survey speed with the Bench.
- :doc:`notebooks/04_tdamm_followup`: Classification depth at R = 300, kilonova time-to-classify and afterglow time-to-redshift.
- :doc:`notebooks/05_exoplanet_atmospheres_cshell`: He 10830 transmission spectroscopy at R = 43,000 and 18,000, and the FSR limits of cross-correlation.

.. toctree::
   :maxdepth: 1
   :hidden:

   notebooks/03_stellar_metallicity
   notebooks/04_tdamm_followup
   notebooks/05_exoplanet_atmospheres_cshell

Science cases versus the array size
-----------------------------------

- :doc:`notebooks/06_array_scaling_overview`: 20 telescopes versus 10 m and 30 m equivalent: the combiner model, limiting magnitudes, bright limits.
- :doc:`notebooks/07_time_domain_lsst_era`: Classification depth, kilonova horizon, SN Ia redshift reach, TDE and AGN monitoring.
- :doc:`notebooks/08_galactic_archaeology`: Single-target metallicity depth, distance reach by tracer, targets per night.
- :doc:`notebooks/09_stellar_activity_flares`: 10 s cadence H-alpha spectroscopy of M dwarfs, flare spectra, host-star monitoring.
- :doc:`notebooks/10_cshell_nir_science`: He 10830 reach, Paschen-beta accretion monitoring, km/s radial velocities.

.. toctree::
   :maxdepth: 1
   :hidden:

   notebooks/06_array_scaling_overview
   notebooks/07_time_domain_lsst_era
   notebooks/08_galactic_archaeology
   notebooks/09_stellar_activity_flares
   notebooks/10_cshell_nir_science
