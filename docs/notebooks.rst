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

.. grid:: 1 2 2 2
   :gutter: 2

   .. grid-item-card:: 01 Quick start
      :link: notebooks/01_quickstart
      :link-type: doc

      Build a calculator, define a source, exposure time and SNR, noise budget, simulated spectrum, sky, overrides, array scaling, your own spectrum, checks and exercises.

   .. grid-item-card:: 02 Sensitivity maps and array scaling
      :link: notebooks/02_sensitivity_maps_and_array_scaling
      :link-type: doc

      SNR over (magnitude, exposure), limiting magnitude versus wavelength and versus the number of telescopes, dark versus grey sky, exposure versus magnitude.

Science cases, nominal array
----------------------------

.. grid:: 1 2 2 2
   :gutter: 2

   .. grid-item-card:: 03 Stellar metallicity
      :link: notebooks/03_stellar_metallicity
      :link-type: doc

      LAMOST-like R~1800 metallicity depth and survey speed with the Bench.

   .. grid-item-card:: 04 Transient classification
      :link: notebooks/04_tdamm_followup
      :link-type: doc

      Classification depth at R = 300, kilonova and afterglow time-to-classify.

   .. grid-item-card:: 05 Exoplanet atmospheres with CSHELL
      :link: notebooks/05_exoplanet_atmospheres_cshell
      :link-type: doc

      He 10830 transmission spectroscopy at R = 43,000 and 18,000, and the FSR limits of cross-correlation.

Science cases versus the array size
-----------------------------------

.. grid:: 1 2 2 2
   :gutter: 2

   .. grid-item-card:: 06 Array scaling overview
      :link: notebooks/06_array_scaling_overview
      :link-type: doc

      20 telescopes versus 10 m and 30 m equivalent: the combiner model, limiting magnitudes, bright limits.

   .. grid-item-card:: 07 Time domain in the LSST era
      :link: notebooks/07_time_domain_lsst_era
      :link-type: doc

      Classification depth, kilonova horizon, SN Ia redshift reach, TDE and AGN monitoring.

   .. grid-item-card:: 08 Galactic archaeology
      :link: notebooks/08_galactic_archaeology
      :link-type: doc

      Single-target metallicity depth, distance reach by tracer, targets per night.

   .. grid-item-card:: 09 Stellar activity and flares
      :link: notebooks/09_stellar_activity_flares
      :link-type: doc

      10 s cadence H-alpha spectroscopy of M dwarfs, flare spectra, host-star monitoring.

   .. grid-item-card:: 10 CSHELL NIR science
      :link: notebooks/10_cshell_nir_science
      :link-type: doc

      He 10830 reach, Paschen-beta accretion monitoring, km/s radial velocities.

.. toctree::
   :hidden:

   notebooks/01_quickstart
   notebooks/02_sensitivity_maps_and_array_scaling
   notebooks/03_stellar_metallicity
   notebooks/04_tdamm_followup
   notebooks/05_exoplanet_atmospheres_cshell
   notebooks/06_array_scaling_overview
   notebooks/07_time_domain_lsst_era
   notebooks/08_galactic_archaeology
   notebooks/09_stellar_activity_flares
   notebooks/10_cshell_nir_science
