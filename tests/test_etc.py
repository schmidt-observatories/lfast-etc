import numpy as np
import pytest

from lfast_etc import LFASTETC, Source, Telescope, Instrument, list_instruments, Sky
from reference_notebook import run as reference_run


def test_regression_against_reference_notebook():
    """Bench, 400-1800 nm resel grid, 5800 K V=15 star with an H-alpha line, SNR=50."""
    wl_ref, t_ref, nstar_ref, nsky_ref = reference_run(spectral_lines=())
    # the reference evaluates the Bench over the whole throughput range, so widen the preset bandpass
    wide = Instrument.from_preset("bench_r1800", wavelength_min_nm=400, wavelength_max_nm=1800,
                                  fiber_image_px=1.0)      # the notebook counts the dispersion direction only
    etc = LFASTETC(wide, n_telescopes=20, lsf=None)     # the oracle is the boxcar per resel
    src = Source.blackbody(5800, mag=15, band="V")
    res = etc.exptime(src, snr=50)
    # the reference grid overshoots lambda_max by one step; compare on the common bins
    n = min(len(wl_ref), len(res.wavelength))
    np.testing.assert_allclose(res.wavelength[:n], wl_ref[:n], rtol=1e-9)
    ok = res.n_star[:n] > 0          # the oracle clamps zero-throughput bins to 1e-20; we return NaN
    assert np.isnan(res.exptime[:n][~ok]).all() and (~ok).sum() > 0
    np.testing.assert_allclose(res.n_star[:n][ok], nstar_ref[:n][ok], rtol=1e-6)
    np.testing.assert_allclose(res.n_sky_pix[:n][ok], nsky_ref[:n][ok], rtol=1e-6)  # cumsum rounding
    np.testing.assert_allclose(res.exptime[:n][ok], t_ref[:n][ok], rtol=1e-6)


def test_snr_exptime_roundtrip():
    etc = LFASTETC("cshell_r18000", wavelength_min_nm=1500, wavelength_max_nm=1560)
    src = Source.blackbody(4000, mag=12, band="J")
    t = etc.exptime(src, snr=30, n_exp=3)
    s = etc.snr(src, exptime=t.exptime, n_exp=3)
    np.testing.assert_allclose(s.snr, 30, rtol=1e-8)


def test_fiber_combine_scales_signal_linearly_and_independent_as_sqrt():
    src = Source.blackbody(5800, mag=20, band="V")
    one = LFASTETC("bench_r1800", telescope=Telescope(1, combine="fiber"))
    t1 = one.rates(src)
    # fiber combine: both star and sky scale by N, detector terms don't
    twenty = LFASTETC("bench_r1800", telescope=Telescope(20, combine="fiber"))
    t20 = twenty.rates(src)
    np.testing.assert_allclose(t20["n_star"], 20 * t1["n_star"])
    np.testing.assert_allclose(t20["n_sky_pix"], 20 * t1["n_sky_pix"])
    assert t20["read_noise"] == t1["read_noise"]
    # independent: SNR = sqrt(N) * SNR_1 exactly
    ind = LFASTETC("bench_r1800", telescope=Telescope(20, combine="independent"))
    s1 = one.snr(src, 100).snr
    s20 = ind.snr(src, 100).snr
    np.testing.assert_allclose(s20, np.sqrt(20) * s1, rtol=1e-10)


def test_multi_array_inputs():
    src = Source.blackbody(5800, mag=20, band="V")
    one = LFASTETC("bench_r1800", telescope=Telescope(20))
    ten = LFASTETC("bench_r1800", telescope=Telescope(200))
    assert ten.telescope.n_inputs == 10 and Telescope(201).n_inputs == 11
    np.testing.assert_allclose(ten.snr(src, 100).snr, np.sqrt(10) * one.snr(src, 100).snr, rtol=1e-10)
    ideal = LFASTETC("bench_r1800", telescope=Telescope(200, telescopes_per_input=None))
    assert np.nanmedian(ideal.snr(src, 100).snr) > np.nanmedian(ten.snr(src, 100).snr)
    assert abs(Telescope(190).equivalent_diameter_m - 10.0) < 0.02


def test_limiting_magnitude_consistent_with_snr():
    etc = LFASTETC("bench_r1800", wavelength_min_nm=500, wavelength_max_nm=600)
    src = Source.blackbody(5800, mag=15, band="V")
    m = etc.limiting_magnitude(src, exptime=600, snr=10, wavelength="all")
    s = etc.snr(src.with_mag(0), 600)  # only to get wavelength grid shape
    assert m.shape == s.snr.shape
    # at the per-bin limiting magnitude, SNR must equal 10 - check one bin end to end
    i = len(m) // 2
    etc_single = LFASTETC("bench_r1800", wavelength_min_nm=500, wavelength_max_nm=600)
    r = etc_single.snr(src.with_mag(m[i]), 600)
    assert abs(r.snr[i] - 10) < 1e-6


def test_snr_grid_matches_direct():
    etc = LFASTETC("bench_r1800")
    src = Source.blackbody(5800, mag=15, band="V")
    grid = etc.snr_grid(src, mags=[12, 18], exptimes=[10, 1000], wavelength=550.0)
    direct = etc.snr(src.with_mag(18), 1000).at(550.0)["snr"]
    assert grid.shape == (2, 2)
    assert abs(grid[1, 1] - direct) / direct < 1e-3   # .at() interpolates, grid uses nearest bin


def test_source_normalisation_and_bands():
    src = Source.blackbody(5800, mag=10, band="V")
    assert abs(src.magnitude("V") - 10) < 1e-9
    # a 5800 K blackbody is bluer in V-J than a 3500 K one
    assert (Source.blackbody(5800, 10).magnitude("J") >
            Source.blackbody(3500, 10).magnitude("J"))
    flat = Source.flat(mag=20, band="V", system="AB")
    assert abs(flat.magnitude("R", "AB") - 20) < 0.01   # box-band F_nu->F_lambda approx
    wl = np.linspace(300, 2000, 500)
    user = Source.from_arrays(wl, np.ones_like(wl), mag=15, band="J")
    assert abs(user.magnitude("J") - 15) < 1e-9


def test_line_equivalent_width_preserved_on_coarse_grid():
    """A 0.05 nm line on a 0.36 nm resel grid keeps its equivalent width (reference notebook point-sampled it)."""
    src = Source.blackbody(5800, 15).add_line(656.3, 0.3, 0.05)
    etc = LFASTETC("bench_r1800")
    ew = np.sum((1 - src.line_factor(etc.wavelength, etc.delta_lambda)) * etc.delta_lambda)
    assert abs(ew - 0.3 * 0.05 * np.sqrt(2 * np.pi)) < 1e-4   # geometric grid: bins overlap slightly
    # fine sampling and bin averaging agree when bins are much narrower than the line
    fine = np.linspace(655, 657.6, 20001)
    np.testing.assert_allclose(src.line_factor(fine, fine[1] - fine[0]), src.line_factor(fine), atol=1e-6)
    assert etc.reads_needed(Source.blackbody(5800, 15), 600) == 1
    assert etc.reads_needed(Source.blackbody(5800, 2), 600) > 1
    with pytest.raises(ValueError):
        src.add_line(500, 1.5, 0.1)


def test_source_errors():
    with pytest.raises(ValueError):
        Source.blackbody(5800, 10, band="Q")
    with pytest.raises(ValueError):
        Source.from_arrays([400, 500, 600], [1, 1, 1], mag=10, band="J")   # no flux in J


def test_instrument_presets_load_and_override():
    for name in list_instruments():
        ins = Instrument.from_preset(name)
        assert ins.resolving_power > 0
    ins = Instrument.from_preset("bench_r1800", resolving_power=5000)
    assert ins.resolving_power == 5000
    assert ins.replace(read_noise=1).read_noise == 1
    with pytest.raises(ValueError):
        Instrument.from_preset("../config/bench_r1800")
    with pytest.raises(ValueError):
        Instrument.from_preset("bench_r1800", resolving_power=-5)


def test_setting_window():
    lo, hi = Instrument.from_preset("cshell_r18000").setting_window(1555.0)
    assert hi - lo == pytest.approx(1555.0 * 2.5e-3)
    assert (lo + hi) / 2 == pytest.approx(1555.0)
    assert Instrument.from_preset("cshell_r18000").setting_window(1080.5)[0] == 1080.0   # clipped
    bench = Instrument.from_preset("bench_r1800")
    assert bench.setting_window(550.0) == (bench.wavelength_min_nm, bench.wavelength_max_nm)
    with pytest.raises(ValueError):
        bench.setting_window(1000.0)


@pytest.mark.parametrize("name", list_instruments())
def test_every_preset_runs_inside_the_throughput_range(name):
    ins = Instrument.from_preset(name)
    lo, hi = Telescope().wavelength_range_nm
    assert lo <= ins.wavelength_min_nm < ins.wavelength_max_nm <= hi
    assert ins.description
    etc = LFASTETC(name)
    assert etc.wavelength.min() >= ins.wavelength_min_nm and etc.wavelength.max() <= ins.wavelength_max_nm
    res = etc.snr(Source.blackbody(5800, mag=12, band="V"), exptime=600)
    assert np.isfinite(np.nanmedian(res.snr)) and np.nanmedian(res.snr) > 0


def test_cshell_presets_follow_the_irtf_detector_and_dispersion():
    """Detector numbers from Brown et al. 2026 (IRTF measurements); pixels per resel scale as 2.5 x 43,000 / R
    (2.7 km/s per pixel), slit losses and efficiency stay placeholders."""
    for name, npix in (("cshell_r18000", 6.0), ("cshell_r21500", 5.0), ("cshell_r43000", 2.5)):
        ins = Instrument.from_preset(name)
        assert ins.pixels_per_resel == npix
        assert ins.pixels_per_resel == pytest.approx(2.5 * 43000 / ins.resolving_power, rel=0.01)
        assert (ins.read_noise, ins.dark_current, ins.min_exptime) == (22.0, 0.5, 8.0)   # Brown et al. 2026, Table 1
        assert ins.full_well == 55000.0        # 1% linear well at the default bias (Brown et al. 2026, Table 1)
        assert ins.free_spectral_range_frac == 2.5e-3 and ins.spectrograph_efficiency == 0.1


def test_bench_modes_follow_the_wiyn_table():
    """R and efficiency scale from the bench_r1800 anchor (600 @ 10.1, 5500 A / 3.35 A, 72%)."""
    bench = [n for n in list_instruments() if n.startswith("bench_")]
    assert len(bench) == 15 and "bench_316_63p4_o14" not in bench   # that mode is < 400 nm
    m = Instrument.from_preset("bench_316_63p4_o9")        # 6563 A / 0.23 A, 50%
    assert abs(m.resolving_power - 1800 * (6563 / 0.23) / (5500 / 3.35)) < 1
    assert abs(m.spectrograph_efficiency - 0.3 * 50 / 72) < 1e-3
    assert (m.wavelength_min_nm, m.wavelength_max_nm) == (643.8, 668.8)
    assert Instrument.from_preset("bench_600_13p9_o2").wavelength_min_nm == 400.0


def test_sky_binning_vectorised_matches_loop():
    sky = Sky()
    wl = np.array([500.0, 800.0, 1500.0, 1799.0])
    dl = np.array([0.01, 1.0, 5.0, 0.001])
    out = sky.radiance_per_bin(wl, dl)
    for i in range(len(wl)):
        m = (sky.wavelength_nm >= wl[i] - dl[i] / 2) & (sky.wavelength_nm <= wl[i] + dl[i] / 2)
        if m.any():
            assert abs(out[i] - sky._flux[m].mean() / 1000) < 1e-12 * max(1, out[i])
        else:
            j = np.argmin(np.abs(sky.wavelength_nm - wl[i]))
            assert abs(out[i] - sky._flux[j] / 1000) < 1e-12


def test_peak_pixel_independent_of_n_telescopes_in_independent_mode():
    src = Source.blackbody(5800, mag=8, band="V")
    p1 = LFASTETC("bench_r1800", telescope=Telescope(1, "50um")).snr(src, 100).peak_e_per_pixel
    p20 = LFASTETC("bench_r1800", telescope=Telescope(20, "50um")).snr(src, 100).peak_e_per_pixel
    np.testing.assert_allclose(p20, p1, rtol=1e-12)
    # fiber mode: all light lands on one detector, so peak scales with N (star-dominated)
    f20 = LFASTETC("bench_r1800", telescope=Telescope(20, "18um")).snr(src, 100).peak_e_per_pixel
    f1 = LFASTETC("bench_r1800", telescope=Telescope(1, "18um")).snr(src, 100).peak_e_per_pixel
    np.testing.assert_allclose(f20, 20 * f1 - 19 * 3e-4 * 100, rtol=1e-9)


def test_min_exptime_flag_and_all_selector():
    etc = LFASTETC("cshell_r18000", wavelength_min_nm=1600, wavelength_max_nm=1605)
    r = etc.exptime(Source.blackbody(5800, mag=4, band="H"), snr=50)
    assert r.below_min_exptime.any()                     # 8 s Fowler minimum
    src = Source.blackbody(5800, 10)
    assert etc.snr_grid(src, [10], [100], wavelength="all").shape == (1, 1, len(etc.wavelength))
    assert etc.limiting_magnitude(src, 100, 10, wavelength="all").shape == etc.wavelength.shape
    with pytest.raises(ValueError):
        etc.snr_grid(src, [10], [100], wavelength="median")


def test_nan_in_zero_throughput_bins_and_at_outside_grid():
    etc = LFASTETC("cshell_r18000")                      # includes the 1355-1420 nm telluric gap
    src = Source.blackbody(5800, mag=10, band="J")
    r = etc.exptime(src, snr=10)
    assert np.isnan(r.exptime).any()
    # outside the telluric cores: K_opt's 1 nm sampling clear, and the resolved transmission the
    # LSF path puts back (water-band lines are far deeper at R=20,000 than their 1 nm average)
    good = (etc.telescope.kopt(etc.wavelength) > 1e-3) & (etc.sky.transmission(etc.wavelength) > 0.5)
    assert np.nanmax(r.exptime[good]) < 1e7     # n_pix is now 36 for R=18,000 (#13)
    m = etc.limiting_magnitude(src, 600, 10, wavelength="all")
    assert np.isnan(m[~(etc.telescope.kopt(etc.wavelength) > 0)]).all()
    assert np.nanmin(m[good]) > 5
    assert np.isnan(r.at(2000.0)["snr"])


def test_saturation_flag():
    etc = LFASTETC("cshell_r18000", wavelength_min_nm=1600, wavelength_max_nm=1620)
    bright = Source.blackbody(5800, mag=2, band="H")
    assert etc.snr(bright, 600).saturated.any()
    assert not etc.snr(bright, 600, n_exp=600).saturated.any()


def test_uniform_grid_and_empty_range():
    etc = LFASTETC("bench_r1800", grid="uniform", step_nm=5.0)
    assert np.allclose(np.diff(etc.wavelength), 5.0)
    with pytest.raises(ValueError):
        LFASTETC("bench_r1800", wavelength_min_nm=900)


def test_validation_rejects_garbage():
    src = Source.blackbody(5800, 15)
    etc = LFASTETC("bench_r1800")
    for bad in (dict(n_exp=0), dict(n_exp=2.5)):
        with pytest.raises(ValueError):
            etc.snr(src, 100, **bad)
    with pytest.raises(ValueError):
        etc.snr(src, -10)
    with pytest.raises(ValueError):
        etc.exptime(src, snr=-5)
    with pytest.raises(ValueError):
        Telescope(2.5)
    with pytest.raises(ValueError):
        Telescope(20, diameter_m=-1)
    with pytest.raises(ValueError):
        Source.blackbody(5800, 10, system="bogus")
    with pytest.raises(ValueError):
        Source.from_arrays([400, 500], [1, 1], flux_unit="jy", mag=10)
    with pytest.raises(ValueError):
        LFASTETC("bench_r1800", telescope=Telescope(5), n_telescopes=20)
    with pytest.raises(ValueError):
        LFASTETC(Instrument.from_preset("bench_r1800", resolving_power=1e6), wavelength_min_nm=1, grid="uniform", step_nm=1e-9)
    with pytest.raises(ValueError):
        etc.snr_grid(src, [10], [10], wavelength=(500, 600, 700))
    with pytest.raises(ValueError):
        LFASTETC("bench_r1800", wavelength_min_nm=0.0)        # 0 is rejected, not treated as "unset"


def test_summary_types_exptime_vs_mag_callable_qe_and_plots():
    import matplotlib
    matplotlib.use("Agg")
    from lfast_etc.plotting import plot_snr_map, plot_limiting_magnitude
    src = Source.blackbody(5800, 15)
    etc = LFASTETC("bench_r1800")
    s = etc.snr(src, 100).summary()
    assert type(s["snr"][50]) is float and set(s["snr"]) == {10, 50, 90}
    t = etc.exptime_vs_mag(src, [10, 15, 20], snr=30, wavelength=(500, 600))
    assert t.shape == (3,) and np.all(np.diff(t) > 0)
    direct = etc.exptime(src.with_mag(15), snr=30)
    m = (etc.wavelength >= 500) & (etc.wavelength <= 600)
    assert abs(t[1] - np.nanmedian(direct.exptime[m])) < 1e-9
    qe_curve = lambda w: 0.9 - 0.5 * (np.asarray(w) - 412) / 300
    e2 = LFASTETC(Instrument.from_preset("bench_r1800", qe=qe_curve))
    assert e2.throughputs()["qe"][0] > e2.throughputs()["qe"][-1]
    ext = LFASTETC("bench_r1800", apply_extinction=True)
    assert np.all(ext.throughputs()["star"] <= etc.throughputs()["star"])
    assert LFASTETC("bench_r1800", sky=Sky(scale=0.5)).sky_rate_per_pixel()[0] == pytest.approx(0.5 * etc.sky_rate_per_pixel()[0])
    assert abs(Source.flat(20, "V", "AB", kind="flam").magnitude("V", "AB") - 20) < 1e-9
    wl = np.linspace(380, 900, 100)
    absolute = Source.from_arrays(wl, np.full_like(wl, 1e-15), flux_unit="W/m2/nm", band="V")
    assert abs(absolute.magnitude("V") - absolute.mag) < 1e-9
    grid = etc.snr_grid(src, [12, 16], [10, 100])
    plot_snr_map([12, 16], [10, 100], grid)
    plot_limiting_magnitude(etc, src, [60, 600], 10)
    direct.plot(); direct.plot(what="n_star"); direct.plot_rates()


def test_new_plot_helpers_and_exports():
    import matplotlib
    matplotlib.use("Agg")
    from lfast_etc.plotting import plot_vs_n_telescopes, sky_surface_brightness_ab
    src = Source.blackbody(5800, 15).add_line(656.3, 0.3, 0.05)
    etc = LFASTETC("bench_r1800")
    res = etc.snr(src, 600)
    res.plot_noise_budget(); res.plot_spectrum(seed=1); etc.plot_throughput(); etc.plot_sky(); src.plot()
    sim = res.simulate(seed=3)
    good = np.isfinite(res.snr)
    assert sim.shape == res.snr.shape and np.all(sim[good] > 0)
    # the realised scatter matches 1/SNR
    rel = (sim[good] - res.n_star[good] * 600) / (res.n_star[good] * 600) * res.snr[good]
    assert abs(np.std(rel) - 1) < 0.1
    t = res.to_table()
    assert len(t) == len(res.wavelength) and "snr" in t.colnames and t.meta["n_exp"] == 1
    sb = sky_surface_brightness_ab(etc)
    m = (etc.wavelength > 540) & (etc.wavelength < 560)
    assert 19.5 < np.nanmedian(sb[m]) < 21.5                 # grey-time V sky ~ 20.3 AB/arcsec^2
    plot_vs_n_telescopes(lambda n: LFASTETC("bench_r1800", n_telescopes=n), src, [5, 20, 80], [60, 600], 10,
                         wavelength=(500, 600))
    with pytest.raises(ValueError):
        LFASTETC("bench_r1800", n_telescopes=0)
    # a zero-throughput bin still collects sky and dark, so the bright-pixel estimate stays finite
    assert LFASTETC("cshell_r18000", wavelength_min_nm=1383, wavelength_max_nm=1383.05).reads_needed(src, 300) == 1
    with pytest.warns(UserWarning):
        Telescope(20, "50um", combine="fiber")


def test_review_fixes_bandpass_clip_selector_and_validation():
    src = Source.blackbody(5800, 15)
    # an explicit range never leaves the preset bandpass
    e = LFASTETC("cshell_r18000", wavelength_min_nm=1000, wavelength_max_nm=1900)
    assert e.wavelength.min() >= 1080 and e.wavelength.max() <= 1750
    with pytest.raises(ValueError):
        LFASTETC("cshell_r18000", wavelength_min_nm=500, wavelength_max_nm=600)
    # a float selector outside the grid raises instead of snapping to the edge bin
    b = LFASTETC("bench_r1800")
    with pytest.raises(ValueError):
        b.snr_grid(src, [15], [100], wavelength=1500.0)
    with pytest.raises(ValueError):
        b.limiting_magnitude(src, 100, 10, wavelength=1500.0)
    # at() interpolates between two valid bins, never across or into an invalid one (#37)
    r = LFASTETC("cshell_r18000", wavelength_min_nm=1368, wavelength_max_nm=1373).exptime(Source.blackbody(5800, 10, "J"), 10)
    i = np.flatnonzero(r.valid & np.roll(r.valid, 1))[1]
    assert np.isfinite(r.at(r.wavelength[i] - 1e-6)["exptime"])
    j = np.flatnonzero(r.valid & ~np.roll(r.valid, 1))[0]       # first valid bin after the telluric gap
    assert np.isnan(r.at(r.wavelength[j] - 1e-6)["exptime"])
    # validation
    for bad in (dict(mag=float("nan")), dict(mag=float("inf"))):
        with pytest.raises(ValueError):
            Source.blackbody(5800, **bad)
    with pytest.raises(ValueError):
        src.add_line(550, 0.3, 0)
    with pytest.raises(ValueError):
        src.add_line(550, float("nan"), 0.1)
    with pytest.raises(ValueError):
        b.snr(src, float("nan"))
    with pytest.raises(ValueError):
        b.reads_needed(src, 600, margin=0)
    with pytest.raises(ValueError):
        Instrument.from_preset("bench_r1800", qe=5)
    with pytest.raises(ValueError):
        Source.from_arrays([400, 500, 600], [1e-15, 1e-15, 1e-15], flux_unit="W/m2/nm", band="Q")
    assert Instrument.from_preset("bench_r1800", name="mine").name == "mine"


def test_peak_pixel_follows_the_brightest_packed_input():
    """#36: N=21 is one full 20-group plus one, so the brightest detector still sees 20 telescopes."""
    src = Source.blackbody(5800, 8, "V")
    T = 54.57          # s: a V=8 star on 20 telescopes fills the STA1 well 1.2x in one read (tophat peak)
    peak = lambda **kw: np.nanmax(LFASTETC("bench_r1800", telescope=Telescope(**kw)).snr(src, T).peak_e_per_pixel)
    full = peak(n_telescopes=20)
    for n in (21, 39, 40, 41):
        assert peak(n_telescopes=n) == pytest.approx(full, rel=1e-12)
    assert peak(n_telescopes=19) < full
    assert peak(n_telescopes=21, combine="independent") == pytest.approx(peak(n_telescopes=1, combine="independent"), rel=1e-12)
    ideal = peak(n_telescopes=40, telescopes_per_input=None)
    assert ideal == pytest.approx(2 * full - 3e-4 * T, rel=1e-9)   # all 40 on one detector (dark once)
    etc = LFASTETC("bench_r1800", telescope=Telescope(21))
    assert etc.snr(src, T).saturated.any() and etc.reads_needed(src, T) == 2
    # the co-added SNR is untouched: N=21 is still noisier than two full arrays
    assert np.nanmedian(etc.snr(src, 100).snr) < np.nanmedian(LFASTETC("bench_r1800", telescope=Telescope(40)).snr(src, 100).snr)


def test_reads_needed_margin_is_a_fraction_of_full_well():
    """#39: margin > 1 would let the brightest pixel exceed the full well."""
    etc = LFASTETC("cshell_r18000", wavelength_min_nm=1081, wavelength_max_nm=1086)
    src = Source.blackbody(4800, 3, "J")
    assert etc.reads_needed(src, 300, margin=1.0) <= etc.reads_needed(src, 300, margin=0.8)
    for bad in (5, 1.01, 0, -1, float("nan")):
        with pytest.raises(ValueError):
            etc.reads_needed(src, 300, margin=bad)


def test_instrument_without_full_well_builds_and_reports_saturation_unknown():
    """#46: the default full_well=inf must pass its own validator; None means the same."""
    ins = Instrument(name="custom", wavelength_min_nm=400, wavelength_max_nm=800,
                     resolving_power=1800, pixels_per_resel=2.2, spectrograph_efficiency=0.3)
    assert ins.full_well == float("inf")
    assert Instrument.from_preset("bench_r1800", full_well=None).full_well == float("inf")
    res = LFASTETC(ins).snr(Source.blackbody(5800, 0, "V"), 3600)
    assert not res.saturated.any() and res.summary()["saturated_fraction"] is None
    assert LFASTETC(ins).reads_needed(Source.blackbody(5800, 0, "V"), 3600) == 1
    assert LFASTETC("bench_r1800").snr(Source.blackbody(5800, 15), 100).summary()["saturated_fraction"] == 0.0
    with pytest.raises(ValueError):
        Instrument.from_preset("bench_r1800", full_well=0)


def test_overflow_is_flagged_not_returned():
    """#47: a finite but absurd SNR overflows the exposure-time quadratic; raise instead of returning inf."""
    etc = LFASTETC("bench_r1800")
    src = Source.blackbody(5800, 15)
    with pytest.raises(ValueError, match="overflow"):
        etc.exptime(src, snr=1e100)
    assert np.isfinite(etc.exptime(src, snr=1e6).exptime[np.isfinite(etc.exptime(src, snr=1e6).snr)]).all()


def test_science_notebook_windows_match_their_text():
    """#34, #35: the integration windows are derived from the widths the prose quotes."""
    import re, sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    from build_science_notebooks import NB
    for name, const, pattern in (("07_time_domain_lsst_era", "WIN_A", r"(\d+) A window"),
                                 ("10_cshell_nir_science", "RV_WINDOW_A", r"(\d+) A window")):
        text = "\n".join(s for k, s in NB[name] if k == "md")
        code = "\n".join(s for k, s in NB[name] if k == "code")
        stated = {int(m) for m in re.findall(pattern, text)}
        assert stated == {int(re.search(rf"{const} = (\d+)", code).group(1))}, (name, stated)
    assert "wavelength_max_nm=1590" not in "\n".join(s for _, s in NB["10_cshell_nir_science"])
    # a 40 A setting at 1550 nm on R=43k is ~111 resels, not ~1100
    assert 100 < len(LFASTETC("cshell_r43000", wavelength_min_nm=1550, wavelength_max_nm=1554).wavelength) < 120


def test_kopt_below_epsilon_is_zero_throughput():
    """#12: a K_opt sample of 4.6e-9 at 687 nm is a telluric core, not a finite throughput."""
    from lfast_etc.telescope import KOPT_EPS
    tel = Telescope()
    assert KOPT_EPS == 1e-6
    assert tel.kopt(687.0) == 0.0 and tel.kopt(687.3) == 0.0        # the sample and the interp ramp next to it
    assert tel.kopt(761.0) == 0.0 and tel.kopt(700.0) > 1e-3
    etc = LFASTETC("bench_316_7p0_o1")
    r = etc.exptime(Source.blackbody(5800, 15), 100)
    m = (r.wavelength > 680) & (r.wavelength < 695)
    assert np.isnan(r.exptime[m]).any()                               # the 687 nm core is invalid ...
    assert np.nanmax(r.exptime[m]) < 2 * np.nanmedian(r.exptime[m])    # ... not a 5x finite spike (B-band wing stays finite)
    clear = etc.telescope.kopt(r.wavelength) > 0
    assert not r.valid[~clear].any()                                    # no valid bin on a zero-K_opt sample ...
    assert r.valid.sum() >= clear.sum() - 10                            # ... and the LSF path only adds the bins touching a core
    assert (LFASTETC("bench_316_7p0_o1", lsf=None).exptime(Source.blackbody(5800, 15), 100).valid == clear).all()


def test_exptime_result_is_nan_in_invalid_bins_and_summary_follows():
    """#37: exptime() must not broadcast the requested SNR into opaque bins."""
    etc = LFASTETC("cshell_r18000", wavelength_min_nm=1368, wavelength_max_nm=1373)
    r = etc.exptime(Source.blackbody(5800, 10, "J"), 10)
    assert (~r.valid).any() and r.valid.any()
    assert np.isnan(r.snr[~r.valid]).all() and np.isnan(r.exptime[~r.valid]).all()
    assert np.allclose(r.snr[r.valid], 10)
    opaque = LFASTETC("cshell_r18000", wavelength_min_nm=1369.5, wavelength_max_nm=1370.5).exptime(Source.blackbody(5800, 10, "J"), 10)
    assert not opaque.valid.any()
    assert all(np.isnan(v) for v in opaque.summary()["snr"].values())


def test_at_exact_centres_and_one_bin_grid():
    """#37: at() returns the sample at an exact centre (valid -> value, invalid -> NaN), also on a one-bin grid."""
    etc = LFASTETC("cshell_r18000", wavelength_min_nm=1368, wavelength_max_nm=1373)
    r = etc.exptime(Source.blackbody(5800, 10, "J"), 10)
    i, j = np.flatnonzero(~r.valid)[0], np.flatnonzero(r.valid)[0]
    bad = r.at(r.wavelength[i])
    assert np.isnan(bad["exptime"]) and np.isnan(bad["snr"]) and bad["n_star"] == 0.0
    good = r.at(r.wavelength[j])
    assert good["exptime"] == r.exptime[j] and good["snr"] == 10.0
    one = LFASTETC("bench_r1800", wavelength_min_nm=550, wavelength_max_nm=550.2).snr(Source.blackbody(5800, 15), 100)
    assert len(one.wavelength) == 1
    assert one.at(one.wavelength[0])["snr"] == one.snr[0] and np.isfinite(one.snr[0])
    assert np.isnan(one.at(one.wavelength[0] + 0.01)["snr"])      # off-centre on a one-bin grid: nothing to interpolate


# ---- #13 / #19: extraction aperture, detector profile, seeing --------------
def test_n_pix_is_the_full_aperture_with_a_circular_fiber_image_by_default():
    b = Instrument.from_preset("bench_r1800")
    assert b.fiber_image_px is None and b.spatial_px == 2.2 and b.n_pix == pytest.approx(2.2 * 2.2)
    assert Instrument.from_preset("bench_r1800", fiber_image_px=1.0).n_pix == 2.2
    assert Instrument.from_preset("bench_r1800", fiber_image_px=3.0).n_pix == pytest.approx(6.6)
    # CSHELL slit modes narrow the dispersion direction only: the spatial height stays the 6 px fiber image
    assert Instrument.from_preset("cshell_r43000").n_pix == pytest.approx(2.5 * 6.0)
    assert Instrument.from_preset("cshell_r18000").n_pix == pytest.approx(6.0 * 6.0)
    with pytest.raises(ValueError):
        Instrument.from_preset("bench_r1800", fiber_image_px=0)


def test_peak_fraction_by_profile():
    from lfast_etc.instrument import peak_pixel_fraction
    u, t, g = (peak_pixel_fraction(2.5, 2.5, p) for p in ("uniform", "tophat", "gaussian"))
    assert u == pytest.approx(1 / 6.25)
    assert u < t < 1 and 0 < g < 1       # a Gaussian of the same FWHM has wings, so its peak can sit below uniform
    assert Instrument.from_preset("bench_r1800").profile == "tophat"
    assert Instrument.from_preset("bench_r1800", profile="uniform").peak_fraction == pytest.approx(1 / 2.2 ** 2)
    # a big disc: the peak pixel holds about 1 / (pi a b) of the light
    assert peak_pixel_fraction(12, 12, "tophat") == pytest.approx(4 / (np.pi * 144), rel=0.15)
    # a wide Gaussian must not be truncated by the render grid: erf^2 of the central pixel
    from scipy.special import erf
    sig = 12 / 2.3548
    assert peak_pixel_fraction(12, 12, "gaussian") == pytest.approx(erf(0.5 / (sig * np.sqrt(2))) ** 2, rel=2e-3)
    assert not Instrument.from_preset("bench_r1800").pixel_image().flags.writeable
    # a sharp Gaussian: nearly everything lands in one pixel
    assert peak_pixel_fraction(0.3, 0.3, "gaussian") > 0.9
    with pytest.raises(ValueError):
        Instrument.from_preset("bench_r1800", profile="airy")


def test_uniform_profile_reproduces_the_uniform_share_peak():
    ins = Instrument.from_preset("bench_r1800", profile="uniform")
    etc = LFASTETC(ins); src = Source.blackbody(5800, 8)
    r = etc.snr(src, 10)
    share = etc.telescope.telescopes_on_brightest_input / etc.telescope.n_telescopes
    expect = ((r.n_star / r.n_pix + r.n_sky_pix) * share + ins.dark_current) * 10
    np.testing.assert_allclose(r.peak_e_per_pixel, expect, rtol=1e-12)
    # top hat peaks higher than the uniform share, so it saturates earlier
    rt = LFASTETC("bench_r1800").snr(src, 10)
    assert np.nanmax(rt.peak_e_per_pixel) > np.nanmax(r.peak_e_per_pixel)
    assert rt.n_pix == r.n_pix


def test_seeing_scales_the_point_source_path_only():
    ref = Telescope()
    assert ref.seeing_arcsec is None and ref.seeing_factor(550.0) == 1.0
    assert Telescope(seeing_arcsec=1.0).seeing_factor(550.0) == pytest.approx(1.0)   # the 1.0" reference
    bad, good = Telescope(seeing_arcsec=2.0), Telescope(seeing_arcsec=0.6)
    assert bad.seeing_factor(550.0) < 1.0 < good.seeing_factor(550.0)
    assert good.seeing_factor(550.0) <= 1 / ref.fiber_coupling(1.0, 550.0)
    assert ref.fiber_coupling(0.01, 550.0) == pytest.approx(1.0, abs=1e-3)   # all the light inside the fiber
    assert ref.fiber_coupling(1.0, 550.0, psf="gaussian") != ref.fiber_coupling(1.0, 550.0, psf="moffat")
    assert bad.seeing_factor(1600.0) > bad.seeing_factor(550.0)                 # seeing improves with wavelength
    s_ref = ref.implied_seeing_ref(550.0)
    assert 1.3 < s_ref < 1.7 and ref.fiber_coupling(s_ref, 550.0) == pytest.approx(ref.eta_fib(550.0), rel=1e-6)
    src = Source.blackbody(5800, 15)
    e0, e2 = LFASTETC("bench_r1800"), LFASTETC("bench_r1800", telescope=bad)
    assert np.nanmedian(e2.snr(src, 100).snr) < np.nanmedian(e0.snr(src, 100).snr)
    np.testing.assert_allclose(e2.sky_rate_per_pixel(), e0.sky_rate_per_pixel())
    with pytest.raises(ValueError):
        Telescope(seeing_arcsec=-1)
    with pytest.raises(ValueError):
        Telescope(seeing_psf="airy")


def test_profile_and_seeing_plots_return_axes():
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ax = Instrument.from_preset("cshell_r43000").plot_profile()
    assert ax.get_images() and "peak" in ax.get_title()
    ax2 = Telescope(seeing_arcsec=1.5).plot_fiber_coupling()
    assert len(ax2.get_lines()) >= 2
    plt.close("all")

def test_bench_presets_use_the_sta1_full_well():
    """WIYN STA1 CCD page: full well 87,500 e-, 1% linearity to 63,000 e-, 3.7 e- read noise at medium gain."""
    for name in (n for n in list_instruments() if n.startswith("bench_")):
        ins = Instrument.from_preset(name)
        assert ins.full_well == 87500.0 and ins.read_noise == 3.7


# ---- Choi et al. 2026 throughput cascade and calibrated seeing ----------------
def test_throughput_cascade_multiplies_to_kopt():
    from lfast_etc.telescope import CHOI_STAGES
    tel = Telescope(seeing_arcsec=1.5)
    w = np.array([450.0, 550.0, 700.0, 1000.0, 1250.0, 1650.0])
    c = tel.throughput_cascade(w)
    assert [s for s, _ in CHOI_STAGES] and all(0 < v <= 1 for _, v in CHOI_STAGES)
    stages = c["stages"]; names = [s["name"] for s in stages]
    assert names[0].startswith("atmosphere") and any("injection" in n for n in names) and any("seeing" in n for n in names)
    for s in stages:
        assert np.all(s["transmission"] >= 0) and np.all(s["transmission"] <= 1.0001)
    np.testing.assert_allclose(c["cumulative"][-1], tel.throughput_point_source(w), rtol=1e-12)
    # the derived residual (mirror + PFC + fiber attenuation) is a plausible 0.6-0.8 in the mid band
    res = [s for s in stages if "mirror" in s["name"]][0]["transmission"]
    assert np.all((res[1:5] > 0.6) & (res[1:5] < 0.8))
    # sky path: no injection, no seeing
    cs = tel.throughput_cascade(w, path="sky")
    np.testing.assert_allclose(cs["cumulative"][-1], tel.throughput_extended(w), rtol=1e-12)
    assert not any("injection" in s["name"] for s in cs["stages"])
    # scalar wavelength returns floats
    c1 = tel.throughput_cascade(550.0)
    assert isinstance(c1["cumulative"][-1], float)


def test_etc_cascade_adds_instrument_and_plots():
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    etc = LFASTETC("bench_r1800", telescope=Telescope(seeing_arcsec=1.2))
    c = etc.throughput_cascade()
    names = [s["name"] for s in c["stages"]]
    assert any("spectrograph" in n for n in names) and any("QE" in n for n in names)
    np.testing.assert_allclose(c["cumulative"][-1] * etc.telescope.area_m2 * 20 * etc.delta_lambda,
                               etc.throughputs()["star"], rtol=1e-12)
    ax = etc.plot_throughput_cascade()
    assert len(ax.get_lines()) >= 8
    plt.close("all")


def test_choi_seeing_model_is_calibrated_and_gentler_than_moffat():
    ref = Telescope()                                   # default seeing_psf is "choi"
    assert ref.seeing_psf == "choi"
    assert Telescope(seeing_arcsec=1.0).seeing_factor(550.0) == pytest.approx(1.0, abs=1e-6)
    t15 = Telescope(seeing_arcsec=1.5)
    f_choi = t15.seeing_factor(550.0)
    f_moffat = Telescope(seeing_arcsec=1.5, seeing_psf="moffat").seeing_factor(550.0)
    assert 0.65 < f_choi < 0.78 and f_choi > f_moffat         # the non-seeing blur softens the seeing dependence
    s = np.array([0.6, 1.0, 1.5, 2.0, 2.5])
    f = np.array([Telescope(seeing_arcsec=x).seeing_factor(800.0) for x in s])
    assert np.all(np.diff(f) < 0)                               # monotonic in seeing
    assert Telescope(seeing_arcsec=2.0).seeing_factor(1650.0) > Telescope(seeing_arcsec=2.0).seeing_factor(550.0)
    # the 50 um fiber barely cares about seeing
    with pytest.warns(UserWarning):
        t50 = Telescope(fiber="50um", seeing_arcsec=2.0, combine="fiber")
    assert 0.75 < t50.seeing_factor(550.0) < 0.95 and t50.seeing_factor(550.0) > Telescope(seeing_arcsec=2.0).seeing_factor(550.0) + 0.25
    # out-of-table seeing is clamped, not extrapolated
    assert Telescope(seeing_arcsec=5.0).seeing_factor(550.0) == pytest.approx(Telescope(seeing_arcsec=3.0).seeing_factor(550.0))
