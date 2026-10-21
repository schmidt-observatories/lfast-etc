import numpy as np
import pytest

from lfast_etc import LFASTETC, Source, Telescope, Instrument, list_instruments, Sky
from reference_notebook import run as reference_run


def test_regression_against_reference_notebook():
    """Bench, 400-1800 nm resel grid, 5800 K V=15 star with an H-alpha line, SNR=50."""
    wl_ref, t_ref, nstar_ref, nsky_ref = reference_run()
    etc = LFASTETC("bench_r1800", n_telescopes=20, wavelength_min_nm=400, wavelength_max_nm=1800)
    src = Source.blackbody(5800, mag=15, band="V").add_line(656.3, 0.3, 0.05)
    res = etc.exptime(src, snr=50)
    # the reference grid overshoots lambda_max by one step; compare on the common bins
    n = min(len(wl_ref), len(res.wavelength))
    np.testing.assert_allclose(res.wavelength[:n], wl_ref[:n], rtol=1e-9)
    np.testing.assert_allclose(res.n_star[:n], nstar_ref[:n], rtol=1e-6, atol=1e-19)
    np.testing.assert_allclose(res.n_sky_pix[:n], nsky_ref[:n], rtol=1e-6, atol=1e-19)  # cumsum rounding; eps clamp
    np.testing.assert_allclose(res.exptime[:n], t_ref[:n], rtol=1e-6)


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
