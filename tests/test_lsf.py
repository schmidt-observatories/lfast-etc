"""Line-spread function and detector-pixel grid (``LFASTETC(lsf=...)``, ``ETCResult.pixel``)."""
import numpy as np
import pytest
from lfast_etc import LFASTETC, Source, Instrument, PixelResult, LSFS


SRC = Source.blackbody(5800.0, 15.0, "V")


def test_continuum_per_resel_numbers_barely_move_and_pixels_sum_to_resels():
    """A smooth source: the convolved, pixel-integrated rates agree with the boxcar reference per
    resel to 1e-3, and the pixels of one resel sum to that resel (integer pixels per resel)."""
    box = LFASTETC("bench_r1800", lsf=None).exptime(SRC, snr=50)
    for lsf in LSFS:
        r = LFASTETC("bench_r1800", lsf=lsf).exptime(SRC, snr=50)
        rel = np.abs(r.n_star / box.n_star - 1)
        assert np.median(rel) < 1e-4 and rel.max() < 3e-2          # the tail: bins on resolved telluric lines
        assert np.median(np.abs(r.exptime / box.exptime - 1)) < 2e-4
        assert isinstance(r.pixel, PixelResult) and len(r.pixel.wavelength) == pytest.approx(len(r.wavelength) * 2.2, abs=3)
    # CSHELL: 6 pixels per resel, pixel 6j centred on resel j, so the resel spans pixels
    # 6j-3 .. 6j+3 with half of each end pixel: that weighted sum is the resel's rate
    etc = LFASTETC("cshell_r18000", wavelength_min_nm=1600, wavelength_max_nm=1604)
    r = etc.snr(SRC.with_mag(10.0), exptime=600)
    px = r.pixel
    j, k = 20, 120
    assert px.wavelength[k] == pytest.approx(etc.wavelength[j], rel=1e-7)   # (1+1/6R)^6 vs 1+1/R
    w = np.array([0.5, 1, 1, 1, 1, 1, 0.5])
    assert (px.n_star[k - 3:k + 4] * w).sum() == pytest.approx(r.n_star[j], rel=1e-3)
    assert (px.n_sky_pix[k - 3:k + 4] * w).sum() / 6 == pytest.approx(r.n_sky_pix[j], rel=1e-3)   # per physical pixel, both
    # the resolved atmosphere gives a continuum source structure *within* a resel (the boxcar cannot)
    assert px.n_star[k - 3:k + 4].max() / px.n_star[k - 3:k + 4].min() > 1.02
    # SNR per extracted pixel x sqrt(pixels per resel) is the resel's SNR, up to that in-resel slope
    assert px.snr[k] * np.sqrt(6) == pytest.approx(r.snr[j], rel=1e-2)
    assert px.exptime[k] == r.exptime[j]
    # on a flat part of the Bench it is exact
    rb = LFASTETC("bench_r1800").snr(SRC, exptime=600)
    jb = int(np.argmin(np.abs(rb.wavelength - 550))); kb = int(np.argmin(np.abs(rb.pixel.wavelength - 550)))
    assert rb.pixel.snr[kb] * np.sqrt(2.2) == pytest.approx(rb.snr[jb], rel=2e-4)


def test_narrow_line_gets_the_instrument_width_and_keeps_its_equivalent_width():
    """A 0.02 nm emission line at 1602 nm on CSHELL R=18,000 (resel 0.089 nm, 6 px): the boxcar
    puts it in one resel; with the Gaussian LSF it spans ~6 pixels at FWHM one resel, and the
    extra electrons above the continuum are the same in both."""
    line = Source.blackbody(5800.0, 10.0, "H").add_line(1602.0, -5.0, 0.02)     # emission, 5x continuum peak
    cont = Source.blackbody(5800.0, 10.0, "H")
    kw = dict(wavelength_min_nm=1600, wavelength_max_nm=1604)
    box_l, box_c = (LFASTETC("cshell_r18000", lsf=None, **kw).snr(s, 600).n_star for s in (line, cont))
    etc = LFASTETC("cshell_r18000", lsf="gaussian", **kw)
    g_l, g_c = (etc.snr(s, 600) for s in (line, cont))
    ew_box, ew_resel = (box_l - box_c).sum(), (g_l.n_star - g_c.n_star).sum()
    excess = g_l.pixel.n_star - g_c.pixel.n_star
    assert excess.sum() == pytest.approx(ew_resel, rel=1e-3)                                      # pixels and resels agree
    # in electrons the LSF path differs from the boxcar by the resolved transmission at the line
    # (the boxcar sees the 1 nm-averaged atmosphere): the ratio is that local factor
    local = etc.sky.transmission(1602.0) / np.clip(etc.sky.transmission_smoothed(np.array([1602.0]))[0], 0.2, 1)
    assert ew_resel / ew_box == pytest.approx(local, rel=0.03)
    above_half = excess > 0.5 * excess.max()
    assert 5 <= above_half.sum() <= 7                                                              # FWHM ~ 6 pixels
    # top hat: same equivalent width, flatter top (more pixels near the peak)
    t_l, t_c = (LFASTETC("cshell_r18000", lsf="tophat", **kw).snr(s, 600) for s in (line, cont))
    ex_t = t_l.pixel.n_star - t_c.pixel.n_star
    assert ex_t.sum() == pytest.approx(excess.sum(), rel=1e-3)
    assert (ex_t > 0.8 * ex_t.max()).sum() > (excess > 0.8 * excess.max()).sum()
    k = int(excess.argmax())
    assert ex_t[k - 5] + ex_t[k + 5] < 0.6 * (excess[k - 5] + excess[k + 5])     # and lighter wings


def test_telluric_cores_stay_invalid_per_resel_and_per_pixel():
    """The 687 nm O2 B band on bench_316_7p0_o1 (#12): zero throughput stays NaN with the LSF."""
    box = LFASTETC("bench_316_7p0_o1", lsf=None).exptime(SRC, snr=50)
    r = LFASTETC("bench_316_7p0_o1").exptime(SRC, snr=50)
    assert np.isnan(box.exptime).sum() > 0
    assert np.isnan(r.exptime)[np.isnan(box.exptime)].all()                 # every boxcar NaN is a NaN here
    assert np.isnan(r.exptime).sum() <= np.isnan(box.exptime).sum() + 6     # plus at most the bins touching a core
    assert np.nanmax(r.exptime) < 10 * np.nanmax(box.exptime)               # and no 1e-6 wing survives
    px = r.pixel
    assert np.isnan(px.exptime).sum() > 0 and np.isnan(px.snr).sum() == np.isnan(px.exptime).sum()
    assert np.all(px.n_star[np.isnan(px.exptime)] == 0)
    assert np.all(np.isfinite(px.snr[px.valid]))


def test_resolved_atmosphere_replaces_the_1nm_one():
    """On the fine grid K_opt's 1 nm atmosphere is divided out and the SkyCalc transmission put
    back resolved: the factor applied has telluric structure (std of several percent at R=43,000
    in the H band) but no net effect (mean 1), and the band total barely moves."""
    kw = dict(wavelength_min_nm=1600, wavelength_max_nm=1604)
    etc = LFASTETC("cshell_r43000", **kw)
    ratio = etc.sky.transmission(etc._fine_w) / np.clip(etc.sky.transmission_smoothed(etc._fine_w), 0.2, 1)
    assert ratio.std() > 0.02 and abs(ratio.mean() - 1) < 0.02
    a = LFASTETC("cshell_r43000", lsf=None, **kw).snr(SRC, 600).n_star
    b = etc.snr(SRC, 600).n_star
    assert b.sum() == pytest.approx(a.sum(), rel=0.02)
    assert (b / a).std() > 0.01                        # redistributed between resels, not just rescaled


def test_pixel_result_simulate_table_and_validation():
    r = LFASTETC("cshell_r18000", wavelength_min_nm=1600, wavelength_max_nm=1604).snr(SRC.with_mag(10.0), 600)
    px = r.pixel
    sim = px.simulate(seed=1)
    assert sim.shape == px.wavelength.shape and np.isfinite(sim[px.valid]).all()
    resid = (sim - px.n_star * px.exptime)[px.valid] / np.sqrt(px.variance()[px.valid])
    assert abs(resid.std() - 1) < 0.2
    assert len(px.to_table()) == len(px.wavelength)
    assert "gaussian LSF" in r.setup
    with pytest.raises(ValueError):
        LFASTETC("bench_r1800", lsf="sinc")
    u = LFASTETC("bench_r1800", grid="uniform")           # pixels need the resel grid: boxcar only
    assert u.lsf is None and u.exptime(SRC, 50).pixel is None


def test_sky_radiance_is_convolved_per_resel_and_per_pixel():
    """The He/OH airglow lines near 1083 nm on CSHELL R=18,000: the raw table (R=20,000, one
    sample per line) is spikier than the LSF-convolved radiance, which spans ~6 pixels, keeps the
    line flux, and is what sky_surface_brightness_ab reports."""
    from lfast_etc.plotting import sky_surface_brightness_ab
    kw = dict(wavelength_min_nm=1082, wavelength_max_nm=1085)
    box, etc = LFASTETC("cshell_r18000", lsf=None, **kw), LFASTETC("cshell_r18000", **kw)
    raw, resel, px = box.sky_radiance(), etc.sky_radiance("resel"), etc.sky_radiance("pixel")
    band = float(etc.sky.radiance_per_bin(np.array([1083.5]), np.array([3.0]))[0])                # the table's own band mean
    assert px.mean() == pytest.approx(band, rel=0.03) and resel.mean() == pytest.approx(band, rel=0.03)   # flux kept
    assert abs(raw.mean() / band - 1) > abs(resel.mean() / band - 1)                               # nearest-sampling was not
    assert px.max() < raw.max() and px.std() / px.mean() < raw.std() / raw.mean()                 # lines spread out
    k = int(px.argmax())
    assert (px[k - 10:k + 11] > 0.5 * px[k]).sum() >= 5                                             # ~one resel wide
    sb = sky_surface_brightness_ab(etc, per="pixel")
    assert sb.shape == etc.pixel_wavelength.shape and np.isfinite(sb).all()
    assert sky_surface_brightness_ab(etc).shape == etc.wavelength.shape
    with pytest.raises(ValueError):
        box.sky_radiance("pixel")
