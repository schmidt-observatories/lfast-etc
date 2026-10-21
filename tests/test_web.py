import pytest

from lfast_etc.web.app import create_app, DEFAULTS


@pytest.fixture
def client():
    app = create_app()
    app.config["TESTING"] = True
    return app.test_client()


def test_get_index_shows_form(client):
    r = client.get("/")
    assert r.status_code == 200
    assert b"<form" in r.data and b'name="mag"' in r.data


def test_post_defaults_shows_summary(client):
    r = client.post("/", data=DEFAULTS)
    assert r.status_code == 200
    assert b"Summary" in r.data and b"bench_r1800" in r.data
    assert b"Input error" not in r.data


def test_post_bad_input_shows_error(client):
    r = client.post("/", data={**DEFAULTS, "mag": "abc"})
    assert r.status_code == 200
    assert b"Input error" in r.data and b"mag" in r.data


def test_post_snr_mode_and_overrides(client):
    r = client.post("/", data={**DEFAULTS, "mode": "snr", "read_noise": "0", "source_type": "powerlaw"})
    assert r.status_code == 200
    assert b"Summary" in r.data and b"RN=0.0" in r.data


def test_api_calc_json(client):
    r = client.get("/api/calc?instrument=bench_r1800&mag=14&snr=20")
    assert r.status_code == 200
    j = r.get_json()
    assert "wavelength" in j and len(j["wavelength"]) == len(j["exptime"]) <= 2000
    assert j["summary"]["n_bins"] > 0


def test_api_calc_bad_input(client):
    r = client.get("/api/calc?mag=abc")
    assert r.status_code == 400 and "error" in r.get_json()


def test_api_nan_bins_are_null(client):
    r = client.get("/api/calc?instrument=cshell_r18000&mag=10&teff=5800&band=J&mode=exptime&snr=10")
    assert r.status_code == 200
    assert b"NaN" not in r.data
    d = r.get_json()
    assert any(v is None for v in d["exptime"])
    assert "below_min_exptime" in d


def test_web_rejects_huge_resolving_power_and_fractional_telescopes(client):
    r = client.get("/api/calc?instrument=cshell_r43000&resolving_power=1e9&mode=snr&exptime=300")
    assert r.status_code == 400
    r = client.get("/api/calc?instrument=bench_r1800&n_telescopes=20.7&mode=snr&exptime=300")
    assert r.status_code == 400 and b"integer" in r.data


def test_web_independent_mode_and_n_exp(client):
    r = client.post("/", data=dict(instrument="cshell_r18000", n_telescopes="20", fiber="50um", combine="independent",
                                   source_type="blackbody", teff="4500", mag="9", band="H", system="Vega",
                                   mode="snr", exptime="300", n_exp="3", wmin="1550", wmax="1560"))
    assert r.status_code == 200 and b"Summary" in r.data and b"Input error" not in r.data
