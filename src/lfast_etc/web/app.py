"""Browser front-end for :mod:`lfast_etc` (Flask + Bokeh).

Run with ``lfast-etc-web [--host H] [--port P]`` or ``python -m lfast_etc.web.app``.
"""
import argparse
import math

import numpy as np
from flask import Flask, jsonify, render_template, request
from bokeh.embed import components
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.plotting import figure
from bokeh.resources import CDN

from lfast_etc import BANDS, LFASTETC, Instrument, Source, Telescope, list_instruments

PURPLE, TEAL, GREY, RED = "#5B2A86", "#00798c", "#6b6b6b", "#c62828"
MAX_PLOT_POINTS = 20000
MAX_API_POINTS = 2000

DEFAULTS = dict(
    instrument="bench_r1800", resolving_power="", pixels_per_resel="",
    spectrograph_efficiency="", slit_throughput="", qe="", read_noise="", dark_current="",
    n_telescopes="20", fiber="18um", combine="auto",
    source_type="blackbody", teff="5800", mag="15", band="V", system="Vega", alpha="0",
    mode="exptime", snr="50", exptime="600", n_exp="1", wmin="", wmax="",
)
OVERRIDE_FIELDS = ("resolving_power", "pixels_per_resel", "spectrograph_efficiency",
                   "slit_throughput", "qe", "read_noise", "dark_current")


# ---- parsing -----------------------------------------------------------------
PLACEHOLDERS = {n: {k: getattr(Instrument.from_preset(n), k) for k in OVERRIDE_FIELDS}
                for n in list_instruments()}


def _num(form, key, lo=None, hi=None, kind=float, required=True):
    raw = (form.get(key) or "").strip()
    if not raw:
        if required:
            raise ValueError(f"'{key}' is required")
        return None
    try:
        v = float(raw)
    except ValueError:
        raise ValueError(f"'{key}' must be a number, got {raw!r}") from None
    if not math.isfinite(v):
        raise ValueError(f"'{key}' must be finite")
    if kind is int:
        if v != int(v):
            raise ValueError(f"'{key}' must be an integer, got {raw!r}")
        v = int(v)
    if lo is not None and v < lo:
        raise ValueError(f"'{key}' must be >= {lo}")
    if hi is not None and v > hi:
        raise ValueError(f"'{key}' must be <= {hi}")
    return v


def _choice(form, key, options):
    v = form.get(key, DEFAULTS[key])
    if v not in options:
        raise ValueError(f"'{key}' must be one of {list(options)}")
    return v


def build(form):
    """Parse a form/query mapping, run the ETC and return (etc, src, res, params)."""
    name = _choice(form, "instrument", list_instruments())
    overrides = {}
    limits = {"resolving_power": 1e6, "pixels_per_resel": 1e4, "spectrograph_efficiency": 1.0,
              "slit_throughput": 1.0, "qe": 1.0, "read_noise": 1e6, "dark_current": 1e6}
    for k in OVERRIDE_FIELDS:
        v = _num(form, k, lo=0, hi=limits.get(k), required=False)
        if v is not None:
            overrides[k] = v
    ins = Instrument.from_preset(name, **overrides)   # Instrument validates > 0 where required

    combine = _choice(form, "combine", ("auto", "fiber", "independent"))
    tel = Telescope(n_telescopes=_num(form, "n_telescopes", lo=1, hi=100000, kind=int),
                    fiber=_choice(form, "fiber", ("18um", "50um")),
                    combine=None if combine == "auto" else combine)

    wmin = _num(form, "wmin", lo=0, required=False)
    wmax = _num(form, "wmax", lo=0, required=False)
    etc = LFASTETC(ins, telescope=tel, wavelength_min_nm=wmin, wavelength_max_nm=wmax)

    band = _choice(form, "band", tuple(BANDS))
    system = _choice(form, "system", ("Vega", "AB"))
    mag = _num(form, "mag", lo=-30, hi=60)
    stype = _choice(form, "source_type", ("blackbody", "flat", "powerlaw"))
    if stype == "blackbody":
        src = Source.blackbody(_num(form, "teff", lo=100, hi=1e6), mag, band, system)
    elif stype == "flat":
        src = Source.flat(mag, band, system)
    else:
        src = Source.powerlaw(_num(form, "alpha", lo=-20, hi=20), mag, band, system)

    n_exp = _num(form, "n_exp", lo=1, hi=100000, kind=int)
    mode = _choice(form, "mode", ("exptime", "snr"))
    if mode == "exptime":
        res = etc.exptime(src, snr=_num(form, "snr", lo=1e-6, hi=1e6), n_exp=n_exp)
    else:
        res = etc.snr(src, exptime=_num(form, "exptime", lo=1e-6, hi=1e9), n_exp=n_exp)
    return etc, src, res, mode


def _finite(v):
    """JSON-safe number(s): every NaN/inf becomes ``null`` (#47)."""
    if v is None:
        return None
    if np.ndim(v):
        return [_finite(x) for x in v]
    return float(v) if np.isfinite(v) else None


def summarize(etc, src, res, mode="snr"):
    s = res.summary()
    lim = etc.limiting_magnitude(src, exptime=3600, snr=10, wavelength=BANDS[src.band][0]) \
        if etc.wavelength.min() <= BANDS[src.band][0] <= etc.wavelength.max() else \
        etc.limiting_magnitude(src, exptime=3600, snr=10)
    solved = res.exptime if mode == "exptime" else res.snr
    return {
        "exptime": {str(k): _finite(v) for k, v in s["exptime"].items()},
        "snr": {str(k): _finite(v) for k, v in s["snr"].items()},
        "saturated_fraction": _finite(s["saturated_fraction"]),
        "below_min_exptime_fraction": _finite(s["below_min_exptime_fraction"]),
        "nan_fraction": float(np.mean(np.isnan(solved))),
        "n_bins": s["n_bins"],
        "limiting_mag_snr10_1hr": _finite(lim),
        "limiting_mag_band": src.band,
        "limiting_mag_system": src.system,
    }


def _decimate(n, limit):
    return max(1, int(np.ceil(n / limit)))


# ---- plots -------------------------------------------------------------------
def _fig(title, ylabel, log=True):
    p = figure(title=title, x_axis_label="Wavelength (nm)", y_axis_label=ylabel,
               y_axis_type="log" if log else "linear", height=300, sizing_mode="stretch_width",
               tools="pan,wheel_zoom,box_zoom,reset,save", toolbar_location="right")
    p.title.text_color = PURPLE
    p.title.text_font_size = "13px"
    p.grid.grid_line_alpha = 0.3
    return p


def make_plots(etc, res, mode):
    step = _decimate(len(res.wavelength), MAX_PLOT_POINTS)
    sl = slice(None, None, step)
    w = res.wavelength[sl]
    y = (res.exptime if mode == "exptime" else res.snr)[sl]
    ylabel = "Total exposure time (s)" if mode == "exptime" else "SNR per resel"
    src1 = ColumnDataSource(dict(w=w, y=y, snr=res.snr[sl], t=res.exptime[sl],
                                 peak=res.peak_e_per_pixel[sl]))
    p1 = _fig(f"{ylabel} vs wavelength", ylabel)
    r1 = p1.line("w", "y", source=src1, color=PURPLE, line_width=1.5)
    p1.add_tools(HoverTool(renderers=[r1], mode="vline", tooltips=[
        ("λ", "@w{0.00} nm"), ("exptime", "@t{0.00} s"), ("SNR", "@snr{0.00}"),
        ("peak", "@peak{0.0} e-/pix")]))
    sat = res.saturated[sl]
    if sat.any():
        rs = p1.scatter(w[sat], y[sat], color=RED, size=4, legend_label="saturated")
        p1.add_tools(HoverTool(renderers=[rs], tooltips=[("saturated at", "$x{0.00} nm")]))
        p1.legend.location = "top_left"

    src2 = ColumnDataSource(dict(w=w, star=res.n_star[sl] / res.n_pix, sky=res.n_sky_pix[sl]))
    p2 = _fig("Electron rates per pixel", "e- / s / pixel")
    rstar = p2.line("w", "star", source=src2, color=PURPLE, line_width=1.5, legend_label="star")
    rsky = p2.line("w", "sky", source=src2, color=TEAL, line_width=1.5, legend_label="sky")
    dark = max(res.dark_pix, 1e-12)
    p2.line([w.min(), w.max()], [dark, dark], color=GREY, line_dash="dashed", legend_label="dark")
    p2.add_tools(HoverTool(renderers=[rstar, rsky], mode="vline", tooltips=[
        ("λ", "@w{0.00} nm"), ("star", "@star{0.000e}"), ("sky", "@sky{0.000e}")]))
    p2.legend.location = "bottom_right"
    p2.legend.click_policy = "hide"

    src3 = ColumnDataSource(dict(w=w, kopt=etc.telescope.kopt(w), eta=etc.telescope.eta_fib(w)))
    p3 = _fig("Telescope throughput", "Throughput", log=False)
    rk = p3.line("w", "kopt", source=src3, color=PURPLE, line_width=1.5, legend_label="K_opt")
    re_ = p3.line("w", "eta", source=src3, color=TEAL, line_width=1.5, legend_label="eta_fib")
    p3.add_tools(HoverTool(renderers=[rk, re_], mode="vline", tooltips=[
        ("λ", "@w{0.00} nm"), ("K_opt", "@kopt{0.000}"), ("eta_fib", "@eta{0.000}")]))
    p3.legend.location = "top_right"

    script, divs = components((p1, p2, p3))
    note = (f"Plots decimated by {step}x ({len(res.wavelength)} bins) for display."
            if step > 1 else "")
    return script, divs, note


# ---- app ---------------------------------------------------------------------
def create_app():
    app = Flask(__name__)

    @app.route("/", methods=["GET", "POST"])
    def index():
        ctx = dict(instruments=list_instruments(), bands=list(BANDS),
                   bokeh_head=CDN.render(), error=None, summary=None, describe=None,
                   plot_script="", plot_divs=(), plot_note="", mode_label="")
        form = dict(DEFAULTS)
        if request.method == "POST":
            form.update({k: v for k, v in request.form.items() if k in DEFAULTS})
            try:
                etc, src, res, mode = build(form)
                ctx["summary"] = summarize(etc, src, res, mode)
                ctx["describe"] = f"{src.describe()}\n{etc.describe()}"
                ctx["plot_script"], ctx["plot_divs"], ctx["plot_note"] = make_plots(etc, res, mode)
                ctx["mode_label"] = mode
            except Exception as e:  # never 500 on user input
                ctx["error"] = str(e) or e.__class__.__name__
        ctx["placeholders"] = PLACEHOLDERS
        return render_template("index.html", form=form, **ctx)

    @app.route("/api/calc")
    def api_calc():
        form = dict(DEFAULTS)
        form.update({k: v for k, v in request.args.items() if k in DEFAULTS})
        try:
            etc, src, res, mode = build(form)
        except Exception as e:
            return jsonify(error=str(e) or e.__class__.__name__), 400
        step = _decimate(len(res.wavelength), MAX_API_POINTS)
        sl = slice(None, None, step)
        return jsonify(
            mode=mode, decimation=step, describe=etc.describe(), source=src.describe(),
            wavelength=res.wavelength[sl].tolist(), exptime=_finite(res.exptime[sl]),
            snr=_finite(res.snr[sl]), saturated=res.saturated[sl].tolist(),
            below_min_exptime=res.below_min_exptime[sl].tolist(),
            summary=summarize(etc, src, res, mode),
        )

    return app


def main():
    ap = argparse.ArgumentParser(description="LFAST ETC web app")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=5000)
    args = ap.parse_args()
    create_app().run(host=args.host, port=args.port, debug=False)


if __name__ == "__main__":
    main()
