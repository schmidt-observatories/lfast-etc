import os, shutil, sys
from pathlib import Path
import importlib.metadata

sys.path.insert(0, os.path.abspath("../src"))
project = "lfast-etc"
copyright = "2026, Gudmundur Stefansson / Schmidt Sciences"
author = "Gudmundur Stefansson"
release = importlib.metadata.version("lfast-etc")

extensions = [
    "sphinx.ext.autodoc", "sphinx.ext.autosummary", "sphinx.ext.napoleon",
    "sphinx.ext.mathjax", "sphinx.ext.viewcode", "sphinx_design",
    "sphinx_copybutton", "nbsphinx", "myst_parser",
    "matplotlib.sphinxext.plot_directive",
]
nbsphinx_execute = "never"          # notebooks are committed with outputs
autodoc_member_order = "bysource"
autosummary_generate = True
autoclass_content = "both"
napoleon_use_param = True
templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints"]
source_suffix = {".rst": "restructuredtext", ".md": "markdown", ".ipynb": "jupyter_notebook"}

# ---- figures generated at build time (.. plot:: directives) ----
plot_include_source = True
plot_html_show_source_link = False
plot_html_show_formats = False
plot_formats = [("png", 110)]
plot_pre_code = """
import numpy as np
import matplotlib.pyplot as plt
from lfast_etc import LFASTETC, Source, Telescope, Instrument, Sky, list_instruments
from lfast_etc.plotting import *
"""

# ---- theme: same as the wcc-etc docs (Read the Docs theme, SOS logo) ----
html_theme = "sphinx_rtd_theme"
html_theme_options = {
    "logo_only": False,
    "navigation_depth": 3,
    "collapse_navigation": False,
    "sticky_navigation": True,
    "style_external_links": True,
}
html_static_path = ["_static"]
html_logo = "logo_schmidt_observatory_system.png"
html_title = f"lfast-etc {release}"
html_css_files = ["custom.css"]

# ---- pull the executed notebooks from ../notebooks into docs/notebooks at build time ----
_root = Path(__file__).parent
_dst = _root / "notebooks"
_dst.mkdir(exist_ok=True)
for nb in (_root.parent / "notebooks").glob("*.ipynb"):
    shutil.copy(nb, _dst / nb.name)

# ---- generate the instrument preset table from the TOML files ----
def _write_preset_table():
    from lfast_etc import list_instruments, Instrument
    rows = []
    for name in list_instruments():
        i = Instrument.from_preset(name)
        rows.append((f"``{name}``", f"{i.wavelength_min_nm:.0f}-{i.wavelength_max_nm:.0f}", f"{i.resolving_power:,.0f}",
                     f"{i.pixels_per_resel:g}", f"{i.spectrograph_efficiency:.3g}", f"{i.slit_throughput:g}",
                     f"{i.qe:g}" if not callable(i.qe) else "curve", f"{i.read_noise:g}", f"{i.dark_current:g}",
                     f"{i.full_well:,.0f}", f"{i.min_exptime:g}", i.description.replace("|", "/")))
    head = ("preset", "nm", "R", "pix/resel", "S_opt", "slit", "QE", "RN e-", "dark e-/s/pix", "full well", "min read s", "description")
    lines = [".. list-table::", "   :header-rows: 1", "   :widths: 18 8 6 6 6 5 5 6 8 8 7 40", ""]
    for r in (head, *rows):
        lines.append("   * - " + r[0])
        lines.extend("     - " + c for c in r[1:])
    (_root / "_generated").mkdir(exist_ok=True)
    (_root / "_generated" / "presets.rst").write_text("\n".join(lines) + "\n")

_write_preset_table()
