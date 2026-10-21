import os, shutil, sys
from pathlib import Path

sys.path.insert(0, os.path.abspath("../src"))
project = "lfast-etc"
copyright = "2026, Gudmundur Stefansson / Schmidt Sciences"
author = "Gudmundur Stefansson"
import importlib.metadata
release = importlib.metadata.version("lfast-etc")

extensions = [
    "sphinx.ext.autodoc", "sphinx.ext.autosummary", "sphinx.ext.napoleon",
    "sphinx.ext.mathjax", "sphinx.ext.viewcode", "sphinx_design",
    "sphinx_copybutton", "nbsphinx", "myst_parser",
]
nbsphinx_execute = "never"          # notebooks are committed with outputs
autodoc_member_order = "bysource"
autosummary_generate = True
autoclass_content = "both"
templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints"]
source_suffix = {".rst": "restructuredtext", ".md": "markdown", ".ipynb": "jupyter_notebook"}
html_theme = "sphinx_book_theme"
html_theme_options = {
    "show_toc_level": 2,
    "repository_url": "https://github.com/schmidt-observatories/lfast-etc",
    "use_repository_button": True,
}
html_static_path = ["_static"]
html_title = "lfast-etc"

# Pull the executed notebooks from ../notebooks into docs/notebooks at build time.
_src = Path(__file__).parent.parent / "notebooks"
_dst = Path(__file__).parent / "notebooks"
_dst.mkdir(exist_ok=True)
for nb in _src.glob("*.ipynb"):
    shutil.copy(nb, _dst / nb.name)
