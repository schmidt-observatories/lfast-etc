"""lfast_etc: Exposure Time Calculator for the LFAST telescope array.

Quick start::

    from lfast_etc import LFASTETC, Source
    etc = LFASTETC("bench_r1800", n_telescopes=20)
    src = Source.blackbody(5800, mag=15, band="V")
    res = etc.exptime(src, snr=50)
    res.plot()
"""
import importlib.metadata

try:
    __version__ = importlib.metadata.version("lfast-etc")
except importlib.metadata.PackageNotFoundError:  # pragma: no cover
    __version__ = "0.0.0"

from .instrument import Instrument, list_instruments
from .telescope import Telescope
from .source import Source, BANDS
from .sky import Sky
from .etc import LFASTETC, ETCResult

__all__ = [
    "__version__", "LFASTETC", "ETCResult", "Instrument", "list_instruments",
    "Telescope", "Source", "Sky", "BANDS",
]
