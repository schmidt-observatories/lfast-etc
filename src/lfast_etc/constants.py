"""Physical constants (SI) and package data paths."""
from pathlib import Path
import numpy as np

H_PLANCK = 6.62607015e-34   # J s
C_LIGHT = 2.99792458e8      # m/s
K_BOLTZMANN = 1.380649e-23  # J/K

DATA_DIR = Path(__file__).parent / "data"


def positive(name, value, allow_zero=False):
    """Validate a scalar parameter at the public boundary; returns it as float."""
    try:
        v = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"{name} must be a number, got {value!r}") from None
    if not np.isfinite(v) or v < 0 or (v == 0 and not allow_zero):
        raise ValueError(f"{name} must be {'>= 0' if allow_zero else '> 0'}, got {value!r}")
    return v
THROUGHPUT_DIR = DATA_DIR / "throughput"
SKY_DIR = DATA_DIR / "sky"
CONFIG_DIR = DATA_DIR / "config"
