"""Physical constants (SI) and package data paths."""
from pathlib import Path

H_PLANCK = 6.62607015e-34   # J s
C_LIGHT = 2.99792458e8      # m/s
K_BOLTZMANN = 1.380649e-23  # J/K

DATA_DIR = Path(__file__).parent / "data"
THROUGHPUT_DIR = DATA_DIR / "throughput"
SKY_DIR = DATA_DIR / "sky"
CONFIG_DIR = DATA_DIR / "config"
