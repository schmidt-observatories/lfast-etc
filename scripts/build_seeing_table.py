"""Build the Choi-calibrated seeing-factor table shipped as data/throughput/seeing_factor_choi.npz.

Choi et al. (2026, SPIE 14147, 1414736) compute the fiber injection efficiency eta_fib as the
encircled energy of PSF_system = seeing (1" Moffat, beta 4.765) (x) mirror figure (x) wind shake
(x) vibration inside the fiber. The shipped eta_fib curves are that product at 1" seeing. To scale
them to other seeing we need the non-seeing part of the blur. We approximate it as a Gaussian
kernel of FWHM g(lambda) plus a grey floor f(lambda), fitted at each wavelength to reproduce BOTH
shipped curves (18 um, 1.48" and 50 um, 4.1" fibers) at 1" seeing; the two fibers give the two
constraints. The table then holds

    factor(s, lambda, fiber) = EE_fiber[Moffat(s * (lambda/500)^-0.2) (x) G(g)] / EE_fiber[Moffat(1" ...) (x) G(g)]

for seeing s at 500 nm from 0.4 to 3.0" and wavelengths 400 to 1800 nm. Runtime a few minutes.
"""
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "src/lfast_etc/data/throughput"
BETA, SEEING_REF = 4.765, 1.0
FIBERS = {"18": 1.48, "50": 4.1}
S = np.round(np.arange(0.4, 3.01, 0.1), 2)
W = np.arange(400.0, 1801.0, 50.0)

n, dx = 768, 0.035                       # +-13.4", 0.035"/px
ax = (np.arange(n) - n / 2) * dx
X, Y = np.meshgrid(ax, ax); R = np.hypot(X, Y)
masks = {k: R <= d / 2 for k, d in FIBERS.items()}


def moffat(fwhm):
    a = fwhm / (2 * np.sqrt(2 ** (1 / BETA) - 1)); p = (1 + (R / a) ** 2) ** (-BETA); return p / p.sum()


def gauss(fwhm):
    s = fwhm / 2.3548; p = np.exp(-0.5 * (R / s) ** 2); return p / p.sum()


def conv(a, b):
    return np.fft.fftshift(np.real(np.fft.ifft2(np.fft.fft2(np.fft.ifftshift(a)) * np.fft.fft2(np.fft.ifftshift(b)))))


def ee(p, k):
    return float(p[masks[k]].sum())


w_eta = np.load(DATA / "Eta_fib_wavelength.npy")
eta = {k: np.load(DATA / f"Eta_fib_value_{k}.npy") for k in FIBERS}

g_tab, f_tab = np.zeros(len(W)), np.zeros(len(W))
factor = {k: np.zeros((len(S), len(W))) for k in FIBERS}
for j, wl in enumerate(W):
    e18, e50 = (float(np.interp(wl, w_eta, eta[k])) for k in ("18", "50"))
    see_ref = moffat(SEEING_REF * (wl / 500.0) ** -0.2)
    def ratio(g):
        p = conv(see_ref, gauss(g)); return ee(p, "18") / ee(p, "50")
    g = brentq(lambda g: ratio(g) - e18 / e50, 0.05, 4.0, xtol=1e-3)
    p0 = conv(see_ref, gauss(g)); f = e50 / ee(p0, "50")
    g_tab[j], f_tab[j] = g, f
    ref = {k: ee(p0, k) for k in FIBERS}
    for i, s in enumerate(S):
        p = conv(moffat(s * (wl / 500.0) ** -0.2), gauss(g))
        for k in FIBERS:
            factor[k][i, j] = ee(p, k) / ref[k]
    print(f"{wl:6.0f} nm  g = {g:.2f}\"  floor = {f:.3f}  factor(1.5\") 18um = {factor['18'][np.argmin(abs(S-1.5)), j]:.3f}")

np.savez(DATA / "seeing_factor_choi.npz", seeing=S, wavelength=W, kernel_fwhm=g_tab, floor=f_tab,
         factor_18=factor["18"], factor_50=factor["50"], beta=BETA, seeing_ref=SEEING_REF)
print("wrote", DATA / "seeing_factor_choi.npz")
