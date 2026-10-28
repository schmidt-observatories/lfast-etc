"""Verbatim port of the reference LFAST ETC notebook (2026-10-12) (Bench configuration).

Used only as the regression oracle for lfast_etc.
"""
import numpy as np
from scipy.integrate import quad
from astropy.io import fits
from lfast_etc.constants import THROUGHPUT_DIR, SKY_DIR

H_PLANCK = 6.62607015e-34
C_LIGHT = 2.99792458e8
K_BOLTZMANN = 1.380649e-23


def run(LAMBDA_MIN_NM=400.0, LAMBDA_MAX_NM=1800.0, T_EFF_K=5800.0, V_MAG=15.0, SNR_TARGET=50.0,
        N_TELESCOPES=20, S_opt=0.3, R_POW=1800, N_PIX=2.2, DARK=3e-4, RN=3.7, QE=0.75,
        FIBER_DIAMETER_ARCSEC=1.48, suffix="18", spectral_lines=({'center_nm': 656.3, 'depth': 0.3, 'sigma_nm': 0.05},)):
    TELESCOPE_AREA_M2 = np.pi * (0.726 / 2.0) ** 2
    OMEGA = np.pi * (FIBER_DIAMETER_ARCSEC / 2) ** 2
    wavelengths = [LAMBDA_MIN_NM]
    while wavelengths[-1] < LAMBDA_MAX_NM:
        wavelengths.append(wavelengths[-1] + wavelengths[-1] / R_POW)
    wl = np.array(wavelengths)
    dl = wl / R_POW

    def I(w, T):
        wm = np.asarray(w, float) * 1e-9
        return 2 * H_PLANCK * C_LIGHT ** 2 / (wm ** 5 * (np.exp(H_PLANCK * C_LIGHT / (wm * K_BOLTZMANN * T)) - 1)) * 1e-9

    f_nu = 3640.0 * 10 ** (-0.4 * V_MAG) * 1e-26
    dnu = C_LIGHT * (90e-9) / (550e-9) ** 2
    F_V = f_nu * dnu
    integral, _ = quad(lambda w: I(w, T_EFF_K), 505, 595)
    alpha = F_V / integral
    F = alpha * I(wl, T_EFF_K)
    for line in spectral_lines:
        F = F * (1 - line['depth'] * np.exp(-0.5 * ((wl - line['center_nm']) / line['sigma_nm']) ** 2))

    with fits.open(SKY_DIR / "skytable_secz1.3.fits") as h:
        sl = np.array(h[1].data['lam']); sf = np.array(h[1].data['flux'])
    o = np.argsort(sl); sl, sf = sl[o], sf[o]
    K = np.interp(wl, np.load(THROUGHPUT_DIR / "Kopt_wavelength.npy"), np.load(THROUGHPUT_DIR / f"Kopt_value_{suffix}.npy"))
    eta = np.interp(wl, np.load(THROUGHPUT_DIR / "Eta_fib_wavelength.npy"), np.load(THROUGHPUT_DIR / f"Eta_fib_value_{suffix}.npy"))

    E = H_PLANCK * C_LIGHT / (wl * 1e-9)
    N_star = F / E * TELESCOPE_AREA_M2 * K * S_opt * dl * QE * N_TELESCOPES

    half = dl / 2
    lo = np.searchsorted(sl, wl - half, 'left'); hi = np.searchsorted(sl, wl + half, 'right')
    rad = np.zeros(len(wl))
    for i in range(len(wl)):
        if hi[i] > lo[i]:
            rad[i] = np.mean(sf[lo[i]:hi[i]])
        else:
            j = np.clip(np.searchsorted(sl, wl[i]), 1, len(sl) - 1)
            rad[i] = sf[j - 1] if abs(wl[i] - sl[j - 1]) <= abs(sl[j] - wl[i]) else sf[j]
    n_sky = rad / 1000.0 * TELESCOPE_AREA_M2 * OMEGA * K / eta * S_opt * dl * QE * N_TELESCOPES / N_PIX

    N_star = np.maximum(N_star, 1e-20); n_sky = np.maximum(n_sky, 1e-20)
    B = N_PIX * (n_sky + DARK)
    a = N_star ** 2; b = -(SNR_TARGET ** 2) * (N_star + B); c = -(SNR_TARGET ** 2) * N_PIX * RN ** 2
    t = (-b + np.sqrt(b ** 2 - 4 * a * c)) / (2 * a)
    return wl, t, N_star, n_sky
