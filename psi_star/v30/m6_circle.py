"""Fourier modes of the full matrix Omega_gamma(a, x) on a circle around a point of (-1/2, 0) and, for comparison, around a point of (-inf, -1/2).
   Upper half of the circle: series.  Lower half: state integrals (lowint.py)."""
import numpy as np
from m5_jump import Om_up, Om_low
for x0, r, n in ((-0.2, 0.1, 16), (-0.75, 0.1, 16)):
    th = 2*np.pi*(np.arange(n) + 0.5)/n
    vals = np.array([Om_up(x0 + r*np.exp(1j*t)) if np.sin(t) > 0 else Om_low(x0 + r*np.exp(1j*t)) for t in th])
    def mode(k): return (vals*np.exp(-1j*k*th)[:, None, None]).mean(axis=0)
    neg = [f"{np.abs(mode(-k)[0]).max():.1e} | {np.abs(mode(-k)[1:]).max():.1e}" for k in (1, 2, 3)]
    pos = [f"{np.abs(mode(k)[0]).max():.1e} | {np.abs(mode(k)[1:]).max():.1e}" for k in (0, 1, 2, 3, 6)]
    print(f"x0 = {x0}, r = {r}, {n} points (first row | block):\n   negative modes k = 1, 2, 3: " + ";  ".join(neg) + "\n   modes k = 0, 1, 2, 3, 6:    " + ";  ".join(pos))
