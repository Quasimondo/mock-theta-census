"""Growing integral I((a-1)/b), 0 < a < 1, at real x (tau = 0): test of a contour through the selected saddle.
V(w) = -w^2/2 - 2 pi i a' w + Li2(-e^w), a' = a - 1; height H = -Im V; saddle w0 = Log W, W the upper root of W^2 + W = e^(-2 pi i a); H(w0) = Vol.
Contour:  G1: ray  i pi - rho e^(-i d1), rho from infinity to 0 (left half-plane, above the strip; Li2(-e^w) with |e^w| < 1)
          G2a: top edge of the strip between i pi and r0 + i pi          (closed form: Im V = -2 pi a r for r >= 0, (eta* - pi) r for r <= 0, eta* = 2 pi (1 - a))
          G2b: vertical line Re w = r0 from Im w = pi down through w0 to 0  (Im V monotone on each side of w0)
          G3: real axis from r0 to +infinity                              (Im V = 2 pi (1 - a) r)
Reports, for each a and several ray angles d1: min over G1 of Im V + Vol (must be > 0), and sup over the contour of (H - H0)/min(d^2, d).
Usage: python3 m1_contour.py"""
import numpy as np
from mpmath import mp, mpf, mpc, polylog, pi, log, exp, sqrt, im, re, arg
mp.dps = 20
def li2(z): return complex(polylog(2, complex(z)))
def setup(a):
    c = np.exp(-2j*np.pi*a); d = np.sqrt(1 + 4*c + 0j); W = (-1 + d)/2
    if W.imag < 0: W = (-1 - d)/2
    w0 = np.log(W); ap = a - 1
    V0 = -w0*w0/2 - 2j*np.pi*ap*w0 + li2(-W)
    return W, w0, ap, V0
def V_any(w, ap): return -w*w/2 - 2j*np.pi*ap*w + li2(-np.exp(w))       # principal Li2: valid in the open strip and for Re w < 0
def scan(a, d1deg):
    W, w0, ap, V0 = setup(a); r0, e0 = w0.real, w0.imag; vol = -V0.imag; etas = 2*np.pi*(1 - a)
    d1 = np.radians(d1deg); out = {}
    # G1
    rho = np.linspace(1e-4, 14, 2801); w = 1j*np.pi - rho*np.exp(-1j*d1)
    imV = np.array([V_any(x, ap).imag for x in w]); out["G1 min ImV+Vol"] = (imV + vol).min(); out["G1 argmin rho"] = rho[(imV + vol).argmin()]
    pts = [(w, imV)]
    # G2a (closed form) and its check slightly inside the strip
    r = np.linspace(min(r0, 0), max(r0, 0), 400); cf = np.where(r >= 0, -2*np.pi*a*r, (etas - np.pi)*r)
    chk = np.array([V_any(x + 1j*(np.pi - 1e-9), ap).imag for x in r]); out["G2a formula defect"] = np.abs(chk - cf).max()
    pts.append((r + 1j*np.pi, cf))
    # G2b
    eta = np.linspace(1e-6, np.pi - 1e-6, 1500); wv = r0 + 1j*eta; imv = np.array([V_any(x, ap).imag for x in wv]); pts.append((wv, imv))
    below = imv[eta < e0]; above = imv[eta > e0]
    out["G2b monotone"] = bool(np.all(np.diff(below) < 1e-12) and np.all(np.diff(above) > -1e-12))
    # G3
    rr = np.linspace(r0, r0 + 14, 800); pts.append((rr + 0j, 2*np.pi*(1 - a)*rr))
    out["G3 formula defect"] = max(abs(V_any(x, ap).imag - 2*np.pi*(1 - a)*x) for x in rr[::80])
    sup = -1e9
    for ww, iv in pts:
        d = np.abs(ww - w0); m = d > 1e-3
        sup = max(sup, ((-iv[m]) - vol).__truediv__(np.minimum(d[m]**2, d[m])).max())
    out["sup ratio"] = sup; out["vol"] = vol; out["r0"] = r0; out["eta0/pi"] = e0/np.pi
    return out
if __name__ == "__main__":
    print("  a     r0     eta0/pi  Vol     | d1:  min over G1 of (Im V + Vol)  [sup over contour of (H-H0)/min(d^2,d)]")
    for a in (0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 0.95, 0.98):
        row = []
        for d1 in (30, 60, 75, 85):
            o = scan(a, d1); row.append(f"{d1}: {o['G1 min ImV+Vol']:+.3f} [{o['sup ratio']:+.3f}]")
        print(f"  {a:.2f}  {o['r0']:+.3f}  {o['eta0/pi']:.3f}   {o['vol']:.4f}  | " + "   ".join(row) + f"   | G2b monotone: {o['G2b monotone']}, formula defects {o['G2a formula defect']:.0e}, {o['G3 formula defect']:.0e}")
