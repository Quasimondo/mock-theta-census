"""Growing integral I((a-1)/b), 0 < a < 1: the contour
   C1: horizontal ray Im w = 2 pi from -infinity to P = i pi - rho_P e^(-i d1)  (rho_P sin d1 = pi)       [d1 = 90 deg: P = 2 pi i]
   C2: segment from P down to i pi                                                                      [d1 = 90 deg: the imaginary axis]
   C3: top edge of the strip between i pi and r0 + i pi
   C4: vertical line Re w = r0 from Im w = pi down through the saddle w0 to 0
   C5: real axis from r0 to +infinity
Height H_tau = -Im V + tau Re V (tau = tan(arg x)); H_tau(w0) is the saddle value.  Reports sup over the contour of (H - H0)/min(d^2, d).
Closed forms at tau = 0: C1 (d1 = 90): Im V = 2 pi a |r|;  C2 (d1 = 90): Im V = Cl2(eta - pi);  C3: -2 pi a r (r >= 0), (eta* - pi) r (r <= 0);  C5: 2 pi (1 - a) r.
Usage: python3 m2_contour.py"""
import numpy as np
from mpmath import mp, polylog, clsin
from m1_contour import setup, li2
mp.dps = 20
def V_left(w, ap):  return -w*w/2 - 2j*np.pi*ap*w + li2(-np.exp(w))                      # Re w <= 0, any Im w (periodic Li2 part)
def V_strip(w, ap): return -w*w/2 - 2j*np.pi*ap*w + li2(-np.exp(w))                      # open strip, principal
def V_top(r, ap):                                                                        # Im w = pi (boundary value from inside the strip)
    w = r + 1j*np.pi; x = np.exp(r)
    L = li2(x) if r <= 0 else complex(np.pi**2/3 - 0.5*r*r - li2(1/x).real, -np.pi*r)    # Li2(x - i0) for x > 1
    return -w*w/2 - 2j*np.pi*ap*w + L
def contour(a, d1deg, n=500):
    W, w0, ap, V0 = setup(a); r0 = w0.real; d1 = np.radians(d1deg)
    P = 1j*np.pi - (np.pi/np.sin(d1))*np.exp(-1j*d1)
    pts = []
    rr = np.linspace(P.real - 12, P.real, n);      pts += [("C1", x + 2j*np.pi, V_left(x + 2j*np.pi, ap)) for x in rr]
    ss = np.linspace(0, 1, n)[1:-1];               pts += [("C2", 1j*np.pi + s*(P - 1j*np.pi), V_left(1j*np.pi + s*(P - 1j*np.pi), ap)) for s in ss]
    tt = np.linspace(min(r0, 0), max(r0, 0), n//3); pts += [("C3", x + 1j*np.pi, V_top(x, ap)) for x in tt]
    ee = np.linspace(1e-6, np.pi - 1e-6, n);       pts += [("C4", r0 + 1j*e, V_strip(r0 + 1j*e, ap)) for e in ee]
    xx = np.linspace(r0, r0 + 12, n//2);           pts += [("C5", x + 0j, V_strip(x + 0j, ap)) for x in xx]
    return w0, V0, pts
def sup_ratio(a, d1deg, tau):
    w0, V0, pts = contour(a, d1deg); H0 = -V0.imag + tau*V0.real; worst = {}
    for name, w, V in pts:
        d = abs(w - w0)
        if d < 1e-3: continue
        val = ((-V.imag + tau*V.real) - H0)/min(d*d, d)
        worst[name] = max(worst.get(name, -1e9), val)
    return worst
if __name__ == "__main__":
    print("closed forms at tau = 0, d1 = 90 deg (a = 0.1): ", end="")
    w0, V0, pts = contour(0.1, 90.0, 200); dev = {}
    for name, w, V in pts:
        if name == "C1": dev[name] = max(dev.get(name, 0), abs(V.imag - 2*np.pi*0.1*abs(w.real)))
        if name == "C2": dev[name] = max(dev.get(name, 0), abs(V.imag - float(clsin(2, w.imag - np.pi))))
        if name == "C3": dev[name] = max(dev.get(name, 0), abs(V.imag + 2*np.pi*0.1*w.real))
        if name == "C5": dev[name] = max(dev.get(name, 0), abs(V.imag - 2*np.pi*0.9*w.real))
    print(", ".join(f"{k}: {v:.0e}" for k, v in dev.items()))
    for d1, tau_deg in ((90.0, 0.0), (85.0, 0.0), (85.0, 3.0), (85.0, 5.0)):
        tau = np.tan(np.radians(tau_deg))
        print(f"ray angle d1 = {d1} deg, arg x = {tau_deg} deg:  sup of (H - H0)/min(d^2, d) per piece  (all must be < 0)")
        for a in (0.01, 0.02, 0.05, 0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9, 0.98, 0.99):
            wst = sup_ratio(a, d1, tau)
            print(f"    a = {a:.2f}: " + "  ".join(f"{k} {v:+.4f}" for k, v in wst.items()) + ("   <-- FAILS" if max(wst.values()) >= 0 else ""))
