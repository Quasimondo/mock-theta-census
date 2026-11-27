"""Addendum v22: F(z) = sum z^n q^(n(n+1)/2) (-q;q)_n at even roots of unity q0 = e^(2 pi i/X), z = e^(2 pi i a), 0 < a < 1.
Theorem 14:  F(z; q0) = i P e^(-pi i a) [ J(a/b) - I((a-1)/b) + z^X I(a/b) ]   at x = 1/X   (finite sum on the left, n < X/2).
The three integrals at real x, double precision, Phi_b from Faddeev's integral (fad.py):
  J(a/b)      over the line L (delta = 10 degrees);
  I((a-1)/b)  over the contour of Theorem 13, tilted version: ray Im w = 2 pi up to P, segment P -> i pi (angle 85 degrees), top edge, vertical line through the saddle, real axis;
  I(a/b)      over the horizontal line through the lower saddle (Theorem 7 of v19).
Usage: python3 n1_zroots.py"""
import numpy as np
from mpmath import mp, mpf, mpc, expjpi
from fad import logPhi_np
XS, WS = np.polynomial.legendre.leggauss(12)

def F_fin(a, X):
    with mp.workdps(30 + int(0.0702*X) + 10):
        q = expjpi(mpf(2)/X); z = expjpi(2*mpf(a)); s = mpc(0); term = mpc(1); qn = mpc(1)
        for n in range(X//2):
            s += term; qn *= q; term *= qn*z*(1 + qn)
        return complex(s)

def roots(a):
    d = np.sqrt(1 + 4*np.exp(-2j*np.pi*a) + 0j); W = (-1 + d)/2
    return (W, -1 - W) if W.imag > 0 else (-1 - W, W)            # (upper, lower)

def integrand(w, X, ap, kernel=False):
    x = 1.0/X; t = -2j*np.pi*x; s = w/(1j*t)
    lg = np.pi**2/(12*t) + t/24 - 0.5*np.log(2) + t*s*s/2 - logPhi_np(w, x) + 2*np.pi*ap*s
    if kernel:
        z = 2*np.pi*s
        lg -= (z + np.log1p(np.exp(-z))) if z.real > 0 else np.log1p(np.exp(z))
    return np.exp(lg)

def seg(f, wfun, dw, a, b, panels):
    tot = 0j; edges = np.linspace(a, b, panels + 1)
    for u, v in zip(edges[:-1], edges[1:]):
        for xi, wi in zip(XS, WS):
            tot += wi*(v - u)/2*f(wfun((u + v)/2 + (v - u)/2*xi))*dw
    return tot

def PI_growing(a, X):
    t = -2j*np.pi/X; T = abs(t); ap = a - 1; Wu, _ = roots(a); w0 = np.log(Wu); r0 = w0.real
    f = lambda w: integrand(w, X, ap)
    d1 = np.radians(85.0); e = np.exp(-1j*d1); rhoP = np.pi/np.sin(d1); P = 1j*np.pi - rhoP*e
    tot = seg(f, lambda r: r + 2j*np.pi, 1.0, P.real - (14*T/(2*np.pi*a) + 2), P.real, 400)      # C1: decay e^(-2 pi a |r|/T)
    tot += seg(f, lambda rho: 1j*np.pi - rho*e, e, 0.0, rhoP, 300)*(1)                            # C2 from P to i pi: rho decreasing, dw = -e d(rho); integrate rho from rhoP to 0
    tot += seg(f, lambda r: r + 1j*np.pi, 1.0, 0.0, r0, 200) if r0 != 0 else 0                    # C3 (signed: from i pi to r0 + i pi)
    tot += seg(f, lambda eta: r0 + 1j*eta, -1j, 0.0, np.pi, 400)                                  # C4 downwards: dw = -i d(eta)
    tot += seg(f, lambda r: r + 0j, 1.0, r0, r0 + 14*T/(2*np.pi*(1 - a)) + 2, 400)                # C5: decay e^(-2 pi (1-a) r/T)
    return tot*np.exp(t/8)/(1j*t)

def PI_decaying(a, X):
    t = -2j*np.pi/X; T = abs(t); _, Wl = roots(a); w1 = np.log(Wl)
    R = 14*T/min(abs(np.angle(1 + Wl)), abs(np.angle(1 + Wl) - np.angle(Wl))) + 2             # the two asymptotic slopes of Im V on the line
    return seg(lambda w: integrand(w, X, a), lambda r: r + 1j*w1.imag, 1.0, w1.real - R, w1.real + R, 800)*np.exp(t/8)/(1j*t)

def PJ(a, X):
    t = -2j*np.pi/X; ga = np.radians(10.0); R = 6.0*np.sqrt(X)
    f = lambda w: integrand(w, X, a, kernel=True)
    return seg(f, lambda r: 1j*t*r*np.exp(-1j*ga), 1j*t*np.exp(-1j*ga), -R, R, 500)*np.exp(t/8)/(1j*t)

if __name__ == "__main__":
    print("Theorem 14: finite sum F(z; e^(2 pi i/X)) against i e^(-pi i a) [P J(a/b) - P I((a-1)/b) + z^X P I(a/b)]   (double precision)")
    for name, a in (("i", 0.25), ("omega", 1/3), ("-omega", 5/6), ("e^(2 pi i 0.3)", 0.3), ("e^(2 pi i 0.7)", 0.7), ("-1", 0.5)):
        for X in (8, 16, 24, 40):
            F = F_fin(a, X); J = PJ(a, X); G = PI_growing(a, X); D = PI_decaying(a, X)
            rhs = 1j*np.exp(-1j*np.pi*a)*(J - G + np.exp(2j*np.pi*a*X)*D)
            print(f"    z = {name:14s} X = {X:2d}: F = {F:.8f}, relative difference {abs(F - rhs)/abs(F):.1e}   (|P J| = {abs(J):.3f}, |P I((a-1)/b)| = {abs(G):.3e}, |P I(a/b)| = {abs(D):.3e})")
