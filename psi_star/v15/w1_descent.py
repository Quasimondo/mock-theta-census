"""Addendum v15, the descent inequality (Proposition) for straight contours.
Contour s = r e^(-i gamma), w = i t s = +- rho e^(i alpha);  beta = (theta + 90 deg)/2,  alpha = beta - delta,  gamma = beta + delta.
Leading exponent:  psi(rho) = Re( e^(-i theta) G(w) ) - pi rho cos(gamma),
    G(w) = -w^2/2 + Li2(-e^w) + pi^2/12  (Re w <= 0),   G(w) = -w^2 - Li2(-e^(-w)) - pi^2/12  (Re w >= 0).
Exact formula (X + iY = int_0^1 log(1 + e^(-lambda |Re-side| ...)) dlambda, see the addendum):
    right half:  psi/rho = -(pi - Y) cos(gamma) - X sin(gamma) - rho sin(2 delta),        X + iY = int_0^1 log(1 + e^(-lambda w)) dlambda
    left half:   psi/rho = -(pi + Y) cos(gamma) + X sin(gamma) - (rho/2) sin(2 delta),    X + iY = int_0^1 log(1 + e^(+lambda w)) dlambda
Bounds:  |Y| <= pi/2,   lambda_alpha <= X <= ln 2,   lambda_alpha = ln(1 - e^(-(pi/2) cot|alpha|)),
hence    psi/rho <= -kappa - c rho sin(2 delta),   kappa = (pi/2) cos(gamma) - max(ln 2, -lambda_alpha) |sin(gamma)|.
The script (1) tests the exact formula and the bounds at random admissible parameters, (2) tabulates kappa for delta = 10 degrees,
(3) compares, for delta = 10 degrees, the actual sup of psi/rho with -kappa.   Floating point (numpy/scipy).  Usage: python3 w1_descent.py"""
import numpy as np
from scipy.special import spence
from scipy.integrate import quad
Li2 = lambda z: spence(1 - z)

def G(w):
    w = np.atleast_1d(np.asarray(w, dtype=complex)); out = np.empty_like(w); m = w.real <= 0
    out[m] = -w[m]**2/2 + Li2(-np.exp(w[m])) + np.pi**2/12
    out[~m] = -w[~m]**2 - Li2(-np.exp(-w[~m])) - np.pi**2/12
    return out

def ell(w, sign):
    f = lambda lam: np.log(1 + np.exp(-sign*lam*w))
    return quad(lambda l: f(l).real, 0, 1, limit=400)[0], quad(lambda l: f(l).imag, 0, 1, limit=400)[0]

def lam_alpha(al):
    return np.log(1 - np.exp(-(np.pi/2)/np.tan(abs(al)))) if abs(al) > 1e-9 else 0.0

def kappa(beta, d):
    al, ga = np.radians(beta - d), np.radians(beta + d)
    return (np.pi/2)*np.cos(ga) - max(np.log(2), -lam_alpha(al))*abs(np.sin(ga))

if __name__ == "__main__":
    rng = np.random.default_rng(1); dev = 0; maxX = -9; minX = 9; maxY = 0; viol = 0; n = 0
    for _ in range(1500):
        beta = rng.uniform(-85, 85); d = rng.uniform(0.5, 90 - abs(beta))
        al, ga = np.radians(beta - d), np.radians(beta + d); S = np.sin(2*np.radians(d)); rho = rng.uniform(0.02, 30)
        for sgn in (1, -1):
            w = sgn*rho*np.exp(1j*al)
            psi = (np.exp(-1j*np.radians(2*beta - 90))*G(w)[0]).real - np.pi*rho*np.cos(ga)
            X, Y = ell(w, sgn)
            form = (-(np.pi - Y)*np.cos(ga) - X*np.sin(ga) - rho*S) if sgn == 1 else (-(np.pi + Y)*np.cos(ga) + X*np.sin(ga) - rho*S/2)
            dev = max(dev, abs(psi/rho - form)); maxX = max(maxX, X - np.log(2)); minX = min(minX, X - lam_alpha(al)); maxY = max(maxY, abs(Y)); n += 1
            if psi/rho > -kappa(beta, d) - (1 if sgn == 1 else 0.5)*rho*S + 1e-9: viol += 1
    print(f"(1) {n} random points (beta in (-85,85), all admissible delta, rho in (0.02,30), both halves):")
    print(f"    max |psi/rho - exact formula| = {dev:.2e}")
    print(f"    max (X - ln 2) = {maxX:.3e} (<= 0),   min (X - lambda_alpha) = {minX:.3e} (>= 0),   max |Y| = {maxY:.4f} (<= pi/2 = 1.5708)")
    print(f"    violations of  psi/rho <= -kappa - c rho sin(2 delta):  {viol}")
    print("(2) kappa for delta = 10 degrees:  beta, theta = arg t, alpha, gamma, kappa")
    worst = 9
    for beta in range(-64, 66, 2):
        k = kappa(beta, 10)
        if -60 <= beta <= 50: worst = min(worst, k)
        if beta % 10 == 0 or beta in (-64, 64): print(f"    beta = {beta:4d}   theta = {2*beta - 90:5d}   alpha = {beta - 10:4d}   gamma = {beta + 10:4d}   kappa = {k:7.4f}")
    fine = min(kappa(b, 10) for b in np.linspace(-60, 50, 2201))
    print(f"    minimum of kappa over -60 <= beta <= 50 (theta from -210 to 10), step 0.05: {fine:.4f}")
    print("(3) delta = 10 degrees: actual sup over rho (both halves) of psi/rho, against -kappa")
    rho = np.concatenate([np.linspace(0.002, 1, 300), np.linspace(1.01, 12, 500), np.linspace(12.1, 60, 300)])
    for beta in (-60, -45, -30, -15, 0, 15, 30, 42.5, 50):
        al, ga = np.radians(beta - 10), np.radians(beta + 10); best = -9
        for sgn in (1, -1):
            w = sgn*rho*np.exp(1j*al)
            psi = (np.exp(-1j*np.radians(2*beta - 90))*G(w)).real - np.pi*rho*np.cos(ga)
            best = max(best, np.max(psi/rho))
        print(f"    beta = {beta:5.1f}  theta = {2*beta - 90:6.1f}:   sup psi/rho = {best:7.4f}    -kappa = {-kappa(beta, 10):7.4f}")
