"""Addendum v23: the kernel integral at real x is bounded by an elementary argument, and the volume rule at even roots of unity follows.
 a  |Phi_b(y)| = 1 for real b and real y (Faddeev's integral): Re log Phi_b at random real points
 b  |P J(a/b)| at x = 1/X: integral over the real line against the integral over the tilted line L, and against the bound 1/(2 sqrt 2 sin(pi a))
 c  Theorem 15: F(z; e^(2 pi i/X)) / [ (1 + W) W^(1/2) (pi/(t (1 + 2W)))^(1/2) e^((V(w0) + pi^2/12)/t) ] - 1, times X, for even X
Usage: python3 p1_unitary.py [a b c]"""
import sys
import numpy as np
from fad import Lam_np
from n1_zroots import integrand, seg, PJ

def part_a():
    rng = np.random.default_rng(1); worst = 0
    for _ in range(200):
        x = 10**rng.uniform(-2.3, -0.6); w = rng.uniform(-3, 3)          # |Im w| = 0 < pi (1 + x): the integral itself
        worst = max(worst, abs(Lam_np(w + 0j, x).real))
    print(f"(a) max |Re log Phi_b| over 200 random real w in [-3, 3], x in [0.005, 0.25]: {worst:.1e}")

def part_b():
    print("(b) |P J(a/b)| at x = 1/X:  real line | line L (10 degrees) | difference | bound 1/(2 sqrt 2 sin(pi a))")
    for a in (0.05, 0.25, 1/3, 0.5, 0.7, 5/6, 0.95):
        for X in (8, 40, 200):
            T = 2*np.pi/X; t = -1j*T; R = 16*T/(2*np.pi*min(a, 1 - a)) + 2*T
            JR = seg(lambda w: integrand(w, X, a, kernel=True), lambda r: r + 0j, 1.0, -R, R, 600)*np.exp(t/8)/(1j*t)
            JL = PJ(a, X) if X <= 40 else None
            bound = 1/(2*np.sqrt(2)*np.sin(np.pi*a))
            print(f"    a = {a:.3f}, X = {X:3d}: {abs(JR):.6f} | " + (f"{abs(JL):.6f} | {abs(JR - JL):.1e}" if JL is not None else "   -     |    -   ") + f" | {bound:.6f}" + ("" if abs(JR) <= bound else "   <-- VIOLATED"))

def part_c():
    from mpmath import mp, mpf, mpc, exp, pi, sqrt, log, polylog, nstr, expjpi, im
    mp.dps = 30
    I = mpc(0, 1)
    def F_root(a, X):
        with mp.workdps(30 + int(0.0702*X) + 10):
            z = expjpi(2*a); q = expjpi(mpf(2)/X); s = mpc(0); term = mpc(1); qn = mpc(1)
            for n in range(X//2):
                s += term; qn *= q; term *= qn*z*(1 + qn)
            return +s
    print("(c) Theorem 15: X (F/prediction - 1) for even X;  the bound of (b) says |F - leading term| <= 1/(2 sqrt 2 sin(pi a)) + (exponentially small)")
    for a in (mpf(1)/4, mpf(1)/3, mpf(5)/6, mpf(3)/10, mpf(1)/20, sqrt(mpf(2))/2):
        z = exp(2*pi*I*a); d = sqrt(1 + 4/z); W = (-1 + d)/2
        if im(W) < 0: W = (-1 - d)/2
        w0 = log(W); V = -w0*w0/2 - 2*pi*I*(a - 1)*w0 + polylog(2, -W)
        row = []
        for X in (400, 800, 1600, 3200):
            t = -2*pi*I/X
            pred = (1 + W)*exp(w0/2)*sqrt(pi/(t*(1 + 2*W)))*exp((V + pi**2/12)/t)
            row.append(nstr(X*(F_root(a, X)/pred - 1), 6))
        print(f"    a = {nstr(a, 5)}: " + "   ".join(row))

if __name__ == "__main__":
    for p_ in (sys.argv[1:] or list("abc")): {"a": part_a, "b": part_b, "c": part_c}[p_]()
