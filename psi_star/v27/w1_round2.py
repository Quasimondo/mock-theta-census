"""Addendum v27: checks for the remarks after the second stress test.
 a  Theorem 15 with kappa = X Vol(a)/2 pi as the parameter (even X):
      |F/Pi - 1| <= C(a)/X + T_J + T_dec (1 + O(1/X)),
      T_J = e^(-kappa) (|1+2W|/|1+W|)^(1/2) / (2 X^(1/2) sin(pi a)),   T_dec = |W| e^(-2 kappa),
    printed: kappa, X |F/Pi - 1|, X T_J, X T_dec, and |F/Pi - 1| / (1.7/X + T_J + 1.1 T_dec)   (1.7: above the measured values of X |F/Pi - 1| once kappa is large);
    also |Pi| against (X/2)^(1/2) |(1+W)/(1+2W)|^(1/2) e^kappa, and |P I(a/b)| against (X/2)^(1/2) (|W|/|1+2W|)^(1/2) e^(-kappa).
 b  Theorem 13 (b): Vol(a), the angle eps(a) measured by the stress test (round 2, T13b.4), eps/Vol, and the tilt allowed by the proof.
 c  Lemma J for x up to 1: |P J(a/b)| against 1/(2 sqrt 2 sin(pi a)), integral over the real line.
Usage: python3 w1_round2.py [a b c]"""
import sys
import numpy as np
from mpmath import mp, mpf, mpc, exp, pi, sqrt, log, polylog, nstr, expjpi, im, fabs, sin

def saddle(a):
    """upper root W, w0 = log W, V(w0) for the growing integral, Vol = -Im V(w0)"""
    a = mpf(a); z = exp(2j*pi*a); d = sqrt(1 + 4/z); W = (-1 + d)/2
    if im(W) < 0: W = (-1 - d)/2
    w0 = log(W); V = -w0*w0/2 - 2j*pi*(a - 1)*w0 + polylog(2, -W)
    return W, w0, V, -im(V)

def F_root(a, X):
    with mp.workdps(30 + int(0.0702*X) + 10):
        z = expjpi(2*mpf(a)); q = expjpi(mpf(2)/X); s = mpc(0); term = mpc(1); qn = mpc(1)
        for n in range(X//2):
            s += term; qn *= q; term *= qn*z*(1 + qn)
        return +s

def part_a():
    mp.dps = 30
    print("(a) Theorem 15 with kappa = X Vol(a)/2 pi")
    print("      a        X     kappa    X|F/Pi-1|     X T_J      X T_dec    |F/Pi-1|/(1.7/X + T_J + 1.1 T_dec)   | |Pi|/closed form - 1 |")
    worst = 0
    for a, Xs in (("0.001", (800, 3200, 6400)), ("0.01", (100, 400, 1600, 3200)), ("0.05", (50, 200, 800)), ("0.3", (10, 20, 50, 200, 800)),
                  ("0.5", (10, 20, 50, 200)), ("0.9", (50, 200, 800)), ("0.99", (400, 1600, 3200)), ("0.999", (3200, 6400))):
        W, w0, V, vol = saddle(a)
        for X in Xs:
            t = -2j*pi/X; kap = X*vol/(2*pi)
            Pi = (1 + W)*exp(w0/2)*sqrt(pi/(t*(1 + 2*W)))*exp((V + pi**2/12)/t)
            lhs = fabs(F_root(a, X)/Pi - 1)
            TJ = exp(-kap)*sqrt(fabs(1 + 2*W)/fabs(1 + W))/(2*sqrt(X)*sin(pi*mpf(a))); Td = fabs(W)*exp(-2*kap)
            ratio = lhs/(mpf("1.7")/X + TJ + mpf("1.1")*Td); worst = max(worst, ratio)
            mod = fabs(fabs(Pi)/(sqrt(mpf(X)/2)*sqrt(fabs((1 + W)/(1 + 2*W)))*exp(kap)) - 1)
            print(f"    {a:>6s} {X:6d} {nstr(kap, 5):>9s} {nstr(X*lhs, 5):>11s} {nstr(X*TJ, 5):>11s} {nstr(X*Td, 5):>11s} {nstr(ratio, 4):>12s} {nstr(mod, 2):>34s}")
    print(f"    largest ratio: {nstr(worst, 4)}")
    from n1_zroots import PI_decaying
    print("    the decaying term: |P I(a/b)| / [(X/2)^(1/2) (|W|/|1+2W|)^(1/2) e^(-kappa)] - 1")
    for a in (0.1, 0.3, 0.5, 0.8):
        W, w0, V, vol = saddle(a); row = []
        for X in (16, 32, 64):
            kap = float(X*vol/(2*pi)); val = abs(PI_decaying(a, X))
            row.append(f"{val/(np.sqrt(X/2)*np.sqrt(float(fabs(W)/fabs(1 + 2*W)))*np.exp(-kap)) - 1:+.4f}")
        print(f"      a = {a}: " + ", ".join(row) + "   (X = 16, 32, 64)")

EPS = {"0.001": 0.08, "0.005": 0.39, "0.01": 0.78, "0.02": 1.55, "0.05": 3.9, "0.1": 7.8, "0.2": 16.2, "0.3": 20.1, "0.5": 17.2,
       "0.7": 9.6, "0.8": 6.2, "0.9": 3.0, "0.95": 1.4, "0.99": 0.25, "0.999": 0.02}       # degrees; docs/psi-star-stress-test-v2.md, T13b.4

def part_b():
    mp.dps = 20
    print("(b) Theorem 13 (b): eps(a) as measured by the stress test, against the volume")
    print("      a      Vol(a)    eps (deg)   eps/Vol (deg)   tilt allowed by the proof, arctan(Vol/(3 pi^2)) (deg)")
    lo = 1e9
    for a, e in EPS.items():
        vol = float(saddle(a)[3]); lo = min(lo, e/vol)
        print(f"    {a:>6s} {vol:9.5f} {e:9.2f} {e/vol:12.1f} {np.degrees(np.arctan(vol/(3*np.pi**2))):22.3f}")
    print(f"    smallest eps/Vol: {lo:.1f} degrees")

def part_c():
    from n1_zroots import integrand, seg
    print("(c) Lemma J for x up to 1: |P J(a/b)| 2 sqrt 2 sin(pi a)  (must be <= 1), integral over the real line")
    worst = 0
    for x in (0.3, 1/3, 0.5, 0.75, 1.0):
        X = 1/x; t = -2j*np.pi*x; T = abs(t); row = []
        for a in (0.05, 0.2, 0.284, 0.5, 0.8, 0.95):
            R = 16*T/(2*np.pi*min(a, 1 - a)) + 2*T
            J = seg(lambda w: integrand(w, X, a, kernel=True), lambda r: r + 0j, 1.0, -R, R, 800)*np.exp(t/8)/(1j*t)
            v = abs(J)*2*np.sqrt(2)*np.sin(np.pi*a); worst = max(worst, v); row.append(f"{v:.4f}")
        print(f"    x = {x:.4f}: " + ", ".join(row) + "   (a = 0.05, 0.2, 0.284, 0.5, 0.8, 0.95)")
    print(f"    largest: {worst:.4f}")

if __name__ == "__main__":
    for p in (sys.argv[1:] or list("abc")):
        {"a": part_a, "b": part_b, "c": part_c}[p]()
