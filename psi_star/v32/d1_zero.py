"""Addendum v32: Omega_gamma(a, x) = S(a, x+1)^-1 S^c(a', x'-1) on all of C minus (-inf, -1/2], in particular at x = 0.
   (a) Im x > 0: against the series J^-1 Lambda^-1 J'  (mpmath quadrature at 30 digits, and double precision)
   (b) the constant: i delta(a, x+1)/delta(a', x'-1) = lambda(a, x)
   (c) Im x < 0 and the interval (-1/2, 0): against Omega^(1) and Omega^(2) of v31
   (d) around x = 0: the value at 0, a circle around 0 (series above, Omega^(1) below, Omega^(3) everywhere), Fourier modes of Omega^(3)
Usage: python3 d1_zero.py [a b c d]"""
import sys
import numpy as np
from mpmath import mp, mpf, mpc, exp, pi, sqrt, nstr, matrix
from common import Setup, I
import frames
from smat import *
np.set_printoptions(linewidth=220)
parts = sys.argv[1:] or ["a", "b", "c", "d"]

if "a" in parts:
    mp.dps = 30
    print("(a) Im x > 0:  max|S(a,x+1)^-1 S^c(a',x'-1) - J(a;x)^-1 Lambda^-1 J(a';x')|")
    for x, a in ((mpc("0.3", "0.2"), mpf(1)/4), (mpc("-0.15", "0.1"), mpc("0.3", "0.1"))):
        xp = x/(1 + 2*x); ap = a/(1 + 2*x)
        S1 = frames.state_matrix(Setup(x + 1), a)
        X2 = xp - 1
        SC = frames.state_matrix(Setup(-X2.conjugate()), -mpc(ap).conjugate()).conjugate()
        O3 = S1**-1*SC
        l = exp(I*pi/12)*exp(-5*pi*I*x*xp/6)*exp(-pi*I*(a - ap) - pi*I*a*ap)/sqrt(1 + 2*x)
        Os = frames.frame_q(Setup(x), a)**-1*matrix([[1, 0, 0], [0, 1/l, 0], [0, 0, -1/l]])*frames.frame_q(Setup(xp), ap)
        print(f"   x = {nstr(x, 5)}, a = {nstr(a, 4)}  (30 digits): {nstr(frames.mnorm(O3 - Os), 3)}   (max entry {nstr(frames.mnorm(Os), 3)})")
    for x, a in ((0.05 + 0.05j, 0.3), (-0.05 + 0.03j, 0.3), (0.4 + 0.5j, 0.3 + 0.1j), (-0.3 + 0.2j, 0.25)):
        print(f"   x = {x}, a = {a}  (double precision): {np.abs(Om3(a, x) - Om_series(a, x)).max():.2e}")

if "b" in parts:
    mp.dps = 30
    print("(b) |i delta(a, x+1)/delta(a', x'-1) / lambda(a, x) - 1|")
    for x, a in ((mpc("0.3", "0.2"), mpf(1)/4), (mpc("-0.4", "0.1"), mpc("0.3", "0.1")), (mpc("-2.3", "0.7"), mpf("1.2"))):
        xp = x/(1 + 2*x); ap = a/(1 + 2*x)
        l = exp(I*pi/12)*exp(-5*pi*I*x*xp/6)*exp(-pi*I*(a - ap) - pi*I*a*ap)/sqrt(1 + 2*x)
        d = frames.delta(Setup(x + 1), a)/frames.delta(Setup(xp - 1), ap)
        print(f"   x = {nstr(x, 5)}, a = {nstr(a, 4)}: {nstr(abs(I*d/l - 1), 3)}")

if "c" in parts:
    a = 0.3
    print("(c) max|Omega^(3) - Omega^(1)| for Im x < 0, a = 0.3")
    for x in (0.5 - 0.15j, 0.05 - 0.05j, -0.2 - 0.3j, -0.25 - 0.05j, -0.45 - 0.05j, -0.7 - 0.3j):
        print(f"   x = {x}: {np.abs(Om3(a, x) - Om1(a, x)).max():.2e}")
    print("    max|Omega^(3) - Omega^(2)| on the interval (-1/2, 0) and max|Omega^(3) - Omega^(1)| on the positive real axis")
    for x in (-0.1, -0.2, -0.3, -0.4):
        print(f"   x = {x}: {np.abs(Om3(a, x) - Om2(a, x)).max():.2e}")
    for x in (0.2, 0.5, 1.0):
        print(f"   x = {x}: {np.abs(Om3(a, x) - Om1(a, x)).max():.2e}")

if "d" in parts:
    a = 0.3
    O0 = Om3(a, 0.0)
    print("(d) x = 0, a = 0.3:  Omega_gamma(a, 0) = S(a, 1)^-1 conj S(-conj a, 1) =")
    print(np.array2string(O0, precision=10))
    for r, n in ((0.05, 16), (0.1, 16), (0.2, 32)):
        th = 2*np.pi*(np.arange(n) + 0.5)/n; du = dl = 0.0; vals = []
        for t in th:
            x = r*np.exp(1j*t); O = Om3(a, x); vals.append(O)
            if np.sin(t) > 0: du = max(du, np.abs(O - Om_series(a, x)).max())
            else: dl = max(dl, np.abs(O - Om1(a, x)).max() if x.real > 0 else np.abs(O - Om2(a, x)).max())
        vals = np.array(vals)
        neg = [np.abs((vals*np.exp(1j*k*th)[:, None, None]).mean(axis=0)).max() for k in (1, 2, 3)]
        print(f"   circle |x| = {r}, {n} points: max|Omega^(3) - series| above = {du:.2e};  max|Omega^(3) - Omega^(1) or Omega^(2)| below = {dl:.2e} (the old formulas lose accuracy as |x| -> 0);"
              f"  |mean - Omega(0)| = {np.abs(vals.mean(axis=0) - O0).max():.2e};  negative Fourier modes k = 1, 2, 3: " + ", ".join(f"{v:.1e}" for v in neg))
