"""Addendum v25: the formulas at roots of unity for X = 1, 2, 3 (x = 1/X; b = 1, 2^(-1/2), 3^(-1/2)).
   odd X (Theorem 11):   A  = P [3 J - I(+) - I(-)],   radial limit of psi1* is A/3;   psi0*: A0 with u = b + 1/(2b), b - 1/(2b) and the factor q^(1/2)
   even X (Theorem 9):   finite sum = P [J + I(+) - I(-)]
   Integrals over the line L (10 degrees), Phi_b from Faddeev's integral (fad.py), double precision.   Also X = 4, 5 for comparison.
 b  Theorem 14 at X = 2:  1 = i P e^(-pi i a) [ J(a/b) - I((a-1)/b) + z^2 I(a/b) ]  at x = 1/2, for several a
 c  the radial limits themselves: psi1*(r) -> 1/3, psi0*(r) -> 1/3, psi1*(r w) -> 1, psi0*(r w) -> 0   (w = e^(2 pi i/3), r = e^(-2 pi eps))
Usage: python3 r1_small.py"""
import numpy as np
from h1_roots import PI_np, fin
from k1_odd import A_odd

print("(a) Theorems 9 and 11 at x = 1/X")
for X in (1, 2, 3, 4, 5):
    row = []
    for case in ("psi1", "psi0"):
        ew = case == "psi0"
        kw = dict(R=7.0*np.sqrt(X) + 3, panels=700)
        J = PI_np(X, "J", ew, **kw); Ip = PI_np(X, "+", ew, **kw); Im = PI_np(X, "-", ew, **kw)
        if X % 2:
            lhs = complex(A_odd(X, case)); rhs = 3*J - Ip - Im; name = "A " if case == "psi1" else "A0"
        else:
            lhs = complex(fin(X, case)); rhs = J + Ip - Im; name = "sum"
        row.append(f"{case}*: {name} = {lhs.real:+.6f}{lhs.imag:+.6f}i, formula - {name} = {abs(rhs - lhs):.1e}  (P J = {J:.6f}, P I(+) = {Ip:.6f}, P I(-) = {Im:.6f})")
    print(f"X = {X}  " + "\n       ".join(row))

from n1_zroots import integrand, seg
from mpmath import mp, mpf, mpc, exp, pi, nstr, expjpi
print("(b) Theorem 14 at X = 2 (the sum is 1)")
X = 2; t = -2j*np.pi/X; ga = np.radians(10.0); R = 12.0
def over_L(ap, kernel):
    return seg(lambda w: integrand(w, X, ap, kernel=kernel), lambda r: 1j*t*r*np.exp(-1j*ga), 1j*t*np.exp(-1j*ga), -R, R, 700)*np.exp(t/8)/(1j*t)
for a in (0.25, 1/3, 0.5, 0.7, 0.9):
    J = over_L(a, True); G = over_L(a - 1, False); D = over_L(a, False)
    rhs = 1j*np.exp(-1j*np.pi*a)*(J - G + np.exp(4j*np.pi*a)*D)
    print(f"    a = {a:.4f}: formula - 1 = {abs(rhs - 1):.1e}   (|P J| = {abs(J):.4f}, |P I((a-1)/b)| = {abs(G):.4f}, |P I(a/b)| = {abs(D):.4f})")
print("(c) radial limits")
mp.dps = 30
def psis(q):
    s1 = s0 = mpc(0); term = mpc(1); qn = mpc(1); n = 0
    while abs(term) > mpf(10)**(-25) or n < 20:
        s1 += term; s0 += term*q**(n + 1); n += 1; qn *= q; term *= -qn*(1 + qn)
    return s1, s0
for name, z0, lim1, lim0 in (("1", mpc(1), mpf(1)/3, mpf(1)/3), ("e^(2 pi i/3)", expjpi(mpf(2)/3), mpf(1), mpf(0))):
    for eps in ("0.01", "0.001", "0.0001"):
        with mp.workdps(30 + int(0.02/float(eps))):          # the terms reach 10^(0.0166/eps) before they cancel
            a, b = psis(z0*exp(-2*pi*mpf(eps)))
        print(f"    root {name:12s} eps = {eps:6s}: |psi1* - {nstr(lim1, 4)}| = {nstr(abs(a - lim1), 3)},  |psi0* - {nstr(lim0, 4)}| = {nstr(abs(b - lim0), 3)}")
