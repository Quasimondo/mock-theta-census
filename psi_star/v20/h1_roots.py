"""Addendum v20: psi1* and psi0* at even roots of unity q0 = e^(2 pi i/X).
 a  Theorem 9: the finite sum  sum_{n < X/2} (-1)^n q0^(n(n+1)/2) (-q0;q0)_n  against  P J(1/(2b)) + P I(1/(2b)) - P I(-1/(2b))  at x = 1/X
    (real x: Phi_b from Faddeev's integral, double precision; the three integrals over the line L with delta = 10 degrees);
    the same for psi0* = sum (-1)^n q0^((n+1)(n+2)/2) (-q0;q0)_n  with the integrals at u = b + 1/(2b), b - 1/(2b) and the factor q^(1/2)
 b  radial limits: psi1*(rho q0) -> finite sum (Lemma R1),  W(p), p W3(p) -> 1 along the induced path (Lemma R2)
 c  Theorem 10: psi1*(q0) / [ 3^(-1/4) (X/2)^(1/2) e^(Cl2 X/2pi) e^(2 pi i/3 - i pi X/6) ] -> 1,  and the same for psi0* with the extra factor e^(2 pi i/3)
Usage: python3 h1_roots.py [a b c]"""
import sys
import numpy as np
from mpmath import mp, mpf, mpc, exp, pi, sqrt, nstr, clsin, expjpi
from fad import logPhi_np
I = mpc(0, 1)

def fin(X, case="psi1"):                       # exact finite sum, mpmath
    with mp.workdps(30 + int(0.0702*X) + 10):
        q = expjpi(mpf(2)/X); s = mpc(0); term = mpc(1); qn = mpc(1)
        for n in range(X//2):
            s += term*(q**(n + 1) if case == "psi0" else 1); qn *= q; term *= -qn*(1 + qn)
        return +s

def PI_np(X, kernel, extra_w=False, R=None, panels=400, order=12):
    """e^(t/8) int C0 e^(t s^2/2) K(s)/Phi_b ds over s = r e^(-i 10deg), x = 1/X;  kernel: 'J' (1/(2 cosh pi s)), '+' (e^(pi s)), '-' (e^(-pi s));
    extra_w: extra factor e^w e^(-t/2)  (the psi0 versions, with q^(1/2) included)."""
    x = 1.0/X; t = -2j*np.pi*x; ga = np.radians(10.0)
    R = R or 6.0*np.sqrt(X)
    xs, ws = np.polynomial.legendre.leggauss(order); tot = 0j
    edges = np.linspace(-R, R, panels + 1)
    for a, b in zip(edges[:-1], edges[1:]):
        for xi, wi in zip(xs, ws):
            r = (a + b)/2 + (b - a)/2*xi; s = r*np.exp(-1j*ga); w = 1j*t*s
            lg = np.pi**2/(12*t) + t/24 - 0.5*np.log(2) + t*s*s/2 - logPhi_np(w, x)
            if kernel == "J":
                z = np.pi*s*(1 if s.real >= 0 else -1); lg -= z + np.log1p(np.exp(-2*z))
            else:
                lg += np.pi*s*(1 if kernel == "+" else -1)
            if extra_w: lg += w - t/2
            tot += wi*(b - a)/2*np.exp(lg)*np.exp(-1j*ga)
    return tot*np.exp(t/8)

def PI_w(X, sign, extra_w=False, order=12):
    """P I(+-1/(2b)) at x = 1/X on the descent contours of v19 (no cancellation): '+': the line Im w = -2 pi/3;
    '-': ray i pi - rho e^(-i pi/6), segment [i pi, 0], positive real axis.   s = w/(i t), ds = dw/(i t)."""
    x = 1.0/X; t = -2j*np.pi*x; T = abs(t)
    xs, ws = np.polynomial.legendre.leggauss(order)
    def f(w):
        sv = w/(1j*t)
        lg = np.pi**2/(12*t) + t/24 - 0.5*np.log(2) + t*sv*sv/2 - logPhi_np(w, x) + np.pi*sv*sign
        if extra_w: lg += w - t/2
        return np.exp(lg)
    def seg(wfun, dfun, a, b, panels):
        tot = 0j; edges = np.linspace(a, b, panels + 1)
        for u, v in zip(edges[:-1], edges[1:]):
            for xi, wi in zip(xs, ws):
                p = (u + v)/2 + (v - u)/2*xi
                tot += wi*(v - u)/2*f(wfun(p))*dfun
        return tot
    if sign > 0:
        R = 14*T + 2
        if extra_w and np.pi/(3*T) > 1: R = max(R, 36/(np.pi/(3*T) - 1))      # the extra e^w slows the decay on the right
        tot = seg(lambda r: r - 2j*np.pi/3, 1.0, -R, R, 1200)
    else:
        e = np.exp(-1j*np.pi/6)
        tot = seg(lambda rho: 1j*np.pi - rho*e, e, 0.0, 7*np.sqrt(T) + 1, 300)
        tot += seg(lambda eta: 1j*eta, -1j, 0.0, np.pi, 300)
        tot += seg(lambda r: r + 0j, 1.0, 0.0, 14*T + 1, 300)
    return tot*np.exp(t/8)/(1j*t)

def part_a():
    print("(a) Theorem 9 at x = 1/X: finite sum against P J + P I(+) - P I(-)   (double precision; J over L, the I over the contours of v19)")
    for X in (4, 6, 8, 12, 16, 24, 40, 60):
        for case in ("psi1", "psi0"):
            ew = case == "psi0"
            J = PI_np(X, "J", ew)
            # psi0, small X: on the horizontal line the extra factor e^w beats the decay e^(-pi r/(3|t|)); use L there (the integral is small, no cancellation)
            Ip = PI_np(X, "+", True) if (ew and X <= 8) else PI_w(X, +1, ew); Im = PI_w(X, -1, ew)
            f = complex(fin(X, case)); rhs = J + Ip - Im
            print(f"    X = {X:2d} {case}*: sum = {f:.10f},  |sum - formula|/|sum| = {abs(f - rhs)/abs(f):.1e}   (P J = {J:.6f}, P I(+) = {Ip:.3e}, P I(-) = {Im:.3e})")

def part_b():
    mp.dps = 40
    print("(b) radial limits at X = 8:  q = q0 e^(-2 pi eps), x = 1/X + i eps")
    X = 8; f = fin(X)
    for eps in ("1e-2", "1e-3", "1e-4"):
        e = mpf(eps); x = mpf(1)/X + I*e; q = exp(2*pi*I*x); p = exp(-pi*I/x)
        s = mpc(0); term = mpc(1); n = 0
        while abs(term) > mpf(10)**(-45) or n < 4*X:
            s += term; n += 1; term *= -q**n*(1 + q**n)
        W = W3 = mpc(0); po = mpc(1); j = 0
        while abs(p)**(j*j) > mpf(10)**(-45):
            W += (-1)**j*p**(j*j + j)*po; W3 += (-1)**j*p**(j*j + 3*j)*po; po *= 1 - p**(2*j + 1); j += 1
        print(f"    eps = {eps}: |psi1*(q) - finite sum| = {nstr(abs(s - f), 3)},  |W(p) - 1| = {nstr(abs(W - 1), 3)},  |p W3(p) - 1| = {nstr(abs(p*W3 - 1), 3)},  |1 - p|/(1 - |p|) = {nstr(abs(1 - p)/(1 - abs(p)), 6)}")

def part_c():
    mp.dps = 30
    Cl2 = clsin(2, pi/3)
    print("(c) Theorem 10: exact finite sum divided by the asymptotic value, minus 1, times X")
    for X in (100, 200, 400, 800, 1600, 3200):
        lead = mpf(3)**(-mpf(1)/4)*sqrt(mpf(X)/2)*exp(Cl2*X/(2*pi))*exp(2*pi*I/3 - I*pi*X/6)
        r1 = fin(X, "psi1")/lead - 1; r0 = fin(X, "psi0")/(lead*exp(2*pi*I/3)) - 1
        print(f"    X = {X:5d}: psi1*: X (ratio - 1) = {nstr(X*r1, 6)};   psi0*: X (ratio - 1) = {nstr(X*r0, 6)}")

if __name__ == "__main__":
    for p_ in (sys.argv[1:] or list("abc")): {"a": part_a, "b": part_b, "c": part_c}[p_]()
