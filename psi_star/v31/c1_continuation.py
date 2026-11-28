"""Addendum v31: the cocycle Omega_gamma continues across (-1/2, 0).
   S^c(a, x) := conj S(-conj a, -conj x)   (holomorphic for x outside [0, inf));   R(a, x) := S^c(a, x) S(a, x)^-1.
   (a) Im x > 0:  Jt(-a; x) S^c(a, x) = Delta(a, x) J(a; x)                                     (mpmath, 30 digits)
   (b) dual shift: S(a+1, x) = C S(a, x),  C = [[1, -1, 0], [0, p, -p/T], [0, -1, 0]]           (Im x > 0: mpmath; real x and Im x < 0: numpy)
   (c) Im x < 0:  R(a', x') = R(a, x)  for x' = x/(1+2x), a' = a/(1+2x)                         (numpy; this is Theorem 22, step 3)
   (d) Omega^(2) = S^c(a,x)^-1 S^c(a',x') against the series above the axis, against Omega^(1) = S(a,x)^-1 S(a',x') below, and on the interval itself
   (e) the algebra of Lemma U: the two scalar equations for (g, u) at random Laurent polynomials; for |q| < 1 an explicit V != 1 (so |q| > 1 is needed)
Usage: python3 c1_continuation.py [a b c d e]    (a: 6 minutes at 30 digits; b, c, d: 3 minutes in double precision; e: seconds)"""
import sys
import numpy as np
from mpmath import mp, mpf, mpc, exp, pi, nstr, matrix
from common import Setup, I
import frames
from npint import NP
from lowint import Low
np.set_printoptions(linewidth=200)
parts = sys.argv[1:] or ["a", "b", "c", "d", "e"]

def tonp(M): return np.array([[complex(M[i, j]) for j in range(M.cols)] for i in range(M.rows)])

def Cmat(a, x):
    p = np.exp(-1j*np.pi/x); T = np.exp(1j*np.pi*(2*a - 1)/x)
    return np.array([[1, -1, 0], [0, p, -p/T], [0, -1, 0]])

def Sraw(a, x):
    um = max(abs(a/np.sqrt(x)), abs((a - 1)/np.sqrt(x)), abs((a + x)/np.sqrt(x)), abs((a + x - 1)/np.sqrt(x))) + 0.2
    if x.imag < 0: return Low(x, umax=um).Smat(a)
    d = min(0.75, 0.9*(np.pi/2 - abs(np.angle(np.sqrt(x)))))          # a steeper line: less cancellation
    g = np.pi*np.sin(2*d); R = (2*np.pi*um + np.sqrt((2*np.pi*um)**2 + 4*g*36.0))/(2*g)
    return NP(x, delta=d, R=R).Smat(a)

def Sany(a, x):
    """S(a, x) in double precision for any x off the negative real axis (Im x >= 0 needs Re x > 0).
       Re a is first reduced to [0, 1) with the dual shift S(a+1) = C(a) S(a), to keep the integrands small."""
    x = complex(x); a = complex(a)
    k = int(np.floor(a.real)); a0 = a - k
    S = Sraw(a0, x)
    for j in range(k): S = Cmat(a0 + j, x) @ S
    for j in range(-k): S = np.linalg.solve(Cmat(a0 - j - 1, x), S)
    return S

def Sc(a, x):
    return np.conj(Sany(-np.conj(a), -np.conj(x)))

def gam(a, x): return a/(1 + 2*x), x/(1 + 2*x)

def lam_np(a, x):
    ap, xp = gam(a, x)
    return np.exp(1j*np.pi/12 - 5j*np.pi*x*xp/6 - 1j*np.pi*(a - ap) - 1j*np.pi*a*ap)/np.sqrt(1 + 2*x)

def Om_series(a, x):
    """Im x > 0: J(a;x)^-1 Lambda^-1 J(a';x') (mpmath series)"""
    ap, xp = gam(a, x); l = lam_np(a, x)
    J = tonp(frames.frame_q(Setup(mpc(x)), mpc(a))); Jp = tonp(frames.frame_q(Setup(mpc(xp)), mpc(ap)))
    return np.linalg.solve(J, np.diag([1, 1/l, -1/l]) @ Jp)

if "a" in parts:
    mp.dps = 30
    print("(a) reflected Theorem 20, Im x > 0:  max|Jt(-a;x) S^c(a,x) - Delta(a,x) J(a;x)|,  S^c(a,x) = conj S(-conj a, -conj x) by quadrature")
    for x, a in ((mpc("0.3", "0.2"), mpf(1)/4), (mpc("-0.2", "0.6"), mpc("0.3", "0.1")), (mpc("-0.35", "0.25"), mpf("0.8"))):
        S = Setup(x); Sm = Setup(-x.conjugate())
        SC = frames.state_matrix(Sm, -a.conjugate() if isinstance(a, mpc) else -a).conjugate()
        dev = frames.mnorm(frames.frame_dual(S, -a)*SC - frames.Delta(S, a)*frames.frame_q(S, a))
        print(f"   x = {nstr(x, 5)}, a = {nstr(a, 4)}: {nstr(dev, 3)}   (max entry {nstr(frames.mnorm(frames.frame_q(S, a)), 3)})")

if "b" in parts:
    print("(b) dual shift  S(a+1, x) = C S(a, x):  max deviation (integrals at a and at a+1, no reduction)")
    for x, a in ((0.3 + 0.2j, 0.25), (0.6 + 0.0j, 0.3 + 0.1j), (0.5 - 0.15j, 0.3), (-0.2 - 0.3j, 0.3 + 0.1j)):
        print(f"   x = {x}, a = {a}: {np.abs(Sraw(a + 1, x) - Cmat(a, x) @ Sraw(a, x)).max():.2e}")

if "c" in parts:
    print("(c) Im x < 0:  R(a, x) = S^c(a, x) S(a, x)^-1 at (a, x) and at (a', x'): max|R' - R|, first row | block   (sizes)")
    for x, a in ((0.5 - 0.15j, 0.3), (0.3 - 0.3j, 0.3 + 0.1j), (-0.2 - 0.3j, 0.3), (-0.3 - 0.1j, 0.25), (-0.7 - 0.2j, 0.3), (0.2 - 0.6j, 0.8)):
        ap, xp = gam(a, x)
        R = Sc(a, x) @ np.linalg.inv(Sany(a, x)); Rp = Sc(ap, xp) @ np.linalg.inv(Sany(ap, xp))
        D = np.abs(Rp - R)
        print(f"   x = {x}, a = {a}: {D[0].max():.2e} | {D[1:].max():.2e}   ({np.abs(R[0]).max():.3g} | {np.abs(R[1:]).max():.3g})")

if "d" in parts:
    a = 0.3
    def Om2(x):
        ap, xp = gam(a, x); return np.linalg.solve(Sc(a, x), Sc(ap, xp))
    def Om1(x):
        ap, xp = gam(a, x); return np.linalg.solve(Sany(a, x), Sany(ap, xp))
    print("(d) Omega^(2) = S^c(a,x)^-1 S^c(a',x'), a = 0.3")
    print("    above the axis (Re x < 0) against the series J^-1 Lambda^-1 J':")
    for x in (-0.25 + 0.05j, -0.15 + 0.02j, -0.4 + 0.1j):
        print(f"      x = {x}: max|Omega^(2) - series| = {np.abs(Om2(x) - Om_series(a, x)).max():.2e}")
    print("    below the axis against Omega^(1) = S(a,x)^-1 S(a',x'):")
    for x in (-0.25 - 0.05j, -0.15 - 0.03j, -0.4 - 0.1j, -0.3 - 0.4j, -0.45 - 0.05j):
        print(f"      x = {x}: max|Omega^(2) - Omega^(1)| = {np.abs(Om2(x) - Om1(x)).max():.2e}")
    print("    on the interval: Omega_gamma(x0) from state integrals at the positive moduli -x0 and -x0/(1+2 x0); the series at x0 + i eps and Omega^(1) at x0 - i eps tend to it")
    for x0 in (-0.1, -0.2, -0.3):
        O0 = Om2(x0); up = []; lo = []
        for eps in (0.02, 0.01, 0.005):
            up.append(f"{np.abs(Om_series(a, x0 + 1j*eps) - O0).max():.3g}")
        for eps in (0.04, 0.02):
            lo.append(f"{np.abs(Om1(x0 - 1j*eps) - O0).max():.3g}")
        print(f"      x0 = {x0}: first row of Omega_gamma(x0) = {np.array2string(O0[0], precision=8)}")
        print(f"            max|series(x0 + i eps) - Omega_gamma(x0)|, eps = 0.02, 0.01, 0.005: " + ", ".join(up) + ";   max|Omega^(1)(x0 - i eps) - Omega_gamma(x0)|, eps = 0.04, 0.02: " + ", ".join(lo))
    x0, r, n = -0.25, 0.1, 16
    th = 2*np.pi*(np.arange(n) + 0.5)/n; du = dl = 0.0
    for t in th:
        x = x0 + r*np.exp(1j*t); O = Om2(x)
        if np.sin(t) > 0: du = max(du, np.abs(O - Om_series(a, x)).max())
        else: dl = max(dl, np.abs(O - Om1(x)).max())
    print(f"    circle |x + 0.25| = 0.1, {n} points: max|Omega^(2) - series| on the upper half = {du:.2e};  max|Omega^(2) - Omega^(1)| on the lower half = {dl:.2e}")

if "e" in parts:
    print("(e) Lemma U")
    rng = np.random.default_rng(1)
    q = 1.3*np.exp(0.7j)
    gc = rng.normal(size=7) + 1j*rng.normal(size=7); uc = rng.normal(size=7) + 1j*rng.normal(size=7)
    lau = lambda c, z: sum(ck*z**k for k, ck in zip(range(-3, 4), c))
    g = lambda z: lau(gc, z); u = lambda z: lau(uc, z)
    r = lambda z: z*q*u(z/q); sf = lambda z: g(z*q) - z*q*u(z)
    A = lambda z: np.array([[0, 1/(z*q*q)], [1, -1/q]])
    V = lambda z: np.array([[g(z), u(z)], [r(z), sf(z)]])
    dev = 0.0
    for z in (0.8 + 0.3j, -1.1 + 0.6j, 0.2 - 1.4j):
        M = A(z) @ V(z*q) - V(z) @ A(z)
        e2 = g(z) - g(z*q*q) - z*q*u(z) + z*q*q*u(z*q)                                   # equation (2')
        e4 = u(z/q) - q*u(z*q) - (g(z*q) - g(z*q*q) - z*q*u(z) + z*q*q*u(z*q))             # equation (4')
        dev = max(dev, abs(M[0, 0]), abs(M[1, 0]), abs(M[0, 1] + e2/(z*q*q)), abs(M[1, 1] + e4/q))
    print(f"   with r(z) = zq u(z/q), s(z) = g(zq) - zq u(z):  A V(zq) - V A = [[0, -(2')/(zq^2)], [0, -(4')/q]],  max deviation at random g, u: {dev:.1e}")
    q = 0.5*np.exp(0.9j)
    def F(z):
        s = 0; t = 1.0 + 0j
        for n in range(1, 80):
            s += t; t *= q**n*z*(1 + q**n)
        return s
    A3 = lambda z: np.array([[1, 0, -1/(z*q*q)], [0, 0, 1/(z*q*q)], [0, 1, -1/q]])
    V3 = lambda z: np.array([[1, -F(z), -F(z*q)], [0, 2, 0], [0, 0, 2]])
    dev = max(np.abs(A3(z) @ V3(z*q) - V3(z) @ A3(z)).max() for z in (0.8 + 0.3j, -1.1 + 0.6j))
    print(f"   |q| < 1: V = [[1, -F(z), -F(zq)], [0, 2, 0], [0, 0, 2]] satisfies A(z) V(zq) = V(z) A(z):  max deviation {dev:.1e}")
