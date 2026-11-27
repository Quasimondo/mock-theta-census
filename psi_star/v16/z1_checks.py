"""Addendum v16 (the wedge around the real x-axis), checks.   x = b^2 = i t/(2 pi), t = |t| e^(i theta); theta = -90 degrees is real x, |q| = 1.
Lam(w) = log Phi_b(w/(2 pi b)) by Faddeev's integral (fad.py).
 (a) Lemma Phi: the integral equals the product formulas for Im x > 0 and Im x < 0; at real x it satisfies the two functional equations
     and the inversion relation.
 (b) Lemma S (strip):  Lam(w) = - sum_{j<K} B_2j(1/2) t^(2j-1) Li_(2-2j)(-e^w)/(2j)! + O(|t|^(2K-1)),  uniformly near real x;
     defect/|t|^(2K-1) printed for |t| halving.
 (c) Lemma F:  |Lam(w) + Li2(-e^w)/t| <= C e^(Re w/2)  for Re w <= 0, |Im w| <= pi;  sup of the left side times e^(-Re w/2) printed.
 (d) uniform bound on the contour (delta = 10 deg):  sup |H| / exp(psi(rho)/|t|),  H = C0 e^(t s^2/2)/(Phi_b(b s) 2 cosh(pi s));
     and the margin of the dual factors (Lemma D).
 (e) the kernel integral at real x: against the Borel sum of Phi, and ( P J - sum_{k<K} a_k t^k )/(a_K t^K).
Parts a, b: mpmath.  Parts c, d, e: numpy double precision (accuracy about 1e-12).     Usage: python3 z1_checks.py [a b c d e]"""
import sys
import numpy as np
from scipy.special import spence
from mpmath import mp, mpf, mpc, exp, log, pi, sqrt, polylog, bernpoly, factorial, nstr
from fad import Lam_mp, Phi_prod_mp, logPhi_np, Lam_np, I
DELTA = 10.0
Li2 = lambda z: spence(1 - z)

def part_a():
    mp.dps = 25
    print("(a) Faddeev integral against the product formulas: |integral/product - 1|")
    for th, T in ((-70, "0.5"), (-80, "0.3"), (-100, "0.3"), (-120, "0.6")):
        t = mpf(T)*exp(I*th*pi/180); x = I*t/(2*pi); out = []
        for w in (mpc("-0.7", "0.4"), mpc("-3", "-2.9"), mpc("0.3", "1.5")):
            out.append(nstr(abs(exp(Lam_mp(w, x))/Phi_prod_mp(w, x) - 1), 3))
        print(f"    theta = {th:5d}, |t| = {T}:  " + "   ".join(out))
    print("    double-precision version (strip integral + functional equation + inversion) against the products, |numpy/product - 1|:")
    for th, T in ((-80, "0.3"), (-100, "0.3"), (-70, "0.15")):
        t = mpf(T)*exp(I*th*pi/180); x = I*t/(2*pi); out = []
        for w in (mpc("-0.7", "0.4"), mpc("-5", "9"), mpc("-12", "-20"), mpc("2", "-1"), mpc("6", "-11")):
            out.append(f"{abs(np.exp(logPhi_np(complex(w), complex(x)))/complex(Phi_prod_mp(w, x)) - 1):.1e}")
        print(f"    theta = {th:5d}, |t| = {T}, w = -0.7+0.4i, -5+9i, -12-20i, 2-i, 6-11i:  " + "  ".join(out))
    print("    at real x (theta = -90): functional equations and inversion relation of the integral")
    for T in ("0.6", "0.2"):
        t = mpf(T)*exp(-I*pi/2); x = I*t/(2*pi); x = mpc(x.real, 0)
        for w in (mpc("-0.8", "0.0"), mpc("0.5", "0.0"), mpc("-2.5", "0.01")):   # w +- i pi must stay inside |Im w| < pi (1 + x)
            d1 = abs(Lam_mp(w - I*pi*x, x) - Lam_mp(w + I*pi*x, x) - log(1 + exp(w)))
            d2 = abs(Lam_mp(w - I*pi, x) - Lam_mp(w + I*pi, x) - log(1 + exp(w/x)))
            d3 = abs(Lam_mp(w, x) + Lam_mp(-w, x) - (-t/24 + pi**2/(6*t) + w*w/(2*t)))
            print(f"    |t| = {T}, w = {nstr(w, 3)}:  shift i b: {nstr(d1, 3)}   shift i/b: {nstr(d2, 3)}   inversion: {nstr(d3, 3)}")
        sys.stdout.flush()

def part_b():
    mp.dps = 30
    print("(b) Lemma S: |Lam(w) + sum_{j<K} B_2j(1/2) t^(2j-1) Li_(2-2j)(-e^w)/(2j)!| / |t|^(2K-1),  |t| = 0.4, 0.2, 0.1, 0.05")
    for th in (-90, -85, -100):
        for w in (mpc("-0.5", "1.2"), mpc("0.3", "-0.8"), mpc("-3", "2.5")):
            for K in (2, 3):
                out = []
                for T in ("0.4", "0.2", "0.1", "0.05"):
                    t = mpf(T)*exp(I*th*pi/180); x = I*t/(2*pi)
                    if th == -90: x = mpc(x.real, 0)
                    ser = -sum(bernpoly(2*j, mpf(1)/2)*t**(2*j - 1)*polylog(2 - 2*j, -exp(w))/factorial(2*j) for j in range(K))
                    out.append(nstr(abs(Lam_mp(w, x) - ser)/abs(t)**(2*K - 1), 3))
                print(f"    theta = {th:5d}, w = {nstr(w, 3):>14}, K = {K}:  " + "  ".join(out))
        sys.stdout.flush()

def part_c():
    print("(c) Lemma F: sup over Re w in [-12, 0], |Im w| <= pi of |Lam(w) + Li2(-e^w)/t| e^(-Re w/2)")
    res = np.linspace(-12, 0, 49); ims = np.linspace(-np.pi, np.pi, 17)
    for th in (-90, -75, -120, -60):
        row = []
        for T in (0.4, 0.2, 0.1, 0.05):
            t = T*np.exp(1j*np.radians(th)); x = 1j*t/(2*np.pi)
            if th == -90: x = complex(x.real, 0.0)
            best = 0; where = None
            for a in res:
                for b in ims:
                    w = complex(a, b)
                    v = abs(Lam_np(w, x) + Li2(-np.exp(w))/t)*np.exp(-a/2)
                    if v > best: best, where = v, w
            row.append(f"|t|={T}: {best:.4f} (at {where.real:.2f}{where.imag:+.2f}i)")
        print(f"    theta = {th:5d}:  " + "   ".join(row)); sys.stdout.flush()

def Gfun(w):
    if w.real <= 0: return -w*w/2 + Li2(-np.exp(w)) + np.pi**2/12
    return -w*w - Li2(-np.exp(-w)) - np.pi**2/12

def logH(r, t, x, case="psi1"):
    """log of the integrand of P J in the variable s = r e^(-i gamma); includes C0 and the extra factor for psi0."""
    beta = np.angle(np.sqrt(x)); ga = beta + np.radians(DELTA)
    s = r*np.exp(-1j*ga); w = 1j*t*s
    z = np.pi*s*(1 if s.real >= 0 else -1)
    l2cosh = z + np.log1p(np.exp(-2*z))
    out = np.pi**2/(12*t) + t/24 - 0.5*np.log(2) + t*s*s/2 - logPhi_np(w, x) - l2cosh
    return (out + w if case == "psi0" else out), w, ga

def part_d():
    print("(d) sup over the contour (delta = 10 deg, |w| <= 30, both halves) of |H| / exp(psi(rho)/|t|)")
    bs = np.radians(np.linspace(-30, 30, 1201)); d = np.radians(DELTA)
    marg = np.cos(bs + d) - np.abs(np.sin(bs - d)*np.sin(2*bs))
    print(f"    margin of the dual factors, min over |beta| <= 30 deg of cos(gamma) - |sin(alpha) sin(2 beta)|: {marg.min():.4f};"
          f"  min of pi/|sin(alpha)|: {(np.pi/np.abs(np.sin(bs - d))).min():.3f}")
    rhos = np.concatenate([np.linspace(0.02, 1, 50), np.linspace(1.1, 30, 290)])
    for th in (-90, -89, -91, -85, -95, -60, -120, -150, -30):
        row = []
        for T in (0.4, 0.2, 0.1):
            t = T*np.exp(1j*np.radians(th)); x = 1j*t/(2*np.pi)
            if th == -90: x = complex(x.real, 0.0)
            best = -1e9; where = None
            for sgn in (1, -1):
                for rho in rhos:
                    lh, w, ga = logH(sgn*rho/T, t, x)
                    psi = (np.exp(-1j*np.radians(th))*Gfun(w)).real - np.pi*rho*np.cos(ga)
                    v = lh.real - psi/T
                    if v > best: best, where = v, sgn*rho
            row.append(f"|t|={T}: {np.exp(best):.4f} (at {where:.2f})")
        print(f"    theta = {th:5d}:  " + "   ".join(row)); sys.stdout.flush()

def PJ_np(t, x, case="psi1", R=16.0, panels=64, order=12):
    xs, ws = np.polynomial.legendre.leggauss(order); tot = 0j
    edges = np.linspace(-R, R, panels + 1)
    for a, b in zip(edges[:-1], edges[1:]):
        for xi, wi in zip(xs, ws):
            r = (a + b)/2 + (b - a)/2*xi
            lh, w, ga = logH(r, t, x, case)
            tot += wi*(b - a)/2*np.exp(lh)*np.exp(-1j*ga)
    return tot*np.exp(t/8 - (t/2 if case == "psi0" else 0))

def part_e():
    from borel import PJ, borel_sum, coeff
    mp.dps = 20
    t = 0.3*np.exp(-1j*np.radians(80)); x = 1j*t/(2*np.pi)
    ref = PJ(mpc(x), delta=mpf(DELTA)*pi/180)
    print("(e) check of the double-precision integrator at theta = -80, |t| = 0.3 against the product-formula integrator of v13:",
          f"{abs(PJ_np(t, x) - complex(ref)):.2e}")
    for case in ("psi1", "psi0"):
        t = 0.3*np.exp(-1j*np.pi/2); x = complex((1j*t/(2*np.pi)).real, 0.0)
        v = PJ_np(t, x, case)
        bs = borel_sum(case, mpc(t), -pi/2, 60)
        print(f"    {case}, real x = {x.real:.6f} (|q| = 1, theta = -90, |t| = 0.3):  P J = {v:.12f}   Borel sum = {complex(bs):.12f}   difference = {abs(v - complex(bs)):.2e}")
        Ks = (1, 3, 6) if case == "psi1" else (2, 3, 6)
        vals = {}
        for T in (0.2, 0.1, 0.05):
            t = T*np.exp(-1j*np.pi/2); x = complex((1j*t/(2*np.pi)).real, 0.0); vals[T] = (t, PJ_np(t, x, case))
        for K in Ks:
            out = []
            for T in (0.2, 0.1, 0.05):
                t, v = vals[T]
                part = sum(complex(coeff(case, k))*t**k for k in range(K))
                out.append(f"{abs((v - part)/(complex(coeff(case, K))*t**K) - 1):.4f}")
            print(f"    {case}, theta = -90, K = {K}: |ratio - 1| at |t| = 0.2, 0.1, 0.05:  " + "   ".join(out))
        sys.stdout.flush()

if __name__ == "__main__":
    for pt in (sys.argv[1:] or ["a", "b", "c", "d", "e"]): {"a": part_a, "b": part_b, "c": part_c, "d": part_d, "e": part_e}[pt]()
