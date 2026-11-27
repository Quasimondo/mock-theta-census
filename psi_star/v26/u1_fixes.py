"""Addendum v26: checks for the corrections after the stress test.
 a  Theorem 7 (F4): on the horizontal line through the lower saddle, Re V has its maximum at the saddle for every 0 < a < 1.
    g(lam) = lam^2 |1 + lam e^(i eta1)|^2 crosses 1 exactly once; sup over the line of Re V - Re V(w1); the value of a where concavity ends.
 b  Theorem 7 at arg x = 30, 45, 60 degrees: (I(a/b)/asymptotic formula - 1)/|t|, I(a/b) from the series of Theorem C (mpmath).
 c  the e^w-version on the horizontal line at real x (F2): slope of log|integrand| on the right against 1 - (pi/3)/|t|;
    I(b + 1/(2b)) over the line against the value over L where the line is admissible.
 d  analyticity of P J across the real x-axis (F3): P J over the fixed line L = e^(-10 deg i) R in the variable y, on a circle in the t-plane
    around a point with real x; negative Fourier modes and the mean value.
Usage: python3 u1_fixes.py [a b c d]"""
import sys
import numpy as np
from mpmath import mp, mpf, mpc, exp, pi, sqrt, log, polylog, nstr
from fad import logPhi_np

def lower_root(a):
    e = np.exp(-2j*np.pi*a); d = np.sqrt(1 + 4*e + 0j)
    r = [(-1 + d)/2, (-1 - d)/2]
    return [W for W in r if W.imag < 0][0]

def part_a():
    print("(a) Re V on the horizontal line through the lower saddle (V = -w^2/2 - 2 pi i a w + Li2(-e^w))")
    lam = np.exp(np.linspace(-30, 30, 600001)); worst = -1e9; multi = 0; n = 0
    for a in list(np.linspace(0.001, 0.999, 999)):
        W = lower_root(a); eta = np.angle(W); c = np.cos(eta)
        g = lam**2*(1 + 2*c*lam + lam**2)
        cross = np.count_nonzero(np.diff(np.sign(g - 1)) != 0); multi += (cross != 1); n += 1
        # Re V' = -(1/2) log g, so Re V(r) - Re V(r1) = -(1/2) int_{r1}^{r} log g
        r = np.log(lam); f = -0.5*np.log(g); F = np.concatenate([[0], np.cumsum((f[1:] + f[:-1])/2*np.diff(r))])
        i1 = np.argmin(np.abs(lam - abs(W))); worst = max(worst, np.max(F - F[i1]))
    print(f"    {n} values of a in [0.001, 0.999]: g = 1 has one solution in all but {multi} cases;  sup (Re V - Re V(w1)) = {worst:.2e}")
    # concavity ends where cos(eta1) = -2 sqrt(2)/3
    lo, hi = 0.5, 0.999
    for _ in range(60):
        mid = (lo + hi)/2
        if np.cos(np.angle(lower_root(mid))) > -2*np.sqrt(2)/3: lo = mid
        else: hi = mid
    print(f"    Re V is concave on the whole line for a <= {lo:.6f} (eta1 = -{np.degrees(-np.angle(lower_root(lo))):.3f} deg); beyond that it is not concave but still maximal at w1")
    for a in (0.9, 0.99):
        W = lower_root(a); c = np.cos(np.angle(W)); lm = (3*abs(c) - np.sqrt(9*c*c - 8))/4
        print(f"    a = {a}: local maximum of g at lam = {lm:.4f}: g = {lm**2*(1 + 2*c*lm + lm**2):.4f} < 1")

def part_b():
    from common import Setup, I
    from a1_formula import Iser
    print("(b) Theorem 7 off the real axis: (I(a/b)/asymptotic - 1)/|t|   (series of Theorem C, 90 digits)")
    mp.dps = 90
    for a in ("0.1", "0.3", "0.5", "0.8", "0.95"):
        a = mpf(a); e = exp(-2j*pi*a); d = sqrt(1 + 4*e); W = [v for v in ((-1 + d)/2, (-1 - d)/2) if v.imag < 0][0]
        w1 = log(W); V = -w1*w1/2 - 2j*pi*a*w1 + polylog(2, -W)
        row = []
        for deg in (30, 45, 60):
            for T in ("0.4", "0.2", "0.1"):
                t = mpf(T)*exp(1j*pi*(deg - 90)/180); x = 1j*t/(2*pi)
                S = Setup(x); b = S.b
                ex = Iser(S, a/b); asy = sqrt(2*pi*t*(1 + W)/(1 + 2*W))*exp(V/t)/(2*pi*b)
                row.append(nstr(abs(ex/asy - 1)/mpf(T), 3))
        print(f"    a = {nstr(a, 3):5s} arg x = 30: {', '.join(row[0:3])};  45: {', '.join(row[3:6])};  60: {', '.join(row[6:9])}   (|t| = 0.4, 0.2, 0.1)")

def part_c():
    from h1_roots import PI_np, PI_w
    print("(c) the e^w-version of the decaying integral on the line Im w = -2 pi/3 at x = 1/X")
    for X in (4, 6, 7, 8, 10, 12):
        x = 1.0/X; t = -2j*np.pi*x
        def lg(r):
            w = r - 2j*np.pi/3; s = w/(1j*t)
            return (t*s*s/2 - logPhi_np(w, x) + np.pi*s + w).real
        slope = (lg(14.0) - lg(10.0))/4.0
        line = f"    X = {X:2d}: slope of log|integrand| on the right = {slope:+.4f}, predicted 1 - (pi/3)/|t| = {1 - X/6:+.4f}"
        if X > 6:
            T = abs(t); Rl = 14*T + 2; Rr = max(Rl, 40/(np.pi/(3*T) - 1))
            f = lambda w: np.exp(np.pi**2/(12*t) + t/24 - 0.5*np.log(2) + t*(w/(1j*t))**2/2 - logPhi_np(w, x) + np.pi*w/(1j*t) + w - t/2)
            xs, ws = np.polynomial.legendre.leggauss(12); a = 0j; edges = np.linspace(-Rl, Rr, 4001)
            for u, v in zip(edges[:-1], edges[1:]):
                for xi, wi in zip(xs, ws):
                    a += wi*(v - u)/2*f((u + v)/2 + (v - u)/2*xi - 2j*np.pi/3)
            a *= np.exp(t/8)/(1j*t); b = PI_np(X, "+", True, R=9*np.sqrt(X), panels=1600)
            line += f";  line against L: relative difference {abs(a - b)/abs(b):.1e}"
        else:
            line += ";  the line is not admissible"
        print(line)

def PJ(x, R=30.0, panels=500, order=12):
    """P J(1/(2b)) over the fixed line y = r e^(-10 deg i); complex x with Re x > 0 (Faddeev's integral)."""
    x = complex(x); t = -2j*np.pi*x; b = np.sqrt(x); e = np.exp(-1j*np.radians(10.0))
    xs, ws = np.polynomial.legendre.leggauss(order); tot = 0j
    edges = np.linspace(-R, R, panels + 1)*abs(b)
    for u, v in zip(edges[:-1], edges[1:]):
        for xi, wi in zip(xs, ws):
            y = ((u + v)/2 + (v - u)/2*xi)*e; s = y/b; w = 1j*t*s
            z = np.pi*s*(1 if s.real >= 0 else -1)
            lgv = np.pi**2/(12*t) + t/24 - 0.5*np.log(2) + t*s*s/2 - logPhi_np(w, x) - z - np.log1p(np.exp(-2*z))
            tot += wi*(v - u)/2*np.exp(lgv)*e/b
    return tot*np.exp(t/8)

def part_d():
    print("(d) P J over the fixed line L on circles |t - t0| = 0.3 |t0| around real x0 (24 points): Fourier modes c_m = mean of P J e^(-i m phi)")
    for x0 in (0.1, 0.04):
        t0 = -2j*np.pi*x0; N = 24; rad = 0.3*abs(t0)
        ph = 2*np.pi*np.arange(N)/N; vals = np.array([PJ(1j*(t0 + rad*np.exp(1j*p))/(2*np.pi)) for p in ph])
        c = lambda m: np.mean(vals*np.exp(-1j*m*ph))
        print(f"    x0 = {x0}: |c_0 - P J(t0)| = {abs(c(0) - PJ(x0)):.1e};  |c_-1|, |c_-2|, |c_-3| = {abs(c(-1)):.1e}, {abs(c(-2)):.1e}, {abs(c(-3)):.1e};"
              f"  |c_1|, |c_2|, |c_3| = {abs(c(1)):.1e}, {abs(c(2)):.1e}, {abs(c(3)):.1e}   (arg x on the circle from {np.degrees(np.arcsin(-0.3)):.1f} to {np.degrees(np.arcsin(0.3)):.1f} deg)")

if __name__ == "__main__":
    for p in (sys.argv[1:] or list("abcd")):
        {"a": part_a, "b": part_b, "c": part_c, "d": part_d}[p]()
