"""Shared code for the second round: the z-family saddles, V with the continuation of v22 §1,
the five-piece contour C of v22 §2 (and the tilted C' of §3), the state integrals over C and over
tilted lines, the finite sums F(z; q0).  Phi_b from the round-1 library psi_star/stress/core.py
(products and Faddeev's integral, own code)."""
import os
import sys
import numpy as np
import mpmath as mp
from math import pi, log, sqrt, sin, cos, tan

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "stress"))
import core  # noqa: E402

DEG = pi / 180
gl_x, gl_w = np.polynomial.legendre.leggauss(40)


# ------------------------------------------------------------------ saddles and volume
def W_grow(a):
    """root of W^2 + W = e^{-2 pi i a} with Im W > 0 (v22)."""
    r = np.roots([1, 1, -np.exp(-2j * pi * a)])
    return complex(r[np.argmax(r.imag)])


def W_dec(a):
    """root with Im W < 0 (v19 Theorem 7)."""
    r = np.roots([1, 1, -np.exp(-2j * pi * a)])
    return complex(r[np.argmin(r.imag)])


def bloch_wigner(z):
    z = mp.mpc(z)
    return float(mp.im(mp.polylog(2, z)) + mp.arg(1 - z) * mp.log(abs(z)))


def vol(a):
    return abs(bloch_wigner(-W_grow(a)))


# ------------------------------------------------------------------ dilogarithm on C (v22 §1)
def li2_real_gt1_minus(xx):
    """Li2(x - i0) for real x > 1: Re Li2(x) - i pi ln x."""
    xx = np.asarray(xx, dtype=float)
    re = pi ** 2 / 3 - 0.5 * np.log(xx) ** 2 - core.li2(1 / xx).real
    return re - 1j * pi * np.log(xx)


def li2_cont(w):
    """Li2(-e^w) as v22 §1 continues it: principal inside the closed strip |Im w| <= pi with boundary
    values from inside, and the 2 pi i-periodic power series for Re w <= 0 (any Im w)."""
    w = np.atleast_1d(np.asarray(w, dtype=complex))
    out = np.empty_like(w)
    left = w.real <= 0
    # Re w <= 0: power series in -e^w, |e^w| <= 1  -> principal Li2 of -e^w (no cut met)
    out[left] = core.li2(-np.exp(w[left]))
    R = ~left
    wr = w[R]
    if np.any(np.abs(wr.imag) > pi + 1e-12):
        raise ValueError("Re w > 0 outside the strip: not defined by v22 §1")
    o = core.li2(-np.exp(wr))
    top = np.abs(np.abs(wr.imag) - pi) < 1e-12
    if np.any(top):
        # -e^w = e^{Re w} > 1 on Im w = +-pi; from inside: Im w -> pi^- gives e^{r} e^{-i0}, -pi^+ gives +i0
        xx = np.exp(wr[top].real)
        v = li2_real_gt1_minus(xx)
        sgn = np.sign(wr[top].imag)
        o[top] = np.where(sgn > 0, v, np.conj(v))
    out[R] = o
    return out


def V(w, ap):
    """V(w) = -w^2/2 - 2 pi i a' w + Li2(-e^w), with a' the shift (a - 1 for the growing integral)."""
    w = np.atleast_1d(np.asarray(w, dtype=complex))
    return -w ** 2 / 2 - 2j * pi * ap * w + li2_cont(w)


# ------------------------------------------------------------------ contour C of v22 §2
def contour_C(a, n=4000, Rl=None, Rr=None, phi1=None):
    """Points of C1..C5 (or C1', C2' tilted by phi1 degrees), returned as list of (name, w, param)
    with clustering near the corners and the saddle."""
    W = W_grow(a)
    w0 = np.log(W)
    r0, e0 = w0.real, w0.imag
    if Rl is None:
        Rl = 60.0
    if Rr is None:
        Rr = 60.0
    g = lambda L, k: np.concatenate([[0.0], np.geomspace(1e-7, L, k)])
    pieces = []
    if phi1 is None:
        pieces.append(("C1", -g(Rl, n)[::-1] + 2j * pi))
        pieces.append(("C2", 1j * (2 * pi - np.concatenate([g(pi / 2, n // 2), pi - g(pi / 2, n // 2)[::-1]]))))
    else:
        P = -pi * tan(phi1 * DEG) + 2j * pi
        pieces.append(("C1'", (-g(Rl, n)[::-1] + P.real) + 2j * pi))
        s = np.concatenate([g(0.5, n // 2), 1 - g(0.5, n // 2)[::-1]])
        pieces.append(("C2'", P + s * (1j * pi - P)))
    s = np.concatenate([g(0.5, n // 2), 1 - g(0.5, n // 2)[::-1]])
    pieces.append(("C3", r0 * s + 1j * pi))
    # C4: eta from pi to 0, clustered at pi, at 0 and at eta0
    e = np.unique(np.concatenate([pi - g((pi - e0) / 2, n // 2), e0 + g((pi - e0) / 2, n // 2),
                                  e0 - g(e0 / 2, n // 2), g(e0 / 2, n // 2)]))[::-1]
    pieces.append(("C4", r0 + 1j * e))
    pieces.append(("C5", r0 + g(Rr, n)))
    return pieces, w0


def dV_near(w, w0):
    """V(w) - V(w0) by the Taylor series at the saddle (V'(w0) = 0), for |w - w0| < 1e-3, where the
    direct difference is dominated by rounding."""
    e = np.exp(w0)
    s = e / (1 + e)
    d = w - w0
    V2 = -1 - s
    V3 = -s * (1 - s)
    V4 = -s * (1 - s) * (1 - 2 * s)
    V5 = -s * (1 - s) * (1 - 6 * s + 6 * s * s)
    return V2 * d ** 2 / 2 + V3 * d ** 3 / 6 + V4 * d ** 4 / 24 + V5 * d ** 5 / 120


def dV(w, w0, ap, v0):
    w = np.atleast_1d(w)
    out = c_V = V(w, ap) - v0
    near = np.abs(w - w0) < 1e-3
    if np.any(near):
        out = np.where(near, dV_near(w, w0), c_V)
    return out


def heights_L1(a, tau=0.0, n=4000, phi1=None):
    """sup over each piece of (H_tau(w) - H_tau(w0)) / min(d^2, d), H_tau = -Im V + tau Re V."""
    pcs, w0 = contour_C(a, n=n, phi1=phi1)
    ap = a - 1
    v0 = complex(V(np.array([w0]), ap)[0])
    H0 = -v0.imag + tau * v0.real
    out = {}
    for name, w in pcs:
        d = np.abs(w - w0)
        m = d > 1e-9
        dv = dV(w[m], w0, ap, v0)
        q = (-dv.imag + tau * dv.real) / np.minimum(d[m] ** 2, d[m])
        i = int(np.argmax(q))
        out[name] = (q[i], w[m][i])
    return out, w0, v0


# ------------------------------------------------------------------ integrals at real x
def panels_osc(a_, b_, c, tabs, base=1.5, cap=60000):
    edges = [a_]
    while edges[-1] < b_:
        r = edges[-1]
        h = base * tabs / (1 + 0.5 * abs(r - c))
        edges.append(min(b_, r + max(h, 1e-7)))
        if len(edges) > cap:
            raise RuntimeError("too many panels")
    e = np.array(edges)
    m = (e[1:] + e[:-1]) / 2
    h = (e[1:] - e[:-1]) / 2
    return (m[:, None] + h[:, None] * gl_x[None, :]).ravel(), (h[:, None] * gl_w[None, :]).ravel()


def lam(w, t, fad):
    return core.lam_any(w, t, fad=fad, cap=40000)


def I_grow_C(a, X=None, x=None, fad=None, extra=0.0, phi1=None):
    """I((a-1)/b) = (1/(2 pi b)) int_C e^{Q/t}/Phi_b dw over C (or C' with tilt phi1), at x (real or
    complex), returned as I * e^{-V(w0)/t} together with t, b, v0.  extra = c adds e^{c w}."""
    if x is None:
        x = 1.0 / X
    x = complex(x)
    t = -2j * pi * x
    tabs = abs(t)
    b = np.sqrt(x)
    if fad is None:
        fad = core.Faddeev(t)
    ap = a - 1
    W = W_grow(a)
    w0 = np.log(W)
    r0, e0 = w0.real, w0.imag
    v0 = complex(V(np.array([w0]), ap)[0])

    def f(w):
        return np.exp((-w ** 2 / 2 - 2j * pi * ap * w) / t - lam(w, t, fad) - v0 / t + extra * w)

    tot = 0j
    # C1: r from -inf to 0 at Im 2 pi (or to Re P); decay e^{-2 pi a |r| / |t|}
    Rl = 70 * tabs / (2 * pi * a) + 2
    if phi1 is None:
        rr, wt = panels_osc(-Rl, 0.0, 0.0, tabs)
        tot += (f(rr + 2j * pi) * wt).sum()
        # C2: w = i eta, eta 2 pi -> pi: dw = i d eta (reversed)
        ee, wt = panels_osc(pi, 2 * pi, 1.5 * pi, tabs)
        tot += -1j * (f(1j * ee) * wt).sum()
    else:
        P = -pi * tan(phi1 * DEG) + 2j * pi
        rr, wt = panels_osc(-Rl + P.real, P.real, P.real, tabs)
        tot += (f(rr + 2j * pi) * wt).sum()
        ss, wt = panels_osc(0.0, 1.0, 0.5, tabs / abs(1j * pi - P))
        tot += (f(P + ss * (1j * pi - P)) * wt).sum() * (1j * pi - P)
    # C3: r from 0 to r0 at Im pi
    if abs(r0) > 0:
        ss, wt = panels_osc(0.0, 1.0, 0.5, tabs / abs(r0))
        tot += (f(ss * r0 + 1j * pi) * wt).sum() * r0
    # C4: eta pi -> 0 at Re r0
    ee, wt = panels_osc(0.0, pi, e0, tabs)
    tot += -1j * (f(r0 + 1j * ee) * wt).sum()
    # C5: r0 -> inf; decay e^{-2 pi (1-a) r / |t|}
    Rr = 70 * tabs / (2 * pi * (1 - a)) + 2 + abs(r0)
    rr, wt = panels_osc(r0, r0 + Rr, r0, tabs)
    tot += (f(rr + 0j) * wt).sum()
    return tot / (2 * pi * b), t, b, v0


def asym13(a, t, b):
    """Theorem 13 (a) leading term divided by e^{V(w0)/t}: (1/(2 pi b)) (2 pi t / (-V''))^{1/2},
    root whose quotient by -i has positive real part."""
    W = W_grow(a)
    Vpp = -(1 + 2 * W) / (1 + W)
    r = np.sqrt(2 * pi * t / (-Vpp))
    if (r / (-1j)).real < 0:
        r = -r
    return r / (2 * pi * b)


def line_integral(x, ap, kernel=False, delta=45.0, fad=None, n=None, extra=0.0, Rmax=None):
    """(1/?) int over the line e^{-i delta} R of e^{-pi i y^2 + 2 pi u y}/Phi_b(y) [/(1 + e^{2 pi y/b})]
    with u = ap/b (+ extra*b), in the w variable: I(u) = (1/(2 pi b)) int e^{Q/t + extra w}/Phi_b dw.
    Returns (value, conditioning = max|integrand| * length / |value|)."""
    x = complex(x)
    t = -2j * pi * x
    b = np.sqrt(x)
    if fad is None and x.real > 0 and (1 / x).real > 0.85:
        fad = core.Faddeev(t)
    e = np.exp(-1j * (delta * DEG - np.angle(b)))     # direction of w = 2 pi b y
    # size of the Gaussian along the line: |e^{-w^2/(2t)}|; choose R where it reaches e^{-60} (left) and
    # use the right-half doubling of the Gaussian (1/Phi_b ~ e^{-w^2/(2t)} there)
    g = (-(e ** 2) / (2 * t)).real
    if g >= 0:
        raise ValueError("line not admissible")
    lin = abs(2 * pi * ap) / abs(t) + abs(extra) + 1
    if kernel:
        lin += 1 / abs(x)          # 1/(1 + e^{w/x}) can grow like e^{|w|/|x|} on one side for complex x
    R = (lin + np.sqrt(lin ** 2 + 4 * (-g) * 70)) / (2 * (-g)) + 2
    if Rmax is not None:
        R = min(R, Rmax)
    rr, wt = panels_osc(-R, R, 0.0, abs(t), base=1.0)
    w = rr * e
    lf = (-w ** 2 / 2 - 2j * pi * ap * w) / t - lam(w, t, fad) + extra * w
    if kernel:
        lf = lf - core.log1pexp(w / x)
    fv = np.exp(lf)
    val = (fv * wt).sum() * e / (2 * pi * b)
    cond = np.abs(fv).max() * 2 * R / (2 * pi * abs(b)) / max(abs(val), 1e-300)
    return val, cond


def Pfac(t, b):
    return 2 ** -0.5 * np.exp(pi ** 2 / (12 * t) + t / 6) / b


def F_finite(a, X, dps=None):
    """F(z; q0) = sum_{n < X/2} z^n q0^{n(n+1)/2} (-q0;q0)_n, z = e^{2 pi i a}, q0 = e^{2 pi i/X}."""
    if dps is None:
        dps = 30 + int(0.08 * X)
    with mp.workdps(dps):
        z = mp.expjpi(2 * mp.mpf(a))
        q0 = mp.expjpi(mp.mpf(2) / X)
        s = mp.mpc(0)
        P = mp.mpc(1)
        for n in range((X + 1) // 2):
            if n > 0:
                P *= 1 + q0 ** n
            s += z ** n * q0 ** (n * (n + 1) // 2) * P
        return s


def F_prediction(a, X):
    """Theorem 15 main term (1 + W) W^{1/2} (pi/(t(1+2W)))^{1/2} e^{(V(w0) + pi^2/12)/t}, as mp."""
    with mp.workdps(40):
        a = mp.mpf(a)
        z = mp.expjpi(2 * a)
        r = mp.polyroots([1, 1, -1 / z])
        W = r[0] if mp.im(r[0]) > 0 else r[1]
        w0 = mp.log(W)
        V0 = -w0 ** 2 / 2 - 2j * mp.pi * (a - 1) * w0 + mp.polylog(2, -W)
        t = -2j * mp.pi / X
        return (1 + W) * mp.exp(w0 / 2) * mp.sqrt(mp.pi / (t * (1 + 2 * W))) * mp.exp((V0 + mp.pi ** 2 / 12) / t)
