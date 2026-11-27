"""S6: Part 6, the state integrals [v19].

 (a) closed formulas for Re V, Im V on Gamma_1, Gamma_2, Gamma_3 and the bounds derived
     from them; (L1) on Gamma for 0 <= tau <= 0.2, then beyond 0.2.
 (b) monotonicity of Im V along the horizontal line through w_1 for random a in (0,1).
 (c) Theorem 6: (exact / asymptotic - 1) / |t| for the growing pair (0 <= arg x <= 11 deg)
     and the decaying pair (0 <= arg x <= 60 deg), |t| from 0.5 down to 0.01; the exact
     integrals by independent quadrature on the descent contours, cross-checked on L.
 (d) Theorem 7: epsilon(a) from (L1) on the line, and (exact/asymptotic - 1)/|t| for random a.
Run:  python3 s6_state.py     (about 8 min)
"""
import time
import numpy as np
import mpmath as mp
from math import pi, log, cos, sin, atan
import core

DEG = pi / 180
rng = np.random.default_rng(6)
CL2 = core.CL2
gl_x, gl_w = np.polynomial.legendre.leggauss(40)


def V(w, ap):
    w = np.asarray(w, dtype=complex)
    return -w ** 2 / 2 - 2j * pi * ap * w + core.li2_neg_exp(w)


def clausen2(th):
    return float(mp.clsin(2, th))


def panels_osc(a, b, c, tabs, base=1.5):
    """Gauss-Legendre panels on [a, b] whose width shrinks like 1/(1 + |r - c|): the phase of
    e^{V/t} turns at a rate ~ (|r| + const)/|t| along the contours used here."""
    edges = [a]
    while edges[-1] < b:
        r = edges[-1]
        h = base * tabs / (1 + 0.5 * abs(r - c))
        edges.append(min(b, r + max(h, 1e-6)))
        if len(edges) > 30000:
            raise RuntimeError("too many panels")
    e = np.array(edges)
    X, Wt = [], []
    for i in range(len(e) - 1):
        m, h = (e[i] + e[i + 1]) / 2, (e[i + 1] - e[i]) / 2
        X.append(m + h * gl_x)
        Wt.append(h * gl_w)
    return np.concatenate(X), np.concatenate(Wt)


def panels(a, b, n):
    """composite Gauss-Legendre nodes/weights on [a, b] with n panels."""
    e = np.linspace(a, b, n + 1)
    X, Wt = [], []
    for i in range(n):
        m, h = (e[i] + e[i + 1]) / 2, (e[i + 1] - e[i]) / 2
        X.append(m + h * gl_x)
        Wt.append(h * gl_w)
    return np.concatenate(X), np.concatenate(Wt)


# ------------------------------------------------------------------ (a)
def part_a():
    print("(a) Gamma_1, Gamma_2, Gamma_3: closed formulas and bounds")
    eta = rng.uniform(0, pi, 2000)
    v2 = V(1j * eta, -0.5)
    e_im = max(abs(v2.imag[i] + clausen2(pi - eta[i])) for i in range(200))
    e_re = np.abs(v2.real - (0.75 * eta ** 2 - pi * eta - pi ** 2 / 12)).max()
    print("   Gamma_2: max |Im V + Cl2(pi - eta)| = %.1e, max |Re V - (3/4 eta^2 - pi eta - pi^2/12)| = %.1e"
          % (e_im, e_re))
    r = np.concatenate([rng.uniform(0, 30, 2000), [0.0]])
    v3 = V(r, -0.5)
    print("   Gamma_3: max |Im V - pi r| = %.1e;  max of Re V + pi^2/12 + r^2/2 = %.2e (claim <= 0)"
          % (np.abs(v3.imag - pi * r).max(), (v3.real + pi ** 2 / 12 + r ** 2 / 2).max()))
    rho = np.concatenate([rng.uniform(0, 40, 4000), np.linspace(0, 2 * pi, 500)])
    w1 = 1j * pi - rho * np.exp(-1j * pi / 6)
    v1 = V(w1, -0.5)
    Q = -w1 ** 2 / 2 + 1j * pi * w1
    print("   Gamma_1: max |Re Q + pi^2/2 + rho^2/4| = %.1e, max |Im Q - rho^2 sin(pi/6)cos(pi/6)| = %.1e"
          % (np.abs(Q.real + pi ** 2 / 2 + rho ** 2 / 4).max(), np.abs(Q.imag - rho ** 2 * sin(pi / 6) * cos(pi / 6)).max()))
    print("            min (Im V - 0.433 rho^2 + 0.01) = %.3e (claim >= 0);  max (Re V + pi^2/3) = %.3e (claim <= 0)"
          % ((v1.imag - 0.433 * rho ** 2 + 0.01).min(), (v1.real + pi ** 2 / 3).max()))
    m = rho <= 2 * pi
    print("            min Im Li2(-e^w) for rho <= 2pi: %.3e (claim >= 0); max |Li2| for rho >= 2pi: %.4f (claim < 0.01)"
          % (core.li2_neg_exp(w1[m]).imag.min(), np.abs(core.li2_neg_exp(w1[~m])).max()))
    # g(eta) on Gamma_2 and the quadratic bound
    for tau in [0.0, 0.1, 0.2]:
        e = np.linspace(0, pi, 20001)
        g = np.array([clausen2(pi - x) for x in e[::20]]) + tau * (0.75 * e[::20] ** 2 - pi * e[::20] - pi ** 2 / 12)
        e2 = e[::20]
        g0 = CL2 + tau * (0.75 * (2 * pi / 3) ** 2 - pi * 2 * pi / 3 - pi ** 2 / 12)
        d = g - g0 + 0.015 * (e2 - 2 * pi / 3) ** 2
        print("   Gamma_2, tau=%.1f: max of g - g(2pi/3) + 0.015 (eta - 2pi/3)^2 = %.3e (claim <= 0)" % (tau, d.max()))
    # (L1) on the whole contour
    print("   (L1) on Gamma: c_eff(tau) = inf over Gamma of (H(w0) - H(w)) / min(|w-w0|^2, |w-w0|)")
    w0 = 2j * pi / 3
    pts = np.concatenate([1j * pi - np.geomspace(1e-4, 40, 800) * np.exp(-1j * pi / 6),
                          1j * np.linspace(0, pi, 4001), np.geomspace(1e-4, 40, 800)])
    pts = pts[np.abs(pts - w0) > 1e-6]
    vv = V(pts, -0.5)
    vw0 = complex(V(np.array([w0]), -0.5)[0])
    print("   H_tau(w0) = Cl2 - tau 5pi^2/12; V(w0) = %.12f%+.12fi  (claim -5pi^2/12 = %.12f, -Cl2 = %.12f)"
          % (vw0.real, vw0.imag, -5 * pi ** 2 / 12, -CL2))
    line = "  "
    tau_fail = None
    for tau in [0.0, 0.05, 0.1, 0.15, 0.2, 0.22, 0.2206, 0.24, 0.2468, 0.26, 0.3, 0.4, 0.5, 0.577, 0.6]:
        H = -vv.imag + tau * vv.real
        H0 = -vw0.imag + tau * vw0.real
        d = np.abs(pts - w0)
        c = ((H0 - H) / np.minimum(d ** 2, d)).min()
        line += " tau=%.4f:%.4f" % (tau, c)
        if c <= 0 and tau_fail is None:
            tau_fail = tau
    print(line)
    # bisection for the first tau with c_eff <= 0
    lo, hi = 0.2, 0.7
    for _ in range(40):
        mid = (lo + hi) / 2
        H = -vv.imag + mid * vv.real
        H0 = -vw0.imag + mid * vw0.real
        d = np.abs(pts - w0)
        if ((H0 - H) / np.minimum(d ** 2, d)).min() > 0:
            lo = mid
        else:
            hi = mid
    print("   (L1) on this Gamma holds for tau < %.4f, i.e. arg x < %.2f deg (claimed: tau <= 0.2, 11.3 deg)"
          % (lo, atan(lo) / DEG))


# ------------------------------------------------------------------ helpers for (b)-(d)
def W_of(a):
    r = np.roots([1, 1, -np.exp(-2j * pi * a)])
    return r[np.argmin(r.imag)]


def line_heights(a, tau, r):
    W = W_of(a)
    w1 = np.log(W)
    w = r + 1j * w1.imag
    v = V(w, a)
    H = -v.imag + tau * v.real
    v1 = complex(V(np.array([w1]), a)[0])
    return H, -v1.imag + tau * v1.real, w1


def part_b():
    print("\n(b) horizontal line through w1: d/dr Im V < 0 for r < r1, > 0 for r > r1")
    bad = 0
    tot = 0
    worst_valley = 1e9
    for a in np.concatenate([rng.uniform(0, 1, 300), [1e-3, 1e-2, 0.99, 0.999, 0.5]]):
        W = W_of(a)
        w1 = np.log(W)
        r1, e1 = w1.real, w1.imag
        r = r1 + np.concatenate([-np.geomspace(1e-6, 60, 400)[::-1], np.geomspace(1e-6, 60, 400)])
        w = r + 1j * e1
        dIm = (-w - 2j * pi * a - np.log1p(np.exp(w))).imag     # Im V'
        sgn_ok = np.all(dIm[r < r1] < 0) and np.all(dIm[r > r1] > 0)
        bad += (not sgn_ok)
        tot += 1
        worst_valley = min(worst_valley, e1 + 2 * pi * a, -pi * a - e1)
    print("   %d values of a (incl. 1e-3, 0.999): sign pattern violated for %d;  min of the valley margins"
          " (eta1 + 2 pi a, -pi a - eta1) = %.3e (claim > 0)" % (tot, bad, worst_valley))


def eps_of_a(a):
    """(L1) on the horizontal line: c0 = inf (ImV - ImV1)/m(r), M = sup (ReV - ReV1)/m(r);
    sufficient tau < c0/M; and the actual largest tau with c_eff > 0."""
    W = W_of(a)
    w1 = np.log(W)
    r = w1.real + np.concatenate([-np.geomspace(1e-5, 80, 1500)[::-1], np.geomspace(1e-5, 80, 1500)])
    w = r + 1j * w1.imag
    v = V(w, a)
    v1 = complex(V(np.array([w1]), a)[0])
    d = np.abs(r - w1.real)
    m = np.minimum(d ** 2, d)
    c0 = ((v.imag - v1.imag) / m).min()
    M = max(((v.real - v1.real) / m).max(), 0)
    lo, hi = 0.0, 50.0
    for _ in range(60):
        mid = (lo + hi) / 2
        c = (((v.imag - v1.imag) - mid * (v.real - v1.real)) / m).min()
        if c > 0:
            lo = mid
        else:
            hi = mid
    return c0, M, (c0 / M if M > 0 else np.inf), lo


# ------------------------------------------------------------------ integrals
def I_decay(tabs, argx, a, extra_ew=False, fad=None):
    """I(a/b) [or I(b + a/b)] at x = |x| e^{i argx}, on the horizontal line through w1.
    Returns I * e^{-V(w1)/t}."""
    x = tabs / (2 * pi) * np.exp(1j * argx * DEG)
    t = -2j * pi * x
    b = np.sqrt(x)
    W = W_of(a)
    w1 = np.log(W)
    e1 = w1.imag
    rate = min(e1 + 2 * pi * a, -(2 * e1 + 2 * pi * a)) * cos(argx * DEG)
    if extra_ew:
        # the factor e^{w} costs e^{r} on the right: decay needs -(2 eta1 + 2 pi a) > |t|
        rate = min(e1 + 2 * pi * a, -(2 * e1 + 2 * pi * a) * cos(argx * DEG) - tabs)
        if rate <= 0:
            return np.nan, t, b, complex(V(np.array([w1]), a)[0])
    R = min(60 * tabs / max(rate, 1e-3), 400.0) + 12 * np.sqrt(tabs)
    # refine R from the actual size of the integrand (the rate above is pessimistic for complex x)
    rc = w1.real + np.linspace(-R, R, 801)
    wc = rc + 1j * e1
    v1c = complex(V(np.array([w1]), a)[0])
    lm = ((-wc ** 2 / 2 - 2j * pi * a * wc) / t - core.lam_any(wc, t, fad=fad, cap=40000) - v1c / t
          + (wc if extra_ew else 0)).real
    big = np.nonzero(lm > lm.max() - 50)[0]
    if big[0] == 0 or big[-1] == len(rc) - 1:
        if R >= 400:
            return np.nan, t, b, v1c
    else:
        R = max(abs(rc[big[0] - 1] - w1.real), abs(rc[big[-1] + 1] - w1.real)) + 1
    r, wt = panels_osc(w1.real - R, w1.real + R, w1.real, tabs)
    w = r + 1j * e1
    v1 = complex(V(np.array([w1]), a)[0])
    lam = core.lam_any(w, t, fad=fad, cap=40000)
    Qw = -w ** 2 / 2 - 2j * pi * a * w
    f = np.exp(Qw / t - lam - v1 / t + (w if extra_ew else 0))
    return (f * wt).sum() / (2 * pi * b), t, b, v1


def I_grow(tabs, argx, extra_ew=False, fad=None):
    """I(-1/(2b)) [or I(b - 1/(2b))] on Gamma_1 + Gamma_2 + Gamma_3, times e^{-V(w0)/t}."""
    x = tabs / (2 * pi) * np.exp(1j * argx * DEG)
    t = -2j * pi * x
    b = np.sqrt(x)
    ap = -0.5
    w0 = 2j * pi / 3
    v0 = complex(V(np.array([w0]), ap)[0])
    tot = 0j
    ce = cos(argx * DEG)
    # Gamma_1: w = i pi - rho e^{-i pi/6}, contributes int_0^inf f e^{-i pi/6} d rho
    R1 = np.sqrt(70 * tabs / (0.433 * ce)) + 1
    rr, wt = panels_osc(0, R1, 0.0, tabs)
    w = 1j * pi - rr * np.exp(-1j * pi / 6)
    tot += (np.exp(-w ** 2 / 2 / t + 1j * pi * w / t - core.lam_any(w, t, fad=fad, cap=40000) - v0 / t
                   + (w if extra_ew else 0)) * wt).sum() * np.exp(-1j * pi / 6)
    # Gamma_2: w = i eta, eta from pi to 0: -i int_0^pi f d eta
    ee, wt = panels(0, pi, int(np.clip(pi / (1.5 * tabs), 30, 3000)))
    w = 1j * ee
    tot += -1j * (np.exp(-w ** 2 / 2 / t + 1j * pi * w / t - core.lam_any(w, t, fad=fad, cap=40000) - v0 / t
                         + (w if extra_ew else 0)) * wt).sum()
    # Gamma_3: w = r, 0..inf
    R3 = 70 * tabs / (pi * ce) + 1
    rr, wt = panels_osc(0, R3, 0.0, tabs)
    w = rr + 0j
    tot += (np.exp(-w ** 2 / 2 / t + 1j * pi * w / t - core.lam_any(w, t, fad=fad, cap=40000) - v0 / t
                   + (w if extra_ew else 0)) * wt).sum()
    return tot / (2 * pi * b), t, b, v0


def I_on_L(tabs, argx, ap, extra_ew=False, delta=10.0):
    """the definition: integral over L = e^{-i delta} R (Im x > 0, moderate |t|), times e^{-V*/t}
    is not applied here (returns the raw value)."""
    x = tabs / (2 * pi) * np.exp(1j * argx * DEG)
    t = -2j * pi * x
    b = np.sqrt(x)
    rr, wt = panels(-12, 12, 3000)
    w = 2 * pi * b * rr * np.exp(-1j * delta * DEG)
    dw = 2 * pi * b * np.exp(-1j * delta * DEG)
    f = np.exp((-w ** 2 / 2 - 2j * pi * ap * w) / t - core.lam_any(w, t, cap=400000) + (w if extra_ew else 0))
    return (f * wt).sum() * dw / (2 * pi * b)


def asym_6(t, b, which):
    """Theorem 6 leading terms divided by exp(V*/t) (V* = -5pi^2/12 +- i Cl2)."""
    base = np.sqrt(t / (2 * pi)) * 3 ** -0.25 / b
    if which == "dec":
        return base * np.exp(1j * pi / 12)
    if which == "dec0":
        return np.exp(-2j * pi / 3) * base * np.exp(1j * pi / 12)
    if which == "grow":
        return base * np.exp(-1j * pi / 12)
    if which == "grow0":
        return np.exp(2j * pi / 3) * base * np.exp(-1j * pi / 12)


def part_c():
    print("\n(c) Theorem 6: (exact/asymptotic - 1)/|t|")
    # cross-check of the contour quadrature against the definition on L (Im x > 0)
    for argx, tabs in [(30, 0.6), (45, 0.5), (11, 0.8)]:
        Id, t, b, v1 = I_decay(tabs, argx, 0.5)
        IL = I_on_L(tabs, argx, 0.5)
        Ig, t, b, v0 = I_grow(tabs, argx)
        IL2 = I_on_L(tabs, argx, -0.5)
        print("   check vs. L at arg x=%d |t|=%.1f: decaying %.1e, growing %.1e (relative)"
              % (argx, tabs, abs(Id * np.exp(v1 / t) / IL - 1), abs(Ig * np.exp(v0 / t) / IL2 - 1)))
    for which, args in [("grow", [0, 2, 5, 8, 11, 11.3]), ("grow0", [0, 5, 11]),
                        ("dec", [0, 10, 30, 45, 60]), ("dec0", [0, 30, 60])]:
        line = "   %-5s" % which
        worst = 0
        for argx in args:
            sub = "  argx=%4.1f:" % argx
            for tabs in [0.5, 0.2, 0.1, 0.05, 0.02, 0.01]:
                t = -2j * pi * tabs / (2 * pi) * np.exp(1j * argx * DEG)
                fad = core.Faddeev(t) if argx < 1e-9 else None
                if which.startswith("grow"):
                    I, t, b, v = I_grow(tabs, argx, which == "grow0", fad)
                else:
                    I, t, b, v = I_decay(tabs, argx, 0.5, which == "dec0", fad)
                rr = abs(I / asym_6(t, b, which) - 1) / tabs
                worst = max(worst, rr)
                sub += " %.3f" % rr
            line += sub
        print(line)
        print("         (|t| = 0.5, 0.2, 0.1, 0.05, 0.02, 0.01); max = %.3f" % worst)
    print("   beyond: growing integral at arg x = 13, 16, 20 deg (the proof's (L1) on Gamma holds below 17.15 deg):")
    for argx in [13.0, 16.0, 20.0]:
        line = "    argx=%.0f:" % argx
        for tabs in [0.2, 0.1, 0.05, 0.02]:
            I, t, b, v = I_grow(tabs, argx)
            line += "  |t|=%.2f: %.3f" % (tabs, abs(I / asym_6(t, b, "grow") - 1) / tabs)
        print(line)


def part_d():
    print("\n(d) Theorem 7: epsilon(a) and (exact/asymptotic - 1)/|t| at real x and at arg x = epsilon(a)/2")
    print("   %7s %9s %9s %10s %10s %12s %s" % ("a", "c0", "M", "eps_suff", "eps_true", "", "ratio at |t| = 0.2, 0.1, 0.05, 0.02"))
    print("   (eps_suff = atan(c0/M), 90 if M = 0; eps_true = largest tau with (L1) on the line, searched up to tau = 50, i.e. 88.85 deg)")
    aa = sorted(list(rng.uniform(0.02, 0.98, 10)) + [0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99])
    for a in aa:
        c0, M, ts, tt = eps_of_a(a)
        es, et = atan(ts) / DEG if np.isfinite(ts) else 90.0, atan(tt) / DEG
        W = w_ = W_of(a)
        line = "   %7.4f %9.4f %9.4f %10.2f %10.2f   " % (a, c0, M, es, et)
        e1 = np.log(W).imag
        if min(e1 + 2 * pi * a, -(2 * e1 + 2 * pi * a)) < 0.05:
            # the integrand decays like e^{-margin |r|/|t|} while its phase turns like r^2/|t|:
            # too slow for this quadrature (a near 0 or 1)
            print(line + " (integral not computed: valley margin %.3f)" % min(e1 + 2 * pi * a, -(2 * e1 + 2 * pi * a)))
            continue
        for argx in [0.0, min(es, 60) / 2]:
            vals = []
            for tabs in [0.2, 0.1, 0.05, 0.02]:
                t = -2j * pi * tabs / (2 * pi) * np.exp(1j * argx * DEG)
                fad = core.Faddeev(t) if argx < 1e-9 else None
                I, t, b, v1 = I_decay(tabs, argx, a, False, fad)
                asym = (1 / (2 * pi * b)) * np.sqrt(2 * pi * t * (1 + W) / (1 + 2 * W))
                vals.append(abs(I / asym - 1) / tabs)
            line += " [argx=%.1f: %s]" % (argx, " ".join("%.2f" % v for v in vals))
        print(line)


def part_e():
    print("\n(e) the forms in section 0.2 (Theorem III): P.I(-+1/(2b)) against -+A e^{+-5 pi i/12} t^{-1/2} e^{-Sbar/t, -S/t},"
          " and q^{1/2} P.I(b +- 1/(2b)) / P.I(+-1/(2b)) against e^{-+2 pi i/3}: (ratio - 1)/|t|")
    S = pi ** 2 / 3 - 1j * CL2
    A = np.sqrt(pi / np.sqrt(3))
    for argx in [0.0, 5.0, 11.0]:
        line = "   arg x=%4.1f:" % argx
        for tabs in [0.2, 0.05, 0.02]:
            x = tabs / (2 * pi) * np.exp(1j * argx * DEG)
            t = -2j * pi * x
            b = np.sqrt(x)
            P = 2 ** -0.5 * np.exp(pi ** 2 / (12 * t) + t / 6) / b
            fad = core.Faddeev(t)
            Ig, _, _, v0 = I_grow(tabs, argx, False, fad)
            Ig0, _, _, _ = I_grow(tabs, argx, True, fad)
            Id, _, _, v1 = I_decay(tabs, argx, 0.5, False, fad)
            Id0, _, _, _ = I_decay(tabs, argx, 0.5, True, fad)
            # compare in logarithmic form to avoid overflow: P I e^{-V/t} vs formula e^{-V/t}
            lg = np.log(P) + np.log(Ig) + v0 / t
            fg = np.log(-A * np.exp(5j * pi / 12) / np.sqrt(t)) - np.conj(S) / t
            ld = np.log(P) + np.log(Id) + v1 / t
            fd = np.log(A * np.exp(-5j * pi / 12) / np.sqrt(t)) - S / t
            r1 = abs(np.exp(lg - fg) - 1) / tabs
            r2 = abs(np.exp(ld - fd) - 1) / tabs
            r3 = abs(np.exp(-t / 2) * Ig0 / Ig / np.exp(2j * pi / 3) - 1) / tabs
            r4 = abs(np.exp(-t / 2) * Id0 / Id / np.exp(-2j * pi / 3) - 1) / tabs
            line += "  |t|=%.2f: %.3f %.3f %.3f %.3f" % (tabs, r1, r2, r3, r4)
        print(line)


if __name__ == "__main__":
    T = time.time()
    part_a()
    part_b()
    part_c()
    part_d()
    part_e()
    print("\nrun time %.0f s" % (time.time() - T))
