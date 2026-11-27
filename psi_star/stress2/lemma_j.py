"""J: Lemmas U and J of v23.

 (a) Lemma U: |Phi_b| = 1 on R for real b in (0, 1]: max |Re log Phi_b(y)| (Faddeev's integral itself);
 (b) the contour turn: sup of log|Phi_b| and of log|1/Phi_b| in the left wedge pi - delta <= arg y <= pi
     (delta = 10 deg as in v23, and 45 deg); the right wedge through the inversion relation;
 (c) |P.J(a/b)| <= 1/(2 sqrt2 sin(pi a)) for a in (0, 1), real x in (0, 1]; sharpness as a -> 0, 1,
     comparison with |1/(1 - 2z)|; J computed over the line e^{-45 deg i} R (equivalent to L) and,
     as a check, over e^{-10 deg i} R and over the real line itself.
Run:  python3 lemma_j.py   (about 6 min)
"""
import time
import numpy as np
from math import pi, sin, log, sqrt
import common2 as c
import core

DEG = pi / 180
rng = np.random.default_rng(23)


def xs():
    X = [1.0, 1.25, 1.5, 2, 3, 4, 6, 8, 12, 16, 24, 40, 100, 200]
    return [1 / v for v in X] + list(rng.uniform(0.005, 1.0, 6))


def part_a():
    print("(a) Lemma U: max over real y of |Re log Phi_b(y)| (Faddeev's integral, w = 2 pi b y in [-60, 60])")
    worst = (0, None)
    for x in xs():
        t = -2j * pi * x
        fad = core.Faddeev(t)
        w = np.concatenate([np.linspace(-60, 60, 601), rng.uniform(-60, 60, 100)]) + 0j
        L = fad.lam(w)          # Re w > 0 through the inversion relation (whose terms are imaginary here)
        L2 = fad.strip(w[np.abs(w) <= 10])       # the integral itself where it is well conditioned
        m = max(np.abs(L.real).max(), np.abs(L2.real).max())
        if m > worst[0]:
            worst = (m, x)
    print("   %d values of x in (0, 1]: max |Re log Phi_b| = %.1e (at x = %.4f)" % (len(xs()), *worst))


def part_b():
    print("\n(b) the left wedge pi - delta <= arg y <= pi: sup log|Phi_b| and sup log|1/Phi_b| (|w| <= 200)")
    for delta in [10.0, 45.0]:
        line = "   delta=%2d:" % delta
        g1 = g2 = -1e9
        for x in xs():
            t = -2j * pi * x
            fad = core.Faddeev(t)
            rho = np.geomspace(1e-3, 200, 400)
            th = np.linspace(pi - delta * DEG, pi, 13)
            w = (rho[:, None] * np.exp(1j * th[None, :])).ravel()
            L = c.lam(w, t, fad).real
            fin = np.isfinite(L)
            L, ww = L[fin], w[fin]
            if L.max() > g1:
                g1, x1, w1 = L.max(), x, ww[np.argmax(L)]
            if (-L).max() > g2:
                g2, x2, w2 = (-L).max(), x, ww[np.argmax(-L)]
        print(line + " sup log|Phi_b| = %.4f (x=%.3f, w=%.2f%+.2fi), sup log|1/Phi_b| = %.4f (x=%.3f, w=%.2f%+.2fi)"
              % (g1, x1, w1.real, w1.imag, g2, x2, w2.real, w2.imag))
    print("   right wedge -delta <= arg y <= 0, delta = 45: sup of log|1/Phi_b| + Re(pi i y^2) (= log|Phi_b(-y)| - log|Phi_b(0)^2|):")
    g = -1e9
    for x in xs():
        t = -2j * pi * x
        fad = core.Faddeev(t)
        rho = np.geomspace(1e-3, 200, 400)
        th = np.linspace(-45 * DEG, 0, 13)
        w = (rho[:, None] * np.exp(1j * th[None, :])).ravel()
        val = (-c.lam(w, t, fad) + w ** 2 / (2 * t)).real      # e^{-pi i y^2} = e^{-w^2/(2t)}
        g = max(g, val[np.isfinite(val)].max())
    print("     %.4f" % g)


def PJ(a, x, delta=45.0, fad=None):
    t = -2j * pi * x
    b = np.sqrt(x)
    J, cond = c.line_integral(x, a, kernel=True, delta=delta, fad=fad)
    return c.Pfac(t, b) * J, cond


def PJ_real_line(a, x, fad):
    """J over the real line itself (Lemma J's form), in s = y/b: b int e^{-pi i x s^2 + 2 pi a s}/(Phi_b 2..)"""
    t = -2j * pi * x
    b = np.sqrt(x)
    S = 40 / (2 * pi * min(a, 1 - a)) + 5
    s, wt = c.panels_osc(-S, S, 0.0, 1.0, base=0.3)
    w = 2 * pi * x * s + 0j                 # w = 2 pi b y = 2 pi b^2 s
    f = np.exp(-1j * pi * x * s ** 2 + 2 * pi * a * s - fad.lam(w) - core.log1pexp(2 * pi * s + 0j))
    return c.Pfac(t, b) * b * (f * wt).sum()


def part_c():
    print("\n(c) |P.J(a/b)| * 2 sqrt2 sin(pi a) (claim <= 1) and |P.J(a/b)| / |1/(1 - 2z)|")
    print("   check of the contour turn: |J(line 45) / J(line 10) - 1| and |J(45) / J(real line) - 1|:")
    for a, x in [(0.25, 1 / 8), (0.5, 1 / 40), (0.7, 0.5), (0.05, 1 / 40)]:
        t = -2j * pi * x
        fad = core.Faddeev(t)
        j45, c45 = PJ(a, x, 45.0, fad)
        j10, c10 = PJ(a, x, 10.0, fad)
        jr = PJ_real_line(a, x, fad)
        print("     a=%.2f x=%.4f: %.1e (cond %.0e), %.1e" % (a, x, abs(j45 / j10 - 1), c10, abs(j45 / jr - 1)))
    A = np.unique(np.concatenate([np.geomspace(1e-3, 0.5, 12), 1 - np.geomspace(1e-3, 0.5, 12)]))
    worst = (0, None)
    print("   %8s " % "a" + " ".join("%8s" % ("X=%g" % (1 / x)) for x in xs()[:14]) + "   lim |1/(1-2z)|*2sqrt2 sin")
    for a in A:
        line = "   %8.4f " % a
        z = np.exp(2j * pi * a)
        for x in xs():
            t = -2j * pi * x
            fad = core.Faddeev(t)
            try:
                v, cond = PJ(a, x, 45.0, fad)
            except RuntimeError:
                if x in xs()[:14]:
                    line += "%8s " % "*"
                continue
            r = abs(v) * 2 * sqrt(2) * sin(pi * a)
            if r > worst[0]:
                worst = (r, (a, x))
            if x in xs()[:14]:
                line += "%8.4f " % r
        print(line + "  %8.4f" % (abs(1 / (1 - 2 * z)) * 2 * sqrt(2) * sin(pi * a)))
    print("   max ratio %.4f at a = %.4f, x = %.4f (claim <= 1)" % (worst[0], *worst[1]))


if __name__ == "__main__":
    T = time.time()
    part_a()
    part_b()
    part_c()
    print("\nrun time %.0f s" % (time.time() - T))
