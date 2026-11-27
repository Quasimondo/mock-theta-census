"""T14: Theorem 14 of v22 (and Theorem 18 (c) of v25 for X = 2).

 (a) F(z; q0) = i P e^{-pi i a} [J(a/b) - I((a-1)/b) + z^X I(a/b)] at x = 1/X, X = 2, 4, ..., 20, 24, 30, 36, 40, 50, four random a and a = 1/2;
     J over e^{-45 deg i} R (kernel), I((a-1)/b) over the contour C of v22 (and over the line where that is
     well conditioned), I(a/b) over the horizontal line through the decaying saddle (and over the line);
 (b) the dual series on the path x = 1/X + i eps: W(T) -> 1, W(p^2 T) -> 1, p T -> z^X;
 (c) odd X: the same right-hand side against F(z; q0) (no claim; what happens).
Run:  python3 t14.py   (about 12 min)
"""
import time
import numpy as np
import mpmath as mp
from math import pi, log
import common2 as c
import core

DEG = pi / 180
rng = np.random.default_rng(14)


def I_dec_line(a, X, fad):
    """I(a/b) at real x = 1/X over the horizontal line through w1 = Log W_dec (v19 §4), times 1."""
    x = 1 / X
    t = -2j * pi * x
    tabs = abs(t)
    b = np.sqrt(x)
    W = c.W_dec(a)
    w1 = np.log(W)
    e1 = w1.imag
    v1 = complex(-w1 ** 2 / 2 - 2j * pi * a * w1 + core.li2(-W))
    rate = min(e1 + 2 * pi * a, -(2 * e1 + 2 * pi * a))
    R = 70 * tabs / rate + 3
    r, wt = c.panels_osc(w1.real - R, w1.real + R, w1.real, tabs)
    w = r + 1j * e1
    f = np.exp((-w ** 2 / 2 - 2j * pi * a * w) / t - c.lam(w, t, fad) - v1 / t)
    return (f * wt).sum() / (2 * pi * b) * np.exp(v1 / t)


def rhs(a, X, use_line=False):
    x = 1 / X
    t = -2j * pi * x
    b = np.sqrt(x)
    fad = core.Faddeev(t)
    J, cj = c.line_integral(x, a, kernel=True, fad=fad)
    if use_line:
        Ig, cg = c.line_integral(x, a - 1, fad=fad)
        Id, cd = c.line_integral(x, a, fad=fad)
    else:
        cg = cd = 1.0
        try:
            Ig, _, _, v0 = c.I_grow_C(a, X, fad=fad)
            Ig = Ig * np.exp(v0 / t)
        except RuntimeError:            # very long C1 or C5 (a near 0 or 1): the 45-deg line instead
            Ig, cg = c.line_integral(x, a - 1, fad=fad)
        try:
            Id = I_dec_line(a, X, fad)
        except RuntimeError:            # slowly decaying horizontal line (a near 0 or 1)
            Id, cd = c.line_integral(x, a, fad=fad)
    z = np.exp(2j * pi * a)
    P = c.Pfac(t, b)
    val = 1j * P * np.exp(-1j * pi * a) * (J - Ig + z ** X * Id)
    return val, max(cj, cg, cd), abs(P * J)


def part_a():
    print("(a) Theorem 14: |F(z; q0) - formula| / |F|  [J over the 45-deg line; I over C and the horizontal line]")
    worst = (0, None)
    allA = []
    for X in list(range(2, 21, 2)) + [24, 30, 36, 40, 50]:
        A = list(rng.uniform(0.02, 0.98, 4)) + [0.5]
        line = "   X=%2d:" % X
        for a in A:
            F = complex(c.F_finite(a, X))
            try:
                v, cond, pj = rhs(a, X)
            except RuntimeError:
                line += "  a=%.3f: *" % a
                continue
            e = abs(F - v) / abs(F)
            allA.append(e)
            if e > worst[0]:
                worst = (e, (X, a))
            line += "  a=%.3f: %.0e" % (a, e)
        print(line)
    print("   worst %.2e at X=%d, a=%.4f; median %.1e" % (worst[0], *worst[1], np.median(allA)))
    print("   the same with all three integrals over the 45-deg line (small X), and the conditioning:")
    for X in [2, 4, 6, 8, 12]:
        line = "   X=%2d:" % X
        for a in [0.1, 0.3, 0.5, 0.8, 0.95]:
            F = complex(c.F_finite(a, X))
            v, cond, pj = rhs(a, X, use_line=True)
            line += "  a=%.2f: %.0e (cond %.0e)" % (a, abs(F - v) / abs(F), cond)
        print(line)


def W_series(T, p):
    s, P, k = mp.mpc(0), mp.mpc(1), 0
    while True:
        if k > 0:
            P *= 1 - p ** (2 * k - 1)
        term = (-T) ** k * p ** (k * k + k) * P
        s += term
        k += 1
        if k > 10 and abs(term) < mp.mpf(10) ** -25:
            break
        if k > 200000:
            break
    return s


def part_b():
    print("\n(b) dual series on x = 1/X + i eps: |W(T) - 1|, |W(p^2 T) - 1|, |p T - z^X|")
    mp.mp.dps = 30
    for X in [4, 8, 20]:
        for a in [0.1, 0.5, 0.9]:
            line = "   X=%2d a=%.1f:" % (X, a)
            for eps in [1e-2, 1e-3, 1e-4]:
                x = mp.mpf(1) / X + 1j * mp.mpf(eps)
                p = mp.exp(-1j * mp.pi / x)
                T = p * mp.exp(2j * mp.pi * a / x)
                z = mp.expjpi(2 * mp.mpf(a))
                w1 = W_series(T, p)
                w2 = W_series(p ** 2 * T, p)
                line += "  eps=%.0e: %.1e %.1e %.1e" % (eps, abs(w1 - 1), abs(w2 - 1), abs(p * T - z ** X))
            print(line)


def part_c():
    print("\n(c) odd X (no claim): |F(z; q0) - formula| / |F| and |F|")
    for X in [3, 5, 7, 9, 15, 21]:
        line = "   X=%2d:" % X
        for a in [0.1, 0.3, 0.5, 0.7]:
            F = complex(c.F_finite(a, X))
            try:
                v, cond, pj = rhs(a, X)
                line += "  a=%.1f: %.1e (|F| %.2e)" % (a, abs(F - v) / abs(F), abs(F))
            except RuntimeError:
                line += "  a=%.1f: *" % a
        print(line)


if __name__ == "__main__":
    T = time.time()
    part_a()
    part_b()
    part_c()
    print("\nrun time %.0f s" % (time.time() - T))
