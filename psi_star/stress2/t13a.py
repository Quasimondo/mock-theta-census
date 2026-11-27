"""T13a: Theorem 13 (a) of v22, the growing integral I((a-1)/b) at real x.

 (a) the closed forms of the heights on C1, C2, C3, C5; V'(w0) = 0; Im V(w0) = -Vol(z);
 (b) (L1) at tau = 0 on each of C1 ... C5: sup (H - H(w0))/min(d^2, d), H = -Im V, for a dense in (0, 1)
     down to 1e-3 and up to 1 - 1e-3, with local zooming at the worst point of each piece;
 (c) the bound |e^{Q/t}/Phi_b| <= C |e^{V/t}| on each piece (log C measured), X = 4 ... 256;
 (d) v22 §5: along the horizontal line through the growing saddle -Im V rises in both directions;
 (e) Theorem 13 (a) itself: (exact/asymptotic - 1)/|t|, the exact integral over C, cross-checked
     against the integral over the line e^{-45 deg i} R (equivalent to L) where that is well conditioned.
Run:  python3 t13a.py   (about 10 min)
"""
import time
import numpy as np
import mpmath as mp
from math import pi, log
import common2 as c
import core

DEG = pi / 180


def a_grid():
    g = np.geomspace(1e-3, 0.5, 40)
    return np.unique(np.round(np.concatenate([g, 1 - g, [0.25, 0.75]]), 12))


def part_a():
    print("(a) saddle and closed forms")
    worst = [0, 0, 0, 0, 0, 0]
    for a in a_grid():
        ap = a - 1
        W = c.W_grow(a)
        w0 = np.log(W)
        worst[0] = max(worst[0], abs(-w0 - 2j * pi * ap - np.log(1 + W)))
        worst[1] = max(worst[1], abs(c.V(np.array([w0]), ap)[0].imag + c.vol(a)))
        r = -np.geomspace(1e-3, 60, 50)
        worst[2] = max(worst[2], np.abs(c.V(r + 2j * pi, ap).imag - 2 * pi * a * np.abs(r)).max())
        eta = np.linspace(pi, 2 * pi, 25)
        worst[3] = max(worst[3], max(abs(c.V(np.array([1j * e]), ap)[0].imag - float(mp.clsin(2, e - pi))) for e in eta))
        rr = w0.real * np.linspace(0, 1, 25)
        v3 = c.V(rr + 1j * pi, ap).imag
        lin = v3[-1] * np.linspace(0, 1, 25)
        worst[4] = max(worst[4], np.abs(v3 - lin).max())
        r5 = w0.real + np.geomspace(1e-3, 60, 50)
        worst[5] = max(worst[5], np.abs(c.V(r5 + 0j, ap).imag - 2 * pi * (1 - a) * r5).max())
    print("   %d values of a in [1e-3, 1 - 1e-3]: |V'(w0)| %.1e; |Im V(w0) + Vol| %.1e; C1 %.1e; C2 %.1e; "
          "C3 linear %.1e; C5 %.1e" % (len(a_grid()), *worst))


def zoom(a, name, wbest, w0, v0, tau=0.0):
    """refine the sup near the worst grid point of a piece (dense local grid)."""
    ap = a - 1
    H0 = -v0.imag + tau * v0.real
    # local direction of the piece
    if name in ("C1", "C3", "C5", "C1'"):
        dirn = 1.0
    else:
        dirn = 1j
    best = -1e9
    span = max(abs(wbest - w0), 1e-6) * 0.2 + 1e-6
    for k in range(3):
        ws = wbest + dirn * np.linspace(-span, span, 2001)
        # stay on the piece
        if name == "C1":
            ws = ws[(ws.real <= 0)]
        if name == "C2":
            ws = ws[(ws.imag >= pi) & (ws.imag <= 2 * pi)]
        if name == "C3":
            lo, hi = sorted([0.0, w0.real])
            ws = ws[(ws.real >= lo) & (ws.real <= hi)]
        if name == "C4":
            ws = ws[(ws.imag >= 0) & (ws.imag <= pi)]
        if name == "C5":
            ws = ws[ws.real >= w0.real]
        d = np.abs(ws - w0)
        ws = ws[d > 1e-9]
        d = d[d > 1e-9]
        if len(ws) == 0:
            break
        dv = c.dV(ws, w0, ap, v0)
        q = (-dv.imag + tau * dv.real) / np.minimum(d ** 2, d)
        i = int(np.argmax(q))
        best = max(best, q[i])
        wbest = ws[i]
        span /= 50
    return best


def part_b():
    print("\n(b) (L1) at tau = 0: sup over each piece of (H - H(w0))/min(d^2, d)  (claim < 0)")
    print("   %10s %9s %10s %10s %10s %10s %10s %10s %10s" % ("a", "Vol", "C1", "C2", "C3", "C4", "C5", "worst", "worst/Vol"))
    rows = []
    for a in a_grid():
        res, w0, v0 = c.heights_L1(a, 0.0, n=3000)
        vals = {}
        for name, (q, wb) in res.items():
            vals[name] = max(q, zoom(a, name, wb, w0, v0))
        wv = max(vals.values())
        V = c.vol(a)
        rows.append((a, V, vals, wv))
        print("   %10.6f %9.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.4f"
              % (a, V, vals["C1"], vals["C2"], vals["C3"], vals["C4"], vals["C5"], wv, wv / V))
    w = max(rows, key=lambda r: r[3])
    print("   largest (least negative) value %.6f at a = %.6f; all negative: %s"
          % (w[3], w[0], all(r[3] < 0 for r in rows)))
    # very close to the ends
    for a in [1e-4, 1e-5, 1 - 1e-4, 1 - 1e-5]:
        res, w0, v0 = c.heights_L1(a, 0.0, n=3000)
        wv = max(max(q, zoom(a, nm, wb, w0, v0)) for nm, (q, wb) in res.items())
        print("   beyond the grid: a = %.6f: worst %.3e, Vol = %.3e, ratio %.3f" % (a, wv, c.vol(a), wv / c.vol(a)))


def part_c():
    print("\n(c) the bound for 1/Phi_b on each piece: sup of log|e^{Q/t}/Phi_b| - Re(V/t) = -Re Lambda + Im Li2_cont/|t|")
    print("    (C1, C2: functional equation once + Lemma F at the shifted point; C3-C5: Lemma B (i))")
    for a in [0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99]:
        line = "   a=%.2f:" % a
        per = {}
        for X in [4, 8, 16, 32, 64, 128, 256]:
            t = -2j * pi / X
            fad = core.Faddeev(t)
            pcs, w0 = c.contour_C(a, n=300, Rl=40.0, Rr=40.0)
            sup = {}
            for name, w in pcs:
                L = c.lam(w, t, fad)
                B = (-L).real + c.li2_cont(w).imag / abs(t)
                B = B[np.isfinite(B)]
                sup[name] = B.max()
            for k, v in sup.items():
                per[k] = max(per.get(k, -1e9), v)
            line += " X=%d:%.3f(%s)" % (X, max(sup.values()), max(sup, key=sup.get))
        print(line)
        print("      per piece, max over X: " + "  ".join("%s %.3f" % (k, v) for k, v in per.items()))


def part_d():
    print("\n(d) v22 §5: along w0 +- s (s > 0), -Im V rises: min over s in (0, 60] of (Im V(w0) - Im V(w0 +- s))/min(s^2, s)")
    worst = (1e9, None)
    for a in np.concatenate([a_grid(), [1e-4, 1 - 1e-4]]):
        W = c.W_grow(a)
        w0 = np.log(W)
        ap = a - 1
        v0 = c.V(np.array([w0]), ap)[0]
        s = np.geomspace(1e-5, 60, 3000)
        for sg in (+1, -1):
            v = c.V(w0 + sg * s, ap)
            q = (v0.imag - v.imag) / np.minimum(s ** 2, s)
            i = int(np.argmin(q))
            if q[i] < worst[0]:
                worst = (q[i], (a, sg, s[i]))
    print("   min = %.4e at a = %.5f, direction %+d, s = %.3g (claim > 0)" % (worst[0], *worst[1]))
    for a in [0.5, 0.2]:
        W = c.W_grow(a)
        w0 = np.log(W)
        v0 = c.V(np.array([w0]), a - 1)[0]
        print("   a=%.1f: rise at s = 2, 8 to the left: %.3f, %.3f (v22: %s)"
              % (a, *(v0.imag - c.V(np.array([w0 - s]), a - 1)[0].imag for s in (2, 8)),
                 "1.2, 7.4" if a == 0.5 else "3.07, 16.4"))


def part_e():
    print("\n(e) Theorem 13 (a): (I_exact/asymptotic - 1)/|t|, I over C")
    print("   cross-check of the integral over C against the line e^{-45 deg i} R:")
    for a, X in [(0.5, 8), (0.25, 8), (0.8, 16), (0.1, 16), (0.5, 24)]:
        x = 1 / X
        t = -2j * pi * x
        fad = core.Faddeev(t)
        Ic, t, b, v0 = c.I_grow_C(a, X, fad=fad)
        Il, cond = c.line_integral(x, a - 1, fad=fad)
        print("     a=%.2f X=%d: |I_C/I_line - 1| = %.1e  (conditioning of the line integral %.1e)"
              % (a, X, abs(Ic * np.exp(v0 / t) / Il - 1), cond))
    print("   %8s %s" % ("a", "  X = 8, 16, 32, 64, 128"))
    for a in [0.02, 0.05, 0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9, 0.95, 0.98]:
        line = "   %8.3f" % a
        for X in [8, 16, 32, 64, 128]:
            tabs = 2 * pi / X
            if 40 * tabs / (2 * pi * min(a, 1 - a)) > 60:
                line += "      -   "
                continue
            t = -2j * pi / X
            fad = core.Faddeev(t)
            try:
                I, t, b, v0 = c.I_grow_C(a, X, fad=fad)
            except RuntimeError:
                line += "      *   "
                continue
            line += " %9.4f" % (abs(I / c.asym13(a, t, b) - 1) / tabs)
        print(line)
    print("   (- : decay length 40|t|/(2 pi min(a, 1-a)) > 60, not computed; * : too many panels)")


if __name__ == "__main__":
    T = time.time()
    part_a()
    part_b()
    part_d()
    part_c()
    part_e()
    print("\nrun time %.0f s" % (time.time() - T))
