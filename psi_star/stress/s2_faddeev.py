"""S2: Part 2, Faddeev's dilogarithm [v16 §1-§3, v19 §2].

 (a) Lemma Phi: Faddeev's integral against both product formulas; the two functional
     equations and the inversion relation at real b; Phi_b(0)^2.
 (b) Lemma S: defect / |t|^{2K-1}, K = 2, 3, |Im w| <= pi - eps1, |arg x| <= 60, |x| <= 1/4.
 (c) Lemma F: sup |Lambda + Li2(-e^w)/t| e^{-Re w/2} on Re w <= 0, |Im w| <= pi (edges incl.).
 (d) Lemma D: the margin 0.383 and the actual size of E on the contour.
 (e) Lemma B (i)-(iii): the constants, 0 <= arg x <= 60.
 (f) accuracy of the double-precision Lambda against mpmath at 40 digits.
Run:  python3 s2_faddeev.py     (about 6 min)
"""
import time
import numpy as np
import mpmath as mp
from math import pi, log, cos, sin
from scipy.optimize import minimize
import core

DEG = pi / 180
rng = np.random.default_rng(2)


def t_of(xabs, argx_deg):
    x = xabs * np.exp(1j * argx_deg * DEG)
    return -2j * pi * x


def lam(w, t, fad=None):
    """Lambda as Faddeev's integral itself (the branch the lemmas are about)."""
    if fad is None:
        fad = core.Faddeev(t)
    return fad.lam(w)


# ------------------------------------------------------------------ (f) reference
def lam_mp(w, t, dps=40):
    """log Phi_b by the product formula in mpmath (Im x != 0)."""
    mp.mp.dps = dps
    t = mp.mpc(t)
    w = mp.mpc(w)
    x = 1j * t / (2 * mp.pi)
    if t.real > 0:
        q, p = mp.exp(-t), mp.exp(-2 * mp.pi ** 2 / t)
        A = mp.qp(-mp.sqrt(q) * mp.exp(w), q, maxterms=10 ** 7)
        B = mp.qp(-p * mp.exp(w / x), p ** 2, maxterms=10 ** 7)
        return complex(mp.log(A) - mp.log(B))
    q1, p1 = mp.exp(t), mp.exp(2 * mp.pi ** 2 / t)
    A = mp.qp(-p1 * mp.exp(w / x), p1 ** 2, maxterms=10 ** 7)
    B = mp.qp(-mp.sqrt(q1) * mp.exp(w), q1, maxterms=10 ** 7)
    return complex(mp.log(A) - mp.log(B))


def part_f():
    print("(f) accuracy of the double-precision Lambda (products, Faddeev) against mpmath (40 digits)")
    worst_p, worst_f = 0, 0
    for tabs, th in [(0.3, -60), (0.1, -120), (0.05, -45), (0.03, -100), (0.5, -150)]:
        t = tabs * np.exp(1j * th * DEG)
        fad = core.Faddeev(t)
        for w in [-0.3 + 1j, -2 - 2.5j, -0.01 + 3.1j, -5 + 0.2j]:
            ref = lam_mp(w, t)
            a = core.lam_prod(np.array([w]), t, nmax_cap=10 ** 7)[0]
            b = fad.lam(np.array([w]))[0]
            worst_p = max(worst_p, abs(np.exp(a - ref) - 1))
            worst_f = max(worst_f, abs(np.exp(b - ref) - 1))
    print("   max |e^{Lambda - Lambda_ref} - 1|: products %.1e, Faddeev %.1e" % (worst_p, worst_f))


# ------------------------------------------------------------------ (a)
def part_a():
    print("\n(a) Lemma Phi")
    worst = 0
    n = 0
    for k in range(60):
        xabs = float(np.exp(rng.uniform(log(0.01), log(0.25))))
        argx = rng.uniform(-89, 89) if k % 3 else rng.choice([-60, 60, -89, 89, -30, 30])
        t = t_of(xabs, argx)
        if abs(t.real) < 2e-3:
            continue
        x = 1j * t / (2 * pi)
        if (1 / x).real < 0.85:
            continue
        fad = core.Faddeev(t)
        w = rng.uniform(-5, 0, 6) + 1j * rng.uniform(-pi, pi, 6)
        a = fad.strip(w)
        b = core.lam_prod(w, t, nmax_cap=10 ** 6)
        if b is None:
            continue
        worst = max(worst, np.abs(np.exp(a - b) - 1).max())
        n += 6
    print("   Faddeev = product formula (Im x > 0) / Phi_{1/b} product (Im x < 0): %d points, "
          "max |ratio - 1| = %.1e" % (n, worst))
    we = [0, 0, 0, 0]
    for k in range(40):
        X = rng.uniform(4, 100) if k > 3 else [4.0, 4.0001, 100.0, 5.5][k]
        x = 1 / X
        t = -2j * pi * x
        fad = core.Faddeev(t)
        w = rng.uniform(-3, 3, 5) + 1j * rng.uniform(-pi + 0.2, pi - 0.2, 5)
        # direction b: Phi(w + t) = (1 + e^{w + t/2}) Phi(w)   [z - ib/2 -> z + ib/2]
        we[0] = max(we[0], np.abs(np.exp(fad.lam(w + t) - fad.lam(w)) / (1 + np.exp(w + t / 2)) - 1).max())
        # inversion
        we[2] = max(we[2], np.abs(np.exp(fad.lam(w) + fad.lam(-w) - (-t / 24 + pi ** 2 / (6 * t) + w ** 2 / (2 * t))) - 1).max())
        # Phi_b(0)^2 = e^{-t/24 + pi^2/(6t)}
        we[3] = max(we[3], abs(np.exp(2 * fad.strip(np.array([0j]))[0] - (-t / 24 + pi ** 2 / (6 * t))) - 1))
    # direction 1/b on the edge of convergence: Faddeev's integral converges for |Im w| < pi (1 + x);
    # compare Lambda(w) for pi < Im w < pi(1+x)... with the shifted value
    w1 = 0
    for X in [4.0, 6.0, 10.0]:
        x = 1 / X
        t = -2j * pi * x
        fad = core.Faddeev(t)
        for wr in [-0.5, -2.0]:
            for a in [0.2, 0.5, 0.8]:
                w = wr + 1j * (pi + a * pi * x)          # inside the convergence strip, outside |Im w| <= pi
                direct = fad.strip(np.array([w]))[0]      # integral itself
                shifted = fad.lam(np.array([w]))[0]       # via the functional equation from w - 2 pi i
                w1 = max(w1, abs(np.exp(direct - shifted) - 1))
    print("   real x: functional eq. (b): %.1e; (1/b) [integral vs functional equation, pi<Im w<pi(1+x)]: %.1e;"
          " inversion: %.1e; Phi_b(0)^2: %.1e" % (we[0], w1, we[2], we[3]))


# ------------------------------------------------------------------ (b) Lemma S
B2 = {0: 1.0, 1: -1 / 12, 2: 7 / 240, 3: -31 / 1344}   # B_{2j}(1/2)
FACT = {0: 1, 1: 2, 2: 24, 3: 720}


def li_neg(n, z):
    """Li_{-n}(z), n >= 0, closed forms."""
    if n == 0:
        return z / (1 - z)
    if n == 2:
        return z * (1 + z) / (1 - z) ** 3
    if n == 4:
        return z * (1 + z) * (1 + 10 * z + z * z) / (1 - z) ** 5
    raise ValueError


def lemma_s_series(w, t, K):
    z = -np.exp(w)
    s = core.li2(z) / t
    for j in range(1, K):
        s = s + B2[j] * t ** (2 * j - 1) * li_neg(2 * j - 2, z) / FACT[j]
    return s


def part_b():
    print("\n(b) Lemma S: sup |defect| / |t|^{2K-1};  defect = -Lambda - sum_{j<K} ...")
    for K in (2, 3):
        for eps1 in (pi / 2, pi / 6, 0.1):
            line = "   K=%d eps1=%.3f:" % (K, eps1)
            for tabs in [1.5, 0.8, 0.4, 0.2, 0.1, 0.05]:
                if K == 3 and tabs < 0.1:
                    continue
                best = (0, None)
                for argx in [-60, -45, -30, -10, 0, 10, 30, 45, 60]:
                    t = tabs * np.exp(1j * (argx - 90) * DEG)
                    fad = core.Faddeev(t)
                    im = np.concatenate([np.linspace(-(pi - eps1), pi - eps1, 15), [pi - eps1, -(pi - eps1)]])
                    re = np.concatenate([-np.geomspace(1e-3, 40, 25), [0.0], np.geomspace(1e-3, 4, 8)])
                    W = (re[:, None] + 1j * im[None, :]).ravel()
                    lw = lam(W, t, fad)
                    d = np.abs(-lw - lemma_s_series(W, t, K)) / tabs ** (2 * K - 1)
                    i = int(np.argmax(d))
                    if d[i] > best[0]:
                        best = (d[i], (argx, W[i]))
                line += "  |t|=%.2f:%.3g" % (tabs, best[0])
            print(line + "   (worst at arg x=%d, w=%.2f%+.2fi)" % (best[1][0], best[1][1].real, best[1][1].imag))


# ------------------------------------------------------------------ (c) Lemma F
def lemma_f_val(w, t, fad):
    return np.abs(lam(w, t, fad) + core.li2_neg_exp(w) / t) * np.exp(-w.real / 2)


def part_c():
    print("\n(c) Lemma F: sup |Lambda + Li2(-e^w)/t| e^{-Re w/2}, Re w <= 0, |Im w| <= pi")
    glob = (0, None)
    for tabs in [1.5708, 1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01]:
        best = (0, None)
        for argx in [-60, -59.9, -45, -20, -5, 0, 5, 20, 45, 59.9, 60]:
            t = tabs * np.exp(1j * (argx - 90) * DEG)
            fad = core.Faddeev(t)
            im = np.concatenate([np.linspace(-pi, pi, 41), [pi - 1e-9, -pi + 1e-9]])
            re = np.concatenate([-np.geomspace(1e-4, 40, 30), [0.0]])
            W = (re[:, None] + 1j * im[None, :]).ravel()
            v = lemma_f_val(W, t, fad)
            i = int(np.argmax(v))
            if v[i] > best[0]:
                best = (v[i], (argx, W[i]))
        # local optimisation from the worst point
        argx, w0 = best[1]
        t = tabs * np.exp(1j * (argx - 90) * DEG)
        fad = core.Faddeev(t)

        def f(v):
            w = complex(min(v[0], 0), np.clip(v[1], -pi, pi))
            return -lemma_f_val(np.array([w]), t, fad)[0]
        res = minimize(f, [w0.real, w0.imag], method="Nelder-Mead", options={"xatol": 1e-9})
        val = -res.fun
        if val > best[0]:
            best = (val, (argx, complex(min(res.x[0], 0), np.clip(res.x[1], -pi, pi))))
        print("   |t|=%.4f (|x|=%.4f): sup = %.5f at arg x=%.1f, w=%.4f%+.4fi"
              % (tabs, tabs / (2 * pi), best[0], best[1][0], best[1][1].real, best[1][1].imag))
        if best[0] > glob[0]:
            glob = (best[0], (tabs, *best[1]))
    print("   overall sup = %.5f  (doc: 0.347; note ln(2)/2 = %.5f)" % (glob[0], log(2) / 2))
    line = "   by direction (|t| = 0.01 and 1.5708):"
    for argx in [-60, -30, -10, 0, 10, 30, 60]:
        vals = []
        for tabs in [0.01, 1.5708]:
            t = tabs * np.exp(1j * (argx - 90) * DEG)
            fad = core.Faddeev(t)
            W = (np.concatenate([-np.geomspace(1e-4, 40, 30), [0.0]])[:, None]
                 + 1j * np.linspace(-pi, pi, 81)[None, :]).ravel()
            vals.append(lemma_f_val(W, t, fad).max())
        line += "  %d: %.4f/%.4f" % (argx, *vals)
    print(line)
    # beyond: |arg x| > 60 and |x| > 1/4
    line = "   beyond: "
    for tabs, argx in [(1.5708, 75), (1.5708, 85), (0.1, 80), (0.1, 89), (2.5, 0), (3.0, 30), (0.1, -85)]:
        t = tabs * np.exp(1j * (argx - 90) * DEG)
        x = 1j * t / (2 * pi)
        if (1 / x).real < 0.85:
            line += " (|t|=%.2f,arg x=%d: line Im sigma=1/2 not admissible)" % (tabs, argx)
            continue
        fad = core.Faddeev(t)
        W = (np.concatenate([-np.geomspace(1e-4, 40, 30), [0.0]])[:, None] + 1j * np.linspace(-pi, pi, 41)[None, :]).ravel()
        line += " (|t|=%.2f,arg x=%d): %.4f;" % (tabs, argx, lemma_f_val(W, t, fad).max())
    print(line)


# ------------------------------------------------------------------ (d) Lemma D
def part_d():
    print("\n(d) Lemma D")
    b = np.linspace(-30, 30, 60001)
    g, a = b + 10, b - 10
    br = np.cos(g * DEG) - np.abs(np.sin(a * DEG) * np.sin(2 * b * DEG))
    print("   min over |beta| <= 30 of cos(gamma) - |sin(alpha) sin(2 beta)| = %.6f at beta = %.3f (claim >= 0.383)"
          % (br.min(), b[np.argmin(br)]))
    print("   min rho with N != 0 on the left half: pi/|sin alpha|max = %.4f (claim >= 4.88)"
          % (pi / np.sin(40 * DEG)))
    print("   actual |E| = |Lambda(w) - Lambda(w - 2 pi i N)| on the left half (products, Im x != 0):")
    for tabs in [0.4, 0.2, 0.1, 0.05]:
        worst = (0, None)
        for th in [-150, -140, -120, -100, -95, -85, -80, -60, -40, -30]:
            beta, al, ga = core.angles(th)
            t = tabs * np.exp(1j * th * DEG)
            rho = np.linspace(0.1, 60, 1200)
            w = -rho * np.exp(1j * al * DEG)
            N = np.round(w.imag / (2 * pi)).astype(int)
            m = N != 0
            if not m.any():
                continue
            lp = core.lam_prod(w[m], t, nmax_cap=10 ** 6)
            lp0 = core.lam_prod(w[m] - 2j * pi * N[m], t, nmax_cap=10 ** 6)
            if lp is None or lp0 is None:
                continue
            E = np.abs(np.exp(lp - lp0) - 1)
            i = int(np.argmax(E))
            if E[i] > worst[0]:
                worst = (E[i], (th, rho[m][i], N[m][i]))
        if worst[1] is None:
            print("   |t|=%.2f: no point with N != 0" % tabs)
            continue
        print("   |t|=%.2f: max |e^E - 1| = %.2e  (theta=%d, rho=%.2f, N=%d);  |t| log(max) = %.2f"
              % (tabs, worst[0], *worst[1], tabs * log(max(worst[0], 1e-300))))


# ------------------------------------------------------------------ (e) Lemma B
def part_e():
    print("\n(e) Lemma B, 0 <= arg x <= 60, |x| <= 1/4")
    print("   (i) sup over the closed strip, all Re w, of log(1/|Phi_b|) - Re(Li2(-e^w)/t):")
    line = "     "
    gi = -1e9
    for tabs in [1.5708, 0.5, 0.1, 0.02]:
        best = -1e9
        for argx in [0, 10, 30, 50, 60]:
            t = tabs * np.exp(1j * (argx - 90) * DEG)
            fad = core.Faddeev(t)
            re = np.concatenate([-np.geomspace(1e-3, 40, 25), [0.0], np.geomspace(1e-3, 40, 25)])
            W = (re[:, None] + 1j * np.linspace(-pi, pi, 31)[None, :]).ravel()
            v = (-lam(W, t, fad) - core.li2_neg_exp(W) / t).real
            best = max(best, v.max())
        gi = max(gi, best)
        line += " |t|=%.2f: %.4f;" % (tabs, best)
    print(line + "  => C = e^{%.4f}" % gi)
    print("   (ii) sup |(1/Phi_b) e^{-Li2/t} - 1| / |t| on |Im w| <= pi - eps1:")
    for eps1 in (pi / 6, pi / 2):
        line = "     eps1=%.3f:" % eps1
        for tabs in [1.5708, 0.5, 0.1, 0.02]:
            best = 0
            for argx in [0, 10, 30, 50, 60]:
                t = tabs * np.exp(1j * (argx - 90) * DEG)
                fad = core.Faddeev(t)
                re = np.concatenate([-np.geomspace(1e-3, 40, 25), [0.0], np.geomspace(1e-3, 40, 25)])
                W = (re[:, None] + 1j * np.linspace(-(pi - eps1), pi - eps1, 21)[None, :]).ravel()
                v = np.abs(np.exp(-lam(W, t, fad) - core.li2_neg_exp(W) / t) - 1) / tabs
                best = max(best, v.max())
            line += " |t|=%.2f: %.4f;" % (tabs, best)
        print(line)
    print("   (iii) ray w = i pi - rho e^{-i delta1}: sup of log(1/|Phi_b|) - Re(Li2/t) - N ln 2, N = 1 + rho/(2pi):")
    for d1 in [30.0, 60.0, 80.0]:
        line = "     delta1=%.0f:" % d1
        for tabs in [1.5708, 0.5, 0.1, 0.05]:
            best = -1e9
            for argx in [0, 5, 10, 30, 60]:
                if d1 + argx > 90:
                    continue
                t = tabs * np.exp(1j * (argx - 90) * DEG)
                fad = core.Faddeev(t)
                rho = np.concatenate([[0.0], np.geomspace(1e-3, 60, 200)])
                W = 1j * pi - rho * np.exp(-1j * d1 * DEG)
                v = (-lam(W, t, fad) - core.li2_neg_exp(W) / t).real - (1 + rho / (2 * pi)) * log(2)
                best = max(best, v.max())
            line += " |t|=%.2f: %.4f;" % (tabs, best)
        print(line)
    print("     without the 2^N factor (is it needed?): sup of log(1/|Phi_b|) - Re(Li2/t), delta1=30:")
    line = "      "
    for tabs in [1.5708, 0.5, 0.1, 0.05]:
        best = -1e9
        for argx in [0, 10, 30, 60]:
            t = tabs * np.exp(1j * (argx - 90) * DEG)
            fad = core.Faddeev(t)
            rho = np.concatenate([[0.0], np.geomspace(1e-3, 60, 200)])
            W = 1j * pi - rho * np.exp(-1j * 30 * DEG)
            v = (-lam(W, t, fad) - core.li2_neg_exp(W) / t).real
            best = max(best, v.max())
        line += " |t|=%.2f: %.4f;" % (tabs, best)
    print(line)


if __name__ == "__main__":
    T = time.time()
    part_f()
    part_a()
    part_b()
    part_c()
    part_d()
    part_e()
    print("\nrun time %.0f s" % (time.time() - T))
