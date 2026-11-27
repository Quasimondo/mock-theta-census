"""P3: what round 1 did not test, and the corrections of v26.

 (a) the closing lemma of v12 §4.1: the exponent bound on the vertical lines Re xi = sigma_N with the
     stated constant c, and the remark "the supremum on the line is this bound times a factor between
     1 and 1.5";
 (b) Lemma L as an abstract statement, on a model family (entire V, analytic g, Gamma = R);
 (c) Theorem 8 / Theorem IV as a whole: psi_1^* at |q| close to 1 in 0 < arg x <= 10 deg against the
     Borel sum + the two exponential terms with W(p), p W3(p);
 (d) Theorem C at multiple precision outside the round-1 sample (|q| < 0.1, |q| > 0.9, arg x > 120 deg,
     delta near 0 and near 90 deg - beta);
 (e) v26: Lemma V (the claims inside its proof), and the growth rate 1 - pi/|t| of the e^w-version of
     the growing integral on Gamma_3 (v26 §2).
Run:  python3 p3.py   (about 15 min)
"""
import os
import sys
import time
import numpy as np
import mpmath as mp
from math import pi, log, sin, cos, tan, sqrt
import common2 as c
import core

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "stress"))
import s5_factorial as S5  # noqa: E402  (round-1 coefficients a_k, own code)

DEG = pi / 180
rng = np.random.default_rng(33)


# ------------------------------------------------------------------ (a) closing lemma
def part_a():
    print("(a) closing lemma: on Re xi = sigma_N (xi = y/b), tau from the point on C upwards:")
    print("    E(tau) = log|Phi_b(b xi)| + Re(-A pi i y^2);  bound -c sigma_N^2, c = pi(1-A) min(|x|^2/Im x, |x| sin 2delta/cos^2(delta-beta))")
    print("    reported: f_N = sup_tau E / (-c sigma_N^2) (the remark says 1 <= f <= 1.5 for the sup), and sup E + c sigma_N^2")
    worst_f = [1e9, -1e9]
    worst_gap = -1e9
    n = 0
    for k in range(40):
        absq = rng.uniform(0.2, 0.8)
        x = complex(rng.uniform(-0.4, 0.4), -log(absq) / (2 * pi))
        b = np.sqrt(x)
        beta = np.angle(b) / DEG
        delta = rng.uniform(0.2, 0.8) * (90 - beta)
        t = -2j * pi * x
        for A in (0.0, 0.5):
            cc = pi * (1 - A) * min(abs(x) ** 2 / x.imag, abs(x) * sin(2 * delta * DEG) / cos((delta - beta) * DEG) ** 2)
            fs = []
            for N in range(2, 9):
                sig = N * x.imag / abs(x) ** 2
                tauC = sig * tan((delta - beta) * DEG)
                tau = tauC + np.concatenate([[0.0], np.geomspace(1e-3, 30 * (1 + sig), 600)])
                xi = sig + 1j * tau
                y = b * xi
                lp = core.lam_prod(2 * pi * b * y, t, nmax_cap=10 ** 6)
                if lp is None:
                    continue
                E = lp.real + (-A * pi * 1j * y ** 2).real
                sup = E.max()
                fs.append(sup / (-cc * sig ** 2))
                worst_gap = max(worst_gap, (sup + cc * sig ** 2) / (sig + 1))
            if fs:
                n += 1
                worst_f[0] = min(worst_f[0], fs[-1])
                worst_f[1] = max(worst_f[1], fs[-1])
    print("   %d (x, delta, A) samples: f at the largest N in [%.3f, %.3f];  max of (sup E + c sigma_N^2)/(sigma_N + 1) = %.3f"
          % (n, worst_f[0], worst_f[1], worst_gap))
    print("   (f >= 1 means the bound holds at leading order; the lower-order terms are bounded by the last number)")
    print("   trend in N for three samples (f_N, N = 2, 4, 8, 16, 32):")
    for x, delta, A in [(0.1 + 0.15j, 30.0, 0.5), (-0.2 + 0.1j, 10.0, 0.0), (0.3 + 0.05j, 40.0, 0.5)]:
        b = np.sqrt(x)
        beta = np.angle(b) / DEG
        t = -2j * pi * x
        cc = pi * (1 - A) * min(abs(x) ** 2 / x.imag, abs(x) * sin(2 * delta * DEG) / cos((delta - beta) * DEG) ** 2)
        which = "|x|^2/Im x" if abs(x) ** 2 / x.imag <= abs(x) * sin(2 * delta * DEG) / cos((delta - beta) * DEG) ** 2 else "point on C"
        fs = []
        for N in [2, 4, 8, 16, 32]:
            sig = N * x.imag / abs(x) ** 2
            tauC = sig * tan((delta - beta) * DEG)
            tau = tauC + np.concatenate([[0.0], np.geomspace(1e-3, 30 * (1 + sig), 1500)])
            y = b * (sig + 1j * tau)
            lp = core.lam_prod(2 * pi * b * y, t, nmax_cap=10 ** 7)
            E = lp.real + (-A * pi * 1j * y ** 2).real
            fs.append(E.max() / (-cc * sig ** 2))
        print("     x=%s delta=%.0f A=%.1f (min attained by %s): %s" % (x, delta, A, which, " ".join("%.4f" % f for f in fs)))


# ------------------------------------------------------------------ (b) Lemma L on a model
def part_b():
    print("\n(b) Lemma L on a model: F(t) = int_R g(w) e^{V(w)/t} dw, V = -w^2/2 - w^4/12 + i w^3/10,"
          " g = e^{iw}/(2 + cos w): (F / (g(0)(2 pi t)^{1/2}) - 1)/|t|")
    gx, gw = np.polynomial.legendre.leggauss(60)
    for argt in [0.0, 30.0, 60.0, 80.0]:
        ww = np.linspace(-6, 6, 4001)
        V = -ww ** 2 / 2 - ww ** 4 / 12 + 1j * ww ** 3 / 10
        H = (V * np.exp(-1j * argt * DEG)).real
        L1 = bool(np.all(H[np.abs(ww) > 1e-9] < 0))
        line = "   arg t = %4.0f ((L1) %s):" % (argt, "holds" if L1 else "fails")
        for tabs in [0.2, 0.1, 0.05, 0.02, 0.01]:
            t = tabs * np.exp(1j * argt * DEG)
            R = 6.0
            e = np.linspace(-R, R, 241)
            ww = ((e[1:] + e[:-1])[:, None] / 2 + (e[1:] - e[:-1])[:, None] / 2 * gx[None, :]).ravel()
            wt = ((e[1:] - e[:-1])[:, None] / 2 * gw[None, :]).ravel()
            V = -ww ** 2 / 2 - ww ** 4 / 12 + 1j * ww ** 3 / 10
            g = np.exp(1j * ww) / (2 + np.cos(ww))
            F = (g * np.exp(V / t) * wt).sum()
            lead = (1 / 3) * np.sqrt(2 * pi * t)
            line += "  %.3f" % (abs(F / lead - 1) / tabs)
        print(line)


# ------------------------------------------------------------------ (c) Theorem 8
def borel_sum(a, t, th):
    """S_theta[Phi](t) = (1/t) int_0^{inf e^{i th}} e^{-xi/t} B(xi) d xi; B from the coefficients by series
    (|xi| < 1.5) and Pade beyond."""
    with mp.workdps(40):
        K = len(a)
        bc = [a[k] / mp.factorial(k) for k in range(K)]
        P, Q = mp.pade(bc, (K - 1) // 2, K - 1 - (K - 1) // 2)
        e = mp.expj(th * mp.pi / 180)
        tt = mp.mpc(t)
        rate = mp.re(e / tt)

        def B(xi):
            if abs(xi) < 1.5:
                return mp.polyval(bc[::-1], xi)
            return mp.polyval(P[::-1], xi) / mp.polyval(Q[::-1], xi)
        f = lambda r: mp.exp(-r * e / tt) * B(r * e) * e
        rmax = 90 / rate
        pts = [0, min(1, rmax)] + ([rmax] if rmax > 1 else [])
        return complex(mp.quad(f, pts) / tt)


def psi1_star(q):
    """direct summation at 40 digits (in double precision the product (-q;q)_n of several hundred factors
    with |q| close to 1 loses ~1 % at arg x = 0.3 deg, |t| = 0.05)."""
    with mp.workdps(40):
        q = mp.mpc(q)
        s, P, n = mp.mpc(0), mp.mpc(1), 0
        while True:
            if n > 0:
                P *= 1 + q ** n
            term = (-1) ** n * q ** (n * (n + 1) // 2) * P
            s += term
            n += 1
            if n > 10 and abs(term) < mp.mpf(10) ** -30 * abs(s):
                break
        return complex(s)


def Wser(p, e3=False):
    s, P, j = 0j, 1 + 0j, 0
    while True:
        if j > 0:
            P *= 1 - p ** (2 * j - 1)
        term = (-1) ** j * p ** (j * j + (3 * j if e3 else j)) * P
        s += term
        j += 1
        if j > 10 and abs(term) < 1e-17:
            break
        if j > 10 ** 6:
            break
    return s


def part_c():
    print("\n(c) Theorem 8 (psi_1^*): R = (psi1* - S[Phi1] - decaying term) / growing term - 1, reported as R/|t|")
    a = [complex(v) for v in S5.coeffs(46, 1)]
    Sb = pi ** 2 / 3 + 1j * c.core.CL2
    S = pi ** 2 / 3 - 1j * c.core.CL2
    A = sqrt(pi / sqrt(3))
    for argx in [10.0, 5.0, 2.0, 1.0, 0.3]:
        line = "   arg x=%4.1f:" % argx
        for tabs in [0.4, 0.2, 0.1, 0.05, 0.03]:
            th = argx - 90
            t = tabs * np.exp(1j * th * DEG)
            q = np.exp(-t)
            p = np.exp(-2 * pi ** 2 / t)
            psi = psi1_star(q)
            Sv = borel_sum(a, t, th)
            grow = A / np.sqrt(t) * np.exp(5j * pi / 12) * np.exp(-Sb / t) * Wser(p)
            dec = A / np.sqrt(t) * np.exp(-5j * pi / 12) * np.exp(-S / t) * p * Wser(p, True)
            R = (psi - Sv - dec) / grow - 1
            line += "  |t|=%.2f: %.3f (|p|=%.2f)" % (tabs, abs(R) / tabs, abs(p))
        print(line)


# ------------------------------------------------------------------ (d) Theorem C at mp
def part_d():
    print("\n(d) Theorem C at 30 digits outside the round-1 sample: |I_quad / (c e^{-pi i u^2/2} D(u)) - 1|")
    mp.mp.dps = 30
    cases = [(0.02 + 0.40j, 10.0, 0.7 + 0.3j, "|q| = 0.08"), (0.05 + 0.45j, 40.0, -0.5 + 0.8j, "|q| = 0.06"),
             (0.10 + 0.014j, 2.0, 0.3 - 0.2j, "|q| = 0.92, arg b ~ 4 deg, delta small"),
             (-0.30 + 0.12j, 3.0, 0.4 + 0.1j, "arg x = 158 deg (beta 79 deg), delta = 3 deg"),
             (0.2 + 0.2j, 66.0, -0.6 - 0.4j, "delta near 90 - beta = 67.5 deg"),
             (0.15 + 0.04j, 1.0, 0.2 + 0.2j, "|q| = 0.78, delta = 1 deg")]
    for x, delta, u, note in cases:
        x = mp.mpc(x)
        b = mp.sqrt(x)
        t = -2j * mp.pi * x
        q = mp.exp(-t)
        lp = -2 * mp.pi ** 2 / t
        p = mp.exp(lp)
        u = mp.mpc(u)

        def logphi(y):
            v = mp.exp(2 * mp.pi * b * y)
            w2 = mp.exp(2 * mp.pi * y / b)
            return mp.log(mp.qp(-mp.exp(-t / 2) * v, q)) - mp.log(mp.qp(-p * w2, p ** 2))

        e = mp.expj(-mp.mpf(delta) * mp.pi / 180)
        f = lambda r: mp.exp(-1j * mp.pi * (r * e) ** 2 + 2 * mp.pi * u * r * e - logphi(r * e)) * e
        s2 = mp.sin(2 * mp.mpf(delta) * mp.pi / 180)
        R = (2 * abs(u) + mp.sqrt(4 * abs(u) ** 2 + 4 * s2 * 80 / mp.pi)) / (2 * s2) + 2
        T0 = time.time()
        Iq = mp.quad(f, mp.linspace(-R, R, 41))
        # D(u)
        D = mp.mpc(0)
        for k in range(0, 400):
            lk = -t * k * k / 4 + 1j * mp.pi * u * b * k
            if k > 10 and mp.re(lk) < -80 - mp.re(lp) * 0:
                pass
            qk = mp.qp(q, q, k) if k > 0 else 1
            row = mp.mpc(0)
            for nn in range(0, 200):
                term = ((-1) ** (k * nn + nn) * mp.exp(lk + 1j * mp.pi * u * nn / b + lp * (mp.mpf(nn * nn) / 2 + nn))
                        / (qk * (mp.qp(p ** 2, p ** 2, nn) if nn > 0 else 1)))
                row += term
                if nn > 5 and abs(term) < mp.mpf(10) ** -40:
                    break
            D += row
            if k > 10 and abs(row) < mp.mpf(10) ** -40:
                break
        cst = 2 ** mp.mpf(-0.5) * mp.expj(-mp.pi / 4) * mp.exp(t / 24) * mp.exp(-mp.pi ** 2 / (6 * t))
        Ic = cst * mp.exp(-1j * mp.pi * u ** 2 / 2) * D
        print("   x=%s delta=%.1f u=%s (%s): %s  (%.0f s)" % (mp.nstr(x, 3), delta, mp.nstr(u, 2), note,
                                                         mp.nstr(abs(Iq / Ic - 1), 3), time.time() - T0))


# ------------------------------------------------------------------ (e) v26
def part_e():
    print("\n(e) v26 §3 Lemma V: g(lambda) = lambda^4 + 2c lambda^3 + lambda^2")
    worst_sol, worst_lm, lm_range = 0, 0, [9, -9]
    thr = None
    for a in np.linspace(0.0005, 0.9995, 1999):
        W = c.W_dec(a)
        e1 = np.angle(W)
        cc = cos(e1)
        lam = np.geomspace(1e-4, 1e3, 200001)
        gl = lam ** 4 + 2 * cc * lam ** 3 + lam ** 2
        nsol = int(np.sum(np.diff(np.sign(gl - 1)) != 0))
        worst_sol = max(worst_sol, nsol)
        if cc < -2 * sqrt(2) / 3:
            lmm = (3 * abs(cc) - sqrt(9 * cc * cc - 8)) / 4
            lm_range = [min(lm_range[0], lmm), max(lm_range[1], lmm)]
            worst_lm = max(worst_lm, lmm ** 4 + 2 * cc * lmm ** 3 + lmm ** 2)
            if thr is None:
                thr = a
    print("   1999 values of a: max number of solutions of g = 1: %d (claim 1); lambda_- in [%.4f, %.4f] (claim (1/2, 0.7071));"
          " max g(lambda_-) = %.4f (claim < 0.16); cos eta1 < -2sqrt2/3 first at a = %.5f (v26: 0.807674)"
          % (worst_sol, lm_range[0], lm_range[1], worst_lm, thr))
    print("   v26 §2, growing e^w-version on Gamma_3 at real x: slope of log|integrand| in r against 1 - pi/|t|")
    line = "  "
    for X in [2, 3, 4, 6]:
        t = -2j * pi / X
        fad = core.Faddeev(t)
        r = np.array([20.0, 30.0]) + 0j
        lf = ((-r ** 2 / 2 + 1j * pi * r) / t - c.lam(r, t, fad) + r).real
        line += "  X=%d: %.4f vs %.4f" % (X, (lf[1] - lf[0]) / 10, 1 - pi / abs(t))
    print(line)


if __name__ == "__main__":
    T = time.time()
    part_a()
    part_b()
    part_e()
    part_c()
    part_d()
    print("\nrun time %.0f s" % (time.time() - T))
