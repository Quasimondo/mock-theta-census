"""S1: Part 1, the kernel identity [v12].  Random x with |q| in [0.1, 0.9], random admissible
delta, random complex u with |u| <= 3.  All integrals by independent quadrature with the
product formula for Phi_b (Im x > 0).

 (a) Theorem C: I(u) on L against c e^{-pi i u^2/2} D(u);  B1 (both recursions) by quadrature.
 (b) Fourier transforms G(w) (over L) and F(w) (over C) on the half-plane H.
 (c) steps (X) and (Y): I(u) and I_1/2(v) against K(lambda) computed by quadrature.
 (d) Theorem A (both routes) with J(u) by quadrature, the residue sums Psi, Psi-check,
     and the dual sums with I from Theorem C (after (a));  (B2), (B3).
 Integrals whose integrand exceeds |result| by more than 1e8 on the contour are skipped
 (double precision), and counted.
Run:  python3 s1_kernel.py     (about 4 min)
"""
import time
import numpy as np
from math import pi, log, cos, sin
import core

DEG = pi / 180
rng = np.random.default_rng(11)
gl_x, gl_w = np.polynomial.legendre.leggauss(40)


def panels(a, b, n):
    e = np.linspace(a, b, n + 1)
    X, Wt = [], []
    for i in range(n):
        m, h = (e[i] + e[i + 1]) / 2, (e[i + 1] - e[i]) / 2
        X.append(m + h * gl_x)
        Wt.append(h * gl_w)
    return np.concatenate(X), np.concatenate(Wt)


def cond_ok(m, v, lim=1e8):
    """the integral is trusted if log max|integrand| - log|value| < log lim (and both finite)."""
    return bool(np.isfinite(m) and np.isfinite(v) and v > 0 and m - log(v) < log(lim))


class Setup:
    def __init__(self, x, delta_deg):
        self.x = complex(x)
        self.t = -2j * pi * self.x
        self.b = np.sqrt(self.x)
        self.q = np.exp(-self.t)
        self.p = np.exp(-2 * pi ** 2 / self.t)
        self.delta = delta_deg * DEG
        self.cb = 1j * (self.b + 1 / self.b) / 2

    def logPhi(self, y):
        """log Phi_b(y) by the product formula (y in the variable of v12)."""
        w = 2 * pi * self.b * np.asarray(y, dtype=complex)
        return core.lam_prod(w, self.t, nmax_cap=10 ** 6)

    def qpoch(self, a, q, n=None):
        """(a; q)_n (n=None: infinite), |q| < 1."""
        if n is None:
            n = int(np.ceil(42 / -log(abs(q)))) + 5
        k = np.arange(n)
        return np.prod(1 - a * q ** k)

    # ----- integrals over L
    def line_integral(self, logf, R=None, h=None, ulin=3.0):
        """int_L e^{logf(y)} dy, L = e^{-i delta} R, trapezoid on a truncated line.
        Returns (value, log max |integrand|); log max = +inf if the ends are not negligible."""
        s2 = sin(2 * self.delta)
        if R is None:
            R = 2 * ulin / s2 + np.sqrt(60 / (pi * s2)) + 2
        R = min(R, 80.0)
        r = np.arange(-R, R + 1e-12, 0.004)
        y = r * np.exp(-1j * self.delta)
        lf = logf(y)
        mx = lf.real.max()
        if max(lf.real[0], lf.real[-1]) > mx - 40:
            return np.nan, np.inf
        val = np.exp(lf).sum() * 0.004 * np.exp(-1j * self.delta)
        return val, mx

    def I_quad(self, u):
        lf = lambda y: -1j * pi * y ** 2 + 2 * pi * u * y - self.logPhi(y)
        return self.line_integral(lf)

    def J_quad(self, u):
        lf = lambda y: (-1j * pi * y ** 2 + 2 * pi * u * y - self.logPhi(y)
                        - core.log1pexp(2 * pi * y / self.b))
        return self.line_integral(lf)

    # ----- Theorem C
    def D(self, u):
        b, q, p, t = self.b, self.q, self.p, self.t
        # ranges: the k-factor has modulus e^{-pi Im(u b) k - Re(t) k^2/4}, the n-factor
        # e^{-pi Im(u/b) n + Re(log p)(n^2/2 + n)}; go 50 e-folds below the maximum of each
        def nterms(lin, quad):
            kk = np.arange(0, 4000)
            e = lin * kk - quad * kk ** 2
            m = e.max()
            last = np.nonzero(e > m - 50)[0][-1]
            return int(last) + 3
        K = nterms(-pi * (u * b).imag, t.real / 4)
        lpr = (-2 * pi ** 2 / t).real
        nmax = nterms(-pi * (u / b).imag + lpr, -lpr / 2)
        k = np.arange(K)
        n = np.arange(nmax)
        lqk = np.concatenate([[0], np.cumsum(np.log(1 - q ** np.arange(1, K)))])         # log (q;q)_k
        lpn = np.concatenate([[0], np.cumsum(np.log(1 - p ** (2 * np.arange(1, nmax))))])  # log (p^2;p^2)_n
        lp = -2 * pi ** 2 / t
        A = (1j * pi * u * b * k[:, None] + 1j * pi * u * n[None, :] / b - t * k[:, None] ** 2 / 4
             + lp * (n[None, :] ** 2 / 2 + n[None, :]) - lqk[:, None] - lpn[None, :])
        sg = (-1.0) ** (k[:, None] * n[None, :] + n[None, :])
        if A.real.max() > 700:
            self.D_cond = np.inf
            self.max_cond = np.inf
            return np.nan
        val = (sg * np.exp(A)).sum()
        self.D_cond = np.exp(A.real.max()) / max(abs(val), 1e-300)
        self.max_cond = max(getattr(self, "max_cond", 0.0), self.D_cond)
        return val

    def c_const(self):
        t = self.t
        return 2 ** -0.5 * np.exp(-1j * pi / 4) * np.exp(t / 24) * np.exp(-pi ** 2 / (6 * t))

    def I_C(self, u):
        return self.c_const() * np.exp(-1j * pi * u ** 2 / 2) * self.D(u)


def random_setup():
    """|q| in [0.1, 0.9], arg x <= 120 deg (beta <= 60), delta in [0.25, 0.75] (90 - beta)."""
    while True:
        absq = rng.uniform(0.1, 0.9)
        imx = -log(absq) / (2 * pi)
        rex = rng.uniform(-0.5, 0.5)
        x = rex + 1j * imx
        beta = np.angle(np.sqrt(x)) / DEG
        if beta <= 60:
            break
    delta = rng.uniform(0.25, 0.75) * (90 - beta)
    return Setup(x, delta), absq


# ------------------------------------------------------------------ (a)
def part_a():
    print("(a) Theorem C (I(u) = c e^{-pi i u^2/2} D(u)) and (B1), by quadrature on L")
    worstC, worstB, n, skip = 0, 0, 0, 0
    loc = None
    while n < 150:
        S, absq = random_setup()
        u = 3 * np.sqrt(rng.random()) * np.exp(2j * pi * rng.random())
        Iq, mx = S.I_quad(u)
        Ic = S.I_C(u)
        if not np.isfinite(mx) or mx - log(abs(Iq) + 1e-300) > log(1e8) or S.D_cond > 1e8:
            skip += 1
            if skip > 2000:
                break
            continue
        e = abs(Iq / Ic - 1)
        if e > worstC:
            worstC, loc = e, (S.x, S.delta / DEG, u)
        # (B1) by quadrature: I(u) = p^{-1} e^{2 pi i u/b} [I(u + 1/b) + p^{-1} I(u + 2/b)]
        I1, m1 = S.I_quad(u + 1 / S.b)
        I2, m2 = S.I_quad(u + 2 / S.b)
        rhs = np.exp(2j * pi * u / S.b) / S.p * (I1 + I2 / S.p)
        J1, n1 = S.I_quad(u + S.b)
        J2, n2 = S.I_quad(u + 2 * S.b)
        rhs2 = np.sqrt(S.q) * np.exp(2j * pi * u * S.b) * (J1 + np.sqrt(S.q) * J2)
        # residual relative to the size of the terms on the right (p^{-1} can be huge)
        sc1 = abs(np.exp(2j * pi * u / S.b) / S.p) * (abs(I1) + abs(I2 / S.p)) + abs(Iq)
        sc2 = abs(np.sqrt(S.q) * np.exp(2j * pi * u * S.b)) * (abs(J1) + abs(np.sqrt(S.q) * J2)) + abs(Iq)
        if all(cond_ok(m, v) for m, v in ((m1, abs(I1)), (m2, abs(I2)), (n1, abs(J1)), (n2, abs(J2)))):
            worstB = max(worstB, abs(rhs - Iq) / sc1, abs(rhs2 - Iq) / sc2)
        n += 1
    print("   %d samples (%d skipped as ill-conditioned): max |I_quad/I_C - 1| = %.2e at x=%s delta=%.1f u=%s"
          % (n, skip, worstC, loc[0], loc[1], loc[2]))
    print("   (B1) both recursions: max residual / (sum of |terms|) = %.2e" % worstB)


# ------------------------------------------------------------------ (b)
def part_b():
    print("\n(b) Fourier transforms on H = {Im(w e^{-i delta}) < 0}")
    wG, wF, n = 0, 0, 0
    while n < 80:
        S, absq = random_setup()
        d = S.delta
        # w in H: w = e^{i delta} (a - i c), c > 0
        w = np.exp(1j * d) * (rng.uniform(-1.5, 1.5) - 1j * rng.uniform(0.05, 1.5))
        # G(w) = int_L e^{2 pi i w y}/Phi_b(y) dy = (i/b)((q~)_inf/(q)_inf) e^{-2 pi i c_b w} Phi_b(w + c_b)
        lf = lambda y: 2j * pi * w * y - S.logPhi(y)
        Gq, mx = S.line_integral(lf, R=30.0)
        qq = S.qpoch(S.q, S.q)
        pp = S.qpoch(S.p ** 2, S.p ** 2)
        Gc = (1j / S.b) * (pp / qq) * np.exp(-2j * pi * S.cb * w) * np.exp(S.logPhi(np.array([w + S.cb]))[0])
        if cond_ok(mx, abs(Gq)):
            wG = max(wG, abs(Gq / Gc - 1))
        # F(w) = int_C Phi_b(y) e^{2 pi i w y} dy = -i b ((q)_inf/(q~)_inf) e^{2 pi i c_b w} / Phi_b(-w - c_b)
        r, wt = panels(0, 30, 600)
        yR = r * np.exp(1j * d)
        yL = -r * np.exp(-1j * d)
        fR = np.exp(S.logPhi(yR) + 2j * pi * w * yR) * np.exp(1j * d)
        fL = np.exp(S.logPhi(yL) + 2j * pi * w * yL) * np.exp(-1j * d)
        Fq = ((fR + fL) * wt).sum()
        Fc = -1j * S.b * (qq / pp) * np.exp(2j * pi * S.cb * w) / np.exp(S.logPhi(np.array([-w - S.cb]))[0])
        ends = max(abs(fR[-1]), abs(fL[-1])) < 1e-16 * max(np.abs(fR).max(), np.abs(fL).max())
        cond = max(np.abs(fR).max(), np.abs(fL).max()) / abs(Fq)
        if ends and cond < 1e8:
            wF = max(wF, abs(Fq / Fc - 1))
        n += 1
    print("   %d samples: max |G_quad/G_formula - 1| = %.2e,  max |F_quad/F_formula - 1| = %.2e" % (n, wG, wF))


# ------------------------------------------------------------------ (c)
def K_quad(S, lam):
    """K(lambda) = int Phi_b(z) e^{pi i z^2 + 2 pi lambda z} dz over e^{i delta} R (through 0, below the poles)."""
    r = np.arange(-14, 14 + 1e-12, 0.004)
    z = r * np.exp(1j * S.delta)
    lf = S.logPhi(z) + 1j * pi * z ** 2 + 2 * pi * lam * z
    mx = lf.real.max()
    if max(lf.real[0], lf.real[-1]) > mx - 40:
        return np.nan, np.inf
    return np.exp(lf).sum() * 0.004 * np.exp(1j * S.delta), mx


def part_c():
    print("\n(c) steps (X) and (Y): the integrals against K(lambda) by quadrature")
    wX, wY, n, skip = 0, 0, 0, 0
    while n < 80 and skip < 3000:
        S, absq = random_setup()
        u = 1.5 * np.sqrt(rng.random()) * np.exp(2j * pi * rng.random())
        cb = S.cb
        qq = S.qpoch(S.q, S.q)
        pp = S.qpoch(S.p ** 2, S.p ** 2)
        Iq, m1 = S.I_quad(u)
        Kq, m2 = K_quad(S, -u - 2j * cb)
        rhsX = np.exp(-1j * pi / 4) * (1j / S.b) * (pp / qq) * np.exp(1j * pi * (1j * u - cb) ** 2 + 2j * pi * cb ** 2) * Kq
        # (Y): I_1/2(v) over C, v = u/2 + i c_b/2
        v = u / 2 + 1j * cb / 2
        r, wt = panels(0, 30, 600)
        yR = r * np.exp(1j * S.delta)
        yL = -r * np.exp(-1j * S.delta)
        fR = np.exp(S.logPhi(yR) - 1j * pi * yR ** 2 / 2 + 2 * pi * v * yR) * np.exp(1j * S.delta)
        fL = np.exp(S.logPhi(yL) - 1j * pi * yL ** 2 / 2 + 2 * pi * v * yL) * np.exp(-1j * S.delta)
        Ih = ((fR + fL) * wt).sum()
        Phi0sq = np.exp(-S.t / 24 + pi ** 2 / (6 * S.t))
        Kq2, m3 = K_quad(S, -2 * v - 1j * cb)
        rhsY = (np.sqrt(2) * np.exp(-1j * pi / 4) * (-1j * S.b) * (qq / pp) / Phi0sq
                * np.exp(2j * pi * (1j * v - cb) ** 2 - 2j * pi * cb ** 2) * Kq2)
        okX = cond_ok(m1, abs(Iq)) and cond_ok(m2, abs(Kq))
        endsY = max(abs(fR[-1]), abs(fL[-1])) < 1e-16 * max(np.abs(fR).max(), np.abs(fL).max())
        okY = endsY and max(np.abs(fR).max(), np.abs(fL).max()) / abs(Ih) < 1e8 and cond_ok(m3, abs(Kq2))
        if not (okX and okY):
            skip += 1
            continue
        wX = max(wX, abs(rhsX / Iq - 1))
        wY = max(wY, abs(rhsY / Ih - 1))
        n += 1
    print("   %d samples (%d skipped): max relative error (X): %.2e, (Y): %.2e" % (n, skip, wX, wY))


# ------------------------------------------------------------------ (d)
def Psi(S, z, check=False):
    q = S.q
    s, P, n = 0j, 1 + 0j, 0
    while True:
        if n > 0:
            P = P * (1 + q ** n) if not check else P / (1 + q ** n)
        term = (q ** (n * (n + 1) // 2) if not check else q ** (n * n + n)) * z ** n * P
        s += term
        n += 1
        if n > 5 and abs(term) < 1e-18 * (abs(s) + 1):
            break
        if n > 5000:
            break
    return s


def part_d():
    print("\n(d) Theorem A (both routes), (B2), (B3); dual sums with I(u -+ m/b) from Theorem C")
    wA1, wA2, wB2a, wB2b, wB3, n, skip = 0, 0, 0, 0, 0, 0, 0
    while n < 120 and skip < 3000:
        S, absq = random_setup()
        u = 3 * np.sqrt(rng.random()) * np.exp(2j * pi * rng.random())
        b, q, p = S.b, S.q, S.p
        Jq, mx = S.J_quad(u)
        if not cond_ok(mx, abs(Jq)):
            skip += 1
            continue
        S.max_cond = 0.0
        pinf = S.qpoch(p, p ** 2)
        mq = np.prod(1 + q ** np.arange(1, int(42 / -log(abs(q))) + 5))
        # dual sums
        s_up, s_dn = 0j, 0j
        Iu = {}
        for m in range(0, 200):
            a = S.I_C(u - m / b) if m >= 1 else S.I_C(u)
            Iu[-m] = a
            if m >= 1:
                s_up += (-1) ** m * a
            bm = S.I_C(u + m / b)
            Iu[m] = bm
            s_dn += (-1) ** m * bm
            if m > 8 and abs(a) + abs(bm) < 1e-17 * (abs(s_up) + abs(s_dn) + abs(Jq)):
                break
        A_up = -1j * b * q ** 0.125 * np.exp(1j * pi * u * b) * (pinf / mq) * Psi(S, np.exp(2j * pi * u * b)) - s_up
        A_dn = 1j * b * q ** 0.125 * np.exp(-1j * pi * u * b) * (pinf / (2 * mq)) * Psi(S, np.exp(-2j * pi * u * b), True) + s_dn
        if not S.max_cond < 1e8:
            skip += 1
            continue
        # the residue sum and the dual sum can be ~1e20 and cancel to |J|: measure against their size
        res_up = A_up + s_up
        res_dn = A_dn - s_dn
        big1 = abs(res_up) + abs(s_up)
        big2 = abs(res_dn) + abs(s_dn)
        if big1 < 1e3 * abs(Jq):
            wA1 = max(wA1, abs(A_up / Jq - 1))
        if big2 < 1e3 * abs(Jq):
            wA2 = max(wA2, abs(A_dn / Jq - 1))
        wA1b = max(globals().get("wA1b", 0), abs(A_up - Jq) / big1)
        wA2b = max(globals().get("wA2b", 0), abs(A_dn - Jq) / big2)
        globals()["wA1b"], globals()["wA2b"] = wA1b, wA2b
        # (B2)
        tB = p * np.exp(2j * pi * u / b)
        def Wcal(T):
            s, P, k = 0j, 1 + 0j, 0
            while True:
                if k > 0:
                    P *= 1 - p ** (2 * k - 1)
                term = (-T) ** k * p ** (k * k + k) * P
                s += term
                k += 1
                if k > 5 and abs(term) < 1e-18:
                    break
            return s
        def Gcal(T, which):
            s, P, k = 0j, 1 + 0j, 0
            while True:
                if which == 0:
                    if k > 0:
                        P *= 1 - p ** (2 * k - 1)
                    term = T ** (-k) * p ** (2 * k * k + k) / P
                else:
                    P *= 1 - p ** (2 * k + 1)
                    term = T ** (-k) * p ** (2 * k * k + 3 * k) / P
                s += term
                k += 1
                if k > 5 and abs(term) < 1e-18 * (abs(s) + 1):
                    break
                if k > 3000:
                    break
            return s
        rhs1 = Iu[0] * (Wcal(tB / p ** 2) - 1) - (tB / p) * S.I_C(u + 1 / b) * Wcal(tB)
        rhs2 = Iu[0] * Gcal(tB, 0) - S.I_C(u + 1 / b) * Gcal(tB, 1)
        sc = abs(s_up) + abs(Iu[0])
        if not S.max_cond < 1e8:
            n += 1
            continue
        wB2a = max(wB2a, abs(rhs1 - s_up) / sc)
        wB2b = max(wB2b, abs(rhs2 - s_dn) / (abs(s_dn) + abs(Iu[0])))
        # (B3) at u = 1/(2b) + k b
        k = int(rng.integers(-2, 3))
        u3 = 1 / (2 * b) + k * b
        tt = p * np.exp(2j * pi * u3 / b)
        def dsum(uu, sign, m0):
            tot = 0j
            for m in range(m0, 200):
                v = S.I_C(uu + sign * m / b)
                tot += (-1) ** m * v
                if m > 8 and abs(v) < 1e-17 * (abs(tot) + 1e-300):
                    break
            return tot
        sup3 = dsum(u3, -1, 1)
        sdn3 = dsum(u3, +1, 0)
        W1 = Wcal(1.0)
        W3 = Wcal(p ** 2)
        Fa = Gfa = sum(p ** (2 * nn * nn - nn) / np.prod(1 - p ** (2 * np.arange(1, nn + 1) - 1)) for nn in range(1, 40))
        Fb = sum(p ** (2 * nn * nn - 3 * nn + 1) / np.prod(1 - p ** (2 * np.arange(1, nn + 1) - 1)) for nn in range(1, 40))
        r3a = S.I_C(u3) * p * W3 - S.I_C(u3 - 1 / b) * W1
        if not S.max_cond < 1e8:
            n += 1
            continue
        r3b = S.I_C(u3) * Fb - S.I_C(u3 - 1 / b) * Fa
        wB3 = max(wB3, abs(tt - 1), abs(r3a - sup3) / (abs(sup3) + abs(S.I_C(u3))), abs(r3b - sdn3) / (abs(sdn3) + abs(S.I_C(u3))))
        n += 1
    print("   %d samples (%d skipped): (A up) %.2e, (A down) %.2e relative to |J| where the two sums cancel by < 1e3;"
          " relative to |residue sum| + |dual sum| in all samples: %.2e, %.2e;" % (n, skip, wA1, wA2,
                                                                                   globals()["wA1b"], globals()["wA2b"]))
    print("   (B2) first %.2e, second %.2e, (B3) %.2e" % (wB2a, wB2b, wB3))


def qseries(S):
    q = S.q
    out = {}
    s1 = s0 = f1 = f0 = 0j
    P = 1 + 0j
    n = 0
    while True:
        if n > 0:
            P *= 1 + q ** n
        a = (-1) ** n * q ** (n * (n + 1) // 2) * P
        s1 += a
        s0 += a * q ** (n + 1)
        c = (-1) ** n * q ** (n * n) / P
        f0 += c
        f1 += c * q ** n
        n += 1
        if n > 5 and abs(a) + abs(c) < 1e-18:
            break
    return s1, s0, f1, f0


def dual(S):
    p = S.p
    W = W3 = Fa = Fb = 0j
    P = 1 + 0j
    j = 0
    while True:
        if j > 0:
            P *= 1 - p ** (2 * j - 1)
        tw = (-1) ** j * p ** (j * j + j) * P
        W += tw
        W3 += tw * p ** (2 * j)
        if j > 0:
            Fa += p ** (2 * j * j - j) / P
            Fb += p ** (2 * j * j - 3 * j + 1) / P
        j += 1
        if j > 5 and abs(tw) + abs(p ** (2 * j * j - 3 * j)) < 1e-18:
            break
    return W, W3, Fa, Fb


def part_e():
    import s5_factorial as S5
    print("\n(e) Theorem I (the four exact formulas) with all integrals by quadrature on L,"
          " and the link P.J(1/(2b)) [v12, over L] = e^{t/8} int H ds [Parts 3-5, over the s-contour]")
    worst = [0, 0, 0, 0]
    wl = [0, 0]
    wP = 0
    n = 0
    skip = 0
    while n < 60 and skip < 2000:
        S, absq = random_setup()
        b, q, p, t = S.b, S.q, S.p, S.t
        P = 2 ** -0.5 * np.exp(pi ** 2 / (12 * t) + t / 6) / b
        mq = np.prod(1 + q ** np.arange(1, int(42 / -log(abs(q))) + 5))
        P2 = mq / (b * q ** 0.125 * S.qpoch(p, p ** 2))
        wP = max(wP, abs(P / P2 - 1))
        J1, m1 = S.J_quad(1 / (2 * b))
        Ip, m2 = S.I_quad(1 / (2 * b))
        Im_, m3 = S.I_quad(-1 / (2 * b))
        J0, m4 = S.J_quad(b + 1 / (2 * b))
        Ip0, m5 = S.I_quad(b + 1 / (2 * b))
        Im0, m6 = S.I_quad(b - 1 / (2 * b))
        vals = [abs(J1), abs(Ip), abs(Im_), abs(J0), abs(Ip0), abs(Im0)]
        if not all(cond_ok(m, v) for m, v in zip([m1, m2, m3, m4, m5, m6], vals)):
            skip += 1
            continue
        s1, s0, f1, f0 = qseries(S)
        W, W3, Fa, Fb = dual(S)
        qh = np.sqrt(q)
        r = [P * (J1 + Ip * p * W3 - Im_ * W),
             2 * P * (J1 - Ip * Fb + Im_ * Fa),
             qh * P * (J0 + Ip0 * p * W3 - Im0 * W),
             2 * qh * P * (J0 - Ip0 * Fb + Im0 * Fa)]
        for k, (lhs, rhs) in enumerate(zip([s1, f1, s0, f0], r)):
            worst[k] = max(worst[k], abs(lhs - rhs) / max(abs(lhs), 1e-300))
        # link to the s-contour integral of Parts 3-5 (direction theta of t)
        th = np.angle(t) / DEG
        wl[0] = max(wl[0], abs(P * J1 / S5.pj_direct(t, 1, th=th) - 1))
        wl[1] = max(wl[1], abs(qh * P * J0 / S5.pj_direct(t, 0, th=th) - 1))
        n += 1
    print("   %d samples (%d skipped): psi1* %.1e, f1* %.1e, psi0* %.1e, f0* %.1e (relative)" % (n, skip, *worst))
    print("   P = 2^{-1/2} e^{pi^2/(12t)+t/6}/b vs (-q;q)/(b q^{1/8} (p;p^2)): %.1e;"
          "  P.J over L vs e^{t/8} int H ds: psi1 %.1e, psi0 %.1e" % (wP, *wl))


if __name__ == "__main__":
    T = time.time()
    part_a()
    part_b()
    part_c()
    part_d()
    part_e()
    print("\nrun time %.0f s" % (time.time() - T))
