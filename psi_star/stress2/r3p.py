"""R3': v25 (Lemma R3', Theorem 18).

 (a) Lemma F at real x in (1/4, 1] and in a complex neighbourhood |x/x0 - 1| <= 0.1:
     sup of |Lambda + Li2(-e^w)/t| e^{-Re w/2} over Re w <= 0, |Im w| <= pi (Faddeev's integral itself);
 (b) Lemma D as used in v25: on the left half of the contour of v16 (|beta| <= 30 deg), the factors
     1 + eta of the functional equation: |eta| <= e^{-0.37 rho}?, |E| <= 1.2?, for |t| <= 2 pi * 1.01;
 (c) Theorems 9, 11, 14 at X = 1, 2, 3 with 25-digit arithmetic: an own multiple-precision Faddeev
     integral (trapezoid rule on Im sigma = 1/2), the integrals over e^{-45 deg i} R by Gauss-Legendre.
Run:  python3 r3p.py   (about 25 min)
"""
import time
import numpy as np
import mpmath as mp
from math import pi, log, cos, sin
import common2 as c
import core

DEG = pi / 180


# ------------------------------------------------------------------ (a)
def part_a():
    print("(a) Lemma F beyond |x| <= 1/4: sup |Lambda + Li2(-e^w)/t| e^{-Re w/2}, Re w <= 0, |Im w| <= pi")
    re = np.concatenate([-np.geomspace(1e-4, 40, 40), [0.0]])
    im = np.concatenate([np.linspace(-pi, pi, 61)])
    W = (re[:, None] + 1j * im[None, :]).ravel()
    glob = (0, None)
    for x0 in [0.25, 0.3, 1 / 3, 0.4, 0.5, 2 / 3, 0.8, 0.9, 1.0]:
        line = "   x0=%.3f:" % x0
        for rho, phi in [(0, 0)] + [(r, ph) for r in (0.05, 0.1) for ph in (0, 90, 180, 270)]:
            x = x0 * (1 + rho * np.exp(1j * phi * DEG))
            if (1 / x).real <= 0.85:
                continue
            t = -2j * pi * x
            fad = core.Faddeev(t)
            v = np.abs(fad.lam(W) + core.li2_neg_exp(W) / t) * np.exp(-W.real / 2)
            m = v.max()
            if m > glob[0]:
                glob = (m, (x, W[np.argmax(v)]))
            if rho == 0:
                line += " real: %.4f;" % m
            elif phi == 0:
                line += " rho=%.2f: %.4f" % (rho, m)
            else:
                line += "/%.4f" % m
        print(line)
    print("   overall sup %.4f at x = %.4f%+.4fi, w = %.3f%+.3fi  (round 1: 0.511 at |x| = 1/4, arg x = 60 deg)"
          % (glob[0], glob[1][0].real, glob[1][0].imag, glob[1][1].real, glob[1][1].imag))


# ------------------------------------------------------------------ (b)
def part_b():
    print("\n(b) Lemma D for |t| <= 2 pi * 1.01 on the left half of the v16 contour, |beta| <= 30 deg")
    w_eta = 0
    w_E = 0
    w_Ebound = 0
    w_cnt = 0
    for tabs in [0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 2 * pi, 2 * pi * 1.01]:
        for beta in np.linspace(-30, 30, 13):
            th = 2 * beta - 90
            t = tabs * np.exp(1j * th * DEG)
            x = 1j * t / (2 * pi)
            al = beta - 10
            rho = np.linspace(0.05, 80, 4000)
            w = -rho * np.exp(1j * al * DEG)
            N = np.round(w.imag / (2 * pi)).astype(int)
            lp = -2 * pi ** 2 / t
            for i in np.nonzero(N)[0]:
                n = N[i]
                lw2 = w[i] / x
                j = np.arange(abs(n))
                e = lw2 + (2 * j + 1) * lp if n > 0 else lw2 - (2 * j + 1) * lp
                eta = np.abs(np.exp(e))
                w_eta = max(w_eta, (eta / np.exp(-0.37 * rho[i])).max())
                E = np.abs(core.log1pexp(e)).sum()
                w_E = max(w_E, E)
                w_Ebound = max(w_Ebound, E / (1.2 * rho[i] * np.exp(-0.37 * rho[i])))
                w_cnt = max(w_cnt, abs(n) / rho[i])
    print("   max |eta| / e^{-0.37 rho} = %.4f;  max sum |log(1 + eta)| = %.4f (claim <= 1.2);"
          "  max of that / (1.2 rho e^{-0.37 rho}) = %.4f;  max (number of factors)/rho = %.3f"
          % (w_eta, w_E, w_Ebound, w_cnt))


# ------------------------------------------------------------------ (c) multiple precision
class MPFaddeev:
    """Lambda(w) = log Phi_b at mp precision: trapezoid on Im sigma = 1/2 for Re w <= 0, |Im w| <= pi,
    the inversion relation for Re w > 0 and the functional equation in the direction 1/b otherwise."""

    def __init__(self, x, dps=30, h=None):
        self.dps = dps
        mp.mp.dps = dps
        self.x = mp.mpf(x) if np.isreal(x) else mp.mpc(x)
        self.t = -2j * mp.pi * self.x
        h = h if h is not None else mp.mpf(2 * pi * 0.5) / (dps * 2.303 + 8)   # e^{-2 pi (1/2)/h} ~ 10^{-dps}
        self.h = h
        self.cache = {}
        self.logp = -2 * mp.pi ** 2 / self.t

    def grid(self, A):
        key = int(A) + 1
        if key not in self.cache:
            n = int(key / self.h) + 1
            sg = [mp.mpf(k) * self.h + 0.5j for k in range(-n, n + 1)]
            x = self.x
            lD = [-mp.log(4 * s) - mp.log(mp.sinh(mp.pi * s)) - mp.log(mp.sinh(mp.pi * x * s)) for s in sg]
            self.cache[key] = (sg, lD)
        return self.cache[key]

    def strip(self, w):
        gap = mp.pi * (1 + mp.re(self.x)) - abs(mp.im(w))
        A = (self.dps * 2.303 + 10 + mp.log(1 / (4 * mp.pi ** 2 * abs(self.x)) + 1)) / gap
        sg, lD = self.grid(A)
        s = mp.fsum(mp.exp(-1j * w * sgk + l) for sgk, l in zip(sg, lD) if abs(mp.re(sgk)) <= A)
        return s * self.h

    def lam(self, w):
        w = mp.mpc(w)
        t = self.t
        if mp.re(w) > 0:
            return -self.lam(-w) - t / 24 + mp.pi ** 2 / (6 * t) + w ** 2 / (2 * t)
        N = int(mp.nint(mp.im(w) / (2 * mp.pi)))
        base = self.strip(w - 2j * mp.pi * N)
        lw2 = w / self.x
        if N > 0:
            corr = -mp.fsum(mp.log(1 + mp.exp(lw2 + (2 * j + 1) * self.logp)) for j in range(N))
        elif N < 0:
            corr = mp.fsum(mp.log(1 + mp.exp(lw2 - (2 * j + 1) * self.logp)) for j in range(-N))
        else:
            corr = 0
        return base + corr


def mp_integral(F, ap, extra, kernel, npan=40, nodes=20, R=None):
    """(1/(2 pi b)) int over w = 2 pi b r e^{-i pi/4}, r real, of e^{Q/t + extra w - Lambda} [/(1 + e^{w/x})]."""
    mp.mp.dps = F.dps
    x, t = F.x, F.t
    b = mp.sqrt(x)
    e = mp.expj(-mp.pi / 4)
    if R is None:
        # |e^{-w^2/(2t)}| along the line: decays like e^{-|w|^2 sin(pi/2)/(2|t|)} on the left
        R = mp.sqrt(2 * abs(t) * (F.dps * 2.303 + 20)) + 4 * abs(2 * mp.pi * ap) + 4 * abs(extra) * abs(t) + 4
    xs, ws = mp.gauss_quadrature(nodes, "legendre") if hasattr(mp, "gauss_quadrature") else (None, None)
    gx, gw = np.polynomial.legendre.leggauss(nodes)
    edges = [-R + 2 * R * k / npan for k in range(npan + 1)]
    tot = mp.mpc(0)
    for k in range(npan):
        a_, b_ = mp.mpf(edges[k]), mp.mpf(edges[k + 1])
        m, hw = (a_ + b_) / 2, (b_ - a_) / 2
        for xi, wi in zip(gx, gw):
            r = m + hw * mp.mpf(xi)
            w = r * e
            lf = (-w ** 2 / 2 - 2j * mp.pi * ap * w) / t - F.lam(w) + extra * w
            if kernel:
                lf -= mp.log(1 + mp.exp(w / x))
            tot += hw * mp.mpf(wi) * mp.exp(lf)
    return tot * e / (2 * mp.pi * b)


def part_c():
    print("\n(c) Theorems 9, 11, 14 at X = 1, 2, 3 with 25-digit arithmetic (own mp Faddeev integral)")
    dps = 25
    for X in [1, 2, 3]:
        T0 = time.time()
        x = mp.mpf(1) / X
        F = MPFaddeev(x, dps)
        mp.mp.dps = dps
        t = F.t
        b = mp.sqrt(F.x)
        P = 2 ** mp.mpf(-0.5) * mp.exp(mp.pi ** 2 / (12 * t) + t / 6) / b
        # check the quadrature: two resolutions of the first integral
        J1 = mp_integral(F, 0.5, 0, True)
        J1b = mp_integral(F, 0.5, 0, True, npan=56, nodes=24)
        Ip = mp_integral(F, 0.5, 0, False)
        Im = mp_integral(F, -0.5, 0, False)
        J0 = mp_integral(F, 0.5, 1, True)
        Ip0 = mp_integral(F, 0.5, 1, False)
        Im0 = mp_integral(F, -0.5, 1, False)
        qh = mp.exp(-t / 2)
        if X == 2:
            f1 = P * (J1 + Ip - Im)
            f0 = qh * P * (J0 + Ip0 - Im0)
            v1, v0 = 1, -1
        else:
            f1 = P * (J1 - Ip / 3 - Im / 3)
            f0 = qh * P * (J0 - Ip0 / 3 - Im0 / 3)
            v1 = mp.mpf(1) / 3 if X == 1 else 1
            v0 = mp.mpf(1) / 3 if X == 1 else 0
        print("   X=%d: quadrature check |J(two resolutions)| rel. %.1e;  psi1* formula - value = %s;  psi0* formula - value = %s  (%.0f s)"
              % (X, float(abs(J1 / J1b - 1)), mp.nstr(abs(f1 - v1), 3), mp.nstr(abs(f0 - v0), 3), time.time() - T0))
        if X == 1:
            print("        P.J(1/2) = %s, P.I(1/2) = %s, P.I(-1/2) = %s  (v25: -0.005723 - 0.239909i, -0.300817 - 0.306142i, -0.716352 - 0.413586i)"
                  % (mp.nstr(P * J1, 8), mp.nstr(P * Ip, 8), mp.nstr(P * Im, 8)))
        if X == 2:
            for a in [mp.mpf(1) / 4, mp.mpf(7) / 10, mp.mpf(1) / 20]:
                Ja = mp_integral(F, a, 0, True)
                Ig = mp_integral(F, a - 1, 0, False)
                Id = mp_integral(F, a, 0, False)
                z = mp.expjpi(2 * a)
                val = 1j * P * mp.expjpi(-a) * (Ja - Ig + z ** 2 * Id)
                print("        Theorem 14 at X = 2, a = %s: formula - 1 = %s" % (mp.nstr(a, 3), mp.nstr(abs(val - 1), 3)))


if __name__ == "__main__":
    T = time.time()
    part_a()
    part_b()
    part_c()
    print("\nrun time %.0f s" % (time.time() - T))
