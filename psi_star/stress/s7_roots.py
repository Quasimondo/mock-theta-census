"""S7: Part 7, roots of unity [v20, v21].

 (a) Lemma R1: block bound |B_m| <= min(1 + C delta, (m+1) C delta) and its two steps, random X, m, delta;
     radial limit of psi_1^* at even roots.
 (b) Lemma R2: each inequality of the proof; the range of delta where the bracket (K = 3/2) is > 0;
     W(p) -> 1, W3(p) -> 1 in the region |1 - p| <= (3/2)(1 - |p|).
 (c) Lemma R2': along the path to p = -1 (X odd) and slightly off it: phi <= delta^2, the lower
     bounds 9/8 and e^{-0.3/delta}, the bound for G(Q), |P_+-|, |psi(p)|, and the limits.
 (d) Lemma R3: continuity of P.J, P.I(+-1/(2b)) at x = 1/X as eps -> 0.
 (e) Theorems 9 and 11 (and 10, 12): the three-integral formulas at x = 1/X, X = 2 ... 40,
     including X = 3 (excluded) and X = 2, 4 (edge), for psi_1^* and psi_0^*.
Run:  python3 s7_roots.py     (about 10 min)
"""
import time
import numpy as np
import mpmath as mp
from math import pi, log, cos, sin, atan, exp, sqrt
import core
import s5_factorial as S5
import s6_state as S6

DEG = pi / 180
rng = np.random.default_rng(7)
CL2 = core.CL2


# ------------------------------------------------------------------ (a)
def part_a():
    print("(a) Lemma R1")
    w1 = w2 = w3 = 0
    loc = None
    for k in range(4000):
        X = int(rng.choice(np.arange(4, 42, 2)))
        m = int(rng.integers(0, 300)) if k % 4 else int(rng.integers(0, 3))
        delta = float(np.exp(rng.uniform(log(1e-7), log(1.0 / X))))
        rho = 1 - delta
        q0 = np.exp(2j * pi / X)
        j = np.arange(1, X + 1)
        cj = rho ** (m * X + j)
        Bm = np.prod(1 + cj * q0 ** j)
        c = rho ** ((m + 1) * X)
        C = X ** 2 * (1 + 3.0 ** (X - 1))
        r1 = abs(Bm) / min(1 + C * delta, (m + 1) * C * delta)
        r2 = abs(Bm - (1 - c ** X)) / (X ** 2 * 3.0 ** (X - 1) * delta)
        r3 = (1 - c ** X) / min(1, (m + 1) * X ** 2 * delta)
        if r1 > w1:
            w1, loc = r1, (X, m, delta)
        w2, w3 = max(w2, r2), max(w3, r3)
    print("   4000 samples: max |B_m|/min(1+C d,(m+1)C d) = %.3e at X=%d m=%d delta=%.2e;"
          " max |B_m-(1-c^X)|/(X^2 3^{X-1} d) = %.3e; max (1-c^X)/min(1,(m+1)X^2 d) = %.4f"
          % (w1, *loc, w2, w3))
    # sharper: what C would do?  sup |B_m| - 1 over delta and m, divided by delta X^2
    sharp = 0
    for X in [4, 6, 10, 20, 40]:
        q0 = np.exp(2j * pi / X)
        j = np.arange(1, X + 1)
        for delta in [1e-6, 1e-4, 1e-2 / X]:
            for m in range(0, 2000, 7):
                cj = (1 - delta) ** (m * X + j)
                sharp = max(sharp, (abs(np.prod(1 + cj * q0 ** j)) - 1) / (delta * X ** 2))
    print("   sup (|B_m| - 1)/(X^2 delta) over X <= 40 = %.4f  (so C ~ X^2 suffices; doc C = X^2 (1 + 3^{X-1}))" % sharp)
    print("   radial limit: |psi_1^*(rho q0) - sum_{n<X/2}| for delta = 1-rho:")
    mp.mp.dps = 30
    for X in [4, 6, 10, 20]:
        q0 = mp.expjpi(mp.mpf(2) / X)
        fin = mp.mpc(0)
        P = mp.mpc(1)
        for n in range(X // 2):
            if n > 0:
                P *= 1 + q0 ** n
            fin += (-1) ** n * q0 ** (n * (n + 1) // 2) * P
        line = "    X=%2d:" % X
        for delta in [1e-2, 1e-3, 1e-4]:
            q = (1 - mp.mpf(delta)) * q0
            s = mp.mpc(0)
            P = mp.mpc(1)
            n = 0
            while True:
                if n > 0:
                    P *= 1 + q ** n
                term = (-1) ** n * q ** (n * (n + 1) // 2) * P
                s += term
                n += 1
                if n > 20 and abs(term) < mp.mpf(10) ** -25:
                    break
            line += "  d=%.0e: %.2e" % (delta, abs(s - fin))
        print(line)


# ------------------------------------------------------------------ (b)
def W_series(p, which=1, dps=30):
    mp.mp.dps = dps
    p = mp.mpc(p)
    s = mp.mpc(0)
    P = mp.mpc(1)
    j = 0
    while True:
        if j > 0:
            P *= 1 - p ** (2 * j - 1)
        e = j * j + j if which == 1 else j * j + 3 * j
        term = (-1) ** j * p ** e * P
        s += term
        j += 1
        if j > 10 and abs(term) < mp.mpf(10) ** (-dps + 5):
            break
        if j > 200000:
            break
    return complex(s)


def part_b():
    print("\n(b) Lemma R2")
    K = 1.5
    f = lambda d: 1 / (2 * K) - np.exp(-1 / K + 2 * d) / (2 - d)
    lo, hi = 0.0, 0.5
    for _ in range(80):
        mid = (lo + hi) / 2
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
    print("   bracket 1/(2K) - e^{-1/K + 2delta}/(2 - delta) at delta -> 0: %.5f (claim 0.077); > 0 for delta < %.5f"
          % (f(0.0), lo))
    # with the sharper exponent e^{-1/K + delta} that the argument actually gives
    f2 = lambda d: 1 / (2 * K) - np.exp(-1 / K + d) / (2 - d)
    lo2, hi2 = 0.0, 0.9
    for _ in range(80):
        mid = (lo2 + hi2) / 2
        if f2(mid) > 0:
            lo2 = mid
        else:
            hi2 = mid
    print("   (with e^{-1/K + delta}: > 0 for delta < %.5f);  largest K with bracket > 0 at delta -> 0: K < %.4f"
          % (lo2, [k for k in np.linspace(1.0, 3, 20001) if 1 / (2 * k) - exp(-1 / k) / 2 > 0][-1]))
    # inequalities: |1 - p^k| <= k K delta, (2j-1)!! (K d)^j <= min(K d, sqrt2 e^{-j}) for j <= j1, tail bound
    worst = [0, 0, 0]
    for k in range(400):
        d = float(np.exp(rng.uniform(log(1e-4), log(0.1))))
        # p with |1 - p| <= K d and |p| = 1 - d: solve for the arg range
        rp = 1 - d
        cphi = (1 + rp ** 2 - (K * d) ** 2) / (2 * rp)
        phimax = np.arccos(min(1.0, cphi))
        phi = rng.uniform(-phimax, phimax) if k % 3 else rng.choice([-1, 1]) * phimax
        p = rp * np.exp(1j * phi)
        kk = np.arange(1, 4000)
        worst[0] = max(worst[0], (np.abs(1 - p ** kk) / (kk * K * d)).max())
        j1 = int(1 / (2 * K * d))
        for j in [1, 2, max(1, j1 // 2), j1]:
            if j < 1:
                continue
            ldf = np.log(np.arange(1, 2 * j, 2, dtype=float) * K * d).sum()
            worst[1] = max(worst[1], exp(ldf - min(log(K * d), 0.5 * log(2) - j)))
        # tail bound for j > j1 (in logarithms)
        lP = 0.0
        mx = -1e300
        lbnd = 0.5 * log(2) - j1 + abs(p) ** (2 * j1 + 1) / (1 - abs(p) ** 2)
        for j in range(1, j1 + 400):
            lP += log(abs(1 - p ** (2 * j - 1)))
            if j > j1:
                mx = max(mx, lP - lbnd)
        worst[2] = max(worst[2], exp(mx))
    print("   max |1-p^k|/(k K delta) = %.4f; max (2j-1)!!(K d)^j / min(K d, sqrt2 e^{-j}) = %.4f;"
          " max |(p;p^2)_j| / tail bound (j > j1) = %.4f   (claims <= 1)" % tuple(worst))
    print("   W(p) - 1 and W3(p) - 1 on the edge |1 - p| = (3/2)(1 - |p|) and on the path of v20 (X = 6):")
    for d in [1e-1, 3e-2, 1e-2, 3e-3, 1e-3]:
        rp = 1 - d
        phimax = np.arccos(min(1.0, (1 + rp ** 2 - (K * d) ** 2) / (2 * rp)))
        p = rp * np.exp(1j * phimax)
        # the path p = e^{-pi i/x}, x = 1/X + i eps with 1 - |p| = d
        X = 6
        eps = -log(1 - d) / (pi * X ** 2)            # first-order inversion of |p|
        xx = 1 / X + 1j * eps
        pp = np.exp(-1j * pi / xx)
        print("    delta=%.0e: edge |W-1|=%.2e |W3-1|=%.2e;  path |1-p|/(1-|p|)=%.4f |W-1|=%.2e |W3-1|=%.2e"
              % (d, abs(W_series(p) - 1), abs(W_series(p, 3) - 1), abs(1 - pp) / (1 - abs(pp)),
                 abs(W_series(pp) - 1), abs(W_series(pp, 3) - 1)))


# ------------------------------------------------------------------ (c)
def path_p(X, eps):
    rho = exp(-pi * eps * X ** 2 / (1 + eps ** 2 * X ** 2))
    phi = pi * eps ** 2 * X ** 3 / (1 + eps ** 2 * X ** 2)
    return rho, phi


def series_vals(p, dps):
    mp.mp.dps = dps
    p = mp.mpc(p)
    # F_a, F_b
    Fa = Fb = mp.mpc(0)
    P = mp.mpc(1)
    n = 1
    minlow, minhigh = mp.inf, mp.inf
    while True:
        P *= 1 - p ** (2 * n - 1)
        ta = p ** (2 * n * n - n) / P
        tb = p ** (2 * n * n - 3 * n + 1) / P
        Fa += ta
        Fb += tb
        if n > 10 and abs(ta) + abs(tb) < mp.mpf(10) ** (-20):
            break
        n += 1
        if n > 100000:
            break
    # psi(p)
    s = mp.mpc(0)
    n = 0
    while True:
        term = p ** (n * (n + 1) // 2)
        s += term
        n += 1
        if n > 10 and abs(term) < mp.mpf(10) ** (-dps + 3):
            break
    # P_+-
    Pp = Pm = mp.mpc(0)
    Q2 = mp.mpc(1)
    k = 0
    while True:
        if k > 0:
            Q2 *= 1 - p ** (2 * k)
        tp = p ** (2 * k * k + k) / Q2
        tm = p ** (2 * k * k - k) / Q2
        Pp += tp
        Pm += tm
        k += 1
        if k > 10 and abs(tp) + abs(tm) < mp.mpf(10) ** (-dps + 3) * (abs(Pp) + 1):
            break
    return complex(Fa), complex(Fb), s, complex(Pp), complex(Pm)


def pochhammer_profile(p, nmax):
    """|(p;p^2)_n| for n = 1..nmax (double)."""
    k = np.arange(1, 2 * nmax, 2)
    f = np.abs(1 - p ** k)
    return np.cumprod(f)


def part_c():
    print("\n(c) Lemma R2' (p -> -1 along x = 1/X + i eps, X odd)")
    for X in [5, 7, 11, 21]:
        line = "   X=%2d:" % X
        for delta in [0.1, 0.03, 0.01, 0.003]:
            eps = -log(1 - delta) / (pi * X ** 2)
            for _ in range(20):        # solve 1 - rho = delta exactly
                rho, phi = path_p(X, eps)
                eps *= log(1 - delta) / log(rho)
            rho, phi = path_p(X, eps)
            line += "  d=%.3f: phi/d^2=%.3f" % (1 - rho, phi / (1 - rho) ** 2)
        print(line + "   (claim phi <= delta^2; -> 1/(pi X) = %.3f)" % (1 / (pi * X)))
    print("   lower bounds for |(p;p^2)_n| (claim: >= (9/8)^n for n <= n1, >= e^{-0.3/delta} for n > n1):")
    for X in [5, 11]:
        for delta in [0.1, 0.03, 0.01, 0.003, 0.001]:
            eps = -log(1 - delta) / (pi * X ** 2)
            for _ in range(20):
                rho, phi = path_p(X, eps)
                eps *= log(1 - delta) / log(rho)
            rho, phi = path_p(X, eps)
            d = 1 - rho
            n1 = len([k for k in range(1, int(1 / d) + 1, 2)])
            out = []
            for tag, ph in [("path", phi), ("phi=4d^2", 4 * d * d), ("phi=d", d), ("phi=3d", 3 * d)]:
                p = -rho * np.exp(1j * ph)
                prof = pochhammer_profile(p, n1 + int(8 / d))
                n = np.arange(1, len(prof) + 1)
                r_lo = (prof[:n1] / (9 / 8) ** n[:n1]).min()
                r_hi = float(np.exp(np.log(prof[n1:]).min() + 0.3 / d))
                out.append("%s: %.3f/%s" % (tag, r_lo, ("%.3g" % r_hi) if np.isfinite(r_hi) else ">1e300"))
            print("    X=%2d d=%.3f n1=%d:  %s" % (X, d, n1, "   ".join(out)))
    print("   G(Q) bound log G <= 2 log(1/(1-Q)) + pi^2/(15 lambda), and |psi(p)| delta^{1/2} e^{pi^2/(8 delta)}, |P+-| / (2 delta^-2 e^{pi^2/(30 delta)}):")
    mp.mp.dps = 30
    gw = -1e9
    for Q in np.concatenate([np.linspace(0.01, 0.99, 99), 1 - np.geomspace(1e-4, 1e-2, 20)]):
        lam = -log(Q)
        nn = np.arange(1, int(60 / lam) + 10)
        nn = nn[(nn % 5 == 1) | (nn % 5 == 4)]
        lg = float(-np.log1p(-Q ** nn).sum())          # log G(Q) = -sum_{n = +-1 mod 5} log(1 - Q^n)
        gw = max(gw, lg - (2 * log(1 / (1 - Q)) + pi ** 2 / (15 * lam)))
    print("    max of log G - bound over Q in (0,1): %.4f (claim <= 0)" % gw)
    for X in [5, 11]:
        for delta in [0.1, 0.05, 0.02, 0.01]:
            eps = -log(1 - delta) / (pi * X ** 2)
            for _ in range(20):
                rho, phi = path_p(X, eps)
                eps *= log(1 - delta) / log(rho)
            rho, phi = path_p(X, eps)
            d = 1 - rho
            p = -rho * np.exp(1j * phi)
            dps = int(25 + (pi ** 2 / (8 * d)) / 2.3)
            Fa, Fb, ps, Pp, Pm = series_vals(p, dps)
            Cpsi = float(abs(ps)) * sqrt(d) * exp(pi ** 2 / (8 * d))
            rP = max(abs(Pp), abs(Pm)) / (2 * d ** -2 * exp(pi ** 2 / (30 * d)))
            W = complex(ps) * Pp - Fa
            pW3 = complex(ps) * Pm - Fb
            print("    X=%2d d=%.3f: C_psi=%.3f  |P+-|/bound=%.2e  F_a+1/3=%.1e F_b-1/3=%.1e  W-1/3=%.1e pW3+1/3=%.1e"
                  % (X, d, Cpsi, rP, abs(Fa + 1 / 3), abs(Fb - 1 / 3), abs(W - 1 / 3), abs(pW3 + 1 / 3)))


# ------------------------------------------------------------------ integrals at x near 1/X
def I_on_L(x, ap, extra, fad, delta=10.0, n=1500):
    """the state integral over L itself (v12's definition), Faddeev's Phi_b; returns (value, cond)."""
    t = -2j * pi * x
    b = np.sqrt(x)
    Rw = np.sqrt(2 * 70 * abs(t) / sin(2 * delta * DEG)) + 8      # Gaussian e^{-|w|^2 sin(2 delta)/(2|t|)}
    R = Rw / (2 * pi * abs(b))
    rr, wt = S6.panels(-R, R, n)
    w = 2 * pi * b * rr * np.exp(-1j * delta * DEG)
    dw = 2 * pi * b * np.exp(-1j * delta * DEG)
    f = np.exp((-w ** 2 / 2 - 2j * pi * ap * w) / t - fad.lam(w) + (w if extra else 0))
    val = (f * wt).sum() * dw / (2 * pi * b)
    return val, np.abs(f).max() * abs(dw) / (2 * pi * abs(b)) / abs(val)


LAST_CHECK = {}


def pieces(X, eps=0.0, version=1):
    """P.J, P.I(1/(2b)), P.I(-1/(2b)) at x = 1/X + i eps (version 0: the b +- 1/(2b) integrals with q^{1/2}).
    For X <= 8 the state integrals are taken over L (their definition; the descent contours of v19 need
    |t| small: for the e^{w} versions the horizontal line is admissible only for |t| < pi/3), with the
    descent-contour values recorded in LAST_CHECK for comparison; for X > 8 over the descent contours
    (over L the cancellation grows like e^{0.16 X})."""
    x = 1 / X + 1j * eps
    t = -2j * pi * x
    tabs, argx = abs(x) * 2 * pi, np.angle(x) / DEG
    b = np.sqrt(x)
    Pf = 2 ** -0.5 * np.exp(pi ** 2 / (12 * t) + t / 6) / b
    fad = core.Faddeev(t)
    PJ = S5.pj_direct(t, version, th=argx - 90)
    qh = np.exp(-t / 2) if version == 0 else 1.0
    try:
        Id, _, _, v1 = S6.I_decay(tabs, argx, 0.5, version == 0, fad)
    except RuntimeError:          # e^{w}-version near |t| = pi/3: tail too long for the panel cap
        if X > 8:
            raise
        Id, v1 = np.nan, 0.0
    Ig, _, _, v0 = S6.I_grow(tabs, argx, version == 0, fad)
    Ipd, Ipg = Id * np.exp(v1 / t), Ig * np.exp(v0 / t)
    if X <= 8:
        IL1, c1 = I_on_L(x, 0.5, version == 0, fad)
        IL2, c2 = I_on_L(x, -0.5, version == 0, fad)
        LAST_CHECK[(X, eps, version)] = (abs(Ipd / IL1 - 1), abs(Ipg / IL2 - 1), c1, c2)
        Ipd, Ipg = IL1, IL2
    return PJ, qh * Pf * Ipd, qh * Pf * Ipg, t


def finite_sums(X):
    mp.mp.dps = 30 + int(0.07 * X)
    z = mp.expjpi(mp.mpf(2) / X)
    out = {}
    for name, nmax in [("half", X // 2), ("full", X)]:
        s1 = s0 = mp.mpc(0)
        P = mp.mpc(1)
        for n in range(nmax):
            if n > 0:
                P *= 1 + z ** n
            s1 += (-1) ** n * z ** (n * (n + 1) // 2) * P
            s0 += (-1) ** n * z ** ((n + 1) * (n + 2) // 2) * P
        out[name] = (complex(s1), complex(s0))
    return out


def part_d():
    print("\n(d) Lemma R3: continuity at x = 1/X: |value(eps) - value(0)| / |value(0)|")
    for X in [4, 5, 8]:
        base = pieces(X, 0.0)
        line = "   X=%d:" % X
        for eps in [1e-2, 1e-3, 1e-4]:
            v = pieces(X, eps)
            line += "  eps=%.0e: J %.1e I+ %.1e I- %.1e" % (eps, *[abs(v[k] - base[k]) / abs(base[k]) for k in range(3)])
        print(line)


def part_e():
    print("\n(e) Theorems 9, 11 at x = 1/X (relative error of the three-integral formula)")
    print("   %3s %12s %12s %12s %12s   %s" % ("X", "psi1 formula", "psi0 formula", "psi1/lead-1", "psi0/lead-1", "note"))
    for X in [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 20, 21, 25, 30, 31, 40, 41]:
        fs = finite_sums(X)
        PJ1, PIp1, PIm1, t = pieces(X, 0.0, 1)
        PJ0, PIp0, PIm0, t = pieces(X, 0.0, 0)
        if X % 2 == 0:
            f1 = PJ1 + PIp1 - PIm1
            f0 = PJ0 + PIp0 - PIm0
            v1, v0 = fs["half"]
            w = 1.0
        else:
            f1 = PJ1 - PIp1 / 3 - PIm1 / 3
            f0 = PJ0 - PIp0 / 3 - PIm0 / 3
            v1, v0 = fs["full"][0] / 3, fs["full"][1] / 3
            w = 1 / 3
        lead = w * 3 ** -0.25 * sqrt(X / 2) * np.exp(CL2 * X / (2 * pi)) * np.exp(2j * pi / 3 - 1j * pi * X / 6)
        note = ""
        if X in (2, 3):
            note = "outside the stated range"
        if X == 4:
            note = "edge |x| = 1/4"
        e0 = abs(f0 - v0) / abs(v0) if abs(v0) > 1e-12 else abs(f0 - v0)
        extra = ""
        if (X, 0.0, 1) in LAST_CHECK:
            c = LAST_CHECK[(X, 0.0, 1)]
            d = LAST_CHECK[(X, 0.0, 0)]
            extra = "  [L vs descent contours: psi1 %.0e/%.0e, psi0 %.0e/%.0e]" % (c[0], c[1], d[0], d[1])
        print("   %3d %12.2e %12.2e %12.3e %12.3e   %-24s X*(psi1/lead - 1) = %.3f%s"
              % (X, abs(f1 - v1) / abs(v1), e0, abs(v1 / lead - 1),
                 abs(v0 / (lead * np.exp(2j * pi / 3)) - 1) if abs(v0) > 1e-12 else float('nan'), note,
                 X * abs(v1 / lead - 1), extra))


if __name__ == "__main__":
    T = time.time()
    part_a()
    part_b()
    part_c()
    part_d()
    part_e()
    print("\nrun time %.0f s" % (time.time() - T))
