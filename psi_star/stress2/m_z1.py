"""M: z = 1 and z = q (v24).

 (a) Proposition M: I(0), I(1/b), I(b), I(b + 1/b) by quadrature over a line equivalent to L, at random x
     with |q| in [0.05, 0.95] (Re x in [-1.5, 1.5], arg b <= 70 deg), against the closed forms;
 (b) Theorem 16: the four laws, with M1, M0 from J(1/b), J(b + 1/b) by quadrature and own q-series;
 (c) Theorem 17: M1, M0 against (60t/pi)^{1/2} L(2/5, 10t), L(1/5, 10t), with an own Mordell integral;
 (d) (R3), (R4) exactly to q^1000 (python-flint, modulo two primes) -- also (R1), (R2);
 (e) Theorem 17 on the unit circle x = 1/X, X = 1 ... 100 (Faddeev; L on the ray at 45 deg), and the
     continuity of J(1/b), J(b + 1/b) as x = 1/X + i eps -> 1/X.
Run:  python3 m_z1.py   (about 10 min)
"""
import time
import numpy as np
import mpmath as mp
from math import pi, sqrt, sin, log
import flint
import common2 as c
import core

DEG = pi / 180
rng = np.random.default_rng(24)
s1, s2 = sin(pi / 5), sin(2 * pi / 5)
EPS = np.exp(-1j * pi / 4)


def qpow(x, al):
    return np.exp(2j * pi * x * al)


def ppow(x, al):
    return np.exp(-1j * pi * al / x)


def closed_M(x):
    k = EPS * 2 / sqrt(5)
    return {"I(0)": k * s2 * qpow(x, -1 / 60) * ppow(x, 1 / 30),
            "I(1/b)": k * s1 * qpow(x, -1 / 60) * ppow(x, 19 / 30),
            "I(b)": k * s1 * qpow(x, -19 / 60) * ppow(x, 1 / 30),
            "I(b+1/b)": -k * s2 * qpow(x, -19 / 60) * ppow(x, 19 / 30)}


def random_x(nearreal=False):
    while True:
        absq = rng.uniform(0.9, 0.95) if nearreal else rng.uniform(0.05, 0.95)
        imx = -log(absq) / (2 * pi)
        x = complex(rng.uniform(-1.5, 1.5), imx)
        if np.angle(np.sqrt(x)) / DEG <= 70:
            return x


def integrals(x):
    """the six integrals by quadrature; u given as (ap, extra) with u = ap/b + extra*b."""
    t = -2j * pi * x
    b = np.sqrt(x)
    beta = np.angle(b) / DEG
    delta = (90 - beta) / 2
    out, cond = {}, {}
    for name, ap, ex, ker in [("I(0)", 0.0, 0.0, False), ("I(1/b)", 1.0, 0.0, False), ("I(b)", 0.0, 1.0, False),
                              ("I(b+1/b)", 1.0, 1.0, False), ("J(1/b)", 1.0, 0.0, True), ("J(b+1/b)", 1.0, 1.0, True)]:
        v, cd = c.line_integral(x, ap, kernel=ker, delta=delta, fad=None, extra=ex)
        out[name], cond[name] = v, cd
    return out, cond


def series(q, kind, N=None):
    """Ramanujan's fifth-order series at a complex q (double, direct summation)."""
    s = 0j
    P = 1 + 0j
    n = 0
    while True:
        if kind == "f0":
            if n > 0:
                P /= 1 + q ** n
            term = q ** (n * n) * P
        elif kind == "f1":
            if n > 0:
                P /= 1 + q ** n
            term = q ** (n * n + n) * P
        elif kind == "psi0":
            if n > 0:
                P *= 1 + q ** n
            term = q ** ((n + 1) * (n + 2) // 2) * P
        elif kind == "psi1":
            if n > 0:
                P *= 1 + q ** n
            term = q ** (n * (n + 1) // 2) * P
        elif kind == "F0":
            if n > 0:
                P /= 1 - q ** (2 * n - 1)
            term = q ** (2 * n * n) * P
        elif kind == "F1":
            P /= 1 - q ** (2 * n + 1)
            term = q ** (2 * n * n + 2 * n) * P
        elif kind == "phi0":
            if n > 0:
                P *= 1 + q ** (2 * n - 1)
            term = q ** (n * n) * P
        elif kind == "phi1":
            if n > 0:
                P *= 1 + q ** (2 * n - 1)
            term = q ** ((n + 1) ** 2) * P
        s += term
        n += 1
        if n > 5 and abs(term) < 1e-18 * (abs(s) + 1e-300):
            break
        if n > 100000:
            break
    return s


def mordell(r, alpha):
    """L(r, alpha) = int_0^inf e^{-3 alpha u^2/2} [cosh((3r-2) alpha u) + cosh((3r-1) alpha u)]/cosh(3 alpha u/2) du,
    on the ray arg u = -arg(alpha)/2 (no pole crossed, see the docstring of the module)."""
    with mp.workdps(30):
        al = mp.mpc(alpha)
        e = mp.expj(-mp.arg(al) / 2)
        f = lambda s: (mp.exp(-3 * al * (s * e) ** 2 / 2) * (mp.cosh((3 * r - 2) * al * s * e) + mp.cosh((3 * r - 1) * al * s * e))
                       / mp.cosh(3 * al * s * e / 2)) * e
        return complex(mp.quad(f, [0, 1, 4, 12, mp.inf]))


def part_abc():
    print("(a)-(c) Proposition M, Theorems 16, 17 at random x (|q| in [0.05, 0.95]; 10 samples with |q| in [0.9, 0.95])")
    wM = {k: 0 for k in ["I(0)", "I(1/b)", "I(b)", "I(b+1/b)"]}
    wL = [0, 0, 0, 0]
    wT = [0, 0]
    wLM = [0, 0, 0, 0]
    worst_cond = 0
    n = skip = 0
    worst_x = None
    while n < 50 and skip < 400:
        x = random_x(nearreal=(n % 5 == 4))
        t = -2j * pi * x
        b = np.sqrt(x)
        try:
            ints, cond = integrals(x)
        except (RuntimeError, ValueError, TypeError):
            skip += 1
            continue
        if any(v > 1e8 for v in cond.values()) or any(not np.isfinite(v) for v in ints.values()):
            skip += 1
            continue
        cm = closed_M(x)
        for k in wM:
            e = abs(ints[k] / cm[k] - 1)
            if e > wM[k]:
                wM[k] = e
                worst_x = x
        q = np.exp(-t)
        p = np.exp(-2 * pi ** 2 / t)
        P = c.Pfac(t, b)
        M1 = 2j * qpow(x, 11 / 60) * P * ints["J(1/b)"]
        M0 = 2j * qpow(x, 29 / 60) * P * ints["J(b+1/b)"]
        cp = np.sqrt(pi * (5 + sqrt(5)) / (10 * t))
        cmn = np.sqrt(pi * (5 - sqrt(5)) / (10 * t))
        F0, F1, f0, f1 = series(p, "F0"), series(p, "F1"), series(q, "f0"), series(q, "f1")
        ph0, ph1 = series(-p, "phi0"), series(-p, "phi1")
        ps0, ps1 = series(q, "psi0"), series(q, "psi1")
        laws = [
            (qpow(x, 11 / 60) * f1, M1 + 2 * cp * ppow(x, -1 / 120) * (F0 - 1) - 2 * cmn * ppow(x, 71 / 120) * F1),
            (2 * qpow(x, 11 / 60) * ps1, -M1 + 2 * cp * ppow(x, -1 / 120) * ph0 + 2 * cmn * ppow(x, -49 / 120) * ph1),
            (qpow(x, -1 / 60) * f0, M0 + 2 * cmn * ppow(x, -1 / 120) * (F0 - 1) + 2 * cp * ppow(x, 71 / 120) * F1),
            (2 * qpow(x, -1 / 60) * ps0, -M0 + 2 * cmn * ppow(x, -1 / 120) * ph0 - 2 * cp * ppow(x, -49 / 120) * ph1)]
        for k, (l, r) in enumerate(laws):
            wL[k] = max(wL[k], abs(l - r) / max(abs(l), 1e-300))
        pre = np.sqrt(60 * t / pi)
        L1, L0 = pre * mordell(0.4, 10 * t), pre * mordell(0.2, 10 * t)
        wT[0] = max(wT[0], abs(M1 / L1 - 1))
        wT[1] = max(wT[1], abs(M0 / L0 - 1))
        # the laws with M1, M0 replaced by the Mordell integrals: q-series and L only, no state integral
        for k, (l, r) in enumerate(laws):
            r2 = r + (L1 - M1) * (1 if k == 0 else -1 if k == 1 else 0) + (L0 - M0) * (1 if k == 2 else -1 if k == 3 else 0)
            wLM[k] = max(wLM[k], abs(l - r2) / max(abs(l), 1e-300))
        worst_cond = max(worst_cond, max(cond.values()))
        n += 1
    print("   %d samples (%d skipped: quadrature conditioning > 1e8 or products too long)" % (n, skip))
    print("   Proposition M, max relative error: " + ", ".join("%s %.1e" % kv for kv in wM.items()) + "  (worst at x = %s)" % worst_x)
    print("   Theorem 16, the four laws (f1, psi1, f0, psi0): " + ", ".join("%.1e" % v for v in wL))
    print("   Theorem 17: M1 vs Mordell %.1e, M0 vs Mordell %.1e  (largest conditioning of a quadrature used: %.0e)"
          % (wT[0], wT[1], worst_cond))
    print("   the four laws with M1, M0 replaced by the Mordell integrals (q-series and L only): " + ", ".join("%.1e" % v for v in wLM))


def part_d(N=1000):
    print("\n(d) series identities to q^%d, exactly modulo two primes (python-flint nmod_poly)" % N)
    for pr in [2 ** 61 - 1, 1000000007]:
        R = lambda coeffs: flint.nmod_poly(coeffs, pr)
        one = R([1])

        def tr(f):
            return flint.nmod_poly(f.coeffs()[:N + 1], pr) if f.degree() > N else f

        def mono(e):
            cf = [0] * (e + 1)
            cf[e] = 1
            return R(cf)

        def inv(f):
            # power series inverse to order N
            g = R([pow(int(f.coeffs()[0]), -1, pr)])
            k = 1
            while k <= N:
                k *= 2
                g = tr((g * (R([2]) - tr(f * g))))
                g = flint.nmod_poly(g.coeffs()[:k], pr)
            return tr(g)

        def onepm(k, sgn):     # 1 + sgn q^k
            cf = [0] * (k + 1)
            cf[0] = 1
            cf[k] = sgn % pr
            return R(cf)

        def sub2(f):           # f(q^2)
            cf = f.coeffs()
            out = [0] * (2 * len(cf))
            for i, v in enumerate(cf):
                out[2 * i] = int(v)
            return tr(R(out))
        # f0, f1, psi0, psi1
        f0 = f1 = ps0 = ps1 = R([0])
        P = one
        Pinv = one
        n = 0
        while n * n <= N:
            if n > 0:
                Pinv = tr(Pinv * inv(onepm(n, 1)))
            f0 = f0 + tr(mono(n * n) * Pinv)
            if n * n + n <= N:
                f1 = f1 + tr(mono(n * n + n) * Pinv)
            n += 1
        n = 0
        while n * (n + 1) // 2 <= N:
            if n > 0:
                P = tr(P * onepm(n, 1))
            ps1 = ps1 + tr(mono(n * (n + 1) // 2) * P)
            if (n + 1) * (n + 2) // 2 <= N:
                ps0 = ps0 + tr(mono((n + 1) * (n + 2) // 2) * P)
            n += 1
        # chi0 = sum q^n/(q^{n+1};q)_n, chi1 = sum q^n/(q^{n+1};q)_{n+1}
        ch0 = ch1 = R([0])
        for n in range(0, N // 2 + 2):
            if 2 * n > N + 2:
                break
            d0 = one
            for j in range(n + 1, 2 * n + 1):
                d0 = tr(d0 * onepm(j, -1))
            d1 = tr(d0 * onepm(2 * n + 1, -1))
            ch0 = ch0 + tr(mono(n) * inv(d0))
            ch1 = ch1 + tr(mono(n) * inv(d1))
        ch0, ch1 = tr(ch0), tr(ch1)
        R3 = tr(ps1 - f1 - tr(mono(1) * sub2(ch1)))
        R4 = tr(ps0 - f0 + R([2]) - sub2(ch0))
        z3 = all(int(v) == 0 for v in R3.coeffs()[:N + 1])
        z4 = all(int(v) == 0 for v in R4.coeffs()[:N + 1])
        print("   mod %d: (R3) psi1 - f1 - q chi1(q^2) = 0: %s;  (R4) psi0 - f0 + 2 - chi0(q^2) = 0: %s" % (pr, z3, z4))


def J_real(x, ap, ex, fad):
    v, cd = c.line_integral(x, ap, kernel=True, delta=45.0, fad=fad, extra=ex)
    return v, cd


def part_e():
    print("\n(e) Theorem 17 on the unit circle: |J/closed form - 1| for J(1/b), J(b+1/b) at x = 1/X")
    rows = []
    for X in list(range(1, 21)) + [25, 30, 40, 50, 60, 70, 80, 90, 100]:
        x = 1.0 / X
        t = -2j * pi * x
        fad = core.Faddeev(t)
        J1, c1 = J_real(x, 1.0, 0.0, fad)
        J0, c0 = J_real(x, 1.0, 1.0, fad)
        pre = EPS * sqrt(15) * t / pi * ppow(x, 1 / 24)
        f1 = pre * qpow(x, -1 / 60) * mordell(0.4, 10 * t)
        f0 = pre * qpow(x, -19 / 60) * mordell(0.2, 10 * t)
        rows.append((X, abs(J1 / f1 - 1), abs(J0 / f0 - 1), c1, c0))
    for i in range(0, len(rows), 5):
        print("   " + "  ".join("X=%3d: %.0e/%.0e" % r[:3] for r in rows[i:i + 5]))
    print("   worst %.1e (conditioning of the line integrals up to %.0e)" % (max(max(r[1], r[2]) for r in rows),
                                                                         max(max(r[3], r[4]) for r in rows)))
    print("   continuity: |J(x0 + i eps) - J(x0)| / |J(x0)|, x0 = 1/X")
    for X in [3, 5, 8, 20]:
        x0 = 1.0 / X
        fad0 = core.Faddeev(-2j * pi * x0)
        base1, _ = J_real(x0, 1.0, 0.0, fad0)
        base0, _ = J_real(x0, 1.0, 1.0, fad0)
        line = "   X=%2d:" % X
        for eps in [1e-2, 1e-3, 1e-4]:
            x = x0 + 1j * eps
            fad = core.Faddeev(-2j * pi * x)
            v1, _ = J_real(x, 1.0, 0.0, fad)
            v0, _ = J_real(x, 1.0, 1.0, fad)
            line += "  eps=%.0e: %.1e %.1e" % (eps, abs(v1 / base1 - 1), abs(v0 / base0 - 1))
        print(line)


if __name__ == "__main__":
    T = time.time()
    part_abc()
    part_d()
    part_e()
    print("\nrun time %.0f s" % (time.time() - T))
