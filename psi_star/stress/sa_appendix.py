"""SA: Appendix [v2 §2, v3 §1, v4 §1].

 (a) identities (I), (II), (V), (VI) and the q-difference equation of G(z) as exact
     integer power series to order q^400;
 (b) odd-root rule, step (c): |(-q;q)_{mJ}| >= 1.4^J when r^{J m^2} >= 1/2, and the
     domination of the terms of f_1^*; step (d): theta_4 R_1 -> 0 and its actual size;
 (c) the odd-root rule itself: psi_1^*(r zeta) -> A/3, psi_0^*(r zeta) -> A_0/3 (m = 3 ... 15,
     all primitive roots of order m <= 9).
Run:  python3 sa_appendix.py     (about 2 min)
"""
import time
import numpy as np
import mpmath as mp
from math import pi, log, sqrt, gcd

NQ = 400


# ------------------------------------------------------------------ exact series in q (Python ints)
def zero():
    return [0] * (NQ + 1)


def mono(e, c=1):
    s = zero()
    if 0 <= e <= NQ:
        s[e] = c
    return s


def add(a, b, sb=1):
    return [x + sb * y for x, y in zip(a, b)]


def mul(a, b):
    out = zero()
    nz = [(i, x) for i, x in enumerate(a) if x]
    for j, y in enumerate(b):
        if y:
            for i, x in nz:
                if i + j > NQ:
                    break
                out[i + j] += x * y
    return out


def times_1pqk(a, k, sgn=1):
    """a * (1 + sgn q^k)"""
    out = a[:]
    for i in range(NQ, k - 1, -1):
        out[i] += sgn * a[i - k]
    return out


def div_1pqk(a, k, sgn=1):
    """a / (1 + sgn q^k)"""
    out = a[:]
    for i in range(k, NQ + 1):
        out[i] -= sgn * out[i - k]
    return out


def shift(a, e, c=1):
    out = zero()
    for i, x in enumerate(a):
        if x and 0 <= i + e <= NQ:
            out[i + e] += c * x
    return out


def psi_stars():
    p1 = zero()
    p0 = zero()
    P = mono(0)            # (-q;q)_n
    n = 0
    while n * (n + 1) // 2 <= NQ:
        if n > 0:
            P = times_1pqk(P, n, +1)
        p1 = add(p1, shift(P, n * (n + 1) // 2, (-1) ** n))
        p0 = add(p0, shift(P, (n + 1) * (n + 2) // 2, (-1) ** n))
        n += 1
    return p1, p0


def f_stars():
    f1, f0 = zero(), zero()
    P = mono(0)            # 1/(-q;q)_n
    n = 0
    while n * n <= NQ:
        if n > 0:
            P = div_1pqk(P, n, +1)
        f1 = add(f1, shift(P, n * n + n, (-1) ** n))
        f0 = add(f0, shift(P, n * n, (-1) ** n))
        n += 1
    return f1, f0


def R_series():
    r1, r0 = zero(), zero()
    P = mono(0)            # 1/(q;q)_n
    n = 0
    while n * n <= NQ:
        if n > 0:
            P = div_1pqk(P, n, -1)
        r1 = add(r1, shift(P, n * n + n, (-1) ** n))
        r0 = add(r0, shift(P, n * n, (-1) ** n))
        n += 1
    return r1, r0


def theta4():
    s = zero()
    n = 0
    while n * n <= NQ:
        s[n * n] += (1 if n == 0 else 2) * (-1) ** n
        n += 1
    return s


def W_series(e3=False):
    s = zero()
    P = mono(0)            # (q;q^2)_j
    j = 0
    while j * j <= NQ:
        if j > 0:
            P = times_1pqk(P, 2 * j - 1, -1)
        s = add(s, shift(P, j * j + (3 * j if e3 else j), (-1) ** j))
        j += 1
    return s


def F_ab(which):
    s = zero()
    P = mono(0)            # 1/(q;q^2)_n
    for n in range(1, 40):
        P = div_1pqk(P, 2 * n - 1, -1)
        e = 2 * n * n - n if which == "a" else 2 * n * n - 3 * n + 1
        if e > NQ:
            break
        s = add(s, shift(P, e))
    return s


def P_pm(sign):
    s = zero()
    P = mono(0)            # 1/(q^2;q^2)_k
    k = 0
    while 2 * k * k - k <= NQ:
        if k > 0:
            P = div_1pqk(P, 2 * k, -1)
        s = add(s, shift(P, 2 * k * k + sign * k))
        k += 1
    return s


def psi_ram():
    s = zero()
    n = 0
    while n * (n + 1) // 2 <= NQ:
        s[n * (n + 1) // 2] += 1
        n += 1
    return s


def part_a():
    print("(a) identities as exact power series to q^%d" % NQ)
    p1, p0 = psi_stars()
    f1, f0 = f_stars()
    r1, r0 = R_series()
    th = theta4()
    I_ = add(add(shift(p0, 0, 2), f0, -1), mul(th, r0))
    II = add(add(shift(p1, 0, 2), f1, -1), mul(th, r1), -1)
    W, W3 = W_series(), W_series(True)
    ps = psi_ram()
    V_ = add(add(W, F_ab("a")), mul(ps, P_pm(+1)), -1)
    VI = add(add(shift(W3, 1), F_ab("b")), mul(ps, P_pm(-1)), -1)
    nz = lambda s: max([i for i, x in enumerate(s) if x] or [-1])
    print("   (I)   2 psi0* - f0* + theta4 R0:  %s" % ("0 to order %d" % NQ if nz(I_) < 0 else "nonzero at q^%d" % min(i for i, x in enumerate(I_) if x)))
    print("   (II)  2 psi1* - f1* - theta4 R1:  %s" % ("0 to order %d" % NQ if nz(II) < 0 else "nonzero at q^%d" % min(i for i, x in enumerate(II) if x)))
    print("   (V)   W + F_a - psi P+:           %s" % ("0 to order %d" % NQ if nz(V_) < 0 else "nonzero at q^%d" % min(i for i, x in enumerate(V_) if x)))
    print("   (VI)  q W3 + F_b - psi P-:        %s" % ("0 to order %d" % NQ if nz(VI) < 0 else "nonzero at q^%d" % min(i for i, x in enumerate(VI) if x)))
    # q-difference equation G(z) = 1 + z q G(z q^2) - z q^2 G(z q^4) at z = -q^{2k+1}: G(-q^5) etc.
    # check the two quoted examples: G(-q^5) = (W - 1 + q^2 W3)/q^3 and G(-1/q) = 1 - W + q W3
    def G_at(e):  # G(-q^e) = sum (-1)^j q^{j^2 + e j} (q;q^2)_j, valid for e >= -1 as a power series
        s = zero()
        P = mono(0)
        j = 0
        while j * j + e * j <= NQ + 50 and j < 60:
            if j > 0:
                P = times_1pqk(P, 2 * j - 1, -1)
            s = add(s, shift(P, j * j + e * j, (-1) ** j))
            j += 1
        return s
    g5 = G_at(5)
    rhs5 = add(add(W, mono(0), -1), shift(W3, 2))     # = q^3 G(-q^5)
    d5 = add(shift(g5, 3), rhs5, -1)
    gm1 = G_at(-1)
    dm1 = add(gm1, add(add(mono(0), W, -1), shift(W3, 1)), -1)
    print("   G(-q^5) = (W - 1 + q^2 W3)/q^3: %s;  G(-1/q) = 1 - W + q W3: %s"
          % ("0 to order %d" % (NQ - 3) if nz(d5[:NQ - 3]) < 0 else "fails", "0 to order %d" % NQ if nz(dm1) < 0 else "fails"))
    # double-sum forms of v2
    lhs1 = zero()
    lhs0 = zero()
    Pk = mono(0)
    for k in range(0, 40):
        if k > 0:
            Pk = div_1pqk(Pk, k, -1)
        for n in range(0, 30):
            e1 = k * (k + 1) // 2 + n * n + n + k * n
            e0 = k * (k + 1) // 2 + n * n + k * n
            if k % 2 == 0 and e1 <= NQ:
                lhs1 = add(lhs1, shift(Pk, e1, (-1) ** n))
            if k % 2 == 1 and e0 <= NQ:
                lhs0 = add(lhs0, shift(Pk, e0, (-1) ** n))
    mq = mono(0)
    for k in range(1, NQ + 1):
        mq = times_1pqk(mq, k, +1)          # (-q;q)_inf
    p1, p0 = psi_stars()
    dd1 = add(mul(mq, p1), lhs1, -1)
    dd0 = add(mul(mq, p0), lhs0, -1)
    print("   double-sum forms (-q)_inf psi1*, (-q)_inf psi0*: %s, %s"
          % ("0 to order %d" % NQ if nz(dd1) < 0 else "nonzero at q^%d" % min(i for i, x in enumerate(dd1) if x),
             "0 to order %d" % NQ if nz(dd0) < 0 else "nonzero at q^%d" % min(i for i, x in enumerate(dd0) if x)))


# ------------------------------------------------------------------ (b), (c) numerics at roots of unity
def part_b():
    print("\n(b) odd-root rule, steps (c) and (d)")
    mp.mp.dps = 30
    worst = 1e9
    for m in [3, 5, 7, 9, 11]:
        zeta = mp.expjpi(mp.mpf(2) / m)
        for d in [1e-2, 1e-3, 1e-4]:
            r = 1 - mp.mpf(d)
            q = r * zeta
            P = mp.mpc(1)
            J = 0
            n = 0
            while True:
                n += 1
                P *= 1 + q ** n
                if n % m == 0:
                    J = n // m
                    if r ** (J * m * m) >= 0.5:
                        worst = min(worst, float(abs(P) / mp.mpf(1.4) ** J))
                    else:
                        break
    print("   min |(-q;q)_{mJ}| / 1.4^J over J with r^{J m^2} >= 1/2 (m <= 11, 1-r >= 1e-4): %.4f (claim >= 1)" % worst)
    print("   theta4 R1 at q = zeta e^{-t}: m^2 t log|theta4 R1| (claimed bound exponent -(pi^2/4 - 1 - eps); "
          "claimed true value -pi^2/3 = %.3f):" % (-pi ** 2 / 3))
    for m in [3, 5, 7]:
        zeta = mp.expjpi(mp.mpf(2) / m)
        line = "    m=%d:" % m
        for m2t in [0.5, 0.25, 0.125, 1 / 16, 1 / 32]:
            t = mp.mpf(m2t) / m ** 2
            mp.mp.dps = int(20 + 4 / (m2t * 2.3))
            q = zeta * mp.exp(-t)
            th4 = mp.jtheta(4, 0, q)
            R1 = mp.mpc(0)
            P = mp.mpc(1)
            n = 0
            while True:
                if n > 0:
                    P *= 1 - q ** n
                term = (-1) ** n * q ** (n * n + n) / P
                R1 += term
                n += 1
                if n > 10 and abs(term) < mp.mpf(10) ** (-mp.mp.dps + 5):
                    break
            line += "  m^2t=%.4f: %.3f" % (m2t, float(m2t * mp.log(abs(th4 * R1))))
        print(line)


def part_c():
    print("\n(c) odd-root rule: |psi1*(r zeta) - A/3|, |psi0*(r zeta) - A0/3| as r -> 1, all primitive roots, m <= 9")
    mp.mp.dps = 60      # the terms reach e^{(ln 2/m)^2/(2(1-r))} before they decay
    for m in [3, 5, 7, 9]:
        line = "   m=%d:" % m
        worst = {1e-2: 0, 1e-3: 0, 1e-4: 0}
        for h in range(1, m):
            if gcd(h, m) != 1:
                continue
            zeta = mp.expjpi(mp.mpf(2 * h) / m)
            A = A0 = mp.mpc(0)
            P = mp.mpc(1)
            for n in range(m):
                if n > 0:
                    P *= 1 + zeta ** n
                A += (-1) ** n * zeta ** (n * (n + 1) // 2) * P
                A0 += (-1) ** n * zeta ** ((n + 1) * (n + 2) // 2) * P
            for d in worst:
                mp.mp.dps = int(30 + (log(2) / m) ** 2 / (2 * d) / 2.3)
                q = (1 - mp.mpf(d)) * zeta
                s1 = s0 = mp.mpc(0)
                P = mp.mpc(1)
                n = 0
                while True:
                    if n > 0:
                        P *= 1 + q ** n
                    t1 = (-1) ** n * q ** (n * (n + 1) // 2) * P
                    s1 += t1
                    s0 += t1 * q ** (n + 1)
                    n += 1
                    if n > 30 and abs(t1) < mp.mpf(10) ** -22 and n * n * d > 120:
                        break
                worst[d] = max(worst[d], float(abs(s1 - A / 3) / max(abs(A), 1)), float(abs(s0 - A0 / 3) / max(abs(A0), 1)))
        for d in worst:
            line += "  1-r=%.0e: %.2e" % (d, worst[d])
        print(line + "   (relative to max(|A|, 1); A0 = 0 for m = 3; should scale like 1 - r)")


if __name__ == "__main__":
    T = time.time()
    part_a()
    part_b()
    part_c()
    print("\nrun time %.0f s" % (time.time() - T))
