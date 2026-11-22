"""Quantum modularity test.
phi(x) = e^{2 pi i alpha x} * F(e^{2 pi i x}) at rationals x (exact terminating sums).
For gamma = [[a,b],[c,d]] the cocycle h(x) = phi(x) - eps*(c x + d)^{-k} * phi(gamma x) should be smooth
(real-analytic) on an interval if phi is quantum modular of weight k.  eps is fitted in closed form;
roughness = sum |Delta h|^2 / sum |Delta phi|^2 over neighbouring rationals (small = smooth)."""
import mpmath as mp
from math import gcd
mp.mp.dps = 40
TWOPI_I = 2j * mp.pi

def kontsevich(h, m):           # F(zeta) = sum_{n<m} (zeta;zeta)_n
    z = mp.e ** (TWOPI_I * mp.mpf(h) / m); s = 0; P = 1
    for n in range(m):
        s += P; P *= (1 - z ** (n + 1))
    return s
def psi1s(h, m):                # sum (-1)^n zeta^{n(n+1)/2} (-zeta;zeta)_n, terminates for even m
    z = mp.e ** (TWOPI_I * mp.mpf(h) / m); s = 0; P = 1
    for n in range(m):
        s += (-1) ** n * z ** (n * (n + 1) // 2) * P
        P *= (1 + z ** (n + 1))
        if abs(P) < mp.mpf(10) ** -30: break
    return s
def psi0s(h, m):                # sum (-1)^n zeta^{(n+1)(n+2)/2} (-zeta;zeta)_n
    z = mp.e ** (TWOPI_I * mp.mpf(h) / m); s = 0; P = 1
    for n in range(m):
        s += (-1) ** n * z ** ((n + 1) * (n + 2) // 2) * P
        P *= (1 + z ** (n + 1))
        if abs(P) < mp.mpf(10) ** -30: break
    return s

def reduce(p, q):
    if q < 0: p, q = -p, -q
    g = gcd(p, q); return p // g, q // g

def roughness(F, alpha, k, gam, pts, cyc):
    a, b, c, d = gam
    D1, D2 = [], []
    vals = []
    for (h, m) in pts:
        x = mp.mpf(h) / m
        p2, q2 = reduce(a * h + b * m, c * h + d * m)
        phi1 = mp.e ** (TWOPI_I * alpha * x) * F(h, m)
        y = mp.mpf(p2) / q2
        phi2 = mp.e ** (TWOPI_I * alpha * y) * F(p2 % q2, q2)
        fac = cyc(x) ** (-k)
        vals.append((x, phi1, fac * phi2))
    vals.sort(key=lambda v: v[0])
    for (x1, f1, g1), (x2, f2, g2) in zip(vals, vals[1:]):
        D1.append(f2 - f1); D2.append(g2 - g1)
    num = sum(u * mp.conj(v) for u, v in zip(D1, D2)); den = sum(abs(v) ** 2 for v in D2)
    eps = num / den
    R = sum(abs(u - eps * v) ** 2 for u, v in zip(D1, D2)) / sum(abs(u) ** 2 for u in D1)
    return R, eps

def points(lo, hi, maxden, even_only):
    out = []
    for m in range(1, maxden + 1):
        if even_only and m % 2: continue
        for h in range(int(lo * m) + 1, int(hi * m) + 1):
            if gcd(h, m) == 1 and lo < h / m < hi: out.append((h, m))
    return out
