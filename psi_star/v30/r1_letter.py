#!/usr/bin/env python3
"""Addendum v30 section 1: (R3), (R4) of v24 from Ramanujan's relations between the fifth-order functions (last letter; proved by Watson).
Exact integer power series to order N (default 300).

 first family    (A0)  f0(q) + 2 F0(q^2) - 2 = phi0(-q^2) + psi0(q)        [ = theta_4(0,-q) G(-q) ]
                 (C0)  chi0(q) = 2 F0(q) - phi0(-q)                                              (R1 of v24)
 second family   (A1)  f1(q) + 2q F1(q^2) = psi1(q) - phi1(-q^2)/q         [ = theta_4(0,-q) H(-q) ]
                 (C1)  q chi1(q) = 2q F1(q) + phi1(-q)                                           (R2 of v24)
 consequences    (R4)  psi0(q) - f0(q) + 2 = chi0(q^2)        = (A0) and (C0) at q^2
                 (R3)  psi1(q) - f1(q)     = q chi1(q^2)      = (A1) and (C1) at q^2
 also checked    2 phi0(-q^2) - f0(q) = theta_4(0,q) G(q),   f1(q) + 2 phi1(-q^2)/q = theta_4(0,q) H(q),
                 psi0(q) - F0(q^2) + 1 = q psi(q^2) H(q^4),  psi1(q) - q F1(q^2) = psi(q^2) G(q^4)    (Lost Notebook, as quoted in arXiv:2007.03664)
"""
import sys
N = int(sys.argv[1]) if len(sys.argv) > 1 else 300

def mul(a, b):
    c = [0]*(N + 1)
    for i, x in enumerate(a):
        if x:
            for j in range(N + 1 - i):
                if b[j]: c[i + j] += x*b[j]
    return c
def inv(a):
    assert a[0] in (1, -1)
    b = [0]*(N + 1); b[0] = a[0]
    for n in range(1, N + 1): b[n] = -a[0]*sum(a[k]*b[n - k] for k in range(1, n + 1))
    return b
def one(): return [1] + [0]*N
def zero(): return [0]*(N + 1)
def mf(a, k, c=1):
    """a * (1 - c q^k)"""
    if k > N: return a[:]
    return [a[i] - (c*a[i - k] if i >= k else 0) for i in range(N + 1)]
def df(a, k):
    """a / (1 - q^k)"""
    b = a[:]
    for i in range(k, N + 1): b[i] += b[i - k]
    return b
def shift(a, k): return ([0]*k + a[:N + 1 - k]) if k <= N else zero()
def add(a, b, cb=1): return [x + cb*y for x, y in zip(a, b)]
def sc(a, c): return [c*x for x in a]
def sub(a, m, sign=1):
    """a(sign q^m)"""
    b = zero()
    for i, x in enumerate(a):
        if i*m > N: break
        b[i*m] = x*(sign**i)
    return b
def unshift(a):
    assert a[0] == 0
    return a[1:] + [0]

f0 = zero(); f1 = zero(); p0 = zero(); p1 = zero(); ph0 = zero(); ph1 = zero(); F0 = zero(); F1 = zero(); c0 = zero(); c1 = zero()
P = one(); Pi = one(); n = 0                       # (-q;q)_n and its inverse
while n*(n + 1)//2 <= N:
    f0 = add(f0, shift(Pi, n*n)); f1 = add(f1, shift(Pi, n*n + n))
    p0 = add(p0, shift(P, (n + 1)*(n + 2)//2)); p1 = add(p1, shift(P, n*(n + 1)//2))
    n += 1; P = mf(P, n, -1); Pi = inv(P)
P = one(); n = 0                                   # (-q;q^2)_n
while n*n <= N:
    ph0 = add(ph0, shift(P, n*n)); ph1 = add(ph1, shift(P, (n + 1)*(n + 1)))
    P = mf(P, 2*n + 1, -1); n += 1
Pi = one(); n = 0                                  # 1/(q;q^2)_n
while 2*n*n <= N:
    F0 = add(F0, shift(Pi, 2*n*n))
    Pi = df(Pi, 2*n + 1); F1 = add(F1, shift(Pi, 2*n*(n + 1))); n += 1
Di = one()                                         # 1/(q^{n+1};q)_n
for n in range(N + 1):
    if n > 0:
        Di = df(df(Di, 2*n - 1), 2*n); Di = mf(Di, n)
    c0 = add(c0, shift(Di, n)); c1 = add(c1, shift(df(Di, 2*n + 1), n))
G = one(); H = one(); mq = one()
for k in range(1, N + 1):
    if k % 5 in (1, 4): G = df(G, k)
    if k % 5 in (2, 3): H = df(H, k)
th = zero(); ps = zero(); n = 0
while n*n <= N:
    th[n*n] += (1 if n == 0 else 2)*(-1)**n; n += 1
n = 0
while n*(n + 1)//2 <= N:
    ps[n*(n + 1)//2] += 1; n += 1
th3 = sub(th, 1, -1)
ph0m = sub(ph0, 2, -1); ph1m = unshift(sub(ph1, 2, -1)); F0s = sub(F0, 2); F1s = sub(F1, 2)

print("order N =", N)
A0l = add(add(f0, sc(F0s, 2)), sc(one(), -2)); A0r = add(ph0m, p0)
A1l = add(f1, shift(sc(F1s, 2), 1)); A1r = add(p1, ph1m, -1)
print("(A0) f0(q) + 2F0(q^2) - 2 = phi0(-q^2) + psi0(q)   :", A0l == A0r, "   both = theta_4(0,-q) G(-q):", A0l == mul(th3, sub(G, 1, -1)))
print("(C0) chi0(q) = 2F0(q) - phi0(-q)                   :", c0 == add(sc(F0, 2), sub(ph0, 1, -1), -1))
print("(A1) f1(q) + 2q F1(q^2) = psi1(q) - phi1(-q^2)/q   :", A1l == A1r, "   both = theta_4(0,-q) H(-q):", A1l == mul(th3, sub(H, 1, -1)))
print("(C1) q chi1(q) = 2q F1(q) + phi1(-q)               :", shift(c1, 1) == add(shift(sc(F1, 2), 1), sub(ph1, 1, -1)))
print("(R4) psi0(q) - f0(q) + 2 = chi0(q^2)               :", add(add(p0, f0, -1), sc(one(), 2)) == sub(c0, 2))
print("(R3) psi1(q) - f1(q) = q chi1(q^2)                 :", add(p1, f1, -1) == shift(sub(c1, 2), 1))
print("further relations of the two families:")
print("     2 phi0(-q^2) - f0(q) = theta_4(0,q) G(q)       :", add(sc(ph0m, 2), f0, -1) == mul(th, G))
print("     f1(q) + 2 phi1(-q^2)/q = theta_4(0,q) H(q)     :", add(f1, sc(ph1m, 2)) == mul(th, H))
print("     psi0(q) - F0(q^2) + 1 = q psi(q^2) H(q^4)      :", add(add(p0, F0s, -1), one()) == shift(mul(sub(ps, 2), sub(H, 4)), 1))
print("     psi1(q) - q F1(q^2) = psi(q^2) G(q^4)          :", add(p1, shift(F1s, 1), -1) == mul(sub(ps, 2), sub(G, 4)))
