#!/usr/bin/env python3
"""Addendum v29, later remark on the Corollary of section 1.

Exact integer power series to order N (default 300).

 (W0) f0 + 2 psi0 = 4q (q^4,q^16,q^20;q^20)_inf/(q^2;q^4)_inf + (q^2,q^3,q^5;q^5)_inf/(-q;q)_inf
 (W1) f1 + 2 psi1 = 4  (q^8,q^12,q^20;q^20)_inf/(q^2;q^4)_inf - (q,q^4,q^5;q^5)_inf/(-q;q)_inf
      (Watson 1937, section 7; in this form Mc Laughlin, arXiv:1906.11997, (4.1), (4.2))
 (L)  f1 psi0 - f0 psi1 = (f1 B0 - f0 B1)/2 = B1 psi0 - B0 psi1,   B_j = f_j + 2 psi_j   (algebra)
 (C)  f1 psi0 - f0 psi1 = -2q sum_{n in Z} (-1)^n q^{n(n+3)/2} s_n,  s_n = sum_{j<=n} (-1)^j (-q;q)_j
"""
import sys

N = int(sys.argv[1]) if len(sys.argv) > 1 else 300


def mul(a, b):
    c = [0]*(N + 1)
    for i, x in enumerate(a):
        if x:
            for j in range(N + 1 - i):
                if b[j]:
                    c[i + j] += x*b[j]
    return c


def inv(a):
    assert a[0] in (1, -1)
    b = [0]*(N + 1)
    b[0] = a[0]
    for n in range(1, N + 1):
        b[n] = -a[0]*sum(a[k]*b[n - k] for k in range(1, n + 1))
    return b


def one():
    return [1] + [0]*N


def fac(k, c=1):
    """1 - c q^k"""
    a = one()
    if k <= N:
        a[k] -= c
    return a


def shift(a, k):
    return [0]*k + a[:N + 1 - k] if k <= N else [0]*(N + 1)


def add(a, b, cb=1):
    return [x + cb*y for x, y in zip(a, b)]


def poch(residues, m, c=1):
    """prod over k >= 1, k mod m in residues, of (1 - c q^k)"""
    a = one()
    for k in range(1, N + 1):
        if k % m in residues:
            a = mul(a, fac(k, c))
    return a


# the four functions
f0 = [0]*(N + 1); f1 = [0]*(N + 1); p0 = [0]*(N + 1); p1 = [0]*(N + 1)
P = one()            # (-q;q)_n
n = 0
while n*(n + 1)//2 <= N:
    Pi = inv(P)
    f0 = add(f0, shift(Pi, n*n)); f1 = add(f1, shift(Pi, n*n + n))
    p0 = add(p0, shift(P, (n + 1)*(n + 2)//2)); p1 = add(p1, shift(P, n*(n + 1)//2))
    n += 1
    P = mul(P, fac(n, -1))

mq = poch(range(1), 1, -1)                       # (-q;q)_inf
d24 = inv(poch([2], 4))                          # 1/(q^2;q^4)_inf
B0 = add(shift([4*x for x in mul(poch([4, 16, 0], 20), d24)], 1), mul(poch([2, 3, 0], 5), inv(mq)))
B1 = add([4*x for x in mul(poch([8, 12, 0], 20), d24)], mul(poch([1, 4, 0], 5), inv(mq)), -1)
print("order N =", N)
print("(W0) f0 + 2 psi0 = theta products :", add(f0, [2*x for x in p0]) == B0)
print("(W1) f1 + 2 psi1 = theta products :", add(f1, [2*x for x in p1]) == B1)

w = add(mul(f1, p0), mul(f0, p1), -1)
h = add(mul(f1, B0), mul(f0, B1), -1)
print("(L)  f1 psi0 - f0 psi1 = (f1 B0 - f0 B1)/2 :", [2*x for x in w] == h)
print("(L)  f1 psi0 - f0 psi1 = B1 psi0 - B0 psi1  :", w == add(mul(B1, p0), mul(B0, p1), -1))

# twice the partial sums: 2 s_n
# negative side: (-q;q)_{-m} = q^{m(m-1)/2} / (2 (-q;q)_{m-1}),  m >= 1
M = 2
while (M - 1)*(M - 1) <= N:
    M += 1
term = {}                                        # 2 (-1)^j (-q;q)_j for j = -m
P = one()
for m in range(1, M + 2):
    term[-m] = [(-1)**m*x for x in shift(inv(P), m*(m - 1)//2)]
    P = mul(P, fac(m, -1))
s2 = {}
acc = [0]*(N + 1)
for m in range(M + 1, 0, -1):                    # tail below -M-1 has order > N after the prefactor
    acc = add(acc, term[-m])
    s2[-m] = acc[:]
rhs = [0]*(N + 1)
for m in range(1, M + 1):
    e = m*(m - 3)//2 + 1
    t = s2[-m]
    if e < 0:
        assert all(x == 0 for x in t[:(-e)])
        t = t[-e:] + [0]*(-e)
    else:
        t = shift(t, e)
    rhs = add(rhs, t, -(-1)**m)
acc = s2[-1][:]
P = one()
n = 0
while n*(n + 3)//2 + 1 <= N:
    acc = add(acc, [2*(-1)**n*x for x in P])
    rhs = add(rhs, shift(acc, n*(n + 3)//2 + 1), -(-1)**n)
    n += 1
    P = mul(P, fac(n, -1))
print("(C)  f1 psi0 - f0 psi1 = -2q sum (-1)^n q^{n(n+3)/2} s_n :", w == rhs)
print("first coefficients of -(f1 psi0 - f0 psi1):", [-x for x in w[:31]])
