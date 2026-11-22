"""Exact power-series checks (integer arithmetic, truncation at q^N).

(I)   psi0*(q) = 1/2 f0*(q) - 1/2 theta4(q) R0(q)
(II)  psi1*(q) = 1/2 f1*(q) + 1/2 theta4(q) R1(q)
with  f0* = sum (-1)^n q^{n^2}/(-q)_n,      f1* = sum (-1)^n q^{n^2+n}/(-q)_n
      R0  = sum (-1)^n q^{n^2}/(q)_n,       R1  = sum (-1)^n q^{n^2+n}/(q)_n
      theta4 = (q)_inf/(-q)_inf = sum (-1)^n q^{n^2}
"""
from flint import fmpz_poly, fmpz_series
import sys

N = int(sys.argv[1]) if len(sys.argv) > 1 else 400
fmpz_series  # noqa


def mul(a, b):
    c = a * b
    return fmpz_poly(c.coeffs()[:N])


def mono(e, c=1):
    if e >= N:
        return fmpz_poly([0])
    return fmpz_poly([0] * e + [c])


def inv(a):
    # power series inverse mod q^N, a(0)=1
    co = a.coeffs() + [0] * N
    b = [0] * N
    b[0] = 1
    for n in range(1, N):
        s = 0
        for k in range(1, n + 1):
            if co[k]:
                s += co[k] * b[n - k]
        b[n] = -s
    return fmpz_poly(b)


def euler_sum(expo, sign, fac, numer=True):
    """sum_n sign^n q^{expo(n)} prod_{k<=n} fac(k)^{+-1}"""
    tot = fmpz_poly([0])
    P = fmpz_poly([1])
    n = 0
    while expo(n) < N:
        if n > 0:
            P = mul(P, fac(n))
        term = P if numer else inv(P)
        tot += mul(mono(expo(n), sign ** n), term)
        n += 1
    return fmpz_poly(tot.coeffs()[:N])


plus = lambda k: fmpz_poly([1] + [0] * (k - 1) + [1]) if k < N else fmpz_poly([1])
minus = lambda k: fmpz_poly([1] + [0] * (k - 1) + [-1]) if k < N else fmpz_poly([1])

psi1s = euler_sum(lambda n: n * (n + 1) // 2, -1, plus)
psi0s = euler_sum(lambda n: (n + 1) * (n + 2) // 2, -1, plus)
psi1 = euler_sum(lambda n: n * (n + 1) // 2, 1, plus)
psi0 = euler_sum(lambda n: (n + 1) * (n + 2) // 2, 1, plus)
f0s = euler_sum(lambda n: n * n, -1, plus, numer=False)
f1s = euler_sum(lambda n: n * n + n, -1, plus, numer=False)
R0 = euler_sum(lambda n: n * n, -1, minus, numer=False)
R1 = euler_sum(lambda n: n * n + n, -1, minus, numer=False)
th4 = fmpz_poly([0])
n = 0
while n * n < N:
    th4 += mono(n * n, (2 if n else 1) * (-1) ** n)
    n += 1

print("psi1* :", psi1s.coeffs()[:30])
print("psi0*/q:", psi0s.coeffs()[1:31])

d1 = 2 * psi1s - f1s - mul(th4, R1)
d0 = 2 * psi0s - f0s + mul(th4, R0)
print("N =", N)
print("(II) 2psi1* - f1* - th4*R1 == 0 :", all(c == 0 for c in d1.coeffs()))
print("(I)  2psi0* - f0* + th4*R0 == 0 :", all(c == 0 for c in d0.coeffs()))

# NOTE (repo edit): the original script also tested "f0 + 2 psi0 = theta4(q^5) G" and
# "f1 + 2 psi1 = theta4(q^5) H". Those are not the classical identities (they print False) and
# have been removed. The correct unsigned statements are B(f0) = f0 + 2 psi0 and
# B(f1) = f1 + 2 psi1 (bilateral series, modular); see Bajpai-Kimport-Liang-Ma-Ricci,
# arXiv:2202.12141, and Mc Laughlin, arXiv:1906.11997, Prop. 4.1.

import json
json.dump({"psi1s": [int(c) for c in psi1s.coeffs()], "psi0s": [int(c) for c in psi0s.coeffs()],
           "f0s": [int(c) for c in f0s.coeffs()], "f1s": [int(c) for c in f1s.coeffs()],
           "R0": [int(c) for c in R0.coeffs()], "R1": [int(c) for c in R1.coeffs()]},
          open("series.json", "w"))
