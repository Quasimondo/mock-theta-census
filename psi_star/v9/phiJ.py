# Exact perturbative series of J(e^{-tau}) = sum_m e^{-tau(2m^2+m)} / prod_{k=1}^m (1 + e^{-2k tau}).
# prod 2/(1+e^{-2k tau}) = e^{tau m(m+1)/2} / prod cosh(k tau), so
# term_m = 2^{-m} exp(-tau (3m^2+m)/2 - sum_k log cosh(k tau)),  log cosh y = sum_{n>=1} 2^{2n}(2^{2n}-1) B_{2n} y^{2n} / (2n (2n)!)
# Coefficients are sums over m of 2^{-m} * polynomial(m): sum_m m^d 2^{-m} = A(d), A(d) = 2[d=0] + sum_{j<d} C(d,j) A(j).
from fractions import Fraction as Fr
from math import comb, factorial
NT = 46
def bern(n):
    B = [Fr(1)]
    for m in range(1, n + 1):
        B.append(-sum(comb(m + 1, k) * B[k] for k in range(m)) / (m + 1))
    return B
B = bern(2 * NT + 2)
def powsum(p):          # S_p(m) = sum_{k=1}^m k^p as a polynomial in m (Faulhaber, B_1 = +1/2 convention)
    c = [Fr(0)] * (p + 2)
    for j in range(p + 1):
        bj = B[j] if j != 1 else Fr(1, 2)
        c[p + 1 - j] += Fr(comb(p + 1, j)) * bj / (p + 1)
    return c
def padd(a, b):
    n = max(len(a), len(b)); return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]
def pmul(a, b):
    c = [Fr(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b): c[i + j] += x * y
    return c
def psc(a, s): return [s * x for x in a]
# exponent E(tau) = sum_n tau^n e_n(m), e_n polynomials in m
E = [[Fr(0)] for _ in range(NT + 1)]
E[1] = [Fr(0), Fr(-1, 2), Fr(-3, 2)]
for n in range(1, NT // 2 + 1):
    a = Fr(2 ** (2 * n) * (2 ** (2 * n) - 1)) * B[2 * n] / (2 * n * factorial(2 * n))
    if 2 * n <= NT: E[2 * n] = padd(E[2 * n], psc(powsum(2 * n), -a))
# exp of the series: Exp' = E' Exp  ->  n X_n = sum_{k=1}^n k E_k X_{n-k}
X = [[Fr(1)]] + [None] * NT
for n in range(1, NT + 1):
    acc = [Fr(0)]
    for k in range(1, n + 1):
        if any(E[k]): acc = padd(acc, psc(pmul(E[k], X[n - k]), Fr(k)))
    X[n] = psc(acc, Fr(1, n))
maxd = max(len(x) for x in X)
A = [Fr(0)] * maxd
for d in range(maxd):
    A[d] = (2 if d == 0 else 0) + sum(comb(d, j) * A[j] for j in range(d))
PHI = [sum(c * A[d] for d, c in enumerate(X[n])) for n in range(NT + 1)]
if __name__ == '__main__':
    print('Phi_J =', ', '.join(str(PHI[n]) for n in range(4)), ' (addendum-v9: 2 - 10 tau + 183 tau^2 - ...)')
