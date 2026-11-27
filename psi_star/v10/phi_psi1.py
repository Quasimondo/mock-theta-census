# Exact perturbative series of psi1* at q = e^{-t}: by addendum-v3 sec. 1 it is that of f1*/2,
# f1*(e^{-t}) = sum_n (-1)^n e^{-t(n^2+n)} / prod_{k<=n}(1+e^{-kt}),  prod 2/(1+e^{-kt}) = e^{t n(n+1)/4}/prod cosh(kt/2)
# term_n = (-1/2)^n exp(-3t(n^2+n)/4 - sum_k log cosh(k t/2));  sum_n n^d x^n with x = -1/2.
from fractions import Fraction as Fr
from math import comb, factorial
import sys
sys.path.insert(0, '/home/claude/v9check')
from phiJ import bern, powsum, padd, pmul, psc, B
NT = 40
E = [[Fr(0)] for _ in range(NT + 1)]
E[1] = [Fr(0), Fr(-3, 4), Fr(-3, 4)]
for n in range(1, NT // 2 + 1):
    a = Fr(2 ** (2 * n) * (2 ** (2 * n) - 1)) * B[2 * n] / (2 * n * factorial(2 * n)) / Fr(2 ** (2 * n))   # log cosh(k t/2)
    E[2 * n] = padd(E[2 * n], psc(powsum(2 * n), -a))
X = [[Fr(1)]] + [None] * NT
for n in range(1, NT + 1):
    acc = [Fr(0)]
    for k in range(1, n + 1):
        if any(E[k]): acc = padd(acc, psc(pmul(E[k], X[n - k]), Fr(k)))
    X[n] = psc(acc, Fr(1, n))
x = Fr(-1, 2); maxd = max(len(v) for v in X)
A = [Fr(0)] * maxd
for d in range(maxd):                       # (1-x) A(d) = [d=0] + x sum_{j<d} C(d,j) A(j)
    A[d] = ((1 if d == 0 else 0) + x * sum(comb(d, j) * A[j] for j in range(d))) / (1 - x)
PHI = [sum(c * A[d] for d, c in enumerate(X[n])) / 2 for n in range(NT + 1)]
if __name__ == '__main__':
    print('Phi =', ', '.join(str(PHI[n]) for n in range(6)), '  (addendum-v3: 1/3, 1/9, 5/162, 5/486, 29/17496, -1937/262440)')
