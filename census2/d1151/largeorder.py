"""Large-order behaviour of the exact expansion d1151(e^{-t}) = sum_k a_k t^k (a_k.json from expand_q1.py).

Fits  a_k ~ K (-1)^k Gamma(k + beta) A^{-k-beta} (1 + s1 A/(k-1/2) + s2 A^2/((k-1/2)(k-3/2)) + ...)
and identifies the constants. Expected output (to about 10-12 digits):
    A = pi^2/60,  beta = 1/2,  pi K^2 = 8(5+sqrt5)/25,  s1 = (35 + 6 sqrt5)/300
The identifications are numerical (PSLQ), not proved.
"""
import json, mpmath as mp
from fractions import Fraction as Fr
mp.mp.dps = 60
a = [Fr(x) for x in json.load(open('a_k.json'))]
af = [mp.mpf(x.numerator) / x.denominator for x in a]
def rich(vals, ks, p):
    M = mp.matrix([[1] + [1 / mp.mpf(k) ** j for j in range(1, p)] for k in ks]); return mp.lu_solve(M, mp.matrix(vals))[0]
r = [af[k + 1] / af[k] for k in range(len(af) - 1)]
d = [r[k + 1] - r[k] for k in range(len(r) - 1)]
ks = list(range(len(d) - 6, len(d)))
M = mp.matrix([[1] + [1 / mp.mpf(k + 1) ** j for j in range(2, 7)] for k in ks])
A_est = -1 / mp.lu_solve(M, mp.matrix([d[k] for k in ks]))[0]
print('A    =', mp.nstr(A_est, 14), '   pi^2/60 =', mp.nstr(mp.pi ** 2 / 60, 14))
A = mp.pi ** 2 / 60
bs = [-A * af[k + 1] / af[k] - k for k in range(len(af) - 1)]
ks = list(range(len(bs) - 6, len(bs)))
print('beta =', mp.nstr(rich([bs[k] for k in ks], [k + 1 for k in ks], 6), 12))
b = mp.mpf(1) / 2
Ks = [af[k] * (-1) ** k * A ** (k + b) / mp.gamma(k + b) for k in range(len(af))]
ks = list(range(len(Ks) - 7, len(Ks)))
K = rich([Ks[k] for k in ks], [k + b for k in ks], 7)
print('pi K^2 =', mp.nstr(mp.pi * K ** 2, 14), '   8(5+sqrt5)/25 =', mp.nstr(8 * (5 + mp.sqrt(5)) / 25, 14))
K = mp.sqrt(8 * (5 + mp.sqrt(5)) / 25 / mp.pi)
ks = list(range(30, 38))
rows = [[1] + [A ** j / mp.fprod([k - b - i for i in range(j)]) for j in range(1, 8)] for k in ks]
sol = mp.lu_solve(mp.matrix(rows), mp.matrix([af[k] * (-1) ** k * A ** (k + b) / (K * mp.gamma(k + b)) for k in ks]))
print('leading =', mp.nstr(sol[0], 12), '  s1 =', mp.nstr(sol[1], 12), '   (35+6 sqrt5)/300 =', mp.nstr((35 + 6 * mp.sqrt(5)) / 300, 12))
print('PSLQ on (s1, 1, sqrt5):', mp.pslq([sol[1], 1, mp.sqrt(5)], maxcoeff=2000, maxsteps=10 ** 6, tol=mp.mpf(10) ** -9))
