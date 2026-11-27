"""(a) exact, to r^299:   psi(r) (r^2;r^4)_inf sum_m r^(m(m+1)/2)/(r^2;r^2)_m  =  JJ(r) + KK(r),
        JJ(r) = sum_{m>=0} r^(2m^2+m)/(-r^2;r^2)_m ,   KK(r) = 2 sum_{n>=0} r^((n+1)^2) (-r^2;r^2)_n = sum_{n>=1} r^(n^2) (-1;r^2)_n
    (JJ + KK is the bilateral sum over all integers n of r^(n^2) (-1;r^2)_n)
(b) the measured limits of fstar.py against JJ, KK at the 20 points, and their first corrections."""
import sys, os, json; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from mpmath import mp, mpf, mpc, exp, pi, nstr, sqrt
N = 300
def mul(a, b):
    c = [0]*N
    for i, x in enumerate(a):
        if x:
            for j in range(N - i):
                if b[j]: c[i + j] += x*b[j]
    return c
def div1p(a, e, sign):            # a / (1 + sign r^e)
    c = [0]*N
    for i in range(N): c[i] = a[i] - (sign*c[i - e] if i >= e else 0)
    return c
def mul1p(a, e, sign): return [a[i] + (sign*a[i - e] if i >= e else 0) for i in range(N)]
one = [1] + [0]*(N - 1)
psi = [0]*N; k = 0
while k*(k + 1)//2 < N: psi[k*(k + 1)//2] = 1; k += 1
pr = one[:]
for e in range(2, N, 4): pr = mul1p(pr, e, -1)
E = [0]*N; den = one[:]; m = 0
while m*(m + 1)//2 < N:
    for i in range(N - m*(m + 1)//2): E[i + m*(m + 1)//2] += den[i]
    m += 1; den = div1p(den, 2*m, -1)
lhs = mul(mul(psi, pr), E)
JJ = [0]*N; KK = [0]*N; dn = one[:]; nu = one[:]; m = 0
while 2*m*m + m < N or (m + 1)**2 < N:
    e = 2*m*m + m
    for i in range(N - e) if e < N else []: JJ[i + e] += dn[i]
    e = (m + 1)**2
    for i in range(N - e) if e < N else []: KK[i + e] += 2*nu[i]
    m += 1; dn = div1p(dn, 2*m, +1); nu = mul1p(nu, 2*m, +1)
print("(a) bilateral identity to r^299:", lhs == [x + y for x, y in zip(JJ, KK)])
print("    JJ =", JJ[:16]); print("    KK =", KK[:16])

mp.dps = 30
res = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "fstar_circle.json")))
dev = {"F0": 0, "G0": 0, "qF1": 0, "mG1": 0}; ratios = {k: [] for k in dev}
for key, v in res.items():
    X0 = Fr(key); XX = mpc(mpf(X0.numerator)/X0.denominator, -mpf(4)/5)
    r = exp(-1j*pi*XX/2)
    J = mpc(0); K = mpc(0); dn = mpc(1); nu = mpc(1)
    for m in range(40):
        J += r**(2*m*m + m)/dn; K += 2*r**((m + 1)**2)*nu
        dn *= 1 + r**(2*m + 2); nu *= 1 + r**(2*m + 2)
    for nm, ex in (("F0", J), ("G0", K), ("qF1", J), ("mG1", K)):
        lim = mpc(mpf(v[nm][0]), mpf(v[nm][1])); c1 = mpc(mpf(v[nm][3]), mpf(v[nm][4]))
        dev[nm] = max(dev[nm], abs(lim - ex)/abs(ex)); ratios[nm].append(c1/ex)
print("(b) largest relative deviation of the measured limits from JJ (F0*, q F1*) and KK (G0, -G1):", {k: nstr(x, 3) for k, x in dev.items()})
for nm in dev:
    mean = sum(ratios[nm])/len(ratios[nm])
    print(f"    first correction / leading function for {nm}: {nstr(mean, 12)}  (spread {nstr(max(abs(x - mean) for x in ratios[nm]), 2)})")
print("    for comparison: 1/12 =", nstr(mpf(1)/12, 12), " 1/(9 sqrt3) =", nstr(1/(9*sqrt(3)), 12), " 1/24 =", nstr(mpf(1)/24, 12))
