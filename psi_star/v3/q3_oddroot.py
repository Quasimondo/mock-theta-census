"""Q3: the odd-root rule  psi1*(zeta) = A/3,  psi0*(zeta) = A0/3  at roots zeta of odd order m.
 (a) finite identity A = B (reversal n -> m-1-n), 40 digits
 (b) f1*, f0* at the root: geometric with ratio -1/2 per m steps, sum = (2/3) B
 (c) size of theta4*R1 and theta4*R0 along the radius, against the elementary bound
 (d) direct radial values of psi* against A/3
"""
import mpmath as mp
from math import gcd
mp.mp.dps = 40
def finite_sums(h, m):
    z = mp.expjpi(mp.mpf(2*h)/m)
    A = B = A0 = B0 = 0; P = mp.mpf(1)
    for n in range(m):
        if n: P *= 1 + z**n
        A  += (-1)**n * z**(n*(n+1)//2) * P;       B  += (-1)**n * z**(n*n+n) / P
        A0 += (-1)**n * z**((n+1)*(n+2)//2) * P;   B0 += (-1)**n * z**(n*n) / P
    return z, A, B, A0, B0
print("(a) |A - B|, |A0 - B0| for all primitive roots of odd order m <= 25:")
worst = 0
for m in range(1, 26, 2):
    for h in range(1, m+1):
        if gcd(h, m) != 1: continue
        z, A, B, A0, B0 = finite_sums(h, m)
        worst = max(worst, abs(A-B), abs(A0-B0))
print("    largest difference:", mp.nstr(worst, 3))

def series(q, expo, sign, kind, tol):
    """sum sign^n q^expo(n) * P_n; kind: '+num' -> (-q)_n, '+den' -> 1/(-q)_n, '-den' -> 1/(q)_n"""
    s = 0; P = mp.mpf(1); n = 0; qn = mp.mpf(1)
    while True:
        if n:
            qn *= q
            P = P*(1+qn) if kind == '+num' else (P/(1+qn) if kind == '+den' else P/(1-qn))
        term = sign**n * q**expo(n) * P
        s += term
        if n > 50 and abs(term) < tol and abs(q)**expo(n) < tol: break
        n += 1
    return s
print("(b),(c),(d) along q = zeta e^(-t):")
print("    m  h   u=m^2 t   |f1* - (2/3)A|   |psi1* - A/3|   |2psi1* - f1* - th4 R1|   u*log|th4 R1|   (bound: < -pi^2/6 = -1.645)")
for (h, m) in [(0, 1), (1, 3), (1, 5), (2, 5), (1, 7), (3, 7), (2, 9)]:
    z, A, B, A0, B0 = finite_sums(h if m > 1 else 1, m)
    for u in (mp.mpf(1)/4, mp.mpf(1)/8, mp.mpf(1)/16, mp.mpf(1)/32):
        t = u/m**2
        mp.mp.dps = 60
        q = z*mp.exp(-t)
        tol = mp.mpf(10)**(-mp.mp.dps+5)
        f1 = series(q, lambda n: n*n+n, -1, '+den', tol)
        R1 = series(q, lambda n: n*n+n, -1, '-den', tol)
        ps = series(q, lambda n: n*(n+1)//2, -1, '+num', tol)
        th4 = mp.jtheta(4, 0, q) if m == 1 else mp.nsum(lambda n: (-1)**int(n)*q**(int(n)**2), [-mp.inf, mp.inf])
        print("   %2d %2d   %-7s   %-14s   %-13s   %-22s   %s" % (m, h, mp.nstr(u, 4), mp.nstr(abs(f1 - mp.mpf(2)/3*A), 3),
              mp.nstr(abs(ps - A/3), 3), mp.nstr(abs(2*ps - f1 - th4*R1), 3), mp.nstr(u*mp.log(abs(th4*R1)), 6)))
