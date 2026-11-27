# Independent check of addendum-v14:
# (a) the facts behind Lemma 3: U(0) = 0, U'(0) = -ln 2, U'' < 0, reflection formula, theta0 = arctan(ln2/pi);
# (b) Theorem 1 numerically: (P*J - sum_{k<K} a_k t^k)/(a_K t^K) - 1 shrinks like |t|, at arg t = -45 deg (own quadrature).
import sys, json
from fractions import Fraction as Fr
import mpmath as mp
sys.path.insert(0, '../v13')
mp.mp.dps = 30
U = lambda tau: -tau**2/2 + mp.polylog(2, -mp.exp(tau)) + mp.pi**2/12
print('(a) U(0) =', mp.nstr(U(0), 3), '  U\'(0) + ln2 =', mp.nstr(mp.diff(U, 0) + mp.log(2), 3),
      '  max U\'\' on [-8, 8] =', mp.nstr(max(mp.diff(U, mp.mpf(s)/4, 2) for s in range(-32, 33)), 4))
print('    reflection defect at tau = 0.7:', mp.nstr(mp.polylog(2, -mp.exp(0.7)) + mp.polylog(2, -mp.exp(-0.7)) + mp.pi**2/6 + mp.mpf(0.7)**2/2, 3),
      '  theta0 =', mp.nstr(mp.degrees(mp.atan(mp.log(2)/mp.pi)), 6), 'deg')
print('    max of U at ln((sqrt5-1)/2):', mp.nstr(U(mp.log((mp.sqrt(5)-1)/2)) - mp.pi**2/60, 3))
coeffs = [Fr(c) for c in json.load(open('../v13/phi_exact.json'))['psi1']]
a = [mp.mpf(c.numerator)/c.denominator for c in coeffs[:12]]
from check_borel import PJ
mp.mp.dps = 25
print('(b) arg t = -45 deg:')
for K in (1, 3):
    row = []
    for T in (mp.mpf('0.2'), mp.mpf('0.1'), mp.mpf('0.05')):
        t = T*mp.expjpi(-mp.mpf(1)/4); x = 1j*t/(2*mp.pi)
        v = PJ(x)
        r = (v - mp.fsum(a[k]*t**k for k in range(K))) / (a[K]*t**K)
        row.append(mp.nstr(abs(r - 1), 3))
    print(f'    K = {K}: |ratio - 1| at |t| = 0.2, 0.1, 0.05: {", ".join(row)}', flush=True)
