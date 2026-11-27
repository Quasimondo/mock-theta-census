# Theorem 2 of addendum-v15, own quadrature (contour as in check_borel.PJ): |ratio - 1| for K = 1 should halve with |t|
import sys, json
from fractions import Fraction as Fr
import mpmath as mp
sys.path.insert(0, '../v13')
from check_borel import PJ
coeffs = [Fr(c) for c in json.load(open('../v13/phi_exact.json'))['psi1']]
a = [mp.mpf(c.numerator)/c.denominator for c in coeffs[:4]]
mp.mp.dps = 25
for th in (0, -150):
    row = []
    for T in (mp.mpf('0.2'), mp.mpf('0.1'), mp.mpf('0.05')):
        t = T*mp.expjpi(mp.mpf(th)/180); x = 1j*t/(2*mp.pi)
        v = PJ(x)
        row.append(mp.nstr(abs((v - a[0])/(a[1]*t) - 1), 3))
    print(f'theta = {th} deg: K = 1, |ratio - 1| at |t| = 0.2, 0.1, 0.05: {", ".join(row)}', flush=True)
