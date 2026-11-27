# Independent check of addendum-v13 sec. 3: P*J(x) = Borel sum of Phi in the direction of t, at |t| = 0.3.
#   own Borel-Pade (mpmath pade + Laplace), own quadrature of J with Phi_b continued to |q| > 1 by inverted products.
import json, sys, glob
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 40
coeffs = [Fr(c) if isinstance(c, str) else Fr(c[0], c[1]) for c in json.load(open('phi_exact.json'))['psi1']]
def borel_sum(t, N=50):
    with mp.workdps(300):
        bc = [mp.mpf(c.numerator)/c.denominator / mp.factorial(k) for k, c in enumerate(coeffs[:2*N+1])]
        P, Q = mp.pade(bc, N, N)
    P = [+c for c in P]; Q = [+c for c in Q]
    Bf = lambda xi: mp.polyval(P[::-1], xi) / mp.polyval(Q[::-1], xi)
    return mp.quad(lambda r: mp.exp(-r) * Bf(t*r), [0, 5, 20, 60, mp.inf])     # direction theta = arg t
def poch(a, q):                       # (a;q)_inf for |q| < 1, continued by 1/(a/q;1/q)_inf for |q| > 1
    if abs(q) > 1: return 1/poch(a/q, 1/q)
    s = mp.mpc(1); k = 0
    while True:
        f = a*q**k; s *= 1 - f
        if abs(f) < mp.mpf(10)**(-mp.mp.dps-5): return s
        k += 1
def PJ(x, c=mp.mpf('0.5'), sign=-1, R=8, pieces=32):
    b = mp.sqrt(x); q = mp.expjpi(2*x); p = mp.expjpi(-1/x); t = -2j*mp.pi*x
    Phi = lambda y: poch(-mp.expjpi(x)*mp.exp(2*mp.pi*b*y), q) / poch(-p*mp.exp(2*mp.pi*y/b), p**2)
    f = lambda y: mp.exp(-mp.pi*1j*y**2 + mp.pi*y/b) / (Phi(y)*(1 + mp.exp(2*mp.pi*y/b)))
    delta = c*(mp.pi/2 - abs(mp.arg(b))); e = mp.expj(sign*delta)
    pts = [-R + 2*R*k/pieces for k in range(pieces+1)]
    J = mp.quad(lambda r: f(r*e)*e, pts)
    return 2**mp.mpf(-0.5) * mp.exp(mp.pi**2/(12*t) + t/6) / b * J
if __name__ == "__main__":
  for phi_deg in (-80, -135):
      t = mp.mpf('0.3') * mp.expjpi(mp.mpf(phi_deg)/180); x = 1j*t/(2*mp.pi)
      S = borel_sum(t)
      v1 = PJ(x, mp.mpf('0.5')); v2 = PJ(x, mp.mpf('0.65'))
      print(f'arg t = {phi_deg} deg, arg x = {mp.nstr(mp.degrees(mp.arg(x)),4)} deg, |q| = {mp.nstr(abs(mp.expjpi(2*x)),4)}:')
      print(f'   Borel sum {mp.nstr(S, 20)}')
      print(f'   P*J       {mp.nstr(v1, 20)}   (contour c=0.5 vs 0.65 differ by {mp.nstr(abs(v1-v2), 2)})   |P*J - Borel| = {mp.nstr(abs(v1-S), 3)}', flush=True)
