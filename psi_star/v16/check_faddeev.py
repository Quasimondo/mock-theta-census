# Independent check of addendum-v16 with own code:
# (a) Lemma Phi (i): Faddeev's integral Lambda(w) = int_{R+i0} e^{-i w s}/(4 s sinh(pi s) sinh(pi x s)) ds equals log of the
#     product formula (|q| < 1: v12 form; |q| > 1: inverted products as in v13), w = 2 pi b y.
# (b) inversion relation at real x: Lambda(w) + Lambda(-w) = -t/24 + pi^2/(6t) + w^2/(2t).
# (c) the Borel sum of Phi at real x (|t| = 0.3) against the value quoted in sec. 6 for P*J there.
import sys, json
import mpmath as mp
from fractions import Fraction as Fr
sys.path.insert(0, '../v13')
mp.mp.dps = 25
def Lam(w, x, eta=mp.mpf('0.5')):
    f = lambda a: (lambda s: mp.exp(-1j*w*s) / (4*s*mp.sinh(mp.pi*s)*mp.sinh(mp.pi*x*s)))(a + 1j*eta)
    return mp.quad(f, [-mp.inf, -20, -5, -1, 0, 1, 5, 20, mp.inf])
def poch(a, q):
    if abs(q) > 1: return 1/poch(a/q, 1/q)
    s = mp.mpc(1); k = 0
    while True:
        f = a*q**k; s *= 1 - f
        if abs(f) < mp.mpf(10)**(-mp.mp.dps - 5): return s
        k += 1
def logPhi_prod(w, x):
    b = mp.sqrt(x); y = w/(2*mp.pi*b); q = mp.expjpi(2*x); p = mp.expjpi(-1/x)
    return mp.log(poch(-mp.expjpi(x)*mp.exp(2*mp.pi*b*y), q) / poch(-p*mp.exp(2*mp.pi*y/b), p**2))
print('(a) Faddeev integral vs product formula (difference mod 2 pi i):')
for xdeg in (20, -20):
    x = mp.mpf('0.2')*mp.expjpi(mp.mpf(xdeg)/180)
    for w in (mp.mpc(0.4, 0.6), mp.mpc(-1.0, -1.5)):
        d = Lam(w, x) - logPhi_prod(w, x); d -= 2j*mp.pi*mp.nint(mp.im(d)/(2*mp.pi))
        print(f'    arg x = {xdeg:+d} deg (|q| {"<" if xdeg > 0 else ">"} 1), w = {mp.nstr(w,3)}: |difference| = {mp.nstr(abs(d), 3)}')
x = mp.mpf('0.3')/(2*mp.pi); t = -2j*mp.pi*x                 # real x, |q| = 1
print('(b) inversion relation at real x = %s:' % mp.nstr(x, 6))
for w in (mp.mpc(0.3, 0.5), mp.mpc(-0.8, 1.9)):
    d = Lam(w, x) + Lam(-w, x) - (-t/24 + mp.pi**2/(6*t) + w**2/(2*t)); d -= 2j*mp.pi*mp.nint(mp.im(d)/(2*mp.pi))
    print(f'    w = {mp.nstr(w,3)}: |defect| = {mp.nstr(abs(d), 3)}')
from check_borel import borel_sum
S = borel_sum(t)
quoted = mp.mpc('0.330580057167', '-0.033046552525')
print(f'(c) Borel sum of Phi at t = -0.3i (own Pade-Laplace): {mp.nstr(S, 15)};  addendum P*J: {mp.nstr(quoted, 12)};  |difference| = {mp.nstr(abs(S - quoted), 3)}')
