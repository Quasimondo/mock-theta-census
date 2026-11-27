# Independent check of addendum-v17 (own code):
# (a) the explicit constants of Lemmas G, T and Corollary G; the infimum of |1 - e^eta| on |eta| >= 1, |Im eta| <= pi;
# (b) Proposition 2 / Corollary R off the real axis: |P*J - sum_{n<N} (-1)^n tau_n| for N = 0..12 at |t| = 0.2, arg t = -60 deg.
import sys
import mpmath as mp
sys.path.insert(0, '../v13')
from check_borel import PJ
mp.mp.dps = 25
c3 = 2*mp.exp(-0.5)*mp.sin(0.5); L = 2 + 1/c3; nu = mp.log(2) - L*mp.mpf('0.04')
mc = 1 + mp.exp(-0.3)*mp.cos(0.3); lam = mp.exp(0.3)/mc
print(f'(a) c3 = {mp.nstr(c3,5)}, L = {mp.nstr(L,5)}, nu = {mp.nstr(nu,4)}, m_c = {mp.nstr(mc,5)}, lambda = {mp.nstr(lam,4)}, '
      f'A0 = {mp.nstr(mp.exp(0.3)/(2*(1-lam)),3)}, B0 = {mp.nstr(1/(0.3*(1-lam)),4)}, B1 = {mp.nstr(2/(0.04*(1-lam)),4)}, 1/(nu c1) = {mp.nstr(1/(nu*0.04),3)}')
inf = min(abs(1 - mp.exp(mp.mpc(a, b))) for a in [x/50 for x in range(-150, 151)] for b in [y/50 for y in range(-157, 158)] if abs(mp.mpc(a, b)) >= 1)
print(f'    inf |1 - e^eta| on the grid = {mp.nstr(inf, 6)} (claimed bound 0.5816; true value 1 - 1/e = {mp.nstr(1 - mp.exp(-1), 6)})')
T = mp.mpf('0.2'); t = T*mp.expjpi(-mp.mpf(1)/3); x = 1j*t/(2*mp.pi); q = mp.exp(-t)
v = PJ(x)
s = 0; P = mp.mpc(1); rows = []
for n in range(13):
    rows.append((n, abs(v - s)))
    s += (-1)**n * q**(n*n+n) / (2*P); P *= 1 + q**(n+1)
print(f'(b) |t| = 0.2, arg t = -60 deg: P*J = {mp.nstr(v, 12)}')
print('    N, |P*J - sum_{n<N}(-1)^n tau_n|, ratio to previous:')
for i, (n, r) in enumerate(rows):
    print(f'      {n:2d}  {mp.nstr(r, 3):>9}  {mp.nstr(r/rows[i-1][1], 3) if i else ""}')
