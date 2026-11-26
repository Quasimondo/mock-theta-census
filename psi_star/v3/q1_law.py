"""Q1: a modularity law of Zagier type for psi1*, psi0* (numerical test).

For X -> infinity through rationals with X = rho mod 2 (rho rational), t = -2 pi i / X:
   psi1*(e^(2 pi i/X)) = Cbar * t^(-1/2) * e^(-Sbar/t) * W(e^(-pi i X)) * (1 + s1bar t + O(t^2)) ,
   psi0*(e^(2 pi i/X)) = e^(2 pi i/3) * (the same, with s1bar replaced by its psi0* value) ,
 Sbar = pi^2/3 + i Cl2(pi/3),  Cbar = sqrt(pi/sqrt3) e^(5 pi i/12),
 W(xi) = sum_j (-1)^j xi^(j^2+j) (xi;xi^2)_j   (finite at odd order; geometric-series value at even order).
Values of psi* at roots of unity are exact finite sums (jval.py); the limit is taken by Neville extrapolation in 1/X.
Usage: python3 q1_law.py [smax]   (denominators of rho up to smax, default 9)
"""
import mpmath as mp, sys
from fractions import Fraction as Fr
from math import gcd
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from jval import J
from wval import W
mp.mp.dps = 60
smax = int(sys.argv[1]) if len(sys.argv) > 1 else 9
Cl2 = mp.clsin(2, mp.pi/3); Sb = mp.pi**2/3 + 1j*Cl2
Cb = mp.sqrt(mp.pi/mp.sqrt(3))*mp.expjpi(mp.mpf(5)/12)
om = mp.expjpi(mp.mpf(2)/3)
s1b = (mp.mpf(1)/4 + 1j*mp.sqrt(3)/54, -mp.mpf(1)/12 + 1j*mp.sqrt(3)/54)      # conjugates of s1 from Q2
Ms = [200, 240, 300, 360, 440, 540, 660, 800, 960, 1160, 1400, 1700, 2000, 2400]
def neville(hs, vs):
    n = len(hs); T = vs[:]
    for k in range(1, n):
        T = [(hs[i+k]*T[i] - hs[i]*T[i+1])/(hs[i+k]-hs[i]) for i in range(n-k)]
    return T[0]
worst = [0, 0, 0, 0]; n = 0
rows = []
for s in range(1, smax+1):
    for r in range(s):
        if gcd(r, s) != 1: continue
        for odd in (0, 1):
            rho = Fr(r, s) + odd
            xi = -rho/2                                   # W at e^(2 pi i * (-rho/2)) = e^(-pi i rho)
            Wv = W(xi.numerator, xi.denominator)
            hs, v1, v0 = [], [], []
            for M in Ms:
                num, den = s, s*(M+odd) + r               # 1/X = s/(s(M+odd)+r)
                j1, j0 = J(num, den)
                X = mp.mpf(den)/s; t = -2j*mp.pi/X
                lead = mp.exp(-Sb/t)/mp.sqrt(t)
                hs.append(1/X); v1.append(j1/lead); v0.append(j0/lead)
            L1, L0 = neville(hs, v1), neville(hs, v0)
            d1, d0 = abs(L1 - Cb*Wv), abs(L0 - om*Cb*Wv)
            # first correction: (psi/(lead*Cb*W) - 1)/t  ->  s1bar
            e1 = e0 = None
            if abs(Wv) > mp.mpf(10)**-10:
                c1 = neville(hs, [(v1[i]/(Cb*Wv) - 1)/(-2j*mp.pi*hs[i]) for i in range(len(hs))])
                c0 = neville(hs, [(v0[i]/(om*Cb*Wv) - 1)/(-2j*mp.pi*hs[i]) for i in range(len(hs))])
                e1, e0 = abs(c1 - s1b[0]), abs(c0 - s1b[1])
                worst[2] = max(worst[2], e1); worst[3] = max(worst[3], e0)
            worst[0] = max(worst[0], d1/max(abs(Wv), 1)); worst[1] = max(worst[1], d0/max(abs(Wv), 1)); n += 1
            rows.append((rho, Wv, d1, d0, e1, e0))
print("rho      W(e^(-pi i rho))                               |L1 - Cbar W|   |L0 - w Cbar W|   |corr1 - s1bar|  |corr0 - s1bar'|")
for rho, Wv, d1, d0, e1, e0 in rows[:16]:
    print("%-7s  %-45s  %-13s   %-15s   %-14s   %s" % (rho, mp.nstr(Wv, 14), mp.nstr(d1, 2), mp.nstr(d0, 2), mp.nstr(e1, 2) if e1 is not None else "-", mp.nstr(e0, 2) if e0 is not None else "-"))
print("...")
print(f"{n} values of rho (denominators <= {smax}, both classes mod 2):")
print("  largest |L1 - Cbar W| / max(|W|,1):", mp.nstr(worst[0], 3), "   largest |L0 - e^(2 pi i/3) Cbar W| / max(|W|,1):", mp.nstr(worst[1], 3))
print("  largest deviation of the first correction from s1bar:", mp.nstr(worst[2], 3), "(psi1*)", mp.nstr(worst[3], 3), "(psi0*)")
