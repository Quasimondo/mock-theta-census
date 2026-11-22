"""Radial limits of psi1*, psi0* at odd-order roots of unity.

Compares high-precision values at q = zeta*exp(-t), extrapolated to t -> 0
by polynomial (Richardson) extrapolation, with the conjectured closed form
    lim = (1/3) * sum_{n<m} a_n(zeta)    (m odd).
Expected output: every diff around 1e-14 (limited by the extrapolation).
Run time: a few minutes.
"""
# Radial limits of psi1*, psi0* at odd-order roots of unity:
# compare high-precision values at q = zeta*exp(-t), extrapolated t->0,
# with the conjectured closed form  lim = (1/3) * sum_{n<m} term_n(zeta).
import mpmath as mp, math
def psi(q, off, N):
    t = q**off; s = t; qn = mp.mpf(1)
    for n in range(1, N):
        qn *= q
        t *= -qn**(1) * (1 + qn) * (q if off else 1); s += t
    return s
def closed(z, off, m):
    t = z**off; s = t; zn = mp.mpf(1)
    for n in range(1, m):
        zn *= z; t *= -zn * (1 + zn) * (z if off else 1); s += t
    return s / 3
roots = [(h, m) for m in (1, 3, 5, 7, 9, 11, 13, 15) for h in range(m) if math.gcd(h, m) == 1 and 2*h <= m]
SC = lambda m: [0.02/(m**1.5*1.5**k) for k in range(9)]
worst = 0
for (h, m) in roots:
    ts = SC(m)
    for off, nm in ((0, 'psi1*'), (1, 'psi0*')):
        vals = []
        for t in ts:
            mp.mp.dps = int(math.pi**2/(12*t*m*m)/2.3) + 50
            z = mp.expjpi(mp.mpf(2*h)/m); q = z*mp.exp(-t)
            N = int(math.sqrt(2*mp.mp.dps*2.31/t)) + 10
            vals.append(psi(q, off, N))
        mp.mp.dps = 30
        # Richardson / polynomial extrapolation in t to t=0
        ext = mp.mpf(0)
        n = len(ts)
        for i in range(n):
            w = mp.mpf(1)
            for j in range(n):
                if j != i: w *= (0 - ts[j]) / (ts[i] - ts[j])
            ext += w * vals[i]
        z = mp.expjpi(mp.mpf(2*h)/m); c = closed(z, off, m)
        d = abs(ext - c); worst = max(worst, d)
        print(f"{h}/{m:<3} {nm}  extrapolated {mp.nstr(ext,12):>34}  closed {mp.nstr(c,12):>34}  diff {mp.nstr(d,2)}")
print("worst diff", mp.nstr(worst, 3))
