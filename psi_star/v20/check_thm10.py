# Independent check of addendum-v20, Theorem 10 (own code): exact finite sums psi1*(q0), psi0*(q0) at
# q0 = e^{2 pi i/X}, X even, against 3^{-1/4} (X/2)^{1/2} e^{Cl2 X/2pi} e^{2pi i/3 - i pi X/6} (psi0*: times e^{2 pi i/3}).
# Prints X (ratio - 1) and a Richardson extrapolation in 1/X.
import mpmath as mp
Cl2 = None
def psis(X):
    q = mp.expjpi(mp.mpf(2)/X); s1 = s0 = 0; P = mp.mpc(1)
    for n in range(X//2):
        s1 += (-1)**n * q**(n*(n+1)//2) * P; s0 += (-1)**n * q**((n+1)*(n+2)//2) * P
        P *= 1 + q**(n+1)
    return s1, s0
print('psi1*(i) =', mp.nstr(psis(4)[0], 15), ' (claim 2 - i; Fable parts sum: ', mp.nstr(mp.mpc(0.276,-0.146)+mp.mpc(-0.539,-0.092)+mp.mpc(2.263,-0.763), 4), ')')
rows = []
for X in (200, 400, 800, 1600, 3200):
    with mp.workdps(30 + X//10):
        Cl2 = mp.clsin(2, mp.pi/3)
        L = mp.mpf(3)**(-0.25)*mp.sqrt(mp.mpf(X)/2)*mp.exp(Cl2*X/(2*mp.pi))*mp.expjpi(mp.mpf(2)/3 - mp.mpf(X)/6)
        s1, s0 = psis(X)
        r1 = X*(s1/L - 1); r0 = X*(s0/(mp.expjpi(mp.mpf(2)/3)*L) - 1)
        rows.append((X, r1, r0))
        print(f'X = {X:5d}:  psi1*: X(ratio-1) = {mp.nstr(r1, 8):>26}   psi0*: X(ratio-1) = {mp.nstr(r0, 8):>26}', flush=True)
for k, name in ((1, 'psi1*'), (2, 'psi0*')):
    v = [r[k] for r in rows]; xs = [mp.mpf(1)/r[0] for r in rows]
    for j in range(1, len(v)):
        for i in range(len(v) - j): v[i] = (xs[i+j]*v[i] - xs[i]*v[i+1])/(xs[i+j] - xs[i])
    print(f'{name}: extrapolated X(ratio-1) -> {mp.nstr(v[0], 10)}   (compare 0.2015 -/+ i pi/2, i pi/6:  pi/2 = {mp.nstr(mp.pi/2, 8)}, pi/6 = {mp.nstr(mp.pi/6, 8)})')
