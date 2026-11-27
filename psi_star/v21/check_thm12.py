# Independent check of addendum-v21 (own code).
# (1) Theorem 12: radial limit A/3 (A = one period of the series, odd X; v3 odd-root rule) against
#     (1/3) 3^{-1/4} (X/2)^{1/2} e^{Cl2 X/2pi} e^{2pi i/3 - i pi X/6}  (psi0*: times e^{2pi i/3}); X(ratio-1), extrapolated.
# (2) Lemma R2': W(p) -> 1/3, p W3(p) -> -1/3 on the path x = 1/X + i eps, X = 5 (W summed directly).
import mpmath as mp
def periods(X):
    q = mp.expjpi(mp.mpf(2)/X); s1 = s0 = 0; P = mp.mpc(1)
    for n in range(X):
        s1 += (-1)**n * q**(n*(n+1)//2) * P; s0 += (-1)**n * q**((n+1)*(n+2)//2) * P
        P *= 1 + q**(n+1)
    return s1/3, s0/3
rows = []
for X in (201, 401, 801, 1601, 3201):
    with mp.workdps(30 + X//8):
        Cl2 = mp.clsin(2, mp.pi/3)
        L = mp.mpf(3)**(-0.25)/3*mp.sqrt(mp.mpf(X)/2)*mp.exp(Cl2*X/(2*mp.pi))*mp.expjpi(mp.mpf(2)/3 - mp.mpf(X)/6)
        a1, a0 = periods(X)
        rows.append((X, X*(a1/L - 1), X*(a0/(mp.expjpi(mp.mpf(2)/3)*L) - 1)))
        print(f'X = {X:5d}:  psi1*: X(ratio-1) = {mp.nstr(rows[-1][1], 8):>26}   psi0*: {mp.nstr(rows[-1][2], 8):>26}', flush=True)
for k, name in ((1, 'psi1*'), (2, 'psi0*')):
    v = [r[k] for r in rows]; xs = [mp.mpf(1)/r[0] for r in rows]
    for j in range(1, len(v)):
        for i in range(len(v) - j): v[i] = (xs[i+j]*v[i] - xs[i]*v[i+1])/(xs[i+j] - xs[i])
    print(f'{name}: extrapolated -> {mp.nstr(v[0], 10)}   (-2 pi i s1 = {mp.nstr(-2j*mp.pi*(mp.mpf(1)/4 + 1j*mp.sqrt(3)/54), 10)})')
print('(2) Lemma R2\', X = 5')
X = 5
for eps in (mp.mpf('1e-3'), mp.mpf('3e-4'), mp.mpf('1e-4')):
    with mp.workdps(150):
        x = mp.mpf(1)/X + 1j*eps; p = mp.expjpi(-1/x)
        W = W3 = 0; po = mp.mpc(1); j = 0
        while abs(p)**(j*j) > mp.mpf(10)**(-60) or j < 5:
            W += (-1)**j * p**(j*j + j) * po; W3 += (-1)**j * p**(j*j + 3*j) * po; po *= 1 - p**(2*j + 1); j += 1
        print(f'  eps = {mp.nstr(eps, 2)}: 1-|p| = {mp.nstr(1-abs(p), 3)}   |W - 1/3| = {mp.nstr(abs(W - mp.mpf(1)/3), 3)}   |p W3 + 1/3| = {mp.nstr(abs(p*W3 + mp.mpf(1)/3), 3)}')
