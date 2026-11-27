# Spot check of addendum-v6 sec. 2: B1 = N0(q) E+'(p) + i q^{-1/4} Y(q) E+(p) behaves like
# c e^{i Cl2/t} (1 + (-1/24 - i/(18 sqrt3)) t + ...) with no dependence on p.
# X = X0 + 2m - i b, q = e^{2 pi i/X}, t = -2 pi i/X, p = e^{-pi i X}
import mpmath as mp
def sum_q(f, q, dps):
    with mp.workdps(dps):
        s = 0; n = 0; P = mp.mpc(1)     # P = (q;q)_{2n} or (q;q)_{2n+1}
        return f(q)
def N0_Y(q, dps):
    with mp.workdps(dps):
        n0 = 0; y = 0; P = mp.mpc(1); n = 0     # P = (q;q)_{2n}
        while True:
            t0 = (-1)**n * q**(n*n+n) / P
            P1 = P * (1 - q**(2*n+1))
            t1 = (-1)**n * q**((n+1)**2) / P1
            n0 += t0; y += t1
            if n > 10 and abs(t0) + abs(t1) < mp.mpf(10)**(-dps+15) * (abs(n0) + abs(y)): break
            P = P1 * (1 - q**(2*n+2)); n += 1
        return +n0, +y
def E(p, sign, a):
    s = 0; P = mp.mpc(1)
    for m in range(200):
        s += sign**m * p**(mp.mpf(m*(m+a))/2) / P; P *= 1 - p**(2*m+2)
    return s
mp.mp.dps = 40
Cl2 = mp.clsin(2, mp.pi/3)
c = mp.mpf(4/3)**0.25 * mp.expjpi(-mp.mpf(1)/12)
corr = -mp.mpf(1)/24 - 1j/(18*mp.sqrt(3))
for X0, b in [(0, 1), (1, mp.mpf('0.7')), (mp.mpf('0.5'), mp.mpf('1.3'))]:
    for m in (100, 200, 400, 800):
        dps = int(60 + abs(mp.mpc(X0 + 2*m, -b))/2)
        with mp.workdps(dps):
            X = mp.mpc(X0 + 2*m, -b); q = mp.exp(2j*mp.pi/X); t = -2j*mp.pi/X
            p = mp.exp(-1j*mp.pi*X)
            n0, y = N0_Y(q, dps)
            B1 = n0 * E(p, 1, 1) + 1j * q**(-mp.mpf(1)/4) * y * E(p, -1, 1)
            r = B1 / (c * mp.exp(1j*Cl2/t) * (1 + corr*t))
        print(f'X0={mp.nstr(X0,2):>4} b={mp.nstr(b,2):>4} |p|={mp.nstr(abs(p),3):>8}  |X|={mp.nstr(abs(X),4):>5}: |ratio - 1| = {mp.nstr(abs(r-1),3):>9}  (|t|^2 = {mp.nstr(abs(t)**2,2)})', flush=True)
