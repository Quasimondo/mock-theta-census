# Spot check of addendum-v8 sec. 2: F0*(q) ~ A(t) J(r) (1 + (1/8 + i/(9 sqrt3)) t), X = X0 + 4m - i b,
# q = e^{2 pi i/X}, t = -2 pi i/X, r = e^{-pi i X/2}, A(t) = sqrt(pi/(2 sqrt3)) e^{-i pi/12} t^{-1/2} exp(-(pi^2/24 + i Cl2/2)/t)
import mpmath as mp
def F0s(q, dps):
    s = 0; P = mp.mpc(1); n = 0
    while True:
        term = (-1)**n * q**(2*n*n) / P
        s += term
        if n > 10 and abs(term) < mp.mpf(10)**(-dps + 15) * abs(s): return s
        P *= 1 - q**(2*n+1); n += 1
def J(r):
    s = 0; P = mp.mpc(1)
    for m in range(80):
        s += r**(2*m*m+m) / P; P *= 1 + r**(2*m+2)
    return s
corr = mp.mpf(1)/8 + 1j/(9*mp.sqrt(3))
for X0, b in [(0, 1), (1, mp.mpf('0.6'))]:
    for m in (50, 100, 200, 400):
        dps = int(60 + (4*m)/3)
        with mp.workdps(dps):
            X = mp.mpc(X0 + 4*m, -b); q = mp.exp(2j*mp.pi/X); t = -2j*mp.pi/X; r = mp.exp(-1j*mp.pi*X/2)
            Cl2 = mp.clsin(2, mp.pi/3)
            A = mp.sqrt(mp.pi/(2*mp.sqrt(3))) * mp.expjpi(-mp.mpf(1)/12) * t**mp.mpf(-0.5) * mp.exp(-(mp.pi**2/24 + 1j*Cl2/2)/t)
            ratio = F0s(q, dps) / (A * J(r) * (1 + corr*t))
            print(f'X0={X0} b={mp.nstr(b,2)} |r|={mp.nstr(abs(r),3)} |X|={mp.nstr(abs(X),5):>6}: |ratio-1| = {mp.nstr(abs(ratio-1),3):>9}   (|t|^2 = {mp.nstr(abs(t)**2,2)})', flush=True)
