# Compare psi1* with the Kashaev invariant of the figure-eight knot 4_1.
# J_N = sum_{n<N} |(q;q)_n|^2, q = e^{2 pi i/N};  known: J_N ~ 3^{-1/4} N^{3/2} e^{Vol N/(2 pi)} (1 + O(1/N)), Vol = 2 Cl2(pi/3).
import mpmath as mp
from fractions import Fraction as Fr
mp.mp.dps = 30
def kashaev(N):
    with mp.workdps(30 + N // 10):
        q = mp.expjpi(mp.mpf(2) / N); s = 0; P = mp.mpc(1)
        for n in range(N):
            s += abs(P)**2; P *= 1 - q**(n+1)
        return +s
def psi1s(X):          # exact value at e^{2 pi i/X}, X integer (addendum-v3 sec. 1)
    with mp.workdps(30 + X // 10):
        z = mp.expjpi(mp.mpf(2) / X); s = 0; P = mp.mpc(1)
        for n in range(X):
            s += (-1)**n * z**(n*(n+1)//2) * P; P *= 1 + z**(n+1)
        return +(s / 3 if X % 2 else s)
Cl2 = mp.clsin(2, mp.pi/3); Vol = 2*Cl2
print('Vol(4_1) = 2 Cl2(pi/3) =', mp.nstr(Vol, 12), ';  psi* grows like exp(Cl2 X/2pi) = exp(Vol X/4pi): half the Kashaev rate exp(Vol N/2pi)')
print('\n1) Kashaev control: J_N / (3^{-1/4} N^{3/2} e^{Vol N/2pi})  -> 1')
for N in (100, 200, 400, 800):
    print(f'   N={N:4d}: {mp.nstr(kashaev(N) / (3**-0.25 * mp.mpf(N)**1.5 * mp.exp(Vol*N/(2*mp.pi))), 10)}')
print('\n2) Same exponential: compare |psi1*(e^{2 pi i/X})| with J_{X/2}  (X even, so W = 1)')
for X in (200, 400, 800, 1600):
    r = abs(psi1s(X)) / kashaev(X // 2)
    print(f'   X={X:5d}: |psi1*| / J_(X/2) = {mp.nstr(r, 8):>12}   * X = {mp.nstr(r*X, 8)}')
# predicted from the laws: |psi1*| ~ |C| sqrt(X/2pi) e^{Vol X/4pi};  J_{X/2} ~ 3^{-1/4} (X/2)^{3/2} e^{Vol X/4pi}
C = mp.sqrt(mp.pi/mp.sqrt(3))
print('   predicted limit of the ratio * X:', mp.nstr(C / mp.sqrt(2*mp.pi) / (3**-0.25 * mp.mpf(2)**-1.5), 8))
print('\n3) Perturbative series in the same variable h = 2 pi i/N with N = X/2, i.e. t = -h/2:')
print('   4_1 (Zagier):   Phi(h)/Phi(0) = 1 + (11/(72 sqrt(-3))) h + ...  -> coefficient of t:', mp.nstr(-11/(72*mp.sqrt(-3))*2, 8))
print('   psi1* (v3/v5):  1 + (1/4 + i/(18 sqrt3)) t + ...                 -> coefficient of t:', mp.nstr(mp.mpf(1)/4 + 1j/(18*mp.sqrt(3)), 8))
