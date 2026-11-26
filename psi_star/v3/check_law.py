# Independent check of addendum-v3 sec. 5.3:
# psi1*(e^{2 pi i/X}) ~ Cbar t^{-1/2} e^{-Sbar/t} W(e^{-pi i X}) (1 + s1bar t),  t = -2 pi i / X
import mpmath as mp
from fractions import Fraction as Fr
mp.mp.dps = 60
def psi_star(fr, shift=0):          # exact value at zeta = e^{2 pi i fr}, rule of sec. 1 at odd order
    m = fr.denominator; z = mp.expjpi(2 * mp.mpf(fr.numerator) / m)
    s = 0; P = mp.mpc(1)
    for n in range(m):
        e = n*(n+1)//2 if shift == 0 else (n+1)*(n+2)//2
        s += (-1)**n * z**e * P
        P *= 1 + z**(n+1)
    if m % 2 == 1: return s / 3
    return s                        # even order: the series terminates within m terms
def W(fr, L=400):                   # W(xi), xi = e^{2 pi i fr}: sum (-1)^j xi^{j^2+j} (xi;xi^2)_j, geometric completion
    xi = mp.expjpi(2 * mp.mpf(fr.numerator) / fr.denominator)
    terms = []; P = mp.mpc(1)
    for j in range(L):
        terms.append((-1)**j * xi**(j*j+j) * P)
        P *= 1 - xi**(2*j+1)
    if abs(terms[-1]) < mp.mpf(10)**-40 * max(1, abs(sum(terms))):
        return sum(terms)
    for p in range(1, L//3):        # find period p with term_{j+p} = c term_j
        c = terms[p] / terms[0] if terms[0] != 0 else None
        if c is not None and all(abs(terms[j+p] - c*terms[j]) < mp.mpf(10)**-30 * (1+abs(terms[j])) for j in range(L - p)):
            return sum(terms[:p]) / (1 - c)
    raise ValueError('no period')
Cl = mp.clsin(2, mp.pi/3)
Sb = mp.pi**2/3 + 1j*Cl
Cb = mp.sqrt(mp.pi/mp.sqrt(3)) * mp.expjpi(mp.mpf(5)/12)
s1b = mp.mpf(1)/4 + 1j*mp.sqrt(3)/54
print('W(1), W(-1), W(-i), W(e^{-i pi/3}):', [mp.nstr(W(Fr(a)), 8) for a in (0, Fr(1, 2), Fr(3, 4), Fr(5, 6))])
for rho in [Fr(0), Fr(1), Fr(1, 2), Fr(1, 3), Fr(3, 2), Fr(2, 5)]:
    rows = []
    for N in (200, 400, 800):
        X = N + rho
        t = -2j * mp.pi * X.denominator / mp.mpf(X.numerator)
        lead = Cb * t**mp.mpf(-0.5) * mp.exp(-Sb / t)
        Wv = W((-X / 2) % 1)
        ratio = psi_star(Fr(1) / X) / (lead * Wv)
        rows.append((N, ratio, (ratio - 1) / t))
    print(f"rho={str(rho):4s} W={mp.nstr(W((-(Fr(400)+rho)/2) % 1), 6):>22}  ratio-1 at X~800: {mp.nstr(abs(rows[-1][1]-1), 3):>8}  (ratio-1)/t -> {mp.nstr(rows[-1][2], 8)}  (s1bar = {mp.nstr(s1b, 8)})")
# Richardson in t for rho = 2/5: (ratio-1)/t = s1 + s2' t + ..., with X = N + rho
print('Richardson, rho = 2/5:')
vals = []
for N in (400, 800, 1600, 3200):
    X = N + Fr(2, 5); t = -2j * mp.pi * X.denominator / mp.mpf(X.numerator)
    r = psi_star(Fr(1) / X) / (Cb * t**mp.mpf(-0.5) * mp.exp(-Sb / t) * W((-X / 2) % 1))
    vals.append((t, (r - 1) / t))
ts = [v[0] for v in vals]; P = [v[1] for v in vals]
for k in range(1, len(P)):
    for i in range(len(P) - k):
        P[i] = (ts[i + k] * P[i] - ts[i] * P[i + 1]) / (ts[i + k] - ts[i])
print('  s1 estimate', mp.nstr(P[0], 12), ' claimed', mp.nstr(s1b, 12), ' diff', mp.nstr(abs(P[0] - s1b), 3))
