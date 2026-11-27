# Independent check of addendum-v4 sec. 4 (converse law) and sec. 5 (mirror law), exact values at roots,
# Neville extrapolation in u = Y^{-1/2} along Y = rho + M0 + 120k (resp. Z = rho + M0 + 240k).
import mpmath as mp
from fractions import Fraction as Fr
mp.mp.dps = 40
def expo(fr): return mp.expjpi(2 * mp.mpf(fr.numerator) / fr.denominator)
def psi_star(fr, shift):            # exact value at e^{2 pi i fr}; rule of addendum-v3 sec. 1
    fr = fr % 1; m = fr.denominator
    with mp.workdps(40 + m // 8):
        return +_psi_star(fr, shift, m)
def _psi_star(fr, shift, m):
    z = expo(fr); s = 0; P = mp.mpc(1)
    for n in range(m):
        e = n*(n+1)//2 if shift == 1 else (n+1)*(n+2)//2
        s += (-1)**n * z**e * P; P *= 1 + z**(n+1)
    return s / 3 if m % 2 else s
def Gsum(fr, a):                    # sum (-1)^j xi^{j^2+a j} (xi;xi^2)_j, a = 1 (W) or 3 (W3), geometric completion
    fr = fr % 1; m = fr.denominator
    with mp.workdps(40 + m // 8):
        return +_Gsum(fr, a, m)
def _Gsum(fr, a, m):
    xi = expo(fr)
    L = 2 * m + 4; terms = []; P = mp.mpc(1)
    for j in range(L):
        terms.append((-1)**j * xi**(j*j + a*j) * P); P *= 1 - xi**(2*j+1)
    if m % 2: return sum(terms[:m+1])            # odd order: terminates
    s = m // 2                                    # even order 2s: block s+j = c * block j
    c = terms[s] / terms[0]
    assert all(abs(terms[j+s] - c*terms[j]) < mp.mpf(10)**-25 * (1 + abs(terms[j])) for j in range(s))
    return sum(terms[:s]) / (1 - c)
phi = (1 + mp.sqrt(5)) / 2
def neville(xs, ys):
    P = list(ys)
    for k in range(1, len(P)):
        for i in range(len(P) - k):
            P[i] = (xs[i+k] * P[i] - xs[i] * P[i+1]) / (xs[i+k] - xs[i])
    return P[0]
c = mp.sqrt(2 * mp.pi / mp.sqrt(5))
print('Converse law (sec. 4): limit of W(e^{2pi i/Y}) sqrt(tau) / prediction')
for rho, M0 in [(Fr(0), 3), (Fr(1, 3), 7), (Fr(2, 5), 1), (Fr(1, 2), 4)]:
    for which in ('W', 'W3'):
        xs, ys = [], []
        for k in (1, 2, 4, 8, 16):
            Y = rho + M0 + 120 * k
            tau = -2j * mp.pi * Y.denominator / Y.numerator
            q1 = -Y / 2
            p0, p1 = psi_star(q1, 0), psi_star(q1, 1)
            e1 = mp.expjpi(mp.mpf(Y.numerator) / (60 * Y.denominator)); e2 = mp.expjpi(-11 * mp.mpf(Y.numerator) / (60 * Y.denominator))
            if which == 'W':
                lhs = Gsum(1 / Y, 1); pred = c * (-e1 * p0 + 1j * e2 * p1)
            else:
                lhs = Gsum(1 / Y, 3); pred = c * (-phi * e1 * p0 - 1j / phi * e2 * p1)
            xs.append(1 / mp.sqrt(mp.mpf(Y.numerator) / Y.denominator)); ys.append(lhs * mp.sqrt(tau) - pred)
        d = neville(xs, ys)
        print(f"   rho={str(rho):4s} {which:3s} |limit - prediction| = {mp.nstr(abs(d), 3)}   (|prediction| ~ {mp.nstr(abs(pred), 3)})")
kap = mp.sqrt(mp.pi / mp.sqrt(5))
print('Mirror law (sec. 5): q = -e^{2 pi i/Z}')
for rho, M0 in [(Fr(0), 5), (Fr(1, 3), 2), (Fr(1, 2), 9)]:
    for which in (1, 0):
        xs, ys = [], []
        for k in (1, 2, 4, 8, 16):
            Z = rho + M0 + 240 * k
            t = -2j * mp.pi * Z.denominator / Z.numerator
            qhat = Fr(1, 2) - Z / 4
            p0, p1 = psi_star(qhat, 0), psi_star(qhat, 1)
            e1 = mp.expjpi(mp.mpf(Z.numerator) / (120 * Z.denominator)); e2 = mp.expjpi(-11 * mp.mpf(Z.numerator) / (120 * Z.denominator))
            lhs = psi_star(Fr(1, 2) + 1 / Z, which)
            pred = kap * ((-1/phi * e1 * p0 + 1j * phi * e2 * p1) if which == 1 else (phi * e1 * p0 - 1j/phi * e2 * p1))
            xs.append(1 / mp.sqrt(mp.mpf(Z.numerator) / Z.denominator)); ys.append(lhs * mp.sqrt(t) - pred)
        d = neville(xs, ys)
        print(f"   rho={str(rho):4s} psi{which}* |limit - prediction| = {mp.nstr(abs(d), 3)}   (|prediction| ~ {mp.nstr(abs(pred), 3)})")
