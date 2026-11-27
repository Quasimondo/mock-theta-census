# Independent check of addendum-v5 sec. 2 (converse law for W inside the disk, with first corrections):
# W(q) = Phi_W(tau) + tau^{-1/2} c [ -u0 psi0*(q') (1 + a0 tau) + i u1 psi1*(q') (1 + a1 tau) ] + O(tau^{3/2})
# q = e^{2 pi i/Y}, Y = Y0 + 120k - i b, tau = -2 pi i/Y, q' = e^{-pi i Y}, u0 = q'^{-1/60}, u1 = q'^{11/60}
from fractions import Fraction as Fr
import mpmath as mp
K = 30
def ser_exp(a):                     # e^{-a t} to degree K
    s = [Fr(1)]; f = Fr(1)
    for k in range(1, K + 1):
        f = f * (-a) / k; s.append(f)
    return s
def smul(x, y):
    z = [Fr(0)] * (K + 1)
    for i, a in enumerate(x):
        if a:
            for j in range(K + 1 - i): z[i + j] += a * y[j]
    return z
Phi = [Fr(0)] * (K + 1); P = [Fr(1)] + [Fr(0)] * K
for j in range(K + 1):
    term = smul(ser_exp(j * j + j), P)
    for i in range(K + 1): Phi[i] += (-1) ** j * term[i]
    one_minus = [Fr(0)] + [-x for x in ser_exp(2 * j + 1)[1:]]     # 1 - e^{-(2j+1)t}
    P = smul(P, one_minus)
import math
print('k! a_k:', [int(Phi[k] * math.factorial(k)) for k in range(8)], '(addendum-v4: 1, -1, 11, -253, 10859, -745981, 75045131, -10399927693)')
def W_at(q, dps):
    with mp.workdps(dps):
        s = 0; P = mp.mpc(1); j = 0
        while True:
            term = (-1) ** j * q ** (j * j + j) * P
            s += term
            if j > 10 and abs(term) < mp.mpf(10) ** (-dps + 20) * (1 + abs(s)): break
            P *= 1 - q ** (2 * j + 1); j += 1
        return +s
def psi_star(q, shift):
    s = 0; P = mp.mpc(1)
    for n in range(400):
        e = n * (n + 1) // 2 if shift == 1 else (n + 1) * (n + 2) // 2
        s += (-1) ** n * q ** e * P; P *= 1 + q ** (n + 1)
    return s
mp.mp.dps = 50
c = mp.sqrt(2 * mp.pi / mp.sqrt(5)); r5 = mp.sqrt(5)
a0 = (20 + 3 * r5) / 120; a1 = (20 - 3 * r5) / 120
for Y0, b in [(3, 1), (7, mp.mpf('0.6'))]:
    print(f'Y0 = {Y0}, b = {b}')
    for k in (2, 4, 8, 16):
        Y = mp.mpc(Y0 + 120 * k, -b)
        q = mp.exp(2j * mp.pi / Y); tau = -2j * mp.pi / Y; qp = mp.exp(-1j * mp.pi * Y)
        dps = int(60 + abs(Y) / 3)
        w1 = W_at(q, dps); w2 = W_at(q, dps + 60)
        # optimal truncation of Phi_W
        terms = [mp.mpf(Phi[i].numerator) / Phi[i].denominator * tau ** i for i in range(K + 1)]
        n_opt = min(range(K + 1), key=lambda i: abs(terms[i]))
        phi_w = sum(terms[:n_opt])
        u0 = mp.exp(1j * mp.pi * Y / 60); u1 = mp.exp(-11j * mp.pi * Y / 60)
        lead = tau ** mp.mpf(-0.5) * c * (-u0 * psi_star(qp, 0) * (1 + a0 * tau) + 1j * u1 * psi_star(qp, 1) * (1 + a1 * tau))
        rel = abs(w1 - phi_w - lead) / abs(lead)
        print(f'   |Y| = {mp.nstr(abs(Y), 5):>6}  relative residual {mp.nstr(rel, 3):>9}   (|tau|^2 = {mp.nstr(abs(tau)**2, 2)}; precision check {mp.nstr(abs(w1 - w2), 2)}; Phi truncated at {n_opt})', flush=True)
