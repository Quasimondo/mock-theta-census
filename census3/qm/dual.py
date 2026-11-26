"""Look for q -> 1/q (zeta -> conj zeta) relations among radial limits:
target(zeta) = zeta^e * sum_i a_i X_i(conj zeta) + sum_j b_j Y_j(zeta), least squares over all roots."""
import json, cmath, itertools
import numpy as np
from math import gcd
L = json.load(open('limits.json'))
def sig(z):
    s, p = 1, 1
    for n in range(400):
        s += (-1) ** n * z ** (n + 1) * p; p *= (1 - z ** (n + 1))
    return s
def sigs(z):
    s, p = 0, 1
    for n in range(400):
        s += z ** (n + 1) * p; p *= (1 - z ** (2 * n + 2))
    return -2 * s
roots = [(0, 1)] + [(h, m) for m in range(2, 16) for h in range(1, m) if gcd(h, m) == 1]
def val(nm, h, m):
    z = cmath.exp(2j * cmath.pi * h / m)
    if nm == 'sigma': return sig(z)
    if nm == 'sigma*': return sigs(z)
    if nm == '1': return 1
    return complex(*L[f"{nm}|{h}/{m}"][:2])
LIB = ['psi1*', 'psi0*', 'sigma', 'sigma*', '1']
import sys
for target in (sys.argv[1:] or ['psi1*', 'psi0*']):
    best = []
    for e in [j / 120 for j in range(-120, 121)]:
        rows, rhs = [], []
        for (h, m) in roots:
            z = cmath.exp(2j * cmath.pi * h / m)
            hb = (m - h) % m
            ze = cmath.exp(2j * cmath.pi * e * h / m)
            row = [ze * val(x, hb, m) for x in LIB] + [val(x, h, m) for x in LIB if x not in (target, '1')]
            rows.append(row); rhs.append(val(target, h, m))
        A = np.array(rows); b = np.array(rhs)
        # real coefficients: stack real and imaginary parts
        AR = np.vstack([A.real, A.imag]); bR = np.concatenate([b.real, b.imag])
        x, res, rk, sv = np.linalg.lstsq(AR, bR, rcond=None)
        r = np.linalg.norm(AR @ x - bR) / np.linalg.norm(bR)
        best.append((r, e, np.round(x, 4)))
    best.sort(key=lambda t: t[0])
    print(target)
    for r, e, x in best[:3]:
        print(f"   e = {e:+.4f}  rel. residual {r:.2e}  coeffs {x}")
