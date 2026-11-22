import numpy as np
from scipy.optimize import minimize
from qmf import psi1s, psi0s, reduce, points
E = lambda al, t: np.exp(2j * np.pi * al * t)
FN = {'psi1*': psi1s, 'psi0*': psi0s}
rec = []
for (h, m) in points(0.05, 2.0, 64, even_only=True):
    if m % 4: continue
    p2, q2 = reduce(h, 2 * h + m)
    rec.append((h / m, p2 / q2, 2 * h / m + 1, {n: complex(f(h, m)) for n, f in FN.items()}, {n: complex(f(p2 % q2, q2)) for n, f in FN.items()}))
rec.sort(key=lambda r: r[0])
X = np.array([r[0] for r in rec]); Y = np.array([r[1] for r in rec]); CX = np.array([r[2] for r in rec])
V = {n: np.array([r[3][n] for r in rec]) for n in FN}; W = {n: np.array([r[4][n] for r in rec]) for n in FN}
def R(p, target):
    a1, a0, k = p
    al = {'psi1*': a1, 'psi0*': a0}
    dphi = np.diff(E(al[target], X) * V[target])
    A = np.stack([np.diff(CX ** (-k) * E(al[n], Y) * W[n]) for n in FN], 1)
    sol, *_ = np.linalg.lstsq(A, dphi, rcond=None); r = dphi - A @ sol
    return np.log10(np.vdot(r, r).real / np.vdot(dphi, dphi).real), sol
for target, start in [('psi1*', (20 / 120, 9 / 120, 0.5)), ('psi0*', (11 / 120, 0.0, 0.5))]:
    best = None
    for s in [start, (start[0] + 0.01, start[1] - 0.01, 0.5), (start[0] - 0.01, start[1] + 0.02, 0.55)]:
        o = minimize(lambda p: R(p, target)[0], s, method='Nelder-Mead', options=dict(xatol=1e-7, fatol=1e-4, maxiter=3000))
        if best is None or o.fun < best.fun: best = o
    lr, sol = R(best.x, target)
    print(f"{target}: best log10 roughness={lr:.2f} at alpha(psi1*)={best.x[0]:.6f} alpha(psi0*)={best.x[1]:.6f} weight={best.x[2]:.6f} coeffs={np.round(sol,5)}")
