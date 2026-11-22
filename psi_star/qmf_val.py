import numpy as np
from qmf import psi1s, psi0s, reduce, points
E = lambda al, t: np.exp(2j * np.pi * al * t)
FN = {'psi1*': psi1s, 'psi0*': psi0s}
def data(cls, lo_den, hi_den):
    rec = []
    for (h, m) in points(0.05, 2.0, hi_den, even_only=True):
        if m % 4 != cls or m < lo_den: continue
        p2, q2 = reduce(h, 2 * h + m)
        rec.append((h / m, p2 / q2, 2 * h / m + 1, {n: complex(f(h, m)) for n, f in FN.items()}, {n: complex(f(p2 % q2, q2)) for n, f in FN.items()}))
    rec.sort(key=lambda r: r[0]); return rec
def mats(rec, p, target):
    a1, a0, k = p; al = {'psi1*': a1, 'psi0*': a0}
    X = np.array([r[0] for r in rec]); Y = np.array([r[1] for r in rec]); CX = np.array([r[2] for r in rec])
    dphi = np.diff(E(al[target], X) * np.array([r[3][target] for r in rec]))
    A = np.stack([np.diff(CX ** (-k) * E(al[n], Y) * np.array([r[4][n] for r in rec])) for n in FN], 1)
    return dphi, A
train = data(0, 2, 64); test2 = data(2, 2, 64); hiden = data(0, 66, 90)
for target, p in [('psi1*', (0.179493, -0.218409, 0.507778)), ('psi1*', (20/120, 9/120, 0.5)), ('psi0*', (0.388115, -0.012957, 0.507942))]:
    d, A = mats(train, p, target); sol, *_ = np.linalg.lstsq(A, d, rcond=None)
    out = []
    for nm, rec in [('train m=0 mod 4, m<=64', train), ('test m=2 mod 4', test2), ('test m in 66..90', hiden)]:
        d2, A2 = mats(rec, p, target); r = d2 - A2 @ sol
        out.append(f"{nm}: {np.vdot(r,r).real/np.vdot(d2,d2).real:.2g}")
    print(target, p, ' | '.join(out))
