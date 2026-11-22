import numpy as np, mpmath as mp
from qmf import psi1s, psi0s, kontsevich, reduce, points
E = lambda al, t: np.exp(2j * np.pi * al * t)
def run(FN, target, alphas, k, maxden, cls=0, fitcoef=True):
    a, b, c, d = 1, 0, 2, 1
    rec = []
    for (h, m) in points(0.05, 2.0, maxden, even_only=True):
        if m % 4 != cls: continue
        p2, q2 = reduce(a * h + b * m, c * h + d * m)
        rec.append((h / m, p2 / q2, 2 * h / m + 1, {n: complex(f(h, m)) for n, f in FN.items()}, {n: complex(f(p2 % q2, q2)) for n, f in FN.items()}))
    rec.sort(key=lambda r: r[0])
    X = np.array([r[0] for r in rec]); Y = np.array([r[1] for r in rec]); CX = np.array([r[2] for r in rec])
    fac = CX ** (-k)
    phi = E(alphas[target], X) * np.array([r[3][target] for r in rec]); dphi = np.diff(phi)
    A = np.stack([np.diff(fac * E(alphas[n], Y) * np.array([r[4][n] for r in rec])) for n in FN], 1)
    sol, *_ = np.linalg.lstsq(A, dphi, rcond=None); res = dphi - A @ sol
    return len(X), np.vdot(res, res).real / np.vdot(dphi, dphi).real, sol
print("Kontsevich control (weight 3/2, alpha 1/24):")
for md in [24, 36, 48, 64]:
    n, R, s = run({'K': kontsevich}, 'K', {'K': 1 / 24}, 1.5, md)
    print(f"   maxden={md:3d}  points={n:4d}  roughness={R:.3g}")
print("psi1* with psi0* (weight 1/2, alpha 20/120 and 9/120):")
for md in [24, 36, 48, 64, 80]:
    n, R, s = run({'psi1*': psi1s, 'psi0*': psi0s}, 'psi1*', {'psi1*': 20 / 120, 'psi0*': 9 / 120}, 0.5, md)
    print(f"   maxden={md:3d}  points={n:4d}  roughness={R:.3g}  coeffs={np.round(s,5)}")
print("wrong-parameter reference (psi1* alone, weight 1/2):")
for md in [24, 48, 80]:
    n, R, s = run({'psi1*': psi1s}, 'psi1*', {'psi1*': 39 / 240}, 0.5, md)
    print(f"   maxden={md:3d}  points={n:4d}  roughness={R:.3g}")
