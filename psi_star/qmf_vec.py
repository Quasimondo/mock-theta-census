import sys, numpy as np, mpmath as mp
from qmf import psi1s, psi0s, kontsevich, reduce, points
FN = {'psi1*': psi1s, 'psi0*': psi0s}
gam = (1, 0, 2, 1); a, b, c, d = gam
pts = points(0.05, 2.0, 40, even_only=True)
X, Y, CX, V, W = [], [], [], {n: [] for n in FN}, {n: [] for n in FN}
for (h, m) in pts:
    p2, q2 = reduce(a * h + b * m, c * h + d * m)
    if q2 % 2: continue
    X.append(h / m); Y.append(p2 / q2); CX.append(c * h / m + d)
    for n, f in FN.items():
        V[n].append(complex(f(h, m))); W[n].append(complex(f(p2 % q2, q2)))
o = np.argsort(X); X = np.array(X)[o]; Y = np.array(Y)[o]; CX = np.array(CX)[o]
for n in FN: V[n] = np.array(V[n])[o]; W[n] = np.array(W[n])[o]
print(len(X), "points")
D = 120
E = lambda al, t: np.exp(2j * np.pi * al * t)
for target in FN:
    best = []
    for k in [0.5, 1.0, 1.5, 2.0]:
        fac = CX ** (-k)
        for j1 in range(D):            # alpha of target component
            phi = E(j1 / D, X) * V[target]; dphi = np.diff(phi); nrm = np.vdot(dphi, dphi).real
            for j0 in range(D):        # alpha of partner component (target uses j1 for its own transformed value)
                cols = []
                for n in FN:
                    al = j1 / D if n == target else j0 / D
                    cols.append(np.diff(fac * E(al, Y) * W[n]))
                A = np.stack(cols, 1)
                sol, *_ = np.linalg.lstsq(A, dphi, rcond=None)
                r = dphi - A @ sol; R = np.vdot(r, r).real / nrm
                best.append((R, k, j1, j0, sol))
    best.sort(key=lambda z: z[0])
    print("== target", target)
    for R, k, j1, j0, sol in best[:5]:
        print(f"  roughness={R:.3g} weight={k} alpha_self={j1}/{D} alpha_partner={j0}/{D} coeffs={np.round(sol,5)} |c|={np.round(np.abs(sol),5)}")
