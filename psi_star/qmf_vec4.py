import numpy as np, mpmath as mp
from qmf import psi1s, psi0s, reduce, points, TWOPI_I
def psi1(h, m):
    z = mp.e ** (TWOPI_I * mp.mpf(h) / m); s = 0; P = 1
    for n in range(m):
        s += z ** (n * (n + 1) // 2) * P; P *= (1 + z ** (n + 1))
        if abs(P) < mp.mpf(10) ** -30: break
    return s
def psi0(h, m):
    z = mp.e ** (TWOPI_I * mp.mpf(h) / m); s = 0; P = 1
    for n in range(m):
        s += z ** ((n + 1) * (n + 2) // 2) * P; P *= (1 + z ** (n + 1))
        if abs(P) < mp.mpf(10) ** -30: break
    return s
FN = {'psi1*': psi1s, 'psi0*': psi0s, 'psi1': psi1, 'psi0': psi0}
a, b, c, d = 1, 0, 2, 1
pts = points(0.05, 2.0, 40, even_only=True)
X, Y, CX = [], [], []; V = {n: [] for n in FN}; W = {n: [] for n in FN}
for (h, m) in pts:
    p2, q2 = reduce(a * h + b * m, c * h + d * m)
    if q2 % 2: continue
    X.append(h / m); Y.append(p2 / q2); CX.append(c * h / m + d)
    for n, f in FN.items():
        V[n].append(complex(f(h, m))); W[n].append(complex(f(p2 % q2, q2)))
o = np.argsort(X); X = np.array(X)[o]; Y = np.array(Y)[o]; CX = np.array(CX)[o]
for n in FN: V[n] = np.array(V[n])[o]; W[n] = np.array(W[n])[o]
E = lambda al, t: np.exp(2j * np.pi * al * t)
D = 120; k = 0.5; fac = CX ** (-k)
# all four components share one alpha per *pair* (psi1 with psi1*, psi0 with psi0*): the sign flip shouldn't change the q-power
for target in ['psi1*', 'psi0*']:
    best = []
    for j1 in range(D):
        for j0 in range(D):
            al = {'psi1*': j1 / D, 'psi1': j1 / D, 'psi0*': j0 / D, 'psi0': j0 / D}
            phi = E(al[target], X) * V[target]; dphi = np.diff(phi); nrm = np.vdot(dphi, dphi).real
            A = np.stack([np.diff(fac * E(al[n], Y) * W[n]) for n in FN], 1)
            sol, *_ = np.linalg.lstsq(A, dphi, rcond=None)
            r = dphi - A @ sol; best.append((np.vdot(r, r).real / nrm, j1, j0, sol))
    best.sort(key=lambda z: z[0])
    print("== target", target, "(components psi1*, psi0*, psi1, psi0; weight 1/2)")
    for R, j1, j0, sol in best[:4]:
        print(f"  roughness={R:.3g}  alpha(psi1 pair)={j1}/{D}  alpha(psi0 pair)={j0}/{D}  coeffs={np.round(sol,5)}  |c|={np.round(np.abs(sol),4)}")
