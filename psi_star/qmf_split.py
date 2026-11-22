import numpy as np, mpmath as mp
from qmf import psi1s, psi0s, kontsevich, reduce, points
FN = {'psi1*': psi1s, 'psi0*': psi0s}
a, b, c, d = 1, 0, 2, 1
pts = points(0.05, 2.0, 60, even_only=True)
rec = []
for (h, m) in pts:
    p2, q2 = reduce(a * h + b * m, c * h + d * m)
    rec.append(dict(x=h / m, y=p2 / q2, cx=c * h / m + d, cls=m % 4,
                    v={n: complex(f(h, m)) for n, f in FN.items()}, w={n: complex(f(p2 % q2, q2)) for n, f in FN.items()}))
E = lambda al, t: np.exp(2j * np.pi * al * t)
D = 120
for cls in (0, 2):
    R0 = sorted([r for r in rec if r['cls'] == cls], key=lambda r: r['x'])
    X = np.array([r['x'] for r in R0]); Y = np.array([r['y'] for r in R0]); CX = np.array([r['cx'] for r in R0])
    V = {n: np.array([r['v'][n] for r in R0]) for n in FN}; W = {n: np.array([r['w'][n] for r in R0]) for n in FN}
    print(f"== denominators m = {cls} mod 4: {len(X)} points")
    for target in FN:
        best = []
        for k in [0.5, 1.0, 1.5]:
            fac = CX ** (-k)
            for j1 in range(D):
                phi = E(j1 / D, X) * V[target]; dphi = np.diff(phi); nrm = np.vdot(dphi, dphi).real
                # single-component and two-component fits
                for j0 in range(D):
                    al = {target: j1 / D}; other = [n for n in FN if n != target][0]; al[other] = j0 / D
                    A1 = np.diff(fac * E(al[target], Y) * W[target])[:, None]
                    A2 = np.stack([np.diff(fac * E(al[n], Y) * W[n]) for n in FN], 1)
                    for tag, A in (('1-comp', A1), ('2-comp', A2)):
                        if tag == '1-comp' and j0: continue
                        sol, *_ = np.linalg.lstsq(A, dphi, rcond=None); r = dphi - A @ sol
                        best.append((np.vdot(r, r).real / nrm, tag, k, j1, j0, sol))
        best.sort(key=lambda z: z[0])
        for R, tag, k, j1, j0, sol in best[:3]:
            print(f"   {target}: roughness={R:.3g} [{tag}] weight={k} alpha={j1}/{D} partner_alpha={j0}/{D} coeffs={np.round(sol,6)} |c|={np.round(np.abs(sol),6)}")
