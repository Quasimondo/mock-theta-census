import numpy as np, mpmath as mp
from scipy.optimize import minimize
from qmf import reduce, points, TWOPI_I
def mk(A, B, sgn):   # sum sgn^n zeta^{(A n^2 + B n)/2} (-zeta;zeta)_n   (terminates at even-order roots)
    def f(h, m):
        z = mp.e ** (TWOPI_I * mp.mpf(h) / m); s = 0; P = 1
        for n in range(m):
            s += sgn ** n * z ** ((A * n * n + B * n) // 2) * P; P *= (1 + z ** (n + 1))
            if abs(P) < mp.mpf(10) ** -30: break
        return s
    return f
E = lambda al, t: np.exp(2j * np.pi * al * t)
def test(FN, target, starts):
    rec = []
    for (h, m) in points(0.05, 2.0, 64, even_only=True):
        if m % 4: continue
        p2, q2 = reduce(h, 2 * h + m)
        rec.append((h / m, p2 / q2, 2 * h / m + 1, {n: complex(f(h, m)) for n, f in FN.items()}, {n: complex(f(p2 % q2, q2)) for n, f in FN.items()}))
    rec.sort(key=lambda r: r[0])
    X = np.array([r[0] for r in rec]); Y = np.array([r[1] for r in rec]); CX = np.array([r[2] for r in rec])
    V = {n: np.array([r[3][n] for r in rec]) for n in FN}; W = {n: np.array([r[4][n] for r in rec]) for n in FN}
    names = list(FN)
    def R(p):
        al = {names[0]: p[0], names[1]: p[1]}
        d = np.diff(E(al[target], X) * V[target])
        A = np.stack([np.diff(CX ** (-p[2]) * E(al[n], Y) * W[n]) for n in FN], 1)
        sol, *_ = np.linalg.lstsq(A, d, rcond=None); r = d - A @ sol
        return np.log10(np.vdot(r, r).real / np.vdot(d, d).real)
    best = min((minimize(R, s, method='Nelder-Mead', options=dict(xatol=1e-7, fatol=1e-4, maxiter=3000)) for s in starts), key=lambda o: o.fun)
    return best.fun, best.x
starts = [(0.15, 0.05, 0.5), (0.3, -0.1, 0.6), (0.05, 0.2, 1.0), (0.4, 0.4, 1.5)]
pairs = {
  'unrelated A: (-1)^n q^{n^2}(-q)_n & (-1)^n q^{n^2+n}(-q)_n': {'a': mk(2, 0, -1), 'b': mk(2, 2, -1)},
  'unrelated B: q^{n(n+1)}(-q)_n & (-1)^n q^{(3n^2+n)/2}(-q)_n': {'a': mk(2, 2, 1), 'b': mk(3, 1, -1)},
  'psi1* & psi0* (for reference)': {'a': mk(1, 1, -1), 'b': mk(1, 3, -1)},
}
for nm, FN in pairs.items():
    lr, x = test(FN, 'a', starts)
    print(f"{nm}:  best log10 roughness = {lr:.2f}  at alpha=({x[0]:.4f}, {x[1]:.4f}) weight={x[2]:.4f}")
