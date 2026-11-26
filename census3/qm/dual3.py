"""Like dual2.py, but uses every library function on the roots where its limit is reliable
(relative extrapolation error < 1e-6), if there are at least MINR such roots."""
import json, cmath, math
import numpy as np
from math import gcd
MINR = 24
ROOTS = [(0, 1)] + [(h, m) for m in range(2, 13) for h in range(1, m) if gcd(h, m) == 1]
L = json.load(open('limits.json'))
def T(nm, h, m): return complex(*L[f"{nm}|{h}/{m}"][:2])
def ok(x): return x is not None and all(map(math.isfinite, x)) and x[2] <= 1e-6 * max(1, abs(complex(x[0], x[1])))
LIB = {}
for l in open('limits2.jsonl'):
    r = json.loads(l); v = r['lim']
    good = {k: complex(x[0], x[1]) for k, x in v.items() if ok(x)}
    LIB[r['id']] = good
ES = [j / 120 for j in range(-120, 121)]
n_used = 0; hits = []
for i, X in LIB.items():
    for conj in (True, False):
        roots = [(h, m) for (h, m) in ROOTS if f"{(m - h) % m if conj else h}/{m}" in X]
        if len(roots) < MINR: continue
        n_used += 1
        for tgt, other in [('psi1*', 'psi0*'), ('psi0*', 'psi1*')]:
            at = (lambda h, m: ((m - h) % m, m)) if conj else (lambda h, m: (h, m))
            tv = np.array([T(tgt, h, m) for h, m in roots])
            xv = np.array([X[f"{a}/{b}"] for a, b in (at(h, m) for h, m in roots)])
            uv = np.array([T(other, *at(h, m)) for h, m in roots])
            ph = np.array([2 * np.pi * h / m for h, m in roots])
            best = (9, 0, None)
            for e in ES:
                z = np.exp(1j * e * ph)
                A = np.stack([z * xv, z * uv, np.ones_like(z)], 1)
                AR = np.vstack([A.real, A.imag]); bR = np.concatenate([tv.real, tv.imag])
                x = np.linalg.lstsq(AR, bR, rcond=None)[0]
                r = np.linalg.norm(AR @ x - bR) / np.linalg.norm(bR)
                if r < best[0]: best = (r, e, x)
            if best[0] < 1e-3:
                hits.append((best[0], tgt, 'conj' if conj else 'same', i, len(roots), best[1], np.round(best[2], 6)))
print(n_used, '(function, orientation) pairs with at least', MINR, 'reliable roots')
for h in sorted(hits, key=lambda t: t[0]): print(h)
