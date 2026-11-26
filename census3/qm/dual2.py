"""q -> 1/q partner search for psi1*, psi0* in the 1,357 bounded census functions.
Fits  T(zeta) = a zeta^e X(zeta') + c zeta^e U(zeta') + b  (a, b, c real; zeta' = conj zeta or zeta;
U = the other psi*), e in j/120, separately on all roots, odd-order roots and even-order roots."""
import json, cmath, sys, math
import numpy as np
from math import gcd
ROOTS = [(0, 1)] + [(h, m) for m in range(2, 13) for h in range(1, m) if gcd(h, m) == 1]
L = json.load(open('limits.json'))
def T(nm, h, m): return complex(*L[f"{nm}|{h}/{m}"][:2])
LIB = {}
for l in open('limits2.jsonl'):
    r = json.loads(l)
    v = r['lim']
    if any(v.get(f"{h}/{m}") is None for h, m in ROOTS): continue
    if any(not all(map(math.isfinite, x)) or x[2] > 1e-6 * max(1, abs(complex(x[0], x[1]))) for x in v.values()): continue
    LIB[r['id']] = {k: complex(x[0], x[1]) for k, x in v.items()}
print(len(LIB), 'library functions with reliable limits at all', len(ROOTS), 'roots', flush=True)
SUBS = {'all': ROOTS, 'odd': [r for r in ROOTS if r[1] % 2], 'even': [r for r in ROOTS if r[1] % 2 == 0]}
ES = [j / 120 for j in range(-120, 121)]
def fit(tv, cols, roots, e):
    A = []; b = []
    for k, (h, m) in enumerate(roots):
        z = cmath.exp(2j * cmath.pi * e * h / m)
        A.append([z * c[k] for c in cols] + [1]); b.append(tv[k])
    A = np.array(A); b = np.array(b)
    AR = np.vstack([A.real, A.imag]); bR = np.concatenate([b.real, b.imag])
    x = np.linalg.lstsq(AR, bR, rcond=None)[0]
    return np.linalg.norm(AR @ x - bR) / (np.linalg.norm(bR) + 1e-300), x
hits = []
for tgt, other in [('psi1*', 'psi0*'), ('psi0*', 'psi1*')]:
    for sub, roots in SUBS.items():
        tv = [T(tgt, h, m) for h, m in roots]
        for conj in (True, False):
            def at(h, m): return ((m - h) % m, m) if conj else (h, m)
            uv = [T(other, *at(h, m)) for h, m in roots]
            for i, X in LIB.items():
                xv = [X[f"{a}/{b}"] for a, b in (at(h, m) for h, m in roots)]
                best = min((fit(tv, [xv, uv], roots, e) + (e,) for e in ES), key=lambda t: t[0])
                if best[0] < 1e-2:
                    hits.append((best[0], tgt, sub, 'conj' if conj else 'same', i, best[2], np.round(best[1], 6)))
hits.sort(key=lambda t: t[0])
json.dump([[float(h[0]), h[1], h[2], h[3], h[4], h[5], list(map(float, h[6]))] for h in hits], open('dual2_hits.json', 'w'))
print(len(hits), 'fits with residual < 1e-3')
for h in hits[:40]: print(h)
