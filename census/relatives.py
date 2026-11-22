import json, collections
from fractions import Fraction as Fr
import search
N = 60
reps = {r['id']: r for r in json.load(open('cands.json'))}
label = {int(k): v for k, v in json.load(open('labels.json')).items()}
def co(i):
    c = dict(reps[i]); c['pat'] = [tuple(x) for x in c['pat']]
    return [int(x) for x in search.coeffs(c, N=N)]
Q = [i for i, v in label.items() if v == 'quantum']
C = {i: co(i) for i in Q}
def variants(a):
    # q -> -q and q-shifts up to 3 (shift right = multiply by q^k)
    out = []
    for sg in (1, -1):
        b = [x * (sg ** i) for i, x in enumerate(a)]
        for k in range(0, 4):
            out.append((sg, k, [0] * k + b[:N - k]))
    return out
def solve(target, basis):
    # exact rational least-solve: target = sum c_j basis_j + c0 ; use Gaussian elimination on all rows
    cols = basis + [[1] + [0] * (N - 1)]
    m = len(cols)
    rows = [[Fr(cols[j][i]) for j in range(m)] + [Fr(target[i])] for i in range(N)]
    # gaussian elimination
    piv = []; r = 0
    for c in range(m):
        p = next((i for i in range(r, N) if rows[i][c] != 0), None)
        if p is None: continue
        rows[r], rows[p] = rows[p], rows[r]
        pv = rows[r][c]; rows[r] = [x / pv for x in rows[r]]
        for i in range(N):
            if i != r and rows[i][c] != 0:
                f = rows[i][c]; rows[i] = [x - f * y for x, y in zip(rows[i], rows[r])]
        piv.append(c); r += 1
    if any(rows[i][m] != 0 for i in range(r, N)): return None
    sol = [Fr(0)] * m
    for i, c in enumerate(piv): sol[c] = rows[i][m]
    return sol
targets = {'psi1*': 830, 'psi0*/q': 1382}
found = collections.defaultdict(list)
for tn, ti in targets.items():
    T = C[ti]
    for i in Q:
        if i == ti: continue
        for sg, k, v in variants(C[i]):
            s = solve(v, [T])
            if s is not None and s[0] != 0:
                found[tn].append((i, sg, k, s)); break
for tn in targets:
    print('==', tn, len(found[tn]), 'relatives of the form g(±q)·q^k = a·T + b')
    for i, sg, k, s in found[tn][:20]:
        print(f"   c{i:5d}  q->{'q' if sg>0 else '-q'} shift {k}:  a={s[0]}, b={s[1]}   {reps[i]['label']}")
# two-term: relatives that need both psi1* and psi0*
both = []
for i in Q:
    if i in targets.values() or any(i == f[0] for l in found.values() for f in l): continue
    for sg, k, v in variants(C[i]):
        for kk in range(0, 3):
            s = solve(v, [C[830], [0] * kk + C[1382][:N - kk]])
            if s is not None and s[0] != 0 and s[1] != 0:
                both.append((i, sg, k, kk, s)); break
        else:
            continue
        break
print('== needing both psi1* and q^kk psi0*/q:', len(both))
for b in both[:15]: print('  ', b[0], b[1], b[2], 'kk', b[3], [str(x) for x in b[4]], reps[b[0]]['label'])
json.dump(dict(found={k: [[f[0], f[1], f[2], [str(x) for x in f[3]]] for f in v] for k, v in found.items()},
               both=[[b[0], b[1], b[2], b[3], [str(x) for x in b[4]]] for b in both]), open('relatives.json', 'w'))
