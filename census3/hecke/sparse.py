"""Step 1: multiply the bounded census functions by eta quotients and look for sparse
(Hecke-type) coefficient patterns, the way f0(q)(q;q)_inf is sparse."""
import itertools, sys
N = 400
def mul(a, b):
    c = [0] * N
    for i, x in enumerate(a):
        if x:
            for j in range(N - i):
                if b[j]: c[i + j] += x * b[j]
    return c
def poch_series(factors):  # product of (1 + s q^k)^{+-1}
    pass
def qpoch(a, step, sign=-1, inv=False, n=None):
    """(sign*q^a; q^step)_n as series (n=None: infinite)."""
    s = [0] * N; s[0] = 1
    k = 0; e = a
    while e < N and (n is None or k < n):
        if inv:  # divide by (1 + sign q^e)
            for i in range(e, N): s[i] -= sign * s[i - e]
        else:
            for i in range(N - 1, e - 1, -1): s[i] += sign * s[i - e]
        k += 1; e += step
    return s
def series_sum(term):
    tot = [0] * N; n = 0
    while True:
        t = term(n)
        if t is None: break
        tot = [x + y for x, y in zip(tot, t)]; n += 1
    return tot
def shift(s, k): return ([0] * k + s)[:N] if k >= 0 else s[-k:] + [0] * (-k)
def psi1s(n):
    e = n * (n + 1) // 2
    if e >= N: return None
    return [(-1) ** n * x for x in shift(qpoch(1, 1, +1, n=n), e)]
def psi0s(n):
    e = (n + 1) * (n + 2) // 2
    if e >= N: return None
    return [(-1) ** n * x for x in shift(qpoch(1, 1, +1, n=n), e)]
def psi3s(n):  # third-order psi twisted: sum_{n>=1} (-1)^n q^{n^2}/(q;q^2)_n
    if n == 0: return [0] * N
    if n * n >= N: return None
    return [(-1) ** n * x for x in shift(qpoch(1, 2, -1, inv=True, n=n), n * n)]
def F0s(n):
    if 2 * n * n >= N: return None
    return [(-1) ** n * x for x in shift(qpoch(1, 2, -1, inv=True, n=n), 2 * n * n)]
def F1s(n):
    if 2 * n * n + 2 * n >= N: return None
    return [(-1) ** n * x for x in shift(qpoch(1, 2, -1, inv=True, n=n + 1), 2 * n * n + 2 * n)]
def T0s(n):
    e = (n + 1) * (n + 2)
    if e >= N: return None
    a = qpoch(2, 2, +1, n=n); b = qpoch(1, 2, +1, inv=True, n=n + 1)
    return [(-1) ** n * x for x in shift(mul(a, b), e)]
T = {'psi1*': series_sum(psi1s), 'psi0*': series_sum(psi0s), 'psi3*': series_sum(psi3s),
     'F0*': series_sum(F0s), 'F1*': series_sum(F1s), 'T0*': series_sum(T0s)}
E = {a: qpoch(a, a) for a in (1, 2, 3, 4, 6)}
EI = {a: qpoch(a, a, inv=True) for a in (1, 2, 3, 4, 6)}
def etaq(ex):
    s = [1] + [0] * (N - 1)
    for a, e in ex.items():
        for _ in range(abs(e)): s = mul(s, E[a] if e > 0 else EI[a])
    return s
res = []
bases = (1, 2, 4)
for ex in itertools.product(range(-2, 3), repeat=3):
    d = dict(zip(bases, ex))
    P = etaq(d)
    for nm, t in T.items():
        c = mul(t, P)
        nz = sum(1 for x in c if x)
        mx = max(abs(x) for x in c)
        res.append((nz, mx, nm, ex))
for nm in T:
    print(nm, 'first coeffs', T[nm][:16])
    best = sorted([r for r in res if r[2] == nm])[:4]
    for nz, mx, _, ex in best:
        print(f"   x (q)^{ex[0]} (q2)^{ex[1]} (q4)^{ex[2]}: {nz}/{N} nonzero, max |c| {mx}")
