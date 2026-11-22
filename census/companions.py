"""Coefficient-only companion search: find Eulerian series g with
   g(±q) * q^k  =  a*psi1* + b*psi0* + c*psi1 + d*psi0 + e      (exactly, to q^59)."""
import itertools, json, sys, numpy as np
from fractions import Fraction as Fr
from multiprocessing import Pool
import search
N = 60
TYPES = {  # factor at index k as {exponent: coeff}
 "(1+q^k)": lambda k: {0: 1, k: 1}, "(1-q^k)": lambda k: {0: 1, k: -1},
 "(1+q^2k)": lambda k: {0: 1, 2*k: 1}, "(1-q^2k)": lambda k: {0: 1, 2*k: -1},
 "(1+q^(2k-1))": lambda k: {0: 1, 2*k-1: 1}, "(1-q^(2k-1))": lambda k: {0: 1, 2*k-1: -1},
 "(1+q^3k)": lambda k: {0: 1, 3*k: 1}, "(1-q^3k)": lambda k: {0: 1, 3*k: -1},
 "(1+q^(3k-1))": lambda k: {0: 1, 3*k-1: 1}, "(1-q^(3k-1))": lambda k: {0: 1, 3*k-1: -1},
 "(1-q^k+q^2k)": lambda k: {0: 1, k: -1, 2*k: 1}, "(1+q^k+q^2k)": lambda k: {0: 1, k: 1, 2*k: 1},
}
POWS = [-2, -1, 1, 2]
def poly(d):
    r = np.zeros(N, dtype=object)
    for k, v in d.items():
        if k < N: r[k] += v
    return r
def mul(a, b):
    return np.convolve(a, b)[:N]
def inv(a):
    r = np.zeros(N, dtype=object); r[0] = 1
    for i in range(1, N):
        r[i] = -np.dot(a[1:i+1], r[i-1::-1][:i])
    return r
def coeffs(A, B, s, pat):
    tot = np.zeros(N, dtype=object); P = poly({0: 1})
    for n in range(0, 400):
        e = n * (A * n + B) // 2
        if n > 0:
            for f, p in pat:
                fp = poly(TYPES[f](n)); g = fp if p > 0 else inv(fp)
                for _ in range(abs(p)): P = mul(P, g)
        if e < N:
            tot[e:] += (s ** n) * P[:N - e]
        if (A > 0 and e >= N and n > 2) or (A == 0 and n >= N + 2): break
    return tot
def label(A, B, s, pat):
    return ("(−1)ⁿ " if s < 0 else "") + f"q^((%dn²%+dn)/2)" % (A, B) + " " + " ".join(f + ("" if p == 1 else f"^{p}") for f, p in pat)

# targets (from the first template's representatives)
psi1s = [int(x) for x in coeffs(1, 1, -1, [("(1+q^k)", 1)])]
psi0s = [int(x) for x in coeffs(1, 3, -1, [("(1+q^k)", 1)])]   # psi0*/q
psi1 = [int(x) for x in coeffs(1, 1, 1, [("(1+q^k)", 1)])]
psi0 = [int(x) for x in coeffs(1, 3, 1, [("(1+q^k)", 1)])]     # psi0/q
one = [1] + [0] * (N - 1)
BASIS = np.array([psi1s, psi0s, psi1, psi0, one], dtype=float).T
def exact(v):
    cols = [psi1s, psi0s, psi1, psi0, one]; m = 5
    rows = [[Fr(cols[j][i]) for j in range(m)] + [Fr(v[i])] for i in range(N)]
    piv, r = [], 0
    for c in range(m):
        p = next((i for i in range(r, N) if rows[i][c] != 0), None)
        if p is None: continue
        rows[r], rows[p] = rows[p], rows[r]; pv = rows[r][c]; rows[r] = [x / pv for x in rows[r]]
        for i in range(N):
            if i != r and rows[i][c] != 0:
                f = rows[i][c]; rows[i] = [x - f * y for x, y in zip(rows[i], rows[r])]
        piv.append(c); r += 1
    if any(rows[i][m] != 0 for i in range(r, N)): return None
    sol = [Fr(0)] * m
    for i, c in enumerate(piv): sol[c] = rows[i][m]
    return sol
def work(job):
    A, B, s, pat = job
    try:
        a = coeffs(A, B, s, pat)
    except Exception:
        return None
    a = [int(x) for x in a]
    if max(abs(x) for x in a) > 10 ** 14: return None
    for sg in (1, -1):
        b = [x * sg ** i for i, x in enumerate(a)]
        for k in range(0, 3):
            v = np.array([0] * k + b[:N - k], dtype=float)
            x, *_ = np.linalg.lstsq(BASIS, v, rcond=None)
            if np.abs(BASIS @ x - v).max() < 1e-6 * max(1, np.abs(v).max()):
                sol = exact([0] * k + b[:N - k])
                if sol and (sol[0] != 0 or sol[1] != 0):
                    return dict(f=label(A, B, s, pat), A=A, B=B, s=s, pat=pat, qsign=sg, shift=k, coef=[str(c) for c in sol])
    return None
if __name__ == "__main__":
    names = list(TYPES)
    pats = [[(a, p)] for a in names for p in POWS] + [[(a, p), (b, r)] for a, b in itertools.combinations(names, 2) for p in POWS for r in POWS]
    AB = [(A, B) for A in range(0, 7) for B in range(-A if A else 1, A + 5) if (A + B) % 2 == 0 and not (A == 0 and B <= 0)]
    jobs = [(A, B, s, pat) for A, B in AB for s in (1, -1) for pat in pats]
    print(len(jobs), "formulas", flush=True)
    hits = []
    with Pool(2) as p:
        for i, r in enumerate(p.imap_unordered(work, jobs, chunksize=64)):
            if r:
                hits.append(r); print("HIT", r["f"], "q->", r["qsign"], "shift", r["shift"], r["coef"], flush=True)
            if i % 10000 == 0: print("..", i, flush=True)
    json.dump(hits, open("companions.json", "w"), indent=1)
    print("done", len(hits), "hits")
