"""Candidate generator for Ramanujan-style (Eulerian) q-series and exact coefficient dedupe.

Candidate:  sum_{n>=0} s^n q^{(A n^2 + B n)/2} prod_{k=1}^{n} prod_j F_j(k)^{p_j}
"""
import itertools, json, sys
from flint import acb, arb, ctx
import catalog

NCOEF = 40
FACT = {  # name: (pretty, coefficient poly of the factor at index k, evaluator with Q(k)=q^k)
    "(1+q^k)":       (lambda k: {0: 1, k: 1},            lambda n, Q: 1 + Q(n)),
    "(1-q^k)":       (lambda k: {0: 1, k: -1},           lambda n, Q: 1 - Q(n)),
    "(1+q^2k)":      (lambda k: {0: 1, 2 * k: 1},        lambda n, Q: 1 + Q(2 * n)),
    "(1-q^2k)":      (lambda k: {0: 1, 2 * k: -1},       lambda n, Q: 1 - Q(2 * n)),
    "(1+q^(2k-1))":  (lambda k: {0: 1, 2 * k - 1: 1},    lambda n, Q: 1 + Q(2 * n - 1)),
    "(1-q^(2k-1))":  (lambda k: {0: 1, 2 * k - 1: -1},   lambda n, Q: 1 - Q(2 * n - 1)),
    "(1-q^k+q^2k)":  (lambda k: {0: 1, k: -1, 2 * k: 1}, lambda n, Q: 1 - Q(n) + Q(2 * n)),
    "(1+q^k+q^2k)":  (lambda k: {0: 1, k: 1, 2 * k: 1},  lambda n, Q: 1 + Q(n) + Q(2 * n)),
}
POWS = [-2, -1, 1]
AB = [(A, B) for A in (1, 2, 3, 4) for B in range(-A, A + 3) if (A + B) % 2 == 0]

def candidates():
    names = list(FACT)
    pats = [((a, p),) for a in names for p in POWS]
    pats += [((a, p), (b, r)) for a, b in itertools.combinations(names, 2) for p in POWS for r in POWS]
    out = []
    for (A, B) in AB:
        for s in (1, -1):
            for pat in pats:
                out.append(dict(A=A, B=B, s=s, pat=pat))
    return out

def label(c):
    e = f"q^{{({c['A']}n²{'+' if c['B'] >= 0 else '−'}{abs(c['B'])}n)/2}}" if c['B'] else f"q^{{{c['A']}n²/2}}"
    fac = " ".join(f"∏{f}" + ("" if p == 1 else f"^{p}") for f, p in c["pat"])
    return ("(−1)ⁿ " if c["s"] < 0 else "") + e + " · " + fac

# ---- exact truncated power series ----
import numpy as np
def mul(a, b, N=NCOEF):
    return list(np.convolve(np.asarray(a, dtype=np.int64), np.asarray(b, dtype=np.int64))[:N])

def poly(d, N=NCOEF):
    r = [0] * N
    for k, v in d.items():
        if k < N:
            r[k] += v
    return r

def inv(a, N=NCOEF):
    r = [0] * N; r[0] = 1
    for i in range(1, N):
        r[i] = -sum(a[j] * r[i - j] for j in range(1, i + 1))
    return r

def coeffs(c, N=NCOEF):
    A, B, s = c["A"], c["B"], c["s"]
    total = [0] * N
    P = [1] + [0] * (N - 1)
    n = 0
    while True:
        e = n * (A * n + B) // 2
        if n > 1 and e >= N:
            break
        if n > 0:
            for f, p in c["pat"]:
                fp = poly(FACT[f][0](n), N)
                g = fp if p > 0 else inv(fp, N)
                for _ in range(abs(p)):
                    P = mul(P, g, N)
        if e < N:
            sg = s ** n
            for i in range(N - e):
                total[e + i] += sg * P[i]
        n += 1
    return total

def canon(cf):
    """canonical form up to additive constant, monomial factor q^v, overall sign, and q -> -q"""
    a = [int(x) for x in cf[1:]]
    nz = [i for i, x in enumerate(a) if x]
    if not nz:
        return None
    v = nz[0]
    a = a[v:v + 25]
    if a[0] < 0:
        a = [-x for x in a]
    b = [x * (-1) ** i for i, x in enumerate(a)]
    if b[0] < 0:
        b = [-x for x in b]
    return tuple(min(a, b))

def numeric_coeffs(fn, N=NCOEF, M=128):
    """coefficients of an analytic function from samples on |q| = 1/2 (rounded to integers)"""
    ctx.prec = 300
    r = arb(1) / 2
    vals = [fn(r * (acb(2 * j) / M).exp_pi_i()) for j in range(M)]
    out = []
    for k in range(N):
        s = acb(0)
        for j in range(M):
            s += vals[j] * (acb(-2 * j * k) / M).exp_pi_i()
        a = (s / M) / r ** k
        out.append(round(float(a.real.mid())))
    return out

def make_fn(c):
    A, B, s, pat = c["A"], c["B"], c["s"], c["pat"]
    def ratio(n, Q):
        v = acb(s)
        for f, p in pat:
            x = FACT[f][1](n, Q)
            v = v * (x ** p if p > 0 else (1 / x) ** (-p))
        return v
    return catalog.eul(lambda n: n * (A * n + B) // 2, ratio)

if __name__ == "__main__":
    C = candidates()
    print(len(C), "candidates")
    known = {}
    for e in catalog.E:
        cf = numeric_coeffs(e["fn"])
        k = canon(cf)
        if k:
            known[k] = e["name"]
    groups = {}
    for i, c in enumerate(C):
        k = canon(coeffs(c))
        if k is None:
            continue
        groups.setdefault(k, []).append(i)
    print(len(groups), "distinct functions (up to constant, q-shift, sign, q->-q)")
    reps = []
    for k, idx in groups.items():
        # simplest representative: fewest factors, smallest |powers|, smallest A, |B|
        best = min(idx, key=lambda i: (len(C[i]["pat"]), sum(abs(p) for _, p in C[i]["pat"]), C[i]["A"], abs(C[i]["B"]), C[i]["s"] < 0))
        reps.append(dict(id=best, n_forms=len(idx), known=known.get(k), coeffs=list(k)[:15], **C[best], label=label(C[best])))
    print(sum(1 for r in reps if r["known"]), "match catalog functions:", sorted({r["known"] for r in reps if r["known"]}))
    json.dump(reps, open("cands.json", "w"))
