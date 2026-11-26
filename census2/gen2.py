"""Candidate generator for rank-2 (double-sum) q-series, with exact dedupe and a q->1 pre-filter.

Candidate:
    S(q) = sum_{m,n>=0} s1^m s2^n q^{(a m^2 + 2b mn + c n^2 + d m + e n)/2}
                       * prod_{k<=m} F1(k) * prod_{k<=n} F2(k) * prod_{k<=m+n} F12(k)

F1, F2 come from FACT (a single factor with a power), F12 from COUPLE.
Usage:  python3 gen2.py [--size small|medium|large]   ->  cands2.json
"""
import itertools, json, math, sys
from fractions import Fraction
import numpy as np
from scipy.special import spence
from scipy.optimize import minimize

NCOEF = 40
# name: (coefficient dict of the factor at index k, evaluator with Q(k) = q^k, dilog decomposition)
# dilog decomposition: list of (d, sigma, weight): factor = prod (1 + sigma q^{d k})^weight asymptotically
FACT = {
    "1/(q)":       (lambda k: {0: 1, k: -1},         lambda k, Q: 1 - Q(k),         [(1, -1, 1)], -1),
    "1/(q^2;q^2)": (lambda k: {0: 1, 2 * k: -1},     lambda k, Q: 1 - Q(2 * k),     [(2, -1, 1)], -1),
    "1/(-q)":      (lambda k: {0: 1, k: 1},          lambda k, Q: 1 + Q(k),         [(1, 1, 1)], -1),
    "(-q)":        (lambda k: {0: 1, k: 1},          lambda k, Q: 1 + Q(k),         [(1, 1, 1)], 1),
    "1/(q;q^2)":   (lambda k: {0: 1, 2 * k - 1: -1}, lambda k, Q: 1 - Q(2 * k - 1), [(2, -1, 1)], -1),
    "(-q;q^2)":    (lambda k: {0: 1, 2 * k - 1: 1},  lambda k, Q: 1 + Q(2 * k - 1), [(2, 1, 1)], 1),
}
FACT["(q)"] = (lambda k: {0: 1, k: -1}, lambda k, Q: 1 - Q(k), [(1, -1, 1)], 1)
COUPLE = {"": None, "1/(q)_{m+n}": "1/(q)", "(q)_{m+n}": "(q)", "(-q)_{m+n}": "(-q)"}

SIZES = {
    "small":  dict(AC=(1, 2, 3), B=(-1, 0, 1, 2), F=("1/(q)", "1/(-q)", "(-q)"), CP=("", "1/(q)_{m+n}"), S=((1, 1), (-1, 1))),
    "medium": dict(AC=(1, 2, 3, 4), B=(-2, -1, 0, 1, 2), F=("1/(q)", "1/(q^2;q^2)", "1/(-q)", "(-q)"),
                   CP=("", "1/(q)_{m+n}", "(-q)_{m+n}"), S=((1, 1), (-1, 1), (-1, -1))),
    "large":  dict(AC=(1, 2, 3, 4), B=(-2, -1, 0, 1, 2), F=tuple(k for k in FACT if k != "(q)"),
                   CP=("", "1/(q)_{m+n}", "(q)_{m+n}", "(-q)_{m+n}"), S=((1, 1), (-1, 1), (1, -1), (-1, -1))),
}

def coupling(name):
    return COUPLE[name]

def expo(c, m, n):
    return (c["a"] * m * m + 2 * c["b"] * m * n + c["c"] * n * n + c["d"] * m + c["e"] * n) // 2

def candidates(size):
    P = SIZES[size]
    out = []
    for a, cc in itertools.product(P["AC"], P["AC"]):
        for b in P["B"]:
            if a * cc - b * b <= 0:
                continue
            for d in range(-a, a + 3):
                if (a + d) % 2:
                    continue
                for e in range(-cc, cc + 3):
                    if (cc + e) % 2:
                        continue
                    base = dict(a=a, b=b, c=cc, d=d, e=e)
                    if min(expo(base, m, n) for m in range(12) for n in range(12) if m + n) < 1:
                        continue
                    for (s1, s2) in P["S"]:
                        for f1, f2 in itertools.product(P["F"], P["F"]):
                            for cp in P["CP"]:
                                if b == 0 and cp == "":
                                    continue  # product of two single sums
                                out.append(dict(base, s1=s1, s2=s2, f1=f1, f2=f2, cp=cp))
    return out

def label(c):
    q = f"q^(({c['a']}m²{2 * c['b']:+d}mn{c['c']:+d}n²{c['d']:+d}m{c['e']:+d}n)/2)"
    sg = ("(−1)^m " if c["s1"] < 0 else "") + ("(−1)^n " if c["s2"] < 0 else "")
    cp = f" · {c['cp']}" if c["cp"] else ""
    return f"Σ {sg}{q} · {c['f1']}_m · {c['f2']}_n{cp}"

# ---------- exact coefficients (int64 is ample for 40 coefficients) ----------
from functools import lru_cache

def apply(P, cf, p, N):
    """multiply (p=1) or divide (p=-1) the series P by a factor given by its coefficient dict"""
    terms = [(k, v) for k, v in cf.items() if k > 0]
    P = P.copy()
    if p > 0:
        for i in range(N - 1, -1, -1):
            for k, v in terms:
                if i - k >= 0:
                    P[i] += v * P[i - k]
    else:
        for i in range(N):
            for k, v in terms:
                if i - k >= 0:
                    P[i] -= v * P[i - k]
    return P

@lru_cache(maxsize=None)
def prods(f, K, N):
    """tuple of truncated series prod_{k<=j} factor(k)^power, j = 0..K"""
    one = np.zeros(N, dtype=np.int64); one[0] = 1
    if f is None:
        return tuple([one] * (K + 1))
    cf, _, _, p = FACT[f]
    out, P = [one], one
    for k in range(1, K + 1):
        P = apply(P, cf(k), p, N)
        out.append(P)
    return tuple(out)

def coeffs(c, N=NCOEF):
    # every (m, n) with E(m, n) < N. With a cross term b != 0 the smallest exponents of a row can
    # lie far from the axes, so the bound K is found by scanning whole rows and columns
    # (an earlier version used only the axes and missed terms when b < 0).
    def row_min(K):
        return min(min(expo(c, K, n) for n in range(0, 6 * K + 6)), min(expo(c, m, K) for m in range(0, 6 * K + 6)))
    K = 2
    while row_min(K) < N:
        K += 1
    K = 6 * K + 6
    P1, P2, P12 = prods(c["f1"], K, N), prods(c["f2"], K, N), prods(coupling(c["cp"]), 2 * K, N)
    tot = np.zeros(N, dtype=np.int64)
    for m in range(K + 1):
        for n in range(K + 1):
            E = expo(c, m, n)
            if E >= N:
                continue
            t = np.convolve(np.convolve(P1[m], P2[n])[:N], P12[m + n])[:N - E]
            tot[E:] += ((c["s1"] ** m) * (c["s2"] ** n)) * t
    return [int(x) for x in tot]

def canon(cf):
    """canonical form up to additive constant, monomial factor, overall sign and q -> -q"""
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

# ---------- q -> 1 pre-filter (non-alternating candidates only) ----------
def Li2(x):
    return spence(1 - x)

def phi(c):
    parts = []
    for f, sel in ((c["f1"], (1, 0)), (c["f2"], (0, 1)), (coupling(c["cp"]), (1, 1))):
        if f is None:
            continue
        _, _, dec, p = FACT[f]
        for d, s, w in dec:
            parts.append((d, s, w * p, sel))
    a, b, cc = c["a"], c["b"], c["c"]
    def F(x):
        u, v = abs(x[0]), abs(x[1])
        val = -(a * u * u + 2 * b * u * v + cc * v * v) / 2
        for d, s, w, (i, j) in parts:
            y = i * u + j * v
            val += w / d * (Li2(-s * math.exp(-d * y)) - Li2(-s))
        return val
    return F

def predict_kappa(c):
    F = phi(c)
    best, arg = -1e9, None
    for u in np.linspace(0, 6, 31):
        for v in np.linspace(0, 6, 31):
            y = F((u, v))
            if y > best:
                best, arg = y, (u, v)
    r = minimize(lambda x: -F(x), arg, method="Nelder-Mead", options=dict(xatol=1e-13, fatol=1e-16, maxiter=4000))
    return max(best, -r.fun, 0.0) / (4 * math.pi ** 2)

def rational(k, D=2000, tol=1e-10):
    f = Fraction(k).limit_denominator(D)
    return abs(float(f) - k) < tol, f

if __name__ == "__main__":
    size = sys.argv[sys.argv.index("--size") + 1] if "--size" in sys.argv else "medium"
    C = candidates(size)
    print(len(C), "candidates", flush=True)
    groups = {}
    for i, c in enumerate(C):
        k = canon(coeffs(c))
        if k is not None:
            groups.setdefault(k, []).append(i)
        if i % 20000 == 0:
            print(" coefficients", i, len(groups), flush=True)
    print(len(groups), "distinct functions", flush=True)
    reps, dropped = [], 0
    for k, idx in groups.items():
        best = min(idx, key=lambda i: (C[i]["s1"] < 0 or C[i]["s2"] < 0, C[i]["cp"] != "", C[i]["a"] + C[i]["c"],
                                       abs(C[i]["b"]), abs(C[i]["d"]) + abs(C[i]["e"])))
        c = C[best]
        rec = dict(id=len(reps), n_forms=len(idx), coeffs=list(k)[:15], **c, label=label(c))
        # pre-filter only applies when some representative has all terms positive
        pos = [i for i in idx if C[i]["s1"] > 0 and C[i]["s2"] > 0]
        if pos:
            kap = predict_kappa(C[pos[0]])
            ok, fr = rational(kap)
            rec["kappa_pred"] = kap
            if kap > 1e-7 and not ok:
                dropped += 1
                continue
        reps.append(rec)
    print(len(reps), "kept after the q->1 dilogarithm pre-filter;", dropped, "dropped", flush=True)
    json.dump(reps, open("cands2.json", "w"))
