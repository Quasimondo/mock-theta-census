"""Q2: perturbative series of psi1*, psi0* at q = -1, q = i (terminating sums, expanded exactly in t)
and at q -> 1 (through f1*/2, f0*/2, which differ from psi* by an exponentially small term; see Q3).
Output: q2_series.json with coefficients a_k of t^k as decimal strings (real, imag).
Usage: python3 q2_series.py K     (default K = 150)
"""
import sys, json
from flint import acb_series, acb, arb, ctx
K = int(sys.argv[1]) if len(sys.argv) > 1 else 150
DIG = 120
def expser(k):                      # exp(-k t) as a series
    return acb_series([0, -k]).exp()
def at_even_root(h, m, shift, prec):
    """psi*(zeta e^(-t)), zeta = e^(2 pi i h/m), m even: sum (-1)^n q^(e(n)) (-q;q)_n,
    e(n) = n(n+1)/2 (shift 0, psi1*) or (n+1)(n+2)/2 (shift 1, psi0*)."""
    ctx.prec = prec; ctx.cap = K + 1
    zeta = (acb(2*h)/m).exp_pi_i()
    total = acb_series([0]); P = acb_series([1]); n = 0
    zero_run = 0
    while True:
        if n:
            P = P * (1 + zeta**n * expser(n))
        e = n*(n+1)//2 if shift == 0 else (n+1)*(n+2)//2
        term = P * (zeta**e * expser(e))
        total = total + (term if n % 2 == 0 else -term)
        # stop when P has no terms below t^(K+1): all coefficients are (numerically) zero
        if n > m and all(abs(P[j]).mid() < arb(10)**(-200) for j in range(0, K+1, max(1, K//12))) and \
           all(abs(P[j]).mid() < arb(10)**(-200) for j in range(K+1)):
            break
        n += 1
    return total, n
def at_one(shift, prec, nmax):
    """(1/2) f1*(e^(-t)) (shift 0) or (1/2) f0*(e^(-t)) (shift 1) as a series: sum (-1)^n q^(n^2+n or n^2)/(-q;q)_n"""
    ctx.prec = prec; ctx.cap = K + 1
    total = acb_series([0]); P = acb_series([1])
    for n in range(nmax):
        if n:
            P = P / (1 + expser(n))
        e = n*n + n if shift == 0 else n*n
        term = P * expser(e)
        total = total + (term if n % 2 == 0 else -term)
    return total / 2
def dump(s):
    out = []
    for k in range(K+1):
        c = s[k]
        rel = max(float(c.real.rad()), float(c.imag.rad())) / max(float(abs(c).mid()), 1e-300) if abs(c).mid() > 0 else 0.0
        out.append([c.real.mid().str(DIG, radius=False), c.imag.mid().str(DIG, radius=False), rel])
    return out
res = {"K": K}
for name, (h, m) in (("-1", (1, 2)), ("i", (1, 4))):
    for shift, fn in ((0, "psi1"), (1, "psi0")):
        s, n = at_even_root(h, m, shift, 4000)
        res[f"{fn}@{name}"] = dump(s)
        print(f"{fn}* at q = {name}: {n} terms of the sum needed; a_0..a_3 =", [s[k].real.mid().str(8, radius=False) + "+" + s[k].imag.mid().str(8, radius=False) + "i" for k in range(4)],
              " worst relative radius:", max(x[2] for x in res[f"{fn}@{name}"]), flush=True)
for shift, fn in ((0, "psi1"), (1, "psi0")):
    s = at_one(shift, 2400, 40*K//2 + 400)
    res[f"{fn}@1"] = dump(s)
    print(f"{fn}* at q -> 1 (as f*/2): a_0..a_4 =", [s[k].real.mid().str(12, radius=False) for k in range(5)], " worst relative radius:", max(x[2] for x in res[f"{fn}@1"]), flush=True)
json.dump(res, open("q2_series.json", "w"))
