"""Radial limits at roots of unity of the bounded functions, by polynomial
extrapolation in t of g(zeta e^{-t}), t = T0*j, j = 1..K.
    python limits.py [mmax]      -> limits.json"""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../census'))
os.chdir(os.path.join(os.path.dirname(__file__), '../../census'))
from math import gcd
import mpmath as mp
from flint import acb, arb, ctx
import catalog, twist
ctx.prec = 200
mp.mp.dps = 50
catalog.add('sigma', 'x', '', catalog.eul(lambda n: n * (n + 1) // 2, lambda n, Q: 1 / (1 + Q(n))))
catalog.BY = {e['name']: e for e in catalog.E}
NAMES = ['psi1*', 'psi0*', 'F0*', 'F1*', 'T0 (8th)*', 'psi*', 'sigma']
MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 13
ROOTS = [(0, 1)] + [(h, m) for m in range(2, MMAX + 1) for h in range(1, m) if gcd(h, m) == 1]
K = 14
def lim(fn, h, m):
    ts = [0.003 / m ** 2 * j for j in range(1, K + 1)]
    vals = []
    for t in ts:
        q = (acb(2 * h) / m).exp_pi_i() * arb(-t).exp()
        v = fn(q)
        vals.append(mp.mpc(float(v.real.mid()), float(v.imag.mid())) if False else mp.mpc(mp.mpf(v.real.mid().str(40, radius=False)), mp.mpf(v.imag.mid().str(40, radius=False))))
    # Neville extrapolation to t = 0
    P = list(vals)
    for k in range(1, K):
        for i in range(K - k):
            P[i] = (ts[i + k] * P[i] - ts[i] * P[i + 1]) / (ts[i + k] - ts[i])
    est = P[0]
    # error estimate: same with K-1 points
    P2 = list(vals[:-1])
    for k in range(1, K - 1):
        for i in range(K - 1 - k):
            P2[i] = (ts[i + k] * P2[i] - ts[i] * P2[i + 1]) / (ts[i + k] - ts[i])
    return est, abs(est - P2[0])
out = {}
if os.path.exists('../census3/qm/limits.json'):
    out = json.load(open('../census3/qm/limits.json'))
for nm in NAMES:
    fn = catalog.BY[nm]['fn']
    for (h, m) in ROOTS:
        key = f"{nm}|{h}/{m}"
        if key in out: continue
        t0 = time.time()
        v, err = lim(fn, h, m)
        out[key] = [float(v.real), float(v.imag), float(err)]
        print(key, mp.nstr(v, 12), 'err', mp.nstr(err, 2), f"{time.time()-t0:.1f}s", flush=True)
        json.dump(out, open('../census3/qm/limits.json', 'w'))
