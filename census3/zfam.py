"""Census 3a: one-parameter families at roots of unity.

    F(z)    = sum_n z^n q^{n(n+1)/2} (-q;q)_n            (psi_1 = F(1), psi_1* = F(-1))
    g2(x,q) = sum_n (-q;q)_n q^{n(n+1)/2} / ((x;q)_{n+1} (q/x;q)_{n+1})
    g3(x,q) = sum_n q^{n(n+1)}            / ((x;q)_{n+1} (q/x;q)_{n+1})

g2, g3 are the universal mock theta functions of Gordon-McIntosh / Hickerson.
z, x run over e^{2 pi i k/N}, N <= NMAX, gcd(k,N)=1 (x = 1 is a pole and is left out).
Each function is probed at the 19 roots h/m with h <= m/2. Coefficients are complex for
non-real x, so the cells of x at h/m and of conj(x) at h/m together give x at all roots.

    python zfam.py --workers 2 [--nmax 12] [--limit 300]
"""
import sys, os, json, math, time, signal
from fractions import Fraction
from math import gcd
from multiprocessing import Pool
from flint import acb
import catalog, probe

NMAX = int(sys.argv[sys.argv.index('--nmax') + 1]) if '--nmax' in sys.argv else 12
LIMIT = int(sys.argv[sys.argv.index('--limit') + 1]) if '--limit' in sys.argv else 300
OUT = 'zfam_results.jsonl'

def unit(k, N):
    return lambda: (acb(2 * k) / N).exp_pi_i() if N > 1 else acb(1)

def F(zf):
    return catalog.eul(lambda n: n * (n + 1) // 2, lambda n, Q: zf() * (1 + Q(n)))

def g(xf, two):
    def ratio(n, Q):
        x = xf()
        r = 1 / ((1 - x * Q(n)) * (1 - Q(n + 1) / x))
        return r * (1 + Q(n)) if two else r
    def init(Q):
        x = xf()
        return 1 / ((1 - x) * (1 - Q(1) / x))
    e = (lambda n: n * (n + 1) // 2) if two else (lambda n: n * (n + 1))
    return catalog.eul(e, ratio, init=init)

FRACS = sorted({Fraction(k, N) for N in range(1, NMAX + 1) for k in range(N) if gcd(k, N) == 1})
NEW = []
for a in FRACS:
    k, N = a.numerator, a.denominator
    tag = f"{k}/{N}"
    catalog.add(f"F[{tag}]", 'F', '', F(unit(k, N))); NEW.append(f"F[{tag}]")
    if a != 0:
        catalog.add(f"g2[{tag}]", 'g2', '', g(unit(k, N), True)); NEW.append(f"g2[{tag}]")
        catalog.add(f"g3[{tag}]", 'g3', '', g(unit(k, N), False)); NEW.append(f"g3[{tag}]")
catalog.BY = {e['name']: e for e in catalog.E}

class Timeout(BaseException):
    pass

def _alarm(*a):
    raise Timeout()

def init_worker():
    try:
        os.nice(19)
    except OSError:
        pass
    signal.signal(signal.SIGALRM, _alarm)

def work(job):
    signal.alarm(LIMIT)
    t0 = time.time()
    try:
        r = probe.probe(job)
    except Timeout:
        r = dict(name=job[0], h=job[1], m=job[2], error='timeout')
    finally:
        signal.alarm(0)
    r['secs'] = round(time.time() - t0, 2)
    return r

def ct(r):
    if 'error' in r:
        return '?'
    if r['resid'] < 1e-25 and r['kappa_err'] < 1e-20:
        return 'E'
    if abs(r['kappa']) < 1e-4:
        return 'Q'
    if r['kappa_err'] < 1e-8 and Fraction(r['kappa_rat']).denominator <= 500:
        return 'R'
    return 'X'

if __name__ == '__main__':
    workers = int(sys.argv[sys.argv.index('--workers') + 1]) if '--workers' in sys.argv else 2
    done = set()
    if os.path.exists(OUT):
        for l in open(OUT):
            d = json.loads(l)
            if 'error' not in d or '--retry' not in sys.argv:
                done.add((d['name'], d['h'], d['m']))
    jobs = [(n, h, m) for (h, m) in probe.ROOTS for n in NEW if (n, h, m) not in done]
    print(f"{len(NEW)} functions, {len(jobs)} probes to do ({len(done)} done)", flush=True)
    t0 = time.time()
    with Pool(workers, initializer=init_worker) as p, open(OUT, 'a') as f:
        for k, r in enumerate(p.imap_unordered(work, jobs), 1):
            f.write(json.dumps(r) + '\n'); f.flush()
            c = ct(r)
            msg = r.get('error') or f"kappa {r['kappa_rat']}"
            el = time.time() - t0
            print(f"  [{k}/{len(jobs)}] {r['name']} at {r['h']}/{r['m']}: {c} {msg}, {r['secs']} s | eta {el / k * (len(jobs) - k) / 3600:.2f} h", flush=True)
