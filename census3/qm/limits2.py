"""Radial limits of the 1,357 census functions bounded at all 19 roots, at all roots of order <= 12.
    python limits2.py --workers 2   -> limits2.jsonl (resumable)"""
import sys, os, json, time, signal
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '../../census'))
os.chdir(os.path.join(HERE, '../../census'))
from math import gcd
from multiprocessing import Pool
import mpmath as mp
from flint import acb, arb, ctx
import catalog, search
ctx.prec = 160
K = 12
ROOTS = [(0, 1)] + [(h, m) for m in range(2, 13) for h in range(1, m) if gcd(h, m) == 1]
labels = json.load(open('labels.json'))
C = {c['id']: c for c in json.load(open('cands.json'))}
IDS = sorted(int(k) for k, v in labels.items() if v == 'quantum')
OUT = os.path.join(HERE, 'limits2.jsonl')
def fn_of(i):
    c = dict(C[i]); c['pat'] = [tuple(x) for x in c['pat']]
    return search.make_fn(c)
def neville(ts, vals):
    P = list(vals)
    for k in range(1, len(P)):
        for i in range(len(P) - k):
            P[i] = (ts[i + k] * P[i] - ts[i] * P[i + 1]) / (ts[i + k] - ts[i])
    return P[0]
class TO(BaseException): pass
def _al(*a): raise TO()
def init():
    os.nice(19); signal.signal(signal.SIGALRM, _al)
def work(i):
    fn = fn_of(i); res = {}
    t0 = time.time()
    for (h, m) in ROOTS:
        ts = [0.003 / m ** 2 * j for j in range(1, K + 1)]
        signal.alarm(60)
        try:
            vals = []
            for t in ts:
                v = fn((acb(2 * h) / m).exp_pi_i() * arb(-t).exp())
                vals.append(complex(float(v.real.mid()), float(v.imag.mid())))
            signal.alarm(0)
        except (TO, Exception) as e:
            signal.alarm(0); res[f"{h}/{m}"] = None; continue
        a = neville(ts, vals); b = neville(ts[:-1], vals[:-1])
        res[f"{h}/{m}"] = [a.real, a.imag, abs(a - b)]
    return dict(id=i, lim=res, secs=round(time.time() - t0, 1))
if __name__ == '__main__':
    w = int(sys.argv[sys.argv.index('--workers') + 1]) if '--workers' in sys.argv else 2
    done = set()
    if os.path.exists(OUT):
        done = {json.loads(l)['id'] for l in open(OUT)}
    todo = [i for i in IDS if i not in done]
    print(len(todo), 'functions to do', flush=True)
    t0 = time.time()
    with Pool(w, initializer=init) as p, open(OUT, 'a') as f:
        for k, r in enumerate(p.imap_unordered(work, todo), 1):
            f.write(json.dumps(r) + '\n'); f.flush()
            if k % 10 == 0:
                el = time.time() - t0
                print(f"{k}/{len(todo)} eta {el / k * (len(todo) - k) / 3600:.2f} h", flush=True)
