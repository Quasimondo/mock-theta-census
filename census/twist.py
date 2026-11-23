"""Sign twists of the catalog mock theta functions, plus F(z) at further z.
Each twisted function is the catalog Eulerian series with (-1)^n inserted.
Runs the 19-root growth test (probe.probe) on every new function."""
import sys, json, re, os
from multiprocessing import Pool
from flint import acb
import catalog, probe

src = open('catalog.py').read()
lines = [l for l in src.splitlines() if l.startswith('M("')]
NEW = []
def tw(e, ratio, **kw):
    return catalog.eul(e, lambda n, Q: -ratio(n, Q), **kw)
def M(name, order, fn):
    catalog.add(name + '*', 'twist', order, fn); NEW.append(name + '*')
env = dict(catalog.__dict__); env.update(M=M, tw=tw)
for l in lines:
    exec(l.replace('eul(', 'tw(', 1), env)

def Fz(zf):
    return catalog.eul(lambda n: n * (n + 1) // 2, lambda n, Q: zf() * (1 + Q(n)))
for nm, zf in [('F(i)', lambda: acb(0, 1)), ('F(omega)', lambda: (acb(2) / 3).exp_pi_i()),
               ('F(-omega)', lambda: -(acb(2) / 3).exp_pi_i())]:
    catalog.add(nm, 'Fz', '', Fz(zf)); NEW.append(nm)
catalog.BY = {e['name']: e for e in catalog.E}

def _nice():
    try:
        os.nice(19)
    except OSError:
        pass

if __name__ == '__main__':
    import os
    _nice()
    workers = int(sys.argv[sys.argv.index('--workers') + 1]) if '--workers' in sys.argv else 2
    names = [a for a in sys.argv[1:] if a in NEW] or NEW
    done = set()
    if os.path.exists('twist_results.jsonl'):
        for l in open('twist_results.jsonl'):
            try:
                d = json.loads(l); done.add((d['name'], d['h'], d['m']))
            except json.JSONDecodeError:
                pass
    jobs = [(n, h, m) for n in names for (h, m) in probe.ROOTS if (n, h, m) not in done]
    print(f"{len(jobs)} probes to do ({len(done)} already done)", flush=True)
    with Pool(workers, initializer=_nice) as p, open('twist_results.jsonl', 'a') as f:
        for k, r in enumerate(p.imap_unordered(probe.probe, jobs), 1):
            f.write(json.dumps(r) + '\n'); f.flush()
            print(f"  [{k}/{len(jobs)}] {r['name']} at {r['h']}/{r['m']}" + (f": error {r['error']}" if 'error' in r else f": kappa {r.get('kappa_rat', r.get('kappa'))}"), flush=True)
