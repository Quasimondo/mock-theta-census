"""Run probe on candidate representatives at given roots, append JSON lines (resumable)."""
import json, sys, signal, os
from multiprocessing import Pool
import catalog, probe, search

reps = json.load(open("cands.json"))
for r in reps:
    r["pat"] = [tuple(x) for x in r["pat"]]
    catalog.BY[f"c{r['id']}"] = dict(fn=search.make_fn(r))

class TO(Exception): pass
def handler(*a): raise TO()

def work(job):
    signal.signal(signal.SIGALRM, handler)
    signal.alarm(job[3])
    try:
        out = probe.probe(job[:3])
    except TO:
        out = dict(name=job[0], h=job[1], m=job[2], error="timeout")
    except Exception as ex:
        out = dict(name=job[0], h=job[1], m=job[2], error=str(ex)[:80])
    finally:
        signal.alarm(0)
    return out

if __name__ == "__main__":
    outfile, roots, ids = sys.argv[1], json.loads(sys.argv[2]), sys.argv[3]
    tlim = int(sys.argv[4]) if len(sys.argv) > 4 else 30
    names = [f"c{r['id']}" for r in reps if not r["known"]] if ids == "unknown" else [f"c{i}" for i in json.load(open(ids))]
    done = set()
    if os.path.exists(outfile):
        for l in open(outfile):
            d = json.loads(l); done.add((d["name"], d["h"], d["m"]))
    jobs = [(n, h, m, tlim) for n in names for h, m in roots if (n, h, m) not in done]
    print(len(jobs), "jobs", flush=True)
    with Pool(2) as p, open(outfile, "a") as f:
        for r in p.imap_unordered(work, jobs, chunksize=4):
            f.write(json.dumps(r) + "\n"); f.flush()
