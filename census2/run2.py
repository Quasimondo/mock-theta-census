"""Staged, resumable, low-priority census of the rank-2 candidates in cands2.json.

Stage A: q -> 1 only.          Drops every function without an exact theta-type law there.
Stage B: roots of order 2-5.   Same test on the survivors.
Stage C: all 19 roots.         Final classification.

Every finished probe is appended to stageX.jsonl at once, so the run can be stopped
(Ctrl-C, reboot) and restarted with the same command; finished probes are skipped.

Usage:
    python3 run2.py --workers 4            # all stages, 4 worker processes, lowest CPU priority
    python3 run2.py --workers 4 --stage A  # one stage only
    python3 run2.py --bench                # time a few probes and estimate the total run time
    python3 run2.py --workers 4 --retry    # also redo earlier errors/timeouts
    python3 run2.py --workers 4 --old-eval # use the original sum2.py evaluator instead of sum2_fast.py

Within a stage the roots are run in order of increasing m (cheap first), and as soon as a function
gets an X or R the rest of its probes in that stage are skipped.

Verdict per probe (ct):
    E  exact theta-type law, rational exponent kappa (modular or mock behaviour at this root)
    Q  bounded at this root
    R  rational kappa, but the law is not exact: power-series corrections (Nahm-sum type, like
       R0, R1 at q = -1 in the psi* identities). Not modular or mock, so dropped from the mock
       search, but listed separately by report2.py.
    X  no rational kappa: dropped
    ?  error or time limit
"""
import json, math, os, sys, time, signal, random
from fractions import Fraction
from multiprocessing import Pool, Array
import catalog, probe
if "--old-eval" in sys.argv:
    import sum2
else:
    import sum2_fast as sum2

STAGES = {
    "A": ([(0, 1)], 120),
    "B": ([(1, 2), (1, 3), (1, 4), (1, 5), (2, 5)], 600),
    "C": ([r for r in probe.ROOTS if r not in [(0, 1), (1, 2), (1, 3), (1, 4), (1, 5), (2, 5)]], 1800),
}
REPS = {r["id"]: r for r in json.load(open("cands2.json"))}
RETRY = "--retry" in sys.argv

LIMITS = (120, 600, 1800)
RDEN = 500  # largest denominator accepted as a 'rational' kappa for the R class

def ct(r):
    if "error" in r:
        return "?"
    if any(abs(r.get("secs", 0) - L) < 1.0 for L in LIMITS):
        return "?"  # cut off by the time limit: the fit used partial data
    if r["resid"] < 1e-25 and r["kappa_err"] < 1e-20:
        return "E"
    if r.get("supp") and r["kappa_err"] < 1e-8 and r["resid"] < 1e3 * math.exp(-r["supp"]):
        return "E"  # exact law up to the bounded part that a rim-limited pass cannot suppress
    if abs(r["kappa"]) < 1e-4:
        return "Q"
    if r["kappa_err"] < 1e-8 and Fraction(r["kappa_rat"]).denominator <= RDEN:
        return "R"
    return "X"

class Timeout(BaseException):
    def __str__(self):
        return "timeout"

def _alarm(*a):
    raise Timeout()

DEAD = None  # shared flags, one per function id: set by the parent when the function got an X

def init_worker(dead=None):
    global DEAD
    DEAD = dead
    try:
        os.nice(19)
    except OSError:
        pass
    signal.signal(signal.SIGINT, signal.SIG_IGN)

def work(job):
    cid, h, m, tlim = job
    name = f"d{cid}"
    if DEAD is not None and DEAD[cid]:
        return dict(name=name, h=h, m=m, skipped=True)  # already out at another root of this stage
    if name not in catalog.BY:
        catalog.BY.clear()  # one function at a time is enough; keeps the worker's memory flat
        catalog.BY[name] = dict(fn=sum2.make_fn(REPS[cid]))
    signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(tlim)
    t0 = time.time()
    try:
        out = probe.probe((name, h, m))
    except Timeout:
        out = dict(name=name, h=h, m=m, error="timeout")
    except Exception as ex:
        out = dict(name=name, h=h, m=m, error=str(ex)[:80])
    finally:
        signal.alarm(0)
    out["secs"] = round(time.time() - t0, 1)
    return out

def load(stage):
    fn = f"stage{stage}.jsonl"
    rows = []
    if os.path.exists(fn):
        for l in open(fn):
            try:
                rows.append(json.loads(l))
            except json.JSONDecodeError:
                pass  # half-written last line after a hard stop
    return rows

def survivors(upto):
    """ids with only E/Q cells in the stages up to and including `upto`.
    Functions with an X are out; functions with a ? (error or time limit) are parked as
    undecided and listed by report2.py, so they can be rerun later with a longer limit."""
    ids = set(REPS)
    for st in "ABC":
        for r in load(st):
            if ct(r) in "XR?":
                ids.discard(int(r["name"][1:]))
        if st == upto:
            break
    return ids

def run_stage(stage, workers):
    roots, tlim = STAGES[stage]
    prev = {"A": None, "B": "A", "C": "B"}[stage]
    ids = sorted(REPS) if prev is None else sorted(survivors(prev))
    if RETRY:
        keep = [r for r in load(stage) if ct(r) != "?"]
        dropped = len(load(stage)) - len(keep)
        with open(f"stage{stage}.jsonl", "w") as f:
            for r in keep:
                f.write(json.dumps(r) + "\n")
        if dropped:
            print(f"  --retry: {dropped} earlier errors/timeouts in stage {stage} will be redone", flush=True)
    done = {(int(r["name"][1:]), r["h"], r["m"]) for r in load(stage)}
    jobs = [(i, h, m, tlim) for i in ids for (h, m) in roots if (i, h, m) not in done]
    out = {int(r["name"][1:]) for r in load(stage) if ct(r) in "XR"}  # an X or R is final; a ? is not
    jobs = [j for j in jobs if j[0] not in out]
    rng = random.Random(0)
    jobs.sort(key=lambda j: (j[2], rng.random()))  # cheap roots (small m) first, so an X there saves the expensive ones
    dead = Array("b", max(REPS) + 1, lock=False)
    print(f"stage {stage}: {len(ids)} functions, {len(jobs)} probes to do ({len(done)} already done)", flush=True)
    if not jobs:
        return
    t0, k, skipped = time.time(), 0, 0
    tally = {"E": 0, "Q": 0, "R": 0, "X": 0, "?": 0}
    WORD = {"E": "exact law", "Q": "bounded", "R": "rational, inexact: out", "X": "out", "?": "error"}
    print(f"  {workers} workers started; one line per finished probe "
          "(E = exact theta-type law, Q = bounded, R = rational exponent but inexact law, X = no rational exponent, ? = error/timeout)", flush=True)
    with Pool(workers, initializer=init_worker, initargs=(dead,)) as pool, open(f"stage{stage}.jsonl", "a") as f:
        try:
            for r in pool.imap_unordered(work, jobs):
                if r.get("skipped"):
                    skipped += 1
                    continue
                f.write(json.dumps(r) + "\n"); f.flush()
                k += 1
                c = ct(r)
                tally[c] += 1
                if c in "XR":
                    dead[int(r["name"][1:])] = 1
                el = time.time() - t0
                kap = f" κ={r['kappa_rat']}" if c in "ER" else ""
                err = f" ({r.get('error', 'hit the time limit')})" if c == "?" else ""
                left = len(jobs) - k - skipped
                print(f"  [{k + skipped}/{len(jobs)}] {r['name']} at {r['h']}/{r['m']}: {c} {WORD[c]}{kap}{err}, {r['secs']} s"
                      f"   | E {tally['E']} Q {tally['Q']} R {tally['R']} X {tally['X']} ? {tally['?']} skipped {skipped}"
                      f" | eta {el/k*left/3600:.2f} h", flush=True)
        except KeyboardInterrupt:
            pool.terminate()
            print("stopped; rerun the same command to continue", flush=True)
            sys.exit(1)

def bench(workers):
    ids = random.Random(1).sample(sorted(REPS), min(12, len(REPS)))
    t0 = time.time()
    with Pool(workers, initializer=init_worker) as pool:
        res = pool.map(work, [(i, 0, 1, STAGES["A"][1]) for i in ids])
    wall = time.time() - t0
    per = wall * workers / len(ids)
    nA = len(REPS)
    print(f"{len(ids)} stage-A probes in {wall:.0f} s with {workers} workers ({per:.1f} core-seconds each)")
    print(f"stage A estimate: {nA} probes, about {nA*per/workers/3600:.1f} h with {workers} workers")
    print("stages B and C depend on how many functions survive A (in the single-sum census about 4 in 10),")
    print("and probes at higher-order roots cost 10-100 times more than at q -> 1.")

if __name__ == "__main__":
    workers = int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 2
    try:
        os.nice(19)
    except OSError:
        pass
    if "--bench" in sys.argv:
        bench(workers)
        sys.exit()
    stages = sys.argv[sys.argv.index("--stage") + 1] if "--stage" in sys.argv else "ABC"
    for st in stages:
        run_stage(st, workers)
    print("done; run  python3 report2.py  for the summary", flush=True)
