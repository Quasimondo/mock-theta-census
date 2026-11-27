"""Summarize the rank-2 census: class of every function that reached a verdict.

Writes report2.json and report2.md. Functions whose coefficients match a function from the
single-sum census (../census/cands.json) are marked with that match.
"""
import json, os, collections
from run2 import REPS, load, ct, STAGES

cells = collections.defaultdict(dict)
for st in "ABC":
    for r in load(st):
        cells[int(r["name"][1:])][(r["h"], r["m"])] = r

single = {}
p = os.path.join("..", "census", "cands.json")
if os.path.exists(p):
    for r in json.load(open(p)):
        single.setdefault(tuple(r["coeffs"]), r.get("known") or f"c{r['id']}: {r.get('label', '')}")

nroots = sum(len(v[0]) for v in STAGES.values())
out = collections.defaultdict(list)
for cid, cs in cells.items():
    t = {k: ct(r) for k, r in cs.items()}
    v = set(t.values())
    if "X" in v:
        cls = "neither"
    elif "R" in v:
        cls = "R type"
    elif len(t) < nroots:
        cls = "pending" if "?" not in v else "pending (errors)"
    elif "?" in v:
        cls = "undecided (errors)"
    elif v == {"E"}:
        cls = "modular"
    elif v == {"Q"}:
        cls = "quantum type"
    else:
        cls = "mock"
    rep = REPS[cid]
    kap = sorted({r["kappa_rat"] for r in cs.values() if ct(r) in "ER"})
    out[cls].append(dict(id=cid, label=rep["label"], coeffs=rep["coeffs"], kappas=kap,
                         single=single.get(tuple(rep["coeffs"])),
                         cells="".join(t.get(r, ".") for r in sorted(t))))
json.dump(out, open("report2.json", "w"), indent=1)
with open("report2.md", "w") as f:
    f.write("# Rank-2 census report\n\n")
    for cls in ["mock", "modular", "quantum type", "R type", "undecided (errors)", "pending", "pending (errors)", "neither"]:
        f.write(f"- {cls}: {len(out[cls])}\n")
    f.write("\nR type = rational exponent but inexact law at some root (power-series corrections, Nahm-sum type).\n")
    byk = collections.Counter(k for r in out["R type"] for k in r["kappas"])
    if byk:
        f.write("R type by exponent: " + ", ".join(f"κ={k}: {n}" for k, n in byk.most_common(15)) + "\n")
    for cls in ["mock", "modular", "undecided (errors)"]:
        f.write(f"\n## {cls}\n\n| id | formula | κ values | first coefficients | single-sum match |\n|---|---|---|---|---|\n")
        for r in sorted(out[cls], key=lambda r: r["id"]):
            f.write(f"| {r['id']} | {r['label']} | {', '.join(r['kappas'])} | {' '.join(map(str, r['coeffs'][:10]))} | {r['single'] or ''} |\n")
print(open("report2.md").read()[:3000])
