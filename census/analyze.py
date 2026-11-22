import json, collections, sys
import catalog

R = json.load(open(sys.argv[1] if len(sys.argv) > 1 else "results.json"))
lab = {e["name"]: e for e in catalog.E}
by = collections.defaultdict(list)
for r in R:
    by[r["name"]].append(r)

def cusp_type(r):
    if "error" in r:
        return "?"
    if r["resid"] < 1e-25 and r["kappa_err"] < 1e-20:
        return "E"          # exact theta-like law at this cusp
    if abs(r["kappa"]) < 1e-4:
        return "Q"          # bounded, with a genuine power series left over (quantum signature)
    return "X"              # growth with an irrational exponent or an inexact law

def predict(types):
    s = set(types)
    if "X" in s: return "neither"
    s.discard("?")
    if not s: return "incomplete"
    if s == {"E"}: return "modular"
    if s == {"Q"}: return "false/quantum"
    return "mock"

rows = []
for e in catalog.E:
    rs = sorted(by[e["name"]], key=lambda r: r["h"] / r["m"])
    types = [cusp_type(r) for r in rs]
    kap = ["%s" % (r.get("kappa_rat") if cusp_type(r) == "E" else ("~0" if cusp_type(r) == "Q" else ("%.6f" % r["kappa"] if "kappa" in r else "?"))) for r in rs]
    off = sorted({r["gamma_rat"] for r in rs if cusp_type(r) == "E"})
    wt = sorted({r["weight_rat"] for r in rs if cusp_type(r) == "E"})
    rows.append(dict(name=e["name"], label=e["label"], order=e["order"], types="".join(types),
                     pred=predict(types), kappas=kap, offsets=off, weights=wt,
                     roots=["%d/%d" % (r["h"], r["m"]) for r in rs],
                     minres=min((r["resid"] for r in rs if "resid" in r), default=None),
                     maxres=max((r["resid"] for r in rs if "resid" in r), default=None)))
json.dump(rows, open("summary.json", "w"), indent=1)

conf = collections.Counter((r["label"], r["pred"]) for r in rows)
labels = ["modular", "mixed", "mock", "false", "other", "neither"]
preds = ["modular", "mock", "false/quantum", "neither", "incomplete"]
print("%-9s" % "label" + "".join("%15s" % p for p in preds))
for l in labels:
    print("%-9s" % l + "".join("%15d" % conf[(l, p)] for p in preds))
print()
for r in rows:
    flag = "" if (r["label"], r["pred"]) in {("modular", "modular"), ("mock", "mock"), ("false", "false/quantum"), ("neither", "neither")} else "   <-- check"
    print(f"{r['name'][:34]:34s} {r['label']:8s} {r['types']}  {r['pred']:13s} kappa={' '.join(r['kappas'])}  w={','.join(r['weights'])} off={','.join(r['offsets'])}{flag}")
