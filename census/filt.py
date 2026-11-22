import json, sys, collections
def ct(r):
    if 'error' in r: return '?'
    if r['resid'] < 1e-25 and r['kappa_err'] < 1e-20: return 'E'
    if abs(r['kappa']) < 1e-4: return 'Q'
    return 'X'
files = sys.argv[1].split(',')
R = {}
for fn in files:
    for l in open(fn):
        d = json.loads(l); R.setdefault(d['name'], []).append(d)
types = {n: ''.join(ct(r) for r in rs) for n, rs in R.items()}
print(collections.Counter('X' in t and 'X' or ('?' in t and '?' or 'ok') for t in types.values()))
surv = [int(n[1:]) for n, t in types.items() if 'X' not in t]
json.dump(surv, open(sys.argv[2], 'w'))
print(len(surv), 'survivors ->', sys.argv[2])
