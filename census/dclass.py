import json, sys, collections
def ct(r):
    if 'error' in r: return '?'
    if r['resid'] < 1e-25 and r['kappa_err'] < 1e-20: return 'E'
    if abs(r['kappa']) < 1e-4: return 'Q'
    return 'X'
ids = json.load(open(sys.argv[1])); files = sys.argv[2].split(',')
R = collections.defaultdict(list)
for f in files:
    for l in open(f):
        d = json.loads(l); R[int(d['name'][1:])].append(d)
out = {'X': [], 'E': [], 'Q': [], '?': []}
for i in ids:
    t = set(ct(r) for r in R[i])
    k = 'X' if 'X' in t else 'E' if 'E' in t else 'Q' if t == {'Q'} else '?'
    out[k].append(i)
print({k: len(v) for k, v in out.items()})
json.dump(out['Q'] + out['?'], open(sys.argv[3], 'w'))
json.dump(out['E'], open(sys.argv[4], 'w'))
