import json, collections
from math import gcd
reps = {r['id']: r for r in json.load(open('cands.json'))}
def ct(r):
    if 'error' in r: return '?'
    if r['resid'] < 1e-25 and r['kappa_err'] < 1e-20: return 'E'
    if abs(r['kappa']) < 1e-4: return 'Q'
    return 'X'
res = collections.defaultdict(dict)
for f in ['stageA.jsonl','stageB1.jsonl','stageB2.jsonl','stageB3.jsonl','stageC.jsonl','quantum_followup.jsonl']:
    for l in open(f):
        d = json.loads(l); res[int(d['name'][1:])][(d['h'], d['m'])] = d
stageC = {o['id']: o['cls'] for o in json.load(open('stageC_summary.json'))}
label = {}
for i, r in reps.items():
    if r['known']: label[i] = 'known(catalog)'; continue
    rs = res[i]
    if i in stageC: label[i] = stageC[i]; continue
    order = [(0,1),(1,2),(1,3),(1,4),(1,6),(1,5),(2,5),(1,7),(2,7),(3,7),(1,8),(3,8),(1,10),(3,10),(1,12),(5,12),(1,14),(3,14),(5,14)]
    lab = None
    for o in order:
        if o in rs and ct(rs[o]) == 'X':
            lab = 'out@%d/%d' % o; break
    label[i] = lab or 'quantum'
json.dump(label, open('labels.json', 'w'))
c = collections.Counter(v if not v.startswith('out@') else v for v in label.values())
print(sorted(c.items(), key=lambda x: -x[1]))
