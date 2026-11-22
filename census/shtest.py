import json, sys, catalog, search, shadow, mpmath as mp
from math import gcd
order=[(h,m) for m in [1,2,3,4,5,6,7,8,10,12,14] for h in range(0,m//2+1) if gcd(h,m)==1]
S={o['id']:o for o in json.load(open('stageC_summary.json'))}
reps={r['id']:r for r in json.load(open('cands.json'))}
def fnc(i):
    c=dict(reps[i]); c['pat']=[tuple(x) for x in c['pat']]; return search.make_fn(c)
cid, kn = int(sys.argv[1]), sys.argv[2]
roots=[order[k] for k,t in enumerate(S[cid]['types']) if t=='Q'][:2]
u=shadow.vector(fnc(cid),roots); v=shadow.vector(catalog.BY[kn]['fn'],roots) if kn in catalog.BY else shadow.vector(fnc(int(kn)),roots)
c,mis=shadow.proportional(u,v)
print(cid,'vs',kn,'roots',roots,'c=',mp.nstr(c,10),'misfit',mp.nstr(mis,3),flush=True)
