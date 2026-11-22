import json, sys, catalog, search, shadow, mpmath as mp
from math import gcd
order=[(h,m) for m in [1,2,3,4,5,6,7,8,10,12,14] for h in range(0,m//2+1) if gcd(h,m)==1]
S={o['id']:o for o in json.load(open('stageC_summary.json'))}
reps={r['id']:r for r in json.load(open('cands.json'))}
def fnc(i):
    c=dict(reps[i]); c['pat']=[tuple(x) for x in c['pat']]; return search.make_fn(c)
cid=int(sys.argv[1]); basis=sys.argv[2].split(',')
roots=[order[k] for k,t in enumerate(S[cid]['types']) if t=='Q'][:2]
def vec(fn):
    out=[]
    for h,m in roots:
        a=shadow.series_at(fn,h,m); out+=a[1:6]
    return out
u=vec(fnc(cid)); V=[vec(catalog.BY[b]['fn']) for b in basis]
# complex least squares u ~ sum c_j V_j  via real embedding
rows=[];rhs=[]
for k in range(len(u)):
    rows.append([V[j][k].real for j in range(len(V))]+[-V[j][k].imag for j in range(len(V))]); rhs.append(u[k].real)
    rows.append([V[j][k].imag for j in range(len(V))]+[V[j][k].real for j in range(len(V))]); rhs.append(u[k].imag)
x,res=mp.qr_solve(mp.matrix(rows),mp.matrix(rhs))
nb=len(V); c=[mp.mpc(x[j],x[nb+j]) for j in range(nb)]
mis=res/mp.sqrt(sum(abs(z)**2 for z in u))
print(cid,'~',' + '.join(f"({mp.nstr(cj,8)})·{b}" for cj,b in zip(c,basis)),'roots',roots,'misfit',mp.nstr(mis,3),flush=True)
