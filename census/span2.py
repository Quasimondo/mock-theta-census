import json, sys, catalog, search, shadow, mpmath as mp
reps={r['id']:r for r in json.load(open('cands.json'))}
def fnc(i):
    c=dict(reps[i]); c['pat']=[tuple(x) for x in c['pat']]; return search.make_fn(c)
target=fnc(8438); basis=[('psi1',catalog.BY['psi1']['fn']),('psi1*',fnc(830))]
roots=[(1,2),(1,4)]
def vec(fn):
    out=[]
    for h,m in roots:
        a=shadow.series_at(fn,h,m); out+=a[1:6]
    return out
u=vec(target); V=[vec(f) for _,f in basis]
rows=[];rhs=[]
for k in range(len(u)):
    rows.append([V[j][k].real for j in range(2)]+[-V[j][k].imag for j in range(2)]); rhs.append(u[k].real)
    rows.append([V[j][k].imag for j in range(2)]+[V[j][k].real for j in range(2)]); rhs.append(u[k].imag)
x,res=mp.qr_solve(mp.matrix(rows),mp.matrix(rhs))
print('coeffs', [mp.nstr(mp.mpc(x[j],x[2+j]),12) for j in range(2)], 'misfit', mp.nstr(res/mp.sqrt(sum(abs(z)**2 for z in u)),3), flush=True)
