exec(open('F3.py').read().split("for tn,t in thetas.items():")[0])
import itertools
P={a:poch_inf(1,-1,a,a) for a in (1,2,4)}
Pinv={a:inv(P[a]) for a in (1,2,4)}
def eta(e):
    r=one()
    for a,k in zip((1,2,4),e):
        for _ in range(abs(k)): r=mul(r, P[a] if k>0 else Pinv[a])
    return r
targets={'B0*':S['B0s']}
# also the twisted bilateral for F1*: negative part of F1* = sum (-1)^n q^{2n^2+2n}/(q;q^2)_{n+1}
F1s=S['F1s']
for tname,B in targets.items():
  for e in itertools.product(range(-2,3),repeat=3):
    if e==(0,0,0): continue
    t=eta(e)
    for c in (1,2,-1,-2):
        Y=mul(sc(c,B), inv(t)); nn=norm(Y)
        if nn and nn[1] in L:
            print(tname,'*',c,'/ eta', e, '=', [ (k,v) for k,v in L[nn[1]]][:3], 'shift', nn[0], flush=True)
print('done')
