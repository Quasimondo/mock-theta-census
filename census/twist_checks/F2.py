from ser import *
from fractions import Fraction as Fr
import itertools
def prodexp(a, K=60):
    if a[0]!=1: return None
    f=a; ivf=inv(f); d=[i*f[i] for i in range(N)]
    lg=[sum(d[j]*ivf[i-j] for j in range(i+1)) for i in range(K)]
    c=[0]*K
    for k in range(1,K):
        c[k]=Fr(lg[k]-sum(n*c[n] for n in range(1,k) if k%n==0),k)
    return c[1:]
def periodic(c):
    for P in range(1,21):
        if all(c[i]==c[i+P] for i in range(len(c)-P)): return P
    return None
odd = lambda k: [(2*k-1,-1,-1)]
F0 = eul(lambda n: 2*n*n, odd, 1)
F1 = mul(eul(lambda n: 2*n*n+2*n, lambda k: [(2*k+1,-1,-1)], 1), inv(add(one(), mono(1,-1))))
phi0 = eul(lambda n: n*n, lambda k: [(2*k-1,1,1)], 1)
phi1 = eul(lambda n: (n+1)**2, lambda k: [(2*k-1,1,1)], 1)
parts = {'phi0(q)':phi0,'phi0(-q)':neg_q(phi0),'phi1(q)':phi1,'phi1(-q)':neg_q(phi1)}
for nameX,X in (('F0',F0),('F1',F1)):
  for (pn,P),c,k in itertools.product(parts.items(), (1,-1,2,-2), (0,1,-1)):
    S = add(X, sc(c,P), mono(0,k))
    # allow leading monomial: shift
    nz=[i for i,x in enumerate(S) if x]
    if not nz: continue
    v=nz[0]; S2=S[v:]+[0]*v
    if S2[0] not in (1,-1,2,-2): continue
    lead=S2[0]
    if any(x % lead for x in S2[:120]): continue
    S3=[x//lead for x in S2]
    ex=prodexp(S3,45)
    if ex is None: continue
    p=periodic(ex)
    if p: print(nameX, '+', c, '*', pn, '+', k, ' = ', lead, '* q^%d *'%v, 'product with periodic exponents, period', p, [str(x) for x in ex[:p]])
