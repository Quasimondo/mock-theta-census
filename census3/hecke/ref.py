from sparse import *
def f0(n):
    if n*n>=N: return None
    return shift(qpoch(1,1,+1,inv=True,n=n), n*n)
def psi1(n):
    e=n*(n+1)//2
    if e>=N: return None
    return shift(qpoch(1,1,+1,n=n), e)
for nm,t in [('f0',series_sum(f0)),('psi1',series_sum(psi1))]:
    for ex in [(1,0,0),(0,1,0),(2,-1,0),(-1,1,0)]:
        c=mul(t,etaq(dict(zip((1,2,4),ex))))
        print(nm,ex,sum(1 for x in c if x),max(abs(x) for x in c))
