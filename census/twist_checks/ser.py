# exact integer power series to N terms
N = 260
def z(): return [0]*N
def one(): a=z(); a[0]=1; return a
def mul(a,b):
    r=z()
    for i,x in enumerate(a):
        if x:
            for j in range(N-i):
                if b[j]: r[i+j]+=x*b[j]
    return r
def add(*xs):
    r=z()
    for a in xs:
        for i in range(N): r[i]+=a[i]
    return r
def sc(c,a): return [c*x for x in a]
def mono(k,c=1):
    a=z()
    if 0<=k<N: a[k]=c
    return a
def mulf(a,k,s):  # a*(1+s q^k)
    r=a[:]
    for i in range(N-1,k-1,-1): r[i]+=s*a[i-k]
    return r
def divf(a,k,s):  # a/(1+s q^k)
    r=a[:]
    for i in range(k,N): r[i]-=s*r[i-k]
    return r
def sub_q(a,m):  # a(q^m)
    r=z()
    for i,x in enumerate(a):
        if i*m<N: r[i*m]=x
    return r
def neg_q(a): return [x*(-1)**i for i,x in enumerate(a)]
def inv(a):
    r=z(); r[0]=1
    assert a[0]==1
    for i in range(1,N): r[i]=-sum(a[j]*r[i-j] for j in range(1,i+1))
    return r
def eul(e, facs, sign=1, n0=0, nmax=200):
    """sum_{n>=n0} sign^n q^{e(n)} prod_{k<=n} facs(k): facs(k) -> list of (kexp, s, power)"""
    tot=z(); P=one()
    for n in range(0, nmax):
        if n>0:
            for (kk,s,p) in facs(n):
                P = mulf(P,kk,s) if p>0 else divf(P,kk,s)
        if n>=n0:
            E=e(n)
            if E>=N: 
                if n>5: break
                continue
            t=[0]*E+P[:N-E]
            tot=add(tot, sc(sign**n, t))
    return tot
