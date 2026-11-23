from ser import *
import itertools, json
S=json.load(open('series.json'))
B0s=S['B0s']
def theta4(m=1):
    a=z(); n=0
    while m*n*n<N:
        a[m*n*n]+= 1 if n==0 else 2*(-1)**n
        n+=1
    return a
def theta3(m=1): return neg_q(theta4(m)) if m==1 else sub_q(neg_q(theta4(1)),m)
def poch_inf(k,s,step=1,start=1):
    a=one(); j=start
    while j<N:
        a=mulf(a,j,s); j+=step
    return a
thetas={'theta4(q)':theta4(1),'theta4(q^2)':theta4(2),'theta3(q)':theta3(1),'theta3(q^2)':sub_q(theta3(1),2),
        '(q;q)inf':poch_inf(1,-1),'(q^2;q^2)inf':poch_inf(1,-1,2,2),'(q;q^2)inf':poch_inf(1,-1,2,1),
        '(-q;q^2)inf':poch_inf(1,1,2,1),'theta4(q)/(-q;q)inf... skip':None}
FAC={'(q;q)':lambda k:[(k,-1,-1)],'(q^2;q^2)':lambda k:[(2*k,-1,-1)],'(q;q^2)':lambda k:[(2*k-1,-1,-1)],
     '(-q;q)':lambda k:[(k,1,-1)],'(-q;q^2)':lambda k:[(2*k-1,1,-1)],'(q^4;q^4)':lambda k:[(4*k,-1,-1)]}
lib={}
for A in range(1,9):
  for B in range(-A,A+5):
    if (A+B)%2: continue
    for fn,fac in FAC.items():
      for sg in (1,-1):
        s=eul(lambda n,A=A,B=B:(A*n*n+B*n)//2, fac, sg)
        lib[(A,B,fn,sg)]=s
def norm(a):
    nz=[i for i,x in enumerate(a) if x]
    if not nz: return None
    v=nz[0]; b=a[v:v+150]
    if b[0]<0: b=[-x for x in b]
    return (v,tuple(b))
L={}
for k,s in lib.items():
    nn=norm(s)
    if nn: L.setdefault(nn[1],[]).append((k,nn[0]))
for tn,t in thetas.items():
    if t is None: continue
    if t[0]!=1: continue
    for c in (1,2,-1,-2):
        Y=mul(sc(c,B0s), inv(t))  # c*B0s/theta
        nn=norm(Y)
        if nn and nn[1] in L:
            print('B0* =', '1/%d'%c if abs(c)!=1 else ('-' if c<0 else ''), tn, '*', [(f"sum {('(-1)^n ' if k[3]<0 else '')}q^(({k[0]}n^2+{k[1]}n)/2)/{k[2]}_n", 'shift q^%d'%(v)) for k,v in L[nn[1]]], '(leading shift %d)'%nn[0])
print('done')
