from ser import *
import json, sympy as sp
S=json.load(open('series.json'))
odd=lambda k:[(2*k-1,-1,-1)]
neg=lambda k:[(k,1,1)]   # (1+q^k)
psi1s=eul(lambda n:n*(n+1)//2, neg, -1); psi0s=eul(lambda n:(n+1)*(n+2)//2, neg, -1)
f0s=eul(lambda n:n*n, lambda k:[(k,1,-1)], -1); f1s=eul(lambda n:n*n+n, lambda k:[(k,1,-1)], -1)
phi0s=S['phi0s']; phi1s=S['phi1s']
basis={'1':one(),'psi1*':psi1s,'psi0*':psi0s,'psi1*(q^2)':sub_q(psi1s,2),'psi0*(q^2)':sub_q(psi0s,2),
       'f0*':f0s,'f1*':f1s,'f0*(q^2)':sub_q(f0s,2),'f1*(q^2)':sub_q(f1s,2),
       'phi0*(q)':phi0s,'phi0*(-q)':neg_q(phi0s),'phi1*(q)':phi1s,'phi1*(-q)':neg_q(phi1s)}
for tn in ('F0s','F1s'):
    T=S[tn]
    names=list(basis); M=sp.Matrix([[basis[b][i] for b in names]+[T[i]] for i in range(200)])
    ns=M.nullspace()
    print(tn, 'relations found:', len(ns))
    for v in ns:
        if v[-1]!=0:
            v=v/v[-1]; print('  ', tn, '=', ' + '.join(f'{-v[i]}*{names[i]}' for i in range(len(names)) if v[i]!=0))
