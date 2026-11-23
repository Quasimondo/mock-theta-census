from ser import *
# fifth-order: F0 = sum q^{2n^2}/(q;q^2)_n ; F1 = sum q^{2n^2+2n}/(q;q^2)_{n+1}
# phi0 = sum q^{n^2}(-q;q^2)_n ; phi1 = sum q^{(n+1)^2}(-q;q^2)_n
odd = lambda k: [(2*k-1,-1,-1)]   # 1/(1-q^{2k-1})
F0s = eul(lambda n: 2*n*n, odd, -1)
F1s = mul(eul(lambda n: 2*n*n+2*n, lambda k: [(2*k+1,-1,-1)], -1), inv(add(one(), mono(1,-1))))  # /(1-q) for (q;q^2)_{n+1}
phi0s = eul(lambda n: n*n, lambda k: [(2*k-1,1,1)], -1)
phi1s = eul(lambda n: (n+1)**2, lambda k: [(2*k-1,1,1)], -1)
F0 = eul(lambda n: 2*n*n, odd, 1)
phi0 = eul(lambda n: n*n, lambda k: [(2*k-1,1,1)], 1)
print('F0*  ', F0s[:20]); print('F1*  ', F1s[:20])
# twist-run coefficients for comparison
print('F0* from twist run: 1, 0, -1, -1, -1, -1, -1, -1, 0, 0, 0, 1, 1, 1, 2, 2')
# classical check: F0 + phi0(-q) - 1 should be modular (theta quotient) -> test sparsity of its product exponents
B0 = add(F0, neg_q(phi0), mono(0,-1))
B0s = add(F0s, neg_q(phi0s), mono(0,-1))
import json; json.dump(dict(F0s=F0s,F1s=F1s,phi0s=phi0s,phi1s=phi1s,B0=B0,B0s=B0s),open('series.json','w'))
def prodexp(a, K=60):
    # a = prod (1-q^n)^{-c_n}; return c_n
    from fractions import Fraction as Fr
    f=a; ivf=inv(f); d=[i*f[i] for i in range(N)]
    lg=[sum(d[j]*ivf[i-j] for j in range(i+1)) for i in range(K)]
    c=[0]*K
    for k in range(1,K):
        c[k]=Fr(lg[k]-sum(n*c[n] for n in range(1,k) if k%n==0),k)
    return [str(x) for x in c[1:]]
print('classical bilateral F0+phi0(-q)-1:', B0[:12], 'product exps', prodexp(B0,30))
print('twisted   bilateral F0*+phi0*(-q)-1:', B0s[:16])
