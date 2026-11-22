# (a) telescoping proof of H(x) = (1-x^2) g2(x,q) + 1 - x, checked symbolically term by term
from sympy import symbols, simplify, Rational
x, q = symbols('x q')
def poch(a, n):
    r = 1
    for k in range(n): r *= (1 - a*q**k)
    return r
ok = True
for n in range(0, 6):
    h  = q**(n*(n-1)//2) * poch(-q, n) / (poch(x*q, n) * poch(q/x, n))
    h1 = q**((n+1)*n//2) * poch(-q, n+1) / (poch(x*q, n+1) * poch(q/x, n+1))
    t  = poch(-q, n) * q**(n*(n+1)//2) / (poch(x, n+1) * poch(q/x, n+1))
    U, U1 = (1 - x*q**n)*h, (1 - x*q**(n+1))*h1
    ok &= simplify(U - U1 - (h - (1 - x**2)*t)) == 0
print("telescoping step U_n - U_{n+1} = h_n - (1-x^2) t_n, n=0..5:", ok, "  (U_0 = 1-x, U_n -> 0)")

# (b) double-sum formulas, exact series to q^N
from flint import fmpz_poly
import json
N = 300
S = json.load(open("series.json"))
def mul(a,b): return fmpz_poly((a*b).coeffs()[:N])
def inv(a):
    co = a.coeffs() + [0]*N; b = [0]*N; b[0] = 1
    for n in range(1, N):
        b[n] = -sum(co[k]*b[n-k] for k in range(1, n+1) if co[k])
    return fmpz_poly(b)
mq = fmpz_poly([1])            # (-q)_inf
for k in range(1, N): mq = mul(mq, fmpz_poly([1]+[0]*(k-1)+[1]))
def dbl(lin, par):
    # sum over k = par mod 2, n>=0 of (-1)^n q^{k(k+1)/2 + n^2 + lin*n + k n}/(q)_k
    tot = fmpz_poly([0]); P = fmpz_poly([1]); k = 0
    while k*(k+1)//2 < N:
        if k > 0: P = mul(P, fmpz_poly([1]+[0]*(k-1)+[-1]))
        if k % 2 == par:
            th = [0]*N; n = 0
            while k*(k+1)//2 + n*n + lin*n + k*n < N:
                th[k*(k+1)//2 + n*n + lin*n + k*n] += (-1)**n; n += 1
            tot += mul(fmpz_poly(th), inv(P))
        k += 1
    return fmpz_poly(tot.coeffs()[:N])
l1 = mul(mq, fmpz_poly(S["psi1s"])); l0 = mul(mq, fmpz_poly(S["psi0s"]))
print("(-q)_inf psi1* = sum_{k even,n} (-1)^n q^{k(k+1)/2+n^2+n+kn}/(q)_k :", l1 == dbl(1, 0))
print("(-q)_inf psi0* = sum_{k odd ,n} (-1)^n q^{k(k+1)/2+n^2+kn}/(q)_k   :", l0 == dbl(0, 1))
