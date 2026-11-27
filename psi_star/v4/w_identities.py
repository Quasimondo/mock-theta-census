"""Identities for W(q) = sum_{j>=0} (-1)^j q^(j^2+j) (q;q^2)_j and its partner, exact to q^(N-1).
 G(z) = sum z^j q^(j^2) (q;q^2)_j ;  W = G(-q),  W3 = G(-q^3)
 (V)   W(q)  + sum_{n>=1} q^(2n^2-n)/(q;q^2)_n   =        psi(q) * sum_k q^(2k^2+k)/(q^2;q^2)_k
 (VI)  q W3(q) + sum_{n>=1} q^(2n^2-3n+1)/(q;q^2)_n =     psi(q) * sum_k q^(2k^2-k)/(q^2;q^2)_k
 psi(q) = sum_{n>=0} q^(n(n+1)/2) = (q^2;q^2)_inf/(q;q^2)_inf
 q-difference equation: G(z) = 1 + z q G(z q^2) - z q^2 G(z q^4)
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qs import *
num = lambda P, n: mfac(P, 2*n-1, -1)            # P_n = P_{n-1} (1 - q^(2n-1))   -> (q;q^2)_n
den = lambda P, n: dfac(P, 2*n-1, -1)            # 1/(q;q^2)_n
den2 = lambda P, n: dfac(P, 2*n, -1)             # 1/(q^2;q^2)_n
def G(shift, sign):                               # sum sign^j q^(j^2 + shift*j) (q;q^2)_j
    return eul(lambda j: j*j + shift*j, num, sign)
W, W3 = G(1, -1), G(3, -1)
psi = zero(); n = 0
while n*(n+1)//2 < N: psi[n*(n+1)//2] += 1; n += 1
assert psi == mul(poch_inf(2, 2, -1), poch_inf(1, 2, -1, inv=True))
def tail(expo):                                   # sum_{n>=1} q^expo(n)/(q;q^2)_n
    s = eul(expo, den); s[expo(0)] -= 1; return s
Fa = tail(lambda n: 2*n*n - n)
Fb = tail(lambda n: 2*n*n - 3*n + 1)
Rp = eul(lambda k: 2*k*k + k, den2)               # sum q^(2k^2+k)/(q^2;q^2)_k
Rm = eul(lambda k: 2*k*k - k, den2)
print("(V)   W + sum_{n>=1} q^(2n^2-n)/(q;q^2)_n = psi(q) sum q^(2k^2+k)/(q^2;q^2)_k :", add(W, Fa) == mul(psi, Rp))
print("(VI)  q W3 + sum_{n>=1} q^(2n^2-3n+1)/(q;q^2)_n = psi(q) sum q^(2k^2-k)/(q^2;q^2)_k :", add(mono(W3, 1), Fb) == mul(psi, Rm))
# q-difference equation at z = -q and z = -1/q
W5 = G(5, -1); Wm = G(-1, -1)
lhs = W; rhs = add(add(one(), mono(W3, 2), -1), mono(W5, 3))
print("G(-q) = 1 - q^2 G(-q^3) + q^3 G(-q^5):", lhs[:N-8] == rhs[:N-8])
rhs = add(add(one(), W, -1), mono(W3, 1))
print("G(-1/q) = 1 - G(-q) + q G(-q^3):", Wm[:N-8] == rhs[:N-8])
print("W  =", W[:16]); print("W3 =", W3[:16])
