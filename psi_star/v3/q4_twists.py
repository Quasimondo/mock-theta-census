"""Q4: identities for the sign-twisted fifth-order F0, F1.  Exact check to q^(N-1)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qs import *
# F0* = sum (-1)^n q^(2n^2)/(q;q^2)_n ,  F1* = sum (-1)^n q^(2n^2+2n)/(q;q^2)_(n+1)
F0s = eul(lambda n: 2*n*n, lambda P, n: dfac(P, 2*n-1, -1), -1)
F1s = eul(lambda n: 2*n*n+2*n, lambda P, n: dfac(P, 2*n+1, -1), -1, start=dfac(one(), 1, -1))
# negative halves of the bilateral sums:  N0 = sum_{n>=1} q^(n^2) (q;q^2)_n ,  N1 = sum_{n>=0} q^((n+1)^2) (q;q^2)_n
phi0m = eul(lambda n: n*n, lambda P, n: mfac(P, 2*n-1, -1))            # sum_{n>=0} q^(n^2)(q;q^2)_n
N0 = phi0m[:]; N0[0] -= 1
N1 = eul(lambda n: (n+1)**2, lambda P, n: mfac(P, 2*n-1, -1))
# X(q) = sum (-1)^i q^(i^2)/(q;q)_(2i),  Y(q) = sum (-1)^i q^((i+1)^2)/(q;q)_(2i+1)
X = eul(lambda n: n*n, lambda P, n: dfac(dfac(P, 2*n-1, -1), 2*n, -1), -1)
Y = eul(lambda n: (n+1)**2, lambda P, n: dfac(dfac(P, 2*n, -1), 2*n+1, -1), -1, start=dfac(one(), 1, -1))
th = theta4(2); mq = poch_inf(1, 1, +1)                                  # theta4(q^2), (-q;q)_inf
lhs0 = add(F0s, N0); rhs0 = mul(mul(th, mq), sub2(X, 2))
print("(III)  F0* + sum_{n>=1} q^(n^2)(q;q^2)_n = theta4(q^2) (-q;q)_inf X(q^2):", lhs0 == rhs0)
lhs1 = mono(add(F1s, mono(N1, -1), -1), 1)                               # q*(F1* - N1/q)
rhs1 = mul(mul(th, mq), sub2(Y, 2))
print("(IV)   q F1* - sum_{n>=0} q^((n+1)^2)(q;q^2)_n = theta4(q^2) (-q;q)_inf Y(q^2):", lhs1 == rhs1)
# Rogers-type reduction of X, Y to R0, R1 ?
R0 = eul(lambda n: n*n, lambda P, n: dfac(P, n, -1), -1)
R1 = eul(lambda n: n*n+n, lambda P, n: dfac(P, n, -1), -1)
print("X = (-q;q)_inf R0 :", X == mul(mq, R0), "   Y = q (-q;q)_inf R1 :", Y == mono(mul(mq, R1), 1))
for nm, T in (("X", X), ("Y", Y)):
    print(nm, T[:14])
print("(-q)R0", mul(mq, R0)[:14]); print("(-q)R1", mul(mq, R1)[:14])
