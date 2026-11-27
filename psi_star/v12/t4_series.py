"""Test 4 (addendum v12): the three q-series identities used in (B3), as formal power series up to p^N (default N = 300).
   (i)   sum_k (-1)^k p^(k^2-k) (p;p^2)_k        = 1 - W + p W3          [ = Wt(p^-2) ]
   (ii)  sum_k p^(2k^2+k) / (p;p^2)_k            = F_b - F_a             [ = G0(1) ]
   (iii) p * sum_k p^(2k^2+3k) / (p;p^2)_(k+1)   = F_a                   [ = p G1(1) ]
   W = sum (-1)^j p^(j^2+j) (p;p^2)_j,  W3 = sum (-1)^j p^(j^2+3j) (p;p^2)_j,
   F_a = sum_{n>=1} p^(2n^2-n)/(p;p^2)_n,  F_b = sum_{n>=1} p^(2n^2-3n+1)/(p;p^2)_n.
Exact integer arithmetic.  Usage: python3 t4_series.py [N]"""
import sys
N = int(sys.argv[1]) if len(sys.argv) > 1 else 300

def mul(a, b):
    c = [0]*(N + 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b[:N + 1 - i]):
                if bj: c[i + j] += ai*bj
    return c

def inv(a):                      # a[0] = 1
    c = [0]*(N + 1); c[0] = 1
    for n in range(1, N + 1):
        c[n] = -sum(a[k]*c[n - k] for k in range(1, n + 1))
    return c

def shift(a, e, s=1):            # s * p^e * a
    c = [0]*(N + 1)
    if e <= N:
        for i in range(N + 1 - e): c[i + e] = s*a[i]
    return c

def add(*xs):
    return [sum(v) for v in zip(*xs)]

one = [1] + [0]*N
poch = [one]                     # poch[k] = (p;p^2)_k
k = 0
while 2*k + 1 <= N + 2 and k < 40:
    f = [0]*(N + 1); f[0] = 1
    if 2*k + 1 <= N: f[2*k + 1] = -1
    poch.append(mul(poch[-1], f)); k += 1
ipoch = [inv(a) for a in poch]
K = len(poch) - 1
zero = [0]*(N + 1)
W = add(zero, *[shift(poch[j], j*j + j, (-1)**j) for j in range(K) if j*j + j <= N])
W3 = add(zero, *[shift(poch[j], j*j + 3*j, (-1)**j) for j in range(K) if j*j + 3*j <= N])
Fa = add(zero, *[shift(ipoch[n], 2*n*n - n) for n in range(1, K) if 2*n*n - n <= N])
Fb = add(zero, *[shift(ipoch[n], 2*n*n - 3*n + 1) for n in range(1, K) if 2*n*n - 3*n + 1 <= N])
lhs1 = add(zero, *[shift(poch[k], k*k - k, (-1)**k) for k in range(K) if k*k - k <= N])
rhs1 = add(one, [-v for v in W], shift(W3, 1))
lhs2 = add(zero, *[shift(ipoch[k], 2*k*k + k) for k in range(K) if 2*k*k + k <= N])
rhs2 = add(Fb, [-v for v in Fa])
lhs3 = add(zero, *[shift(ipoch[k + 1], 2*k*k + 3*k + 1) for k in range(K - 1) if 2*k*k + 3*k + 1 <= N])
print(f"formal power series up to p^{N}:")
print("   (i)   Wt(p^-2) = 1 - W + p W3 :", lhs1 == rhs1, "   first coefficients:", lhs1[:12])
print("   (ii)  G0(1) = F_b - F_a       :", lhs2 == rhs2, "   first coefficients:", lhs2[:12])
print("   (iii) p G1(1) = F_a           :", lhs3 == Fa, "   first coefficients:", Fa[:12])
