"""Checks for  D = theta4 * f + 4R  with the closed form of R.   Run:  python3 verify_R.py
Needs only the standard library, plus mpmath for the last (numerical) part.

  D = sum_{m,n>=0} q^(m^2-mn+n^2) / ((-q)_m (-q)_n)
  S_d = sum_{n>=0} q^(n^2+dn) / ((-q)_n (-q)_(n+d))          (the diagonal m - n = d of D; S_0 = f)
  Lemma:  S_d + S_(d+1) = 2/(-q)_d
  R = sum_{d>k>=0} (-1)^(d+k+1) q^(d^2) / (-q)_k
"""
from fractions import Fraction as Fr
N = 320   # exact checks hold to q^(N-1)

def mul(a, b):
    r = [0] * N
    for i, x in enumerate(a):
        if x:
            for j in range(N - i):
                if b[j]:
                    r[i + j] += x * b[j]
    return r
def mono(a, e, c=1):
    r = [0] * N
    for i, x in enumerate(a):
        if x and i + e < N:
            r[i + e] += c * x
    return r
def add(a, b, c=1):
    return [x + c * y for x, y in zip(a, b)]
def div1p(a, k):            # a / (1 + q^k)
    r = a[:]
    for i in range(k, N):
        r[i] -= r[i - k]
    return r

INV = [[1] + [0] * (N - 1)]                 # INV[n] = 1/(-q;q)_n
for n in range(1, 64):
    INV.append(div1p(INV[-1], n))

def S(d):
    r, n = [0] * N, 0
    while n * n + d * n < N:
        r = add(r, mono(mul(INV[n], INV[n + d]), n * n + d * n)); n += 1
    return r

# D summed over ALL (m,n) with m^2-mn+n^2 < N  (the form is >= 3/4 max(m,n)^2, so m,n < 40 is enough for N = 320)
D = [0] * N
for m in range(40):
    for n in range(40):
        e = m * m - m * n + n * n
        if e < N:
            D = add(D, mono(mul(INV[m], INV[n]), e))
th4 = [0] * N
n = 0
while n * n < N:
    th4[n * n] += (1 if n == 0 else 2) * (-1) ** n; n += 1
f = S(0)
R4 = add(D, mul(th4, f), -1)
assert all(x % 4 == 0 for x in R4)
R = [x // 4 for x in R4]                      # R as defined in the brief, computed independently of the formula
print("R, q^0..q^20:", R[:21])

print("1. Lemma S_d + S_(d+1) = 2/(-q)_d for d = 0..15:",
      all(add(S(d), S(d + 1)) == [2 * x for x in INV[d]] for d in range(16)))

Rf, d = [0] * N, 1
while d * d < N:
    P = [0] * N
    for k in range(d):
        P = add(P, INV[k], (-1) ** k)
    Rf = add(Rf, mono(P, d * d, (-1) ** (d - 1))); d += 1
print(f"2. R = sum_(d>k>=0) (-1)^(d+k+1) q^(d^2)/(-q)_k, exact to q^{N-1}:", Rf == R)

Rg, i = [0] * N, 1
while i * i < N:
    G, n = [0] * N, 0
    while n * n + 2 * i * n < N:
        G = add(G, mono(INV[n], n * n + 2 * i * n)); n += 1
    Rg = add(Rg, mono(G, i * i, (-1) ** (i - 1))); i += 1
print(f"3. R = sum_(i>=1) (-1)^(i-1) q^(i^2) Phi(q^(2i)), Phi(z) = sum z^n q^(n^2)/(-q)_n, exact to q^{N-1}:", Rg == R)

# 4. general z:  D_z = sum q^(m^2-mn+n^2)/((zq)_m (q/z)_n)
#    = theta(z) * rank(z) - (1-1/z) sum_{d>k>=0} z^(d-k) q^(d^2)/(zq)_k - (1-z) sum_{d>k>=0} z^(k-d) q^(d^2)/(q/z)_k
M = 28
def fmul(a, b):
    r = [Fr(0)] * M
    for i, x in enumerate(a):
        if x:
            for j in range(M - i):
                r[i + j] += x * b[j]
    return r
def fmono(a, e, c=1):
    r = [Fr(0)] * M
    for i, x in enumerate(a):
        if i + e < M:
            r[i + e] += c * x
    return r
def fadd(a, b, c=1):
    return [x + c * y for x, y in zip(a, b)]
def pinv(c, n):             # 1/(c q;q)_n
    r = [Fr(1)] + [Fr(0)] * (M - 1)
    for k in range(1, n + 1):
        for i in range(k, M):
            r[i] += c * r[i - k]
    return r
for z in (Fr(2), Fr(-3), Fr(5, 7)):
    A = [pinv(z, n) for n in range(M)]
    B = [pinv(1 / z, n) for n in range(M)]
    Dz = [Fr(0)] * M
    for m in range(12):
        for n in range(12):
            e = m * m - m * n + n * n
            if e < M:
                Dz = fadd(Dz, fmono(fmul(A[m], B[n]), e))
    rank, n = [Fr(0)] * M, 0
    while n * n < M:
        rank = fadd(rank, fmono(fmul(A[n], B[n]), n * n)); n += 1
    th = [Fr(0)] * M; th[0] = Fr(1); d = 1
    while d * d < M:
        th[d * d] += z ** d + z ** (-d); d += 1
    rhs, d = fmul(rank, th), 1
    while d * d < M:
        for k in range(d):
            rhs = fadd(rhs, fmono(A[k], d * d), -(1 - 1 / z) * z ** (d - k))
            rhs = fadd(rhs, fmono(B[k], d * d), -(1 - z) * z ** (k - d))
        d += 1
    print(f"4. general z = {z}: D_z = theta(z)*rank(z) - tails, exact to q^{M-1}:", rhs == Dz)

# 5. numerical: growth of R along q = -e^(-t)
try:
    import mpmath as mp
except ImportError:
    mp = None
if mp:
    mp.mp.dps = 50
    def Rnum(t):
        q = -mp.exp(-t)
        terms = [(-1) ** (d - 1) * q ** (d * d) for d in range(1, int(3 / t) + 301)]
        T = [0] * (len(terms) + 1)
        for k in range(len(terms) - 1, -1, -1):
            T[k] = T[k + 1] + terms[k]
        s, den, qk = 0, mp.mpf(1), mp.mpf(1)
        for k in range(len(terms)):
            if k:
                qk *= q; den *= 1 + qk
            s += (-1) ** k * T[k] / den
        return s
    C = mp.sqrt(2 / (5 - mp.sqrt(5)))
    ts = [mp.mpf(1) / (100 * 2 ** i) for i in range(7)]
    a = [Rnum(t) / (C * mp.jtheta(3, 0, mp.exp(-t)) * mp.exp(mp.pi ** 2 / (60 * t))) for t in ts]
    def rich(v):
        tab = [v[:]]
        for j in range(1, len(v)):
            p = tab[-1]
            tab.append([(2 ** j * p[i + 1] - p[i]) / (2 ** j - 1) for i in range(len(p) - 1)])
        return tab[-1][0]
    print("5. R(-e^-t) / ( sqrt(2/(5-sqrt5)) * theta3(e^-t) * e^(pi^2/60t) ) -> ", mp.nstr(rich(a), 12),
          "   first correction sigma_1 =", mp.nstr(rich([(a[i] - 1) / ts[i] for i in range(len(a))]), 11))
