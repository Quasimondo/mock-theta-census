"""Addendum v24: exact power-series identities used at z = 1, to order N (default 300), integer arithmetic.
  (R1) chi0(q) = 2 F0(q) - phi0(-q)                     (Ramanujan; Watson 1937)
  (R2) q chi1(q) = 2 q F1(q) + phi1(-q)                 (Ramanujan; Watson 1937)
  (R3) psi1(q) - f1(q) = q chi1(q^2)
  (R4) psi0(q) - f0(q) + 2 = chi0(q^2)
  (R5) G(q) H(-q) + G(-q) H(q) = 2/(q^2;q^4)_inf^2      (Ramanujan, one of the forty identities)
  (S1) sum q^(n^2)/(q^4;q^4)_n = G(q)/(-q^2;q^2)_inf          (S2) sum q^(n^2+2n)/(q^4;q^4)_n = H(q)/(-q^2;q^2)_inf
  (S3) sum q^(m^2)/(q;q)_(2m) = G(q^4)/(q;q^2)_inf            (S4) sum q^(m^2+m)/(q;q)_(2m+1) = H(-q)/(q;q^2)_inf
  (S5) sum q^(m^2+m)/(q;q)_(2m) = G(-q)/(q;q^2)_inf           (S6) sum q^(m^2+2m)/(q;q)_(2m+1) = H(q^4)/(q;q^2)_inf
       (Rogers; in Slater's list)
Usage: python3 s1_series.py [N]"""
import sys
N = int(sys.argv[1]) if len(sys.argv) > 1 else 300

def mul(a, b):
    c = [0]*(N + 1)
    for i, x in enumerate(a):
        if x:
            for j in range(0, N + 1 - i):
                if b[j]: c[i + j] += x*b[j]
    return c
def add(a, b, s=1): return [x + s*y for x, y in zip(a, b)]
def mono(k, c=1):
    a = [0]*(N + 1)
    if k <= N: a[k] = c
    return a
def one_minus(k, c=1):          # 1 - c q^k
    a = mono(0); 
    if k <= N: a[k] -= c
    return a
def inv(a):                     # a[0] = +-1
    b = [0]*(N + 1); b[0] = a[0]
    for n in range(1, N + 1):
        b[n] = -a[0]*sum(a[k]*b[n - k] for k in range(1, n + 1))
    return b
def sub2(a, m=2):               # q -> q^m
    b = [0]*(N + 1)
    for i, x in enumerate(a):
        if i*m <= N: b[i*m] = x
    return b
def neg(a): return [x*(-1)**i for i, x in enumerate(a)]      # q -> -q
def shift(a, k): return [0]*k + a[:N + 1 - k]

def series(expo, numer_step=None, denom_step=None, n0=0):
    """sum_{n>=0} q^expo(n) * prod_{j<=n} numer_step(j) / prod_{j<=n} denom_step(j)"""
    s = [0]*(N + 1); num = mono(0); n = 0
    while expo(n) <= N:
        if n >= 1 or n0:
            pass
        s = add(s, shift(num, expo(n))); n += 1
        if numer_step: num = mul(num, numer_step(n))
        if denom_step: num = mul(num, inv(denom_step(n)))
    return s

f0 = series(lambda n: n*n, denom_step=lambda n: one_minus(n, -1))
f1 = series(lambda n: n*n + n, denom_step=lambda n: one_minus(n, -1))
psi0 = series(lambda n: (n + 1)*(n + 2)//2, numer_step=lambda n: one_minus(n, -1))
psi1 = series(lambda n: n*(n + 1)//2, numer_step=lambda n: one_minus(n, -1))
phi0 = series(lambda n: n*n, numer_step=lambda n: one_minus(2*n - 1, -1))
phi1 = series(lambda n: (n + 1)**2, numer_step=lambda n: one_minus(2*n - 1, -1))
F0 = series(lambda n: 2*n*n, denom_step=lambda n: one_minus(2*n - 1))
# F1 = sum q^(2n^2+2n)/(q;q^2)_(n+1)
F1 = [0]*(N + 1); den = inv(one_minus(1)); n = 0
while 2*n*n + 2*n <= N:
    F1 = add(F1, shift(den, 2*n*n + 2*n)); n += 1; den = mul(den, inv(one_minus(2*n + 1)))
# chi0 = sum q^n/(q^(n+1);q)_n,  chi1 = sum q^n/(q^(n+1);q)_(n+1)
chi0 = [0]*(N + 1); chi1 = [0]*(N + 1)
for n in range(N + 1):
    d = mono(0)
    for j in range(n + 1, 2*n + 1): d = mul(d, one_minus(j))
    chi0 = add(chi0, shift(inv(d), n))
    d = mul(d, one_minus(2*n + 1))
    chi1 = add(chi1, shift(inv(d), n))
G = series(lambda n: n*n, denom_step=lambda n: one_minus(n))
H = series(lambda n: n*n + n, denom_step=lambda n: one_minus(n))

def report(name, a, b):
    d = [i for i in range(N + 1) if a[i] != b[i]]
    print(f"{name}: " + (f"holds to order q^{N}" if not d else f"FAILS first at q^{d[0]}"))

report("(R1) chi0(q) = 2 F0(q) - phi0(-q)", chi0, add([2*x for x in F0], neg(phi0), -1))
report("(R2) q chi1(q) = 2 q F1(q) + phi1(-q)", shift(chi1, 1), add(shift([2*x for x in F1], 1), neg(phi1)))
report("(R3) psi1(q) - f1(q) = q chi1(q^2)", add(psi1, f1, -1), shift(sub2(chi1), 1))
report("(R4) psi0(q) - f0(q) + 2 = chi0(q^2)", add(add(psi0, f0, -1), mono(0, 2)), sub2(chi0))
P = mono(0)
for k in range(2, N + 1, 4): P = mul(P, one_minus(k))
P2 = inv(mul(P, P))
report("(R5) G(q) H(-q) + G(-q) H(q) = 2/(q^2;q^4)_inf^2", add(mul(G, neg(H)), mul(neg(G), H)), [2*x for x in P2])

def qpoch(start, step, c=1):     # (c q^start; q^step)_inf
    a = mono(0); k = start
    while k <= N: a = mul(a, one_minus(k, c)); k += step
    return a
def qfac_series(expo, upto):     # sum_m q^expo(m) / (q;q)_(upto(m))
    s = [0]*(N + 1); d = mono(0); j = 0; m = 0
    while expo(m) <= N:
        while j < upto(m): j += 1; d = mul(d, inv(one_minus(j)))
        s = add(s, shift(d, expo(m))); m += 1
    return s
e2 = inv(qpoch(2, 2, -1)); o = inv(qpoch(1, 2))
report("(S1) sum q^(n^2)/(q^4;q^4)_n = G(q)/(-q^2;q^2)_inf", series(lambda n: n*n, denom_step=lambda n: one_minus(4*n)), mul(G, e2))
report("(S2) sum q^(n^2+2n)/(q^4;q^4)_n = H(q)/(-q^2;q^2)_inf", series(lambda n: n*n + 2*n, denom_step=lambda n: one_minus(4*n)), mul(H, e2))
report("(S3) sum q^(m^2)/(q;q)_(2m) = G(q^4)/(q;q^2)_inf", qfac_series(lambda m: m*m, lambda m: 2*m), mul(sub2(G, 4), o))
report("(S4) sum q^(m^2+m)/(q;q)_(2m+1) = H(-q)/(q;q^2)_inf", qfac_series(lambda m: m*m + m, lambda m: 2*m + 1), mul(neg(H), o))
report("(S5) sum q^(m^2+m)/(q;q)_(2m) = G(-q)/(q;q^2)_inf", qfac_series(lambda m: m*m + m, lambda m: 2*m), mul(neg(G), o))
report("(S6) sum q^(m^2+2m)/(q;q)_(2m+1) = H(q^4)/(q;q^2)_inf", qfac_series(lambda m: m*m + 2*m, lambda m: 2*m + 1), mul(sub2(H, 4), o))
