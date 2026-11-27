"""Addendum v18, task 3: exact power series checks to order N (default 300).
 (1) W(1/p) = phi0(-p),  p W(p) = -phi1(-p),  G0(p) = F0(p),  G1(p) = F1(p)       (fifth-order mock theta functions)
       W(t) = sum (-t)^k p^(k^2+k) (p;p^2)_k,  G0(t) = sum t^-k p^(2k^2+k)/(p;p^2)_k,  G1(t) = sum t^-k p^(2k^2+3k)/(p;p^2)_(k+1)
       phi0(q) = sum q^(n^2) (-q;q^2)_n, phi1(q) = sum q^((n+1)^2) (-q;q^2)_n, F0(q) = sum q^(2n^2)/(q;q^2)_n, F1(q) = sum q^(2n^2+2n)/(q;q^2)_(n+1)
 (2) the four Rogers-Ramanujan type evaluations behind I(0), I(1/b)  (Slater's list):
       sum q^(m^2)/(q;q)_2m = G(q^4)/(q;q^2)_inf,       sum q^(m^2+m)/(q;q)_(2m+1) = H(-q)/(q;q^2)_inf,
       sum r^(n^2)/(r^4;r^4)_n = G(r)/(-r^2;r^2)_inf,   sum r^(n^2+2n)/(r^4;r^4)_n = H(r)/(-r^2;r^2)_inf
Usage: python3 c1_series.py [N]"""
import sys
N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
M = N + 1

def mul(a, b):
    c = [0]*M
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b[:M - i]):
                if y: c[i + j] += x*y
    return c
def mono(e, c=1):
    v = [0]*M
    if e < M: v[e] = c
    return v
def add(a, b, s=1): return [x + s*y for x, y in zip(a, b)]
def inv(a):
    b = [0]*M; b[0] = 1
    for k in range(1, M): b[k] = -sum(a[i]*b[k - i] for i in range(1, k + 1))
    return b
def shift(a, e): return ([0]*e + a)[:M]
def one_minus(e, c=1):          # 1 - c q^e
    v = [0]*M; v[0] = 1
    if e < M: v[e] -= c
    return v
def prod_inf(f):                # product over j >= 0 of f(j) (lists), until f returns None
    out = mono(0); j = 0
    while True:
        g = f(j)
        if g is None: return out
        out = mul(out, g); j += 1
def subs(a, k, sign=1):         # a(sign * q^k)
    v = [0]*M
    for i, x in enumerate(a):
        if i*k < M: v[i*k] = x*(sign**i)
    return v
def report(name, a, b):
    d = [i for i in range(M) if a[i] != b[i]]
    print(f"    {name}: " + ("exact to order %d" % N if not d else "FAILS at order %d" % d[0]))

zero = [0]*M
# (1)
phi0 = zero[:]; phi1 = zero[:]; po = mono(0); n = 0          # po = (-q;q^2)_n
while n*n < M:
    phi0 = add(phi0, shift(po, n*n)); phi1 = add(phi1, shift(po, (n + 1)**2))
    po = mul(po, one_minus(2*n + 1, -1)); n += 1
W_inv = zero[:]; pW = zero[:]; po = mono(0); k = 0           # po = (p;p^2)_k
while k*k < M:
    W_inv = add(W_inv, shift(po, k*k), (-1)**k)              # W(1/p) = sum (-1)^k p^(k^2) (p;p^2)_k
    pW = add(pW, shift(po, k*k + 2*k + 1), (-1)**k)          # p W(p) = sum (-1)^k p^(k^2+2k+1) (p;p^2)_k
    po = mul(po, one_minus(2*k + 1)); k += 1
print("(1) dual functions at z = 1")
report("W(1/p) = phi0(-p)", W_inv, subs(phi0, 1, -1))
report("p W(p) = -phi1(-p)", pW, [-x for x in subs(phi1, 1, -1)])
F0 = zero[:]; F1 = zero[:]; ipo = mono(0); k = 0             # ipo = 1/(q;q^2)_k
while 2*k*k < M:
    F0 = add(F0, shift(ipo, 2*k*k)); ipo = mul(ipo, inv(one_minus(2*k + 1))); F1 = add(F1, shift(ipo, 2*k*k + 2*k)); k += 1
print("    G0(p) = F0(p), G1(p) = F1(p): the same series term by term (t = p: p^(-k) p^(2k^2+k) = p^(2k^2), p^(-k) p^(2k^2+3k) = p^(2k^2+2k))")
# (2)
G = inv(prod_inf(lambda j: mul(one_minus(5*j + 1), one_minus(5*j + 4)) if 5*j + 1 < M else None))
H = inv(prod_inf(lambda j: mul(one_minus(5*j + 2), one_minus(5*j + 3)) if 5*j + 2 < M else None))
qq2 = prod_inf(lambda j: one_minus(2*j + 1) if 2*j + 1 < M else None)          # (q;q^2)_inf
mq2 = prod_inf(lambda j: one_minus(2*j + 2, -1) if 2*j + 2 < M else None)      # (-q^2;q^2)_inf
Xe = zero[:]; Xo = zero[:]; iq = mono(0); m = 0                                # iq = 1/(q;q)_2m
while m*m < M:
    Xe = add(Xe, shift(iq, m*m)); iq = mul(iq, inv(one_minus(2*m + 1))); Xo = add(Xo, shift(iq, m*m + m)); iq = mul(iq, inv(one_minus(2*m + 2))); m += 1
print("(2) Rogers-Ramanujan type evaluations")
report("sum q^(m^2)/(q)_2m = G(q^4)/(q;q^2)_inf", mul(Xe, qq2), subs(G, 4))
report("sum q^(m^2+m)/(q)_(2m+1) = H(-q)/(q;q^2)_inf", mul(Xo, qq2), subs(H, 1, -1))
S0 = zero[:]; S1 = zero[:]; i4 = mono(0); n = 0                                # i4 = 1/(r^4;r^4)_n
while n*n < M:
    S0 = add(S0, shift(i4, n*n)); S1 = add(S1, shift(i4, n*n + 2*n)); n += 1; i4 = mul(i4, inv(one_minus(4*n)))
report("sum r^(n^2)/(r^4;r^4)_n = G(r)/(-r^2;r^2)_inf", mul(S0, mq2), G)
report("sum r^(n^2+2n)/(r^4;r^4)_n = H(r)/(-r^2;r^2)_inf", mul(S1, mq2), H)
