"""Exact checks to q^309 (integer power series).
 (a) (-q;q)_inf R0 = N0,  (-q;q)_inf R1 = N1         (R0, R1: sign-twisted Rogers-Ramanujan sums)
 (b) (q;q^2)_inf P+ = E+,  (q;q^2)_inf P- = E-       (P+-: half-shifted sums; E+- as in hfam.py, nome q)
 (c) E+' E- + E-' E+ = 2                             (Casoratian on the dual side)
 (d) X N0 + N1 Y = 1                                 (Casoratian on the q side)"""
N = 310
def mul(a, b):
    c = [0]*N
    for i, x in enumerate(a):
        if x:
            for j in range(N - i):
                if b[j]: c[i + j] += x*b[j]
    return c
def inv(a):
    b = [0]*N; b[0] = 1
    for n in range(1, N): b[n] = -sum(a[k]*b[n - k] for k in range(1, n + 1) if a[k])
    return b
def one(): return [1] + [0]*(N - 1)
def times1m(a, e, sign=-1):          # a * (1 + sign q^e)
    c = a[:]
    for i in range(N - e): c[i + e] += sign*a[i]
    return c
def shift(a, e): return [0]*e + a[:N - e] if e < N else [0]*N
def add(acc, a, s=1):
    for i in range(N): acc[i] += s*a[i]

# 1/(q;q)_k for all needed k, by repeated division
invq = [one()]
k = 0
while k < 40:
    k += 1
    d = one(); d[k] = -1 if k < N else 0
    invq.append(mul(invq[-1], inv(d)))
mq = one(); oq2 = one()                 # (-q;q)_inf and (q;q^2)_inf
for i in range(1, N): mq = times1m(mq, i, +1)
for i in range(1, N, 2): oq2 = times1m(oq2, i, -1)
R0 = [0]*N; R1 = [0]*N; N0 = [0]*N; N1 = [0]*N; X = [0]*N; Y = [0]*N
for n in range(18):
    s = -1 if n & 1 else 1
    add(R0, shift(invq[n], n*n), s); add(R1, shift(invq[n], n*n + n), s)
    add(N0, shift(invq[2*n], n*n + n), s); add(N1, shift(invq[2*n + 1], n*n + n), s)
    add(X, shift(invq[2*n], n*n), s); add(Y, shift(invq[2*n + 1], (n + 1)**2), s)
# base q^2 denominators 1/(q^2;q^2)_m
invq2 = [one()]
for m in range(1, 30):
    d = one()
    if 2*m < N: d[2*m] = -1
    invq2.append(mul(invq2[-1], inv(d)))
Pp = [0]*N; Pm = [0]*N; Ep = [0]*N; Em = [0]*N; Epp = [0]*N; Emp = [0]*N
for m in range(30):
    if 2*m*m + m < N: add(Pp, shift(invq2[m], 2*m*m + m))
    if m == 0: add(Pm, invq2[0])
    elif 2*m*m - m < N: add(Pm, shift(invq2[m], 2*m*m - m))
    s = -1 if m & 1 else 1
    a = shift(invq2[m], m*(m + 1)//2); b = shift(invq2[m], m*(m + 3)//2)
    add(Ep, a, s); add(Em, b, s); add(Epp, a); add(Emp, b)
print("(a)", mul(mq, R0) == N0, mul(mq, R1) == N1)
print("(b)", mul(oq2, Pp) == Ep, mul(oq2, Pm) == Em)
c = mul(Epp, Em); add(c, mul(Emp, Ep)); print("(c)", c == [2] + [0]*(N - 1))
d = mul(X, N0); add(d, mul(N1, Y)); print("(d)", d == one())
print("N0 =", N0[:12]); print("N1 =", N1[:12]); print("E+ =", Ep[:12]); print("E+' =", Epp[:12])
