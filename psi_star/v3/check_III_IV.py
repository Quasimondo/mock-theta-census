# Independent check of (III), (IV) as integer q-series to q^N
N = 330
def mul(a, b):
    c = [0]*N
    for i, x in enumerate(a):
        if x:
            for j in range(N-i):
                if b[j]: c[i+j] += x*b[j]
    return c
def poch(a, step, sign, n=None, inv=False):   # prod_{k<n} (1 + sign q^{a+k step})^{+-1}
    s = [0]*N; s[0] = 1; k, e = 0, a
    while e < N and (n is None or k < n):
        if inv:
            for i in range(e, N): s[i] -= sign*s[i-e]
        else:
            for i in range(N-1, e-1, -1): s[i] += sign*s[i-e]
        k += 1; e += step
    return s
def sh(s, k): return ([0]*k + s)[:N]
def tot(terms):
    out = [0]*N
    for t in terms:
        out = [x+y for x, y in zip(out, t)]
    return out
def sub_q2(s):   # f(q) -> f(q^2)
    o = [0]*N
    for i, x in enumerate(s):
        if 2*i < N: o[2*i] = x
    return o
F0s = tot([[(-1)**n*x for x in sh(poch(1, 2, -1, n, True), 2*n*n)] for n in range(13) if 2*n*n < N])
F1s = tot([[(-1)**n*x for x in sh(poch(1, 2, -1, n+1, True), 2*n*n+2*n)] for n in range(13) if 2*n*n+2*n < N])
L3 = tot([sh(poch(1, 2, -1, n), n*n) for n in range(1, 19) if n*n < N])
L4 = tot([sh(poch(1, 2, -1, n), (n+1)**2) for n in range(0, 19) if (n+1)**2 < N])
X = tot([[(-1)**n*x for x in sh(poch(1, 1, -1, 2*n, True), n*n)] for n in range(19) if n*n < N])
Y = tot([[(-1)**n*x for x in sh(poch(1, 1, -1, 2*n+1, True), (n+1)**2)] for n in range(19) if (n+1)**2 < N])
th4 = [0]*N
for n in range(-20, 21):
    if n*n < N: th4[n*n] += (-1)**n
pref = mul(sub_q2(th4), poch(1, 1, +1))         # theta4(q^2) (-q;q)_inf
lhs3 = [a+b for a, b in zip(F0s, L3)]; rhs3 = mul(pref, sub_q2(X))
lhs4 = [a-b for a, b in zip(sh(F1s, 1), L4)]; rhs4 = mul(pref, sub_q2(Y))
print('(III) holds to q^%d:' % (N-1), lhs3 == rhs3, ' (IV):', lhs4 == rhs4)
print('lhs3', lhs3[:12]); print('X', X[:12])
