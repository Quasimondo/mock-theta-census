# Independent check of (V), (VI) and the G recursion as exact integer q-series to q^(N-1)
N = 330
def mul(a, b):
    c = [0]*N
    for i, x in enumerate(a):
        if x:
            for j in range(N-i):
                if b[j]: c[i+j] += x*b[j]
    return c
def poch(a, step, sign, n=None, inv=False):    # prod_{k<n} (1 + sign q^{a+k step})^{+-1}
    s = [0]*N; s[0] = 1; k, e = 0, a
    while e < N and (n is None or k < n):
        if inv:
            for i in range(e, N): s[i] -= sign*s[i-e]
        else:
            for i in range(N-1, e-1, -1): s[i] += sign*s[i-e]
        k += 1; e += step
    return s
def sh(s, k): return ([0]*k + s)[:N] if k >= 0 else (s[-k:] + [0]*(-k))
def tot(ts):
    o = [0]*N
    for t in ts: o = [a+b for a, b in zip(o, t)]
    return o
def sc(s, c): return [c*x for x in s]
W  = tot([sc(sh(poch(1, 2, -1, j), j*j+j), (-1)**j) for j in range(20) if j*j+j < N])
W3 = tot([sc(sh(poch(1, 2, -1, j), j*j+3*j), (-1)**j) for j in range(20) if j*j+3*j < N])
L5 = tot([sh(poch(1, 2, -1, n, True), 2*n*n-n) for n in range(1, 14) if 2*n*n-n < N])
L6 = tot([sh(poch(1, 2, -1, n, True), 2*n*n-3*n+1) for n in range(1, 14) if 2*n*n-3*n+1 < N])
Pp = tot([sh(poch(2, 2, -1, k, True), 2*k*k+k) for k in range(14) if 2*k*k+k < N])
Pm = tot([sh(poch(2, 2, -1, k, True), 2*k*k-k) for k in range(14) if 2*k*k-k < N])
psi = [0]*N
for n in range(30):
    if n*(n+1)//2 < N: psi[n*(n+1)//2] = 1
print('(V)  holds to q^%d:' % (N-1), [a+b for a, b in zip(W, L5)] == mul(psi, Pp))
print('(VI) holds to q^%d:' % (N-1), [a+b for a, b in zip(sh(W3, 1), L6)] == mul(psi, Pm))
# G(z) = sum z^j q^{j^2} (q;q^2)_j ; check G(-q^5) = (W - 1 + q^2 W3)/q^3 and G(-1/q) = 1 - W + q W3
G5 = tot([sc(sh(poch(1, 2, -1, j), j*j+5*j), (-1)**j) for j in range(20) if j*j+5*j < N])
rhs = [a-(1 if i == 0 else 0)+b for i, (a, b) in enumerate(zip(W, sh(W3, 2)))]
print('G(-q^5) relation:', sh(G5, 3) == rhs)
Gm = tot([sc(sh(poch(1, 2, -1, j), j*j-j), (-1)**j) for j in range(20) if j*j-j < N])
print('G(-1/q) relation:', Gm == [(1 if i == 0 else 0)-a+b for i, (a, b) in enumerate(zip(W, sh(W3, 1)))])
