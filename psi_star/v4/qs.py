"""Exact truncated integer q-series (shared helpers)."""
N = 320
def zero(): return [0]*N
def one():
    a = [0]*N; a[0] = 1; return a
def mul(a, b):
    r = [0]*N
    for i, x in enumerate(a):
        if x:
            for j in range(N-i):
                y = b[j]
                if y: r[i+j] += x*y
    return r
def mono(a, e, c=1):
    r = [0]*N
    if e < 0:
        assert all(x == 0 for x in a[:-e]), "negative shift would drop terms"
    for i, x in enumerate(a):
        if x and 0 <= i+e < N: r[i+e] += c*x
    return r
def add(a, b, c=1): return [x + c*y for x, y in zip(a, b)]
def mfac(a, k, s):          # a * (1 + s q^k)
    r = a[:]
    for i in range(N-1, k-1, -1): r[i] += s*a[i-k]
    return r
def dfac(a, k, s):          # a / (1 + s q^k)
    r = a[:]
    for i in range(k, N): r[i] -= s*r[i-k]
    return r
def poch_inf(a, step, s, inv=False):     # (−s q^a; q^step)_inf, i.e. prod (1 + s q^(a+step*j))
    r = one(); e = a
    while e < N:
        r = dfac(r, e, s) if inv else mfac(r, e, s); e += step
    return r
def eul(expo, fac_step, sign=1, nmax=10**6, start=None):
    """sum_n sign^n q^expo(n) P_n, with P_0 = start (default 1) and P_n = fac_step(P_{n-1}, n)."""
    out = zero(); P = start if start is not None else one(); n = 0
    while n < nmax:
        if n: P = fac_step(P, n)
        e = expo(n)
        if e >= N: 
            if n > 3: break
        else:
            out = add(out, mono(P, e), sign**n)
        n += 1
    return out
def theta4(k=1):            # sum (-1)^n q^(k n^2)
    a = zero(); n = 0
    while k*n*n < N:
        a[k*n*n] += (1 if n == 0 else 2)*(-1)**n; n += 1
    return a
def sub2(a, k):             # a(q^k)
    r = zero()
    for i, x in enumerate(a):
        if i*k < N: r[i*k] = x
    return r
