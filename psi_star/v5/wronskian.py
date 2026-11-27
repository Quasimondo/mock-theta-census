"""Exact q-series (to q^N) of the W family and of the Wronskian-type combination  Wr = q W3 Fa - W Fb."""
N = 40
def mul(a, b):
    c = [0]*N
    for i, x in enumerate(a):
        if x:
            for j in range(N - i): c[i + j] += x*b[j]
    return c
def inv(a):          # 1/a, a[0] = 1
    b = [0]*N; b[0] = 1
    for n in range(1, N): b[n] = -sum(a[k]*b[n - k] for k in range(1, n + 1))
    return b
def mono(e): 
    m = [0]*N
    if e < N: m[e] = 1
    return m
def fam():
    W = [0]*N; W3 = [0]*N; Fa = [0]*N; Fb = [0]*N
    poch = mono(0)               # (q;q^2)_j
    j = 0
    while j*j + j < N or (j >= 1 and 2*j*j - 3*j + 1 < N):
        s = -1 if j & 1 else 1
        for i in range(N):
            if poch[i]:
                if i + j*j + j < N: W[i + j*j + j] += s*poch[i]
                if i + j*j + 3*j < N: W3[i + j*j + 3*j] += s*poch[i]
        if j >= 1:
            ip = inv(poch)
            for i in range(N):
                if ip[i]:
                    if i + 2*j*j - j < N: Fa[i + 2*j*j - j] += ip[i]
                    if i + 2*j*j - 3*j + 1 < N: Fb[i + 2*j*j - 3*j + 1] += ip[i]
        f = mono(0); 
        if 2*j + 1 < N: f[2*j + 1] = -1
        poch = mul(poch, f); j += 1
    return W, W3, Fa, Fb
if __name__ == "__main__":
    W, W3, Fa, Fb = fam()
    qW3 = [0] + W3[:-1]
    Wr = [x - y for x, y in zip(mul(qW3, Fa), mul(W, Fb))]
    print("W   ", W[:16]); print("W3  ", W3[:16]); print("Fa  ", Fa[:16]); print("Fb  ", Fb[:16])
    print("Wr = q W3 Fa - W Fb:", Wr[:24])
    print("2 Wr:", [2*x for x in Wr[:20]])
