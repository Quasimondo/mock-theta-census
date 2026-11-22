import sys, mpmath as mp
from qmf import *
F = {'psi1*': psi1s, 'psi0*': psi0s, 'kontsevich': kontsevich}[sys.argv[1]]
gam = tuple(int(x) for x in sys.argv[2].split(','))       # a,b,c,d
lo, hi, maxden = float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5])
D = int(sys.argv[6])
a, b, c, d = gam
pts = points(lo, hi, maxden, even_only=True)
# precompute
data = []
for (h, m) in pts:
    p2, q2 = reduce(a * h + b * m, c * h + d * m)
    if q2 % 2: continue
    x = mp.mpf(h) / m; y = mp.mpf(p2) / q2
    data.append((x, y, F(h, m), F(p2 % q2, q2), c * x + d))
data.sort(key=lambda v: v[0])
print(len(data), "points", flush=True)
res = []
for k in [0, mp.mpf(1)/2, 1, mp.mpf(3)/2, 2, mp.mpf(5)/2]:
    for j in range(D):
        al = mp.mpf(j) / D
        P1 = [mp.e ** (TWOPI_I * al * x) * f for x, y, f, g, cx in data]
        P2 = [cx ** (-k) * mp.e ** (TWOPI_I * al * y) * g for x, y, f, g, cx in data]
        D1 = [u - v for u, v in zip(P1[1:], P1)]; D2 = [u - v for u, v in zip(P2[1:], P2)]
        num = sum(u * mp.conj(v) for u, v in zip(D1, D2)); den = sum(abs(v) ** 2 for v in D2)
        eps = num / den
        R = sum(abs(u - eps * v) ** 2 for u, v in zip(D1, D2)) / sum(abs(u) ** 2 for u in D1)
        res.append((R, float(k), j, eps))
res.sort(key=lambda r: r[0])
for R, k, j, eps in res[:6]:
    print(f"  roughness={mp.nstr(R,3):>10}  weight={k}  alpha={j}/{D}  eps={mp.nstr(eps,6)} |eps|={mp.nstr(abs(eps),6)}", flush=True)
print("  median:", mp.nstr(res[len(res)//2][0], 3))
