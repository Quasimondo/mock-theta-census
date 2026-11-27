# Independent check of addendum-v28 (own code).
# (1) Proposition N, (N1)-(N5) and the Q0 versions, as exact integer power series to q^N.
# (2) Theorem 19 (a), (b): |t| log(|psi1*| / (A |t|^{-1/2})) and |t| log(|f1*| / (2A |t|^{-1/2})) on rays arg x = eps,
#     psi1*, f1* summed directly, against lambda1, lambda2.
import mpmath as mp, sys
N = 200
def mul(a, b):
    c = [0]*(N+1)
    for i, x in enumerate(a):
        if x:
            for j in range(N+1-i): c[i+j] += x*b[j]
    return c
def inv(a):                                   # a[0] = 1
    b = [0]*(N+1); b[0] = 1
    for n in range(1, N+1): b[n] = -sum(a[k]*b[n-k] for k in range(1, n+1))
    return b
def mono(e, c=1):
    s = [0]*(N+1)
    if 0 <= e <= N: s[e] = c
    return s
def add(*xs): return [sum(v) for v in zip(*xs)]
def poch(sign, step, start, n):               # prod_{j=0}^{n-1} (1 + sign q^{start + step j})
    p = mono(0)
    for j in range(n):
        f = mono(0); e = start + step*j
        if e <= N: f[e] += sign
        p = mul(p, f)
    return p
def poch_inf(sign, step, start): return poch(sign, step, start, (N - start)//step + 1 if start <= N else 0)
M = poch_inf(+1, 1, 1)                        # (-q;q)_inf
def series(term, nmax):                       # term(n) -> list
    s = [0]*(N+1)
    for n in range(nmax): s = add(s, term(n))
    return s
f1s = series(lambda n: [(-1)**n*x for x in mul(mono(n*n+n), inv(poch(+1, 1, 1, n)))] if n*n+n <= N else [0]*(N+1), 20)
f0s = series(lambda n: [(-1)**n*x for x in mul(mono(n*n), inv(poch(+1, 1, 1, n)))] if n*n <= N else [0]*(N+1), 20)
psi1s = series(lambda n: [(-1)**n*x for x in mul(mono(n*(n+1)//2), poch(+1, 1, 1, n))] if n*(n+1)//2 <= N else [0]*(N+1), 25)
psi0s = series(lambda n: [(-1)**n*x for x in mul(mono((n+1)*(n+2)//2), poch(+1, 1, 1, n))] if (n+1)*(n+2)//2 <= N else [0]*(N+1), 25)
th4 = add(mono(0), *[mono(n*n, 2*(-1)**n) for n in range(1, 16)])
R1 = series(lambda n: [(-1)**n*x for x in mul(mono(n*n+n), inv(poch(-1, 1, 1, n)))] if n*n+n <= N else [0]*(N+1), 20)
R0 = series(lambda n: [(-1)**n*x for x in mul(mono(n*n), inv(poch(-1, 1, 1, n)))] if n*n <= N else [0]*(N+1), 20)
qk = [inv(poch(-1, 1, 1, k)) for k in range(61)]
def dsum(Q, ch, nrange, krange):
    s = [0]*(N+1)
    for n in nrange:
        for k in krange:
            e = Q(n, k)
            if e > N or e < 0: continue
            s = add(s, [ch(n, k)*x for x in mul(mono(e), qk[k])])
    return s
Q1 = lambda n, k: n*n + n + n*k + k*(k+1)//2
Q0 = lambda n, k: n*n + n*k + k*(k+1)//2
ok = lambda a, b: all(x == y for x, y in zip(a, b))
print('(1) Proposition N, exact to q^%d' % N)
print('  (N1)  M f1* :', ok(mul(M, f1s), dsum(Q1, lambda n, k: (-1)**n, range(15), range(21))))
print('  (N1\') M f0* :', ok(mul(M, f0s), dsum(Q0, lambda n, k: (-1)**n, range(15), range(21))))
print('  (N2)  M th4 R1:', ok(mul(M, mul(th4, R1)), dsum(Q1, lambda n, k: (-1)**(n+k), range(15), range(21))))
print('  (N2\') M th4 R0:', ok(mul(M, mul(th4, R0)), dsum(Q0, lambda n, k: (-1)**(n+k), range(15), range(21))))
print('  (N3)  M psi1* = k even part:', ok(mul(M, psi1s), dsum(Q1, lambda n, k: (-1)**n, range(15), range(0, 21, 2))))
print('  (N3\') M psi0* = k odd part (Q0):', ok(mul(M, psi0s), dsum(Q0, lambda n, k: (-1)**n, range(15), range(1, 21, 2))))
rhs4 = [-x for x in mul(th4, series(lambda m: [(-1)**m*x for x in mul(mono(m*m+m), qk[2*m+1])] if m*m+m <= N and 2*m+1 < 61 else [0]*(N+1), 15))]
rhs4b = mul(th4, series(lambda m: [(-1)**m*x for x in mul(mono(m*m+m), qk[2*m])] if m*m+m <= N and 2*m < 61 else [0]*(N+1), 15))
print('  (N4)  n over Z, Q1:', ok(dsum(Q1, lambda n, k: (-1)**n, range(-90, 15), range(61)), rhs4))
print('  (N4\') n over Z, Q0:', ok(dsum(Q0, lambda n, k: (-1)**n, range(-90, 15), range(61)), rhs4b))
q2inf = poch_inf(-1, 2, 2); q2k = [inv(poch(-1, 2, 2, j)) for j in range(101)]
s5 = [0]*(N+1)
for n in range(25):
    for j in range(101):
        e = n*(n+1)//2 + 2*n*j + 2*j
        if e > N: continue
        s5 = add(s5, [(-1)**n*x for x in mul(mono(e), mul(qk[n] if n < 25 else [0]*(N+1), q2k[j]))])
print('  (N5)  psi1* = (q^2;q^2)_inf sum ...:', ok(psi1s, mul(q2inf, s5)))
if len(sys.argv) > 1 and sys.argv[1] == 'series': sys.exit()
print('(2) Theorem 19 on rays arg x = eps, |t| = 0.05')
A = mp.sqrt(mp.pi/mp.sqrt(3)); Cl2 = mp.clsin(2, mp.pi/3)
for eps in (1, 2, 5, 10):
    with mp.workdps(60):
        e = mp.radians(eps); T = mp.mpf('0.05'); t = T*mp.expj(e - mp.pi/2); q = mp.exp(-t)
        s1 = s2 = 0; P = mp.mpc(1); n = 0
        while True:
            a = (-1)**n*q**(n*(n+1)//2)*P; b = (-1)**n*q**(n*n+n)/P
            s1 += a; s2 += b
            if n > 20 and abs(a) < mp.mpf(10)**-50 and abs(b) < mp.mpf(10)**-50: break
            n += 1; P *= 1 + q**n
        l1 = Cl2*mp.cos(e) - mp.pi**2/3*mp.sin(e); l2 = Cl2*mp.cos(e) - 7*mp.pi**2/3*mp.sin(e)
        m1 = T*mp.log(abs(s1)/(A*T**-0.5)); m2 = T*mp.log(abs(s2)/(2*A*T**-0.5))
        print(f'  eps = {eps:2d}: psi1* {mp.nstr(m1, 6)} vs lambda1 {mp.nstr(l1, 6)};   f1* ' + (f'{mp.nstr(m2, 6)} vs lambda2 {mp.nstr(l2, 6)}' if l2 > 0 else f'|f1*| = {mp.nstr(abs(s2), 5)} (lambda2 < 0)'))
