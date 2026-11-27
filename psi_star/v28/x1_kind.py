"""Addendum v28: what kind of object is f*?
 a  Exact series identities to order N (default 300), integer arithmetic.  M = (-q;q)_inf, th4 = (q;q)_inf/(-q;q)_inf,
    Q1(n,k) = n^2 + n + nk + k(k+1)/2,  Q0(n,k) = n^2 + nk + k(k+1)/2:
      (N1)  M f1*      = sum_{n,k>=0} (-1)^n     q^Q1/(q;q)_k        (N1') M f0*      = the same with Q0
      (N2)  M th4 R1   = sum_{n,k>=0} (-1)^(n+k) q^Q1/(q;q)_k        (N2') M th4 R0   = the same with Q0
      (N3)  M psi1*    = sum over k even of (N1)                      (N3') M psi0*    = sum over k odd of (N1')
      (N4)  sum_{n in Z, k>=0} (-1)^n q^Q1/(q;q)_k = -th4 sum_m (-1)^m q^(m^2+m)/(q;q)_(2m+1)
      (N4') sum_{n in Z, k>=0} (-1)^n q^Q0/(q;q)_k =  th4 sum_m (-1)^m q^(m^2+m)/(q;q)_(2m)
      (N5)  psi1* = (q^2;q^2)_inf sum_{n,j>=0} (-1)^n q^(n(n+1)/2 + 2nj + 2j)/((q;q)_n (q^2;q^2)_j)
 b  Growth at the cusp x -> 0 along rays arg x = eps:  |t| log(|f|/(c |t|^(-1/2)))  against
      lambda1 = Cl2 cos(eps) - (pi^2/3) sin(eps)   (psi1*, c = A)     lambda2 = Cl2 cos(eps) - (7 pi^2/3) sin(eps)   (f1*, c = 2A),
    and for Ramanujan's psi1, f1 for comparison:  |t| log|psi1| against (pi^2/60) sin(eps), |f1| bounded.
 c  The q-difference equations in z:  L y = y(z) - z q y(zq) - z q^2 y(z q^2);  L F = 1,  L [Fcheck(1/z)/z] = -2.
Usage: python3 x1_kind.py [a b c] [N]"""
import sys

def part_a(N):
    def mul(a, b):
        c = [0]*(N + 1)
        for i, x in enumerate(a):
            if x:
                for j in range(N + 1 - i):
                    if b[j]: c[i + j] += x*b[j]
        return c
    def add(a, b, s=1): return [x + s*y for x, y in zip(a, b)]
    def mono(k, c=1):
        a = [0]*(N + 1)
        if 0 <= k <= N: a[k] = c
        return a
    def omq(k, c=1):                      # 1 - c q^k
        a = mono(0)
        if k <= N: a[k] -= c
        return a
    def inv(a):
        b = [0]*(N + 1); b[0] = a[0]
        for n in range(1, N + 1): b[n] = -a[0]*sum(a[k]*b[n - k] for k in range(1, n + 1))
        return b
    def shift(a, k): return ([0]*k + a[:N + 1 - k]) if k <= N else [0]*(N + 1)
    def report(name, a, b):
        d = [i for i in range(N + 1) if a[i] != b[i]]
        print(f"    {name}: " + (f"holds to order q^{N}" if not d else f"FAILS first at q^{d[0]}"))
    # single sums
    def onesided(expo, sign, step, start=None):
        s = [0]*(N + 1); cur = start or mono(0); n = 0
        while expo(n) <= N:
            s = add(s, shift(cur, expo(n)), sign(n)); n += 1; cur = mul(cur, step(n))
        return s
    alt = lambda n: (-1)**n
    f1 = onesided(lambda n: n*n + n, alt, lambda n: inv(omq(n, -1)))
    f0 = onesided(lambda n: n*n, alt, lambda n: inv(omq(n, -1)))
    psi1 = onesided(lambda n: n*(n + 1)//2, alt, lambda n: omq(n, -1))
    psi0 = onesided(lambda n: (n + 1)*(n + 2)//2, alt, lambda n: omq(n, -1))
    R1 = onesided(lambda n: n*n + n, alt, lambda n: inv(omq(n)))
    R0 = onesided(lambda n: n*n, alt, lambda n: inv(omq(n)))
    M = mono(0); Pq = mono(0); P2 = mono(0)
    for k in range(1, N + 1):
        M = mul(M, omq(k, -1)); Pq = mul(Pq, omq(k))
        if k % 2 == 0: P2 = mul(P2, omq(k))
    th4 = mul(Pq, inv(M))
    # double sums: for each k, q^(k(k+1)/2)/(q)_k times the n-sum
    def dsum(delta, chi, bilateral=False):
        tot = [0]*(N + 1); ik = mono(0); k = 0
        while k*(k + 1)//2 - ((k + delta)**2)//4 - 1 <= N:
            inner = [0]*(N + 1); rng = range(-2*N, 2*N) if bilateral else range(0, 2*N)
            for n in rng:
                e = n*n + (delta + k)*n + k*(k + 1)//2
                if 0 <= e <= N: inner[e] += (-1)**(n % 2)
            if any(inner): tot = add(tot, mul(inner, ik), chi(k))
            k += 1; ik = mul(ik, inv(omq(k)))
            if k > 2*N: break
        return tot
    one = lambda k: 1; par = lambda k: (-1)**k; ev = lambda k: 1 - k % 2; od = lambda k: k % 2
    report("(N1)  M f1* = sum (-1)^n q^Q1/(q)_k", mul(M, f1), dsum(1, one))
    report("(N2)  M th4 R1 = sum (-1)^(n+k) q^Q1/(q)_k", mul(M, mul(th4, R1)), dsum(1, par))
    report("(N3)  M psi1* = the part k even", mul(M, psi1), dsum(1, ev))
    report("(N1') M f0* = sum (-1)^n q^Q0/(q)_k", mul(M, f0), dsum(0, one))
    report("(N2') M th4 R0 = sum (-1)^(n+k) q^Q0/(q)_k", mul(M, mul(th4, R0)), dsum(0, par))
    report("(N3') M psi0* = the part k odd", mul(M, psi0), dsum(0, od))
    def qfac(expo, upto):
        s = [0]*(N + 1); d = mono(0); j = 0; m = 0
        while expo(m) <= N:
            while j < upto(m): j += 1; d = mul(d, inv(omq(j)))
            s = add(s, shift(d, expo(m)), (-1)**m); m += 1
        return s
    report("(N4)  n over Z: = -th4 sum (-1)^m q^(m^2+m)/(q)_(2m+1)", dsum(1, one, True), [-x for x in mul(th4, qfac(lambda m: m*m + m, lambda m: 2*m + 1))])
    report("(N4') n over Z: =  th4 sum (-1)^m q^(m^2+m)/(q)_(2m)", dsum(0, one, True), mul(th4, qfac(lambda m: m*m + m, lambda m: 2*m)))
    tot = [0]*(N + 1); i2 = mono(0); j = 0
    while 2*j <= N:
        inner = [0]*(N + 1); iq = mono(0); n = 0
        while n*(n + 1)//2 + 2*n*j + 2*j <= N:
            inner = add(inner, shift(iq, n*(n + 1)//2 + 2*n*j + 2*j), (-1)**n); n += 1; iq = mul(iq, inv(omq(n)))
        tot = add(tot, mul(inner, i2)); j += 1; i2 = mul(i2, inv(omq(2*j)))
    report("(N5)  psi1* = (q^2;q^2)_inf sum (-1)^n q^(n(n+1)/2+2nj+2j)/((q)_n (q^2;q^2)_j)", psi1, mul(P2, tot))

def part_b():
    from mpmath import mp, mpf, mpc, exp, pi, sqrt, log, sin, cos, clsin, nstr, radians, fabs
    cl2 = clsin(2, pi/3); A = sqrt(pi/sqrt(3))
    def sums(q):
        """psi1*, f1*, psi1, f1 at q"""
        ps = fs = p = f = mpc(0); tp = tps = mpc(1); tf = tfs = mpc(1); qn = mpc(1); n = 0; eps = mpf(10)**(-mp.dps + 5)
        while True:
            ps += tps; p += tp; fs += tfs; f += tf; n += 1; qn *= q
            tp *= qn*(1 + qn); tps *= -qn*(1 + qn); tf *= qn*qn/(1 + qn); tfs *= -qn*qn/(1 + qn)
            if fabs(tp) < eps*max(1, fabs(p)) and n > 30 and fabs(tp) < 1: break
        return ps, fs, p, f
    print("    eps   |t|    psi1*: measured, lambda1       f1*: measured, lambda2 (or |f1*| if lambda2 < 0)      psi1: |t| log|psi1|, (pi^2/60) sin eps      |f1|")
    for deg in (1, 2, 5, 10):
        e = radians(deg); l1 = cl2*cos(e) - pi**2/3*sin(e); l2 = cl2*cos(e) - 7*pi**2/3*sin(e)
        for T in ("0.2", "0.1", "0.05"):
            T = mpf(T)
            with mp.workdps(40 + int(0.11/(T*sin(e)))):
                t = T*exp(1j*(e - pi/2)); ps, fs, p, f = sums(exp(-t))
                m1 = T*log(fabs(ps)/(A/sqrt(T))); m2 = T*log(fabs(fs)/(2*A/sqrt(T)))
                col2 = f"{nstr(m2, 6):>10s} {nstr(l2, 6):>10s}" if l2 > 0 else f"   |f1*| = {nstr(fabs(fs), 6):>10s}  "
                print(f"    {deg:3d}  {nstr(T, 3):>5s}   {nstr(m1, 6):>10s} {nstr(l1, 6):>10s}     {col2}          {nstr(T*log(fabs(p)), 5):>10s} {nstr(pi**2/60*sin(e), 5):>10s}     {nstr(fabs(f), 6)}")

def part_c():
    from mpmath import mp, mpf, mpc, nstr, fabs
    mp.dps = 30
    def F(z, q):
        s = mpc(0); term = mpc(1); qn = mpc(1)
        for n in range(400):
            s += term; qn *= q; term *= z*qn*(1 + qn)
        return s
    def Fc(v, q):
        s = mpc(0); term = mpc(1); qn = mpc(1)
        for n in range(400):
            s += term; qn *= q; term *= v*qn*qn/(1 + qn)
        return s
    worst = [mpf(0), mpf(0), mpf(0)]
    for z, q in ((mpc("0.7", "0.4"), mpc("0.3", "0.2")), (mpc(-1), mpc("0.5", "-0.3")), (mpc("-2.2", "1.1"), mpc("-0.4", "0.45")), (mpc("0.31", "-2.5"), mpc("0.6", "0.1"))):
        L = lambda y: y(z) - z*q*y(z*q) - z*q*q*y(z*q*q)
        y2 = lambda w: Fc(1/w, q)/w
        B = lambda w: F(w, q) + y2(w)/2
        worst = [max(worst[0], fabs(L(lambda w: F(w, q)) - 1)), max(worst[1], fabs(L(y2) + 2)), max(worst[2], fabs(L(B)))]
    print(f"    |L F - 1| <= {nstr(worst[0], 3)},  |L[Fcheck(1/z)/z] + 2| <= {nstr(worst[1], 3)},  |L[F + Fcheck(1/z)/(2z)]| <= {nstr(worst[2], 3)}   (4 pairs z, q)")

if __name__ == "__main__":
    parts = [a for a in sys.argv[1:] if a in "abc"] or list("abc")
    num = [int(a) for a in sys.argv[1:] if a.isdigit()]
    if "a" in parts: print("(a) series identities"); part_a(num[0] if num else 300)
    if "b" in parts: print("(b) growth along rays arg x = eps (degrees)"); part_b()
    if "c" in parts: print("(c) q-difference equations in z"); part_c()
