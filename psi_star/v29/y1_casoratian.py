"""Addendum v29: Casoratians of the trivial series with the homogeneous solutions, on the q-side and on the dual side.
 L y = y(z) - z q y(zq) - z q^2 y(z q^2);  F(z) = sum z^n q^(n(n+1)/2) (-q;q)_n,  G(z) = Fc(1/z)/(2z),  B = F + G (bilateral, L B = 0).
 a  q-side.  D(z) = B(z) F(zq) - B(zq) F(z) = G(z) F(zq) - G(zq) F(z):
      (1) D(z) + z q^2 D(zq) = -B(zq)
      (2) D(z) = -sum_{n in Z} (-1)^n q^(n(n+3)/2) s_n z^n,   s_n = sum_{j<=n} (-1)^j (-q;q)_j
      (3) D(z) = -sum_{k>=0} (-1)^k z^k q^(k(k+3)/2) B(z q^(k+1))
      (4) z = -1:  f0* psi1* - f1* psi0* = -2q sum_n q^(n(n+3)/2) s_n;     z = 1:  f1 psi0 - f0 psi1 = -2q sum_n (-1)^n q^(n(n+3)/2) s_n
      (5) B is not (multiplier) x r(1/z)/z, r(v) = sum q^(n^2+n) v^n/(q;q)_n:  B(z) y(zq) + B(zq) y(z) for y = r(1/z)/z, relative size
 b  series, exact to order N:  s_(-1) = -(1/2) sigma_-,  sigma_- = sum (-1)^k q^(k(k+1)/2)/(-q;q)_k = sum_{n>=0} q^(n(3n+1)/2) (1 - q^(2n+1))   (Rogers)
 c  dual side.  I_m = I((a+m)/b), w_m = W(T p^(-2m)), K_m = -I_(m-1) w_m + T p^(1-2m) I_m w_(m-1)   (K_0 is the state-integral part of Theorem Z):
      (1) W(T) + T p^2 W(p^2 T) - T p^3 W(p^4 T) = 1
      (2) I_(m-1) = T p^(-2m) (I_m + I_(m+1)/p)                     ((B1) of v12)
      (3) K_m + K_(m+1) = -I_m,   and   J((a+m)/b) + J((a+m+1)/b) = I_m
      (4) K_0 = sum_{j>=1} (-1)^j I_(-j)
      (5) Theorem Z: F(z) = i P e^(-pi i a) [J(a/b) + K_0]
      (6) z = 1 (a = 0): K_0 = I(0) (phi0(-p) - 1) + I(1/b) phi1(-p)/p
Usage: python3 y1_casoratian.py [a b c] [N]"""
import sys
from mpmath import mp, mpf, mpc, exp, pi, nstr, fabs, qp

def F(z, q):
    s = mpc(0); term = mpc(1); qn = mpc(1); n = 0
    while n < 30 or fabs(term) > mpf(10)**(-mp.dps - 10):
        s += term; n += 1; qn *= q; term *= z*qn*(1 + qn)
    return s
def Fc(v, q):
    s = mpc(0); term = mpc(1); qn = mpc(1); n = 0
    while n < 30 or fabs(term) > mpf(10)**(-mp.dps - 10):
        s += term; n += 1; qn *= q; term *= v*qn*qn/(1 + qn)
    return s
G = lambda z, q: Fc(1/z, q)/(2*z)
B = lambda z, q: F(z, q) + G(z, q)
def s_table(q, n0):
    """s_n for -n0 <= n <= n0: s_n = sum_{j<=n} (-1)^j (-q;q)_j, with (-q;q)_(-m) = q^(m(m-1)/2)/(2 (-q;q)_(m-1))"""
    neg = {}; tot = mpc(0); po = mpc(1)                     # tail sums from below
    terms = []
    for m in range(1, 400):
        terms.append((-1)**m*q**(m*(m - 1)//2)/(2*po)); po *= 1 + q**m
        if fabs(terms[-1]) < mpf(10)**(-mp.dps - 10) and m > n0 + 5: break
    s = {}; acc = mpc(0)
    for m in range(len(terms), 0, -1):
        acc += terms[m - 1]; s[-m] = acc
    po = mpc(1); cur = s[-1]
    for n in range(0, n0 + 1):
        if n > 0: po *= 1 + q**n
        cur += (-1)**n*po; s[n] = cur
    return s

def part_a():
    mp.dps = 40
    print("(a) q-side")
    w = [mpf(0)]*5
    for z, q in ((mpc("0.7", "0.4"), mpc("0.3", "0.2")), (mpc("-1.3", "0.2"), mpc("0.45", "-0.3")), (mpc("0.31", "-2.5"), mpc("-0.4", "0.35")), (mpc(-1), mpc("0.5", "0.1")), (mpc(1), mpc("0.2", "0.4"))):
        D = lambda zz: G(zz, q)*F(zz*q, q) - G(zz*q, q)*F(zz, q)
        d0 = D(z)
        w[0] = max(w[0], fabs(B(z, q)*F(z*q, q) - B(z*q, q)*F(z, q) - d0)/fabs(d0))
        w[1] = max(w[1], fabs(d0 + z*q*q*D(z*q) + B(z*q, q))/fabs(d0))
        s = s_table(q, 60)
        ser = -sum((-1)**n*q**(n*(n + 3)//2)*s[n]*z**n for n in range(-60, 61) if n in s)
        w[2] = max(w[2], fabs(ser - d0)/fabs(d0))
        it = -sum((-1)**k*z**k*q**(k*(k + 3)//2)*B(z*q**(k + 1), q) for k in range(0, 40))
        w[3] = max(w[3], fabs(it - d0)/fabs(d0))
        r = lambda v: sum(q**(n*n + n)*v**n/qp(q, q, n) for n in range(0, 40))
        y = lambda zz: r(1/zz)/zz
        w[4] = max(w[4], 0)
        print(f"    z = {nstr(z, 4)}, q = {nstr(q, 4)}: D = {nstr(d0, 15)};  (5) |B(z)y(zq) + B(zq)y(z)|/|B(z)y(zq)| = {nstr(fabs(B(z, q)*y(z*q) + B(z*q, q)*y(z))/fabs(B(z, q)*y(z*q)), 4)}")
    print(f"    largest relative deviations: two forms of D {nstr(w[0], 3)};  (1) {nstr(w[1], 3)};  (2) {nstr(w[2], 3)};  (3) {nstr(w[3], 3)}")
    for q in (mpc("0.3", "0.2"), mpc("-0.4", "0.3"), mpf("0.6")):
        s = s_table(q, 60)
        def fam(sgn):               # f1, f0, psi1, psi0 with sign sgn in the n-th term
            f1 = Fc(sgn, q); f0 = Fc(sgn/q, q); psi1 = F(sgn, q); psi0 = q*F(sgn*q, q)
            return f1, f0, psi1, psi0
        f1, f0, p1, p0 = fam(-1); lhs = f0*p1 - f1*p0; rhs = -2*q*sum(q**(n*(n + 3)//2)*s[n] for n in s)
        g1, g0, r1, r0 = fam(1); lhs2 = g1*r0 - g0*r1; rhs2 = -2*q*sum((-1)**n*q**(n*(n + 3)//2)*s[n] for n in s)
        print(f"    (4) q = {nstr(q, 3)}: f0* psi1* - f1* psi0* = {nstr(lhs, 15)}, deviation {nstr(fabs(lhs - rhs), 3)};   f1 psi0 - f0 psi1 = {nstr(lhs2, 15)}, deviation {nstr(fabs(lhs2 - rhs2), 3)}")

def part_b(N):
    print("(b) series")
    def mul(a, b):
        c = [0]*(N + 1)
        for i, x in enumerate(a):
            if x:
                for j in range(N + 1 - i):
                    if b[j]: c[i + j] += x*b[j]
        return c
    def inv(a):
        b = [0]*(N + 1); b[0] = a[0]
        for n in range(1, N + 1): b[n] = -a[0]*sum(a[k]*b[n - k] for k in range(1, n + 1))
        return b
    sig = [0]*(N + 1); den = [1] + [0]*N; k = 0
    while k*(k + 1)//2 <= N:
        idn = inv(den)
        for i in range(N + 1 - k*(k + 1)//2): sig[i + k*(k + 1)//2] += (-1)**k*idn[i]
        k += 1; f = [1] + [0]*N
        if k <= N: f[k] = 1
        den = mul(den, f)
    ft = [0]*(N + 1); n = 0
    while n*(3*n + 1)//2 <= N:
        ft[n*(3*n + 1)//2] += 1
        if n*(3*n + 1)//2 + 2*n + 1 <= N: ft[n*(3*n + 1)//2 + 2*n + 1] -= 1
        n += 1
    d = [i for i in range(N + 1) if sig[i] != ft[i]]
    print("    sigma_- = sum q^(n(3n+1)/2)(1 - q^(2n+1)): " + (f"holds to order q^{N}" if not d else f"FAILS first at q^{d[0]}"))

def part_c():
    from common import Setup, I
    from a1_formula import Iser, Wt
    from c2_z1 import mock5
    mp.dps = 30
    print("(c) dual side")
    w = [mpf(0)]*6
    for x in (mpc("0.3", "0.2"), mpc("0.2", "0.5"), mpc("0.11", "0.06")):
        S = Setup(x); b = S.b; p = S.p
        for a in (mpf(1)/2, mpf(1)/4, mpf("0.8"), mpc("0.3", "0.1")):
            T = p*exp(2*pi*I*a/x); z = exp(2*pi*I*a)
            Im = lambda m: Iser(S, (a + m)/b); wm = lambda m: Wt(S, T*p**(-2*m))
            K = lambda m: -Im(m - 1)*wm(m) + T*p**(1 - 2*m)*Im(m)*wm(m - 1)
            w[0] = max(w[0], fabs(Wt(S, T) + T*p*p*Wt(S, p*p*T) - T*p**3*Wt(S, p**4*T) - 1))
            for m in (-1, 0, 1):
                w[1] = max(w[1], fabs(Im(m - 1) - T*p**(-2*m)*(Im(m) + Im(m + 1)/p))/fabs(Im(m - 1)))
            for m in (-1, 0):
                w[2] = max(w[2], fabs(K(m) + K(m + 1) + Im(m))/fabs(Im(m)))
            w[2] = max(w[2], fabs(S.Ju(a/b) + S.Ju((a + 1)/b) - Im(0))/fabs(Im(0)))
            ser = sum((-1)**j*Im(-j) for j in range(1, 40))
            w[3] = max(w[3], fabs(K(0) - ser)/fabs(K(0)))
            w[4] = max(w[4], fabs(S.Psi(z) - I*S.P()*exp(-pi*I*a)*(S.Ju(a/b) + K(0))))
        ph0, ph1, F0, F1 = mock5(p); a = mpf(0); T = p
        K0 = -Iser(S, -1/b)*Wt(S, T) + T*p*Iser(S, 0)*Wt(S, T*p*p)
        w[5] = max(w[5], fabs(K0 - (Iser(S, 0)*(ph0 - 1) + Iser(S, 1/b)*ph1/p))/fabs(K0))
    print(f"    largest deviations at 3 points x and a = 1/2, 1/4, 0.8, 0.3 + 0.1i:  (1) {nstr(w[0], 3)};  (2) {nstr(w[1], 3)};  (3) {nstr(w[2], 3)};  (4) {nstr(w[3], 3)};  (5) {nstr(w[4], 3)};  (6) {nstr(w[5], 3)}")

if __name__ == "__main__":
    parts = [a for a in sys.argv[1:] if a in "abc"] or list("abc")
    num = [int(a) for a in sys.argv[1:] if a.isdigit()]
    if "a" in parts: part_a()
    if "b" in parts: part_b(num[0] if num else 300)
    if "c" in parts: part_c()
