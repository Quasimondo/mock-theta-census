# Independent check of addendum-v29 (own code).
# (1) Proposition C2: D = B(z)F(zq) - B(zq)F(z) = G(z)F(zq) - G(zq)F(z); D(z) + zq^2 D(zq) = -B(zq); the Laurent
#     series (ii) with s_n; the sum over shifts (iii). B(-1) = theta4 R1/2.
# (2) The Corollary: f0* psi1* - f1* psi0* and f1 psi0 - f0 psi1 against -2q sum (+-1)^n q^{n(n+3)/2} s_n,
#     numerically at three q and the second also as an exact power series (Fractions) to q^120.
# (3) Rogers' false theta identity sigma_- = sum q^{n(3n+1)/2}(1 - q^{2n+1}), exact to q^200.
# (4) Proposition C3: (i) the equation for calW; (iv) K0 = sum_{j>=1} (-1)^j I((a-j)/b) against its definition, with
#     I from the double series of Theorem C (own code); (v) J(a/b) + J((a+1)/b) = I(a/b) by quadrature.
import mpmath as mp
from fractions import Fraction
mp.mp.dps = 40
def pm(q, n):                      # (-q;q)_n for all integers n, with (-q;q)_{-m} = q^{m(m-1)/2}/(2(-q;q)_{m-1})
    if n >= 0:
        p = mp.mpc(1)
        for j in range(1, n+1): p *= 1 + q**j
        return p
    m = -n
    return q**(m*(m-1)//2)/(2*pm(q, m-1))
def F(z, q, N=120): return sum(z**n*q**(n*(n+1)//2)*pm(q, n) for n in range(N))
def B(z, q, N=120): return sum(z**n*q**(n*(n+1)//2)*pm(q, n) for n in range(-N, N))
def G(z, q): return B(z, q) - F(z, q)
def Fc(v, q, N=120):               # F-check(v) = sum q^{n^2+n} v^n/(-q;q)_n
    return sum(q**(n*n+n)*v**n/pm(q, n) for n in range(N))
print('(1) Proposition C2')
for z, q in ((mp.mpc(0.7, 0.4), mp.mpc(0.3, 0.2)), (mp.mpc(-1.3, 0.2), mp.mpc(0.45, -0.3)), (mp.mpc(1, 0), mp.mpc(0.2, 0.4))):
    D1 = B(z, q)*F(z*q, q) - B(z*q, q)*F(z, q); D2 = G(z, q)*F(z*q, q) - G(z*q, q)*F(z, q)
    rec = D1 + z*q*q*(B(z*q, q)*F(z*q*q, q) - B(z*q*q, q)*F(z*q, q)) + B(z*q, q)
    def s(n): return sum((-1)**j*pm(q, j) for j in range(-80, n+1))
    lau = -sum((-1)**n*q**(n*(n+3)//2)*s(n)*z**n for n in range(-40, 40))
    shifts = -sum((-1)**k*z**k*q**(k*(k+3)//2)*B(z*q**(k+1), q) for k in range(60))
    g_ok = abs(G(z, q) - Fc(1/z, q)/(2*z))
    print(f'  z = {mp.nstr(z,3)}, q = {mp.nstr(q,3)}: |D1-D2| {mp.nstr(abs(D1-D2),2)}  recursion {mp.nstr(abs(rec),2)}  Laurent {mp.nstr(abs(lau-D1),2)}  shifts {mp.nstr(abs(shifts-D1),2)}  G = Fcheck(1/z)/(2z): {mp.nstr(g_ok,2)}   (|D| = {mp.nstr(abs(D1),3)})')
q = mp.mpc(0.5, 0.1)
th4 = 1 + 2*sum((-1)**n*q**(n*n) for n in range(1, 40))
R1 = sum((-1)**n*q**(n*n+n)/mp.qp(q, q, n) for n in range(60))
print('  B(-1) - theta4 R1/2 =', mp.nstr(abs(B(-1, q) - th4*R1/2), 3))
print('(2) the Corollary')
for q in (mp.mpc(0.3, 0.2), mp.mpc(-0.4, 0.3), mp.mpf(0.6)):
    def s(n): return sum((-1)**j*pm(q, j) for j in range(-80, n+1))
    psi1s = sum((-1)**n*q**(n*(n+1)//2)*pm(q, n) for n in range(150)); psi0s = sum((-1)**n*q**((n+1)*(n+2)//2)*pm(q, n) for n in range(150))
    f1s = sum((-1)**n*q**(n*n+n)/pm(q, n) for n in range(80)); f0s = sum((-1)**n*q**(n*n)/pm(q, n) for n in range(80))
    psi1 = sum(q**(n*(n+1)//2)*pm(q, n) for n in range(150)); psi0 = sum(q**((n+1)*(n+2)//2)*pm(q, n) for n in range(150))
    f1 = sum(q**(n*n+n)/pm(q, n) for n in range(80)); f0 = sum(q**(n*n)/pm(q, n) for n in range(80))
    S1 = -2*q*sum(q**(n*(n+3)//2)*s(n) for n in range(-40, 40)); S2 = -2*q*sum((-1)**n*q**(n*(n+3)//2)*s(n) for n in range(-40, 40))
    print(f'  q = {mp.nstr(q,3)}: f0*psi1* - f1*psi0* - rhs = {mp.nstr(abs(f0s*psi1s - f1s*psi0s - S1),2)};  f1 psi0 - f0 psi1 - rhs = {mp.nstr(abs(f1*psi0 - f0*psi1 - S2),2)}')
N = 120
def pmul(a, b):
    c = [Fraction(0)]*(N+1)
    for i, x in enumerate(a):
        if x:
            for j in range(N+1-i): c[i+j] += x*b[j]
    return c
def pinv(a):
    b = [Fraction(0)]*(N+1); b[0] = Fraction(1)/a[0]
    for n in range(1, N+1): b[n] = -sum(a[k]*b[n-k] for k in range(1, n+1))/a[0]
    return b
def mono(e, c=1):
    s = [Fraction(0)]*(N+1)
    if 0 <= e <= N: s[e] = Fraction(c)
    return s
def padd(*xs): return [sum(v) for v in zip(*xs)]
def pp(n):                         # (-q;q)_n as power series, all n
    if n >= 0:
        p = mono(0)
        for j in range(1, n+1):
            f = mono(0); f[j] += 1 if j <= N else 0; p = pmul(p, f)
        return p
    m = -n
    return pmul(mono(m*(m-1)//2, Fraction(1, 2)), pinv(pp(m-1)))
ppc = {n: pp(n) for n in range(-20, 20)}
def ser(sgnexp, num, nmax):
    s = mono(0, 0)
    for n in range(nmax):
        e = sgnexp(n)
        if e > N: break
        s = padd(s, pmul(mono(e), ppc[n] if num else pinv(ppc[n])))
    return s
psi1 = ser(lambda n: n*(n+1)//2, True, 20); psi0 = ser(lambda n: (n+1)*(n+2)//2, True, 20)
f1 = ser(lambda n: n*n+n, False, 20); f0 = ser(lambda n: n*n, False, 20)
lhs = padd(pmul(f1, psi0), [-x for x in pmul(f0, psi1)])
sn = {}
acc = mono(0, 0)
for j in range(-19, 19):
    acc = padd(acc, [(-1)**j*x for x in ppc[j]]); sn[j] = acc[:]
rhs = mono(0, 0)
for n in range(-19, 18):
    e = n*(n+3)//2 + 1
    if e <= N: rhs = padd(rhs, [-2*(-1)**n*x for x in pmul(mono(e), sn[n])])
print('  f1 psi0 - f0 psi1 = -2q sum (-1)^n q^{n(n+3)/2} s_n, exact to q^%d:' % N, all(a == b for a, b in zip(lhs, rhs)))
sig = mono(0, 0)
for k in range(20):
    if k*(k+1)//2 > N: break
    sig = padd(sig, pmul(mono(k*(k+1)//2, (-1)**k), pinv(ppc[k])))
rog = mono(0, 0)
for n in range(20): rog = padd(rog, mono(n*(3*n+1)//2), mono(n*(3*n+1)//2 + 2*n + 1, -1))
print('(3) Rogers false theta identity, exact to q^%d:' % N, all(a == b for a, b in zip(sig, rog)))
print('(4) Proposition C3')
mp.mp.dps = 30
def calW(T, p, K=60):
    s = 0; po = mp.mpc(1)
    for k in range(K): s += (-T)**k*p**(k*k+k)*po; po *= 1 - p**(2*k+1)
    return s
for x in (mp.mpc(0.3, 0.2), mp.mpc(-0.2, 0.25)):
    p = mp.expjpi(-1/x)
    for T in (mp.mpc(0.4, 0.3), mp.mpc(-1.2, 0.5)):
        print(f'  (i) x = {mp.nstr(x,3)}, T = {mp.nstr(T,3)}: {mp.nstr(abs(calW(T,p) + T*p**2*calW(p**2*T,p) - T*p**3*calW(p**4*T,p) - 1), 2)}')
def pinf(a, q):
    s = mp.mpc(1); k = 0
    while True:
        f = a*q**k; s *= 1 + f
        if abs(f) < mp.mpf(10)**(-mp.mp.dps-5): return s
        k += 1
for x in (mp.mpc(0.3, 0.2), mp.mpc(0.15, 0.25)):
    b = mp.sqrt(x); q = mp.expjpi(2*x); p = mp.expjpi(-1/x)
    c = 1/mp.sqrt(2)*mp.expjpi(-mp.mpf(1)/4)*mp.expjpi(-x/12)*mp.expjpi(-1/(12*x))
    def I(u, K=50):
        D = 0; qk = mp.mpc(1)
        for k in range(K):
            if k: qk *= 1 - q**k
            pn = mp.mpc(1)
            for n in range(K):
                if n: pn *= 1 - p**(2*n)
                D += (-1)**(k*n+n)*mp.exp(mp.pi*1j*u*(b*k + n/b))*mp.expjpi(x*k*k/2)*mp.expjpi(-(n*n/mp.mpf(2) + n)/x)/(qk*pn)
        return c*mp.exp(-mp.pi*1j*u**2/2)*D
    Phi = lambda y: pinf(mp.expjpi(x)*mp.exp(2*mp.pi*b*y), q)/pinf(p*mp.exp(2*mp.pi*y/b), p**2)
    e = mp.expj(-(mp.pi/2 - mp.arg(b))/2)
    J = lambda u: mp.quad(lambda r: mp.exp(-mp.pi*1j*(r*e)**2 + 2*mp.pi*u*r*e)/(Phi(r*e)*(1 + mp.exp(2*mp.pi*r*e/b)))*e, mp.linspace(-6, 6, 25))
    for a in (mp.mpf(1)/4, mp.mpc(0.3, 0.1)):
        T = mp.expjpi((2*a - 1)/x)
        K0 = -I((a-1)/b)*calW(T, p) + p*T*I(a/b)*calW(p*p*T, p)
        Ksum = sum((-1)**j*I((a-j)/b) for j in range(1, 14))
        Jv = J(a/b) + J((a+1)/b) - I(a/b)
        P = 2**mp.mpf(-0.5)*mp.exp(mp.pi**2/(12*(-2j*mp.pi*x)) + (-2j*mp.pi*x)/6)/b
        z = mp.expj(2*mp.pi*a); Fz = sum(z**n*q**(n*(n+1)//2)*pm(q, n) for n in range(200))
        thz = abs(Fz - 1j*P*mp.expjpi(-a)*(J(a/b) + K0))
        print(f'  x = {mp.nstr(x,3)}, a = {mp.nstr(a,3)}: (iv) |K0 - sum| = {mp.nstr(abs(K0 - Ksum),2)}  (v) |J0 + J1 - I0| = {mp.nstr(abs(Jv),2)}  Theorem Z |F - iPe^(-pi i a)(J0 + K0)| = {mp.nstr(thz,2)}', flush=True)
