"""The four bilinear combinations B1..B4 of addendum v6 against state integrals of Faddeev's quantum dilogarithm.

For b^2 = x in the upper half-plane:  q = e^(2 pi i x),  qt = e^(-2 pi i/x),  p = qt^(1/2) = e^(-pi i/x),  c_b = i (b + 1/b)/2,
    Phi_b(y) = (-q^(1/2) e^(2 pi b y); q)_inf / (-qt^(1/2) e^(2 pi y/b); qt)_inf ,
    I(u) = int_C Phi_b(y) exp(-pi i y^2/2 + 2 pi u y) dy ,     C: from infinity e^(i(pi - delta)) through 0 to infinity e^(i delta).
Claim (residues at the poles y = c_b + i b k + i n/b, k, n >= 0):
    I(u_j) = -i b (q;q)_inf/(qt;qt)_inf exp(-pi i c_b^2/2 + 2 pi u_j c_b) B_j ,
    u_1 = b/4,  u_2 = -b/4,  u_3 = b/4 - 1/(2b),  u_4 = -b/4 - 1/(2b).
Usage: python3 state_integral.py [dps]"""
import sys
from mpmath import mp, mpf, mpc, exp, pi, sqrt, qp, quad, nstr
mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 30
I = mpc(0, 1)

def bilinear(x):
    q = exp(2*pi*I*x); p = exp(-pi*I/x); q4 = exp(2*pi*I*x/4)
    N0 = N1 = X = Y = mpc(0)
    n = 0; pe = mpc(1)                       # (q)_(2n)
    while True:
        po = pe*(1 - q**(2*n + 1)); s = (-1)**n
        a = q**(n*n)/pe; c = q**(n*n)/po
        X += s*a; N0 += s*a*q**n; N1 += s*c*q**n; Y += s*c*q**(2*n + 1)
        if abs(q)**(n*n) < mpf(10)**(-mp.dps - 25): break
        n += 1; pe = po*(1 - q**(2*n))
    Ep = Em = Epp = Emp = mpc(0); m = 0; den = mpc(1)
    while True:
        a = p**(m*(m + 1)//2)/den
        Epp += a; Emp += a*p**m; Ep += (-1)**m*a; Em += (-1)**m*a*p**m
        if abs(a) < mpf(10)**(-mp.dps - 25): break
        m += 1; den *= 1 - p**(2*m)
    return [N0*Epp + I/q4*Y*Ep, X*Ep + I*q4*N1*Epp, N0*Emp - I/q4*Y*Em, X*Em - I*q4*N1*Emp]

def integrals(x, delta=mpf(2)/5, L=12):
    b = sqrt(x); q = exp(2*pi*I*x); qt = exp(-2*pi*I/x); qh = exp(pi*I*x); qth = exp(-pi*I/x)
    def Phi(y): return qp(-qh*exp(2*pi*b*y), q)/qp(-qth*exp(2*pi*y/b), qt)
    out = []
    for u in (b/4, -b/4, b/4 - 1/(2*b), -b/4 - 1/(2*b)):
        f = lambda y: Phi(y)*exp(-pi*I*y*y/2 + 2*pi*u*y)
        eL = exp(-I*delta); eR = exp(I*delta)
        v = quad(lambda s: f(s*eL)*eL, [-L, -8, -5, -3, -1.5, 0]) + quad(lambda s: f(s*eR)*eR, [0, 1.5, 3, 5, 8, L])
        cb = I*(b + 1/b)/2
        pref = -I*b*qp(q, q)/qp(qt, qt)*exp(-pi*I*cb*cb/2 + 2*pi*u*cb)
        out.append(v/pref)
    return out

if __name__ == "__main__":
    for x in (1/mpc(3, -mpf(2)/5), mpc("0.3", "0.2"), mpc("1.2", "0.5"), mpc("0.11", "0.03")):
        B = bilinear(x); J = integrals(x)
        print("x =", nstr(x, 8))
        for j in range(4):
            print(f"   B{j+1} = {nstr(B[j], 20)}    integral/prefactor = {nstr(J[j], 20)}    |ratio - 1| = {nstr(abs(J[j]/B[j] - 1), 3)}")
        sys.stdout.flush()
