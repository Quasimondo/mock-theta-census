"""Finite-x test of the three-sector formula for psi1*.
XX = K - i b,  q = e^(2 pi i/XX) = e^(-t),  t = -2 pi i/XX,  p = e^(-pi i XX).   Conjecture:
    psi1*(q) = S[Phi](t) + (kappa/2) [ B2 * p W3(p) - B4 * W(p) ]
with  S[Phi] = Borel sum of the perturbative series Phi of psi1* at q -> 1, taken along the direction of t,
      B2 = X(q) E+(p) + i q^(1/4) N1(q) E+'(p),   B4 = X(q) E-(p) - i q^(1/4) N1(q) E-'(p)      (the state integrals of v6/v7),
      kappa = theta4(q) (p;p^2)_inf / (2 i q^(1/4) (-q;q)_inf psi(p)).
The script prints  T = psi1* - (kappa/2)[...]  against the Borel-Pade-Laplace value of S[Phi], for two Pade orders."""
import sys, os, json
from mpmath import mp, mpf, mpc, exp, pi, sqrt, nstr, quad, pade, polyval, factorial
HERE = os.path.dirname(os.path.abspath(__file__))
I = mpc(0, 1)
SER = json.load(open(os.path.join(HERE, "q2_series_v3.json")))["psi1@1"]

def borel_sum(t, N):
    old = mp.dps; mp.dps = 150
    c = [mpf(SER[k][0])/factorial(k) for k in range(2*N + 1)]
    P, Q = pade(c, N, N)
    P = P[::-1]; Q = Q[::-1]
    mp.dps = old
    f = lambda u: exp(-u)*polyval(P, u*t)/polyval(Q, u*t)
    return quad(f, [0, 5, 15, 30, 60, 120, 250])

def run(K, b):
    XX = mpc(K, -b); t = -2*pi*I/XX; q = exp(-t); p = exp(-pi*I*XX); q4 = exp(-t/4)
    tiny = mpf(10)**(-mp.dps - 5)
    # q side
    psi1 = mpc(0); tri = mpc(1); poch = mpc(1); n = 0
    while abs(q)**(n*(n + 1)//2)*mpf(10)**12 > tiny:
        psi1 += (-1)**n*tri*poch; n += 1; tri *= q**n; poch *= 1 + q**n
    X = N1 = mpc(0); pe = mpc(1); n = 0
    while abs(q)**(n*n)*mpf(10)**12 > tiny:
        po = pe*(1 - q**(2*n + 1))
        X += (-1)**n*q**(n*n)/pe; N1 += (-1)**n*q**(n*n + n)/po
        n += 1; pe = po*(1 - q**(2*n))
    # dual side
    Ep = Em = Epp = Emp = mpc(0); den = mpc(1); m = 0
    while abs(p)**(m*(m + 1)//2) > tiny:
        a = p**(m*(m + 1)//2)/den
        Epp += a; Emp += a*p**m; Ep += (-1)**m*a; Em += (-1)**m*a*p**m
        m += 1; den *= 1 - p**(2*m)
    W = W3 = mpc(0); po = mpc(1); j = 0
    while abs(p)**(j*j) > tiny:
        W += (-1)**j*p**(j*j + j)*po; W3 += (-1)**j*p**(j*j + 3*j)*po
        po *= 1 - p**(2*j + 1); j += 1
    psip = sum(p**(k*(k + 1)//2) for k in range(60))
    th4 = sqrt(pi/t)*sum(exp(-pi**2*(n + mpf(1)/2)**2/t) for n in range(-40, 40))
    kappa = th4/(2*I*q4*exp(pi**2/(12*t) + t/24)/sqrt(2)*psip)
    B2 = X*Ep + I*q4*N1*Epp; B4 = X*Em - I*q4*N1*Emp
    chk = abs(I*q4*N1 - (B2*Em - B4*Ep)/2)/abs(N1)
    dec = kappa/2*B2*p*W3; gro = -kappa/2*B4*W
    T = psi1 - dec - gro
    S60 = borel_sum(t, 60); S50 = borel_sum(t, 50)
    print(f"XX = {K} - {b} i:  |psi1*| = {nstr(abs(psi1), 4)}   |growing term| = {nstr(abs(gro), 4)}   |decaying term| = {nstr(abs(dec), 4)}   (inversion check {nstr(chk, 2)})")
    print(f"     T        = {nstr(T, 22)}")
    print(f"     S[Phi]   = {nstr(S60, 22)}     (Pade 60 vs 50: {nstr(abs(S60 - S50), 2)})")
    d = T - S60
    print(f"     |T - S|  = {nstr(abs(d), 3)}    (T - S)/(kappa B2 W/2) = {nstr(d/(kappa/2*B2*W), 10)}    (T - S)/decaying term = {nstr(d/dec, 10)}")
    sys.stdout.flush()

if __name__ == "__main__":
    mp.dps = 60
    for K, b in ((50, "0.5"), ("50.7", "0.5"), ("60.3", "0.4"), ("71.1", "0.6"), ("80.4", "0.5")):
        run(mpf(K), mpf(b))
