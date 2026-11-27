"""Two integrals for the one-sided sum psi1*.   b^2 = x (Im x > 0), q = e^(2 pi i x) = e^(-t), p = e^(-pi i/x), qt = p^2,
    Phi_b(y) = (-q^(1/2) e^(2 pi b y); q)_inf / (-p e^(2 pi y/b); qt)_inf.
(A) full line with a kernel:   JA = int_L  e^(-pi i y^2 + pi y/b) / ( Phi_b(y) (1 + e^(2 pi y/b)) ) dy ,   L = e^(-i delta) R.
    Residues at y = i b (n + 1/2), n >= 0, give   JA = b q^(1/8) (p;qt)_inf/(-q;q)_inf * psi1*(q).
(B) half line, no kernel:      JB = int from infinity e^(i(pi - delta)) to 0 of  e^(-pi i y^2 + pi y/b) / Phi_b(y) dy ,
    compared, with the same prefactor, to  T = psi1* - (kappa/2)[B2 p W3 - B4 W]  (v10: the Borel sum of the perturbative series)."""
import sys
from mpmath import mp, mpf, mpc, exp, pi, sqrt, qp, quad, nstr
I = mpc(0, 1)

def pieces(x):
    t = -2*pi*I*x; q = exp(-t); p = exp(-pi*I/x); q4 = exp(-t/4)
    tiny = mpf(10)**(-mp.dps - 8)
    psi1 = mpc(0); tri = mpc(1); poch = mpc(1); n = 0
    while abs(q)**(n*(n + 1)//2)*mpf(10)**6 > tiny:
        psi1 += (-1)**n*tri*poch; n += 1; tri *= q**n; poch *= 1 + q**n
    X = N1 = mpc(0); pe = mpc(1); n = 0
    while abs(q)**(n*n)*mpf(10)**6 > tiny:
        po = pe*(1 - q**(2*n + 1)); X += (-1)**n*q**(n*n)/pe; N1 += (-1)**n*q**(n*n + n)/po
        n += 1; pe = po*(1 - q**(2*n))
    Ep = Em = Epp = Emp = mpc(0); den = mpc(1); m = 0
    while abs(p)**(m*(m + 1)//2) > tiny:
        a = p**(m*(m + 1)//2)/den
        Epp += a; Emp += a*p**m; Ep += (-1)**m*a; Em += (-1)**m*a*p**m
        m += 1; den *= 1 - p**(2*m)
    W = W3 = mpc(0); po = mpc(1); j = 0
    while abs(p)**(j*j) > tiny:
        W += (-1)**j*p**(j*j + j)*po; W3 += (-1)**j*p**(j*j + 3*j)*po; po *= 1 - p**(2*j + 1); j += 1
    psip = sum(p**(k*(k + 1)//2) for k in range(80))
    kappa = qp(q, q)/qp(-q, q)*qp(p, p*p)/(2*I*q4*qp(-q, q)*psip)        # theta4(q) = (q)_inf/(-q)_inf
    B2 = X*Ep + I*q4*N1*Epp; B4 = X*Em - I*q4*N1*Emp
    dec = kappa/2*B2*p*W3; gro = -kappa/2*B4*W
    return psi1, dec, gro, W, kappa*B2/2

def integrals(x, delta=mpf(2)/5, L=9):
    b = sqrt(x); q = exp(2*pi*I*x); p = exp(-pi*I/x); qt = p*p; qh = exp(pi*I*x)
    Phi = lambda y: qp(-qh*exp(2*pi*b*y), q)/qp(-p*exp(2*pi*y/b), qt)
    g = lambda y: exp(-pi*I*y*y + pi*y/b)/Phi(y)
    e = exp(-I*delta)
    JA = quad(lambda s: g(s*e)/(1 + exp(2*pi*s*e/b))*e, [-L, -6, -4, -2, -1, 0, 1, 2, 4, 6, L])
    JB = quad(lambda s: g(s*e)*e, [-L, -6, -4, -2, -1, 0])
    pref = b*exp(2*pi*I*x/8)*qp(p, qt)/qp(-q, q)
    return JA/pref, JB/pref

if __name__ == "__main__":
    mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 25
    for x in (mpc("0.3", "0.2"), 1/mpc(3, -mpf(2)/5), mpc("0.11", "0.03")):
        psi1, dec, gro, W, kB2 = pieces(x); T = psi1 - dec - gro
        A, Bv = integrals(x)
        print("x =", nstr(x, 8), "  arg t (deg) =", nstr(mp.arg(-2*pi*I*x)*180/pi, 5), "  |growing| =", nstr(abs(gro), 3), " |decaying| =", nstr(abs(dec), 3))
        print("   psi1*          =", nstr(psi1, 20))
        print("   (A)/prefactor  =", nstr(A, 20), "   |A/psi1* - 1| =", nstr(abs(A/psi1 - 1), 3))
        print("   T              =", nstr(T, 20))
        print("   (B)/prefactor  =", nstr(Bv, 20), "   |B - T| =", nstr(abs(Bv - T), 3), "   (B - T)/gro =", nstr((Bv - T)/gro, 8), "  (B - T)/dec =", nstr((Bv - T)/dec, 8), "  (B - T)/(kappa B2 W/2) =", nstr((Bv - T)/(kB2*W), 8))
        sys.stdout.flush()
