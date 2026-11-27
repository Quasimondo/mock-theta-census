"""Addendum v18: (a) the formula in the fundamental strip 0 < Re a < 1,
        F(z) = i P e^(-pi i a) [ J(a/b) - I((a-1)/b) W(T) + p T I(a/b) W(p^2 T) ],   T = p e^(2 pi i a/x) = e^(pi i (2a-1)/x);
    (b) z = 1 (a = 0): the fifth-order mock theta functions
        psi1(q) =    i P [ J(0) + I(0) (phi0(-p) - 1) + I(1/b) phi1(-p)/p ]
        f1(q)   = -2 i P [ J(0) - I(0) F0(p) + I(1/b) F1(p) ]
        with I(0) = c D(0), I(1/b) = c p^(1/2) D(1/b) and, r = p^(1/2) = e^(-pi i/(2x)):
        D(0)   = [ G(q^4) H(-r) + q^(1/4) H(-q) H(r) ] / ((q;q^2)_inf (-p;p)_inf)
        D(1/b) = [ G(q^4) G(-r) - q^(1/4) H(-q) G(r) ] / ((q;q^2)_inf (-p;p)_inf)        (G, H the Rogers-Ramanujan functions)
Usage: python3 c2_z1.py [digits]"""
import sys
from mpmath import mp, mpf, mpc, exp, pi, sqrt, qp, nstr
from common import Setup, I
from a1_formula import Wt, Iser

def RR(y):                     # G(y) = sum y^(n^2)/(y;y)_n,  H(y) = sum y^(n^2+n)/(y;y)_n
    g = h = mpc(0); den = mpc(1); n = 0
    while abs(y)**(n*n) > mpf(10)**(-mp.dps - 8):
        g += y**(n*n)/den; h += y**(n*n + n)/den; n += 1; den *= 1 - y**n
    return g, h

def mock5(p):                  # phi0(-p), phi1(-p), F0(p), F1(p)
    ph0 = ph1 = F0 = F1 = mpc(0); po = mpc(1); ipo = mpc(1); n = 0
    while abs(p)**(n*n) > mpf(10)**(-mp.dps - 8):
        ph0 += (-1)**n*p**(n*n)*po; ph1 += (-1)**(n + 1)*p**((n + 1)**2)*po     # phi(-p): q^(n^2) -> (-1)^n p^(n^2), (-q;q^2)_n -> (p;p^2)_n
        F0 += p**(2*n*n)/ipo; ipo *= 1 - p**(2*n + 1); F1 += p**(2*n*n + 2*n)/ipo
        po *= 1 - p**(2*n + 1); n += 1
    return ph0, ph1, F0, F1

if __name__ == "__main__":
    mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    for x in (mpc("0.3", "0.2"), mpc("0.11", "0.03")):
        S = Setup(x); b = S.b; p = S.p; q = S.q
        print(f"x = {nstr(x, 6)}")
        for name, a in (("i", mpf(1)/4), ("-omega", mpf(5)/6)):
            z = exp(2*pi*I*a); T = p*exp(2*pi*I*a/x)
            rhs = I*S.P()*exp(-pi*I*a)*(S.Ju(a/b) - Iser(S, (a - 1)/b)*Wt(S, T) + p*T*Iser(S, a/b)*Wt(S, p*p*T))
            print(f"  (a) z = {name}: |F - strip formula| = {nstr(abs(S.Psi(z) - rhs), 3)};  |T| = {nstr(abs(T), 4)}, |W(T) - 1| = {nstr(abs(Wt(S, T) - 1), 3)}, |W(p^2 T) - 1| = {nstr(abs(Wt(S, p*p*T) - 1), 3)}")
        J0 = S.Ju(0); I0 = Iser(S, 0); I1 = Iser(S, 1/b)
        ph0, ph1, F0, F1 = mock5(p)
        psi1 = S.Psi(1); f1 = S.fstar(1)
        print(f"  (b) psi1(q) = {nstr(psi1, 18)}:  |psi1 - formula| = {nstr(abs(psi1 - I*S.P()*(J0 + I0*(ph0 - 1) + I1*ph1/p)), 3)}")
        print(f"      f1(q)   = {nstr(f1, 18)}:  |f1 - formula|   = {nstr(abs(f1 + 2*I*S.P()*(J0 - I0*F0 + I1*F1)), 3)}")
        r = S.ppow(mpf(1)/2); q4 = S.qpow(mpf(1)/4)
        G4, _ = RR(q**4); _, Hm = RR(-q); Gr, Hr = RR(r); Gmr, Hmr = RR(-r)
        den = qp(q, q*q)*qp(-p, p)
        D0 = (G4*Hmr + q4*Hm*Hr)/den; D1 = (G4*Gmr - q4*Hm*Gr)/den
        print(f"      |D(0) - RR form|/|D(0)| = {nstr(abs(S.D(0) - D0)/abs(D0), 3)},  |D(1/b) - RR form|/|D(1/b)| = {nstr(abs(S.D(1/b) - D1)/abs(D1), 3)}")
