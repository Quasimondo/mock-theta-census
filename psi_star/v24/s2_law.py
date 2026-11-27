"""Addendum v24: the family at z = 1 and z = q.  t = -2 pi i x, q = e^(-t), p = e^(-2 pi^2/t), r = p^(1/2),
   s1 = sin(pi/5), s2 = sin(2 pi/5), c- = (2/sqrt 5) s1 sqrt(pi/t), c+ = (2/sqrt 5) s2 sqrt(pi/t), e = e^(-pi i/4).
  a  Proposition M: the four state integrals are monomials
        I(0)   =  e (2/sqrt 5) s2 q^(-1/60)  p^(1/30)        I(1/b)     =  e (2/sqrt 5) s1 q^(-1/60)  p^(19/30)
        I(b)   =  e (2/sqrt 5) s1 q^(-19/60) p^(1/30)        I(b + 1/b) = -e (2/sqrt 5) s2 q^(-19/60) p^(19/30)
     (series of Theorem C of v12, and for I(0) the integral itself)
  b  Theorem 16: with M1 = 2 i q^(11/60) P J(1/b), M0 = 2 i q^(29/60) P J(b + 1/b)
        q^(11/60) f1     =  M1 + 2 c+ p^(-1/120) (F0(p) - 1) - 2 c- p^(71/120) F1(p)
        2 q^(11/60) psi1 = -M1 + 2 c+ p^(-1/120) phi0(-p)    + 2 c- p^(-49/120) phi1(-p)
        q^(-1/60) f0     =  M0 + 2 c- p^(-1/120) (F0(p) - 1) + 2 c+ p^(71/120) F1(p)
        2 q^(-1/60) psi0 = -M0 + 2 c- p^(-1/120) phi0(-p)    - 2 c+ p^(-49/120) phi1(-p)
  c  Theorem 17: M1 = sqrt(60 t/pi) L(2/5, 10 t),  M0 = sqrt(60 t/pi) L(1/5, 10 t),
        L(r, a) = int_0^inf e^(-3 a x^2/2) [cosh((3r-2) a x) + cosh((3r-1) a x)]/cosh(3 a x/2) dx
  d  the quoted law for chi_0, chi_1 (Gordon-McIntosh, as printed in arXiv:2604.20034, (6a), (6b)), on its own, at alpha = 2t
  e  the transformation formulas of the Rogers-Ramanujan functions used in the proof of Proposition M, and the eta identity
Usage: python3 s2_law.py [a b c d e] [digits]"""
import sys
from mpmath import mp, mpf, mpc, exp, pi, sqrt, sin, cosh, quad, nstr, arg, linspace, fabs, qp
from common import Setup, I
from a1_formula import Iser
from c2_z1 import mock5, RR

def L(r, a):                      # continued to complex a (Re a > 0) by turning the ray: x = e^(-i arg(a)/2) rho
    e = exp(-I*arg(a)/2); A = fabs(a)
    f = lambda rho: exp(-1.5*A*rho*rho)*(cosh((3*r - 2)*a*e*rho) + cosh((3*r - 1)*a*e*rho))/cosh(1.5*a*e*rho)
    return e*quad(f, linspace(0, 14/sqrt(A), 29))

def chis(Q):                      # chi0(Q), chi1(Q)
    c0 = c1 = mpc(0); n = 0
    while abs(Q)**n > mpf(10)**(-mp.dps - 8):
        d = mpc(1)
        for j in range(n + 1, 2*n + 1): d *= 1 - Q**j
        c0 += Q**n/d; c1 += Q**n/(d*(1 - Q**(2*n + 1))); n += 1
    return c0, c1

PTS = [mpc("0.3", "0.2"), mpc("0.11", "0.03"), mpc("0.2", "0.5"), mpc(0, "0.3"), mpc("0.05", "0.4"), mpc("0.6", "0.25")]

if __name__ == "__main__":
    parts = [a for a in sys.argv[1:] if a in "abcde"] or list("abcde")
    dig = [int(a) for a in sys.argv[1:] if a.isdigit()]
    mp.dps = dig[0] if dig else 30
    s5 = sqrt(mpf(5)); s1 = sin(pi/5); s2 = sin(2*pi/5); e4 = exp(-I*pi/4); k = 2/s5
    for x in PTS:
        S = Setup(x); t = -2*pi*I*x; b = S.b; p = S.p; q = S.q; P = S.P()
        cm = k*s1*sqrt(pi/t); cp = k*s2*sqrt(pi/t)
        qp_ = lambda a: S.qpow(mpf(a)); pp = lambda a: S.ppow(mpf(a))
        print(f"x = {nstr(x, 5)}   (t = {nstr(t, 5)})")
        if "a" in parts:
            ex = [(0, e4*k*s2*S.qpow(-mpf(1)/60)*S.ppow(mpf(1)/30), "I(0)"),
                  (1/b, e4*k*s1*S.qpow(-mpf(1)/60)*S.ppow(mpf(19)/30), "I(1/b)"),
                  (b, e4*k*s1*S.qpow(-mpf(19)/60)*S.ppow(mpf(1)/30), "I(b)"),
                  (b + 1/b, -e4*k*s2*S.qpow(-mpf(19)/60)*S.ppow(mpf(19)/30), "I(b+1/b)")]
            out = [f"{nm}: {nstr(abs(Iser(S, u)/v - 1), 3)}" for u, v, nm in ex]
            print("  (a) |series/monomial - 1|:  " + ",  ".join(out) + f";   integral I(0): {nstr(abs(S.Iu(0)/ex[0][1] - 1), 3)}")
        if set(parts) & set("bc"):
            M1 = 2*I*S.qpow(mpf(11)/60)*P*S.Ju(1/b); M0 = 2*I*S.qpow(mpf(29)/60)*P*S.Ju(b + 1/b)
        if "b" in parts:
            ph0, ph1, F0, F1 = mock5(p)
            f1 = S.fstar(1); psi1 = S.Psi(1); f0 = S.fstar(1/q); psi0 = q*S.Psi(q)
            A = S.ppow(-mpf(1)/120); B = S.ppow(mpf(71)/120); C = S.ppow(-mpf(49)/120)
            d = [abs(S.qpow(mpf(11)/60)*f1 - (M1 + 2*cp*A*(F0 - 1) - 2*cm*B*F1)),
                 abs(2*S.qpow(mpf(11)/60)*psi1 - (-M1 + 2*cp*A*ph0 + 2*cm*C*ph1)),
                 abs(S.qpow(-mpf(1)/60)*f0 - (M0 + 2*cm*A*(F0 - 1) + 2*cp*B*F1)),
                 abs(2*S.qpow(-mpf(1)/60)*psi0 - (-M0 + 2*cm*A*ph0 - 2*cp*C*ph1))]
            print("  (b) |law - series|:  f1 " + nstr(d[0], 3) + ",  psi1 " + nstr(d[1], 3) + ",  f0 " + nstr(d[2], 3) + ",  psi0 " + nstr(d[3], 3)
                  + f";   f1 = {nstr(f1, 12)}, psi0 = {nstr(psi0, 12)}")
        if "c" in parts:
            m1 = sqrt(60*t/pi)*L(mpf(2)/5, 10*t); m0 = sqrt(60*t/pi)*L(mpf(1)/5, 10*t)
            print(f"  (c) M1 = {nstr(M1, 20)},  |M1 - sqrt(60 t/pi) L(2/5, 10 t)| = {nstr(abs(M1 - m1), 3)}")
            print(f"      M0 = {nstr(M0, 20)},  |M0 - sqrt(60 t/pi) L(1/5, 10 t)| = {nstr(abs(M0 - m0), 3)}")
            j1 = e4*sqrt(mpf(15))*t/pi*S.qpow(-mpf(1)/60)*S.ppow(mpf(1)/24)*L(mpf(2)/5, 10*t)
            j0 = e4*sqrt(mpf(15))*t/pi*S.qpow(-mpf(19)/60)*S.ppow(mpf(1)/24)*L(mpf(1)/5, 10*t)
            print(f"      the integrals themselves: |J(1/b)/closed form - 1| = {nstr(abs(S.Ju(1/b)/j1 - 1), 3)},  |J(b + 1/b)/closed form - 1| = {nstr(abs(S.Ju(b + 1/b)/j0 - 1), 3)}")
        if "d" in parts:
            al = 2*t; Q = exp(-al)                     # alpha = 2t: nome q^2, and q_1^4 = p
            c0, c1 = chis(Q); d0, d1 = chis(p)
            A = sqrt(pi*(5 - s5)/(5*al)); B = sqrt(pi*(5 + s5)/(5*al)); K = sqrt(135*al/(2*pi))
            l0 = -A*S.ppow(-mpf(1)/120)*(d0 - 2) - B*S.ppow(mpf(71)/120)*d1 - K*L(mpf(1)/5, 5*al)
            l1 = -B*S.ppow(-mpf(1)/120)*(d0 - 2) + A*S.ppow(mpf(71)/120)*d1 - K*L(mpf(2)/5, 5*al)
            print(f"  (d) quoted chi law: |(6a)| = {nstr(abs(exp(al/120)*(c0 - 2) - l0), 3)},  |(6b)| = {nstr(abs(exp(-71*al/120)*c1 - l1), 3)}")
        if "e" in parts:
            r = S.ppow(mpf(1)/2); r1 = S.ppow(-mpf(1)/120); r11 = S.ppow(mpf(11)/120)      # r^(-1/60), r^(11/60)
            G4, H4 = RR(q**4); Gm, Hm = RR(-q); Gr, Hr = RR(r); Gmr, Hmr = RR(-r)
            d = [abs(G4 - S.qpow(mpf(1)/15)*k*(s2*r1*Gr + s1*r11*Hr)), abs(H4 - S.qpow(-mpf(11)/15)*k*(s1*r1*Gr - s2*r11*Hr)),
                 abs(Gm - S.qpow(mpf(1)/60)*k*(s1*r1*Gmr + s2*r11*Hmr)), abs(Hm - S.qpow(-mpf(11)/60)*k*(s2*r1*Gmr - s1*r11*Hmr)),
                 abs(qp(q, q*q)*qp(p, p*p) - sqrt(2)*S.qpow(mpf(1)/24)*S.ppow(mpf(1)/24))]
            z = exp(2*pi*I/60); Sm = [[k*s2, k*s1], [k*s1, -k*s2]]; D2 = [z**-2, z**22]
            SDS = [[sum(Sm[i][m]*D2[m]*Sm[m][j] for m in range(2)) for j in range(2)] for i in range(2)]
            Nm = [[z**-1*SDS[0][0], z**-1*SDS[0][1]], [z**11*SDS[1][0], z**11*SDS[1][1]]]
            Ne = [[k*s1, k*s2*z**-6], [k*s2*z**6, -k*s1]]
            dm = max(abs(Nm[i][j] - Ne[i][j]) for i in range(2) for j in range(2))
            N0 = G4*Hmr + S.qpow(mpf(1)/4)*Hm*Hr; N1 = G4*Gmr - S.qpow(mpf(1)/4)*Hm*Gr
            Nb = Gm*Hr + S.qpow(mpf(3)/4)*H4*Hmr; Nb1 = Gm*Gr - S.qpow(mpf(3)/4)*H4*Gmr
            den = qp(q, q*q)*qp(-p, p)
            dD = [abs(S.D(u)*den/v - 1) for u, v in ((0, N0), (1/b, N1), (b, Nb), (b + 1/b, Nb1))]
            print("  (e) the four bilinear forms of D(u) (step 1 of the proof): " + ", ".join(nstr(v, 3) for v in dD) + f";   matrix T S T^2 S: {nstr(dm, 3)}")
            print("      G(q^4), H(q^4), G(-q), H(-q) in terms of G(+-r), H(+-r):  " + ", ".join(nstr(v, 3) for v in d[:4]) + f";   (q;q^2)(p;p^2) = sqrt 2 (qp)^(1/24): {nstr(d[4], 3)}")
