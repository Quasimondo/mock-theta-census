"""Test 3 (addendum v12): the individual steps of the proof of Theorem C, and the dual sums of Lemma B.
 (a) inversion  Phi_b(y) Phi_b(-y) = Phi_b(0)^2 e^(pi i y^2),  Phi_b(0)^2 = q^(1/24) qt^(-1/24);  (q)_inf/(qt)_inf = q^(-1/24) qt^(1/24)/sqrt(-i x)
 (b) F(w) = int_C Phi_b(y) e^(2 pi i w y) dy = -i b (q)_inf/(qt)_inf e^(2 pi i c_b w) / Phi_b(-w - c_b)
     G(w) = int_L e^(2 pi i w y)/Phi_b(y) dy = (i/b) (qt)_inf/(q)_inf e^(-2 pi i c_b w) Phi_b(w + c_b)
 (c) K(lam) = int_Sigma Phi_b(z) e^(pi i z^2 + 2 pi lam z) dz,  Sigma = e^(i delta) R:
     I(u)      = e^(-i pi/4) (i/b) (qt)_inf/(q)_inf e^(pi i (i u - c_b)^2 + 2 pi i c_b^2) K(-u - 2 i c_b)
     I_half(v) = sqrt(2) e^(-i pi/4) (-i b) (q)_inf/(qt)_inf Phi_b(0)^(-2) e^(2 pi i (i v - c_b)^2 - 2 pi i c_b^2) K(-2 v - i c_b)
     where I_half(v) = int_C Phi_b(s) e^(-pi i s^2/2 + 2 pi v s) ds  (the integral of v7).
 (d) Lemma B for general u, with I by its series:  t = p e^(2 pi i u/b),  Wt(t) = sum_k (-t)^k p^(k^2+k) (p;p^2)_k,
     sum_{m>=1} (-1)^m I(u - m/b) = I(u) [Wt(t/p^2) - 1] - (t/p) I(u + 1/b) Wt(t);     Wt(1) = W, Wt(p^2) = W3, Wt(p^-2) = 1 - W + p W3.
     lower route:  sum_{m>=0} (-1)^m I(u + m/b) = I(u) G0(t) - I(u + 1/b) G1(t),
         G0(t) = sum_k t^(-k) p^(2k^2+k)/(p;p^2)_k,   G1(t) = sum_k t^(-k) p^(2k^2+3k)/(p;p^2)_(k+1);   G0(1) = F_b - F_a,  G1(1) = F_a/p.
Usage: python3 t3_steps.py [digits]"""
import sys
from mpmath import mp, mpf, mpc, exp, pi, sqrt, qp, quad, nstr
from common import Setup, POINTS, I


def Iser(S, u):
    return S.cconst()*exp(-pi*I*u*u/2)*S.D(u)


def Wt(S, t):
    p = S.p; tot = mpc(0); po = mpc(1); k = 0
    while True:
        a = (-t)**k*p**(k*k + k)*po; tot += a
        if abs(a) < S.tiny and k > 6: break
        po *= 1 - p**(2*k + 1); k += 1
    return tot


def G01(S, t):
    p = S.p; g0 = g1 = mpc(0); po = mpc(1); k = 0
    while True:
        a = t**(-k)*p**(2*k*k + k)/po; po *= 1 - p**(2*k + 1)
        c = t**(-k)*p**(2*k*k + 3*k)/po
        g0 += a; g1 += c
        if abs(a) < S.tiny and abs(c) < S.tiny and k > 6: break
        k += 1
    return g0, g1


if __name__ == "__main__":
    mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    for x in POINTS[:3]:
        S = Setup(x); b = S.b; q = S.q; p = S.p; qt = p*p
        cb = I*(b + 1/b)/2; delta = mpf(2)/5
        print("x =", nstr(x, 8))
        # (a)
        Phi0sq = S.qpow(mpf(1)/24)*S.ppow(-mpf(1)/12)
        y = mpc("0.31", "0.17")
        print("   (a) |Phi(y)Phi(-y)/(Phi(0)^2 e^(pi i y^2)) - 1| =", nstr(abs(S.Phi(y)*S.Phi(-y)/(Phi0sq*exp(pi*I*y*y)) - 1), 3),
              "   |Phi(0)^2/(q^(1/24) qt^(-1/24)) - 1| =", nstr(abs(S.Phi(0)**2/Phi0sq - 1), 3),
              "   |(q)/(qt) sqrt(-ix) q^(1/24) qt^(-1/24) - 1| =", nstr(abs(qp(q, q)/qp(qt, qt)*sqrt(-I*x)/Phi0sq**(-1) - 1), 3))
        # (b)
        eL = exp(-I*delta); eR = exp(I*delta)
        def onC(f, Lm=10):
            return quad(lambda s: f(s*eL)*eL, [-Lm, -6, -4, -2.5, -1.2, 0]) + quad(lambda s: f(s*eR)*eR, [0, 1.2, 2.5, 4, 6, Lm])
        w = exp(I*delta)*mpc("0.4", "-0.9")
        F = onC(lambda yy: S.Phi(yy)*exp(2*pi*I*w*yy), Lm=14)
        Fth = -I*b*qp(q, q)/qp(qt, qt)*exp(2*pi*I*cb*w)/S.Phi(-w - cb)
        G = S.line(lambda yy: exp(2*pi*I*w*yy)/S.Phi(yy), Lmax=14)
        Gth = I/b*qp(qt, qt)/qp(q, q)*exp(-2*pi*I*cb*w)*S.Phi(w + cb)
        print("   (b) w =", nstr(w, 6), "  |F/F_formula - 1| =", nstr(abs(F/Fth - 1), 3), "   |G/G_formula - 1| =", nstr(abs(G/Gth - 1), 3))
        # (c)
        eS = exp(I*delta)
        def K(lam):
            return quad(lambda s: S.Phi(s*eS)*exp(pi*I*(s*eS)**2 + 2*pi*lam*s*eS)*eS, [-10, -6, -4, -2.5, -1.2, 0, 1.2, 2.5, 4, 6, 10])
        for u in (1/(2*b), mpc("0.37", "-0.21")):
            v = u/2 + I*cb/2
            Iu = S.Iu(u)
            Ih = onC(lambda s: S.Phi(s)*exp(-pi*I*s*s/2 + 2*pi*v*s))
            Kv = K(-u - 2*I*cb)
            r1 = exp(-I*pi/4)*(I/b)*qp(qt, qt)/qp(q, q)*exp(pi*I*(I*u - cb)**2 + 2*pi*I*cb*cb)*Kv
            r2 = sqrt(2)*exp(-I*pi/4)*(-I*b)*qp(q, q)/qp(qt, qt)/Phi0sq*exp(2*pi*I*(I*v - cb)**2 - 2*pi*I*cb*cb)*Kv
            print("   (c) u =", nstr(u, 6), "  |I(u)/(c1 K) - 1| =", nstr(abs(Iu/r1 - 1), 3), "   |I_half(u/2 + i c_b/2)/(c2 K) - 1| =", nstr(abs(Ih/r2 - 1), 3))
        # (d)
        print("   (d) |Wt(1)-W| =", nstr(abs(Wt(S, 1) - S.pside()[4]), 3), "  |Wt(p^2)-W3| =", nstr(abs(Wt(S, p*p) - S.pside()[5]), 3),
              "  |Wt(p^-2) - (1 - W + p W3)| =", nstr(abs(Wt(S, 1/(p*p)) - (1 - S.pside()[4] + p*S.pside()[5])), 3))
        Fa, Fb = S.Fab(); g0, g1 = G01(S, 1)
        print("       |G0(1) - (F_b - F_a)| =", nstr(abs(g0 - (Fb - Fa)), 3), "  |G1(1) - F_a/p| =", nstr(abs(g1 - Fa/p), 3))
        for u in (1/(2*b), b + 1/(2*b), mpc("0.37", "-0.21"), mpc("-0.3", "0.5")):
            t = p*exp(2*pi*I*u/b)
            lhs = mpc(0); m = 1; small = 0
            while small < 3:
                a = (-1)**m*Iser(S, u - m/b); lhs += a; m += 1
                small = small + 1 if abs(a) < S.tiny*(1 + abs(lhs)) else 0
            rhs = Iser(S, u)*(Wt(S, t/(p*p)) - 1) - t/p*Iser(S, u + 1/b)*Wt(S, t)
            rec = max(abs(Iser(S, u - mm/b) - t*p**(2*mm - 2)*Iser(S, u - (mm - 1)/b) - t*p**(2*mm - 3)*Iser(S, u - (mm - 2)/b))/abs(Iser(S, u - mm/b)) for mm in (0, 1, 2, 3))
            lo = mpc(0); m = 0; small = 0
            while small < 3:
                a = (-1)**m*Iser(S, u + m/b); lo += a; m += 1
                small = small + 1 if abs(a) < S.tiny*(1 + abs(lo)) else 0
            g0, g1 = G01(S, t)
            lo_rhs = Iser(S, u)*g0 - Iser(S, u + 1/b)*g1
            print("       u =", nstr(u, 6), "  lower: |sum_m (-1)^m I(u+m/b)| =", nstr(abs(lo), 3), "  |lhs - rhs| =", nstr(abs(lo - lo_rhs), 3))
            print("       u =", nstr(u, 6), "  |sum_m (-1)^m I(u-m/b)| =", nstr(abs(lhs), 3), "  |lhs - rhs| =", nstr(abs(lhs - rhs), 3), "   recursion (series), max rel. defect m=0..3:", nstr(rec, 3))
        sys.stdout.flush()
