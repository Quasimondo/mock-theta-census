"""Test 2 (addendum v12, Theorem A and the corollaries).
For general u, with I(u) evaluated by its series (Theorem C, tested in t1):
   upper route:  J(u) = -i b q^(1/8) e^(pi i u b) (p;p^2)_inf/(-q;q)_inf * Psi(e^(2 pi i u b))  -  sum_{m>=1} (-1)^m I(u - m/b)
   lower route:  J(u) = +i b q^(1/8) e^(-pi i u b) (p;p^2)_inf/(2(-q;q)_inf) * fstar(e^(-2 pi i u b))  +  sum_{m>=0} (-1)^m I(u + m/b)
with Psi(z) = sum q^(n(n+1)/2) z^n (-q;q)_n,  fstar(z) = sum q^(n^2+n) z^n/(-q;q)_n.
Corollaries at u = 1/(2b) and u = b + 1/(2b):
   psi1* = P J(1/(2b))           + (kappa/2) [ B2 p W3(p) - B4 W(p) ]
   psi0* = P q^(1/2) J(b+1/(2b)) - i q^(1/4) (kappa/2) [ B1 p W3(p) + B3 W(p) ]
   f1*   = 2 P J(1/(2b))           - kappa [ B2 F_b(p) - B4 F_a(p) ]
   f0*   = 2 P q^(1/2) J(b+1/(2b)) + i q^(1/4) kappa [ B1 F_b(p) + B3 F_a(p) ]
Usage: python3 t2_kernel.py [digits]"""
import sys
from mpmath import mp, mpf, mpc, exp, pi, qp, nstr
from common import Setup, POINTS, I


def Iser(S, u):
    return S.cconst()*exp(-pi*I*u*u/2)*S.D(u)


def dual_sum(S, u, sign):          # sign=-1: sum_{m>=1} (-1)^m I(u - m/b);  sign=+1: sum_{m>=0} (-1)^m I(u + m/b)
    tot = mpc(0); m = 1 if sign < 0 else 0; small = 0
    while small < 3:
        t = (-1)**m*Iser(S, u + sign*m/S.b); tot += t; m += 1
        small = small + 1 if abs(t) < S.tiny*(1 + abs(tot)) else 0
        if m > 400: raise RuntimeError("dual sum does not converge")
    return tot, m


if __name__ == "__main__":
    mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    for x in POINTS:
        S = Setup(x); b = S.b; q = S.q; p = S.p
        q8 = S.qpow(mpf(1)/8); A0 = qp(p, p*p)/qp(-q, q)
        print("x =", nstr(x, 8), "  |q| =", nstr(abs(q), 4), " |p| =", nstr(abs(p), 4))
        us = [("1/(2b)", 1/(2*b)), ("b+1/(2b)", b + 1/(2*b)), ("0.37-0.21i", mpc("0.37", "-0.21")), ("-0.3+0.5i", mpc("-0.3", "0.5"))]
        Jv = {}
        for name, u in us:
            J = S.Ju(u); Jv[name] = J
            up, mu = dual_sum(S, u, -1); lo, ml = dual_sum(S, u, +1)
            resU = -I*b*q8*exp(pi*I*u*b)*A0*S.Psi(exp(2*pi*I*u*b))
            resL = I*b*q8*exp(-pi*I*u*b)*A0/2*S.fstar(exp(-2*pi*I*u*b))
            print(f"   u = {name}:  J(u) = {nstr(J, 20)}")
            print(f"       upper route: |residues| = {nstr(abs(resU), 3)}, |dual sum| = {nstr(abs(up), 3)} ({mu} terms),  |J - (res - sum)| = {nstr(abs(J - (resU - up)), 3)}")
            print(f"       lower route: |residues| = {nstr(abs(resL), 3)}, |dual sum| = {nstr(abs(lo), 3)} ({ml} terms),  |J - (res + sum)| = {nstr(abs(J - (resL + lo)), 3)}")
            sys.stdout.flush()
        B = S.B(); Ep, Em, Epp, Emp, W, W3 = S.pside(); Fa, Fb = S.Fab()
        P = S.P(); kap = S.kappa(); q2 = S.qpow(mpf(1)/2); q4 = S.qpow(mpf(1)/4)
        psi1 = S.Psi(-1); psi0 = q*S.Psi(-q); f1 = S.fstar(-1)
        f0 = mpc(0); term = mpc(1); n = 0
        while abs(term) > S.tiny or n < 5:
            f0 += term; n += 1; term *= -q**(2*n - 1)/(1 + q**n)
        J1 = Jv["1/(2b)"]; J0 = Jv["b+1/(2b)"]
        print("   psi1* =", nstr(psi1, 20), "  |psi1* - rhs| =", nstr(abs(psi1 - (P*J1 + kap/2*(B[1]*p*W3 - B[3]*W))), 3))
        print("   psi0* =", nstr(psi0, 20), "  |psi0* - rhs| =", nstr(abs(psi0 - (P*q2*J0 - I*q4*kap/2*(B[0]*p*W3 + B[2]*W))), 3))
        print("   f1*   =", nstr(f1, 20), "  |f1* - rhs|   =", nstr(abs(f1 - (2*P*J1 - kap*(B[1]*Fb - B[3]*Fa))), 3))
        print("   f0*   =", nstr(f0, 20), "  |f0* - rhs|   =", nstr(abs(f0 - (2*P*q2*J0 + I*q4*kap*(B[0]*Fb + B[2]*Fa))), 3))
        sys.stdout.flush()
