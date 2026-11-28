"""Addendum v30, Theorem 20: the matrix of state integrals S(a, x) = Jt(zetat; p)^-1 Delta(a, x) J(zeta; q).
   (a) Jt S = Delta J with S by quadrature (30 digits), det S, det J, P elementary.
   (b) the Casoratian of H (v6) and the equation of H, on both sides.
   (c) the Corollary: the dual trivial series from F(z), F(zq) and the six state integrals."""
import sys
from mpmath import mp, mpf, mpc, exp, pi, sqrt, matrix, nstr, det
from common import Setup, I
from frames import *
mp.dps = 30
print("(a) Theorem 20: Jt(zetat;p) S(a,x) = Delta(a,x) J(zeta;q), S by quadrature; det S; P elementary")
for x in (mpc("0.3", "0.2"), mpc("0.45", "0.3"), mpc("-0.2", "0.6")):
    S = Setup(x)
    print(f"  x = {nstr(x, 6)}:  |P - P_elem|/|P| = {nstr(abs(S.P() - Pelem(S))/abs(S.P()), 3)}")
    for a in (mpf(1)/4, mpf(1)/2, mpc("0.3", "0.1"), mpf("0.8"), mpf("1.7")):
        SM = state_matrix(S, a); J = frame_q(S, a); Jt = frame_dual(S, a); D = Delta(S, a)
        lhs = Jt*SM; rhs = D*J
        dS = det(SM); dform = -I/(2*x)*S.qpow(-mpf(7)/6)*S.ppow(mpf(7)/12)*exp(-3*pi*I*a + pi*I*a/x - pi*I*a*a/x)
        print(f"     a = {nstr(a, 4):<14s} max|Jt S - Delta J| = {nstr(mnorm(lhs - rhs), 3)} (max entry {nstr(mnorm(rhs), 3)});  |det S/formula - 1| = {nstr(abs(dS/dform - 1), 3)};"
              f"  det J q^(3/4) zeta/(-2) - 1 = {nstr(abs(det(J)*S.qpow(mpf(3)/4)*exp(pi*I*a)/(-2) - 1), 3)}")
print("(b) H_e(y) H_e(y r^(1/2)) - H_o(y) H_o(y r^(1/2)) = 1 and H(y) - H(yr) = y r^(1/4) H(y r^(1/2)), for r = q and r = p^2")
from frames import Hq, Hp
for x in (mpc("0.3", "0.2"), mpc("-0.2", "0.6")):
    S = Setup(x)
    for y in (mpc("0.7", "0.4"), mpc("-1.3", "0.2")):
        out = []
        for Hf, r, rh, r4 in ((Hq, S.q, S.qpow(mpf(1)/2), S.qpow(mpf(1)/4)), (Hp, S.p**2, S.p, S.ppow(mpf(1)/2))):
            e = lambda v: (Hf(S, v) + Hf(S, -v))/2
            o = lambda v: (Hf(S, v) - Hf(S, -v))/2
            out.append(nstr(abs(e(y)*e(y*rh) - o(y)*o(y*rh) - 1), 3) + ", " + nstr(abs(Hf(S, y) - Hf(S, y*r) - y*r4*Hf(S, y*rh)), 3))
        print(f"   x = {nstr(x, 5)}, y = {nstr(y, 4)}:  r = q: {out[0]};   r = p^2: {out[1]}")
print("(c) Corollary: W(T) = {[F(z) - Jc(a)] i0(a+x) - [F(zq) - Jc(a+x)] i0(a)}/det S,  pT W(p^2 T) = {[F(z) - Jc(a)] i1(a+x) - [F(zq) - Jc(a+x)] i1(a)}/det S")
from a1_formula import Wt
for x, a in ((mpc("0.3", "0.2"), mpf(1)/4), (mpc("0.45", "0.3"), mpc("0.3", "0.1")), (mpc("-0.2", "0.6"), mpf("0.8"))):
    S = Setup(x); SM = state_matrix(S, a); z = exp(2*pi*I*a); p = S.p; T = p*exp(2*pi*I*a/x)
    F0 = S.Psi(z); F1 = S.Psi(z*S.q); d = det(SM)
    w1 = ((F0 - SM[0, 1])*SM[1, 2] - (F1 - SM[0, 2])*SM[1, 1])/d
    w2 = ((F0 - SM[0, 1])*SM[2, 2] - (F1 - SM[0, 2])*SM[2, 1])/d
    print(f"   x = {nstr(x, 5)}, a = {nstr(a, 4)}: W(T) = {nstr(Wt(S, T), 12)}, deviation {nstr(abs(w1 - Wt(S, T)), 3)};   pT W(p^2 T): deviation {nstr(abs(w2 - p*T*Wt(S, p*p*T)), 3)}")
