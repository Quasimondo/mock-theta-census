"""Addendum v30: the q-side frame under Gamma_0(2).
   T: x -> x+1 (a fixed):                J(a; x+1) = C_T J(a; x)
   gamma = [[1,0],[2,1]]: x' = x/(1+2x), a' = a/(1+2x):   J(a'; x') = Lambda J(a; x) Omega,  Omega = S(a,x)^-1 S(a',x'),  Lambda = diag(1, lam, -lam)
   lam = (1+2x)^(-1/2) e^(pi i/12) e^(-5 pi i x x'/6) e^(-pi i (a - a')) e^(-pi i a a')"""
import sys
from mpmath import mp, mpf, mpc, exp, pi, sqrt, matrix, nstr, det
from common import Setup, I
from frames import *
mp.dps = 30
CT = matrix([[1, 0, 0], [0, (1 + I)/2, (1 - I)/2], [0, (1 - I)/2, (1 + I)/2]])

def lam(a, x):
    xp = x/(1 + 2*x); ap = a/(1 + 2*x)
    return exp(I*pi/12)*exp(-5*pi*I*x*xp/6)*exp(-pi*I*(a - ap) - pi*I*a*ap)/sqrt(1 + 2*x)

def Lam(a, x):
    l = lam(a, x); return matrix([[1, 0, 0], [0, l, 0], [0, 0, -l]])

def Omega(a, x):
    xp = x/(1 + 2*x); ap = a/(1 + 2*x)
    return state_matrix(Setup(x), a)**-1*state_matrix(Setup(xp), ap)

if __name__ == "__main__":
    print("(a) T:  J(a; x+1) = C_T J(a; x)")
    for x, a in ((mpc("0.3", "0.2"), mpf(1)/4), (mpc("-0.2", "0.6"), mpc("0.3", "0.1"))):
        print(f"   x = {nstr(x, 5)}, a = {nstr(a, 4)}: max deviation {nstr(mnorm(frame_q(Setup(x + 1), a) - CT*frame_q(Setup(x), a)), 3)}")
    print("(b) gamma:  J(a'; x') = Lambda J(a; x) Omega,  Omega from four columns of state integrals (quadrature)")
    for x, a in ((mpc("0.3", "0.2"), mpf(1)/4), (mpc("0.45", "0.3"), mpc("0.3", "0.1")), (mpc("-0.2", "0.6"), mpf("0.8")), (mpc("0.9", "0.5"), mpf(1)/2)):
        xp = x/(1 + 2*x); ap = a/(1 + 2*x); S = Setup(x); Sp = Setup(xp)
        d1 = delta(S, a)/delta(Sp, ap)
        lhs = frame_q(Sp, ap); rhs = Lam(a, x)*frame_q(S, a)*Omega(a, x)
        print(f"   x = {nstr(x, 5)}, a = {nstr(a, 4)}: x' = {nstr(xp, 5)}; |lam/(delta/delta') - 1| = {nstr(abs(lam(a, x)/d1 - 1), 3)};  max|J' - Lambda J Omega| = {nstr(mnorm(lhs - rhs), 3)} (max entry {nstr(mnorm(lhs), 3)})")
    print("(c) the relation E^2 = -1, E = gamma T^-1:  J(-a; x) = Lam(a1,x1-1) C_T^-1 Lam(a,x-1) C_T^-1 J(a; x) Omega(a,x-1) Omega(a1,x1-1),  x1 = (x-1)/(2x-1), a1 = a/(2x-1)")
    for x, a in ((mpc("0.3", "0.4"), mpf(1)/4), (mpc("0.6", "0.5"), mpc("0.3", "0.1"))):
        x1 = (x - 1)/(2*x - 1); a1 = a/(2*x - 1); Ci = CT**-1
        lhs = frame_q(Setup(x), -a)
        Om = Omega(a, x - 1)*Omega(a1, x1 - 1)
        rhs = Lam(a1, x1 - 1)*Ci*Lam(a, x - 1)*Ci*frame_q(Setup(x), a)*Om
        print(f"   x = {nstr(x, 5)}, a = {nstr(a, 4)}: max deviation {nstr(mnorm(lhs - rhs), 3)} (max entry {nstr(mnorm(lhs), 3)});  first row of the product of the two Omegas: {[nstr(Om[0, j], 8) for j in range(3)]}")
        Fm = Setup(x).Psi(exp(-2*pi*I*a)); J = frame_q(Setup(x), a)
        print(f"        F(1/z) = {nstr(Fm, 15)},  (1, F(z), F(zq)) . column 2 of the product = {nstr(Om[0, 1] + J[0, 1]*Om[1, 1] + J[0, 2]*Om[2, 1], 15)}")
