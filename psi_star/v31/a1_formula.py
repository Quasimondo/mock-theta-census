"""Addendum v18, task 1: the exact formula for F(z) = sum z^n q^(n(n+1)/2) (-q;q)_n from Theorems A-C of v12 at u = a/b.
   z = e^(2 pi i a),  zt = e^(2 pi i a/x),  t = p zt,
   (Z-up)   F(z)      =  i P e^(-pi i a) [ J(a/b) + I(a/b) (W(t/p^2) - 1) - (t/p) I((a+1)/b) W(t) ]
   (Z-down) Fc(1/z)   = -2 i P e^(pi i a) [ J(a/b) - I(a/b) G0(t) + I((a+1)/b) G1(t) ],     Fc(v) = sum q^(n^2+n) v^n/(-q;q)_n
   I(u) = c e^(-pi i u^2/2) D(u),   D((a)/b) = H_even(zeta;q) H(-zetat p; p^2) + H_odd(zeta;q) H(zetat p; p^2),  zeta = e^(pi i a), zetat = e^(pi i a/x)
Usage: python3 a1_formula.py [digits]"""
import sys
from mpmath import mp, mpf, mpc, exp, pi, sqrt, nstr
from common import Setup, I

def Wt(S, t):                 # sum (-t)^k p^(k^2+k) (p;p^2)_k
    p = S.p; s = mpc(0); po = mpc(1); k = 0
    while True:
        term = (-t)**k*p**(k*k + k)*po
        s += term
        if abs(term) < S.tiny and k > 3: break
        po *= 1 - p**(2*k + 1); k += 1
    return s

def G01(S, t):                # G0 = sum t^-k p^(2k^2+k)/(p;p^2)_k,  G1 = sum t^-k p^(2k^2+3k)/(p;p^2)_(k+1)
    p = S.p; g0 = g1 = mpc(0); po = mpc(1); k = 0
    while True:
        a0 = t**(-k)*p**(2*k*k + k)/po; po1 = po*(1 - p**(2*k + 1)); a1 = t**(-k)*p**(2*k*k + 3*k)/po1
        g0 += a0; g1 += a1
        if abs(a0) < S.tiny and abs(a1) < S.tiny and k > 3: break
        po = po1; k += 1
    return g0, g1

def Hfun(zz, r, half):        # H(zz; r) = sum zz^m r^(m^2/4)/(r;r)_m, with r^(1/4) = half given explicitly
    s = mpc(0); den = mpc(1); m = 0
    while True:
        term = zz**m*half**(m*m)/den
        s += term
        if abs(term) < mpf(10)**(-mp.dps - 8) and m > 4: break
        m += 1; den *= 1 - r**m
    return s

def Iser(S, u):               # Theorem C
    return S.cconst()*exp(-pi*I*u*u/2)*S.D(u)

def check(x, a, label):
    S = Setup(x); b = S.b; p = S.p; u = a/b
    z = exp(2*pi*I*a); t = p*exp(2*pi*I*a/x)
    J = S.Ju(u); I0 = Iser(S, u); I1 = Iser(S, u + 1/b)
    up = I*S.P()*exp(-pi*I*a)*(J + I0*(Wt(S, t/p**2) - 1) - t/p*I1*Wt(S, t))
    g0, g1 = G01(S, t)
    dn = -2*I*S.P()*exp(pi*I*a)*(J - I0*g0 + I1*g1)
    F = S.Psi(z); Fc = S.fstar(1/z)
    # D(u) through H
    ze = exp(pi*I*a); zt = exp(pi*I*a/x)
    He = (Hfun(ze, S.q, S.qpow(mpf(1)/4)) + Hfun(-ze, S.q, S.qpow(mpf(1)/4)))/2
    Ho = (Hfun(ze, S.q, S.qpow(mpf(1)/4)) - Hfun(-ze, S.q, S.qpow(mpf(1)/4)))/2
    Dh = He*Hfun(-zt*p, p*p, S.ppow(mpf(1)/2)) + Ho*Hfun(zt*p, p*p, S.ppow(mpf(1)/2))
    # quadrature of I(u) as an independent check of Theorem C at this u
    Iq = S.Iu(u)
    print(f"  {label:<26s} F(z) = {nstr(F, 18)}")
    print(f"      |F - (Z-up)| = {nstr(abs(F - up), 3)}   |Fc(1/z) - (Z-down)| = {nstr(abs(Fc - dn), 3)}   |D - D via H|/|D| = {nstr(abs(S.D(u) - Dh)/abs(Dh), 3)}   |I quad - I series|/|I| = {nstr(abs(Iq - I0)/abs(I0), 3)}")
    print(f"      sizes: |i P e^(-pi i a) J| = {nstr(abs(S.P()*exp(-pi*I*a)*J), 5)}, |t| = {nstr(abs(t), 5)}, |I(a/b)| = {nstr(abs(I0), 5)}, |I((a+1)/b)| = {nstr(abs(I1), 5)}")

if __name__ == "__main__":
    mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    cases = [(mpf(1)/4, "z = i (a = 1/4)"), (-mpf(3)/4, "z = i (a = -3/4)"), (mpf(5)/6, "z = -omega (a = 5/6)"), (-mpf(1)/6, "z = -omega (a = -1/6)"),
             (mpf(1)/2, "z = -1 (a = 1/2)"), (mpf(0), "z = 1 (a = 0)"), (mpf(1), "z = 1 (a = 1)"), (mpf(1)/3, "z = omega (a = 1/3)"),
             (mpc("0.3", "0.1"), "a = 0.3 + 0.1i")]
    for x in (mpc("0.3", "0.2"), 1/mpc(3, -mpf(2)/5), mpc("0.11", "0.03")):
        print(f"x = {nstr(x, 8)}   |q| = {nstr(abs(exp(2*pi*I*x)), 4)}, |p| = {nstr(abs(exp(-pi*I/x)), 4)}")
        for a, lab in cases: check(x, a, lab)
