"""Addendum v21: psi1*, psi0* at odd roots of unity q0 = e^(2 pi i/X).
 a  Theorem 11:  A = sum_{n < X} (-1)^n q0^(n(n+1)/2) (-q0;q0)_n  against  P [3 J - I(+) - I(-)]  at x = 1/X  (radial limit of psi1* is A/3);
    the same for psi0* (A0 with q0^((n+1)(n+2)/2), integrals at u = b + 1/(2b), b - 1/(2b) with the factor q^(1/2))
 b  Lemma R2': along x = 1/X + i eps (X = 5):  W(p) -> 1/3,  p W3(p) -> -1/3,  F_a(p) -> -1/3,  F_b(p) -> 1/3,  psi(p) P+-(p) -> 0,  and psi1*(q) -> A/3
 c  Theorem 12:  (A/3) / [ (1/3) 3^(-1/4) (X/2)^(1/2) e^(Cl2 X/2pi) e^(2 pi i/3 - i pi X/6) ] -> 1,  psi0* with the extra factor e^(2 pi i/3)
Usage: python3 k1_odd.py [a b c]"""
import sys
from mpmath import mp, mpf, mpc, exp, pi, sqrt, nstr, clsin, expjpi, qp
from h1_roots import PI_np, PI_w
I = mpc(0, 1)

def A_odd(X, case="psi1"):
    with mp.workdps(40 + int(0.11*X)):
        q = expjpi(mpf(2)/X); s = mpc(0); term = mpc(1); qn = mpc(1)
        for n in range(X):
            s += term*(q**(n + 1) if case == "psi0" else 1); qn *= q; term *= -qn*(1 + qn)
        return +s

def part_a():
    print("(a) Theorem 11 at x = 1/X, X odd: A against P [3 J - I(+) - I(-)]   (double precision)")
    for X in (5, 7, 9, 15, 25, 41, 61):
        row = []
        for case in ("psi1", "psi0"):
            ew = case == "psi0"
            J = PI_np(X, "J", ew); Ip = PI_np(X, "+", True) if (ew and X <= 8) else PI_w(X, +1, ew); Im = PI_w(X, -1, ew)
            A = complex(A_odd(X, case)); rhs = 3*J - Ip - Im
            row.append(f"{case}*: A = {A:.8f}, relative difference {abs(A - rhs)/abs(A):.1e}")
        print(f"    X = {X:2d}  " + ";   ".join(row))

def part_b():
    X = 5
    print(f"(b) Lemma R2' along x = 1/X + i eps, X = {X}")
    A3 = A_odd(X)/3
    for eps, dps in (("1e-2", 40), ("1e-3", 60), ("3e-4", 120), ("1e-4", 300)):
        mp.dps = dps
        e = mpf(eps); x = mpf(1)/X + I*e; q = exp(2*pi*I*x); p = exp(-pi*I/x); tiny = mpf(10)**(-dps - 5)
        W = W3 = mpc(0); po = mpc(1); j = 0
        while j < 20 or abs(p)**(j*j)*abs(po) > tiny:
            W += (-1)**j*p**(j*j + j)*po; W3 += (-1)**j*p**(j*j + 3*j)*po; po *= 1 - p**(2*j + 1); j += 1
        Fa = Fb = mpc(0); po = mpc(1); n = 1
        while n < 20 or abs(p)**(2*n*n - 3*n + 1)/abs(po) > tiny:
            po *= 1 - p**(2*n - 1); Fa += p**(2*n*n - n)/po; Fb += p**(2*n*n - 3*n + 1)/po; n += 1
        psi = sum(p**(k*(k + 1)//2) for k in range(int(sqrt(2*(dps + 10)*mp.log(10)/(1 - abs(p)))) + 5))
        Pp = Pm = mpc(0); den = mpc(1); k = 0
        while k < 20 or abs(p)**(2*k*k - k)/abs(den) > tiny:
            Pp += p**(2*k*k + k)/den; Pm += p**(2*k*k - k)/den; k += 1; den *= 1 - p**(2*k)
        s = mpc(0); term = mpc(1); n = 0
        while abs(term) > tiny or n < 4*X:
            s += term; n += 1; term *= -q**n*(1 + q**n)
        d = 1 - abs(p)
        print(f"    eps = {eps} (1 - |p| = {nstr(d, 3)}, arg(-p)/(1-|p|)^2 = {nstr(mp.arg(-p)/d**2, 3)}): |W - 1/3| = {nstr(abs(W - mpf(1)/3), 3)}, |p W3 + 1/3| = {nstr(abs(p*W3 + mpf(1)/3), 3)}, "
              f"|F_a + 1/3| = {nstr(abs(Fa + mpf(1)/3), 3)}, |F_b - 1/3| = {nstr(abs(Fb - mpf(1)/3), 3)}, |psi P+| = {nstr(abs(psi*Pp), 3)}, |psi P-| = {nstr(abs(psi*Pm), 3)}, "
              f"identities (V), (VI): {nstr(abs(W + Fa - psi*Pp), 2)}, {nstr(abs(p*W3 + Fb - psi*Pm), 2)};  |psi1*(q) - A/3| = {nstr(abs(s - A3), 3)}")

def part_c():
    mp.dps = 30
    Cl2 = clsin(2, pi/3)
    print("(c) Theorem 12: X (exact/asymptotic - 1);   -2 pi i (1/4 + i sqrt(3)/54) = " + nstr(-2*pi*I*(mpf(1)/4 + I*sqrt(3)/54), 6))
    for X in (101, 201, 401, 801, 1601, 3201):
        lead = mpf(3)**(-mpf(1)/4)*sqrt(mpf(X)/2)*exp(Cl2*X/(2*pi))*exp(2*pi*I/3 - I*pi*X/6)/3
        r1 = A_odd(X, "psi1")/3/lead - 1; r0 = A_odd(X, "psi0")/3/(lead*exp(2*pi*I/3)) - 1
        print(f"    X = {X:5d}: psi1*: {nstr(X*r1, 6)};   psi0*: {nstr(X*r0, 6)}")

if __name__ == "__main__":
    for p_ in (sys.argv[1:] or list("abc")): {"a": part_a, "b": part_b, "c": part_c}[p_]()
