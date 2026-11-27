"""Addendum v18, comparison with census3/zfamily of the repository (commit bbac2c2, 27 Nov 2026).
There: V(a, x) = 2 pi i a x - x^2/2 + Li2(-e^(-x)) + pi^2/12,  x = -log(w) + 2 pi i k,  w^2 + w = 1/z;  selected (root sign, k) per z
in volumes.py:  a = 1/2: (-1, -3),  1/4: (+1, 0),  5/6: (-1, 0),  1/3: (-1, +1),  0: (+1, 0)   (root = (-1 + sgn sqrt(1 + 4/z))/2).
Here: W = their w, my w = -x, V_k'(w) = -w^2/2 - 2 pi i (a + k') w + Li2(-W).
 (1) each selected point translated: half-plane of the root, the branch k' for which it is a true critical point, signed Im V, D(-W);
 (2) the same list of |Im V(a, x)| for k = -3..3 as branches.py uses, to show which candidates are distinguishable by |Im V| alone;
 (3) |psi1*(e^(2 pi i/X))| / J_41(X/2) * X/2 -> 1  (the figure-eight comparison of psi_star/knot), and the constant 3^(-1/4) sqrt(X/2) from v18.
Usage: python3 d1_compare.py"""
from mpmath import mp, mpf, mpc, exp, pi, sqrt, log, polylog, nstr, im, re, arg, expjpi
from b1_roots import F_root, rule, bloch_wigner
mp.dps = 30
I = mpc(0, 1)
sel = [(mpf(1)/2, -1, -3), (mpf(1)/4, 1, 0), (mpf(5)/6, -1, 0), (mpf(1)/3, -1, 1), (mpf(0), 1, 0)]
print("(1) the selected points of census3/zfamily/volumes.py in the notation of v18")
for a, sgn, k in sel:
    Wt = (-1 + sgn*sqrt(1 + 4*expjpi(-2*a)))/2
    x = -log(Wt) + 2*pi*I*k
    Vt = 2*pi*I*a*x - x*x/2 + polylog(2, -exp(-x)) + pi**2/12
    w = -x
    kp = (-w - log(1 + Wt))/(2*pi*I) - a                  # branch on which (W, w) is a critical point
    Wr, wr, kr, _, Vr = rule(a) if a != 0 else (None, None, None, None, None)
    same = "-" if a == 0 else ("yes" if abs(Wt - Wr) < mpf(10)**(-20) and abs(w - wr) < mpf(10)**(-20) else ("same root, w shifted by %s * 2 pi i" % nstr((w - wr)/(2*pi*I), 3) if abs(Wt - Wr) < mpf(10)**(-20) else "other root"))
    print(f"    a = {nstr(a, 4):8s} W = {nstr(Wt, 8)} ({'upper' if im(Wt) > 0 else ('lower' if im(Wt) < 0 else 'real')}),  Im w/pi = {nstr(im(w)/pi, 5)},  critical on branch k' = {nstr(kp, 5)},  Im V = {nstr(im(Vt), 12)},  D(-W) = {nstr(bloch_wigner(-Wt), 12)},  saddle of the v18 rule: {same}")
print("(2) |Im V(a, x)|/(2 pi) for k = -3..3 (their fixed-a potential), root +1 | root -1;  * marks true critical points (k' = 0)")
for a, sgn0, k0 in sel[:4]:
    row = []
    for sgn in (1, -1):
        Wt = (-1 + sgn*sqrt(1 + 4*expjpi(-2*a)))/2; part = []
        for k in range(-3, 4):
            x = -log(Wt) + 2*pi*I*k
            Vt = 2*pi*I*a*x - x*x/2 + polylog(2, -exp(-x)) + pi**2/12
            crit = abs((x - log(1 + Wt))/(2*pi*I) - a) < mpf(10)**(-20)     # stationarity of V(a, .) with j = 0
            part.append(nstr(abs(im(Vt))/(2*pi), 6) + ("*" if crit else ""))
        row.append(" ".join(part))
    print(f"    a = {nstr(a, 4):8s} " + "  |  ".join(row))
print("(3) figure-eight comparison:  J_41(N) = sum_{k<N} prod_{j<=k} |1 - e^(2 pi i j/N)|^2")
for X in (802, 1602, 3202):
    N = X//2
    with mp.workdps(30 + int(0.15*X)):
        zeta = expjpi(mpf(2)/N); s = mpf(0); pr = mpf(1)
        for kk in range(N):
            s += pr; pr *= abs(1 - zeta**(kk + 1))**2
        Fv = abs(F_root(mpf(1)/2, X))
        ratio = +(Fv/s*X/2)
        c = +(Fv/(mpf(3)**(-mpf(1)/4)*sqrt(mpf(X)/2)*exp(bloch_wigner(expjpi(mpf(1)/3))*X/(2*pi))))
    print(f"    X = {X}: |psi1*|/J_41(X/2) * X/2 = {nstr(ratio, 10)};   |psi1*| / (3^(-1/4) sqrt(X/2) e^(Cl2 X/2pi)) = {nstr(c, 10)}")
