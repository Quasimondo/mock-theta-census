"""Q1, second sector (exploratory): subtract the growing sector from psi1*(e^(2 pi i/X)) and look at what is left.
growing sector:  Cbar t^(-1/2) e^(-Sbar/t) W(e^(-pi i X)) * conj-series( sum s_j t^j ),  s_j from the large-order fit (Q2)."""
import json, mpmath as mp, sys
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from jval import J
from wval import W
from fractions import Fraction as Fr
mp.mp.dps = 120
half = mp.mpf(1)/2
D = json.load(open("q2_series.json")); K = D["K"]
a1 = [mp.mpf(r) for r, i, _ in D["psi1@1"]]
Cl2 = mp.clsin(2, mp.pi/3); S = mp.pi**2/3 - 1j*Cl2; Sb = mp.conj(S)
Cb = mp.sqrt(mp.pi/mp.sqrt(3))*mp.expjpi(mp.mpf(5)/12)
def fit_one(a, S, J, k0):
    rows, rhs = [], []
    for k in range(k0, K+1):
        sc = mp.gamma(k+half) * abs(S)**(-(k+half))
        row = []
        for j in range(J):
            g = mp.gamma(k+half-j) * S**(-(k+half-j)) / sc
            row += [2*g.real, -2*g.imag]
        rows.append(row); rhs.append(a[k]/sc)
    sol, _ = mp.qr_solve(mp.matrix(rows), mp.matrix(rhs))
    return [mp.mpc(sol[2*j], sol[2*j+1]) for j in range(J)]
cA = fit_one(a1, S, 18, 60); cB = fit_one(a1, S, 16, 70)
s = [cA[j]/cA[0] for j in range(18)]
print("attached series coefficients s_j (psi1*), and agreement between two fits:")
for j in range(12):
    print("  s_%d = %s   (diff %s)" % (j, mp.nstr(s[j], 12), mp.nstr(abs(s[j] - cB[j]/cB[0]), 2)))
def pert(t, kmax):          # optimal truncation of sum a_k t^k
    tot = 0; best = None
    for k in range(kmax):
        term = a1[k]*t**k
        tot += term
    return tot
for rho in (Fr(0), Fr(1), Fr(1, 2), Fr(1, 3), Fr(2, 3)):
    print("rho =", rho)
    for M in (30, 40, 50, 60, 80):
        X = M + rho
        num, den = X.denominator, X.numerator
        j1, j0 = J(num, den)
        t = -2j*mp.pi*X.denominator/X.numerator
        xi = -rho/2; Wv = W(xi.numerator, xi.denominator)
        outs = []
        for jm in (6, 8, 10):
            grow = Cb*mp.exp(-Sb/t)/mp.sqrt(t)*Wv*sum(mp.conj(s[j])*t**j for j in range(jm))
            outs.append(j1 - grow)
        # perturbative series at this t, truncated at its smallest term
        terms = [a1[k]*t**k for k in range(60)]
        kopt = min(range(5, 60), key=lambda k: abs(terms[k]))
        P = sum(terms[:kopt])
        print("   X = %-6s remainder (6, 8, 10 terms of the growing series): %s | %s | %s    perturbative series: %s" % (X, mp.nstr(outs[0], 8), mp.nstr(outs[1], 8), mp.nstr(outs[2], 8), mp.nstr(P, 8)))
