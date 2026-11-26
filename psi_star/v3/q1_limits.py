"""Q1: leading behaviour of psi*(e^(2 pi i / X)) for X -> infinity with X = rho mod 2 fixed (rho rational).
L_j(rho) := lim  psi_j*(e^(2 pi i/X)) * sqrt(t) * e^(Sbar/t),   t = -2 pi i / X,  Sbar = pi^2/3 + i Cl2(pi/3).
Neville extrapolation in h = 1/X.  Output: q1_limits.json"""
import mpmath as mp, sys, json
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from jval import J
from math import gcd
mp.mp.dps = 40
Cl2 = mp.clsin(2, mp.pi/3); Sb = mp.pi**2/3 + 1j*Cl2
Ms = [120, 160, 200, 260, 320, 400, 520, 640, 800, 1000]
def neville(hs, vs):
    n = len(hs); T = vs[:]
    for k in range(1, n):
        T = [(hs[i+k]*T[i] - hs[i]*T[i+1])/(hs[i+k]-hs[i]) for i in range(n-k)]
    return T[0]
def limits(r, s, odd):
    hs, v1, v0 = [], [], []
    for M in Ms:
        M += odd
        num, den = s, s*M + r                      # 1/X = s/(sM+r)
        j1, j0 = J(num, den)
        X = mp.mpf(den)/s; t = -2j*mp.pi/X
        lead = mp.exp(-Sb/t)/mp.sqrt(t)
        hs.append(1/X); v1.append(j1/lead); v0.append(j0/lead)
    a1, a0 = neville(hs, v1), neville(hs, v0)
    b1, b0 = neville(hs[1:], v1[1:]), neville(hs[1:], v0[1:])
    return a1, a0, max(abs(a1-b1), abs(a0-b0))
if __name__ == "__main__":
    smax = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    res = {}
    for s in range(1, smax+1):
        for r in range(s):
            if gcd(r, s) != 1: continue
            for odd in (0, 1):
                a1, a0, err = limits(r, s, odd)
                res[f"{r + odd*s}/{s}"] = [str(a1.real), str(a1.imag), str(a0.real), str(a0.imag), float(err)]
    json.dump(res, open("q1_limits.json", "w"))
    print(len(res), "values of rho; largest extrapolation uncertainty", max(v[4] for v in res.values()))
