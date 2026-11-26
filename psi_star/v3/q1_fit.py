"""Q1 (negative result): is L_j(rho) a constant linear combination of 1, psi1*, psi0* at a related point?
Needs q1_limits.json (python3 q1_limits.py 9)."""
import mpmath as mp, sys, json
from fractions import Fraction as Fr
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from jval import J
mp.mp.dps = 40
D = json.load(open("q1_limits.json"))
rhos = [Fr(k) for k in D]
L1 = {Fr(k): mp.mpc(mp.mpf(v[0]), mp.mpf(v[1])) for k, v in D.items()}
L0 = {Fr(k): mp.mpc(mp.mpf(v[2]), mp.mpf(v[3])) for k, v in D.items()}
def Jv(y, hat=False):
    j1, j0 = J(y.numerator, y.denominator)
    if hat and y.denominator % 2: j1, j0 = 3*j1, 3*j0
    return j1, j0
cands = {"rho/2": lambda r: r/2, "(rho+1)/2": lambda r: (r+1)/2, "rho": lambda r: r, "2 rho": lambda r: 2*r, "rho/4": lambda r: r/4,
         "-rho/2": lambda r: -r/2, "-(rho+1)/2": lambda r: -(r+1)/2, "-rho": lambda r: -r, "-2rho": lambda r: -2*r}
for nm, f in cands.items():
    for hat in (False, True):
        rows = []
        for r in rhos:
            j1, j0 = Jv(f(r), hat)
            rows.append([mp.mpc(1), j1, j0])
        M = mp.matrix(rows)
        out = []
        for tgt in (L1, L0):
            b = mp.matrix([tgt[r] for r in rhos])
            x, res = mp.qr_solve(M, b)
            out.append((x, res/mp.norm(b)))
        print("%-12s hat=%d  rel residual L1: %s  L0: %s" % (nm, hat, mp.nstr(out[0][1], 3), mp.nstr(out[1][1], 3)))
        if out[0][1] < 1e-6:
            print("     L1 coefficients (1, J1, J0):", [mp.nstr(v, 12) for v in out[0][0]])
            print("     L0 coefficients (1, J1, J0):", [mp.nstr(v, 12) for v in out[1][0]])
