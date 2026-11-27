"""Fit the forward horocycle limits (e2_forward.json) on W-family functions at the dual nome qt = e^(-pi i X)."""
import sys, os, json; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from flint import acb, ctx
from mpmath import mp, mpf, mpc, matrix, qr_solve, nstr, norm
import horo
from e1_limits import tomp
mp.dps = 40
D = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "e2_forward.json")))
pts = []
for key, v in D.items():
    X0, b = [Fr(x) for x in key.split("|")]
    ctx.prec = 400
    qt = (-acb.pi()*acb(0, 1)*horo.cx(X0, b)).exp()
    d = {k: tomp(x) for k, x in horo.near_q(qt, 200).items()}
    d["one"] = mpc(1)
    d["lim"] = {k: mpc(mpf(x[0]), mpf(x[1])) for k, x in v.items()}
    d["c1"] = {k: mpc(mpf(x[3]), mpf(x[4])) for k, x in v.items()}
    pts.append(d)
print(len(pts), "points; largest extrapolation uncertainty:", max(x[2] for v in D.values() for x in v.values()))
names = ["W", "W3", "Fa", "Fb", "one"]
def fit(name, which="lim", names=names):
    A = matrix(len(pts), len(names)); y = matrix(len(pts), 1)
    for i, d in enumerate(pts):
        for j, nm in enumerate(names): A[i, j] = d[nm]
        y[i] = d[which][name]
    x, res = qr_solve(A, y)
    print(f"{name} ({which}) on {names}: relative residual {nstr(res/norm(y), 3)}")
    print("     ", "  ".join(f"{nm}: {nstr(c, 15)}" for nm, c in zip(names, x)))
for nm in ("psi1", "psi0", "f1", "f0"): fit(nm)
print("first corrections (coefficient of t):")
for nm in ("psi1", "psi0", "f1", "f0"): fit(nm, "c1")
