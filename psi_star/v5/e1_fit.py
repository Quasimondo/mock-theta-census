"""Fit the horocycle limits of e1_limits.json on functions of the dual nome q' = e^(-pi i Y)."""
import sys, os, json; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from mpmath import mp, mpf, mpc, matrix, lu_solve, qr_solve, pi, sqrt, exp, nstr, norm
import horo
mp.dps = 40
D = json.load(open(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "e1_limits.json")))
def tomp(z):
    z = z.mid(); return mpc(mpf(z.real.str(45, radius=False)), mpf(z.imag.str(45, radius=False)))
pts = []
for key, v in D.items():
    Y0, b = [Fr(x) for x in key.split("|")]
    d = {k: tomp(x) for k, x in horo.dual((Y0, b), 300).items()}
    Y = mpc(mpf(Y0.numerator)/Y0.denominator, -mpf(b.numerator)/b.denominator)
    d["u0"] = exp(-1j*pi*Y*(-mpf(1)/60)); d["u1"] = exp(-1j*pi*Y*(mpf(11)/60))
    d["lim"] = {k: mpc(mpf(x[0]), mpf(x[1])) for k, x in v.items()}; d["key"] = key
    pts.append(d)
print(len(pts), "points; largest extrapolation uncertainty:", max(x[2] for v in D.values() for x in v.values()))

def fit(name, basis, show=True):
    A = matrix(len(pts), len(basis)); y = matrix(len(pts), 1)
    for i, d in enumerate(pts):
        for j, (lab, f) in enumerate(basis): A[i, j] = f(d)
        y[i] = d["lim"][name]
    x, res = qr_solve(A, y)
    rel = res/norm(y)
    if show:
        print(f"{name} on [{', '.join(l for l, _ in basis)}]: relative residual {nstr(rel, 3)}")
        for (lab, _), c in zip(basis, x): print(f"      {lab}: {nstr(c, 22)}")
    return x, rel

B_R = [("q'^(-1/60) R0", lambda d: d["u0"]*d["R0"]), ("q'^(11/60) R1", lambda d: d["u1"]*d["R1"])]
B_4 = [("q'^(-1/60) psi0*", lambda d: d["u0"]*d["psi0"]), ("q'^(11/60) psi1*", lambda d: d["u1"]*d["psi1"]),
       ("q'^(-1/60) f0*", lambda d: d["u0"]*d["f0"]), ("q'^(11/60) f1*", lambda d: d["u1"]*d["f1"])]
if __name__ == "__main__":
    for nm in ("Pp", "Pm"): fit(nm, B_R)
    for nm in ("NW", "NW3"): fit(nm, B_4)
    print("first corrections (coefficient of tau):")
    for d in pts: d["lim"] = {k: mpc(mpf(x[3]), mpf(x[4])) for k, x in D[d["key"]].items()}
    for nm in ("Pp", "Pm"): fit(nm, B_R)
    for nm in ("NW", "NW3"): fit(nm, B_4)
