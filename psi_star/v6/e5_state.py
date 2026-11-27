"""The four bilinear combinations ("state integrals") along horocycles.
X = X0 + 2 m - i b,  q = e^(2 pi i/X),  t = -2 pi i/X,  p = e^(-pi i X) fixed.
   B1 = N0(q) Ep'(p) + i q^(-1/4) Y(q) Ep(p)        B3 = N0(q) Em'(p) - i q^(-1/4) Y(q) Em(p)
   B2 = X(q) Ep(p)  + i q^(1/4) N1(q) Ep'(p)        B4 = X(q) Em(p)  - i q^(1/4) N1(q) Em'(p)
Claim: B1 e^(-i Cl2/t), B2 e^(-i Cl2/t), B3 e^(+i Cl2/t), B4 e^(+i Cl2/t) tend to constants that do not depend on p,
with corrections in powers of t whose coefficients do not depend on p either.   Cl2 = Cl2(pi/3)."""
import sys, os, json; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from flint import acb, ctx
from mpmath import mp, mpf, mpc, nstr
import horo, hfam
from e1_limits import extrap
HERE = os.path.dirname(os.path.abspath(__file__))
MS = [200, 225, 255, 290, 330, 375, 425, 480, 500]
def tomp(z):
    z = z.mid(); return mpc(mpf(z.real.str(60, radius=False)), mpf(z.imag.str(60, radius=False)))

def point(X0, b):
    rows = []
    I = acb(0, 1)
    for m in MS:
        Xp = (X0 + 2*m, b); K = float(X0 + 2*m)
        r = hfam.qside(Xp, int(0.9*K) + 400, target=int(0.5*K) + 150)
        Xc = horo.cx(*Xp); pi = acb.pi()
        t = -2*pi*I/Xc
        d = hfam.dside((-pi*I*Xc).exp())
        q4 = (-t/4).exp()
        cl2 = acb((pi*I/3).exp().polylog(2).imag)
        up = (I*cl2/t).exp(); dn = (-I*cl2/t).exp()
        B1 = r["N0"]*d["Epp"] + I/q4*r["Y"]*d["Ep"]
        B3 = r["N0"]*d["Emp"] - I/q4*r["Y"]*d["Em"]
        B2 = r["X"]*d["Ep"] + I*q4*r["N1"]*d["Epp"]
        B4 = r["X"]*d["Em"] - I*q4*r["N1"]*d["Emp"]
        rows.append((tomp(t), [tomp(B1*dn), tomp(B2*dn), tomp(B3*up), tomp(B4*up)]))
    return extrap(rows, ["B1", "B2", "B3", "B4"])

if __name__ == "__main__":
    pts = [(Fr(a), Fr(b)) for a, b in (("0", "2/5"), ("3/10", "2/5"), ("11/10", "2/5"), ("17/10", "2/5"), ("7/10", "7/10"), ("13/10", "7/10"))]
    path = os.path.join(HERE, "e5_state.json")
    res = json.load(open(path)) if os.path.exists(path) else {}
    for X0, b in pts:
        key = f"{X0}|{b}"
        if key not in res:
            res[key] = point(X0, b); json.dump(res, open(path, "w"), indent=0)
    mp.dps = 25
    print("extrapolation uncertainty:", max(x[2] for v in res.values() for x in v.values()))
    for nm in ("B1", "B2", "B3", "B4"):
        print(nm, " limit and coefficient of t at each point (X0|b):")
        for key, v in res.items():
            print(f"   {key:10s} {v[nm][0][:22]:>24s} {v[nm][1][:22]:>24s} i     c1: {v[nm][3][:18]:>20s} {v[nm][4][:18]:>20s} i")
