"""Forward direction inside the disk.  X = X0 + 2 m - i b,  q = e^(2 pi i/X),  t = -2 pi i/X,  dual nome qt = e^(-pi i X) fixed.
Normalised values  g = psi * sqrt(t) e^(Sbar/t) / Cbar  (psi0*, f0* also divided by e^(2 pi i/3)) are extrapolated to m -> infinity
and written to e2_forward.json.  Sbar = pi^2/3 + i Cl2(pi/3),  Cbar = sqrt(pi/sqrt3) e^(5 pi i/12)."""
import sys, os, json, math; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from flint import acb, arb, ctx
from mpmath import mp
import horo
from e1_limits import tomp, extrap
MS = [400, 460, 530, 610, 700, 800, 920, 1060, 1220, 1400, 1600, 1840]

def point(X0, b):
    rows = []
    for m in MS:
        Xp = (X0 + 2*m, b); K = float(X0 + 2*m)
        r = horo.fwd(Xp, 300 + int(0.25*K))
        X = horo.cx(*Xp); I = acb(0, 1); pi = acb.pi()
        t = -2*pi*I/X
        Sbar = pi*pi/3 + I*CL2()
        Cbar = (pi/acb(3).sqrt()).sqrt()*(5*pi*I/12).exp()
        norm = t.sqrt()*(Sbar/t).exp()/Cbar
        w3 = (2*pi*I/3).exp()
        rows.append((tomp(t), [tomp(r["psi1"]*norm), tomp(r["psi0"]*norm/w3), tomp(r["f1"]*norm), tomp(r["f0"]*norm/w3)]))
    return extrap(rows, ["psi1", "psi0", "f1", "f0"])

def CL2():
    # Cl2(pi/3) = Im Li2(e^(i pi/3))
    I = acb(0, 1)
    return acb((acb.pi()*I/3).exp().polylog(2).imag)

if __name__ == "__main__":
    pts = [(Fr(a), Fr(b)) for b in ("2/5", "7/10") for a in ("0", "3/10", "7/10", "11/10", "8/5", "19/10")]
    if len(sys.argv) > 1: pts = pts[:int(sys.argv[1])]
    res = {}
    for X0, b in pts:
        res[f"{X0}|{b}"] = point(X0, b)
        print(X0, b, {k: v[2] for k, v in res[f"{X0}|{b}"].items()}, flush=True)
        json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "e2_forward.json"), "w"), indent=0)
