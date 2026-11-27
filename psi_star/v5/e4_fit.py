"""Compare the decaying-sector function M (e4_decay.json, e4_circle_*.json) with 2 Wr(qt),  Wr = q W3 Fa - W Fb."""
import sys, os, json, glob; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from flint import acb, ctx
from mpmath import mp, mpf, mpc, nstr, sqrt
import horo
from e1_limits import tomp
HERE = os.path.dirname(os.path.abspath(__file__))
mp.dps = 30
data = []
for key, v in json.load(open(os.path.join(HERE, "e4_decay.json"))).items():
    X0, b = [Fr(x) for x in key.split("|")]; data.append((X0, b, v))
for path in sorted(glob.glob(os.path.join(HERE, "e4_circle_*.json"))):
    num, den = os.path.basename(path)[:-5].split("_")[2:]
    for key, v in json.load(open(path)).items(): data.append((Fr(key), Fr(int(num), int(den)), v))
worst = 0; r1 = []; r0 = []
for X0, b, v in data:
    ctx.prec = 300
    qt = (-acb.pi()*acb(0, 1)*horo.cx(X0, b)).exp()
    d = {k: tomp(x) for k, x in horo.near_q(qt, 150).items()}
    Wr2 = 2*(d["q"]*d["W3"]*d["Fa"] - d["W"]*d["Fb"])
    for nm in ("M1", "M0"):
        M = mpc(mpf(v[nm][0]), mpf(v[nm][1])); worst = max(worst, abs(M - Wr2)/abs(Wr2))
    r1.append(mpc(mpf(v["M1"][3]), mpf(v["M1"][4]))/Wr2); r0.append(mpc(mpf(v["M0"][3]), mpf(v["M0"][4]))/Wr2)
print(len(data), "points (|qt| = e^(-pi b), b in", sorted({str(b) for _, b, _ in data}), ")")
print("largest relative deviation |M - 2 Wr| / |2 Wr| (M1 and M0):", nstr(worst, 3))
for lab, r, ex in (("psi1*", r1, mpf(1)/4 - 1j/(18*sqrt(3))), ("psi0*", r0, -mpf(1)/12 - 1j/(18*sqrt(3)))):
    m = sum(r)/len(r)
    print(f"first correction / (2 Wr) for {lab}: mean {nstr(m, 12)}, spread {nstr(max(abs(x - m) for x in r), 3)};  conjugate of the growing-sector value: {nstr(ex, 12)}")
