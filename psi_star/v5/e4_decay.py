"""The decaying sector of the forward law, by exact cancellation of the growing one.
X = X0 + 2 m - i b,  q = e^(2 pi i/X),  t = -2 pi i/X,  qt = e^(-pi i X) fixed (|qt| = e^(-pi b)).
Growing sector (e2):  psi1*(q) ~ A(t) W(qt),  f1*(q) ~ -2 A(t) Fa(qt)  with the same A(t).  Hence
    D1 = 2 Fa(qt) psi1*(q) + W(qt) f1*(q) - 2 (W + Fa)(qt) Phi(t)
has no growing sector; Phi = perturbative series of psi1* at q -> 1 (that of f1* is 2 Phi), 140 terms from v3.
If the trivial sector is exactly Phi(t) times 1, D1 is the decaying sector:  D1 ~ C t^(-1/2) e^(-S/t) M1(qt) (1 + ...),
S = pi^2/3 - i Cl2(pi/3),  C = sqrt(pi/sqrt3) e^(-5 pi i/12).  Same for psi0*, f0* (D0, M0; extra factor e^(-2 pi i/3)).
Writes e4_decay.json."""
import sys, os, json, math; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from flint import acb, arb, ctx
from mpmath import mp
import horo
from e1_limits import extrap
from mpmath import mpf, mpc
HERE = os.path.dirname(os.path.abspath(__file__))
MS = [200, 225, 255, 290, 330, 375, 425, 480, 500]
NT = 140
V3 = json.load(open(os.path.join(HERE, "q2_series_v3.json")))

def tomp(z):
    z = z.mid(); return mpc(mpf(z.real.str(60, radius=False)), mpf(z.imag.str(60, radius=False)))

def phi(name, t):
    s = acb(0)
    for k in range(NT, -1, -1): s = s*t + arb(V3[name][k][0])
    return s

def point(X0, b, verbose=False):
    rows = []
    for m in MS:
        Xp = (X0 + 2*m, b); K = float(X0 + 2*m)
        target = int(0.47*K) + 160
        prec = target + int(0.25*K) + 200
        r = horo.fwd(Xp, prec, target)
        X = horo.cx(*Xp); I = acb(0, 1); pi = acb.pi()
        qt = (-pi*I*X).exp()
        d = horo.near_q(qt, int(prec*0.6932/(math.pi*float(b))**0.5) + 60)
        t = -2*pi*I/X
        cl2 = acb((pi*I/3).exp().polylog(2).imag)
        S = pi*pi/3 - I*cl2
        C = (pi/acb(3).sqrt()).sqrt()*(-5*pi*I/12).exp()
        norm = t.sqrt()*(S/t).exp()/C
        w3 = (-2*pi*I/3).exp()
        B = d["W"] + d["Fa"]
        D1 = 2*d["Fa"]*r["psi1"] + d["W"]*r["f1"] - 2*B*phi("psi1@1", t)
        D0 = 2*d["Fa"]*r["psi0"] + d["W"]*r["f0"] - 2*B*phi("psi0@1", t)
        if verbose:
            print(f"   K = {K:7.1f}: |psi1*| = {float(abs(r['psi1']).mid()):.3e}  |D1| = {float(abs(D1).mid()):.3e}  |D1 norm| = {float(abs(D1*norm).mid()):.6f}", flush=True)
        rows.append((tomp(t), [tomp(D1*norm), tomp(D0*norm/w3)]))
    return extrap(rows, ["M1", "M0"])

if __name__ == "__main__":
    pts = [(Fr(a), Fr(b)) for b in ("2/5", "7/10") for a in ("0", "3/10", "7/10", "11/10", "8/5", "19/10")]
    n = int(sys.argv[1]) if len(sys.argv) > 1 else len(pts)
    res = {}
    for i, (X0, b) in enumerate(pts[:n]):
        res[f"{X0}|{b}"] = point(X0, b, verbose=(i == 0))
        print(X0, b, {k: (v[0][:22], v[1][:22], v[2]) for k, v in res[f"{X0}|{b}"].items()}, flush=True)
        json.dump(res, open(os.path.join(HERE, "e4_decay.json"), "w"), indent=0)
