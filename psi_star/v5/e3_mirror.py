"""Mirror law inside the disk.  Z = Z0 + 240 k - i b,  q = -e^(2 pi i/Z) = -e^(-t),  t = -2 pi i/Z,  qh = -e^(-pi i Z/2) fixed.
   N1 = (psi1*(q) - Phi_1(t)) sqrt(t),  N0 = (psi0*(q) - Phi_0(t)) sqrt(t),   Phi = perturbative series at q = -1 (55 terms).
Extrapolates to k -> infinity and fits on psi*, f* at qh.   python3 e3_mirror.py [number of points]"""
import sys, os, json; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from flint import acb, arb, fmpq, ctx
from mpmath import mp, mpf, mpc, matrix, qr_solve, nstr, norm, pi, exp, sqrt
import horo, pert
from e1_limits import tomp, extrap
KS = [13, 14, 15, 17, 19, 21, 23, 25, 28, 31]
NT = 55
HERE = os.path.dirname(os.path.abspath(__file__))
ser = pert.series_m1()

def phi(coeffs, t):
    s = acb(0)
    for k in range(NT, -1, -1):
        f = Fr(coeffs[k]); s = s*t + arb(fmpq(f.numerator, f.denominator))
    return s

def point(Z0, b):
    rows = []
    for k in KS:
        Zp = (Z0 + 240*k, b); K = float(Z0 + 240*k)
        r = horo.mfwd(Zp, int(0.2*K) + 250)
        t = -2*acb.pi()*acb(0, 1)/horo.cx(*Zp); st = t.sqrt()
        rows.append((tomp(t), [tomp((r["psi1"] - phi(ser["psi1"], t))*st), tomp((r["psi0"] - phi(ser["psi0"], t))*st)]))
    return extrap(rows, ["N1", "N0"])

if __name__ == "__main__":
    pts = [(Fr(a), Fr(b)) for b in ("1/2", "9/10") for a in ("0", "3/10", "7/10", "8/5")]
    if len(sys.argv) > 1: pts = pts[:int(sys.argv[1])]
    path = os.path.join(HERE, "e3_mirror.json")
    res = json.load(open(path)) if os.path.exists(path) else {}
    for Z0, b in pts:
        key = f"{Z0}|{b}"
        if key not in res:
            res[key] = point(Z0, b)
            json.dump(res, open(path, "w"), indent=0)
        print(key, {k: v[2] for k, v in res[key].items()}, flush=True)
    # fit
    mp.dps = 40
    P = []
    for key, v in res.items():
        Z0, b = [Fr(x) for x in key.split("|")]
        ctx.prec = 400
        qh = -(-acb.pi()*acb(0, 1)*horo.cx(Z0, b)/2).exp()
        d = {k: tomp(x) for k, x in horo.dual_q(qh, 400).items()}
        Z = mpc(mpf(Z0.numerator)/Z0.denominator, -mpf(b.numerator)/b.denominator)
        d["u0"] = exp(1j*pi*Z/120); d["u1"] = exp(-11j*pi*Z/120); d["v"] = v
        P.append(d)
    basis = [("u0 psi0*", lambda d: d["u0"]*d["psi0"]), ("u1 psi1*", lambda d: d["u1"]*d["psi1"]),
             ("u0 f0*", lambda d: d["u0"]*d["f0"]), ("u1 f1*", lambda d: d["u1"]*d["f1"])]
    kap = sqrt(pi/sqrt(5)); phi_ = (1 + sqrt(5))/2
    print("kappa/phi =", nstr(kap/phi_, 20), " kappa*phi =", nstr(kap*phi_, 20), "  [u0 = e^(i pi Z/120), u1 = e^(-11 i pi Z/120)]")
    for which, lab in ((0, "limit"), (3, "coefficient of t")):
        for nm in ("N1", "N0"):
            A = matrix(len(P), 4); y = matrix(len(P), 1)
            for i, d in enumerate(P):
                for j, (_, f) in enumerate(basis): A[i, j] = f(d)
                y[i] = mpc(mpf(d["v"][nm][which]), mpf(d["v"][nm][which + 1]))
            x, r = qr_solve(A, y)
            print(f"{nm} {lab}: relative residual {nstr(r/norm(y), 3)}")
            for (l, _), c in zip(basis, x): print(f"      {l}: {nstr(c, 20)}")
