"""Horocycle limits for the W family.  Y = Y0 + 120 k - i b,  q = e^(2 pi i/Y),  tau = -2 pi i/Y,  q' = e^(-pi i Y) fixed.
For each point (Y0, b) extrapolate to k -> infinity (polynomial extrapolation in tau):
   Pp = P_+(q),  Pm = P_-(q),
   NW = (W(q) - Phi0(tau)) sqrt(tau),  NW3 = (W3(q) - Phi0_3(tau)) sqrt(tau),
with Phi0, Phi0_3 the perturbative series of W, W3 truncated at 50 terms.  Writes e1_limits.json."""
import sys, os, json, math; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from flint import acb, arb, fmpq, ctx
from mpmath import mp, mpf, mpc
import horo, pert
mp.dps = 45
KS = [10, 12, 14, 17, 20, 24, 28, 33, 39, 46]
NT = 50
ser = pert.series()

def tomp(z):
    z = z.mid(); return mpc(mpf(z.real.str(50, radius=False)), mpf(z.imag.str(50, radius=False)))

def phi0(coeffs, tau):
    s = acb(0)
    for k in range(NT, -1, -1):
        f = Fr(coeffs[k]); s = s*tau + arb(fmpq(f.numerator, f.denominator))
    return s

def neville(xs, ys):
    """value at x = 0 of the interpolating polynomial; also the value without the last node (error estimate)"""
    def val(xs, ys):
        p = list(ys); n = len(xs)
        for m in range(1, n):
            for i in range(n - m):
                p[i] = (xs[i + m]*p[i] - xs[i]*p[i + 1])/(xs[i + m] - xs[i])
        return p[0]
    a = val(xs, ys); b = val(xs[1:], ys[1:])
    return a, abs(a - b)

def extrap(rows, names):
    """limit (Neville, with error estimate) and the coefficient of tau^1 of the interpolating polynomial"""
    out = {}
    xs = [t for t, _ in rows]
    for i, name in enumerate(names):
        ys = [y[i] for _, y in rows]
        v, err = neville(xs, ys)
        mp.dps = 120
        n = len(xs); V = mp.matrix(n, n)
        for r in range(n):
            for c in range(n): V[r, c] = xs[r]**c
        co = mp.lu_solve(V, mp.matrix(ys))
        mp.dps = 45
        out[name] = [mp.nstr(v.real, 40), mp.nstr(v.imag, 40), float(err), mp.nstr(co[1].real, 30), mp.nstr(co[1].imag, 30)]
    return out

def point(Y0, b):
    rows = []
    for k in KS:
        Yp = (Y0 + 120*k, b)
        K = float(Y0 + 120*k)
        prec = int(0.0808*K/0.6932*1.15) + 200
        r = horo.near(Yp, prec)
        Y = horo.cx(*Yp); tau = -2*acb.pi()*acb(0, 1)/Y; st = tau.sqrt()
        NW = (r["W"] - phi0(ser["W"], tau))*st
        NW3 = (r["W3"] - phi0(ser["W3"], tau))*st
        rows.append((tomp(tau), [tomp(r["Pp"]), tomp(r["Pm"]), tomp(NW), tomp(NW3)]))
    return extrap(rows, ["Pp", "Pm", "NW", "NW3"])

if __name__ == "__main__":
    pts = [(Fr(a), Fr(b)) for b in ("2/5", "7/10") for a in ("0", "3/10", "7/10", "11/10", "8/5", "37/10")]
    if len(sys.argv) > 1: pts = pts[:int(sys.argv[1])]
    res = {}
    for Y0, b in pts:
        res[f"{Y0}|{b}"] = point(Y0, b)
        print(Y0, b, {k: v[2] for k, v in res[f"{Y0}|{b}"].items()}, flush=True)
        json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "e1_limits.json"), "w"), indent=0)
