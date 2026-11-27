"""Converse direction: W(e^(2 pi i/Y)) and W3(e^(2 pi i/Y)) for Y -> infinity.
 M(Y) = W(e^(2 pi i/Y)) * sqrt(tau), tau = -2 pi i/Y.  For Y = Y0 + 120 k the limit k -> infinity is taken by Neville
 extrapolation in Y^(-1/2).  Writes conv_limits.json: for each class "rho|M0" the limits for W and W3.
Usage: python3 conv_limits.py            (8 values of rho, Y up to about 10,800; 1.5 minutes)
       python3 conv_limits.py 0,1/2 long (chosen rho, Y up to about 21,700; writes conv_limits_hi.json)"""
import sys, os, json; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from fractions import Fraction as Fr
from wflint import W_both
mp.mp.dps = 40
KS = [8, 11, 16, 23, 32, 45, 64, 90]
if len(sys.argv) > 2 and sys.argv[2] == "long":
    KS = [8, 11, 16, 23, 32, 45, 64, 90, 128, 181]
def Wm(p, n):
    prec = int(0.26*n/max(p, 1)) + 300
    a, b = W_both(p, n, prec)
    f = lambda z: mp.mpc(mp.mpf(z.real.mid().str(40, radius=False)), mp.mpf(z.imag.mid().str(40, radius=False)))
    return f(a), f(b)
def neville(xs, vs):
    n = len(xs); T = vs[:]
    for k in range(1, n):
        T = [(xs[i+k]*T[i] - xs[i]*T[i+1])/(xs[i+k]-xs[i]) for i in range(n-k)]
    return T[0]
def limit(rho, M0):
    xs, v1, v3 = [], [], []
    for k in KS:
        Y = Fr(M0 + 120*k) + rho
        w, w3 = Wm(Y.denominator, Y.numerator)             # e^(2 pi i/Y) = e^(2 pi i den/num)
        tau = -2j*mp.pi*Y.denominator/Y.numerator
        xs.append(mp.sqrt(mp.mpf(Y.denominator)/Y.numerator)); v1.append(w*mp.sqrt(tau)); v3.append(w3*mp.sqrt(tau))
    return neville(xs, v1), neville(xs, v3), abs(neville(xs, v1) - neville(xs[1:], v1[1:]))
if __name__ == "__main__":
    rhos = [Fr(0), Fr(1, 2), Fr(1, 3), Fr(2, 3), Fr(1, 4), Fr(3, 4), Fr(1, 5), Fr(2, 5)]
    outname = "conv_limits.json"
    if len(sys.argv) > 1 and sys.argv[1] != "all":
        rhos = [Fr(x) for x in sys.argv[1].split(",")]; outname = "conv_limits_hi.json"
    out = {}
    for rho in rhos:
        for M0 in range(10):
            a, b, err = limit(rho, M0)
            out[f"{rho}|{M0}"] = [str(a.real), str(a.imag), str(b.real), str(b.imag), float(err)]
        print(rho, "done; extrapolation uncertainty about", max(out[f"{rho}|{M0}"][4] for M0 in range(10)), flush=True)
    json.dump(out, open(outname, "w"))
