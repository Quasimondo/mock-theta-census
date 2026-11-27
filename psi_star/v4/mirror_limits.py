"""Mirror law at q = -1: psi1*, psi0* at x = 1/2 + 1/Z, i.e. q = -e^(2 pi i/Z), for Z -> infinity.
 M_j(Z) = psi_j*(-e^(2 pi i/Z)) * sqrt(t), t = -2 pi i/Z.  For Z = Z0 + 240 k the limit is taken by Neville
 extrapolation in Z^(-1/2).  Writes mirror_limits.json for Z0 = rho + M0, rho in a list, M0 = 0..15."""
import sys, os, json; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from fractions import Fraction as Fr
from jflint import J_both
mp.mp.dps = 40
KS = [5, 7, 10, 14, 20, 28, 40]
NM0 = 16
def Jm(x):
    p, n = x.numerator, x.denominator
    Zapprox = 1/abs(float(x) - 0.5)
    prec = int(1.1*Zapprox) + 400
    while True:                                   # raise the precision until the ball is tight
        a, b = J_both(p, n, prec)
        if float(abs(a).rad()) < 1e-30 and float(abs(b).rad()) < 1e-30: break
        prec *= 2
    f = lambda z: mp.mpc(mp.mpf(z.real.mid().str(40, radius=False)), mp.mpf(z.imag.mid().str(40, radius=False)))
    return f(a), f(b)
def neville(xs, vs):
    n = len(xs); T = vs[:]
    for k in range(1, n):
        T = [(xs[i+k]*T[i] - xs[i]*T[i+1])/(xs[i+k]-xs[i]) for i in range(n-k)]
    return T[0]
def limit(Z0):
    xs, v1, v0 = [], [], []
    for k in KS:
        Z = Fr(Z0) + 240*k
        x = Fr(1, 2) + 1/Z
        j1, j0 = Jm(x)
        t = -2j*mp.pi*Z.denominator/Z.numerator
        xs.append(mp.sqrt(mp.mpf(Z.denominator)/Z.numerator)); v1.append(j1*mp.sqrt(t)); v0.append(j0*mp.sqrt(t))
    return neville(xs, v1), neville(xs, v0), abs(neville(xs, v1) - neville(xs[1:], v1[1:]))
if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "probe"
    out = {}
    if which == "probe":
        for Z0 in (0, 8, 16, 120, 2, 10, 4, 12, 1, 9, 3, 5, 7, 6):
            a, b, err = limit(Z0)
            ph = mp.expjpi(-mp.mpf(Z0)/120)
            print("Z0 = %-4d  L1 e^(-i pi Z0/120) = %s   L0 e^(-i pi Z0/120) = %s   (unc. %s)" % (Z0, mp.nstr(a*ph, 10), mp.nstr(b*ph, 10), mp.nstr(err, 2)), flush=True)
    else:
        rhos = [Fr(x) for x in which.split(",")]
        for rho in rhos:
            for M0 in range(NM0):
                a, b, err = limit(rho + M0)
                out[f"{rho}|{M0}"] = [str(a.real), str(a.imag), str(b.real), str(b.imag), float(err)]
            print(rho, "done; extrapolation uncertainty about", max(out[f"{rho}|{M0}"][4] for M0 in range(NM0)), flush=True)
            json.dump(out, open("mirror_limits.json", "w"))
