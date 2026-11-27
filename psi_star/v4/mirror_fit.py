"""Mirror law at q = -1, analysis of mirror_limits.json.
 L_j(Z) = lim psi_j*(-e^(2 pi i/Z)) sqrt(t)  is fitted, over all classes of Z, by
     e^(i pi Z/120) * (c1 + c2 U1 + c3 U2)  +  e^(-11 i pi Z/120) * (c4 + c5 U1 + c6 U2)
 with (U1, U2) = (psi1*, psi0*) or (W, W3) at a dual root e^(2 pi i d(Z)), for several candidate maps d."""
import sys, os, json; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from fractions import Fraction as Fr
from jval import J
from wval import W as Wmp
mp.mp.dps = 30
D = json.load(open(sys.argv[1] if len(sys.argv) > 1 else "mirror_limits.json"))
data = []
for key, v in D.items():
    rho, M0 = key.split("|"); Z = Fr(rho) + int(M0)
    data.append((Z, mp.mpc(mp.mpf(v[0]), mp.mpf(v[1])), mp.mpc(mp.mpf(v[2]), mp.mpf(v[3]))))
def W3mp(p, n):
    # W3 = (G(-1/q) - 1 + W)/q  is avoided; direct finite sum with the same rule as W
    from math import gcd
    g = gcd(p, n); p //= g; n //= g
    xi = mp.expjpi(mp.mpf(2*p)/n)
    last = (n+1)//2 if n % 2 else n//2
    s = mp.mpc(0); P = mp.mpc(1)
    for j in range(last):
        if j: P *= 1 - xi**(2*j-1)
        s += (-1)**j*xi**(j*j+3*j)*P
    return s/3 if n % 2 == 0 else s
fam = {"psi*": lambda y: J(y.numerator, y.denominator), "W": lambda y: (Wmp(y.numerator, y.denominator), W3mp(y.numerator, y.denominator))}
maps = {"1/2 - Z/4": lambda Z: Fr(1, 2) - Z/4, "1/2 + Z/4": lambda Z: Fr(1, 2) + Z/4, "-Z/4": lambda Z: -Z/4, "Z/4": lambda Z: Z/4,
        "-Z/2": lambda Z: -Z/2, "Z/2": lambda Z: Z/2, "-Z/8": lambda Z: -Z/8, "Z/8": lambda Z: Z/8}
best = None
for fn, f in fam.items():
    for mn, d in maps.items():
        rows = []
        for Z, l1, l0 in data:
            u1, u2 = f(d(Z))
            Zf = mp.mpf(Z.numerator)/Z.denominator
            p1, p2 = mp.expjpi(Zf/120), mp.expjpi(-11*Zf/120)
            rows.append([p1, p1*u1, p1*u2, p2, p2*u1, p2*u2])
        M = mp.matrix(rows)
        res = []
        sols = []
        for idx in (1, 2):
            b = mp.matrix([x[idx] for x in data])
            x, r = mp.qr_solve(M, b)
            res.append(r/mp.norm(b)); sols.append(x)
        print("%-5s at e^(2 pi i (%s)):  relative residual psi1*: %s  psi0*: %s" % (fn, mn, mp.nstr(res[0], 3), mp.nstr(res[1], 3)))
        if max(res) < 1e-5:
            best = (fn, mn, sols)
if best:
    fn, mn, sols = best
    names = ["1", "U1", "U2"]
    for idx, tgt in enumerate(("psi1*", "psi0*")):
        x = sols[idx]
        print("  %s:  e^(i pi Z/120) [ %s ]" % (tgt, " + ".join("(%s) %s" % (mp.nstr(x[i], 10), names[i]) for i in range(3))))
        print("         + e^(-11 i pi Z/120) [ %s ]" % " + ".join("(%s) %s" % (mp.nstr(x[3+i], 10), names[i]) for i in range(3)))
    print("  with (U1, U2) =", "(psi1*, psi0*)" if fn == "psi*" else "(W, W3)", "at e^(2 pi i (%s))" % mn)

if best and best[0] == "psi*":
    kap = mp.sqrt(mp.pi/mp.sqrt(5)); phi = (1+mp.sqrt(5))/2
    x1, x0 = best[2]
    exp1 = [0, 0, -kap/phi, 0, 1j*kap*phi, 0]
    exp0 = [0, 0, kap*phi, 0, -1j*kap/phi, 0]
    dev = max(max(abs(x1[i]-exp1[i]) for i in range(6)), max(abs(x0[i]-exp0[i]) for i in range(6)))
    print("  closed form: psi1* -> kappa [ -(1/phi) e^(i pi Z/120) psi0*(qh) + i phi e^(-11 i pi Z/120) psi1*(qh) ],")
    print("               psi0* -> kappa [   phi   e^(i pi Z/120) psi0*(qh) - (i/phi) e^(-11 i pi Z/120) psi1*(qh) ],  kappa = sqrt(pi/sqrt5), qh = -e^(-i pi Z/2)")
    print("  largest deviation of the 12 fitted coefficients from this:", mp.nstr(dev, 3))
