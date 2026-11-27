"""M1, M0 of e4_decay on a circle |qt| = e^(-pi b): X0 = j/10, j = 0..19, then Fourier coefficients in qt."""
import sys, os, json; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from mpmath import mp, mpf, mpc, exp, pi, nstr
import e4_decay
HERE = os.path.dirname(os.path.abspath(__file__))
b = Fr(sys.argv[1]) if len(sys.argv) > 1 else Fr(2, 5)
path = os.path.join(HERE, f"e4_circle_{b.numerator}_{b.denominator}.json")
res = json.load(open(path)) if os.path.exists(path) else {}
for j in range(20):
    key = str(Fr(j, 10))
    if key not in res:
        res[key] = e4_decay.point(Fr(j, 10), b)
        json.dump(res, open(path, "w"), indent=0)
mp.dps = 30
print("b =", b, " largest extrapolation uncertainty:", max(x[2] for v in res.values() for x in v.values()))
r = exp(-pi*mpf(b.numerator)/b.denominator)
for nm in ("M1", "M0"):
    vals = [mpc(mpf(res[str(Fr(j, 10))][nm][0]), mpf(res[str(Fr(j, 10))][nm][1])) for j in range(20)]
    # qt_j = r e^(-pi i j/10);  M = sum c_n qt^n  ->  c_n r^n = (1/20) sum_j M_j e^(+pi i j n/10)
    print(nm, "Taylor coefficients in qt (n = 0..11):")
    for n in range(12):
        c = sum(v*exp(1j*pi*j*n/10) for j, v in enumerate(vals))/20/r**n
        print(f"   n={n:2d}  {nstr(c.real, 16):>22s}  {nstr(c.imag, 16):>22s} i")
