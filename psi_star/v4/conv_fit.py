"""Converse direction, analysis of conv_limits.json.
 Step 1: for each class (rho, parity of M0):  L(Y) = e^(i pi Y/60) * ( a + b e^(-i pi Y/5) )   (two frequencies)
 Step 2: a and b as constant combinations of 1, psi1*, psi0* at a dual root of unity."""
import sys, os, json; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from fractions import Fraction as Fr
from jval import J
mp.mp.dps = 30
D = json.load(open(sys.argv[1] if len(sys.argv) > 1 else "conv_limits.json"))
cls = {}
for key, v in D.items():
    rho, M0 = key.split("|"); rho = Fr(rho); M0 = int(M0)
    cls.setdefault((rho, M0 % 2), []).append((M0, mp.mpc(mp.mpf(v[0]), mp.mpf(v[1])), mp.mpc(mp.mpf(v[2]), mp.mpf(v[3]))))
AB = {}
worst = 0
for (rho, par), lst in sorted(cls.items()):
    rows = []; r1 = []; r3 = []
    for M0, l1, l3 in lst:
        Y = mp.mpf(M0) + mp.mpf(rho.numerator)/rho.denominator
        ph = mp.expjpi(Y/60)
        rows.append([ph, ph*mp.expjpi(-Y/5)]); r1.append(l1); r3.append(l3)
    M = mp.matrix(rows)
    x1, res1 = mp.qr_solve(M, mp.matrix(r1)); x3, res3 = mp.qr_solve(M, mp.matrix(r3))
    AB[(rho, par)] = (x1[0], x1[1], x3[0], x3[1])
    worst = max(worst, res1, res3)
print("step 1: two-frequency form holds in all", len(AB), "classes; largest residual", mp.nstr(worst, 3))
keys = sorted(AB)
def hat(y):
    j1, j0 = J(y.numerator, y.denominator)
    return j1, j0
points = {"-Y/2": lambda r: -r/2, "Y/2": lambda r: r/2, "-(Y+1)/2": lambda r: -(r+1)/2, "(Y+1)/2": lambda r: (r+1)/2, "-Y": lambda r: -r, "Y": lambda r: r}
names = ["a(W)", "b(W)", "a(W3)", "b(W3)"]
for pn, f in points.items():
    rows = []
    for (rho, par) in keys:
        j1, j0 = hat(f(rho + par))
        rows.append([mp.mpc(1), j1, j0])
    M = mp.matrix(rows)
    line = []
    sols = []
    for i in range(4):
        b = mp.matrix([AB[k][i] for k in keys])
        x, res = mp.qr_solve(M, b)
        line.append(mp.nstr(res/mp.norm(b), 3)); sols.append(x)
    print("step 2, psi* at e^(2 pi i * %s):  relative residuals" % pn, dict(zip(names, line)))
    if all(float(v) < 1e-5 for v in line):
        for i in range(4):
            print("     %s = (%s) + (%s) psi1* + (%s) psi0*" % (names[i], mp.nstr(sols[i][0], 9), mp.nstr(sols[i][1], 9), mp.nstr(sols[i][2], 9)))
# closed forms
c = mp.sqrt(2*mp.pi/mp.sqrt(5)); phi = (1+mp.sqrt(5))/2
dev = 0
for (rho, par) in keys:
    y = -(rho + par)/2
    j1, j0 = J(y.numerator, y.denominator)
    aW, bW, aW3, bW3 = AB[(rho, par)]
    dev = max(dev, abs(aW + c*j0), abs(bW - 1j*c*j1), abs(aW3 + c*phi*j0), abs(bW3 + 1j*c/phi*j1))
print("closed forms  a(W) = -c psi0*,  b(W) = i c psi1*,  a(W3) = -c phi psi0*,  b(W3) = -(i c/phi) psi1*,  c = sqrt(2 pi/sqrt5), phi = golden ratio:")
print("   largest deviation over the", len(keys), "classes:", mp.nstr(dev, 3))
json.dump({f"{k[0]}|{k[1]}": [str(v.real) + "," + str(v.imag) for v in AB[k]] for k in keys}, open("conv_ab.json", "w"))
