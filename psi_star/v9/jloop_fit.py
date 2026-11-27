"""Second step of jloop.py with the exact exponents -1/120 and 71/120 and a constant in the basis.
   L(Y) = q'^(-1/120) A_0(q') + q'^(71/120) A_1(q'),   q'^a := exp(-pi i Y a/2),  Y = Y0 - i b."""
import sys, os, json; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from mpmath import mp, mpf, mpc, nstr, exp, pi, matrix, qr_solve, norm, sqrt
import jloop
mp.dps = 30
res = json.load(open(os.path.join(jloop.HERE, "jloop.json")))
val = lambda Y0: mpc(mpf(res[str(Y0)][0]), mpf(res[str(Y0)][1]))
bases = [Fr(0), Fr(1, 2), Fr(11, 10), Fr(17, 10), Fr(12, 5), Fr(33, 10)]
ex = [-mpf(1)/120, mpf(71)/120]
amps = []
for Y0p in bases:
    A = matrix(6, 2); y = matrix(6, 1)
    for n in range(6):
        Y = mpc(mpf((Y0p + 4*n).numerator)/(Y0p + 4*n).denominator, -mpf(4)/5)
        for j in range(2): A[n, j] = exp(-1j*pi*Y*ex[j]/2)
        y[n] = val(Y0p + 4*n)
    x, rs = qr_solve(A, y)
    Y = mpc(mpf(Y0p.numerator)/Y0p.denominator, -mpf(4)/5)
    qd = exp(-1j*pi*Y/2); d = jloop.duals(qd); d["1"] = mpc(1); d["F1*"] = d["q F1*"]/qd
    amps.append((d, x, rs))
print("two-sector fit, largest residual:", nstr(max(a[2] for a in amps), 3))
names = ["F0*", "1", "F1*", "G0", "G1"]
for j, lab in ((0, "A_0 (sector q'^(-1/120))"), (1, "A_1 (sector q'^(71/120))")):
    A = matrix(6, 5); y = matrix(6, 1)
    for i, (d, x, _) in enumerate(amps):
        for c, nm in enumerate(names): A[i, c] = d[nm]
        y[i] = x[j]
    sol, rs = qr_solve(A, y)
    print(lab, " relative residual", nstr(rs/norm(y), 3))
    for c, nm in zip(sol, names): print(f"      {nm:4s} {nstr(c, 20)}")
c = sqrt(2*pi/sqrt(5)); phi = (1 + sqrt(5))/2
print("c = sqrt(2 pi/sqrt5) =", nstr(c, 20), "  c/2 =", nstr(c/2, 20), "  c phi =", nstr(c*phi, 20), "  c/phi =", nstr(c/phi, 20))
