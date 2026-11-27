"""Addendum v13, step 3: the kernel integral P J against the Borel sum of Phi, on a full turn of directions.
For t = |t| e^(i phi):  x = i t/(2 pi)  (so arg x = phi + 90 degrees; Im x < 0, |q| > 1, for phi < -90),
   P J(x)  by quadrature (two contour angles, to show the accuracy of the quadrature),
   S_theta[Phi](t)  by Pade-Laplace with theta = phi clamped to [-170, -25] degrees (inside the sector free of Borel singularities),
   for two Pade orders.
Usage: python3 s3_sectors.py [digits] [|t|] [case] [phi ...]"""
import sys
from mpmath import mp, mpf, mpc, exp, pi, sqrt, nstr, arg
from borel import PJ, borel_sum, I

if __name__ == "__main__":
    mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 25
    T = mpf(sys.argv[2]) if len(sys.argv) > 2 else mpf("0.3")
    case = sys.argv[3] if len(sys.argv) > 3 else "psi1"
    phis = [int(a) for a in sys.argv[4:]] or [-250, -220, -190, -160, -135, -110, -100, -80, -60, -30, 0, 30, 60]
    print(f"{case}, |t| = {T}, {mp.dps} digits")
    for phi in phis:
        t = T*exp(I*phi*pi/180); x = I*t/(2*pi); beta = arg(sqrt(x))
        dmax = pi/2 - abs(beta)
        v1 = PJ(x, delta=dmax*mpf("0.5"), case=case); v2 = PJ(x, delta=dmax*mpf("0.65"), case=case)
        theta = max(-170, min(-25, phi))*pi/180
        s1 = borel_sum(case, t, theta, 60); s2 = borel_sum(case, t, theta, 90)
        print(f"phi = {phi:5d}  arg x = {phi + 90:5d}  |q| = {nstr(abs(exp(-t)), 5)}   P J = {nstr(v1, 20)}   quadrature: {nstr(abs(v1 - v2), 2)}   "
              f"Borel (theta = {nstr(theta*180/pi, 4)}): Pade 90 vs 60: {nstr(abs(s1 - s2), 2)}   |P J - Borel| = {nstr(abs(v1 - s2), 3)}")
        sys.stdout.flush()
