"""Test 5 (addendum v12, closing lemma of section 4.1): the two bounds on the lines Re(y/b) = sigma_N = N Im x/|x|^2.
 (a) sup over the line of |1/(-p w; p^2)_inf|  against  |p|^(N^2) / (|p|;|p|^2)_inf^2          (must be <=)
 (b) sup over the part of the line above C of |Phi_b(y) e^(-pi i y^2/2)|  (A = 1/2)  against  exp(-c sigma_N^2),
     c = (pi/2) min( |x|^2/Im x ,  |x| sin(2 delta)/cos^2(delta - beta) )   (the two endpoint values of Im(y^2), see the addendum).
     The lemma says the sup is at most a constant times this; the ratio is printed and must stay bounded in N.
Usage: python3 t5_bounds.py [digits]"""
import sys
from mpmath import mp, mpf, mpc, exp, pi, qp, nstr, arg, tan, sin, cos
from common import Setup, POINTS, I

if __name__ == "__main__":
    mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    delta = mpf(2)/5
    for x in POINTS[:2] + POINTS[3:]:
        S = Setup(x); b = S.b; p = S.p; ap = abs(p)
        print("x =", nstr(x, 8))
        for N in (1, 2, 3, 4, 5, 6):
            sig = N*x.imag/abs(x)**2
            tauC = sig*tan(delta - arg(b))                 # where the line meets the right ray of C
            taus = [tauC + mpf(j)/40 for j in range(0, 1200)]
            supA = max(abs(1/qp(-p*exp(2*pi*(sig + I*(mpf(j)/200))), p*p)) for j in range(200))   # periodic in tau with period 1
            bndA = ap**(N*N)/qp(ap, ap*ap)**2
            supB = max(abs(S.Phi(b*(sig + I*t))*exp(-pi*I*(b*(sig + I*t))**2/2)) for t in taus)
            cc = pi/2*min(abs(x)**2/x.imag, abs(x)*sin(2*delta)/cos(delta - arg(b))**2); bndB = exp(-cc*sig**2)
            print(f"   N={N}:  (a) sup = {nstr(supA, 4)}  bound = {nstr(bndA, 4)}  ok = {supA <= bndA}     (b) sup = {nstr(supB, 4)}  exp(-c sigma^2) = {nstr(bndB, 4)}  ratio = {nstr(supB/bndB, 4)}")
        sys.stdout.flush()
