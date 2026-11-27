"""Test 1 (addendum v12, Theorem C and Lemma B).
 (a) I(1/(2b)) = c p^(1/8) B2,  I(-1/(2b)) = c p^(1/8) B4,   c = 2^(-1/2) e^(-i pi/4) q^(-1/24) qt^(1/24);
     and c p^(1/8) = kappa/(2P) (kappa of v10, P the prefactor of v11).
 (b) general u:  I(u) = c e^(-pi i u^2/2) D(u).
 (c) recursion  I(u - m/b) = p^(2m-2) I(u-(m-1)/b) + p^(2m-3) I(u-(m-2)/b)  at u = 1/(2b).
Usage: python3 t1_kernel_free.py [digits]"""
import sys
from mpmath import mp, mpf, mpc, exp, pi, nstr
from common import Setup, POINTS, I

if __name__ == "__main__":
    mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    for x in POINTS:
        S = Setup(x); b = S.b; p = S.p
        B = S.B(); c = S.cconst(); p8 = S.ppow(mpf(1)/8)
        print("x =", nstr(x, 8), "  |q| =", nstr(abs(S.q), 4), " |p| =", nstr(abs(p), 4))
        print("   |c p^(1/8) / (kappa/(2P)) - 1| =", nstr(abs(c*p8/(S.kappa()/(2*S.P())) - 1), 3))
        # for m >= 1 the line is moved to L + i(m-1/2)/b (no poles crossed), which centres the Gaussian
        Im = {m: S.Iu((mpf(1)/2 - m)/b, y0=(I*(m - mpf(1)/2)/b if m >= 1 else 0)) for m in (-1, 0, 1, 2)}
        print("   I(1/(2b))          =", nstr(Im[0], 20), "  |I/(c p^(1/8) B2) - 1| =", nstr(abs(Im[0]/(c*p8*B[1]) - 1), 3))
        I1 = S.Iu(-1/(2*b))        # on L itself, for (a)
        print("   I(-1/(2b))         =", nstr(I1, 20), "  |I/(c p^(1/8) B4) - 1| =", nstr(abs(I1/(c*p8*B[3]) - 1), 3))
        for m in (1, 2):
            rhs = p**(2*m - 2)*Im[m - 1] + p**(2*m - 3)*Im[m - 2]
            print(f"   recursion m={m}:  |I_m/(p^(2m-2) I_(m-1) + p^(2m-3) I_(m-2)) - 1| =", nstr(abs(Im[m]/rhs - 1), 3))
        for u in (mpc("0.37", "-0.21"), mpc("-0.8", "0.45"), b + 1/(2*b), b - 1/(2*b)):
            Iu = S.Iu(u); rhs = c*exp(-pi*I*u*u/2)*S.D(u)
            print("   u =", nstr(u, 6), "  I(u) =", nstr(Iu, 20), "  |I/(c e^(-pi i u^2/2) D(u)) - 1| =", nstr(abs(Iu/rhs - 1), 3))
        print("   D(b+1/(2b))/B1 - 1 =", nstr(abs(S.D(b + 1/(2*b))/B[0] - 1), 3), "   D(b-1/(2b))/B3 - 1 =", nstr(abs(S.D(b - 1/(2*b))/B[2] - 1), 3),
              "   D(1/(2b))/B2 - 1 =", nstr(abs(S.D(1/(2*b))/B[1] - 1), 3), "   D(-1/(2b))/B4 - 1 =", nstr(abs(S.D(-1/(2*b))/B[3] - 1), 3))
        sys.stdout.flush()
